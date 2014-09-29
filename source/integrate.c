#include "../include/standard_include.h"

#include "../include/typedefs.h"
#include "../include/proto_integrate.h"
#include "../include/proto_ext.h"
#include "../include/proto_communicate_wrappers.h"
#include "../include/defines.h"

void chain_force(CHAIN* chain, BEADS *beads,int myid,MPI_Comm world,int numprocs)
{
#include "../include/typ_mask.h"
  double omp2;
  double my_vchain;
  int ibead;

  omp2 = chain->omp*chain->omp;
  if(myid == 0){
   for(ibead=2;ibead <= beads->my_np;ibead++) {
    beads->fu[ibead].x = -beads->m_stage[ibead]*omp2*beads->u[ibead].x;
   } /* endfor */
  } else {
   for(ibead=1;ibead <= beads->my_np;ibead++) {
    beads->fu[ibead].x = -beads->m_stage[ibead]*omp2*beads->u[ibead].x;
   } /* endfor */
  }

/* Zero force on 1st bead only on master node */
  if(myid == 0) {
   beads->fu[1].x = 0.0;
  }

  my_vchain  = 0.0;
  if(myid == 0) {
   for(ibead=2;ibead <= beads->my_np;ibead++) {
    my_vchain += 0.5*beads->m_stage[ibead]*omp2*
                         (beads->u[ibead].x*beads->u[ibead].x);
   }/* endfor */
  } else {
   for(ibead=1;ibead <= beads->my_np;ibead++) {
    my_vchain += 0.5*beads->m_stage[ibead]*omp2*
                         (beads->u[ibead].x*beads->u[ibead].x);
   }/* endfor */
  }/* endif */

/* Collect up all the local vchains */
 
  if(numprocs > 1) {
    Reduce(&my_vchain,&(chain->v_chain),1,MPI_DOUBLE,MPI_SUM,0,world);
  } else {
    chain->v_chain = my_vchain;
  }
}/* end function */
/* -------------------------------------------------------------------------*/

void integrate(CHAIN *chain,BEADS* beads,THERM *therm,T_MATRICES *t_matrices,
               double temp,double dt,
               int trans_opt,int trans_type,
               int myid,int numprocs,MPI_Datatype mpi_vector,
                MPI_Comm world)
{
#include "../include/typ_mask.h"
   double dts,dt2;
   int ibead,ires,icalls;

   dt2 = 0.5*dt;

   for(ires=1;ires <= therm->nres;ires++){
    for(icalls=1;icalls <= therm->ncalls;icalls++){
     dts = therm->dtsuz[icalls];
     suz_int(therm,beads,temp,dts);
    }
   }  

   for(ibead=1;ibead <= beads->my_np; ibead++){
    beads->v[ibead].x += dt2*(beads->fu[ibead].x + beads->fu_ext[ibead].x) 
                       /beads->mp_stage[ibead];
   }/* endfor */

   for(ibead=1;ibead <= beads->my_np; ibead++){
    beads->u[ibead].x += dt*beads->v[ibead].x;
   } /* endfor */

   chain_force(chain,beads,myid,world,numprocs);

  if(trans_opt == 0)
    if(numprocs > 1){
      trans_u_to_x_rec_par(beads,myid,numprocs,mpi_vector,world);
    } else {
      trans_u_to_x_rec(beads);
    } /* endif */
  else
   trans_u_to_x_mat(beads,t_matrices,myid,numprocs,mpi_vector,world);
  ext_force(beads,world,numprocs);
  if(trans_opt == 0)
    if(numprocs > 1){
     trans_fx_to_fu_rec_par(beads,myid,numprocs,mpi_vector,world);
    } else {
     trans_fx_to_fu_rec(beads);
    }
  else
   trans_fx_to_fu_mat(beads,t_matrices,myid,numprocs,trans_type,
                      mpi_vector,world);

   for(ibead=1;ibead <= beads->my_np; ibead++){
    beads->v[ibead].x += dt2*(beads->fu[ibead].x + beads->fu_ext[ibead].x) 
                       /beads->mp_stage[ibead];
   } /* endfor */

   for(ires=1;ires <= therm->nres;ires++){ 
    for(icalls=1;icalls <= therm->ncalls;icalls++){
     dts = therm->dtsuz[icalls];
     suz_int(therm,beads,temp,dts);
    }
   }   

} /* end function */

 void suz_int(THERM *therm,BEADS *beads,double temp,double dts)
{
#include "../include/typ_mask.h"
  double xkin,ckine,ckewant,aa,f1,f2;
  int l,ibead,nch;
  int start_proc = therm->start_proc;
  double rbeta2;
  VECTOR **eta    = therm->eta;
  VECTOR **etadot = therm->etadot;
  VECTOR **feta   = therm->feta;
  double *q       = therm->q;
  VECTOR *v       = beads->v;
  double *mass    = beads->mp_stage;
  nch             = therm->nch;

  rbeta2 = 0.5*temp/RBOLTZ;
  for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
    xkin = 0.5*mass[ibead]*v[ibead].x*v[ibead].x;
    feta[1][ibead].x = 2.0*(xkin - rbeta2);

  }

  for(l=2; l<=nch; l++) {
   for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
     ckine = 0.5*q[ibead]*etadot[l-1][ibead].x*etadot[l-1][ibead].x;
     ckewant = rbeta2;
     feta[l][ibead].x = 2.0*(ckine - ckewant);

    }
  }
   for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
     aa = exp(-0.125*dts*etadot[nch-1][ibead].x);
     etadot[nch][ibead].x = etadot[nch][ibead].x*aa*aa 
                        + 0.25*dts*feta[nch][ibead].x*aa/q[ibead];
     ckine = 0.5*q[ibead]*etadot[nch][ibead].x*etadot[nch][ibead].x;
     feta[nch+1][ibead].x = 2.0*(ckine - ckewant);
     f1 = feta[nch-1][ibead].x;
     f2 = feta[nch+1][ibead].x;
     feta[nch-1][ibead].x = f1 + f2;

   }
/*-----------------------------------------------------------------------*/
  for(l=1; l<=nch-1; l++) {
   for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
      aa = exp(-0.125*dts*etadot[nch+1-l][ibead].x);
      etadot[nch-l][ibead].x = etadot[nch-l][ibead].x*aa*aa + 
                      0.25*dts*feta[nch-l][ibead].x*aa/q[ibead];

   }
 }
/*-----------------------------------------------------------------------*/
   for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
     aa = exp(-0.5*dts*etadot[1][ibead].x);
     v[ibead].x *= aa;

   }

  for(l=1; l<=nch; l++){
   for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
        eta[l][ibead].x += 0.5*dts*etadot[l][ibead].x;
        feta[l][ibead].x = 0.0;

   }
  }

   for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
     xkin = 0.5*mass[ibead]*v[ibead].x*v[ibead].x;
     feta[1][ibead].x = 2.0*(xkin - rbeta2);
     ckine = 0.5*q[ibead]*etadot[nch][ibead].x*etadot[nch][ibead].x;
     feta[nch-1][ibead].x = 2.0*(ckine - ckewant);

  }
/*-----------------------------------------------------------------------*/
     for(l=1; l<=nch-1; l++) {
      for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
       aa = exp(-0.125*dts*etadot[l+1][ibead].x);
       etadot[l][ibead].x = etadot[l][ibead].x*aa*aa 
                          + 0.25*dts*feta[l][ibead].x*aa/q[ibead];
       ckine = 0.5*q[ibead]*etadot[l][ibead].x*etadot[l][ibead].x;
       ckewant = rbeta2;
       feta[l+1][ibead].x += 2.0*(ckine - ckewant);

     }
    }

    for(ibead=start_proc;ibead <= beads->my_np; ibead++){  
      aa = exp(-0.125*dts*etadot[nch-1][ibead].x);
      etadot[nch][ibead].x = etadot[nch][ibead].x*aa*aa 
                           + 0.25*dts*feta[nch][ibead].x*aa/q[ibead];

    }

} /* end function */
