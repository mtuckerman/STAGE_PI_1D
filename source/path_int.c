/* 3D path integral code (serial and parallel) */

#include "../include/standard_include.h"

#include "../include/typedefs.h"
#include "../include/proto.h"
#include "../include/proto_ext.h"
#include "../include/proto_integrate.h"
#include "../include/proto_observables.h"
#include "../include/proto_utils.h"
#include "../include/proto_init.h"
#include "../include/proto_communicate_wrappers.h"
#include "../include/defines.h"


int main(int argc, char *argv[])
{
#include "../include/typ_mask.h"
  BEADS *beads;
  CHAIN *chain;
  THERM *therm;
  SIM_PARAMS *sim_params;
  T_MATRICES *t_matrices;


/* Averages, etc. */
  double etot0,etot,econs=0.0;      /* Total energy */
  double tinst,tbar;                /* Instantaneous and avg. temperature */
  double tcent,tcbar;               /* Instantaneous and avg. temperature */
  double vir,vir_bar;               /* Instantaneous and avg. vir est */
  double prim,prim_bar;             /* Instantaneous and avg. prim est */

/* Coord file information */
  char coord_name[MAXWORD];         /* Input file name */
  FILE *coord_file;                 /* coordinate file pointer */
  FILE *cent_file;                 /* coordinate file pointer */
  int nskip,nwrite_cent;

/* Counter variables */
  int istep,nch1,iii;    
  int ibead;

/* Centroid MD option */
  int cmd_opt;
  double gamma;

/* Random number generator seed */
  double qseed=10985755.0;          /* random number seed */

/* MPI communication variables */
   int numprocs,myid;
   MPI_Datatype mpi_vector;
   MPI_Comm world;

/* -------------------------------------------------------------------------*/
/*  Begin main program  */
/* -------------------------------------------------------------------------*/
/* Initialize communication */

  Init(&argc,&argv,&world);
  Comm_rank(world,&myid);
  Comm_size(world,&numprocs);

/* Define the vector type */
 if(numprocs > 1){
  Type_contiguous(1,MPI_DOUBLE,&mpi_vector);
  Type_commit(&mpi_vector);
 }/* endif */

/* Malloc the memory */
  beads           = (BEADS *) malloc(sizeof(BEADS));
  chain           = (CHAIN *) malloc(sizeof(CHAIN));
  therm           = (THERM *) malloc(sizeof(THERM));
  sim_params      = (SIM_PARAMS *) malloc(sizeof(SIM_PARAMS));
  t_matrices      = (T_MATRICES *) malloc(sizeof(T_MATRICES));

/* Set simulation parameters and broadcast information */
 if(myid == 0) {
  set_sim(sim_params,chain,beads,therm,argv[1],coord_name,numprocs,&nskip,
          &nwrite_cent,&cmd_opt,&gamma);
 }
  if(myid == 0 && cmd_opt == 1){
    cent_file = fopen(argv[2],"w");
  }

  if(numprocs > 1){
   Bcast(&(sim_params->nsteps),1,MPI_INT,0,world);
   Bcast(&(beads->np),1,MPI_INT,0,world);
   Bcast(&(chain->np),1,MPI_INT,0,world);
   Bcast(&(beads->my_np),1,MPI_INT,0,world);
   Bcast(&(chain->my_np),1,MPI_INT,0,world);
   Bcast(&(sim_params->dt),1,MPI_DOUBLE,0,world);
   Bcast(&(sim_params->temp),1,MPI_DOUBLE,0,world);
   Bcast(&(sim_params->temp_p),1,MPI_DOUBLE,0,world);
   Bcast(&(beads->mass),1,MPI_DOUBLE,0,world);
   Bcast(&(beads->mass_const),1,MPI_DOUBLE,0,world);
   Bcast(&(therm->nch),1,MPI_INT,0,world);
   Bcast(&(therm->msuz),1,MPI_INT,0,world);
   Bcast(&(therm->nres),1,MPI_INT,0,world);
   Bcast(&(therm->ncalls),1,MPI_INT,0,world);
   Bcast(&(chain->omp),1,MPI_DOUBLE,0,world);
   Bcast(coord_name,MAXWORD,MPI_CHAR,0,world);
   Bcast(&(sim_params->trans_opt),1,MPI_INT,0,world);
   Bcast(&(sim_params->trans_type),1,MPI_INT,0,world);
   Bcast(&nskip,1,MPI_INT,0,world);
   Bcast(&cmd_opt,1,MPI_INT,0,world);
  }/* endif */

 therm->start_proc = 1;
 if(cmd_opt == 1) therm->start_proc = (myid == 0 ? 2:1);


  beads->u        = (VECTOR *)malloc(chain->my_np*sizeof(VECTOR))-1;
  beads->fu       = (VECTOR *)malloc(chain->my_np*sizeof(VECTOR))-1;
  beads->fu_ext   = (VECTOR *)malloc(chain->my_np*sizeof(VECTOR))-1;
  beads->r        = (VECTOR *)malloc(chain->my_np*sizeof(VECTOR))-1;
  beads->fr       = (VECTOR *)malloc(chain->my_np*sizeof(VECTOR))-1;
  beads->v        = (VECTOR *)malloc(chain->my_np*sizeof(VECTOR))-1;
  beads->m_stage  = (double *)malloc(chain->my_np*sizeof(double))-1;
  beads->mp_stage = (double *)malloc(chain->my_np*sizeof(double))-1;
  therm->q        = (double *)malloc(chain->my_np*sizeof(double))-1;
  therm->eta      = cmall_mat(1,therm->nch,1,chain->my_np);
  therm->etadot   = cmall_mat(1,therm->nch,1,chain->my_np);
  nch1 = therm->nch+1;
  therm->feta     = cmall_mat(1,nch1,1,chain->my_np);
  therm->dtsuz    = (double *)malloc(therm->ncalls*sizeof(double))-1;

/* Initialize Suzuki Yoshida time steps and broadcast */
  if(myid == 0) {
   init_suz(therm,sim_params->dt);
  }
  if(numprocs > 1)
   Bcast(&(therm->dtsuz[1]),therm->ncalls,MPI_DOUBLE,0,
             world);

/* Zero external forces */

  for(ibead=1;ibead <= beads->my_np; ibead++){
   beads->fu_ext[ibead].x = 0.0;
  }

/* If using matrix transformation option, assign matrices */

  if(sim_params->trans_opt == 1){
   t_matrices->u_to_x   = mall_mat(1,beads->my_np,1,beads->np);
   t_matrices->fx_to_fu = mall_mat(1,beads->my_np,1,beads->np);
   if(sim_params->trans_type != 1){

     t_matrices->eig = (double *) malloc(beads->my_np*sizeof(double))-1;
     beads->eig      = (double *) malloc(beads->my_np*sizeof(double))-1;

   }
   assign_t_mats(t_matrices,beads->my_np,beads->np,myid,sim_params->trans_type,
                 world,numprocs);

   if(sim_params->trans_type != 1){
    for(ibead=1;ibead<=beads->my_np;ibead++){
     beads->eig[ibead] = t_matrices->eig[ibead];
    }
   }/* endif */
  }/* endif trans_opt */

/* Read in or set up initial configuration of beads */

   set_coords(beads,therm,coord_name,&qseed,chain->omp,sim_params->temp,
              sim_params->temp_p,
              sim_params->trans_type,world,mpi_vector,numprocs,myid,
              gamma);

   write_coords(coord_file,coord_name,beads,therm,
                myid,numprocs,mpi_vector,world); 

  chain_force(chain,beads,myid,world,numprocs);


  if(sim_params->trans_opt == 0)
    if(numprocs > 1){
     trans_u_to_x_rec_par(beads,myid,numprocs,mpi_vector,world);
    } else {
     trans_u_to_x_rec(beads);
    }
  else
   trans_u_to_x_mat(beads,t_matrices,myid,numprocs,mpi_vector,world);
  ext_force(beads,world,numprocs);
  if(sim_params->trans_opt == 0)
    if(numprocs > 1){
      trans_fx_to_fu_rec_par(beads,myid,numprocs,mpi_vector,world); 
    } else {
      trans_fx_to_fu_rec(beads); 
    } /* endif */
  else
   trans_fx_to_fu_mat(beads,t_matrices,myid,numprocs,sim_params->trans_type,
                      mpi_vector,world); 

  etot0 = fict_kinet(beads,world,numprocs) + chain->v_chain + beads->v_ext
        + therm_energy(therm,chain->my_np,sim_params->temp,world,numprocs); 



  tinst = fict_kinet(beads,world,numprocs);
  tinst *= 2.0*RBOLTZ/(chain->np);
  tbar = tinst;

 
  if(myid == 0){
   tcent = cent_kinet(beads,world,numprocs);
   tcent *= 2.0*RBOLTZ;
   tcbar = tinst;
  }/* endif */

  vir = vir_est(beads,sim_params->temp,myid,world,numprocs);
  prim = prim_est(beads,chain,sim_params->temp,myid,world,numprocs);
  if(myid == 0){
   vir_bar = vir;
   prim_bar = prim;
  }

  for(istep=1;istep <= sim_params->nsteps;istep++){

   integrate(chain,beads,therm,t_matrices,sim_params->temp,sim_params->dt,
             sim_params->trans_opt,sim_params->trans_type,
             myid,numprocs,mpi_vector,world);

   etot = fict_kinet(beads,world,numprocs) + chain->v_chain + beads->v_ext
        + therm_energy(therm,chain->my_np,sim_params->temp,world,numprocs);

   tinst = fict_kinet(beads,world,numprocs);

   if(myid == 0) {
    tinst *= 2.0*RBOLTZ/(chain->np);
    tbar += tinst;
   }

   if(myid == 0) {
    tcent = cent_kinet(beads,world,numprocs);
    tcent *= 2.0*RBOLTZ;
    tcbar += tcent;
   }

    vir = vir_est(beads,sim_params->temp,myid,world,numprocs);
    prim = prim_est(beads,chain,sim_params->temp,myid,world,numprocs);
    if(myid == 0){
      vir_bar += vir;
      prim_bar += prim;
      econs += fabs(etot-etot0)/fabs(etot0);
    }/* endif myid */

    if((istep % nskip) == 0) {
     if(myid == 0){
      printf("******************************************************\n");
      printf("Step %d\n",istep);
      printf("Econs   %.8g\n",econs/((double) istep));
      printf("Temp    %.8g     %.8g\n",tinst,tbar/((double) (istep+1)));
      if(cmd_opt == 1){
       printf("Tcent   %.8g     %.8g\n",tinst,tbar/((double) (istep+1)));
      }/* endif */
      printf("Vir     %.8g     %.8g\n",vir,vir_bar/((double) (istep+1)));
      printf("Prim    %.8g     %.8g\n",prim,prim_bar/((double) (istep+1)));

      putchar('\n');
      fflush(stdout);
     }/* endif myid */


/* Write out configuration */

     write_coords(coord_file,coord_name,beads,therm,
                  myid,numprocs,mpi_vector,world); 

    }/* endif write to screen */

     if(cmd_opt == 1 && myid == 0 && ((istep % nwrite_cent) == 0)){
       fprintf(cent_file," %d     %.8g \n",istep,beads->u[1].x);
#ifdef SAVE_DISK_SPACE_OFF
       fprintf(cent_file," %d     %.8g \n",istep,beads->v[1].x);
#endif
       fflush(cent_file);
     }

  }/* endfor istep loop */

   Finalize();

 return 0;
}/* end main program */
/* -------------------------------------------------------------------------*/
double fict_kinet(BEADS *beads,MPI_Comm world,int numprocs)
{
#include "../include/typ_mask.h"
  double kinetic,my_kinet;
  int ibead;

  my_kinet = 0.0;
   for(ibead=1;ibead <= beads->my_np; ibead++){
     my_kinet += 0.5*beads->mp_stage[ibead]*
                 (beads->v[ibead].x*beads->v[ibead].x);
   }
  if(numprocs > 1) {
   Reduce(&my_kinet,&kinetic,1,MPI_DOUBLE,MPI_SUM,0,world);
  } else {
    kinetic = my_kinet;
  }/* endif */

  return kinetic;
}

/* -------------------------------------------------------------------------*/
double cent_kinet(BEADS *beads,MPI_Comm world,int numprocs)
{
#include "../include/typ_mask.h"
  double kinetic,my_kinet;
  int ibead;

  my_kinet = 0.0;
  my_kinet += 0.5*beads->mp_stage[1]*
                 (beads->v[1].x*beads->v[1].x);
  kinetic = my_kinet;

  return kinetic;
}

/* -------------------------------------------------------------------------*/
double therm_energy(THERM *therm,int my_np,double temp,MPI_Comm world,int numprocs)
{ 
#include "../include/typ_mask.h"
  double etherm,my_etherm;
  double beta;
  int ibead,ich,nch;
  int start_proc = therm->start_proc;

  VECTOR **eta    = therm->eta;
  VECTOR **etadot = therm->etadot;
  double *q       = therm->q;
 
  nch = therm->nch; 
  beta = RBOLTZ/temp;
  my_etherm = 0.0;
  for(ich=1;ich <= nch; ich++){
   for(ibead=start_proc;ibead <= my_np; ibead++){
    my_etherm += 0.5*q[ibead]*(
                 etadot[ich][ibead].x*etadot[ich][ibead].x)
               + eta[ich][ibead].x/beta;
   }
  }
  for(ibead=start_proc;ibead <= my_np; ibead++){
    my_etherm += 
      eta[nch-1][ibead].x/beta;
  }

  if(numprocs > 1){
   Reduce(&my_etherm,&etherm,1,MPI_DOUBLE,MPI_SUM,0,world);
  } else {
    etherm = my_etherm;
  }/* endif */

 return etherm;
} /* end function */

/* -------------------------------------------------------------------------*/
void write_coords(FILE *coord_file,char *coord_name,BEADS *beads,THERM *therm,
                  int myid,int numprocs,MPI_Datatype mpi_vector,
                  MPI_Comm world)
{
#include "../include/typ_mask.h"
   int ibead,ich,iproc;
   int istart,indb;
   static int master=0;
   double rrnp;
   VECTOR *vscratch1,*vscratch2;              /* Vector scratch arrays */
   MPI_Status stat;

   if(myid == 0){

/* malloc some vector scratch memory */
    vscratch1 = (VECTOR *)malloc(beads->np*sizeof(VECTOR))-1;

/* open coordinate file */
    coord_file = fopen(coord_name,"w");
   }/* endif myid */


/* Particle positions */

   if(numprocs > 1){
     Gather(&(beads->u[1]),beads->my_np,mpi_vector,
                &(vscratch1[1]),beads->my_np,mpi_vector,
                0,world);
   } else {
    for(ibead=1;ibead <= beads->np;ibead++){
      vscratch1[ibead].x = beads->u[ibead].x;
    }/* endfor */
   }/* endif numprocs */

    if(myid == 0){
     for(ibead=1; ibead <= beads->np; ibead++){
      fprintf(coord_file,"%lf\n",
             vscratch1[ibead].x);
     }/* endfor */
   } /* endif */

/* Particle velocities */

   if(numprocs > 1){
     Gather(&(beads->v[1]),beads->my_np,mpi_vector,
                &(vscratch1[1]),beads->my_np,mpi_vector,
                0,world);
   } else {
    for(ibead=1;ibead <= beads->np;ibead++){
      vscratch1[ibead].x = beads->v[ibead].x;
    }/* endfor */
   }/* endif */

    if(myid == 0){
     for(ibead=1; ibead <= beads->np; ibead++){
      fprintf(coord_file,"%lf\n",
             vscratch1[ibead].x);
    }/* endfor */
   } /* endif */

/* Thermostat positions */
    for(ich=1;ich <= therm->nch; ich++){
      if(numprocs > 1){
         Gather(&(therm->eta[ich][1]),beads->my_np,mpi_vector,
                    &(vscratch1[1]),beads->my_np,mpi_vector,
                    0,world);
      } else {
       for(ibead=1;ibead <= beads->np;ibead++){
         vscratch1[ibead].x = therm->eta[ich][ibead].x;
       }/* endfor */
      }/* endif */

    if(myid == 0){
     for(ibead=1; ibead <= beads->np; ibead++){
      fprintf(coord_file,"%lf\n",
             vscratch1[ibead].x);
     }/* endfor */
    }/* endif */
   }/* endfor ich */

/* Thermostat velocities */
    for(ich=1;ich <= therm->nch; ich++){
      if(numprocs > 1){
        Gather(&(therm->etadot[ich][1]),beads->my_np,mpi_vector,
                   &(vscratch1[1]),beads->my_np,mpi_vector,
                   0,world);
      } else {
        for(ibead=1;ibead <= beads->np;ibead++){
          vscratch1[ibead].x = therm->etadot[ich][ibead].x;
        }/* endfor */
      } /* endif */

    if(myid == 0){
     for(ibead=1; ibead <= beads->np; ibead++){
      fprintf(coord_file,"%lf\n",
             vscratch1[ibead].x);
     }/* endfor */
    }/* endif */
   }/* endfor ich */

   if(myid == 0){
    free(&(vscratch1[1]));
    fflush(coord_file);
    fclose(coord_file);
   }

}/* end function */

void assign_t_mats(T_MATRICES *t_matrices,int my_np,int np,int myid,
                   int trans_type,MPI_Comm world,int numprocs)
{
#include "../include/typ_mask.h"
   int i,j,k,np2,shift;
   double **a,*ap,**u,**ui,*eig,*fv1,*fv2;/* Malloc and use only if rs_ is called */
   double tol=1.0e-8;
   char c1='V',c2='L';
   double rrnp;
   double temp1,temp2;
   int job=1,ierr;                        /* Use only if rs_ is called */

 if(trans_type == 1){
 
   for(i=1;i<=my_np;i++) {
    for(j=1;j<=np;j++) {
      t_matrices->u_to_x[i][j] = 0.0;
      t_matrices->fx_to_fu[i][j] = 0.0;
    }
   }

   if(myid == 0){
    for(i=1;i<=np;i++){
      t_matrices->fx_to_fu[1][i] = 1.0;
    } 
    for(i=1;i<=my_np;i++){
      t_matrices->u_to_x[i][1] = 1.0;
    } 
    for(i=2;i<=my_np;i++) {
     for(j=i;j<=np;j++) {
       t_matrices->u_to_x[i][j] = ((double) (i-1)/(double) (j-1));
     }
    }
    for(i=2;i<=my_np;i++) {
     for(j=2;j<=i;j++) {
       t_matrices->fx_to_fu[i][j] = ((double) (j-1)/(double) (i-1));
     }
    }
   } else {
    for(i=1;i<=my_np;i++){
      t_matrices->u_to_x[i][1] = 1.0;
    } 
    for(i=1;i<=my_np;i++) {
     shift = myid*my_np + i;
     for(j=shift;j<=np;j++) {
      t_matrices->u_to_x[i][j] = ((double) (shift-1)/(double) (j-1));
     }
    }
    for(i=1;i<=my_np;i++) {
     shift = myid*my_np + i;
     for(j=2;j<=shift;j++) {
       t_matrices->fx_to_fu[i][j] = ((double) (j-1)/(double) (shift-1));
     }
    }
   } /* endif */

 } else {

   np2 = np*(np+1)/2;
   u  = mall_mat(1,np,1,np);
   ui = mall_mat(1,np,1,np);
   eig = (double *) malloc(np*sizeof(double))-1;
   a  = mall_mat(1,np,1,np);
   ap = (double *) malloc(np2*sizeof(double))-1;
   fv1 = (double *) malloc(3*np*sizeof(double))-1;
   fv2 = (double *) malloc(np*sizeof(double))-1;

  if(myid == 0){

/* Fill matrix */

   for(i=1;i<=np;i++){
    for(j=1;j<=np;j++){
     a[i][j] = 0.0;
    }
   }
   for(i=2;i<=np-1;i++){
    a[i][i-1] = -1.0;
    a[i][i] = 2.0;
    a[i][i+1] = -1.0;
   }
  a[1][1] = 2.0;
  a[1][2] = -1.0;
  a[np][np] = 2.0;
  a[np][np-1] = -1.0;
  a[1][np] = -1.0;
  a[np][1] = -1.0; 

    k=0;
    for(i=1;i<=np;i++){
      for(j=i;j<=np;j++){
       ++k;
       ap[k] = a[i][j];
      }
    }

   k=0;
   job=1;
   np2 = 2*np;
#define JAC
#ifdef JAC
   JACOBI(&np,&np,&(a[1][1]),&(eig[1]),&(u[1][1]),&job);
#endif
#ifdef DS
   DSPEV(&c1,&c2,&np,&(ap[1]),&(eig[1]),&(u[1][1]),&np,&(fv1[1]),&k);
   DSPEV(&job,&(ap[1]),&(eig[1]),&(u[1][1]),&np,&(fv1[1]),&np2);
#endif


   for(i=1;i<=np;i++){
    if(eig[i] < tol) k=i;
   }
   temp1 = eig[1];
   temp2 = eig[k];
   eig[1] = temp2;
   eig[k] = temp1;
   for(i=1;i<=np;i++){
    temp1 = u[1][i];
    temp2 = u[k][i];
    u[1][i] = temp2;
    u[k][i] = temp1;
   }


   rrnp = sqrt((double)np);
   for(i=1;i<=np;i++){
    for(j=1;j<=np;j++){
/*     u[i][j] *= -rrnp; */
       u[i][j] *= rrnp; 
    }
   }
   for(i=1;i<=np;i++){
    for(j=1;j<=np;j++){
     ui[i][j] = u[j][i];
    }
   }
   for(i=1;i<=np;i++){
     eig[i] *=(double)np;
   } 
  }/* endif myid */

/* Next divide up matrices among processors */

    np2 = np*np;
    if(numprocs > 1){
      Bcast(&(u[1][1]),np2,MPI_DOUBLE,0,world);
      Bcast(&(ui[1][1]),np2,MPI_DOUBLE,0,world);
      Bcast(&(eig[1]),np,MPI_DOUBLE,0,world);
    }/* endif */
   for(i=1;i<=my_np;i++){
     shift = myid*my_np + i;
     for(j=1;j<=np;j++){
      t_matrices->u_to_x[i][j] = ui[shift][j];
      t_matrices->fx_to_fu[i][j] = u[shift][j];
     }   
    t_matrices->eig[i] = eig[shift];
   }/* endfor */
   

/* Free allocated scratch memory */

     free((char *)(a[1]+np-NR_END));
     free((char *)(a+1-NR_END)); 
     free(&(fv1[1]));
     free(&(fv2[1]));
     free(&(ap[1]));
     free(&(eig[1]));
     free((char *)(u[1]+np-NR_END));
     free((char *)(u+1-NR_END));
     free((char *)(ui[1]+np-NR_END));
     free((char *)(ui+1-NR_END));

 } /* endif trans_type */

} /* end function */
