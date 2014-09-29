#include "../include/standard_include.h"

#include "../include/typedefs.h"
#include "../include/proto_ext.h"
#include "../include/proto_communicate_wrappers.h"

void trans_u_to_x_rec(BEADS *beads)
{
#include "../include/typ_mask.h"
  int ibead,ibeadp,ibeadm,np;
  int iii;

  np = beads->np;
  beads->r[1].x = beads->u[1].x;

  beads->r[np].x = beads->u[np].x + beads->u[1].x;

  for(ibead=beads->np-1;ibead >=2;ibead--){
   ibeadp = ibead+1;
   ibeadm = ibead-1;
   beads->r[ibead].x = beads->u[ibead].x + 
          (((double) ibeadm)*beads->r[ibeadp].x + beads->r[1].x)
         /((double) ibead);
  }/* endfor */
} /* end function */

void trans_u_to_x_rec_par(BEADS *beads,int myid,int numprocs,
                          MPI_Datatype mpi_vector,MPI_Comm world)
{
#include "../include/typ_mask.h"
  int ibead,ibeadp,ibeadm,np;
  int ibeadm_r,ibead_r;
  int my_np,iproc;
  VECTOR utemp,xtemp;
  MPI_Status stat;

  np = beads->np;
  my_np = beads->my_np;

/* store the first bead in a temporary */
  if(myid == 0) {
    utemp.x = beads->u[1].x;
  }

/* Broadcast to other processors */
  if(numprocs > 1) Bcast(&utemp,1,mpi_vector,0,world);

/* Handle endpoint beads 1st and last ones */
  if(myid == 0) {
    beads->r[1].x = utemp.x;
  }
  if(myid == numprocs-1){ 
    beads->r[my_np].x = beads->u[my_np].x + utemp.x;
  }

/* Do staging transformation */
/* Code for first and last processors looks different so do them
   separately                                                     */

 if(myid == numprocs-1){
  iproc = numprocs-1;
  for(ibead=beads->my_np-1;ibead >=1;ibead--){
   ibeadp = ibead+1;
   ibeadm = ibead-1;
   ibeadm_r = iproc*beads->my_np + ibeadm;
   ibead_r  = iproc*beads->my_np + ibead;
   beads->r[ibead].x = beads->u[ibead].x + 
          (((double) ibeadm_r)*beads->r[ibeadp].x + utemp.x)
         /((double) ibead_r);
  }/* endfor */
  xtemp.x = beads->r[1].x;

  Ssend(&xtemp,1,mpi_vector,iproc-1,iproc-1,world);

 }/* endif myid */
 if(myid == 0){

  Recv(&xtemp,1,mpi_vector,1,MPI_ANY_TAG,world);


   if(beads->my_np > 1){
    ibead = beads->my_np;
    ibeadm = ibead-1;
    beads->r[ibead].x = beads->u[ibead].x + 
           (((double) ibeadm)*xtemp.x + utemp.x)
          /((double) ibead);
  }
  for(ibead=beads->my_np-1;ibead >=2;ibead--){
   ibeadp = ibead+1;
   ibeadm = ibead-1;
   beads->r[ibead].x = beads->u[ibead].x + 
          (((double) ibeadm)*beads->r[ibeadp].x + utemp.x)
         /((double) ibead);
  }/* endfor */
 }

/* All other processors */
 for(iproc=numprocs-2;iproc >=1;iproc--){
  if(myid == iproc){

   Recv(&xtemp,1,mpi_vector,iproc+1,MPI_ANY_TAG,world);

    ibead  = beads->my_np;
    ibeadm = ibead-1;
    ibeadm_r = iproc*beads->my_np + ibeadm;
    ibead_r  = iproc*beads->my_np + ibead;
    beads->r[ibead].x = beads->u[ibead].x + 
           (((double) ibeadm_r)*xtemp.x + utemp.x)
          /((double) ibead_r);
 
   for(ibead=beads->my_np-1;ibead >=1;ibead--){
    ibeadp = ibead+1;
    ibeadm = ibead-1;
    ibeadm_r = iproc*beads->my_np + ibeadm;
    ibead_r  = iproc*beads->my_np + ibead;
    beads->r[ibead].x = beads->u[ibead].x + 
           (((double) ibeadm_r)*beads->r[ibeadp].x + utemp.x)
          /((double) ibead_r);
   }/* endfor */
   xtemp.x = beads->r[1].x;

   Ssend(&xtemp,1,mpi_vector,iproc-1,iproc-1,world);

  }
 }

} /* end function */

void trans_fx_to_fu_rec(BEADS *beads)
{
#include "../include/typ_mask.h"
  int ibead,ibeadp,ibeadm,ibeadm2,np;
  double sum_x=0.0;

  np = beads->np;
  for(ibead=1;ibead <= np;ibead++){
   beads->fr[ibead].x /= ((double) np);
  }
 beads->v_ext /=((double) np);

  for(ibead=1;ibead <= np;ibead++){
   sum_x += beads->fr[ibead].x;
  }/* endfor */
  beads->fu_ext[1].x = sum_x;

  for(ibead=2;ibead <= np;ibead++){
   ibeadm2 = ibead-2;
   ibeadm  = ibead-1;
   beads->fu_ext[ibead].x = 
          ((double) ibeadm2)*beads->fu_ext[ibeadm].x/((double) ibeadm)
       +  beads->fr[ibead].x;
 } /* endfor */ 
} /* end function */


void trans_fx_to_fu_rec_par(BEADS *beads,int myid,int numprocs,
                            MPI_Datatype mpi_vector,MPI_Comm world)
{
#include "../include/typ_mask.h"
  int ibead,ibeadp,ibeadm,ibeadm2,np;
  int ibeadm_r,ibeadm2_r;
  int iproc;
  double sum_x=0.0;
  MPI_Status stat;
  VECTOR fxtemp;

  np = beads->np;
  for(ibead=1;ibead <= beads->my_np;ibead++){
   beads->fr[ibead].x /= ((double) np);
  }
  if(myid == 0) beads->v_ext /=((double) np);

  for(ibead=1;ibead <= beads->my_np;ibead++){
   sum_x += beads->fr[ibead].x;
  }/* endfor */

  Reduce(&sum_x,&(beads->fu_ext[1].x),1,MPI_DOUBLE,MPI_SUM,0,world);

  if(myid == numprocs-1){

   Recv(&fxtemp,1,mpi_vector,numprocs-2,MPI_ANY_TAG,world);

    iproc = numprocs-1;
    ibead = 1;
    ibeadm2_r = iproc*beads->my_np + ibead - 2;
    ibeadm  = ibead-1;
    ibeadm_r = iproc*beads->my_np + ibead - 1;
    beads->fu_ext[ibead].x = 
           ((double) ibeadm2_r)*fxtemp.x/((double) ibeadm_r)
        +  beads->fr[ibead].x;
   for(ibead=2;ibead <= beads->my_np;ibead++){
    ibeadm2 = ibead-2;
    ibeadm2_r = iproc*beads->my_np + ibead - 2;
    ibeadm  = ibead-1;
    ibeadm_r = iproc*beads->my_np + ibead - 1;
    beads->fu_ext[ibead].x = 
           ((double) ibeadm2_r)*beads->fu_ext[ibeadm].x/((double) ibeadm_r)
        +  beads->fr[ibead].x;
   } /* endfor */ 
   fxtemp = beads->fu_ext[beads->my_np];
  }
  if(myid == 0) {
   for(ibead=2;ibead <= beads->my_np;ibead++){
    ibeadm2 = ibead-2;
    ibeadm  = ibead-1;
    beads->fu_ext[ibead].x = 
           ((double) ibeadm2)*beads->fu_ext[ibeadm].x/((double) ibeadm)
        +  beads->fr[ibead].x;
   } /* endfor */ 
   fxtemp.x = beads->fu_ext[beads->my_np].x;

   Ssend(&fxtemp,1,mpi_vector,1,1,world);

  }
  for(iproc=1;iproc <= numprocs-2; iproc++){
    if(myid == iproc){

     Recv(&fxtemp,1,mpi_vector,iproc-1,MPI_ANY_TAG,world);

      ibead = 1;
      ibeadm2_r = iproc*beads->my_np + ibead - 2;
      ibeadm_r = iproc*beads->my_np + ibead -1;
      beads->fu_ext[ibead].x = 
             ((double) ibeadm2_r)*fxtemp.x/((double) ibeadm_r)
          +  beads->fr[ibead].x;
     for(ibead=2;ibead <= beads->my_np;ibead++){
      ibeadm2 = ibead-2;
      ibeadm2_r = iproc*beads->my_np + ibead - 2;
      ibeadm  = ibead-1;
      ibeadm_r = iproc*beads->my_np + ibead - 1;
      beads->fu_ext[ibead].x = 
             ((double) ibeadm2_r)*beads->fu_ext[ibeadm].x/((double) ibeadm_r)
          +  beads->fr[ibead].x;
     } /* endfor */ 
     fxtemp.x = beads->fu_ext[beads->my_np].x;

     Ssend(&fxtemp,1,mpi_vector,iproc+1,iproc+1,world);

    }
  }


} /* end function */

void trans_u_to_x_mat(BEADS *beads,T_MATRICES *t_matrices,
                  int myid,int numprocs,
                  MPI_Datatype mpi_vector,MPI_Comm world)
{
#include "../include/typ_mask.h"
  VECTOR *utemp;
  int i,j;

/* allocate some scratch space */
  utemp = (VECTOR *) malloc(beads->np*sizeof(VECTOR))-1;

/* Gather from all the processors the u values into one vector */

  if(numprocs > 1){
    Allgather(&(beads->u[1]),beads->my_np,mpi_vector,
                  &(utemp[1]),beads->my_np,mpi_vector,0,world);
  } else {
    for(i=1;i <= beads->np;i++){
      utemp[i].x = beads->u[i].x;
    }/* endfor */
  }/* endfor */

/* Each processor transforms to get its r values */
  for(i=1;i<=beads->my_np;i++){
    beads->r[i].x = 0.0;
  }
  for(i=1;i<=beads->my_np;i++){
   for(j=1;j<=beads->np;j++){
     beads->r[i].x += t_matrices->u_to_x[i][j]*utemp[j].x;
   }
  }

  free(&(utemp[1]));
}
void trans_fx_to_fu_mat(BEADS *beads,T_MATRICES *t_matrices,
                    int myid,int numprocs,int trans_type,
                    MPI_Datatype mpi_vector,MPI_Comm world)
{
#include "../include/typ_mask.h"
   VECTOR *frtemp;
   int i,j;
   int np,ibead;

/* allocate some scratch space */
   frtemp = (VECTOR *) malloc(beads->np*sizeof(VECTOR))-1;

  np = beads->np;
   for(ibead=1;ibead <= beads->my_np;ibead++){
    beads->fr[ibead].x /= ((double) np);
   }
  if(myid == 0) beads->v_ext /=((double) np);

/* Gather from all the processors the fr values into one vector */
  if(numprocs > 1){
    Allgather(&(beads->fr[1]),beads->my_np,mpi_vector,
                  &(frtemp[1]),beads->my_np,mpi_vector,0,world);

  } else { 

    for(ibead=1;ibead <= beads->np;ibead++){
      frtemp[ibead].x = beads->fr[ibead].x;
    }/* endfor */

  }/* endif */
/* Each processor transforms to get its fu values */
  for(i=1;i<=beads->my_np;i++){
    beads->fu_ext[i].x = 0.0;
  }
  for(i=1;i<=beads->my_np;i++){
   for(j=1;j<=beads->np;j++){
     beads->fu_ext[i].x += t_matrices->fx_to_fu[i][j]*frtemp[j].x;
   }
  }

   free(&(frtemp[1]));
}
