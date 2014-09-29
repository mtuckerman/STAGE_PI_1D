#include "../include/standard_include.h"

#include "../include/typedefs.h"
#include "../include/proto_observables.h"
#include "../include/proto_communicate_wrappers.h"
#include "../include/defines.h"


double vir_est(BEADS *beads,double temp,int myid,MPI_Comm world,int numprocs)
{
#include "../include/typ_mask.h"
 double my_est,est=0.0;
 double my_cent,cent=0.0;
 int ibead;

 my_cent = 0.0;
 for(ibead=1;ibead <= beads->my_np; ibead++){
   my_cent += beads->r[ibead].x;
 }  
 if(numprocs > 1){
   Reduce(&my_cent,&cent,1,MPI_DOUBLE,MPI_SUM,0,world);
 } else {
   cent = my_cent;
 }
 if(myid == 0) cent /= ((double) beads->np);

  my_est = 0.0;
  for(ibead=1;ibead <= beads->my_np; ibead++) {
   my_est -= 
   (  (beads->r[ibead].x-cent)*beads->fr[ibead].x);
  }

  if(numprocs > 1){
   Reduce(&my_est,&est,1,MPI_DOUBLE,MPI_SUM,0,world);
  } else {
    est = my_est;
  }/* endif */

 if(myid == 0){
  est /= 2.0;
  est += beads->v_ext + 0.5*temp;
 }

 return est;
} /* end function */
/*--------------------------------------------------------------------------*/
double prim_est(BEADS *beads,CHAIN *chain,double temp,int myid,
                MPI_Comm world,int numprocs)
{
#include "../include/typ_mask.h"
 double omp2;
 double beta;
 double my_est,est=0.0;
 int ibead;

  omp2 = chain->omp*chain->omp;
  beta = RBOLTZ/temp;
  my_est = 0.0;
  if(myid == 0) {
   for(ibead=2;ibead <= beads->my_np; ibead++) {
    my_est -= 0.5*beads->m_stage[ibead]*omp2*
    (  beads->u[ibead].x*beads->u[ibead].x);
   }
  } else {
   for(ibead=1;ibead <= beads->my_np; ibead++) {
    my_est -= 0.5*beads->m_stage[ibead]*omp2*
    (  beads->u[ibead].x*beads->u[ibead].x);
   }
 }

  if(numprocs > 1){
    Reduce(&my_est,&est,1,MPI_DOUBLE,MPI_SUM,0,world);
  } else {
    est = my_est;
  }/* endif */

 if(myid == 0){
  est += beads->v_ext + 0.5*((double) beads->np)/beta;
 }

 return est;
} /* end function */
