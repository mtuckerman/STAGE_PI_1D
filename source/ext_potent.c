#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#ifdef PARALLEL
#include "mpi.h"
#else
#include "../include/mpi_f.h"
#endif

#include "../include/typedefs.h"
#include "../include/proto_ext.h"
#include "../include/proto_communicate_wrappers.h"

void ext_force(BEADS *beads,MPI_Comm world,int numprocs)
{
#include "../include/typ_mask.h"
/* Local parameters of potential */
/* Harmonic oscillator potential */
/* Put your own potential in here if you want */

 static double omega=10.0;
 double om2;
 double my_vext;
 int ibead,iii;

 om2 = omega*omega;
 my_vext = 0.0;
 for(ibead=1;ibead <= beads->my_np; ibead++){
  beads->fr[ibead].x = -beads->mass*om2*beads->r[ibead].x;
  my_vext += 0.5*beads->mass*om2*
   (  beads->r[ibead].x*beads->r[ibead].x);
 }/* endfor */

 if(numprocs > 1) {
   Reduce(&my_vext,&beads->v_ext,1,MPI_DOUBLE,MPI_SUM,0,world);
 } else {
   beads->v_ext = my_vext;
 }/* endif */

} /* end function */

