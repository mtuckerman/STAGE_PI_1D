#include "../include/standard_include.h"

#include "../include/typedefs.h"
#include "../include/proto_init.h"
#include "../include/proto_utils.h"
#include "../include/proto_communicate_wrappers.h"
#include "../include/defines.h"

void *cmalloc(size_t);

void set_sim(SIM_PARAMS *sim_params,CHAIN *chain,BEADS *beads,THERM *therm,
             char *in_name, char *coord_name,int numprocs,int *nskip,
             int *nwrite_cent,int *cmd_opt,double *gamma)
           
{
#include "../include/typ_mask.h"
  FILE *in_file;

  if((in_file = fopen(in_name,"r")) == NULL){
    printf("File %s not found\n",in_name);
    exit(1);
  }
   fscanf(in_file,"%d",&(sim_params->nsteps));
   fscanf(in_file,"%d",&(beads->np));
   chain->np = beads->np;
   beads->my_np = beads->np/numprocs;
   chain->my_np = beads->my_np;
   fscanf(in_file,"%lf",&(sim_params->dt));
   fscanf(in_file,"%lf",&(sim_params->temp));
   fscanf(in_file,"%lf",&(sim_params->temp_p));
   fscanf(in_file,"%lf",&(beads->mass));
   fscanf(in_file,"%lf",&(beads->mass_const));
   fscanf(in_file,"%d",&(therm->nch));
   fscanf(in_file,"%d",&(therm->msuz));
   fscanf(in_file,"%d",&(therm->nres));
   fscanf(in_file,"%s",coord_name);
   fscanf(in_file,"%d %d",&(sim_params->trans_opt),&(sim_params->trans_type));
   if(sim_params->trans_type == 0 && sim_params->trans_opt != 1){
     printf("Cannot do a recursive transformation for normal modes\n");
     Finalize();
     exit(1);
   }
   fscanf(in_file,"%d %d",nskip,nwrite_cent);
   fscanf(in_file,"%lf",&(beads->x0));
   fscanf(in_file,"%d",cmd_opt);
   if(sim_params->trans_opt == 0 && *cmd_opt==1){
     printf("Cannot do a recursive transformation for centroids\n");
     Finalize();
     exit(1);
   }
   if(sim_params->trans_type != 0 && *cmd_opt==1){
     printf("Must do normal modes for centroids\n");
     Finalize();
     exit(1);
   }
   fscanf(in_file,"%lf",gamma);
   fclose(in_file);

  chain->omp = sim_params->temp*sqrt((double) chain->np)/RBOLTZ;

  switch(therm->msuz)
  {
    case 2:
     therm->ncalls = 3; break;
    case 3:
     therm->ncalls = 5; break;
    case 4:
     therm->ncalls = 7; break;
    case 5:
     therm->ncalls = 15; break;
    case 6:
     therm->ncalls = 25; break;
    case 7:
     therm->ncalls = 125; break;
    case 8:
     therm->ncalls = 625; break;
   default:
    printf("%d suzuki yoshida order not programmed\n", therm->msuz);
    exit(1);
  } /* end switch */

} /* end function */

/* -------------------------------------------------------------------------*/
void set_coords(BEADS *beads,THERM *therm,char *coord_name,double *qseed,
                double omp,double temp,double temp_p,int trans_type,
                MPI_Comm world,MPI_Datatype mpi_vector,
                int numprocs,int myid,double gamma)
{
#include "../include/typ_mask.h"
  FILE *coord_file;                 /* Input file pointer */
  double *scratch1;                 /* Scratch arrays */
  double *scratch2;
  double *scratch3;
  VECTOR *vscratch1;                /* Vector scratch */
  double beta;
  double omp2;
  double beta_p,omp_p,omp_p2;
  int ibead,ich;
  int iproc,indb;
  char command[MAXWORD];
  MPI_Status stat;
  int master=0;
  int np = beads->np;

  beta = RBOLTZ/temp;
  omp2 = omp*omp;
  beta_p = RBOLTZ/(temp_p);
  omp_p = sqrt(((double) np))*temp_p/RBOLTZ;
  omp_p2 = omp_p*omp_p;


 if(myid == master){
/* malloc some scratch */
  scratch1 = (double *)cmalloc(np*sizeof(double))-1;
  scratch2 = (double *)cmalloc(np*sizeof(double))-1;
  scratch3 = (double *)cmalloc(np*sizeof(double))-1;
  vscratch1 = (VECTOR *)cmalloc(np*sizeof(VECTOR))-1; 
 } /* endif myid */

 Barrier(world);


 if(trans_type != 1){
   if(numprocs > 1) {
    Gather(&(beads->eig[1]),beads->my_np,MPI_DOUBLE,
                &(scratch3[1]),beads->my_np,MPI_DOUBLE,
                0,world);
   } else {
     for(ibead=1;ibead <= beads->np;ibead++){
        scratch3[ibead] = beads->eig[ibead];
     }/* endfor */
   }/* endif */

 }


/* set staging masses */
 if(myid == master){
  scratch1[1] = beads->mass;
  scratch2[1] = beads->mass*beads->mass_const;
  if(trans_type == 1){
   for(ibead=2; ibead <= np; ibead++){
    scratch1[ibead] = ((double) ibead)/((double) (ibead-1))*beads->mass;
    scratch2[ibead] = scratch1[ibead]*beads->mass_const*gamma;
   } /* endfor */
  } else {
   for(ibead=2; ibead <= np; ibead++){
    scratch1[ibead] = scratch3[ibead]*beads->mass;
    scratch2[ibead] = scratch3[ibead]*beads->mass*beads->mass_const*gamma;
   } /* endfor */
  }/* endif trans_opt */
 }/* endif myid */

/* Divide masses among processors */

 if(numprocs > 1){
  Scatter(&(scratch1[1]),beads->my_np,MPI_DOUBLE,
              &(beads->m_stage[1]),beads->my_np,MPI_DOUBLE,
              0,world);

  Scatter(&(scratch2[1]),beads->my_np,MPI_DOUBLE,
              &(beads->mp_stage[1]),beads->my_np,MPI_DOUBLE,
              0,world);
 } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        beads->m_stage[ibead] = scratch1[ibead];
        beads->mp_stage[ibead] = scratch2[ibead];
     }/* endfor */

 }/* endif */
/* Set coordinates */

 if((coord_file = fopen(coord_name,"r")) == NULL) {

/* positions */
   if(myid == master){
    for(ibead=2;ibead <= np;ibead++){
     scratch2[ibead] = 1.0/sqrt(beta_p*omp_p2*scratch1[ibead]);
    } /* endfor */
    scratch2[1] = 1.0;
    gauss_v(vscratch1,scratch2,qseed,np);
    vscratch1[1].x = beads->x0;
   }/* endif myid */

/* Divide particle positions among processors */

   if(numprocs > 1) {
    Scatter(&(vscratch1[1]),beads->my_np,mpi_vector,
              &(beads->u[1]),beads->my_np,mpi_vector,
              0,world);
   } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        beads->u[ibead].x = vscratch1[ibead].x;
     }/* endfor */

   }/* endif */

/* velocities */

 if(myid == master) {
   for(ibead=1;ibead <= np;ibead++){
    scratch2[ibead] = 1.0/sqrt(beta*scratch1[ibead]);
   } /* endfor */
   gauss_v(vscratch1,scratch2,qseed,np);
  }/* endif myid */

/* Divide particle velocities among processors */

 if(numprocs > 1){
     Scatter(&(vscratch1[1]),beads->my_np,mpi_vector,
              &(beads->v[1]),beads->my_np,mpi_vector,
              0,world);
 } else {
     for(ibead=1;ibead <= beads->np;ibead++){
        beads->v[ibead].x = vscratch1[ibead].x;
     }/* endfor */

 }/* endif */

/* Thermostats */

  for(ich=1;ich <= therm->nch; ich++){
   for(ibead=1;ibead <= beads->my_np;ibead++){
    therm->eta[ich][ibead].x = 0.0;
   }
  }

/* Thermostat mass parameters */

 if(myid == master) {
  for(ibead=1;ibead <= np;ibead++){
   scratch1[ibead] = 1.0/(beta*omp2);
  }
 }/* endif myid */
/* Divide thermostat masses among processors */

 if(numprocs > 1) {
    Scatter(&(scratch1[1]),beads->my_np,MPI_DOUBLE,
              &(therm->q[1]),beads->my_np,MPI_DOUBLE,
              0,world);
 } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        therm->q[ibead] = scratch1[ibead];
     }/* endfor */

 }/* endif */

/* Thermostat velocities */

 if(myid == master) {
   for(ibead=1;ibead <= np;ibead++){
    scratch2[ibead] = 1.0/sqrt(beta*scratch1[ibead]);
   } /* endfor */
 }/* endif myid */
  for(ich=1;ich <= therm->nch; ich++){

   if(myid == master) gauss_v(vscratch1,scratch2,qseed,np);

/* Divide thermostat velocities among processors */

   if(numprocs > 1){
     Scatter(&(vscratch1[1]),beads->my_np,mpi_vector,
              &(therm->etadot[ich][1]),beads->my_np,mpi_vector,
              0,world);
   } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        therm->etadot[ich][ibead].x = vscratch1[ibead].x;
     }/* endfor */
   }/* endif */
  } /* endfor ich */

/* Touch coord file to make sure it is there */

  if(myid == master){
   sprintf(command,"touch %s",coord_name);
   system(command);
  }

 } else /* if coord file is found */ {
/* Set thermostat mass parameters */

  if(myid == master) {
   for(ibead=1;ibead <= np;ibead++){
    scratch1[ibead] = 1.0/(beta*omp2);
   }
  }/* endif */
/* Divide thermostat masses among processors */

  if(numprocs > 1){
     Scatter(&(scratch1[1]),beads->my_np,MPI_DOUBLE,
              &(therm->q[1]),beads->my_np,MPI_DOUBLE,
              0,world);
  } else {
     for(ibead=1;ibead <= beads->np;ibead++){
        therm->q[ibead] = scratch1[ibead];
     }/* endfor */
  } /* endif */

/* Particle positions */
 if(myid == master){
  for(ibead=1;ibead <= np;ibead++){
   if(
      (fscanf(coord_file,"%lf",
         &(vscratch1[ibead].x)))
     < 0)  {printf("File empty:  %s\n",coord_name); exit(1);}
  }/* endfor */
 }/* endif */
/* Divide particle positions among processors */

 if(numprocs > 1){
     Scatter(&(vscratch1[1]),beads->my_np,mpi_vector,
              &(beads->u[1]),beads->my_np,mpi_vector,
              0,world);

 } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        beads->u[ibead].x = vscratch1[ibead].x;
     }/* endfor */
 } /* endif */

/* Particle velocities */
 if(myid == master) {
  for(ibead=1;ibead <= np;ibead++){
   fscanf(coord_file,"%lf",
         &(vscratch1[ibead].x));
  }/* endfor */
 }/* endif */
/* Divide particle velocities among processors */

 if(numprocs > 1){
     Scatter(&(vscratch1[1]),beads->my_np,mpi_vector,
              &(beads->v[1]),beads->my_np,mpi_vector,
              0,world);
 } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        beads->v[ibead].x = vscratch1[ibead].x;

     }/* endfor */
 } /* endif */
/* Thermostat positions */
 for(ich=1;ich <= therm->nch; ich++){
  if(myid == master){ 
   for(ibead=1;ibead <= np;ibead++){
    fscanf(coord_file,"%lf",
          &(vscratch1[ibead].x));
   }
  } /* endif myid */
/* Divide thermostat positions among processors */

  if(numprocs > 1){
    Scatter(&(vscratch1[1]),beads->my_np,mpi_vector,
               &(therm->eta[ich][1]),beads->my_np,mpi_vector,
               0,world);

  } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        therm->eta[ich][ibead].x = vscratch1[ibead].x;
     }/* endfor */
  } /* endif */
 } /* endfor ich */

/* Thermostat velocities */
 for(ich=1;ich <= therm->nch; ich++){
  if(myid == master){
   for(ibead=1;ibead <= np;ibead++){
    fscanf(coord_file,"%lf",
          &(vscratch1[ibead].x));
   }
  }/* endif myid */
/* Divide thermostat velocities among processors */
  
  if(numprocs > 1){
    Scatter(&(vscratch1[1]),beads->my_np,mpi_vector,
               &(therm->etadot[ich][1]),beads->my_np,mpi_vector,
               0,world);
  } else {

     for(ibead=1;ibead <= beads->np;ibead++){
        therm->etadot[ich][ibead].x = vscratch1[ibead].x;
     }/* endfor */
  } /* endif */
 } /* endfor ich */

  if(myid == master) fclose(coord_file);
}/* endif */

 if(myid == master) {
  free(&(scratch1[1]));
  free(&(scratch2[1]));
  free(&(scratch3[1]));
  free(&(vscratch1[1]));
 }

} /* end function */

/* -------------------------------------------------------------------------*/
void init_suz(THERM *therm,double dt)
{
#include "../include/typ_mask.h"
  double psuz[6][6];
  double xnt,dno,dnit;
  double pyosh;
  double w_yosh6_1,w_yosh6_2,w_yosh6_3,w_yosh6_0;
  double w_yosh8_1,w_yosh8_2,w_yosh8_3,w_yosh8_4,w_yosh8_5,w_yosh8_6,w_yosh8_7,w_yosh8_0;
  int no,i,j,k,l,m,iii;

  dnit = (double) therm->nres;
  for(i=0; i<= 4; i++)
   for(j=0; j<= 5; j++)
    psuz[i][j] = 0.0;

  for(no=2; no<=therm->msuz; no++) {
    dno = (double) no;
    xnt = 1.0/(2.0*dno - 1.0);
    for(i=1; i <=5; i++) {
      psuz[no][i] = 1.0/(4.0 - pow(4.0,xnt));
      if((i % 3) == 0)
        psuz[no][i] = 1.0 - 4.0/(4.0 - pow(4.0,xnt));
     } /* endfor */
  } /* endfor */

 switch(therm->msuz) {
   case 2:
     pyosh = 1.0/(2.0 - pow(2.0,(1.0/3.0)));
     therm->dtsuz[1] = pyosh*dt/dnit;
     therm->dtsuz[2] = -(pow(2.0,(1.0/3.0))*pyosh*dt/dnit);
     therm->dtsuz[3] = pyosh*dt/dnit;
     break;

   case 3:
    k=0;
    for(i=1; i<=5; i++) {
      k++;
      therm->dtsuz[k] = psuz[2][i]*dt/dnit;  
    }
    break;

    case 4:
       w_yosh6_1 =  0.784513610477560;
       w_yosh6_2 =  0.235573213359357;
       w_yosh6_3 = -1.177679984178870;
       w_yosh6_0 = 1.00 - 2.00*(w_yosh6_1+w_yosh6_2+w_yosh6_3);
       therm->dtsuz[1] = w_yosh6_3*dt/dnit;
       therm->dtsuz[2] = w_yosh6_2*dt/dnit;
       therm->dtsuz[3] = w_yosh6_1*dt/dnit;
       therm->dtsuz[4] = w_yosh6_0*dt/dnit;
       therm->dtsuz[5] = w_yosh6_1*dt/dnit;
       therm->dtsuz[6] = w_yosh6_2*dt/dnit;
       therm->dtsuz[7] = w_yosh6_3*dt/dnit;
       break;

    case 5:
      w_yosh8_1 =  0.102799849391985;
      w_yosh8_2 = -1.96061023297549;
      w_yosh8_3 =  1.93813913762276;
      w_yosh8_4 = -0.158240635368243;
      w_yosh8_5 = -1.44485223686048;
      w_yosh8_6 =  0.253693336566229;
      w_yosh8_7 =  0.914844246229740;
      w_yosh8_0 = 1.0 - 2.0*(w_yosh8_1+w_yosh8_2+w_yosh8_3
                           +w_yosh8_4+w_yosh8_5+w_yosh8_6
                           +w_yosh8_7);
       therm->dtsuz[1] = w_yosh8_7*dt/dnit;
       therm->dtsuz[2] = w_yosh8_6*dt/dnit;
       therm->dtsuz[3] = w_yosh8_5*dt/dnit;
       therm->dtsuz[4] = w_yosh8_4*dt/dnit;
       therm->dtsuz[5] = w_yosh8_3*dt/dnit;
       therm->dtsuz[6] = w_yosh8_2*dt/dnit;
       therm->dtsuz[7] = w_yosh8_1*dt/dnit;
       therm->dtsuz[8] = w_yosh8_0*dt/dnit;
       therm->dtsuz[9] = w_yosh8_1*dt/dnit;
       therm->dtsuz[10] = w_yosh8_2*dt/dnit;
       therm->dtsuz[11] = w_yosh8_3*dt/dnit;
       therm->dtsuz[12] = w_yosh8_4*dt/dnit;
       therm->dtsuz[13] = w_yosh8_5*dt/dnit;
       therm->dtsuz[14] = w_yosh8_6*dt/dnit;
       therm->dtsuz[15] = w_yosh8_7*dt/dnit;
       break;

    case 6:
    k=0;
    for(j=1; j<=5; j++) {
     for(i=1; i<=5; i++) {
      k++;
      therm->dtsuz[k] = psuz[2][i]*psuz[3][j]*dt/dnit;
     }
    }
   break;
    case 7:
    k=0;
    for(j=1; j<=5; j++) {
     for(i=1; i<=5; i++) {
      for(l=1; l<=5; l++) {
       k++;
       therm->dtsuz[k] = psuz[2][i]*psuz[3][j]*psuz[4][l]*dt/dnit;
      }
     }
    }
    break;
    case 8:
     k=0;
     for(j=1; j<=5; j++) {
      for(i=1; i<=5; i++) {
       for(l=1; l<=5; l++) {
        for(m=1; m<=5; m++) {
         k++;
         therm->dtsuz[k] = psuz[2][i]*psuz[3][j]*psuz[4][l]* 
                           psuz[5][m]*dt/dnit;
       }
      }
     }
    }
   break;
  } /* end switch */
}

void *cmalloc(size_t len)
{ /* begin routine */
  void *mem_ptr;
  double request;

  if(len == 0) return NULL;
  mem_ptr = malloc(len);
  if(mem_ptr == NULL) {
   request = ((double) len)*1.0e-6;
   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   printf("Dude, like you've just requested %g MBytes\n",request);
   printf("of memory -- get real\n\n");
   printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
   fflush(stdout);
   exit(1);
  }

  return mem_ptr;

} /* end routine */



