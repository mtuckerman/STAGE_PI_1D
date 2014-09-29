#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXWORD 50

int main(int argc, char *argv[])
{
  FILE *ifp,*ofp,*wfp;
  int npts;
  int iii,i,j,ind,step;
  int ncorr;
  double *x,*y,*z,*r;
  double *vx,*vy,*vz;
  double *cqq,*cvv;
  double time,dt;
  double rbar = 0.0;
  char command[MAXWORD];

  sprintf(command,"wc -l %s > wc.out",argv[1]);
  system(command);
  wfp = fopen("wc.out","r");
  fscanf(wfp,"%d %*s",&npts);
  fclose(wfp);
  system("rm wc.out");

  if((ifp = fopen(argv[1],"r")) == NULL){
    printf("File %s not found\n",argv[1]);
    exit(1);
  }
  ofp = fopen(argv[2],"w"); 

  npts /=2;
  x = (double *) malloc(npts*sizeof(double))-1;
  y = (double *) malloc(npts*sizeof(double))-1;
  z = (double *) malloc(npts*sizeof(double))-1;
  r = (double *) malloc(npts*sizeof(double))-1;
  vx = (double *) malloc(npts*sizeof(double))-1;
  vy = (double *) malloc(npts*sizeof(double))-1;
  vz = (double *) malloc(npts*sizeof(double))-1;

  printf("Number of points is %d\n",npts);
  printf("How many points in correlation function?\n");
  scanf("%d",&ncorr);
  printf("what is the time step?\n");
  scanf("%lf",&dt);

  cqq = (double *) malloc(ncorr*sizeof(double));
  cvv = (double *) malloc(ncorr*sizeof(double));

/* Read in data */

  for(i=1;i<=npts;i++){
   fscanf(ifp,"%d %lf %lf %lf\n",&step,x+i,y+i,z+i);
#ifdef SAVE_DISK_SPACE_OFF
   fscanf(ifp,"%d %lf %lf %lf\n",&step,vx+i,vy+i,vz+i);
#endif
  }
  for(i=1;i<=npts;i++){
    r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]); 
  }
  for(i=1;i<=npts;i++){
    rbar += r[i];
  }
  rbar /=npts;
  printf("rbar %.12g \n",rbar);
  for(i=1;i<=npts;i++){
    r[i] -= rbar;
  }

/*Compute centroid correlation functions */

  for(i=0;i<ncorr;i++){
    cqq[i] = 0.0;   cvv[i] = 0.0;
    for(j=1;j<npts-i;j++){
      ind = i+j;
      cqq[i] += r[j]*r[ind];
    }/* endfor j */
   cqq[i] /= (npts-i);
  }/* endfor i */

  for(i=0;i<ncorr;i++){
   time = i*dt;
   fprintf(ofp,"%.8g %.8g \n",time,cqq[i]);
  }
  
}
