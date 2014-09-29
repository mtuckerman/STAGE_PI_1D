#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/typedefs.h"
#include "../include/proto_utils.h"
#include "../include/defines.h"

void gauss_v(VECTOR *v,double *width,double *qseed,int num)
{
  int i,iii;
  double phi[3],chi[3],twopi;

  twopi = 2.0*M_PI;
  for(i=1; i<= num;i++){
    phi[0] = twopi*ran_essl(qseed);
    chi[0] = ran_essl(qseed);
    v[i].x = sqrt(log(chi[0])*(-2.0))*cos(phi[0])*width[i];
   }
} /* end function */

/* -------------------------------------------------------------------------*/
double ran_essl(double *qseed)
{
  int n=1;
  double x;
  DURAND(qseed,&n,&x);
  return x;
}

/* -------------------------------------------------------------------------*/

VECTOR **cmall_mat(long nrl, long nrh, long ncl, long nch)
/* allocate a VECTOR matrix with subscript range m[nrl...nrh][ncl...nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  VECTOR **m;

  if(nrow <= 0 || ncol <= 0) return (VECTOR **)NULL;
/* allocate pointers to rows */
  m=(VECTOR **) malloc((size_t)((nrow+NR_END)*sizeof(VECTOR*)));
  if(!m) 
    {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("allocation failure of row pointers in cmall_mat\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(VECTOR *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(VECTOR)));
  if(!m[nrl]) 
    {
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("allocation failure of rows in cmall_mat\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;

/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */

/* -------------------------------------------------------------------------*/

double **mall_mat(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl...nrh][ncl...nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  if(nrow <= 0 || ncol <= 0) return (double **)NULL;
/* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if(!m) 
    {
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("allocation failure of row pointers in mall_mat\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }
  m += NR_END;
  m -= nrl;

/* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if(!m[nrl]) 
    {
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("allocation failure of rows in mall_mat\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;

/* return pointer to array of pointers to rows */
   return m;  
}/* end routine */






