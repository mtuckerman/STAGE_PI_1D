 typedef struct {
   double x;
 } VECTOR;

 typedef struct {
  VECTOR *r,*fr;                    /* Primitive Variables */
  VECTOR *u,*fu,*fu_ext;            /* Staging Variables */
  VECTOR *v;                        /* Fictitious Velocities */
  double *m_stage;                  /* Staging masses */
  double *mp_stage;                 /* Kinetic masses = const x m_stage */
  double *eig;
  int np,my_np;                     /* Number of path integral beads 
                                       and number on a processor */
  double mass;                      /* Mass of the quantum particle */
  double mass_const;                /* Const fact between real mass and  */
  double v_ext;                     /* External potential */
  double x0;                        /* Initial condition on centroid */
 } BEADS;

typedef struct {
  double omp;                       /* Chain Frequency */
  double v_chain;                   /* Harmonic potential of chain */
  int np,my_np;                     /* Number of path integral beads 
                                       and number on a processor */
} CHAIN;

typedef struct{
  double dt;                        /* Time step */
  int nsteps;                       /* Number of time steps */
  double temp;                      /* Simulation temperature */
  double temp_p;                    /* Initialization Temperature for paths */
  int trans_opt;                    /* Use recursive (0) or matrix (1) transformation */
  int trans_type;                   /* Use staging or normal mode transformation*/
} SIM_PARAMS;

typedef struct {
 VECTOR **eta;                     /* Thermostat position */
 VECTOR **etadot;                  /* Thermostat velocity */
 VECTOR **feta;                    /* Thermostat force */
 double *q;                        /* Thermostat mass parameter */
 double *dtsuz;                    /* Suzuki/Yoshida time steps */
 int msuz,ncalls,nch,nres;         /* Thermostat integration parameters */
 int start_proc;                   /* Bead starting index for centroids */
} THERM;

typedef struct{
  double **u_to_x;
  double **fx_to_fu;
  double *eig;
} T_MATRICES;
