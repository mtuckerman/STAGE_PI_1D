void set_sim(SIM_PARAMS *,CHAIN *,BEADS *,THERM *,char *,char *,int,int *,
             int *,int *,double *);

void set_coords(BEADS *,THERM *,char *,double *,double,double,double,
                int,MPI_Comm,MPI_Datatype,int,int,double);

void init_suz(THERM *,double );
