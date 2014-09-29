void chain_force(CHAIN *, BEADS *,int,MPI_Comm,int);

void integrate(CHAIN *,BEADS* ,THERM *,T_MATRICES *,
               double,double,int,int,int,int,
               MPI_Datatype,MPI_Comm);

void suz_int(THERM *,BEADS *,double ,double );
