double fict_kinet(BEADS *,MPI_Comm,int);

double cent_kinet(BEADS *,MPI_Comm,int);

double therm_energy(THERM *,int,double,MPI_Comm,int);

void assign_t_mats(T_MATRICES *,int,int,int,int,MPI_Comm,int);

void write_coords(FILE *,char *,BEADS *,THERM *,int,int,MPI_Datatype,MPI_Comm);

