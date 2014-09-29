void ext_force(BEADS *,MPI_Comm,int);
void trans_u_to_x_rec_par(BEADS *,int ,int ,MPI_Datatype ,MPI_Comm );
void trans_u_to_x_rec(BEADS *);
void trans_fx_to_fu_rec_par(BEADS *,int ,int ,MPI_Datatype,MPI_Comm);
void trans_fx_to_fu_rec(BEADS *);
void trans_u_to_x_mat(BEADS *,T_MATRICES *,int ,int ,
                      MPI_Datatype ,MPI_Comm );
void trans_fx_to_fu_mat(BEADS *,T_MATRICES *,int ,int ,int ,
                        MPI_Datatype,MPI_Comm);
