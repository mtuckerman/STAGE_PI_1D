           implicit real*8 (a-h,o-z)
           open(unit=1,file='cent_file_tmp',status='old')
           open(unit=2,file='trace.out')
C
           WRITE(6,*) 'Enter time step'
           READ(5,*) dt
 10        CONTINUE
           READ(1,*,END=20) istep,x1,y1,z1
           dx = x1
           dy = y1
           dz = z1
           r = sqrt(dx*dx + dy*dy + dz*dz)
           write(2,*) dfloat(istep)*dt*0.024188/1000.0,r
           GOTO 10
 20        CONTINUE
           stop
           end
           
