         implicit real*8 (a-h,o-z)
         open(unit=1,file='Prim_128')
         open(unit=2,file='Prim_128s')
         npts = 20000
         do i=1,npts
           read(1,*) x,y
           write(2,*) 50*i,x,y
         enddo
         stop
         end
