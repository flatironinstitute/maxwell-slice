      implicit real *8 (a-h,o-z)
      x = 1.0d0
      y = 5.0d0
      z = 7.0d0
      u = 4.0d0
      open(unit=20, file='data.m',status='unknown')
      open(unit=22, file='dd.bin',status='new',form='unformatted',
     1      access='direct', recl=8*3)
      write(20,*) x,y,z,u
      write(22,rec=1) x, y, z
      write(22,rec=2) x, y, z
c   write(22,rec=2) y
c   write(22,rec=3) z
c   write(22,rec=4) u
      rewind(20)
      read(20,*) x,y,z,u
      write(6,*) x,y,z,u
c   rewind(22)
c   read(22) x,y,z,u
      write(6,*) x,y,z,u
      stop
      end