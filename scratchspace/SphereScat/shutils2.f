C  Copyright (C) 2007: Leslie Greengard
C  Contact: greengard@cims.nyu.edu
C
C  All rights reserved.
C  This copy is not for distribution without permission of the author.
C
C**********************************************************************C
C     Miscelllaneous utility functions for spherical harmonic
C     expansions.
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine unrolltoshxp(nterms,sol,sh1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine converts an "unrolled" spherical harmonic
c     expansion to a doubly indexed array where (i,j) refers to 
c     the Y_{i,j} coefficient.
c
c     INPUT:  
c             sol  (unrolled)
c     OUTPUT: 
c             sh1  (matrix format)
c
      implicit real *8 (a-h,o-z)
      complex *16 sol(1)
      complex *16 sh1(0:nterms,-nterms:nterms)
c
      next = 1
      do i = 0,nterms
         do j = -i,i
            sh1(i,j) = sol(next)
	    next = next+1
         enddo
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shxptounroll(nterms,sol,sh1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine converts a double indexed spherical harmonic
c     expansion to an "unrolled" one. 
c
c     INPUT: 
c              sh1  (matrix format)
c     OUTPUT:  
c              sol  (unrolled)
c
      implicit real *8 (a-h,o-z)
      complex *16 sol(1)
      complex *16 sh1(0:nterms,-nterms:nterms)
c
      next = 1
      do i = 0,nterms
         do j = -i,i
            sol(next) = sh1(i,j) 
	    next = next+1
         enddo
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sourcetoshxp(nterms,x0y0z0,source,
     1          str,zk0,sphererad,scale,pot,dpotdn,wrk,lw,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real *8 (a-h,o-z)
c
c------------------------------------------------------------------
c     WARNING!!!!  f77 static arrays set to length 1000. Thus,
c     expansion lengths limited. => nterms < 900 or so. 
c                                     
c------------------------------------------------------------------
c     utility function:
c
c     takes source at (source(1), source(2), source(3)) of strength str
c     and frequency zk0 
c     and computes induced potential and outward normal derivative 
c     of induced potential ON sphere centered at x0y0z0 of 
c     radius sphererad. These functions are tabulated as spherical
c     harmonic expansions.
c
c     zk0 is Helmholtz parameter.
c   
c     Method:
c
c     If source is inside sphere, H expansion is generated.
c     If source is outside sphere, J expansion is generated.
c
c
c     INPUT:
c
c     nterms        length of expansion
c     x0y0z0        sphere center
c     source        source location
c     str           source strength
c     zk0           Helmholtz parameter
c     sphererad     sphere radius
c     scale         expansion scaling parameter
c     wrk           workspace
c     lw            length of workspace 
c                   at least (nterms+1)**2 + 6*nterms + 10000
c
c     OUTPUT:
c
c     pot           spherical harmonic expansion of potential phi
c     dpotdn        spherical harmonic expansion of dphi/dn
c     ier           error rreturn code 
c
c
      integer *4 iw(0:1000)
      integer *4 lw
      dimension source(3)
      dimension x0y0z0(3)
      complex *16 pot(0:nterms,-nterms:nterms)
      complex *16 dpotdn(0:nterms,-nterms:nterms)
      complex *16 z,zk0,wronsk,str
      complex *16 htemp(0:1000),hptemp(0:1000)
      complex *16 wrk(1)
c
      pi = 4*datan(1.0d0)
cc      call prin2(' in sourcetoshxp source is *',source,3)
cc      call prin2(' in sourcetoshxp x0y0z0 is *',x0y0z0,3)
      dx = source(1)-x0y0z0(1)
      dy = source(2)-x0y0z0(2)
      dz = source(3)-x0y0z0(3)
      r = dsqrt(dx*dx+dy*dy+dz*dz)
ccc      call prin2(' in sourcetoshxp r is *',r,1)
ccc      call prin2(' in sourcetoshxp sphererad is *',sphererad,1)
      ns = 1
      ifder = 1
      if (r.lt.sphererad) then
	 call h3dformmp(ier,zk0,scale,source,str,ns,x0y0z0,nterms,
     1        pot,wrk,lw,lused)
         z = zk0*sphererad
         call h3dall(nterms,z,scale,htemp,ifder,hptemp)
ccc	 call prin2(' in sourcetoshxp htemp is *',htemp,2*nterms+2)
         do n = 0,nterms
            do m = -n,n
	       dpotdn(n,m) = pot(n,m)*zk0*hptemp(n)
	       pot(n,m) = pot(n,m)*htemp(n)
            enddo
         enddo
      else
	 call h3dformta(ier,zk0,scale,source,str,ns,x0y0z0,nterms,
     1        pot,wrk,lw,lused)
         z = zk0*sphererad
	 ljtemp = 1000
         call jfuns3d(ier,nterms,z,scale,htemp,ifder,hptemp,
     1        ljtemp,iw,ntop)
         do n = 0,nterms
            do m = -n,n
	       dpotdn(n,m) = pot(n,m)*zk0*hptemp(n)
	       pot(n,m) = pot(n,m)*htemp(n)
            enddo
         enddo
      endif
      return
      end
c
