c
c
c***********************************************************************
      subroutine rotviaproj(alpha,nterms,m1,m2,mpole,lmp,marray2,lmpn)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for applying rotation operator about
c	the y-axis determined by angle alpha.
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c       If the rotation operator is to be used multiple times, and
c       memory is available, one can precompute and store the 
c       multipliers used in evalall (see below). This has not yet been
c       implemented.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       alpha:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
c
c       m1    : max m index for first expansion   NOT IMPLEMENTED
c       m2    : max m index for second expansion   NOT IMPLEMENTED
c
c               PARAMETERS CURRENTLY IGNORED.
c
C       mpole   coefficients of original multiple expansion
C       lmp     leading dim for mpole (must exceed nterms)
C       lmpn    leading dim for marray2 (must exceed nterms)
c       w     :  work array 
c       lw    :  length of work array 
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c       lused    amount of workspace used.
c       ier      error return flag
c                0 successful execution
c                1 insufficient memory
c
C---------------------------------------------------------------------
c
c
c
      implicit real *8 (a-h,o-z)
      integer nquad
      real *8 w
      allocatable w(:)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
c
      nquad = 2*nterms+2
      ictheta = 1
      istheta = ictheta+nquad
      icphi = istheta+nquad
      isphi = icphi+nquad
      iynm = isphi+nquad
      iynmd = iynm + (nterms+1)**2
      irat1 = iynmd + (nterms+1)**2
      irat2 = irat1 + (nterms+1)**2
      iuval = irat2 + (nterms+1)**2
      iuder = iuval + 2*nquad*(nterms+1)
      iephi = iuder + 2*nquad*(nterms+1)
      iwsave = iephi + 2*(2*nterms+1)
      iavec = iwsave + 4*nquad+20
      ibvec = iavec + 2*nquad
      lused = ibvec + 2*nquad
      allocate (w(lused))
c
      call rotviaproj0(alpha,nquad,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),w(iuval),w(iuder),
     1           w(iephi),w(iwsave),w(iavec),w(ibvec))
      deallocate(w)
      return
      end
c
c
c***********************************************************************
      subroutine rotviaproj0(alpha,nquad,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn,cthetas,sthetas,cphis,sphis,ynm,ynmd,
     1           rat1,rat2,uval,uder,ephis,wsave,avec,bvec)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad,nterms
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
      complex *16 avec(nquad)
      complex *16 bvec(nquad)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 ephis(-nterms:nterms)
      real *8 wsave(4*nquad+20)
C
C
      call getmeridian(alpha,nquad,cthetas,sthetas,cphis,sphis)    
      call evalall(alpha,nquad,cthetas,sthetas,cphis,sphis,mpole,
     2           lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
      call projectonynm(nquad,uval,uder,ynm,ynmd,marray2,lmpn,nterms,
     2           m2,wsave,avec,bvec)
      return
      end
C
C
C***********************************************************************
      subroutine getmeridian(alpha,nquad,cthetas,sthetas,cphis,sphis)
C***********************************************************************
C     Purpose:
C
C           For a rotation of angle ALPHA about the y-axis, this
C           subroutine returns the NQUAD equispaced nodes on the 
C           rotated equator in the original coordinate system.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     alpha  = angle of rotation
C     nquad = number of quadrature nodes in equator.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     cthetas = cos(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     sthetas = sin(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     cphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C     sphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad
      real *8 cthetas(nquad)
      real *8 sthetas(nquad)
      real *8 cphis(nquad)
      real *8 sphis(nquad)
C
      pi = 4.0d0*datan(1.0d0)
C
      ca = cos(alpha)
      sa = sin(alpha)
      do i = 1,nquad
	 im1 = i-1
         phi = 2*pi*im1/nquad
	 theta = pi/2.0d0
         xp = cos(phi)*sin(theta)
         yp = sin(phi)*sin(theta)
         zp = cos(theta)
         x = ca*xp + sa*zp
         y = yp
         z = -sa*xp + ca*zp
         proj = sqrt(x**2+y**2)
	 if (proj.le.1.0d-16) then
	    cphis(i) = 1.0d0
	    sphis(i) = 0.0d0
	 else
	    cphis(i) = x/proj
	    sphis(i) = y/proj
	 endif
	 cthetas(i) = z
	 sthetas(i) = proj
      enddo
      return
      end
C
C***********************************************************************
      subroutine evalall(alpha,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     alpha    : angle of rotation about y-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
C     nterms   : order of multipole expansion
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     ephis    : work array for exp(i m phi) values
C     rat1     : work array for accelerating ynm calculation.
C     rat2     : work array for accelerating ynm calculation.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     uval(i,j) : contribution to potential 
C                 of multipole terms of order j at ith quad node.
C     uder(i,j) : contributions to theta derivative of potential
C                 of multipole terms of order j at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer ndeg,morder, nquad
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
ccc      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 ephi1,ephis(-nterms:nterms)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 uv,utheta,uphi,ztmp1,ztmp2,ztsum
      complex *16 ux,uy,uz,imag
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data imag/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      calpha = cos(alpha)
      salpha = -sin(alpha)
      call ylgndrini(nterms,rat1,rat2)
      do jj=1,nquad
	 ctheta = cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = cphis(jj)
	 sphi = sphis(jj)
         dir1 = -salpha
         dir2 = 0
         dir3 = calpha
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
ccc	 call ylgndr2s(nterms,ctheta,ynm,ynmd)
	 call ylgndr2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
         ephi1 = dcmplx(cphis(jj),sphis(jj))
	 ephis(1) = ephi1
	 ephis(-1) = dconjg(ephi1)
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	    ephis(-i) = ephis(-i+1)*dconjg(ephi1)
	 enddo
c
	 do ndeg = 0,nterms
	    uv = ynm(ndeg,0)*mpole(ndeg,0)
	    utheta=ynmd(ndeg,0)*stheta*mpole(ndeg,0)
	    uphi = 0.0d0
	    do morder = 1,ndeg
               ztmp1 = ephis(morder)*mpole(ndeg,morder)
               ztmp2 = ephis(-morder)*mpole(ndeg,-morder)
	       ztsum = ztmp1+ztmp2
	       uv = uv + stheta*ynm(ndeg,morder)*ztsum
	       utheta = utheta + ynmd(ndeg,morder)*ztsum
	       uphi = uphi - ynm(ndeg,morder)*morder*(ztmp1-ztmp2)
	    enddo
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = utheta*proj2+uphi*imag*proj1
	 enddo
      enddo
      return
      end
C
C
C
C
C
C
C***********************************************************************
      subroutine projectonynm(nquad,uval,uder,
     1           ynm,ynmd,marray,lmpn,nterms,m2,wsave,avec,bvec)
C***********************************************************************
C
C     This subroutine projects from values on equator for each multipole
C     order (uval, uder = dudthteta) 
C     onto spherical harmonics
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nquad    : number of points on equator
C     uval     : F values on equator
C     uder     : dFdtheta values on equator
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     lmpn     : leading dim of marray (must exceed nterms)
C     nterms   : order of expansion
C     m2       : NOT YET IMPLEMENTED (for reduced number of degrees in 
C                expansion (second index)
C     wsave    : work array for FFT (dimension at least 4*nquad+20)
C     avec     : work array of length nquad for FFT (complex)
C     bvec     : work array of length nquad for FFT (complex)
C---------------------------------------------------------------------
C     OUTPUT:
C
C     marray   : rotated expansion 
C
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad, norder
      complex *16 ephi,ephi1,uval(nquad,0:1)
      complex *16 uder(nquad,0:1)
      complex *16 utheta,uphi,ztmp1,ztmp2
      complex *16 alpha,beta,imag
      complex *16 marray(0:lmpn,-lmpn:lmpn)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 wsave(4*nquad+20)
      complex *16 avec(nquad)
      complex *16 bvec(nquad)
C
      data imag/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      theta = pi/2.0d0
      ctheta = cos(theta)
      stheta = sin(theta)
      h = stheta/nquad
      h = 1.0d0/nquad
      call ylgndr2s(nterms,ctheta,ynm,ynmd)
      call zffti(nquad,wsave)
      do norder=0,nterms
	 do ii = 1,nquad
	    avec(ii) = uval(ii,norder)
	    bvec(ii) = uder(ii,norder)
         enddo
	 call zfftf(nquad,avec,wsave)
	 call zfftf(nquad,bvec,wsave)
         do m = -norder,norder
	    if (m.ge.0)  alpha = avec(m+1)*h
	    if (m.lt.0)  alpha = avec(nquad+m+1)*h
	    if (m.ge.0)  beta = bvec(m+1)*h
	    if (m.lt.0)  beta = bvec(nquad+m+1)*h
            marray(norder,m) = (alpha*ynm(norder,abs(m)) -
     1        beta*ynmd(norder,abs(m))*stheta)/
     1        (ynm(norder,abs(m))**2 + (stheta*ynmd(norder,abs(m)))**2)
         enddo
      enddo
      return
      end
c
