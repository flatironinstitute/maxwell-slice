C***********************************************************************
      subroutine h3dmplocquadu(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,
     2           radius,xnodes,wts,nquad)
C***********************************************************************
C
C     Memory management wrapper for subroutine h3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           sc1     = scaling parameter for mpole expansion
C           x0y0z0 = center of original multiple expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           sc2     = scaling parameter for local expansion
C           xnynzn = center of shifted local expansion
C           nterms2 = order of local expansion
C           radius  = radius of sphere on which local expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta direction.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of shifted local expansion
C
C           lused = amount of workspace w actually used.
C           ier   = error return flag
C
C                   1 insufficient workspace for dc, marray, rd, 
C                                  ephi arrays in this routine
C                                  before calling h3dmplocquad0.
C                   8   lwfjs is insufficient in 
C                       h3dmploczshiftstab/h3drescalestab
C                       for jfuns3d.
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C
C     WORK ESTIMATE:
C
C     Let nmax = max(nterms,nterms2). Then
C     w should be at least 14*(nmax)^2 + 43*(nmax) 
C                          + 16*nquad*nmax + 16*nquad + 4000
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      real *8 w
      allocatable w(:)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 imag,wavek
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      ird1 = imarray+lmarray
      lrd = (ldc+1)**2 
      ird2 = ird1+lrd
      iephi = ird2+lrd
      lephi = 2*(2*ldc+3) + 3 
      iephi2 = iephi+lephi
      iynm = iephi2+lephi
      lynm = (ldc+1)**2
      iynmd = iynm+lynm
      imp2 = iynmd+lynm
      iphitemp = imp2+(ldc+1)*(2*ldc+1)*2
      lphitemp = nq*(2*ldc+1)*2
      iphitempn = iphitemp+lphitemp
      ifhs = iphitempn+lphitemp
      ifhder = ifhs+ 2*(nterms+1) + 3
      ifjs = ifhder+ 2*(nterms+1) + 3
      lwfjs = nterms2+1000
      lfjs = 2*(lwfjs+1) + 3
      ifjder = ifjs+lfjs
      lfjder = 2*(nterms2+1)+3
      iiscale = ifjder+lfjder
      liscale = (lwfjs+1)+3
      ictheta = iiscale+liscale
      istheta = ictheta+ nq
      icphi = istheta+ nq
      isphi = icphi+ nq
      iwsave = isphi+ nq
      iavec = iwsave+ 4*nq+20
      ibvec = iavec+ 2*nq
      lused = ibvec+ 2*nq
      allocate (w(lused)) 
c
      call h3dmplocquad0(wavek,sc1,x0y0z0,mpole,nterms,sc2,xnynzn,
     1         local,nterms2,w(imarray),ldc,w(ird1),w(ird2),
     2         w(iephi),w(iephi2),radius,xnodes,wts,nquad,nq,
     3         w(iynm),w(iynmd),w(ictheta),w(istheta),
     3         w(icphi),w(isphi),w(iwsave),w(iavec),w(ibvec),w(imp2),
     4         w(iphitemp),w(iphitempn),w(ifhs),w(ifhder),
     4         w(ifjs),w(ifjder),w(iiscale),lwfjs,ier)
      deallocate(w) 
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmplocquad0(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,marray,ldc,rd1,rd2,ephi,ephi2,
     2           radius,xnodes,wts,nquad,nq,ynm,ynmd,
     3           cthetas,sthetas,cphis,sphis,wsave,avec,bvec,mp2,
     4           phitemp,phitempn,fhs,fhder,fjs,fjder,iscale,lwfjs,ier)

C***********************************************************************

C     USAGE:
C
C           Convert multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted local expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of local expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for local expansion
C           radius  = radius of sphere on which local expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes used (really nquad**2)
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of shifted local expansion
c           ier      : error return code
c              8      lwfjs insufficient for jfuns3d in h3drescale.
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           ldc      must exceed max(nterms,nterms2).
C           rd1,rd2  work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C           w       = work array 
C
C           LOTS MORE
C
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      integer  iscale(0:lwfjs)
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1),rvec(3)
      real *8 rd1(0:ldc,0:ldc)
      real *8 rd2(0:ldc,0:ldc)
      real *8 zshift
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      real *8 cthetas(nq),sthetas(nq)
      real *8 cphis(nq),sphis(nq)
      real *8 wsave(4*nq+20)
      complex *16 avec(nq),bvec(nq)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitempn(nq,-ldc:ldc)
      complex *16 mp2(0:ldc,-ldc:ldc)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
      complex *16 ephi(-ldc-1:ldc+1),imag
      complex *16 ephi2(-ldc-1:ldc+1)
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polar(rvec,d,theta,phi)
c
      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c     create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c     a rotation of THETA radians about the Yprime axis after PHI
c     radians about the z-axis.
c     The PHI rotation is carried out on the fly by multiplying 
c     mpole and ephi inside the following loop. 
c
      do l=0,nterms
         do mp=-l,l
            mpole(l,mp)  = mpole(l,mp)*ephi(mp)
         enddo
      enddo
c
      nquse = 2*nterms+2
      call rotviaproj0(theta,nquse,nterms,nterms,nterms,mpole,nterms,
     1     marray,ldc,cthetas,sthetas,cphis,sphis,ynm,ynmd,rd1,rd2,
     2     phitemp,phitempn,ephi2,wsave,avec,bvec)
c
c      Undo PHI rotation of mpole, so that it is unchanged on return.
c
      do l=0,nterms
         do mp=-l,l
            mpole(l,mp)  = mpole(l,mp)*ephi(-mp)
         enddo
      enddo
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmploczshiftstab(wavek,marray,sc1,ldc,nterms,local,
     1      sc2,nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder,
     3      iscale,lwfjs,ier)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      nquse = 2*nterms2+2
      call rotviaproj0(-theta,nquse,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc,cthetas,sthetas,cphis,sphis,ynm,ynmd,
     2     rd1,rd2,phitemp,phitempn,ephi2,wsave,avec,bvec)
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c***********************************************************************
      subroutine h3dmploczshiftstab(zk,mpole,scale,lmp,nterms,local,
     1      scale2,lmpn,nterms2,radius,zshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder,
     3      iscale,lwfjs,ier)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
C---------------------------------------------------------------------
c     INPUT:
c
c     zk       : Helmholtz parameter
c     mpole    : coefficients of original multipole exp.
c     scale    : scale parameter for mpole
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in original expansion
c
c     scale2   : scale parameter for local
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     radius   : radius of sphere about new center on which field
c                is evaluated
c     zshift   : shifting distance along z-axis
c                             (always assumed positive)
C     xnodes  = Legendre nodes (precomputed)
C     wts     = Legendre weights (precomputed)
C     nquad   = number of quadrature nodes in theta direction.
c     w(lw)    : workspace
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
c     ier      : error return code
c              8      lwfjs insufficient for jfuns3d in h3drescale.
c
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nterms2,nquad,ier
      integer l,lw,m,jnew,knew
      integer iscale(0:lwfjs)
      real *8 zshift
      real *8 xnodes(1),wts(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 phitempn(nquad,-nterms:nterms)
      complex *16 mp2(0:lmpn,-lmpn:lmpn)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 mpole(0:lmp,-lmp:lmp),zk
      complex *16 local(0:lmpn,-lmpn:lmpn)
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
      call h3dmpevalspherenmstab(mpole,zk,scale,zshift,radius,
     2     nterms,lmp,ynm,ynmd,phitemp,phitempn,nquad,xnodes,
     3     fhs,fhder)
      call h3dprojlocsepstab(nterms2,lmpn,nquad,nterms,xnodes,wts,
     1     phitemp,phitempn,local,mp2,ynm)
      call h3drescalestab(nterms2,lmpn,local,mp2,radius,zk,scale2,
     2     fjs,fjder,iscale,lwfjs,ier)
      return
      end
C
C
