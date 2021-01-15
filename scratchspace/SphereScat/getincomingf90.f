C  Copyright (C) 2007: Leslie Greengard
C  Contact: greengard@cims.nyu.edu
C
C  All rights reserved.
C  This copy is not for distribution without permission of the author.
C
c***********************************************************************
      subroutine phi_incoming(str,kvec,targ,phi,gradphi)
c***********************************************************************
c
c     This subroutine defines incoming field as 
c  
c     phi = str * exp( i kvec \cdot targ)
c
c     INPUT:
c
c     str        (complex *16) : strength parameter
c     kvec       (complex *16) : vector parameter
c     targ           (real *8) : coordinates of sample point
c
c     OUTPUT:
c
c     phi        (complex *16) : potential
c     gradphi    (complex *16) : gradient of potential
c
c***********************************************************************
c
c
      implicit real *8 (a-h,o-z)
      real *8 dot,kvec(3),targ(3)
      complex *16 phi,gradphi(3),eye, str
      data eye/(0.0d0,1.0d0)/
C
      dot = kvec(1)*targ(1)+kvec(2)*targ(2)+kvec(3)*targ(3)
      phi = str*cdexp( eye*dot)
      gradphi(1) = phi*eye*kvec(1)
      gradphi(2) = phi*eye*kvec(2)
      gradphi(3) = phi*eye*kvec(3)
      return
      end

c***********************************************************************
      subroutine getincoming(str,kvec,radius,center,nterms,
ccc     1           sh1,sh2,w,lw,ier)
     1           sh1,sh2)
c***********************************************************************
c
c     This subroutine expands an arbitrary incoming field
c     defined by subroutine 
c             
c        phi_incoming(str,kvec,targ,phi,gradphi)
c
c     on the sphere of radius RADIUS as a spherical harmonic expansion.
c
c     INPUT:
c
c     str        (complex *16) : strength parameter
c     kvec       (complex *16) : vector parameter
c     radius      (real *8 )   : radius of sphere about (0,0,zshift)
c                                where phival, phivaln are computed.
c     center      (real *8)    : coordinates of sphere center
c     nterms     (integer *4 ) : order of expansion
c     w          (complex *16) : workspace
c     lw          (integer *4) : length of workspace
c
c     OUTPUT:
c
c     sh1  (complex *16)  :    sph. harmonic expansion of potential
c     sh2  (complex *16)  :    sph. harmonic expansion of normal 
c                              derivative of potential
c     ier                 :    error return code
c                              not yet implemented
c
c
c     WARNING:  xnodes, wts, bb are dimensioned as local variables
c               at length 2000  => nterms < 1000.
c
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer *4 nterms
ccc      real *8 w(lw),xnodes(2000),wts(2000),bb(2000)
      real *8 targ(3), center(3)
      real *8 kvec(3)
      complex *16 sh1(0:nterms,-nterms:nterms)
      complex *16 sh2(0:nterms,-nterms:nterms)
      real *8 :: w
      allocatable w(:)
C
      integer *4 l,m,jnew,knew
C
      pi = 4.0d0*datan(1.0d0)
      nquad = 2*nterms
      iphival = 1
      iphivaln = iphival + 2*nquad*nquad
      iynm = iphivaln + 2*nquad*nquad
      ixnodes = iynm + (nterms+1)**2
      iwts = ixnodes+ nquad
      ibb = iwts+ nquad
      iused = ibb+nquad
      lw = iused+ 100000
ccc      if (iused. lt. lw) then
ccc         ier = 1
ccc	 return
ccc      endif
      allocate (w(lw))
      call gaussq(0,nquad,0,0,0,0,w(ibb),w(ixnodes),w(iwts))
c
      call sampleincoming(str,kvec,radius,center,nquad,w(ixnodes),
     1                    w(iphival),w(iphivaln))
c
      call h3dprojloc(nterms,nterms,nquad,w(ixnodes),w(iwts),
ccc     1                w(iphival),sh1,w(iynm),w(iused),lw-iused,ier)
     1                w(iphival),sh1,w(iynm))
ccc      call prin2(' after projloc ier is *',ier,1)
      call h3dprojloc(nterms,nterms,nquad,w(ixnodes),w(iwts),
ccc     1                w(iphivaln),sh2,w(iynm),w(iused),lw-iused,ier)
     1                w(iphivaln),sh2,w(iynm))
ccc      call prin2(' after second projloc ier is *',ier,1)
      deallocate (w)
      return
      end
c
c***********************************************************************
      subroutine sampleincoming(str,kvec,radius,center,nquad,xnodes,
     1           phival,phivaln)
c***********************************************************************
c
c     This subroutine samples an incoming field
c     defined by subroutine phi_incoming(x,y,z)
c     on the sphere of radius RADIUS centered at CENTER.
c
c     INPUT:
c
c     str        (complex *16) : complex parameter to send to 
c                                subroutine phi_incoming.
c     kvec        (complex *16) : complex 3-vector to send to 
c                                subroutine phi_incoming.
c     radius      (real *8 ) : radius of sphere
c                              where phival, phivaln are computed.
c     center      (real *8 ) : sphere center
c     nquad    (integer *4 ) : number of quadrature nodes
c                              on target sphere is nquad*nquad.
c     xnodes      (real *8 ) : Legendre nodes x_j = cos theta_j.
c
c     OUTPUT:
c
c     phival  (complex *16)  : value of potential on tensor product
c                              mesh on target sphere.
c     phivaln  (complex *16) : value of normal derivative of potential 
c                              on tensor product mesh on target sphere.
c
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer *4 nterms
      real *8 targ(3), center(3)
      real *8 xnodes(1),kvec(3)
      complex *16 phival(nquad,nquad)
      complex *16 phivaln(nquad,nquad)
      complex *16 pot,gradphi(3),str,eye
      data eye/(0.0d0,1.0d0)/
C
      integer *4 l,m,jnew,knew
C
      pi = 4.0d0*datan(1.0d0)
      iffld = 0
      do jj=1,nquad
      do kk=1,nquad
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
	 phi = 2*pi*kk/nquad
	 cosphi = dcos(phi)
	 sinphi = dsin(phi)
	 rx = stheta*cosphi
	 ry = stheta*sinphi
	 rz = ctheta
	 targ(1) = center(1) + rx*radius
	 targ(2) = center(2) + ry*radius
	 targ(3) = center(3) + rz*radius
         call phi_incoming(str,kvec,targ,pot,gradphi)
c
c        scale by radius/eye = -eye*radius
c        has something to do with local expansion evaluation.
c        may ttrack this down more carefully (or not).
c
         phival(jj,kk) = -eye*pot*radius
         phivaln(jj,kk) = -eye*radius*
     1   	 (gradphi(1)*rx+gradphi(2)*ry+gradphi(3)*rz)
      enddo
      enddo
      return
      end
c
