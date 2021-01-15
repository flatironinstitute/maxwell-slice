c
c**********************************************************************
      subroutine h3dpwdata(ier,kvec,rscale,center,
     &		radius,nterms,pot,dpotdn,w,lw,lused)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c PURPOSE:
c
c     This subroutine creates the j-expansion at the center
c     center(3), due to a plane wave defined by kvec.
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformta0 below.
c
c INPUT:
c
c     wavek   (complex *16)     : the Helmholtz coefficient
c     rscale   (real *8)        : scaling parameter
c                                   should be less than one in magnitude.
c                                   Needed for low frequency regime only
c                                  with rsclale abs(wavek) recommended.
c     center(3)   (real *8)     : coordinates of the expansion center
c     nterms   (integer *4)     : order of the j-expansion
c     w           (real *8)     : workspace
c     lw       (integer *4)     : workspace length
c
c
c OUTPUT:
c
c     ier      (integer *4)     : error return code
c		ier=0	returned successfully;
c		ier=8	not enough memory supplied to this subroutine
c	 	ier=16  d is out of range in hfuns3d
c     pot(0:nterms,0:nterms)    : potential in s.h. basis
c              (complex *16)    :
c     dpotdn(0:nterms,0:nterms) : normal deriv of potential in s.h. basis
c              (complex *16)    :
c     lused    (integer *4)     : amount of work space "w" used
c
c
      real *8 center(3),kvec(3)
      complex *16 wavek,pot(0:nterms,-nterms:nterms)
      complex *16 dpotdn(0:nterms,-nterms:nterms)
      real *8 w(1)
c
c ... assign workspaces:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      lused=iephi+lephi
c
      if (lused.gt.lw) then
         ier=8
         return
      endif
c
      do l = 0,nterms
         do m = -l,l
            pot(l,m) = 0.0d0
            dpotdn(l,m) = 0.0d0
         enddo
      enddo
c
      call prinf(' calling pwtasc0 *',ier,0)
      call h3dpwdata0(jer,kvec,rscale,center,radius,
     &		nterms,pot,dpotdn,w(ipp),w(iephi))
      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dpwdata0(ier,kvec,rscale,
     &		center,radius,nterms,pot,dpotdn,pp,ephi)
      implicit real *8 (a-h,o-z)
c
c     See h3dpwta for comments
c
c**********************************************************************
c
      integer *4 iw(0:2000)
      real *8 kvec(3),center(3),zdiff(3)
      complex *16 wavek,pot(0:nterms,-nterms:nterms)
      complex *16 dpotdn(0:nterms,-nterms:nterms)
c
      real *8 pp(0:nterms,0:nterms)
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 ztmp,z,zconst
      complex *16 fjs(0:2000)
      complex *16 fjder(0:2000)
      complex *16 cscale,eye
c
      data eye/(0.0d0,1.0d0)/
      data thresh/1.0d-15/
c
c ... Initializing...
c
      ier=0
c
      zdiff(1) = kvec(1)
      zdiff(2) = kvec(2)
      zdiff(3) = kvec(3)
      wavek= dsqrt( kvec(1)**2 + kvec(2)**2 + kvec(3)**2 )
      cdot = kvec(1)*center(1) + kvec(2)*center(2) + kvec(3)*center(3)
      zconst = cdexp(eye*cdot)/(eye*wavek)
      call prin2(' zconst is *',zconst,2)
      call cart2polar(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1=dcmplx(cphi,sphi)
c
c ... get all the e^{eye*m*phi}:
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      do 1040 i=2,nterms
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
 1040 continue
c
c ... get the associated Legendre functions:
c
      call ylgndr(nterms,ctheta,pp)
c
c ... finally, calculate the coeff locexp:
c
      lwfjs = 2000
      ifder = 1
      z = wavek*radius
      call jfuns3d(jer,nterms,z,rscale,fjs,ifder,fjder,
     1             lwfjs,iw,ntop)
      cscale = rscale
      pot(0,0)=pot(0,0) + zconst*cscale*fjs(0)
      dpotdn(0,0)=dpotdn(0,0) + zconst*cscale*fjder(0)*wavek
      do 2040 n=1,nterms
         cscale = eye*cscale*rscale
         pot(n,0)=pot(n,0) + pp(n,0)*zconst*cscale*fjs(n)
         dpotdn(n,0)=dpotdn(n,0) + pp(n,0)*zconst*cscale*fjder(n)*wavek
         do 2020 m=1,n
            ztmp=pp(n,m)*zconst
	    pot(n,m)=pot(n,m) + ztmp*ephi(-m)*cscale*fjs(n)
	    pot(n,-m)=pot(n,-m) + ztmp*ephi(m)*cscale*fjs(n)
	    dpotdn(n,m)=dpotdn(n,m) + 
     1 	    ztmp*ephi(-m)*cscale*fjder(n)*wavek
	    dpotdn(n,-m)=dpotdn(n,-m) + 
     1	    ztmp*ephi(m)*cscale*fjder(n)*wavek
 2020    continue
 2040 continue
      return
      end
c
c
