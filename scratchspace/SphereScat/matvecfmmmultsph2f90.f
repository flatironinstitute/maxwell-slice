C  Copyright (C) 2007: Leslie Greengard
C  Contact: greengard@cims.nyu.edu
C
C  All rights reserved.
C  This copy is not for distribution without permission of the author.
C
c**************************************************************
c
      subroutine multini3dt(nsp,nterms,nquad,x0y0z0,sphererads,
     1           scale,zk0,zk,beta0,delta0,beta,delta,ier)
c
c     initialization code for matrix vector multiply.
c     Allows subsequent call with simple calling sequence.
c     Stores problem-specific data:
c
c**************************************************************
c**************************************************************
c     WARNING: Static arrays assume no more than 1000 spheres.....
c*****************************************************************
c*****************************************************************
c
c     nsp     number of spheres
c     nterms    length of expansions
c     x0y0z0    sphere positions
c     sphererads  sphere radii
c     beta(i)   coefficient of Pot(in sphere i)
c     beta0     coefficient of Pot(in exterior)
c     delta(i)  coefficient of dPotdn(in sphere i)
c     delta0    coefficient of dPotdn(in exterior)
c     zk(i)     Helmholtz coefficient of sphere i
c     zk0       Helmholtz coefficient in ext region
c
c
      implicit none
      integer *4 nsp,nsp7,nterms,nterms7,nsmax,nquad,nquad7
      integer *4 i,iflag,ier,ierloc,lenw,n
      parameter(nsmax=1000)
      real *8    bb(2000), xnodes(2000), wts(2000)
      real *8    sphererads(nsp), sphererads7(nsmax)
      real *8    scale(nsp), scale7(nsmax)
      real *8 x0y0z0(3,nsp),x0y0z07(3,nsmax)
      complex *16    beta(1), beta7(nsmax)
      complex *16    beta0, beta07
      complex *16    delta(1), delta7(nsmax)
      complex *16    delta0, delta07
      complex *16    zk(1), zk7(nsmax)
      complex *16    zk0, zk07
      complex *16    vecin(1), vecout(1)
      complex *16    work(1)
      save  xnodes,wts
      save  nsp7,nterms7,nquad7
      save  sphererads7
      save  x0y0z07,scale7
      save  beta7,beta07,delta7,delta07,zk7,zk07
c
      nquad7 = nquad
      nsp7 = nsp
      nterms7 = nterms
      zk07 = zk0
      beta07 = beta0
      delta07 = delta0
      call gaussq(0,nquad,0,0,0,0,bb,xnodes,wts)
      do i = 1,nsp
         x0y0z07(1,i) = x0y0z0(1,i)
         x0y0z07(2,i) = x0y0z0(2,i)
         x0y0z07(3,i) = x0y0z0(3,i)
         sphererads7(i) = sphererads(i)
         scale7(i) = scale(i)
         beta7(i) = beta(i)
         delta7(i) = delta(i)
         zk7(i) = zk(i)
      enddo
      return
c
      entry multfmm3dt(vecin,vecout,n,work,lenw)
c
c       actual matrix vector product, using fmm routine
c       followed by (left) block diagonal preconditioner
c
ccc	call prinf(' in multfmm3dt nquad7 is *',nquad7,1)
        call matvecfmmpc3t(vecin,vecout,nsp7,nterms7,x0y0z07,
     1       sphererads7,scale7,beta7,beta07,delta7,delta07,
     2       zk7,zk07,work,lenw,xnodes,wts,nquad7,ierloc)
        iflag = 0
        call diagprec3dt(vecout,nsp7,nterms7,x0y0z07,
     1       sphererads7,scale7,zk07,zk7,beta07,delta07,beta7,delta7,
     2       work,lenw,iflag)
      return
      end
c
c
c
c
c
c
      subroutine matvecfmmpc3t(vecin,vecout,nsp,nterms,x0y0z0,
     1           sphererads,scale,beta,beta0,delta,delta0,zk,zk0,
     2           work,lenw,xnodes,wts,nquad,ier)
      implicit none
      integer *4 iw(0:1000)
      integer *4 iused,impole,ilocal,nterms,jstart,nblock,nquad
      integer *4 ivectmp,j1,j1start,j2,j2start
      integer *4 ljtemp,ifder,ntop,next
      integer *4 lenw,nmpole,nmloc,nmwork,lenleft,iprec
      integer *4 inform(10),ier, ier2
      integer *4 kk,j,k,nsp,nn,nn2,nsys
      real *8    xnodes(1),wts(1)
      real *8    x0y0z0(3,nsp)
      real *8    sphererads(nsp)
      real *8    scale(nsp)
      complex *16 vecin(1),vecout(1)
      complex *16 work(lenw)
      complex *16 z,zk(1),zk0,wronsk
      complex *16 beta(1),delta(1),beta0,delta0
      complex *16 jtemp(0:1000)
      complex *16 jptemp(0:1000)
      complex *16 htemp(0:1000)
      complex *16 hptemp(0:1000)
      complex *16 jtemp0(0:1000)
      complex *16 jptemp0(0:1000)
      complex *16 htemp0(0:1000)
      complex *16 hptemp0(0:1000)
c
c***************************************************************
c     WARNING: lots of static arrays here. 
c              expansions assumed to be much less than 1000th order.
c*****************************************************************
c
c     This routine applies the matrix discretizing the acoustic
c     scattering integral equation based on H expansions.
c     The variables in VECIN are ordered 
c   
c        [sigma_1, sigma_2, .., sigma_NSP ]
c
c     where 
c          sigma is the vector of spherical harmonic expansion 
c     coefficient for the exterior of the corresponding sphere.
c
c     INPUT:
c
c     vecin         (complex *16)    input vector:
c                                   dimension nsp*(nterms+1)**2
c     nsp           (integer *4)     number of spheres
c     nterms        (integer *4)    order of SH expansions 
c     x0y0z0(3,i)   (real *8)       hole center coordinates
c     sphererads(i) (real *8)       radius of ith sphere
c     beta(i)       (complex *16)  jump coefficient for E on inside
c                                   of sphere i  
c     beta0         (complex *16)  jump coefficient for E in exterior
c                                   domain
c     delta(i)      (complex *16)  jump coefficient for dEdn on inside
c                                   of sphere i  
c     delta0         (complex *16)  jump coefficient for dEdn in exterior
c                                   domain
c
c     zk(i)          (complex *16)  Helmholtz frequency parameter for
c                                   sphere i  
c     zk0            (complex *16)  Helmholtz frequency parameter for
c                                   exterior  
c     work           (complex *16)  work array.
c     lenw           (integer *4)   length of work array.
c
c     OUTPUT:
c
c     vecout       (complex *16)    output vector:
c                                   dimension nsp*(nterms+1)**2
c
c-------------------------------------------------------------------
c
      ier = 0
      nn = (nterms+1)**2
      nblock = (nterms+1)*(2*nterms+1)
      impole = 1
      ilocal = impole+nblock
      ivectmp = ilocal+nblock
      iused = ivectmp+nn
      if (iused.gt.lenw) then
         ier = 1
	 return
      endif
      nn2 = 2*nn
c
c
c     First compute "self-interaction" - i.e. diagonal block.
c     Because of mathematical formalism, this block is itself
c     a diagonal matrix.
c
c     Thus, we initialize the output vector (vecout) as
c     relevant diagonal matrix times vecin.
c     
c     The subroutine getdiag3dt computes the appropriate diagonal
c     operator.
c
      do j = 1,nsp
         call getdiag3dt(work(1),nterms,sphererads(j),scale(j),
     1        zk0,zk(j),beta0,beta(j),delta0,delta(j),
     2        work(nn+1),lenw-nn)
         jstart = (j-1)*nn
         do k = 1,nn
            vecout(jstart+k) = work(k)*vecin(jstart+k)
         enddo
      enddo
ccc      call prin2(' after diag vecout *',vecout,2*nn)
c
c     use multiple scattering to get contribution from offdiagonal blocks.
c
c     Each sphere (j2) has a multipole expansion which induces a local
c     expansion on sphere (j1).
c
c     The algorithm is: 
c            
c     1) unroll input vector (h expansion coefficients for sphere j2).
c     2) translate to j-expansion on sphere j1.
c     3) apply formula (see notes).
c
      do j1 = 1,nsp
         j1start = (j1-1)*nn
         do j2 = 1,nsp
            if (j2.eq.j1) goto 333
            j2start = (j2-1)*nn
            call unrolltoshxp(nterms,vecin(j2start+1),work(impole))
ccc	    call prinm(work(1),nterms)
            call h3dmplocquadu(zk0,scale(j2),x0y0z0(1,j2),
     1           work(impole),nterms,scale(j1),x0y0z0(1,j1),
     1           work(ilocal),nterms,sphererads(j1),xnodes,wts,nquad)
ccc	    call prinm(local,nterms)
            call shxptounroll(nterms,work(ivectmp),work(ilocal))
	    ljtemp = 1000
	    ifder = 1
	    z = zk0*sphererads(j1)
            call jfuns3d(ier2,nterms,z,scale(j1),jtemp0,ifder,jptemp0,
     1           ljtemp,iw,ntop)
	    z = zk(j1)*sphererads(j1)
            call jfuns3d(ier2,nterms,z,scale(j1),jtemp,ifder,jptemp,
     1           ljtemp,iw,ntop)
	    next = 1
            do j = 0,nterms
               do k = -j,j
                  wronsk = delta0*zk0*beta(j1)*jptemp0(j)*jtemp(j)-
     1                  delta(j1)*zk(j1)*beta0*jptemp(j)*jtemp0(j)
                  vecout(j1start+next) = vecout(j1start+next) +
     1            work(ivectmp+next-1)*wronsk
                  next = next+1
               enddo
            enddo
333      continue
         enddo
      enddo
ccc      call prin2(' after off diag vecout *',vecout,2*nn)
444      continue
      return
      end
c
c
      subroutine getdiag3dt(diagvec,nterms,sphererad,scale,zk0,zk,
     1             beta0,beta,delta0,delta,work,lenw)
      implicit none
      integer *4 iw(0:1000),ier,j,k,ntop,ifexpon
      integer *4 ljtemp,nn,nn2,nterms,lenw
      integer *4 ifder,next
      real *8 sphererad,scale
      complex *16 diagvec(1)
      complex *16 z,zk,zk0,wronsk
      complex *16 wronsk0
      complex *16 beta,beta0
      complex *16 delta,delta0
      complex *16 jtemp(0:1000)
      complex *16 htemp(0:1000)
      complex *16 jptemp(0:1000)
      complex *16 hptemp(0:1000)
      complex *16 jtemp0(0:1000)
      complex *16 htemp0(0:1000)
      complex *16 jptemp0(0:1000)
      complex *16 hptemp0(0:1000)
      complex *16 work(lenw)
c
c     This routine is used to get the diagonal self-interaction matrix
c     (as a vector). 
c
c     INPUT:
c
c     nterms      (integer *4)      number of terms in expansions
c     sphererad(i)   (real *8)     radius of ith sphere
c     scale          (real *8)      scaling parameter
c     zk0            (complex *16)  ext Helmholtz frequency parameter
c     zk             (complex *16)  int Helmholtz frequency parameter
c     beta0          (complex *16)  coeff of jump in E parameter for ext.
c     beta           (complex *16)  coeff of jump in E parameter for int.
c     delta0         (complex *16)  coeff of jump in dEdn parameter for ext.
c     delta          (complex *16)  coeff of jump in dEdn parameter for int.
c     work(lenw)     (complex *16)  workspace
c
c     OUTPUT:
c
c     diagvec        (complex *16)  vector of diagonal entries
c
c-------------------------------------------------------------------
c
c**************************************************************
c     WARNING: Static arrays assume expansions less than 1000th order
c*****************************************************************
c
      ljtemp = 1000
      nn = (nterms+1)**2
c
      z = zk*sphererad
c
      ifder = 1
      call jfuns3d(ier,nterms,z,scale,jtemp,ifder,jptemp,
     1     ljtemp,iw,ntop)
c
      z = zk0*sphererad
      call h3dall(nterms,z,scale,htemp0,ifder,hptemp0)
c
      next = 1
      do j = 0,nterms
         do k = -j,j
            wronsk = delta0*zk0*beta*hptemp0(j)*jtemp(j)-
     1            delta*zk*beta0*jptemp(j)*htemp0(j)
            diagvec(next) = wronsk
            next = next+1
         enddo
      enddo
      return
      end
c
c
      subroutine diagprec3dt(vecin,nsp,nterms,x0y0z0,
     1           sphererads,scale,zk0,zk,beta0,delta0,beta,delta,
     2           work,lenw,iflag)
cc      implicit real *8 (a-h,o-z)
      implicit none
      integer *4 nterms,jstart,iflag
      integer *4 lenw,nmpole,nmloc,nmwork,lenleft,iprec
      integer *4 inform(10),ier
      integer *4 kk,j,k,nsp,nn,nsys
      real *8    x0y0z0(3,nsp)
      real *8    sphererads(nsp)
      real *8    scale(nsp)
      complex *16 det
      complex *16 vecin(1),a1,a2
      complex *16 work(lenw)
      complex *16 z,zk0,zk(nsp)
      complex *16 beta(nsp),beta0
      complex *16 delta(nsp),delta0
c
c     This routine applies the diagonal preconditioner for the
c     acoustic scattering integral equation based on H-expansions.
c   
c
c     INPUT:
c
c     vecin        (complex *16)    input vector:
c                                   dimension nsp*(nterms+1)**2
c     nsp         (integer *4)      number of spheres
c     nterms        (integer *4)    number of terms in sphere
c                                   expansions 
c     x0y0z0(3,i)     (real *8)     sphere center coordinates
c     sphererads(i)   (real *8)     radius of ith sphere
c     scale           (real *8)     scaling parameter
c
c     zk0            (complex *16)  Helmholtz frequency parameter for
c                                   exterior  
c     zk(i)          (complex *16)  Helmholtz frequency parameter for
c                                   interior of sphere i  
c     beta0          (complex *16)  coeff of jump in E parameter for ext.
c     delta0         (complex *16)  coeff of jump in dEdn parameter for ext.
c     beta           (complex *16)  coeff of jump in E parameter for int.
c     delta          (complex *16)  coeff of jump in dEdn parameter for int.
c     work           (complex *16)  work array.
c     lenw           (integer *4)   length of work array.
c     iflag          (integer *4)   0 ->  invert D, 1 ->  invert D*
c                    not yet implemented - easy here, harder for multfmm...
c
c     OUTPUT
c
c     vecin       (complex *16)     vecin is overwritten by preconditioner
c                                   applied to  vecin
c
c-------------------------------------------------------------------
c
      nn = (nterms+1)**2
c
c     set vecout =  diagonal matrix times vecin 
c
      do j = 1,nsp
         call getdiag3dt(work(1),nterms,sphererads(j),scale(j),zk0,
     1        zk(j),beta0,beta(j),delta0,delta(j),
     2        work(nn+1),lenw)
         jstart = (j-1)*nn
         do k = 1,nn
            a1 = vecin(jstart+k)
            vecin(jstart+k) = a1/work(k)
         enddo
      enddo
      return
      end
c
