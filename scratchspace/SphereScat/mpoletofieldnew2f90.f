C  Copyright (C) 2007: Leslie Greengard
C  Contact: greengard@cims.nyu.edu
C
C  All rights reserved.
C  This copy is not for distribution without permission of the author.
C
C**********************************************************************C
      subroutine mpoletofield(vecin,vecout,rhs,beta0,beta,nsp,nterms,
     1           x0y0z0,sphererads,scale,zk0,zk,work,lenw,ier)
      implicit none
      integer *4 iw(0:1000)
      integer *4 nquad,nterms,jstart,ldc,nblock,ldcblock,nsp
      integer *4 ivectmp,impole,ilocal,iused
      integer *4 j1,j1start,j2,j2start
      integer *4 ljtemp,ifder,ntop,next
      integer *4 lenw,nmpole,nmloc,nmwork,lenleft,iprec
      integer *4 inform(10),ier, ier2
      integer *4 kk,j,k,nn,nn2,nsys
      real *8    x0y0z0(3,nsp)
      real *8    sphererads(nsp)
      real *8    scale(nsp)
      real *8    bb(2000),xnodes(2000),wts(2000)
      complex *16 rhs(1)
      complex *16 vecin(1)
      complex *16 vecout(1)
      complex *16 work(lenw)
      complex *16 beta0,beta(nsp)
      complex *16 z,zk0,wronsk,zk(1)
      complex *16 jtemp0(0:1000)
      complex *16 jptemp0(0:1000)
      complex *16 jtemp(0:1000)
      complex *16 jptemp(0:1000)
cc      complex *16 vectmp(100 000)
cc      complex *16 local(100 000)
c
c**********************************************************************c
c     WARNING: various one-dimensional vectors are allocated here 
c     irresponsibly and dimensioned 1000 - code will break for 
c     very high order expansion lengths. To fix, change 1000 to 
c     a suitably larger number or do more memory allocation work...
c**********************************************************************c
c
c     This routine takes the multipole expansions on all spheres and
c     computes the induced local expansions on all spheres.
c     I.e. H expansions  -> local expansions of fields
c     The variables in VECIN are ordered 
c   
c        [mpole_1, mpole_2, .., mpole_NSP ]
c
c     where the multipole coeffs are ordered in compressed (linear) 
c     fashion and the output local expansions are also stored in 
c     compressed (linear) order.
c
c     INPUT:
c
c     vecin         (complex *16)    input vector:
c                                   dimension nsp*(nterms+1)**2
c     nsp           (integer *4)     number of spheres
c     nterms        (integer *4)    order of SH expansions 
c     x0y0z0(3,i)   (real *8)       hole center coordinates
c     sphererads(i) (real *8)       radius of ith sphere
c     zk(i)          (complex *16)  Helmholtz frequency parameter int of
c                                   sphere i  
c     zk0            (complex *16)  Helmholtz frequency parameter for
c                                   exterior  
c     work           (complex *16)  work array.
c     lenw           (integer *4)   length of work array.
c
c     OUTPUT:
c
c     vecout       (complex *16)    local exp of potential on each sphere
c                                   this is NOT the field - must be
c                                   multiplied by j_n to get correct value.
c                                   (i.e. can subsequrntly use taeval for 
c                                    evaluation inside)
c                                   dimension nsp*(nterms+1)**2
c
c-------------------------------------------------------------------
c
      nquad = 2*nterms
      ljtemp = 1000
      call gaussq(0,nquad,0,0,0,0,bb,xnodes,wts)
      ier = 0
      impole = 1
      nn = (nterms+1)**2
      nblock = (nterms+1)*(2*nterms+1)
      ilocal = impole + nblock
      ivectmp = ilocal + nblock
      iused = ivectmp + nn
      if (iused.gt.lenw) then
         ier = 2
	 return
      endif
c
      nn2 = 2*nn
c
c     set vecout =  diagonal matrix times vecin 
c
      do j = 1,nsp
         jstart = (j-1)*nn
cc         call getselffield(vecin(jstart+1),rhs(jstart+1),beta0,beta(j),
cc     1        nterms,sphererads(j),scale(j),zk0,zk(j),work(1))
         call getselffield(vecin(jstart+1),
     1        nterms,sphererads(j),scale(j),zk0,work(1))
         do k = 1,nn
            vecout(jstart+k) = work(k)
         enddo
      enddo
ccc      call prin2(' after diag vecout *',vecout,2*nn)
c
c     use translation operators to get contribution from offdiagonal blocks.
c     and to create j-expansion on each sphere corresponding to 
c     EXTERIOR Helmholtz parameter.
c
      do j1 = 1,nsp
         j1start = (j1-1)*nn
         do j2 = 1,nsp
            if (j2.eq.j1) goto 333
            j2start = (j2-1)*nn
            call unrolltoshxp(nterms,vecin(j2start+1),work(1))
cc            call h3dmplocquad(zk0,x0y0z0(1,j2),x0y0z0(1,j1),
cc     1           work(1),work(ilocal),nterms,nterms,
cc     1           scale(j2),scale(j1),sphererads(j1),xnodes,wts,
cc     1           nquad)
            call h3dmplocquadu(zk0,scale(j2),x0y0z0(1,j2),
     1           work(1),nterms,scale(j1),x0y0z0(1,j1),work(ilocal),
     1           nterms,sphererads(j1),xnodes,wts,nquad)
            call shxptounroll(nterms,work(ivectmp),work(ilocal))
cc            call shxptounroll(nterms,vectmp,local)
	    ifder = 1
	    z = zk0*sphererads(j1)
	 ljtemp = 1000
            call jfuns3d(ier,nterms,z,scale(j1),jtemp0,ifder,jptemp0,
     1           ljtemp,iw,ntop)
	    next = 1
            do j = 0,nterms
               do k = -j,j
                  vecout(j1start+next) = vecout(j1start+next) +
     1            work(ivectmp+next-1)*jtemp0(j)
ccc     1            vectmp(next)*jtemp0(j)
                  next = next+1
               enddo
            enddo
333      continue
         enddo
      enddo
c
c     Now convert from field values ON each sphere to j-expansion
c     coefficients in interior using continuity conditions (see Notes).
c
ccc      call prin2(' after off diag vecout *',vecout,2*nn)
      do j = 1,nsp
         z = zk(j)*sphererads(j)
	 ljtemp = 1000
         call jfuns3d(ier,nterms,z,scale(j),jtemp,ifder,jptemp,
     1           ljtemp,iw,ntop)
         jstart = (j-1)*nn
         next = 1
         do j1 = 0,nterms
            do k = -j1,j1
               vecout(jstart+next)=
     1 	        (beta0*vecout(jstart+next)-rhs(jstart+next))/
     2              (beta(j)*jtemp(j1))
               next = next+1
            enddo
         enddo
      enddo
      return
      end
c

      subroutine getselffield(diagvec,nterms,
     1             sphererad,scale,zk0,localvec)
      implicit none
      integer *4 iw(0:1000),ier,j,k,ntop,ifexpon
      integer *4 ljtemp,nn,nn2,nterms,lenw
      integer *4 ifder,next
      real *8 sphererad,scale
      complex *16 diagvec(1)
      complex *16 rhs(1)
      complex *16 z,zk0,zk,wronsk
      complex *16 htemp0(0:1000)
      complex *16 hptemp0(0:1000)
ccc      complex *16 jtemp(0:1000)
ccc      complex *16 jptemp(0:1000)
      complex *16 localvec(1)
c
c     This routine is used to compute the field on a sphere from its
c     own H-expansion (as a spherical harmonic expansion). 
c
c     INPUT:
c
c     diagvec     (complex *16)  coefficients of H expansion on sphere
c     nterms      (integer *4)      number of terms in expansions
c     sphererad   (real *8)     radius of sphere
c     scale       (real *8)     scale parameter for expansion
c     zk0            (complex *16)  ext Helmholtz frequency parameter
c
c     OUTPUT:
c
c     localvec    (complex *16)  value of field ON sphere
c
c-------------------------------------------------------------------
c
c**********************************************************************c
c     WARNING: various one-dimensional vectors are allocated here 
c     irresponsibly and dimensioned 1000 - code will break for 
c     very high order expansion lengths. To fix, change 1000 to 
c     a suitably larger number or do more memory allocation work...
c**********************************************************************c
c
      z = zk0*sphererad
      ifder = 1
      call h3dall(nterms,z,scale,htemp0,ifder,hptemp0)
c
      next = 1
      do j = 0,nterms
         do k = -j,j
            localvec(next) = diagvec(next)*htemp0(j)
            next = next+1
         enddo
      enddo
      return
      end
c
