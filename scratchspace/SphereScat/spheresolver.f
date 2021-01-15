C  Copyright (C) 2007: Leslie Greengard
C  Contact: greengard@cims.nyu.edu
C
C  All rights reserved.
C  This copy is not for distribution without permission of the author.
C
C**********************************************************************C
      subroutine spheresolver(nsphere,nterms,nquad,nn1,nn,x0y0z0,
     1     sphererads,scale,zk0,zk,beta0,delta0,beta,delta,
     1     multfmm3dt,rhs,rhs2,rhsjoin,eps,numit,sol,
     1     niter,errs,ngmrec,cw,www,lcw,ierr)
c
c     See Notes for a discussion of mathematical approach.
c    
c     Step 1: create rhs for linear system by combining the 
c     jump condition data in rhs and rhs2 into a single vector.
c
c     INPUT:
c 
c     nsphere  =    number of spheres
c     nterms   =    order of expansions (currently the same for all
c                     spheres)
c     nquad    =    number of quadrature nodes used in projection
c                   2*nterms is generally sufficient.
c     nn1      =    number of unknowns associated with a single sphere
c     nn       =    total number of unknowns 
c     x0y0z0   =    sphere centers
c     sphererads   =  sphere radii
c     scale   =       scale parameters for each sphere
c     zk0     =       exterior Helmholtz coefficient
c     zk      =       vector of interior Helmholtz coefficients for
c                     all spheres
c     beta0, delta0  =  material coefficients for exterior region.
c     beta, delta  =  material coefficients for interior regions.
c
c                      [beta phi] = rhs
c                      [delta phi_n] = rhs2
c
c     multfmm3dt =  user-supplied subroutine to carry out matrix-vector
c                   products
c     rhs    =    first jump condition [beta phi]
c     rhs2   =    second jump condition [delta phi_n]
c     eps    =    tolerance for GMRES scheme.
c     numit  =    maximum number of iterations allowed for GMRES.
c     ngmrec =    maximum number of iterations before restarting GMRES.
c     cw     =    workspace
c     www    =    second workspace
c     lcw    =    length of BOTH workspaces.
c
c     OUTPUT:
c
c     rhsjoin = right-hand side for scattering solver - a linear 
c               combination of rhs and rhs2.
c
c-------------------------------------------------------                  
c     sol    =     solution vector.
c
c                  i.e. compact (unrolled) H-expansion coefficients
c                  for spheres 1,2,...,nsphere
c-------------------------------------------------------                  
c
c     niter  =     number of iterations used.
c     errs  =      errors at each iteration.
c     ierr   =     error return code, not fully implemented
c
c
      implicit none
      integer *4 lcw,nterms,nsphere,nn,nn1,ier,ierr
      integer *4 numit,ngmrec,iflag,lw,niter,i,iffld,lused
      integer *4 nquad,nsp,ltemp,ifder,ntop,next,j,k
      real *8 sphererads(nsphere)
      real *8 x0y0z0(3,nsphere)
      complex *16 zk(nsphere)
      complex *16 beta(nsphere),delta(nsphere)
c
      complex *16 sol(nn)
      complex *16 rhs(nn),rhs2(nn),rhsjoin(nn)
      complex *16 beta0,delta0,zk0,z
      complex *16 cw(lcw),eye
      complex *16 www(lcw),hout
      real *8 targ(3),scale(nsphere),eps
      real *8 errs(1000)
      external multfmm3dt
      data eye/(0.0d0,1.0d0)/
c
      ierr = 0
c
c     Combine f and g according to formula to construct
c     right-hand side for integral equation and store in vector
c     rhs.
c
      do nsp = 1,nsphere
         z = zk(nsp)*sphererads(nsp)
         ltemp = (lcw-2)/3
         ifder = 1
         call jfuns3d(ier,nterms,z,scale(nsp),cw(1),ifder,
     1             cw(ltemp+1),ltemp,cw(2*ltemp+1),ntop)
         next = (nsp-1)*nn1+1
         do j = 0,nterms
            do k = -j,j
               hout  = beta(nsp)*cw(j+1)*rhs2(next) -
     1                 delta(nsp)*zk(nsp)*cw(ltemp+1+j)*rhs(next)
               rhsjoin(next)  = hout
	       next = next+1
            enddo
         enddo
      enddo
c
c    now solve iteratively using GMRES with input parameters.
c
      iflag = 0
      call diagprec3dt(rhsjoin,nsphere,nterms,x0y0z0,sphererads,
     1     scale,zk0,zk,beta0,delta0,beta,delta,cw,lw,iflag)
ccc      call prin2('after diagprec, rhs=*',rhs,2*nn)
      call multini3dt(nsphere,nterms,nquad,x0y0z0,sphererads,
     1     scale,zk0,zk,beta0,delta0,beta,delta,ier)
      call cgmrespc(ier,nn,multfmm3dt,rhsjoin,eps,numit,sol,
     1     niter,errs,ngmrec,cw,www,lcw)
ccc      call prin2(' errs is *',errs,niter)
      return
      end
