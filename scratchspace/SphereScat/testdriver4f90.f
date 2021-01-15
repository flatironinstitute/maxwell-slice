C  Copyright (C) 2007: Leslie Greengard
C  Contact: greengard@cims.nyu.edu
C
C  All rights reserved.
C  This copy is not for distribution without permission of the author.
C
C**********************************************************************C
c
c       Driver for acoustic transmission problem in system of spheres.
c       The governing equation is the Helmholtz equation for a 
c       scalar phi with a different Helmholtz coefficient in each 
c       subdomain and the jump/continuity
c       conditions
c
c       (1A)   [phi] = f
c       (1B)   [phi_n] = g
c
c       where phi_n means the normal derivative of phi and [phi] is the
c       jump in phi across an interface.
c
c       In the true scattering problem, an incoming field phi^in is 
c       specified in the exterior domain and we seek to solve
c    
c              [phi^scat] = -[phi^in]
c              [phi_n^scat]    = -[phi_n^in]
c
c       (This is not the most general boundary condition, but an easy
c       model to specify.)
c       There may be some parameters beta, delta scatteered throughout
c       left over from a more general code with 
c              [beta phi] = f
c              [delta phi_n] = g
c       as boundary conditions. They're set to one here in this driver.
c
c       A test solution to verify accuracy/resolution is obtained
c       using the method of "artificial sources".
c       For this, phi is defined in each distinct domain in terms of
c       artificial sources. In particular,
c       the field inside each sphere j is that due to source numbered
c        "0" which is located outside all spheres (with the sphere-specific
c       Helmholtz parameter), and the field in the 
c       exterior region is that due to a set of sources numbered 
c       j = 1,...,nsource with source j inside sphere j and the exterior 
c       Helmholtz parameter. This global solution is denoted by phi_test.
c       
c       We can then compute the
c       right hand sides f and g corresponding to phi_test, namely
c
c              f = [phi^test]
c              g = [phi_n^test]
c
c       which is a special case of (1A), (1B)).
c       We can then solve for phi with this boundary data 
c       and see if it matches the known (hidden) solution at arbitrary 
c       test locations.
c 
c       The algorithm is based on defining phi_j in sphere j by 
c
c         phi_j = (phi_0 - f)
c
c       and finding phi_0 (the exterior solution) which satisfies: 
c
c       (d/dn) phi_0 - (d/dn) phi_j = g 
c
c*****************************************************************
c       Convention [.]  means outer value minus inner value.
c*****************************************************************
c
c       This code is limited in its capabilities,
c       max of 20 spheres set below and expansions of order
c       at most 200 also set below. Lots to clean up...
c
        implicit none
        integer *4 lcw,lsol,maxsphere,ntmax,ntargs,n1
        integer *4 maxit,nterms,nsphere,nn1,nn,nsource,iffld,ier
        integer *4 numit,ngmrec,nquad,niter,i,j,k,next,lused,ierr
        integer *4 nnx,nny,nnz,ifin,ifpr_geom,iscat
        parameter(lcw=2000 000)
        parameter(maxsphere=20)
        parameter(ntmax=200)
        parameter(lsol=(ntmax+1)*(ntmax+1)*maxsphere)
c
        real *8 sphererads(maxsphere)
        real *8 source(3,0:maxsphere)
        real *8 centers(3,maxsphere)
        real *8 scale(maxsphere)
        complex *16 str(0:maxsphere)
        complex *16 zk(maxsphere)
        complex *16 beta0,delta0
        complex *16 beta(maxsphere), delta(maxsphere)
c
        real *8 kvec(3),eps,t0,t1,dot,xxx(3)
        real *8 second,hhx,hhy,hhz,xlow,ylow,zlow
        real *8 xhigh,yhigh,zhigh
        integer *4 iw(lcw)
        complex *16 sol(lsol)
        complex *16 rhs(lsol), rhs2(lsol), rhsjoin(lsol)
        complex *16 sh1((ntmax+1)*(ntmax+1))
        complex *16 sh2((ntmax+1)*(ntmax+1)*maxsphere)
        complex *16 cw(lcw),www(lcw),eye
        complex *16 z,zk0,fld2(3,maxsphere+1)
        complex *16 hout
        real *8, allocatable :: errs(:)
        real *8, allocatable :: targ(:,:)
        complex *16, allocatable :: pot(:)
        complex *16, allocatable :: fld(:,:)
        complex *16, allocatable :: htot(:)
        external multfmm3dt
        data eye/(0.0d0,1.0d0)/
        call prini(6,13)
C
C       Some parameters
C
c     nterms = order of multipole expansions
c     nsphere = number of spheres
c     nn1 = number of unknowns on single sphere.
c     nn = total number of unknowns 
c     zk0,zk(j) = Helmholtz parameters
c     scale = scaling parameter for numerical stability, used
c             in scaling j_n and h_n functions - documented 
c             elsewhere. 
c             Strong suggestion: set scale(j) = abs(zk0*sphererads(j))
c                                if that product is << 1.
c     sphererads(i) = radius of ith sphere
c     targ = output points 
c     centers = sphere centers
c     source = source locations   0 - exterior, others inside spheres
c     ns = number of sources
c     str = strength of sources
c
c**********************************************************************
c     For now, make the workspaces iw, cw, www large - the memory 
c     managment is still under contruction.
c**********************************************************************
c
c     ifin not yet used
c     
c
      open(unit = 55,file = 'controls.dat',status='unknown')
      read(55,*) ifin
      read(55,*) iscat
      read(55,*) xxx(1)
      read(55,*) xxx(2)
      read(55,*) xxx(3)
      read(55,*) ifpr_geom
      read(55,*) maxit
      read(55,*) eps
      read(55,*) xlow
      read(55,*) ylow
      read(55,*) zlow
      read(55,*) xhigh
      read(55,*) yhigh
      read(55,*) zhigh
      read(55,*) nnx
      read(55,*) nny
      read(55,*) nnz
c
      allocate(targ(3,nnx*nny*nnz))
      allocate(pot(nnx*nny*nnz))
      allocate(fld(3,nnx*nny*nnz))
      allocate(htot(nnx*nny*nnz))
      allocate(errs(maxit))
c
      nterms = 10
      nn1 = (nterms+1)**2
      call getspheredat(nsphere,sphererads,centers,kvec,zk)
c
      zk0 = dsqrt( kvec(1)**2 + kvec(2)**2 + kvec(3)**2 )
      nn = nsphere*nn1
      call prinf('systems size is *',nn,1 )
      do i = 1,nsphere
         if (cdabs(zk0*sphererads(i)).lt.1.0d0) then
            scale(i) = cdabs(zk0*sphererads(i))
         else
            scale(i) = 1.0d0
         endif
      enddo
c
c     nsource is number of sources used in defining exterior field
c     for artificial test problem.
c
      nsource = nsphere
      source(1,0)=xxx(1)
      source(2,0)=xxx(2)
      source(3,0)=xxx(3)
      str(0) = 1.1d2
c
c     these sources i = 1, nsphere are only used in creating
c     artificial test problem.
c
      do i = 1,nsphere
         source(1,i)=centers(1,i) -0.2*sphererads(i)
         source(2,i)=centers(2,i) +0.1*sphererads(i)
         source(3,i)=centers(3,i) -0.12*sphererads(i)
         str(i) = dcmplx(1.2d0+i,0.5d0)
      enddo
c
c     beta,delta parameters set to one (for simplicity in this
c     demo code)
c
      beta0 = 1.0d0
      delta0 = 1.0d0
      do i = 1,nsphere
         beta(i) = 1.0d0
         delta(i) = 1.0d0
      enddo
c
c     get boundary data on spheres in Fourier (spherical
c     harmonic) domain.
c     (cw, iw are workspaces used to compute j expansions.)
c
c
      if (iscat.eq.0) then
c
c     create artificial data for known solution
c
         call createrhst(nterms,nsource,nsphere,centers,
     1          source,str,zk0,zk,beta0,beta,delta0,delta,
     1          sphererads,scale,rhs,rhs2,sh1,sh2,cw,lcw,iw,ier)
      else
c
c     create data corresponding to an incoming plane wave.
c
         call createrhspw(nterms,str(0),kvec,nsphere,centers,
     1          zk,beta0,beta,delta0,delta,
     2          sphererads,scale,rhs,rhs2,sh1,sh2,cw,lcw,iw,ier)
      endif
c
      call prin2('created right-hand side*',rhs,2*nn)
      call prin2('created right-hand side 2*',rhs2,2*nn)
c
c     solve the discretized system, obtaining the surface density
c
c     eps is GMRES tolerance 
c     numit is maximum number of iterations allowed
c     ngmrec is "restart parameter" for GMRES
c     nquad is number of points used for quadrature rules on sphere.
c           2*nterms is sufficient to yield extremely high accuracy.
c
      numit = maxit
      ngmrec = 50
      nquad = 2*nterms
      t0 = second()
      call spheresolver(nsphere,nterms,nquad,nn1,nn,centers,sphererads,
     1     scale,zk0,zk,beta0,delta0,beta,delta,
     1     multfmm3dt,rhs,rhs2,rhsjoin,eps,numit,sol,
     1     niter,errs,ngmrec,cw,www,lcw,ierr)
c
      t1 = second()
      call prin2('time for solve *',t1-t0,1)
      call prin2('finished iterative solution*',sol,0)
      call prin2(' errs is *',errs,niter)
ccc      call prin2('sol 1 is *',sol,2*nn1)
ccc      call prin2('sol 2 is *',sol(nn1+1),2*nn1)
c
c     check solution at exterior point by brute force:
c     i.e. evaluate contribution from each multipole expansion.
c
c     This is a useful check when iscat=0 where
c     createrhst constructs artificial data corresponding to 
c     known exact solution.
c     Ohterwise the error here makes no sense....
c 
      targ(1,1)=10
      targ(2,1)=5
      targ(3,1)=6
      iffld = 1
c
c     evaluate the solution of the Helmholtz equation at the
c     exterior test point directly
c
      call hpotfld3dall(iffld,source(1,1),str(1),nsource,
     1       targ,zk0,pot(1),fld2(1,1))
ccc      call prin2('exterior pt: direct is*',pot(1),2)
ccc      call prin2('diff  =*',htot-pot,2)
ccc      call prin2('ratio =*',htot/pot,2)
c
c    generate interior fields for each sphere and compute solution
c    at test point inside (for convenience, at the  location 
c    source(*,j)).
c
c     This is a useful check when using iscat = 0 where
c     createrhst constructs artificial data corresponding to 
c     known exact solution.
c     Otherwise the error here makes no sense....
c 
c     mpoletofield converts the solution vector sol to another 
c     sequence of expansion coefficients sh2 - the local (j)-expansion
c     coefficients on each of the spheres representing the interior
c     solution.
c
      call mpoletofield(sol,sh2,rhs,beta0,beta,nsphere,nterms,
     1     centers,sphererads,scale,zk0,zk,cw,lcw,ier)
      call prinf('after mpolefield ier =*',ier,1)
c
c     loop over all spheres and compute solution at test points.
c     targ(*,1) is in exterior.
c     targ(*,i+1) is in sphere i for i=1,nsphere.
c
      n1 = 1
      do i = 1,nsphere
         targ(1,i+1) = source(1,i)
         targ(2,i+1) = source(2,i)
         targ(3,i+1) = source(3,i)
c
c     compute field due to point source
c
         call hpotfld3dall(iffld,source(1,0),str(0),n1,
     1       source(1,i),zk(i),pot(i+1),fld2(1,i+1))
      enddo
      ntargs = nsphere+1
      call evalsolgrid(iffld,sol,zk0,zk,scale,nsphere,sphererads,
     1    centers,nterms,sh1,sh2,targ,ntargs,htot,
     2    fld,cw,lcw,lused,ier)
c
c     compute (incoming) plane wave at interior target
c
      call prinf('in exterior *',i,0)
      call prin2('computed sol is*',htot(1),2)
      call prin2('direct sol is *',pot(1),2)
      call prin2('computed fld is*',fld(1,1),6)
      call prin2('direct fld is *',fld2(1,1),6)
      call prin2('diff  =*',htot(1)-pot(1),2)
      call prin2('ratio =*',htot(1)/pot(1),2)
      do i = 2,nsphere+1
         call prinf('sphere is number *',i-1,1)
         call prin2('inside sphere: hout is*',htot(i),2)
         call prin2('inside sphere: direct is*',pot(i),2)
         call prin2('computed fld is*',fld(1,i),6)
         call prin2('direct fld is *',fld2(1,i),6)
         call prin2('diff  =*',htot(i)-pot(i),2)
         call prin2('ratio =*',htot(i)/pot(i),2)
      enddo
c
c     TO DO: Write out sol, sh2 to disk along with sphere data.
c     Build postprocessing routine that reads in this data,
c     gets target grid info, and evaluates as needed.
c
      hhx = (xhigh-xlow)/nnx
      hhy = (yhigh-ylow)/nny
      hhz = (zhigh-zlow)/nnz
      next = 0
      do k = 1,nnz
      do j = 1,nny
      do i = 1,nnx
         next = next+1
         targ(1,next) = xlow + i*hhx
         targ(2,next) = ylow + j*hhy
         targ(3,next) = zlow + k*hhz
      enddo
      enddo
      enddo
      write(6,*) 'next = ',next
      ntargs = nnx*nny*nnz
      call evalsolgrid(iffld,sol,zk0,zk,scale,nsphere,sphererads,
     1    centers,nterms,sh1,sh2,targ,ntargs,htot,
     2    fld,cw,lcw,lused,ier)
      open(unit=20, file='data.m',status='unknown') 
      open(unit=21, file='field.m',status='unknown') 
      write(20,*) 'edens = ['
      write(21,*) 'cfield = ['
      do i = 1,ntargs
         write(20,1001) targ(1,i), targ(2,i), targ(3,i),
     1               real(htot(i))
         write(21,1001) fld(1,i),fld(2,i),fld(3,i)
      enddo
      write(20,*) '];'
      write(21,*) '];'
      write(20,*) 'X = reshape(edens(:,1),',nnx,',',nny,',',nnz,');'
      write(20,*) 'Y = reshape(edens(:,2),',nnx,',',nny,',',nnz,');'
      write(20,*) 'Z = reshape(edens(:,3),',nnx,',',nny,',',nnz,');'
      write(20,*) 'U = reshape(edens(:,4),',nnx,',',nny,',',nnz,');'
      write(20,*) '[X,Y,Z] = meshgrid(-1.9:.1:8);'
      write(20,*) '[X,Y,Z] = meshgrid(-1.9:.1:8);'
      write(20,*) 'xslice = [-1.2,0.8,2];'
      write(20,*) 'yslice = [];'
      write(20,*) 'zslice = 0;'
      write(20,*) 'slice(X,Y,Z,U,xslice,yslice,zslice)'
c
c
c     call plotting routine that generates MATLAB output of sphere
c     geometry with field computed on surface of each sphere.
c     sphereplot.m  - spheres, sphereplot2.m  - plot of solution on
c     surface of spheres.
c
      call prin2('plotting solution *',sol,0)
      call plotspheres(nsphere,centers,sphererads,nterms,nn1,
     1             zk0,scale,sol)
1001  format(6D12.3)
      stop
      end
c
c
      subroutine createrhst(nterms,nsource,nsphere,centers,
     1          source,str,zk0,zk,beta0,beta,delta0,delta,
     2          sphererads,scale,rhs,rhs2,sh1,sh2,cw,lcw,iw,ier)
      implicit real *8 (a-h,o-z)
      integer *4 nsource,nsphere,nterms, iw(lcw), ier
      dimension sphererads(1)
      dimension scale(1)
      dimension source(3,0:1)
      dimension centers(3,1)
      complex *16 str(0:1)
      complex *16 z,zk0,eye
      complex *16 zk(1),beta0,beta(1),delta0,delta(1)
      complex *16 cw(lcw)
      complex *16 hout,houtx,houty,houtz
      complex *16 rhs(1)
      complex *16 rhs2(1)
      complex *16 sh1(0:nterms,-nterms:nterms)
      complex *16 sh2(0:nterms,-nterms:nterms)
      data eye/(0.0d0,1.0d0)/
c
c   
c     create artificial data corresponding to known
c                  exact solution. That solution is:
c
c                  phi (outside) is the field due to the 
c                  sources source(*,i),  i = 1,nsources (obviously 
c                  using the exterior Helmholtz parameter)
c
c                  phi (inside) is the field due to the single
c                  source source(*,0) (obviously using the Helmholtz 
c                  parameter specific to that sphere)
c
c
c     first get f and g (i.e. jumps induced by artificial
c     sources) f stored in rhs, g in rhs2.
c
      nn = (nterms+1)**2
      nsys = nsphere*nn
      do ii = 1,nsys
         rhs(ii) = 0d0
         rhs2(ii) = 0d0
      enddo 
      do nsp = 1,nsphere
         do ii = 1,nsource
            call sourcetoshxp(nterms,centers(1,nsp),source(1,ii),
     1          str(ii),zk0,sphererads(nsp),scale(nsp),sh1,sh2,
     2          cw,lcw,ier)
c
            next = (nsp-1)*nn+1
            do j = 0,nterms
            do k = -j,j
               rhs(next) = rhs(next) + beta0*sh1(j,k)
               rhs2(next) = rhs2(next) + delta0*sh2(j,k)
               next = next+1
            enddo
            enddo
         enddo
      enddo
      do nsp = 1,nsphere
         call sourcetoshxp(nterms,centers(1,nsp),source(1,0),
     1        str(0),zk(nsp),sphererads(nsp),scale(nsp),sh1,sh2,
     2        cw,lcw,ier)
c
         next = (nsp-1)*nn+1
         do j = 0,nterms
         do k = -j,j
            rhs(next)  = rhs(next) - beta(nsp)*sh1(j,k)
            rhs2(next) = rhs2(next) - delta(nsp)*sh2(j,k)
	    next = next+1
         enddo
         enddo
      enddo
      return
      end
c
c
c
      subroutine getspheredat(nsphere,
     1     sphererads,centers,kvec,zk)
      implicit none
      integer nsphere,i
      real *8 sphererads(nsphere),kvec(3)
      real *8 centers(3,nsphere)
      real *8 x1,x2,x3,x4,x5,x6,x7
      complex *16 zk(nsphere)
c
      open (unit = 4, file = 'spinsimp.dat')
      read(4,*) nsphere
      read(4,*) x1,x2,x3
      kvec(1) = x1
      kvec(2) = x2
      kvec(3) = x3
      do i=1,nsphere
         read(4,*) sphererads(i), centers(1,i), centers(2,i),
     1          centers(3,i), x1,x2
          zk(i) = dcmplx(x1,x2)
          call prin2('zk is *',zk(i),2)
          call prin2('centers is *',centers(1,i),3)
      enddo 
      return
      end
c
      subroutine createrhspw(nterms,str,kvec,nsphere,centers,
     1          zk,beta0,beta,delta0,delta,
     2          sphererads,scale,rhs,rhs2,sh1,sh2,cw,lcw,iw,ier)
      implicit real *8 (a-h,o-z)
      integer *4 nsource,nsphere,nterms, iw(lcw), ier
      dimension sphererads(1)
      dimension kvec(1)
      dimension scale(1)
      dimension centers(3,1)
      complex *16 str
      complex *16 z,eye
      complex *16 zk(1),beta0,beta(1),delta0,delta(1)
      complex *16 cw(lcw)
      complex *16 hout,houtx,houty,houtz
      complex *16 rhs(1)
      complex *16 rhs2(1)
      complex *16 sh1(0:nterms,-nterms:nterms)
      complex *16 sh2(0:nterms,-nterms:nterms)
      data eye/(0.0d0,1.0d0)/
c
c     creates right-hand side for scattering from an incoming plane wave
c     defined by kvec, with complex strength str. (The complex strength
c     allows for an arbitrary phase.
c
c     phi_in  = str * exp( i * ( kvec dot r) )
c     
c     get f and g induced by plane wave source.
c     f stored in rhs, g in rhs2.
c
      nn = (nterms+1)**2
      nsys = nsphere*nn
      do ii = 1,nsys
         rhs(ii) = 0d0
         rhs2(ii) = 0d0
      enddo 
      do nsp = 1,nsphere
         call h3dpwdata(ier,kvec,scale(nsp),centers(1,nsp),
     1          sphererads(nsp),nterms,sh1,sh2,cw,lcw,lused)
         next = (nsp-1)*nn+1
         do j = 0,nterms
         do k = -j,j
            rhs(next)  = -str*beta0*sh1(j,k)
            rhs2(next) = -str*delta0*sh2(j,k)
	    next = next+1
         enddo
         enddo
      enddo
      return
      end

