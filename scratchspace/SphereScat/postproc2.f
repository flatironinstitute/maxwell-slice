C  Copyright (C) 2007: Leslie Greengard
C  Contact: greengard@cims.nyu.edu
C
C  All rights reserved.
C  This copy is not for distribution without permission of the author.
C
C**********************************************************************C
      subroutine evalextscat(iffld,sol,zk0,scale,nspheres,x0y0z0,
     1           nterms,sh1,targ,htot,fld,cw,lcw,lused,ier)
c***********************************************************************
c
c     evaluate solution at exterior point by direct evaluation
c     of all multipole expansions
c
c     INPUT: 
c
c     iffld = 0 -> only potential computed, 1 -> field also computed.
c     sol = solution vector from scattering integral equation.
c     zk0 = Helmholtz parameter
c     scale = scaling vector for multipole expansions (see Solver)
c     nspheres = number of spheres 
c     x0y0z0 = vector of sphere center locations
c     nterms = length of multipole expansions
c     sh1    = workspace of length at least (nterms+1)*(2*nterms+1)
c     targ = target location
c     cw =  workspace
c     lcw =  workspace
c  
c     OUTPUT: 
c
c     hout = potential at target
c     fld =  negative gradient of potential at target
c     lused = workspace used
c     ier = error flag, 0 -> no error, 1 -> insufficient workspace
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      dimension scale(nspheres)
      dimension x0y0z0(3,nspheres)
      dimension targ(3)
      integer *4 lcw,i
      complex *16 sol(*),fld(3),fldloc(3)
      complex *16 zk0
      complex *16 sh1(*)
      complex *16 cw(lcw),eye
      complex *16 hout,htot
      data eye/(0.0d0,1.0d0)/   
c
      nn1 = (nterms+1)**2
      htot = 0.0d0
      fld(1) = 0.0d0
      fld(2) = 0.0d0
      fld(3) = 0.0d0
      do i = 1,nspheres
         call unrolltoshxp(nterms,sol((i-1)*nn1+1),sh1)
         call h3dmpeval(zk0,scale(i),x0y0z0(1,i),sh1,
     1       nterms,targ,hout,iffld,fldloc,cw,lcw,lused,ier)
         htot = htot + hout
         if (iffld.eq.1) then 
            fld(1) = fld(1) + fldloc(1)
            fld(2) = fld(2) + fldloc(2)
            fld(3) = fld(3) + fldloc(3)
         endif 
      enddo
c
      return
      end
c
c***********************************************************************
      subroutine evalintscat(iffld,sh2,zk,scale,x0y0z0,
     1           nterms,sh1,targ,hout,fld,cw,lcw,lused,ier)
c***********************************************************************
c
c     evaluate solution at interior point
c
c     INPUT: 
c
c     iffld = 0 -> only potential computed, 1 -> field also computed.
c     sh2 = local expansion in unrolled (compressed) format.
c     zk = Helmholtz parameter
c     scale = scaling vector for multipole expansions (see Solver)
c     x0y0z0 = vector of sphere center locations
c     nterms = length of multipole expansions
c     sh1    = workspace of length at least (nterms+1)*(2*nterms+1)
c     targ = target location
c     cw =  workspace
c     lcw =  workspace
c  
c     OUTPUT: 
c
c     hout = potential at target
c     fld =  negative gradient of potential at target
c     lused = workspace used
c     ier = error flag, 0 -> no error, 1 -> insufficient workspace
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      dimension x0y0z0(3)
      dimension targ(3)
      integer *4 lcw,i
      real *8 scale
      complex *16 fld(3)
      complex *16 zk
      complex *16 sh1(*), sh2(*)
      complex *16 cw(lcw),eye
      complex *16 hout
      data eye/(0.0d0,1.0d0)/   
c
      call unrolltoshxp(nterms,sh2,sh1)
      call h3dtaeval(zk,scale,x0y0z0(1),sh1,
     1       nterms,targ,hout,iffld,fld,cw,lcw,lused,ier)
c
      return
      end
c
      subroutine plotspheres(nspheres,xcs,sphererads,nterms,nn1,
     1             zk0,scale,sol)
      implicit real *8 (a-h,o-z)
      integer *4 nspheres, nn1, lcw, iffld
      parameter(lcw=100000)
      real *8 xcs(3,nspheres)
      real *8 sphererads(nspheres)
      real *8 targ(3), scale(1)
      real *8 x1(100,100)
      real *8 y1(100,100)
      real *8 z1(100,100)
      real *8 cc(100,100)
      real *8 theta,phi
      complex *16 zk0,sh1(10 000), sol(1), cw(lcw)
      complex *16 htot, hout,fld(3)
c
      open (unit = 3, file = 'sphereplot.m')
      open (unit = 4, file = 'sphereplot2.m')
      pi = 4*datan(1.0d0)
      write(3,*)'[X,Y,Z] = sphere(20)'
      do i = 1,nspheres
        write(3,*)' X1=',sphererads(i),'*X + ',xcs(1,i),'*ones(21,21);'
        write(3,*)' Y1=',sphererads(i),'*Y + ',xcs(2,i),'*ones(21,21);'
        write(3,*)' Z1=',sphererads(i),'*Z + ',xcs(3,i),'*ones(21,21);'
        write(3,*)' surf(X1,Y1,Z1)'
        if (i.eq.1) write(3,*)' hold on'
        if (i.eq.1) write(3,*)' axis equal'
      enddo
      write(3,*)' axis image'
      do i = 1,nspheres
         npts = 40
         do j = 1,npts
            theta = pi*(j-1)/(npts-1)
            do k = 1,npts
               phi = 2*pi*(k-1)/(npts-1)
               x1(j,k) = xcs(1,i) + 
     1            sphererads(i)*sin(theta)*cos(phi)
               y1(j,k) = xcs(2,i) + 
     2            sphererads(i)*sin(theta)*sin(phi)
               z1(j,k) = xcs(3,i) + sphererads(i)*cos(theta)
               targ(1) = x1(j,k)
               targ(2) = y1(j,k)
               targ(3) = z1(j,k)
               htot = 0.0d0
               do ii = 1,nspheres
                  call unrolltoshxp(nterms,sol((ii-1)*nn1+1),sh1)
                  iffld = 0
                  call h3dmpeval(zk0,scale(ii),xcs(1,ii),sh1,
     1               nterms,targ,hout,iffld,fld,cw,lcw,lused,ier)
                  htot = htot + hout
               enddo
               cc(j,k) = dreal(htot)
            enddo
         enddo
         write(4,*)' XX1=['
         do j = 1,npts
            write(4,*)(x1(j,k),k=1,npts)
         enddo
         write(4,*)' ];'
         write(4,*)' YY1=['
         do j = 1,npts
            write(4,*)(y1(j,k),k=1,npts)
         enddo
         write(4,*)' ];'
         write(4,*)' ZZ1=['
         do j = 1,npts
            write(4,*)(z1(j,k),k=1,npts)
         enddo
         write(4,*)' ];'
         write(4,*)' CC=['
         do j = 1,npts
            write(4,*)(cc(j,k),k=1,npts)
         enddo
         write(4,*)' ];'
         write(4,*)' surf(XX1,YY1,ZZ1,CC)'
         if (i.eq.1) write(4,*)' hold on'
         if (i.eq.1) write(4,*)' axis equal'
      enddo
      write(4,*)' axis image'
      write(4,*)' figure(1)'
      return
      end
c
C**********************************************************************C
      subroutine evalsolgrid(iffld,sol,zk0,zk,scale,nspheres,rads,
     1    x0y0z0,nterms,sh1,sh2,targ,ntargs,htot,fld,cw,lcw,lused,ier)
c***********************************************************************
c
c     evaluate solution at exterior point by direct evaluation
c     of all multipole expansions
c
c     INPUT: 
c
c     iffld = 0 -> only potential computed, 1 -> field also computed.
c     sol = solution vector from scattering integral equation.
c     zk0 = exterior Helmholtz parameter
c     zk  = interior Helmholtz parameters
c     scale = scaling vector for multipole expansions (see Solver)
c     nspheres = number of spheres 
c     rads   = sphere radii 
c     x0y0z0 = vector of sphere center locations
c     nterms = length of multipole expansions
c     sh1    = workspace of length at least (nterms+1)*(2*nterms+1)
c     sh2    = unrolled solution
c     targ = target locations
c     ntarg = number of targets
c     cw =  workspace
c     lcw =  workspace
c  
c     OUTPUT: 
c
c     htot = potential at targets
c     fld =  negative gradient of potential at targets
c     lused = workspace used
c     ier = error flag, 0 -> no error, 1 -> insufficient workspace
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 scale(nspheres)
      real *8 rads(nspheres)
      real *8 x0y0z0(3,nspheres)
      real *8 targ(3,ntargs)
      integer *4 lcw,i
      complex *16 sol(*)
      complex *16 zk0
      complex *16 zk(nspheres)
      complex *16 sh1(*)
      complex *16 sh2(*)
      complex *16 cw(lcw),eye
      complex *16 hout,htot(ntargs),fld(3,ntargs)
      data eye/(0.0d0,1.0d0)/   
c
      nn1 = (nterms+1)**2
      do itarg = 1,ntargs
ccc         write(17,*) ' targ ',targ(1,itarg),targ(2,itarg),targ(3,itarg)
         do i = 1,nspheres
            r1 = targ(1,itarg) - x0y0z0(1,i)
            r2 = targ(2,itarg) - x0y0z0(2,i)
            r3 = targ(3,itarg) - x0y0z0(3,i)
            rr = dsqrt(r1*r1+r2*r2+r3*r3)
            if (rr.lt. rads(i)) then
ccc              write(17,*) ' inside sphere ',i
              call evalintscat(iffld,sh2((i-1)*nn1+1),zk(i),scale(i),
     1           x0y0z0(1,i),nterms,sh1,targ(1,itarg),htot(itarg),
     2           fld(1,itarg),cw,lcw,lused,ier)
ccc              write(17,*) ' htot is ',htot(itarg)
              goto 100
            endif
         enddo
         call evalextscat(iffld,sol,zk0,scale,nspheres,x0y0z0,
     1           nterms,sh1,targ(1,itarg),htot(itarg),
     2           fld(1,itarg),cw,lcw,lused,ier)

100      continue
      enddo
      return
      end
