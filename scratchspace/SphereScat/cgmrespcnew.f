ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine cgmrespc(ier,n,mult,y,eps,maxiter,
     1     x,niter,errs,nrec,w,wfmm,lenwfmm,par1,par2)
        implicit real *8 (a-h,o-z)
        dimension errs(1)
        complex *16 x(1),y(1),w(1)
        complex *16 wfmm(1)
c
        external mult
c
c     This subroutine solves a complex linear system Ax=y by means
c     of GMRES algorithm. This is a memory management routine for 
c     cgmres1 which performs the actual work.
c
c     Input parameters:
c
c     n - the dimensionality of the linear system
c     y - the right hand side
c     a - the matrix of the system (or whatever other parameter)
c     eps - the required accuracy
c     mult - the user-defined matrix-vector multiplication subroutine
c
c     the calling sequence for mult must be
c
c     mult(x,y,n,wfmm,lenwfmm,par1,par2)               (1)
c
c     in (1), a is a matrix the system or whatever other parameter, x is
c     an input vector, y is a product Ax, and n is the dimensionality of
c     a, x, and y.
c
c     maxiter - the maximum number of iteration permitted
c     nrec - the maximum number of iteration 
c            which GMRES algorithm needs to be restarted
c
c     w - must be at least (nrec*2+4)*n complex *16 elements long
c
c     Output parameters:
c
c     ier - error return code
c        ier=0 normal execution of the subroutine
c        ier=4 means that the maximum number iterations maxiter
c           has been reached without achieving the required accuracy eps
c        ier=8 means that the errors failed to decrease before the maximum
c           number of iterations has been reached or the required accuracy
c           eps has benn reached. 
c
c     x - the solution of the system
c     niter - the number of iterations performed 
c     errs - the array of errors produced by the algorithm. 
c        errs(i)=||y-Ax_i||, where x_i is a solution obtained on the i-th
c        GMRES iteration.
c
ccc	call prinf(' in cgmrespc lenwfmm is *',lenwfmm,1)
        ie=1
        lie=n*nrec
c
        iae=ie+lie
        liae=n*nrec
c
        iz=iae+liae
        lz=n
c
        iw=iz+lz
        lw=n
c
        ixr=iw+lw
        lxr=n
c
        izr=ixr+lxr
        lzr=n
c
        call cgmres1(ier,n,mult,y,eps,maxiter,x,niter,errs,
     1     nrec,w(ie),w(iae),w(iz),w(iw),w(ixr),w(izr),wfmm,lenwfmm,
     2     par1,par2)

        return
        end
c
c
c
c
c
        subroutine cgmres1(ier,n,mult,y,eps,maxiter,
     1     x,niter,errs,nrec,e,ae,z,w,xr,zr,wfmm,lenwfmm,
     2     par1,par2)
        implicit real *8 (a-h,o-z)
c
        dimension errs(1)
        complex *16 x(1),y(1),e(n,1),ae(n,1),
     1       z(1),w(1),xr(1),zr(1),d
        complex *16 wfmm(1),par1,par2
c
        external mult
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        ier=0
        niter=0
c
        do i=1,n
ccc           x(i)=0
c
c     set x=y for cgmrespp  
c
           x(i)=y(i)
        enddo
c
        dold=1d+100
c
        m=0
c
        do 4000 iter=1,maxiter
c     
c     ... restart GMRES
c
        if( m .eq. nrec ) then
           m=0
        endif
c
        m=m+1
        niter=niter+1
c
cccc        call prinf('niter=*',niter,1)
ccc        call prinf('m=*',m,1)
c
        if( m .ne. 1 ) goto 1600
c
c     ... this is the first iteration, calculate the residual z=y-Ax
c
        call mult(x,z,n,wfmm,lenwfmm,par1,par2)
c
        do i=1,n
           z(i)=y(i)-z(i)
        enddo
c
c       store initial residual norm as rnorm0
c
        rnorm0=0
        do i=1,n
           rnorm0=rnorm0+dconjg(z(i))*z(i)
        enddo
        rnorm0=sqrt(rnorm0)
c
c     ... store the residual z for future use
c
        do i=1,n
           zr(i)=z(i)
        enddo
c
 1600 continue
c
c     ... retrieve the residual z
c
        do i=1,n
           z(i)=zr(i)
        enddo
c
c     ... compute the error
c
        d=0
        do i=1,n
           d=d+dconjg(z(i))*z(i)
        enddo
        d=sqrt(d)
ccc        errs(niter)=dreal(d)
        errs(niter)=dreal(d)/rnorm0
c
        if( abs(d) .lt. eps ) return
c
        if( abs(d) .gt. dold ) then
c
c     ... the errors stopped decreasing, abort
c
           niter=niter-1
           do i=1,n
              x(i)=xr(i)
           enddo
           ier=8
           return
        endif
        dold=d
c
c     ... compute the new direction w=Az
c
        call mult(z,w,n,wfmm,lenwfmm,par1,par2)
c
c     ... orthogonalize w to all preceeding ae
c
        do j=1,m-1
        d=0
        do i=1,n
           d=d+dconjg(ae(i,j))*w(i)
        enddo
        do i=1,n
           w(i)=w(i)-d*ae(i,j)
           z(i)=z(i)-d*e(i,j)
        enddo
        enddo
c
c     ... normalize the current direction w
c
        d=0
        do i=1,n
           d=d+dconjg(w(i))*w(i)
        enddo
        d=1/sqrt(d)
c
c     ... store e and ae 
c
        do i=1,n
           ae(i,m)=w(i)*d
           e(i,m)=z(i)*d
        enddo
c
c     ... update the solution x and the residual z
c
        d=0
        do i=1,n
           d=d+dconjg(ae(i,m))*zr(i)
        enddo
        do i=1,n
           zr(i)=zr(i)-d*ae(i,m)
        enddo
        do i=1,n
           xr(i)=x(i)
           x(i)=x(i)+d*e(i,m)
        enddo
c
 4000   continue
c
c     ... the maximum number of iterations has been reached, abort
c
        ier=4
c
        return
        end
c
c
c
c
c
        subroutine cgmscapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        complex *16 x(1),y(1),prod
c
        prod=0
        do 1200 i=1,n
           prod=prod+dconjg(x(i))*y(i)
 1200   continue
        return
        end
c
c
c
c
c
        subroutine cgmrespp(w,nrec,niter,n,b,a,aw)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,1), b(n,n), a(n,n), cd, aw(1), alpha
c
c     This subroutine finds the inverse (or an approximation of the
c     inverse) of the matrix a by using the rank-1 updates contructed
c     via a preceding call to cgmres.
c
c     This subroutine assumes that we will start with the unit matrix,
c     which is to say that, in cgmres, the initial approximation of the
c     solution x MUST be set to the right hand side y.
c
c     input parameters:
c
c     a - the matrix to be inverted
c     w - the rank-1 update vectors from cgmres routine
c     n - the dimensionality of the problem
c     nrec - the maximum number of iteration 
c            before GMRES algorithm needs to be restarted
c     niter - the actual number iteration performed 
c     
c     aw - work array. must be at least n complex *16 elements long
c
c     output parameters:
c
c     b - the inverse of a
c     

        alpha=1
c        
        do i=1,n
           do j=1,n
              b(i,j)=0
           enddo
           b(i,i)=alpha
        enddo
c
        do k=1,niter
c
           call cgmmult(a,w(1,k+nrec),aw,n)
c
           do i=1,n
              do j=1,n
                 b(i,j)=b(i,j)+dconjg(w(i,k))
     $                *(w(j,k+nrec)-aw(j)*alpha) 
              enddo
           enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine cgmmult(a,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(1),y(1),d
c
        do 1400 i=1,n
           d=0
           do 1200 j=1,n
              d=d+dconjg(a(j,i))*x(j)
 1200      continue
           y(i)=d
 1400   continue
c     
        return
        end
c     
c
c
c
c
