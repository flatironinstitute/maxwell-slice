
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       This is the end of the debugging code and the beginning of the
c       evaluation code for the Legendre functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c       This file contains a set of subroutines for the handling of
c       Legendre functions. It contains 10 subroutines that
c       are user-callable. Following is a brief description of these
c       subroutines.
c
c
c
c       ylgndr - evaluate normalized Legendre functions
c
c       ylgndr2 - evaluate normalized Legendre functions and their derivatives
c
c       ylgndr2s - evaluate normalized Legendre functions and their
c            derivatives. Ynm(x), m>0 values are scaled by 1/sqrt(1-x^2),
c            Ynm(x), m>0 derivatives with respect to are scaled by
c            sqrt(1-x^2).
c
c       ylgndrini - precompute the recursion coefficients.
c
c       ylgndrf - evaluate normalized Legendre functions. Fast version,
c            the recursion coefficients are precomputed with ylgndrini.
c
c       ylgndr2f - evaluate normalized Legendre functions and their
c            derivatives. Fast version, the recursion coefficients are
c            precomputed with ylgndrini.
c
c       ylgndr2sf - evaluate normalized Legendre functions and their
c            derivatives. Ynm(x), m>0 values are scaled by 1/sqrt(1-x^2),
c            Ynm(x), m>0 derivatives with respect to are scaled by
c            sqrt(1-x^2).  Fast version, the recursion coefficients are
c            precomputed with ylgndrini.
c
c       zylgndr - evaluate normalized Legendre functions of a complex argument
c
c       zylgndr2 - evaluate normalized Legendre functions of a complex argument
c            and their derivatives
c
c       zylgndrsc - evaluate normalized Legendre functions of a complex
c            argument. Scaled version to prevent overflows for large z values.
c
c       Added 03/08/09
c
c       ylgndr2s_trunc - evaluate normalized Legendre functions and their
c            derivatives (with scaling). Same as ylgndr2s, 
c            but recursion carried out only to m = m2 rather than nmax.
c
c       ylgndrf_trunc - evaluate normalized Legendre functions. Same as ylgndrf 
c            but recursion carried out only to m = m2 rather than nmax.
c
c       ylgndr2f_trunc - evaluate normalized Legendre functions and their
c            derivatives. Same as ylgndr2f, 
c            but recursion carried out only to m = m2 rather than nmax.
c
c       ylgndr2sf_trunc - evaluate normalized Legendre functions and their
c            derivatives (with scaling). Same as ylgndr2sf,
c            but recursion carried out only to m = m2 rather than nmax.
c
c
      subroutine ylgndr(nmax, x, y)
      implicit real *8 (a-h,o-z)
c
c     Evaluate normalized Legendre function 
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x)
c     for 0 <= n <= nmax  and  0 <= m <= n.  Other elements of y
c     will contain undefined values.
c
      integer nmax, m
      real *8 x, y(0:nmax,0:nmax), u, v, w
      u=-sqrt(1-x*x)
      y(0,0)=1
      do 10 m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 do 20 n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
 20      continue
 10   continue
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndr2(nmax, x, y, d)
      implicit real *8 (a-h,o-z)
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
      integer nmax, m
      real *8 x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du, v, w
      u=-sqrt(1-x*x)
      du=x/sqrt(1-x*x)
      y(0,0)=1
      d(0,0)=0
      do 10 m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*x)/u**2
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
	 do 20 n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
 20      continue
 10   continue
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndr2s(nmax, x, y, d)
      implicit real *8 (a-h,o-z)
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
      integer nmax, m
      real *8 x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, v, w
      u=-sqrt(1-x*x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      do 120 n=m+2, nmax
        y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
120   continue
c
c       ... then, evaluate scaled associated Legendre functions
c
      do 210 m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*x)
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+(1-x**2)*y(m,m))*sqrt(2*m+1.0d0)
	 do 220 n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+(1-x**2)*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
220      continue
210   continue
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       fast version of real valued Legendre functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine ylgndrini(nmax, rat1, rat2)
      implicit real *8 (a-h,o-z)
c
c     Precompute the recurrence coefficients for the fast
c     evaluation of normalized Legendre functions and their derivatives
c    
c     Parameters:
c       nmax                      must be non-negative
c       rat1(0:nmax,0:nmax)       recurence coefficient
c       rat2(0:nmax,0:nmax)       recurence coefficient
c
      integer nmax, m
      real *8 x, rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      rat1(0,0)=1
      rat2(0,0)=1
      do 10 m=0, nmax
	 if (m.gt.0)  rat1(m,m)=sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  rat2(m,m)=1
	 if (m.lt.nmax)  rat1(m+1,m)=sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  rat2(m+1,m)=1
	 do 20 n=m+2, nmax
	    rat1(n,m)=(2*n-1)
            rat2(n,m)=sqrt((n+m-1.0d0)*(n-m-1.0d0))
	    rat1(n,m)=rat1(n,m)/sqrt(dble(n-m)*(n+m))
	    rat2(n,m)=rat2(n,m)/sqrt(dble(n-m)*(n+m))
 20      continue
 10   continue
c
      return
      end
c
c
c
c
c
      subroutine ylgndrf(nmax, x, y, rat1, rat2)
      implicit real *8 (a-h,o-z)
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
      integer nmax, m
      real *8 x, y(0:nmax,0:nmax), u, v, w
      real *8 rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-x*x)
      y(0,0)=1
      do 10 m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 do 20 n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
 20      continue
 10   continue
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
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
      subroutine ylgndr2f(nmax, x, y, d, rat1, rat2)
      implicit real *8 (a-h,o-z)
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
      integer nmax, m
      real *8 x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du, v, w
      real *8 rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-x*x)
      du=x/sqrt(1-x*x)
      y(0,0)=1
      d(0,0)=0
      do 10 m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*x)/u**2
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
	 do 20 n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
 20      continue
 10   continue
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
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
      subroutine ylgndr2sf(nmax, x, y, d, rat1, rat2)
      implicit real *8 (a-h,o-z)
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
      integer nmax, m
      real *8 x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, v, w
      real *8 rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-x*x)
      u2 = 1-x*x
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      do 120 n=m+2, nmax
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
120   continue
c
c       ... then, evaluate scaled associated Legendre functions
c
      do 210 m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*x)
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  
ccc     $      d(m+1,m)=(x*d(m,m)+(1-x**2)*y(m,m))*rat1(m+1,m)
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 do 220 n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
ccc	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+(1-x**2)*y(n-1,m))-
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
220      continue
210   continue
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       complex valued Legendre functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine zylgndr(nmax, z, y)
      implicit none
c
c     Evaluate normalized Legendre function for complex argument
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     z                         complex*16
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(z)
c     for 0 <= n <= nmax  and  0 <= m <= n.  Other elements of y
c     will contain undefined values.
c
      integer nmax
      complex*16 z, y(0:nmax,0:nmax)
c
      integer m,n
      complex*16 u, v, w
c
      u=-sqrt(1-z*z)
      y(0,0)=1
      do 10 m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*sqrt(2*m+1.0d0)
	 do 20 n=m+2, nmax
	    y(n,m)=((2*n-1)*z*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
 20      continue
 10   continue
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine zylgndr2(nmax, z, y, d)
      implicit real *8 (a-h,o-z)
c
c     Evaluate normalized Legendre functions and their derivatives
c       for complex argument
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values, complex *16
c     d(0:nmax,0:nmax)      resulting derivative values, complex *16
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention is valid for the derivative.
c
      integer nmax, m
      complex *16 z, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du, v, w
      u=-sqrt(1-z*z)
      du=z/sqrt(1-z*z)
      y(0,0)=1
      d(0,0)=0
      do 10 m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*z)/u**2
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  d(m+1,m)=(z*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
	 do 20 n=m+2, nmax
	    y(n,m)=((2*n-1)*z*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(z*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
 20      continue
 10   continue
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine zylgndrsc(nmax, z,scale, ysc)
      implicit none
c
c     Evaluate scaled versions of zylgndr.  zylgndr is the complex version
c     of the normalized Legendre function.  Scaling removes the possible
c     conditioning errors from zylgndr evaluated at large arguments.
c
c     y_nm = scale**(-n) * ysc_nm
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     z (complex*16)            the argument
c     scale (real*8)
c          looks like scale wants to be O(1/z) and scale<1.
c     ysc(0:nmax,0:nmax)        resulting function values
c
      integer nmax
      complex*16 z, ysc(0:nmax,0:nmax)
      real*8 scale
c
      integer m,n
      complex*16 u, v, w
c
      u=-sqrt(1-z*z)
      ysc(0,0)=1
      do 10 m=0, nmax
	 if (m.gt.0)  then
            ysc(m,m)=ysc(m-1,m-1)*scale*u*sqrt((2*m-1.0d0)/(2*m))
c           call prinf('m=*',m,1)
c           call prin2('ysc(m,m)=*',ysc(m,m),2)
         endif
	 if (m.lt.nmax)  then
            ysc(m+1,m)=z*scale*ysc(m,m)*sqrt(2*m+1.0d0)
         endif
	 do 20 n=m+2, nmax
	    ysc(n,m)=((2*n-1)*scale*z*ysc(n-1,m) - 
     1           sqrt((n+m-1.0d0)*(n-m-1.0d0))*scale**2*ysc(n-2,m))
     2           /sqrt((n-m+0.0d0)*(n+m))
 20      continue
 10   continue
      do n=0, nmax
	 do m=0, n
	    ysc(n,m)=ysc(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       truncated recurrences for Legendre functions:
c    
c       Ynm(theta) for n = 0,nmax  but m = -m2,...,m2  with m2 < nmax.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
      subroutine ylgndr2s_trunc(nmax, m2, x, y, d)
      implicit real *8 (a-h,o-z)
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
      integer nmax, m, m2
      real *8 x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, v, w
      u=-sqrt(1-x*x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      do 120 n=m+2, nmax
        y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
120   continue
c
c       ... then, evaluate scaled associated Legendre functions
c
      do 210 m=1, m2
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*x)
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+(1-x**2)*y(m,m))*sqrt(2*m+1.0d0)
	 do 220 n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+(1-x**2)*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
220      continue
210   continue
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
c
      subroutine ylgndrf_trunc(nmax, m2, x, y, rat1, rat2)
      implicit real *8 (a-h,o-z)
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., min(m2,n).
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
      integer nmax, m, m2
      real *8 x, y(0:nmax,0:nmax), u, v, w
      real *8 rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-x*x)
      y(0,0)=1
      do 10 m=0, m2
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 do 20 n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
 20      continue
 10   continue
c
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
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
      subroutine ylgndr2f_trunc(nmax, m2, x, y, d, rat1, rat2)
      implicit real *8 (a-h,o-z)
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., min(m2,n).
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
      integer nmax, m, m2
      real *8 x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du, v, w
      real *8 rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-x*x)
      du=x/sqrt(1-x*x)
      y(0,0)=1
      d(0,0)=0
      do 10 m=0, m2
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*x)/u**2
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
	 do 20 n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
 20      continue
 10   continue
c
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
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
      subroutine ylgndr2sf_trunc(nmax, m2, x, y, d, rat1, rat2)
      implicit real *8 (a-h,o-z)
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., min(n,m2).
c
c     Parameters:
c     nmax                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
      integer nmax, m, m2
      real *8 x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, v, w
      real *8 rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-x*x)
      u2 = 1-x*x
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      do 120 n=m+2, nmax
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
120   continue
c
c       ... then, evaluate scaled associated Legendre functions
c
      do 210 m=1, m2
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m*x)
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 do 220 n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
220      continue
210   continue
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
