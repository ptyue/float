c***********************************************************************
      subroutine shap1d ( shg1,shg1x,we1 )
c
c     this routine calculates 1D shape functions and their derivatives
c 
c     OUTPUT:
c       shg1(nd,intg,nel) = shape functions
c       shg1x(nd,intg,nel) = kxi derivative of the shape functions
c       w(intg,nel) = weights for integration rules
c       where : 
c               nd = number number
c               intg = numerical integration points
c               nel = index for elemets
c                   = 1 : linear element
c                     2 : quadratic element
c***********************************************************************
      implicit none
      real*8 shg1(3,3,2),shg1x(3,3,2),we1(3,2)
c
      real*8 xi(3),we(3), x
      integer i
c
c     3 points integration rule
c     -------------------------
      xi(3) = dsqrt(0.6d0)
      xi(1) =-xi(3)
      xi(2) = 0.d0
      we(1) = 5.d0/9.d0
      we(3) = we(1)
      we(2) = 8.d0/9.d0
c
c     shape functions and their derivatives
c     -------------------------------------
      do i=1,3
        x = xi(i)
c
c       linear element
        shg1(1,i,1)  = 0.5d0*(1.d0-x)
        shg1(2,i,1)  = 0.5d0*(1.d0+x)
        shg1x(1,i,1) = -0.5d0
        shg1x(2,i,1) =  0.5d0
        we1(i,1) = we(i)
c
c       quadratic element
        shg1(1,i,2)  = -0.5d0*x*(1.d0-x)
        shg1(2,i,2)  =  1.d0-x*x
        shg1(3,i,2)  =  0.5d0*x*(1.d0+x)
        shg1x(1,i,2) =  x-0.5d0
        shg1x(2,i,2) = -2.d0*x
        shg1x(3,i,2) =  x+0.5d0
        we1(i,2) = we(i)
c
      enddo
c
      return
      end
