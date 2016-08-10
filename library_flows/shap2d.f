c***********************************************************************
      subroutine shap2d ( shg2,shg2x,shg2y,we2 )
c
c     this routine calculates 2D shape functions and their derivatives
c 
c     OUTPUT:
c       shg2(nd,intg,nel) = shape functions
c       shg2x(nd,intg,nel) = kxi derivative of the shape functions
c       shg2y(nd,intg,nel) = eta derivative of the shape functions
c       we2(intg,nel) = weights for integration rules
c       where : 
c               nd = number number
c               intg = numerical integration points
c               nel = index for elemets
c                   = 1 : P1 (linear 3-node triangular element)
c                     2 : P2 (quadratic 6-node triangular element)
c                     3 : Q1 (bi-linear 4-node quadrilateral element)
c                     4 : Q2 (bi-quadratic 9-node quadrilateral element)
c                     5 : 8-node serendipity element
c                     6 : P0 (piece-wise constant triangular element)
c***********************************************************************
      implicit none
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
c
      real*8 xi(9),et(9),we(9), x,y,z,a,b,c
      integer i
c
c     7 points integration rule for triangles : nodes and weights
c     -----------------------------------------------------------
      xi(1) = 1.d0/3.d0
      et(1) = xi(1)
      we(1) = 0.112515000150000d0
      xi(2) = 0.797426985353087d0
      et(2) = 0.101286507323456d0
      we(2) = 0.0629695902724135d0
      xi(3) = et(2)
      et(3) = xi(2)
      we(3) = we(2)
      xi(4) = xi(3)
      et(4) = et(2)
      we(4) = we(2)
      xi(5) = 0.470142064105115d0
      et(5) = 0.059715871789770d0
      we(5) = 0.066197076394253d0
      xi(6) = xi(5)
      et(6) = xi(5)
      we(6) = we(5)
      xi(7) = et(5)
      et(7) = xi(5)
      we(7) = we(5)
c
c     shape functions and their derivatives
c     -------------------------------------
      do i=1,7
        x = xi(i)
        y = et(i)
        z = 1.d0 - x - y
c
c       constant element
        shg2(1,i,6)  = 1.d0
        shg2x(1,i,6) = 0.d0
        shg2y(1,i,6) = 0.d0
        we2(i,6)     = we(i)
c
c       linear element
        shg2(1,i,1)  =  z
        shg2(2,i,1)  =  x
        shg2(3,i,1)  =  y
        shg2x(1,i,1) = -1.d0
        shg2x(2,i,1) =  1.d0
        shg2x(3,i,1) =  0.d0
        shg2y(1,i,1) = -1.d0
        shg2y(2,i,1) =  0.d0
        shg2y(3,i,1) =  1.d0
        we2(i,1)     =  we(i)
c
c       quadratic element
        shg2(1,i,2)  = z*(2.d0*z-1.d0)
        shg2(2,i,2)  = x*(2.d0*x-1.d0)
        shg2(3,i,2)  = y*(2.d0*y-1.d0)
        shg2(4,i,2)  = 4.d0*z*x
        shg2(5,i,2)  = 4.d0*x*y
        shg2(6,i,2)  = 4.d0*y*z
        shg2x(1,i,2) =-4.d0*z+1.d0
        shg2x(2,i,2) = 4.d0*x-1.d0
        shg2x(3,i,2) = 0.d0
        shg2x(4,i,2) = 4.d0*(z-x)
        shg2x(5,i,2) = 4.d0*y
        shg2x(6,i,2) =-4.d0*y
        shg2y(1,i,2) =-4.d0*z+1.d0
        shg2y(2,i,2) = 0.d0
        shg2y(3,i,2) = 4.d0*y-1.d0
        shg2y(4,i,2) =-4.d0*x
        shg2y(5,i,2) = 4.d0*x
        shg2y(6,i,2) = 4.d0*(z-y)
        we2(i,2) = we(i)
c
      enddo
c
c     9 points integration rule for quadrilaterals : nodes and weights
c     ----------------------------------------------------------------
      a = dsqrt(0.6d0)
      b = 5.d0/9.d0
      c = 8.d0/9.d0
      xi(1) = a
      et(1) = a
      we(1) = b*b
      xi(2) = 0.d0
      et(2) = a
      we(2) = b*c
      xi(3) =-a
      et(3) = a
      we(3) = b*b
      xi(4) = a
      et(4) = 0.d0
      we(4) = b*c
      xi(5) = 0.d0
      et(5) = 0.d0
      we(5) = c*c
      xi(6) =-a
      et(6) = 0.d0
      we(6) = b*c
      xi(7) = a
      et(7) =-a
      we(7) = b*b
      xi(8) = 0.d0
      et(8) =-a
      we(8) = b*c
      xi(9) =-a
      et(9) =-a
      we(9) = b*b
c
c     shape functions and their derivatives
c     -------------------------------------
      do i=1,9
        x = xi(i)
        y = et(i)
c
c       linear element
        shg2(1,i,3)  = (1.d0+x)*(1.d0+y)*0.25d0
        shg2(2,i,3)  = (1.d0-x)*(1.d0+y)*0.25d0
        shg2(3,i,3)  = (1.d0-x)*(1.d0-y)*0.25d0
        shg2(4,i,3)  = (1.d0+x)*(1.d0-y)*0.25d0
        shg2x(1,i,3) = (1.d0+y)*0.25d0
        shg2x(2,i,3) =-(1.d0+y)*0.25d0
        shg2x(3,i,3) =-(1.d0-y)*0.25d0
        shg2x(4,i,3) = (1.d0-y)*0.25d0
        shg2y(1,i,3) = (1.d0+x)*0.25d0
        shg2y(2,i,3) = (1.d0-x)*0.25d0
        shg2y(3,i,3) =-(1.d0-x)*0.25d0
        shg2y(4,i,3) =-(1.d0+x)*0.25d0
        we2(i,3) = we(i)
c
c       quadratic element
        shg2(1,i,4)  =  x*(1.d0+x)*y*(1.d0+y)*0.25d0
        shg2(2,i,4)  = -x*(1.d0-x)*y*(1.d0+y)*0.25d0
        shg2(3,i,4)  =  x*(1.d0-x)*y*(1.d0-y)*0.25d0
        shg2(4,i,4)  = -x*(1.d0+x)*y*(1.d0-y)*0.25d0
        shg2(5,i,4)  = (1.d0-x*x)*y*(1.d0+y)*0.5d0
        shg2(6,i,4)  = -x*(1.d0-x)*(1.d0-y*y)*0.5d0
        shg2(7,i,4)  = -(1.d0-x*x)*y*(1.d0-y)*0.5d0
        shg2(8,i,4)  =  x*(1.d0+x)*(1.d0-y*y)*0.5d0
        shg2(9,i,4)  =  (1.d0-x*x)*(1.d0-y*y)
        shg2x(1,i,4) =  (1.d0+2.d0*x)*y*(1.d0+y)*0.25d0
        shg2x(2,i,4) = -(1.d0-2.d0*x)*y*(1.d0+y)*0.25d0
        shg2x(3,i,4) =  (1.d0-2.d0*x)*y*(1.d0-y)*0.25d0
        shg2x(4,i,4) = -(1.d0+2.d0*x)*y*(1.d0-y)*0.25d0
        shg2x(5,i,4) = -x*y*(1.d0+y)
        shg2x(6,i,4) = -(1.d0-2.d0*x)*(1.d0-y*y)*0.5d0
        shg2x(7,i,4) =  x*y*(1.d0-y)
        shg2x(8,i,4) =  (1.d0+2.d0*x)*(1.d0-y*y)*0.5d0
        shg2x(9,i,4) = -2.d0*x*(1.d0-y*y)
        shg2y(1,i,4) =  x*(1.d0+x)*(1.d0+2.d0*y)*0.25d0
        shg2y(2,i,4) = -x*(1.d0-x)*(1.d0+2.d0*y)*0.25d0
        shg2y(3,i,4) =  x*(1.d0-x)*(1.d0-2.d0*y)*0.25d0
        shg2y(4,i,4) = -x*(1.d0+x)*(1.d0-2.d0*y)*0.25d0
        shg2y(5,i,4) =  (1.d0-x*x)*(1.d0+2.d0*y)*0.5d0
        shg2y(6,i,4) =  x*(1.d0-x)*y
        shg2y(7,i,4) = -(1.d0-x*x)*(1.d0-2.d0*y)*0.5d0
        shg2y(8,i,4) = -x*(1.d0+x)*y
        shg2y(9,i,4) = -2.d0*y*(1.d0-x*x)
        we2(i,4) = we(i)
c
c       serendipity element
        shg2(1,i,5)  = (1.d0+x)*(1.d0+y)*(x+y-1.d0)*0.25d0
        shg2(2,i,5)  = (1.d0-x)*(1.d0+y)*(y-x-1.d0)*0.25d0
        shg2(3,i,5)  = (1.d0-x)*(1.d0-y)*(-x-y-1.d0)*0.25d0
        shg2(4,i,5)  = (1.d0+x)*(1.d0-y)*(x-y-1.d0)*0.25d0
        shg2(5,i,5)  = (1.d0-x*x)*(1.d0+y)*0.5d0
        shg2(6,i,5)  = (1.d0-x)*(1.d0-y*y)*0.5d0
        shg2(7,i,5)  = (1.d0-x*x)*(1.d0-y)*0.5d0
        shg2(8,i,5)  = (1.d0+x)*(1.d0-y*y)*0.5d0
        shg2(1,i,5)  = (1.d0+y)*(2.d0*x+y)*0.25d0
        shg2x(2,i,5) = (1.d0+y)*(2.d0*x-y)*0.25d0
        shg2x(3,i,5) = (1.d0-y)*(2.d0*x+y)*0.25d0
        shg2x(4,i,5) = (1.d0-y)*(2.d0*x-y)*0.25d0
        shg2x(5,i,5) = -x*(1.d0+y)
        shg2x(6,i,5) = -(1.d0-y*y)*0.5d0
        shg2x(7,i,5) = -x*(1.d0-y)
        shg2x(8,i,5) = (1.d0-y*y)*0.5d0
        shg2y(1,i,5) = (1.d0+x)*(2.d0*y+x)*0.25d0
        shg2y(2,i,5) = (1.d0-x)*(2.d0*y-x)*0.25d0
        shg2y(3,i,5) = (1.d0-x)*(2.d0*y+x)*0.25d0
        shg2y(4,i,5) = (1.d0+x)*(2.d0*y-x)*0.25d0
        shg2y(5,i,5) = (1.d0-x*x)*0.5d0
        shg2y(6,i,5) = -(1.d0-x)*y
        shg2y(7,i,5) = -(1.d0-x*x)*0.5d0
        shg2y(8,i,5) = -(1.d0+x)*y
        we2(i,5) = we(i)
c
      enddo
c
      return
      end
