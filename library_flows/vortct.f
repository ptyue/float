c***********************************************************************
      subroutine VORTCT ( vort, u,dudx,dudy,ndgvelo,ngmx,ncpg,
     &                  nbrigid,nbfluid,nbelast,nvert,nnode,nelem,inod,
     &                  reft,x,y,selas,ncps,ndgelas,pelas,icoe,dt,
     &                  nbd,ibdnod,xcv,ycv,strsx,strsy, iwk,rwk )

c
c     This routine calculates vorticity, W = dV/dx - dU/dy
c     INPUT:
c        nvert = dimension of vorticity
c        u(2,ndgvelo) = velocity
c        x(nnode), y(nnode), inod(6,nelem) = mesh information
c
c     OUTPUT:
c       vort(nvert) = vorticity
c       dudx(2,ndgvelo) = velocity gradient dudx
c       dudy(2,ndgvelo) = velocity gradient dudy
c
c        iwk,rwk(nnode) = working array
c***********************************************************************
      implicit none
      integer nvert,nnode,nelem,ndgvelo,ngmx,inod(ngmx,nelem), ncpg
     &        ,nbrigid,nbfluid,nbelast,reft(nelem),iwk(nnode)
      real*8 vort(nvert),x(nnode),y(nnode),
     &       u(ncpg,ndgvelo),dudx(ncpg,ndgvelo),dudy(ncpg,ndgvelo),
     &       rwk(nnode)
      integer ndgelas,ncps,nbd,ibdnod(nbd), icoe
      double precision selas(ncps,ndgelas),pelas(10,icoe),dt
      real*8 xcv(*),ycv(*),strsx(*),strsy(*)
c
      integer i,k
      double precision tmp,t1
c
c     velocity gradient
c     -----------------
      call GRADU ( u,dudx,dudy,ndgvelo, x,y,inod,nnode,nelem, 
     &             reft, iwk,rwk )
c
c     vorticity
c     ---------
      do i=1,nvert
        vort(i) = dudx(2,i)-dudy(1,i)
      enddo
c
c     stress
c     ------
      if ( nbelast.gt.0 ) then
         tmp = dt*0.5*pelas(1,icoe)/(1.d0+pelas(2,icoe)) 
         do i=1,ndgelas
            if ( iwk(i).gt.1+nbrigid+nbfluid ) then
              t1 = 0.5*(dudx(1,i)+dudy(2,i))
              selas(1,i) = selas(1,i) + 2.0*(dudx(1,i)-t1)*tmp
              selas(2,i) = selas(2,i) + (dudy(1,i)+dudx(2,i))*tmp
              selas(3,i) = selas(3,i) + 2.0*(dudy(2,i)-t1)*tmp
            endif
         enddo 
      endif
c
      do i=1,nbd
         k = ibdnod(i)
         xcv(i) = x(k)
         ycv(i) = y(k)
         strsx(i) = dudx(1,k)
         strsy(i) = dudy(1,k) + dudx(2,k)
      enddo
c
      return
      end
c
c***********************************************************************
      subroutine GRADU ( f, dfdx, dfdy, n, x, y,inod, nnode,nelem, 
     &                   reft, indx,rwk )
c
c     calculates the gradient for input f
c     INPUT:
c        n = dimension of vector f
c        f(2,n) = input function
c        x(nnode), y(nnode), inod(6,nelem) = mesh information
c     OUTPUT:
c        dfdx(2,n), dfdx(2,n)
c
c        indx(nnode), rwk(nnode) = working array
c
c     F = F(i)*Psi(i), Fx = F(i)*Px(i)  &  Fy = F(i)*Py(i)
c     Px = Pxi*XIx + Pet*ETx , Py = Pxi*XIy + Pet*ETy
c     Xxi = X(i)*Pxi(i), Xet = X(i)*Pet(i),
c     Yxi = Y(i)*Pxi(i), Yet = Y(i)*Pet(i)
c     -> XIx = Yet/J, XIy = -Xet/J, ETx = -Yxi/J, ETy = Xxi/J
c      and  J = Xxi*Yet - Xet*Yxi
c***********************************************************************
      implicit none
      integer nnode,nelem,n,inod(6,nelem), reft(nelem)
     &        ,indx(nnode)
      real*8 x(nnode),y(nnode),f(2,n),dfdx(2,n),dfdy(2,n), rwk(nnode)
c
      integer ng,nd,k,i,ne,nn,j,ndloc, ref
      real*8 pg(9),pgxi(9),pget(9),pgx(9,9),pgy(9,9),puxi(9,9),puet(9,9)
      real*8 xi(6),et(6),xn(6),yn(6), fn(2,9)
      real*8 x0,dxxi,dyet,dxet,dyxi,ajac,t1,t2,t3,t4,xcoor, tmp
      data xi/0.,1.,0.,0.5,0.5,0./, et/0.,0.,1.,0.,0.5,0.5/
c
      ng = ndloc(0,2)
      nd = ng
c
c     calculate the shape functions
c     -----------------------------
      do k=1,nd
        call SHAPE2(ng,xi(k),et(k),pg,pgxi,pget,1)
        do i=1,ng
          pgx(i,k) = pgxi(i)
          pgy(i,k) = pget(i)
        enddo
        call SHAPE2(nd,xi(k),et(k),pg,pgxi,pget,1)
        do i=1,nd
          puxi(i,k) = pgxi(i)
          puet(i,k) = pget(i)
        enddo
      enddo
c
c     initialize
c     ----------
      do i=1,n
        dfdx(1,i) = 0.
        dfdx(2,i) = 0.
        dfdy(1,i) = 0.
        dfdy(2,i) = 0.
        rwk(i) = 0.0
        indx(i) = 0
      enddo
c
      do ne=1,nelem
         ref = reft(ne)
         do j=1,nd
            nn = inod(j,ne)
            if ( nn.le.n) indx(nn) = max(ref,indx(nn))
         enddo
      enddo
c
c     for each element
c     ----------------
      do ne =1,nelem
        ref = reft(ne)
c
        x0 = x(inod(1,ne))
        do k=1,ng
          nn = inod(k,ne)
          xn(k) = xcoor(x0,x(nn))
          yn(k) = y(nn)
        enddo
        do j=1,2
          do k=1,nd
            fn(j,k) = f(j,inod(k,ne))
          enddo
        enddo
c
        do j=1,nd
          nn = inod(j,ne)
          if ( indx(nn).eq.ref ) then
            dxxi = 0.d0
            dxet = 0.d0
            dyxi = 0.d0
            dyet = 0.d0
            do k=1,ng
               dxxi = dxxi + xn(k)*pgx(k,j)
               dxet = dxet + xn(k)*pgy(k,j)
               dyxi = dyxi + yn(k)*pgx(k,j)
               dyet = dyet + yn(k)*pgy(k,j)
            enddo
            ajac = 1.0/(dxxi*dyet-dxet*dyxi+1.d-40)
            dxxi = dxxi*ajac
            dyet = dyet*ajac
            dxet = dxet*ajac
            dyxi = dyxi*ajac
c
            t1 = 0.d0
            t2 = 0.d0
            t3 = 0.d0
            t4 = 0.d0
            do k=1,nd
               tmp = puxi(k,j)*dyet-puet(k,j)*dyxi
               t1 = t1 + fn(1,k)*tmp
               t3 = t3 + fn(2,k)*tmp
               tmp = -puxi(k,j)*dxet+puet(k,j)*dxxi
               t2 = t2 + fn(1,k)*tmp
               t4 = t4 + fn(2,k)*tmp
            enddo
c
c            rwk(nn)    = rwk(nn)    + ajac
c            dfdx(1,nn) = dfdx(1,nn) + t1*ajac
c            dfdy(1,nn) = dfdy(1,nn) + t2*ajac
c            dfdx(2,nn) = dfdx(2,nn) + t3*ajac
c            dfdy(2,nn) = dfdy(2,nn) + t4*ajac
            rwk(nn)    = rwk(nn)    + 1
            dfdx(1,nn) = dfdx(1,nn) + t1
            dfdy(1,nn) = dfdy(1,nn) + t2
            dfdx(2,nn) = dfdx(2,nn) + t3
            dfdy(2,nn) = dfdy(2,nn) + t4

          endif
c
        enddo
      enddo
c
      do i=1,nnode
         tmp = 1.d0/rwk(i)
         dfdx(1,i) = dfdx(1,i)*tmp
         dfdy(1,i) = dfdy(1,i)*tmp
         dfdx(2,i) = dfdx(2,i)*tmp
         dfdy(2,i) = dfdy(2,i)*tmp
c	write(20,'(2i5,4e13.5)') i,indx(i),dfdx(1,i),
c     &       dfdy(1,i)+dfdx(2,i),dfdy(2,i),dfdx(1,i)+dfdy(2,i)
      enddo
c
      return
      end
c
c
c***********************************************************************
      subroutine SHAPE2 ( nd, x,y, psi,psix,psiy, indx )
c
c     this routine calculates the shape functions
c        nd=6: for triangular element with 6 nodes (area coor.)
c        nd=9: for quadrilateral element with 9 nodes
c***********************************************************************
      implicit none
      integer nd, indx
      double precision psi(9),psix(9),psiy(9),x,y
c
      if ( nd.eq.3 ) then
        psi(1) = 1.d0-x-y
        psi(2) = x
        psi(3) = y
        if ( indx.eq.0 ) return
        psix(1) =-1.d0
        psix(2) = 1.d0
        psix(3) = 0.d0
        psiy(1) =-1.d0
        psiy(2) = 0.d0
        psiy(3) = 1.d0
      endif
c
      if ( nd.eq.6 ) then
        psi(1)  = 1.d0-3.d0*(x+y)+2.d0*(x+y)*(x+y)
        psi(2)  = x*(2.d0*x-1.d0)
        psi(3)  = y*(2.d0*y-1.d0)
        psi(4)  = 4.d0*x*(1.d0-x-y)
        psi(5)  = 4.d0*x*y
        psi(6)  = 4.d0*y*(1.d0-x-y)
        if ( indx.eq.0 ) return
        psix(1) =-3.d0+4.d0*(x+y)
        psix(2) = 4.d0*x-1.d0
        psix(3) = 0.d0
        psix(4) = 4.d0-8.d0*x-4.d0*y
        psix(5) = 4.d0*y
        psix(6) =-4.d0*y
        psiy(1) =-3.d0+4.d0*(x+y)
        psiy(2) = 0.d0
        psiy(3) = 4.d0*y-1.d0
        psiy(4) =-4.d0*x
        psiy(5) = 4.d0*x
        psiy(6) = 4.d0-4.d0*x-8.d0*y
      endif
c
      if ( nd.eq.9 ) then
        psi(1)  = x*(1.d0+x)*y*(1.d0+y)*0.25d0
        psi(2)  =-x*(1.d0-x)*y*(1.d0+y)*0.25d0
        psi(3)  = x*(1.d0-x)*y*(1.d0-y)*0.25d0
        psi(4)  =-x*(1.d0+x)*y*(1.d0-y)*0.25d0
        psi(5)  = (1.d0-x*x)*y*(1.d0+y)*0.5d0
        psi(6)  =-x*(1.d0-x)*(1.d0-y*y)*0.5d0
        psi(7)  =-(1.d0-x*x)*y*(1.d0-y)*0.5d0
        psi(8)  = x*(1.d0+x)*(1.d0-y*y)*0.5d0
        psi(9)  = (1.d0-x*x)*(1.d0-y*y)
        if ( indx.eq.0 ) return
        psix(1) = (1.d0+2.d0*x)*y*(1.d0+y)*0.25d0
        psix(2) =-(1.d0-2.d0*x)*y*(1.d0+y)*0.25d0
        psix(3) = (1.d0-2.d0*x)*y*(1.d0-y)*0.25d0
        psix(4) =-(1.d0+2.d0*x)*y*(1.d0-y)*0.25d0
        psix(5) =-x*y*(1.d0+y)
        psix(6) =-(1.d0-2.d0*x)*(1.d0-y*y)*0.5d0
        psix(7) = x*y*(1.d0-y)
        psix(8) = (1.d0+2.d0*x)*(1.d0-y*y)*0.5d0
        psix(9) =-2.d0*x*(1.d0-y*y)
        psiy(1) = x*(1.d0+x)*(1.d0+2.d0*y)*0.25d0
        psiy(2) =-x*(1.d0-x)*(1.d0+2.d0*y)*0.25d0
        psiy(3) = x*(1.d0-x)*(1.d0-2.d0*y)*0.25d0
        psiy(4) =-x*(1.d0+x)*(1.d0-2.d0*y)*0.25d0
        psiy(5) = (1.d0-x*x)*(1.d0+2.d0*y)*0.5d0
        psiy(6) = x*(1.d0-x)*y
        psiy(7) =-(1.d0-x*x)*(1.d0-2.d0*y)*0.5d0
        psiy(8) =-x*(1.d0+x)*y
        psiy(9) =-2.d0*y*(1.d0-x*x)
      endif
c
      return
      end
