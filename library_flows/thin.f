c***********************************************************************
      subroutine pthin ( Nth,ivisc,pvis,upmax,wdth, ueta)
c
c     This routine computes the Poiseuille velocity profile of a shear
c       rate dependent fluid
c
c     Input:
c       ivisc       = index for viscosity function
c       pvis(10)    = viscoelastic parameters
c                     (for the continuous phase)
c       upmax, wdth = inflow conditions
c
c     Output:
c       ueta(1,2*Nth+1) = position of inflow section (y)
c       ueta(2,2*Nth+1) = velocity profile at inflow (u)
c       ueta(3,2*Nth+1) = shear rate at inflow (du/dy)
c       ueta(4,2*Nth+1) = viscosity profile at inflow (eta)
c***********************************************************************
      implicit none
      integer Nth,ivisc
      real*8 pvis(10),ueta(4,2*Nth+1),upmax,wdth
c
      integer ite,i,k
      real*8 dpdx,h,eta22,power,sum,zerr
c
      h = 0.5*wdth/Nth
c
      if ( ivisc.eq.2 ) then
        power = 0.5*(pvis(3)-1.d0)
      elseif ( ivisc.eq.3 ) then
        power = pvis(3)-1.d0
      else
        write(*,*) 'fully-developed flow is not implemented for',
     &             ' this fluid', ivisc
        stop
      endif
c
c     ueta(2,k) is temperary used as ueta(3,k) at former iteration
c
      ite = 0
      do k=1,Nth+1
        ueta(1,k) = h*(k-1)
        ueta(3,k) = 4.*upmax*(1.-2.*ueta(1,k)/wdth)/wdth
        ueta(4,k) = pvis(1)
        ueta(2,k) = ueta(3,k)
      enddo
c
50    continue
c
c     calculating the viscosity funtion
      if ( ivisc.eq.2 ) then
        do k=1,Nth+1
          eta22 = (pvis(2)*ueta(3,k))**2
          ueta(4,k) = pvis(4)+(pvis(1)-pvis(4))*(1.+eta22)**power
        enddo
      elseif ( ivisc.eq.3 ) then
        do k=1,Nth+1
          ueta(4,k) = pvis(1)*(dabs(ueta(3,k)+1.d-20))**power
        enddo
      endif
c
c     relating the velocity at the center to dpdx
      sum = -0.25*wdth/ueta(4,1)
      do i=2,Nth
        sum = sum + (ueta(1,i)-0.5*wdth)/ueta(4,i)
      enddo
      dpdx = upmax/(h*sum)
c      write(*,*) 'ite=',ite,' dpdx=',dpdx
c
c     calculating the shear-rate
      do k=1,Nth+1
        ueta(3,k) = dpdx*(ueta(1,k)-0.5*wdth)/ueta(4,k)
      enddo
c
c     check the convergence
      sum = 0.d0
      do k=1,Nth+1
        sum = sum + (ueta(3,k)-ueta(2,k))**2
      enddo
      zerr = dsqrt(sum)
      if ( zerr.le.1.d-10 ) go to 200
c
      do k=1,Nth+1
        ueta(2,k) = ueta(3,k)
      enddo
      ite = ite+1
      go to 50
c
200   continue
c
c     Trapezoidal integration
      ueta(2,1)=0.d0
      do k=2,Nth+1
        ueta(2,k) = ueta(2,k-1) + 0.5*h*(ueta(3,k-1)+ueta(3,k))
      enddo
c
      do k=2,Nth+1
        ueta(1,k+Nth) =  wdth-ueta(1,Nth+2-k)
        ueta(2,k+Nth) =  ueta(2,Nth+2-k)
        ueta(3,k+Nth) = -ueta(3,Nth+2-k)
        ueta(4,k+Nth) =  ueta(4,Nth+2-k)
      enddo
c
      return
      end
c
c***********************************************************************
      subroutine thin ( ycur,ueta,Nth,uthin,duthin,etathn )
c
c     find the values of velocity, its derivative and viscosity
c       using interpolation scheme
c
c     Input:
c       ycur    = current coordinate in y direction
c       ueta    = position, velocity, shear-rate and viscosity profiles
c     Output:
c       uthin   = velocity u(ycur) in the inflow section
c       duthin  = du/dy (ycur) in the inflow section
c       etathn  = viscosity (ycur) in the inflow section
c***********************************************************************
      implicit none
      integer Nth
      real*8 ueta(4,2*Nth+1), ycur,uthin,duthin,etathn
c
      integer k
      real*8 yrel,h,x1,x2
c
      h = ueta(1,2)-ueta(1,1)
      do k=1,2*Nth+1
        yrel = ycur-ueta(1,k)
        if ( dabs(yrel).le.1.d-8 ) then
          uthin = ueta(2,k)
          duthin= ueta(3,k)
          etathn= ueta(4,k)
          return
        endif
        if ( yrel.lt.0.d0 ) then
          x1 = yrel/h
          x2 = 1.d0-x1
          uthin  = x1*ueta(2,k-1) + x2*ueta(2,k)
          duthin = x1*ueta(3,k-1) + x2*ueta(3,k)
          etathn = x1*ueta(4,k-1) + x2*ueta(4,k)
          return
        endif
      enddo
c
      return
      end
