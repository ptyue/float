c***********************************************************************
      subroutine masbal ( iou, kaxi, nrdgeom,ndggeom,x,y, ndgvelo,u,v, 
     &                   nbd,ibdnod,nic,ic,nbound,nside, shg1,shg1x,we1)
c
c     this subroutine calculates the mass balance 
c     May 14, 1996, Howard Hu
c***********************************************************************
      implicit none
      integer iou,kaxi, nrdgeom,ndggeom,ndgvelo,
     &        nbd,ibdnod(nbd),nic,ic(nic+1),nbound,nside(nbound)
      real*8 u(ndgvelo),v(ndgvelo),x(ndggeom),y(ndggeom)
      real*8 shg1(3,3,2),shg1x(3,3,2),we1(3,2)
c
      integer nidx,nd,is,ica1,ica,icb,nb,kk,nel,m,k,n,j,i,ib
      real*8 xn(3),yn(3),b1dx(3,3),b1dy(3,3),b1ds(3,3),
     &       fact,c,fp,fn, xcoor
c
c     preparation
c     -----------
      nidx = nrdgeom
      nd = nidx + 1
      fact = 1.d0
      if ( kaxi.eq.1 ) fact = 2.d0*3.1415926536d0
c
      write(iou,*) ' '
      write(iou,*) '******************************'
      write(iou,*) '*                            *'
      write(iou,*) '*   table of mass balance    *'
      write(iou,*) '*                            *'
      write(iou,*) '******************************'
      write(iou,1311)
      write(iou,1310)
      write(iou,1303)
      write(iou,1310)
      write(iou,1311)
      write(iou,1310)
c
c     calculation of flow rates
c     -------------------------
      fp = 0.d0
      fn = 0.d0
      is = 0
      do nb=1,nbound
        ica1 = ic(is+1)
        do kk=1,nside(nb)
          c = 0.d0
          is = is + 1
          ica = ic(is)
          icb = ic(is+1)
          nel = (icb-ica)/nidx
          do m=1,nel
            ib = ica + nidx*(m-1)
c
            n = ibdnod(ib)
            xn(1) = x(n)
            yn(1) = y(n)
            do i=2,nd
              k = ib + i - 1
              if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
              n = ibdnod(k)
              xn(i) = xcoor(xn(1),x(n))
              yn(i) = y(n)
            enddo
c
            call mat1ds (xn,yn,nidx,kaxi,b1dx,b1dy,b1ds,shg1,shg1x,we1)
            do i=1,nd
              do j=2,nd
                b1dx(i,1) = b1dx(i,1) + b1dx(i,j)
                b1dy(i,1) = b1dy(i,1) + b1dy(i,j)
              enddo
            enddo
c
            do i=1,nd
              k = ib + i - 1
              if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
              n = ibdnod(k)
              c = c + u(n)*b1dy(i,1) - v(n)*b1dx(i,1)
            enddo
c
          enddo
c
          c = c*fact
          if ( c.lt.0.d0 ) then
            c = dabs(c)
            fn = fn + c
            write(iou,1306) is,ica,icb,c
          else
            fp = fp + c
            write(iou,1304) is,ica,icb,c
          endif
c
        enddo
      enddo
c
      write(iou,1310)
      write(iou,1314) fn,fp
      write(iou,1310)
      write(iou,1311)
      write(iou,1307)
      write(iou,1316) fp-fn
      write(iou,1307)
      write(iou,1311)
c
1303  format('+ side +  1st corner  +  2nd corner  +',
     &       '     inflow     +    outflow     +')
1304  format('+ ',i4,' + ',2(4x,i4,4x,' + '),6x,'-',7x,' + ',
     &       1pe14.7,' +')
1306  format('+ ',i4,' + ',2(4x,i4,4x,' + '),1pe14.7,' + ',6x,'-',
     &       7x,' +')
1307  format('+',70x,'+')
1310  format('+ ',4x,' + ',2(13x,'+ '),14x,' + ',14x,' +')
1311  format(72('-'))
1314  format('+ ',4x,' + ',2(13x,'+ '),14('-'),' + ',14('-'),' +'/
     &       '+ ',4x,' + ',2(13x,'+ '),14x,    ' + ',14x,    ' +'/
     &       '+ ',4x,' + ',2(13x,'+ '),1pe14.7,' + ',1pe14.7,' +')
1316  format('+ absolute error on mass balance :',e10.3,27x,'+')
c
      return
      end
c
c***********************************************************************
      subroutine mat1ds ( xn,yn,nrd,kaxi,b1dx,b1dy,b1ds,shg1,shg1x,we1 )
c
c     this routine performs s-integration along an arbitrary boundary
c     January 26, 1996, by Howard H. Hu
c***********************************************************************
      implicit none
      integer nrd, kaxi
      real*8 xn(nrd+1),yn(nrd+1)
      real*8 b1dx(3,3),b1dy(3,3),b1ds(3,3),
     &       shg1(3,3,2),shg1x(3,3,2),we1(3,2)
c
      real*8 xx,dxxic,dyxic,dsxic,alpha,awx,awy,aws,awxj,awyj,awsj
      integer nd,i,j,l
c
c     0. preparations
c     ---------------
      nd = nrd+1
c
c     2. initialization
c     -----------------
      do j=1,3
        do i=1,3
          b1dx(i,j) = 0.d0
          b1dy(i,j) = 0.d0
          b1ds(i,j) = 0.d0
        enddo
      enddo
c
c     3. loop over integration points
c     -------------------------------
      do l=1,3
c
c       calculating derivatives
c
        xx    = 0.
        dxxic = 0.
        dyxic = 0.
        do i=1,nd
          xx    = xx    + xn(i)*shg1(i,l,nrd)
          dxxic = dxxic + xn(i)*shg1x(i,l,nrd)
          dyxic = dyxic + yn(i)*shg1x(i,l,nrd)
        enddo
        dsxic = dsqrt(dxxic*dxxic + dyxic*dyxic)
        alpha = (1.d0 + kaxi*(xx-1.d0)) * we1(l,nrd)
        awx = alpha*dxxic
        awy = alpha*dyxic
        aws = alpha*dsxic
c
c       evaluating line integrals
c
        do j=1,nd
          xx = shg1(j,l,nrd)
          awxj = awx*xx
          awyj = awy*xx
          awsj = aws*xx
          do i=1,nd
            xx = shg1(i,l,nrd)
            b1dx(i,j) = b1dx(i,j) + awxj*xx
            b1dy(i,j) = b1dy(i,j) + awyj*xx
            b1ds(i,j) = b1ds(i,j) + awsj*xx
          enddo
        enddo
c
      enddo
c
      return
      end
