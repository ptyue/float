c***********************************************************************
      subroutine dropspread(t,ntime,npfmt,io,iout,ishm,
     &                  nvert,nelem,ngmx,ncmx,
     &                  inod,nec,nrdgeom,nrdphi,ndggeom,ndgphi,x,y,phi,
     &                  nic,ic,nbd,ibdnod,iwork,rwork)
c***********************************************************************
c   This subroutine output the interface information for drop spreading
c
c     Input:
c        t:       time
c        ntime:   the n-th time step
c        npfmt:   inter?????.dat is output every npfmt steps
c        io:      the file is opened and closed within this subroutine
c                 (input as io5)
c        iout:    i/o unit for error messages(input as io6)
c        ishm:    0 for triangular cells and 1 for quadrilateral cells
c        phi:     phase-field variables
c        ....
c     Output:
c        two data files
c        dropspread.dat: n, time, drop radius, drop height
c        inter?????.dat: the position of interfacial points             
c     Working:
c        iwork,rwork
c     Written Feb 23, 2005 by Pengtao Yue
c***********************************************************************
      implicit none
      integer io,iout,ishm,nvert,nelem,ngmx,ncmx,
     &        nrdgeom,nrdphi,ndggeom,ndgphi,ntime,npfmt,nic,nbd
      integer inod(ngmx,nelem),nec(ncmx,nelem),iwork(nelem),ic(nic+1),
     &        ibdnod(nbd)
      real(8) t
      real(8) x(ndggeom),y(ndggeom),phi(ndgphi),rwork(nelem*2+nvert)
c
      integer npt,lrintr,lrphi,i
      real(8) f0,r,h
      character buffer*5,filename*14
c           
      f0=0.d0
      lrintr=1
      lrphi=lrintr+2*nelem
c     


      call dropspreadrh(nrdgeom,ndggeom,x,y,nrdphi,ndgphi,phi,
     &                     nic,ic,nbd,ibdnod,f0,r,h)

      open(io,file='dropspread.dat',access='append')
      write(io,100)ntime, t,r,h
100   format(1x,i5, 3(1x,g12.5))
	close(io)
      if(mod(ntime,npfmt)/=0)return
      write(buffer,'(i5.5)')ntime
      filename='inter'//buffer//'.dat'
      open(io,file=filename)
      call getintrfcbd(iout,ishm,nvert,nelem,ngmx,inod,ncmx,nec,
     &                 x,y,phi,f0,npt,rwork(lrintr),iwork,rwork(lrphi))
      do i=1,npt
         write(io,'(2(1x,g14.6))')rwork(lrintr+2*(i-1):lrintr+2*(i-1)+1)
      enddo
!	write(*,*)filename
      close(io)
      end


c***********************************************************************
      subroutine dropspreadrh(nrdgeom,ndggeom,x,y,nrdphi,ndgphi,phi,
     &                     nic,ic,nbd,ibdnod,f0,r,h)
c
c     calcaute drop spreading radius and drop height, 
c     the drop is centered at (0,0)
c     Input:
c     Output:  
c        r:    radius of drop
c        h:    height of drop
c     Note:  applicable range: 2>=nrdgeom>=nrdphi
c            center of the drop is (0,0)
c
c     Mar 31, 2009, Pengtao Yue
c        modified from subroutine DropDfExt
c***********************************************************************
      implicit none
      integer nrdgeom,nrdphi,ndggeom,ndgphi,nic,nbd
      integer ic(nic+1),ibdnod(nbd)
      real*8  x(ndggeom),y(ndggeom),phi(ndgphi),f0,r,h
c
      real(8) xr,yr,f1,f2,c1,c2
      integer i,n1,n2,k
c
c      f0=0.d0
      xr=1.
      yr=1.
      if(nrdgeom.eq.2.and.nrdphi.eq.1)then
         k=2
      else
         k=1
      endif
c     along x axis
      do i=ic(1),ic(2)-k-1,k
         n1=ibdnod(i)
         n2=ibdnod(i+k)
         f1=phi(n1)
         f2=phi(n2)
         if(min(f1,f2).gt.f0.or.max(f1,f2).lt.f0)cycle
         if(f1.eq.f0)then
            xr=x(n1)
            exit
         elseif(f2.eq.f0)then
            xr=x(n2)
            exit
         else
            c1=(f2-f0)/(f2-f1)
            c2=(f0-f1)/(f2-f1)
            xr=c1*x(n1)+c2*x(n2)
            exit
         endif
      enddo
c     along y axis
      do i=ic(4),ic(5)-k-1,k
         n1=ibdnod(i)
         n2=ibdnod(i+k)
         f1=phi(n1)
         f2=phi(n2)
         if(min(f1,f2).gt.f0.or.max(f1,f2).lt.f0)cycle
         if(f1.eq.f0)then
            yr=y(n1)
            exit
         elseif(f2.eq.f0)then
            yr=y(n2)
            exit
         else
            c1=(f2-f0)/(f2-f1)
            c2=(f0-f1)/(f2-f1)
            yr=c1*y(n1)+c2*y(n2)
            exit
         endif
      enddo
      r=xr
      h=yr
c      open(io,file='dropspread.dat',access='append')
c      write(io,100)t,r,h
c100   format(8(1x,g12.5))
c      close(io)
      end
