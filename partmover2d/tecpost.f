c***********************************************************************
      subroutine tecpost ( it,t, nbrigid,xpos, xcv,ycv, strsx,strsy,
     &                     nnode,nvert,nelem,inod,x,y,
     &                     nbd,ibdnod,nic,ic, nbound,nside,bdseg,pwall,
     &                     ngmx,ncpg,ncps,ncpb,ncpvelo,ncpelas,ncpphi,
     &                     u,p,strm,vort,selas,phi,psi,iou)
c     Jan, 2005, Pengtao Yue
c     This routine writes the results into a file 'tecxxxx.dat' for tecplot
c
c     Input :
c       it  : iteration courante
c       nnode,nvert,nelem,inod,x,y, nbd,ibdnod,nic,ic: mesh inf.
c       u,p : velocities and pressure
c       strm, vort: streamfunction and vorticity
c***********************************************************************
      implicit none
      integer ngmx,ncpg,ncps,ncpb
      integer it,ncpvelo,ncpelas,nbrigid,nnode,nvert,nelem,
     &        inod(ngmx,nelem),nbd,ibdnod(nbd),nic,ic(nic+1), iou
     &       ,nbound,nside(nbound),bdseg
      real*8 x(nnode),y(nnode),u(ncpg,nvert),p(nvert),
     &      strm(nvert),vort(nvert),selas(ncps,nvert),
     &      xpos(ncpb,nbrigid),t,xcv(nbd),ycv(nbd),strsx(nbd),strsy(nbd)
     &     ,pwall(5,bdseg)
c*CH start
      integer ncpphi
      real*8 phi(nvert),psi(nvert)
c*CH end
      real*8 xlength
      logical xperid
      common /xperid/ xlength,xperid
c
      integer i,k,lenstr,j, nn
      character*12 auxnam
      character*5 strval
      real*8 dx
c
      close(unit=iou)
      call cvnust(it,5,strval,lenstr)
      auxnam = 'tec'//strval(1:lenstr)//'.dat'
      if(xperid)then         
      nn=0
      do k=1,nelem
         dx=max(abs(x(inod(1,k))-x(inod(2,k))),
     &          abs(x(inod(2,k))-x(inod(3,k))),
     &          abs(x(inod(3,k))-x(inod(1,k))))
         if(dx.lt.0.5*xlength)then
            nn=nn+1
         endif
      enddo
      else
      nn=nelem
      endif

      open(unit=iou,file=auxnam)
        write(iou,'(a,1pe13.5,a)')'Title="partflow 2D, t=',t,'"'
        write(iou,'(a)')'Variables="x", "y", "u", "v", "p", "strm"'
        write(iou,'(a)')'"vort"'
        do j=1,ncpelas
        write(iou,'(5H"taup,I1,1H")')j
        enddo
        if(ncpphi.eq.1) write(iou,'(a)')'"phi", "psi"'
        write(iou,'(a,2(1x,a,i7))')'ZONE', 'N=', nvert, 'E=', nn
        write(iou,'(a)')'F=FEBLOCK, ET=TRIANGLE'
        
        write(iou,1015) (x(i),i=1,nvert)
        write(iou,1015) (y(i),i=1,nvert)
        write(iou,1015) (u(1,i),i=1,nvert)
        write(iou,1015) (u(2,i),i=1,nvert)
        write(iou,1015) (p(i),i=1,nvert)
        write(iou,1015) (strm(i),i=1,nvert)
        write(iou,1015) (vort(i),i=1,nvert)
        do j=1,ncpelas
          write(iou,1015) (selas(j,i),i=1,nvert)
        enddo
        if(ncpphi.eq.1) then
          write(iou,1015) phi(1:nvert)
          write(iou,1015) psi(1:nvert)
        endif
        if(xperid)then   
        do k=1,nelem
         dx=max(abs(x(inod(1,k))-x(inod(2,k))),
     &          abs(x(inod(2,k))-x(inod(3,k))),
     &          abs(x(inod(3,k))-x(inod(1,k))))
         if(dx.lt.0.5*xlength)then
          write(iou,1013) (inod(i,k),i=1,3)
         endif
        enddo
        else
         write(iou,1013) ((inod(i,k),i=1,3),k=1,nelem)
        endif
c
      close(unit=iou)
c
      return
c
1013  format(3i6)
1015  format(6(1pe13.5))
c
      end
