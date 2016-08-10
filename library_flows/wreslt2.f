c***********************************************************************
      subroutine WRESLT2 ( it,t, nbrigid,xpos, xcv,ycv, strsx,strsy,
     &                     nnode,nvert,nelem,inod,x,y,
     &                     nbd,ibdnod,nic,ic, nbound,nside,bdseg,pwall,
     &                     ngmx,ncpg,ncps,ncpb,ncpvelo,ncpelas,
     &                     u,p,strm,vort,selas, iou)
c
c     This routine writes the results into a file 'hhxxxx.fmt'
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
c
      integer i,k,lenstr,j
      character*10 auxnam
      character*4 strval
c
      close(unit=iou)
      call cvnust(it,4,strval,lenstr)
      auxnam = 'hh'//strval(1:lenstr)//'.fmt'
      open(unit=iou,file=auxnam)
        write(iou,1010) it,t,nbrigid
        write(iou,1011) nnode,nvert,nelem,nbd,nic
        write(iou,1013) ((inod(i,k),i=1,3),k=1,nelem)
        write(iou,1012) (ibdnod(k),k=1,nbd)
        write(iou,1012) (ic(k),k=1,nic+1)
        write(iou,1015) (x(i),i=1,nvert)
        write(iou,1015) (y(i),i=1,nvert)
        write(iou,1015) (u(1,i),i=1,nvert)
        write(iou,1015) (u(2,i),i=1,nvert)
        write(iou,1015) (p(i),i=1,nvert)
        write(iou,1015) (strm(i),i=1,nvert)
        write(iou,1015) (vort(i),i=1,nvert)
        write(iou,1015) ((xpos(i,k),i=1,3),k=1,nbrigid)
        write(iou,1011) ncpelas
        do j=1,ncpelas
          write(iou,1015) (selas(j,i),i=1,nvert)
        enddo
        write(iou,1011) 0
        write(iou,1011) nbd
        do i=1,nbd
          write(iou,1014) xcv(i),ycv(i),strsx(i),strsy(i)
        enddo
c
        write(iou,1012) nbound 
        write(iou,1012) (nside(i),i=1,nbound)
        write(iou,1012) bdseg
        do i=1,bdseg
           write(iou,1015) (pwall(j,i),j=1,5)
        enddo
c
      close(unit=iou)
c
      return
c
1010  format(i5,1pe13.5,i5)
1011  format(5i7)
1012  format(13i6)
1013  format(3i6)
1014  format(4(1pe13.5))
1015  format(6(1pe13.5))
c
      end
