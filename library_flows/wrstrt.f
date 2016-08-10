c***********************************************************************
      subroutine WRSTRT ( iou, it,t,dt, nbrigid,xpos,uvom,duvot,fxym,
     &                    ndgvelo,ndgpres,ndgelas,ndgphi,ncpvelo,
     &                    ncpelas,nbelast,ncpphi,
     &                    nvert,nnode,nelem,inod,x,y,nec,area,aspr,reft,
     &                    nbd,ibdnod,nic,ic,nbound,nside,pwall,bdseg,
     &                    u,p,dudt,umesh,selas,desdt,phi,psi,dphidt,
     &                    ngmx,ncmx,ncpg,ncps,ncpb,
     &                    noflow,ncpveln,ncpelan,
     &                    dist,meshgen)


c
c     this routine writes a restart file 'hh.stt'
c***********************************************************************
      implicit none
      integer nnode,nvert,nelem,ngmx,ncmx,inod(ngmx,nelem),
     &        nbd,ibdnod(nbd),nic,ic(nic+1),nbound,nside(nbound),
     &        nec(ncmx,nelem),reft(nelem),ndgvelo,ndgpres,ndgelas,
     &        ncpvelo,ncpelas, it,nbrigid, iou,
     &        ncpg,ncps,ncpb, bdseg, nbelast
      real*8 t,dt,x(nnode),y(nnode),area(nelem),aspr(nelem),
     &       p(ndgpres),u(ncpg,ndgvelo),dudt(ncpg,ndgvelo),
     &       umesh(ncpg,ndgvelo),selas(ncps,ndgelas),
     &       desdt(ncps,ndgelas),
     &       xpos(ncpb,nbrigid),uvom(ncpb,nbrigid),fxym(ncpb,nbrigid),
     &       duvot(ncpb,nbrigid),pwall(5,bdseg)
c     for initialization      
      logical noflow
      integer ncpveln,ncpelan
c*CH start
      integer ndgphi, ncpphi
      real*8 phi(ndgphi), psi(ndgphi), dphidt(ndgphi)
c*CH end
      real(8) dist(nvert)
      integer meshgen
c
      integer i,j,k
c
      open(unit=iou,file='hh.stt',form='unformatted')
        write(iou) it,t,dt
        write(iou) ((xpos(i,k),i=1,ncpb),k=1,nbrigid)
        write(iou) ((uvom(i,k),i=1,ncpb),k=1,nbrigid)
        write(iou) ((fxym(i,k),i=1,ncpb),k=1,nbrigid)
        write(iou) ((duvot(i,k),i=1,ncpb),k=1,nbrigid)
        write(iou) nvert,nnode,nelem
        write(iou) ((inod(k,i),k=1,ngmx),i=1,nelem)
        write(iou) ((nec(k,i),k=1,ncmx),i=1,nelem)
        write(iou) (x(i),y(i),i=1,nnode)
        write(iou) (area(i),aspr(i),i=1,nelem)
        write(iou) (reft(i),i=1,nelem)
        write(iou) nbd,nic,nbound
        write(iou) (ibdnod(i),i=1,nbd)
        write(iou) (ic(i),i=1,nic+1)
        write(iou) (nside(i),i=1,nbound)
        write(iou) ndgvelo
        if(noflow)then
        write(iou) ((u(j,i),i=1,ndgvelo),j=1,ncpveln)
        write(iou) ((dudt(j,i),i=1,ndgvelo),j=1,ncpveln)
        else
        write(iou) ((u(j,i),i=1,ndgvelo),j=1,ncpvelo)
        write(iou) ((dudt(j,i),i=1,ndgvelo),j=1,ncpvelo)
        endif
        write(iou) ((umesh(j,i),i=1,ndgvelo),j=1,ncpg)
        write(iou) ndgpres
        write(iou) (p(i),i=1,ndgpres)
c
        write(iou) ndgelas
        if(noflow)then
        do j=1,ncpelan
          write(iou) (selas(j,i),i=1,ndgelas)
          write(iou) (desdt(j,i),i=1,ndgelas)
        enddo
        else
        do j=1,ncpelas
          write(iou) (selas(j,i),i=1,ndgelas)
          write(iou) (desdt(j,i),i=1,ndgelas)
        enddo
        endif
c
        if ( nbelast.gt.0 ) then
          do j=1,3
            write(iou) (selas(j,i),i=1,ndgelas)
          enddo
        endif
c*CH start     modifications related with C-H equation
        if(ncpphi.eq.1)then
          write(iou) ndgphi
          write(iou) phi(1:ndgphi)
          write(iou) psi(1:ndgphi)
          write(iou) dphidt(1:ndgphi)
        endif
c*CH end
        write(iou) bdseg
        do j=1,bdseg
          write(iou) (pwall(i,j),i=1,5)
        enddo
        if(meshgen==2)write(iou) dist(1:nvert)
c
      close(iou)
c
      return
      end
