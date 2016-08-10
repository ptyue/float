c***********************************************************************
      subroutine MSHOLD ( nnode,nvert,nelem,inod,x,y,nec, ngmx,ncmx,
     &                    nbd,ibdnod,nic,ic,
     &                    nnodd,nvrtd,nelmd,inodd,xold,yold,necold,
     &                    nbdd,ibddd,icold )
c
c     This routine makes a copy of the mesh
c
c     INPUT:
c       nnode,nvert,nelem,inod,x,y,nec = mesh information
c
c     OUTPUT:
c       nnodd,nvrtd,nelmd,inodd,xold,yold,necold = copy of the above
c***********************************************************************
      implicit none
      integer ngmx,ncmx
      integer nnode,nvert,nelem,inod(ngmx,nelem),nec(ncmx,nelem),
     &        nnodd,nvrtd,nelmd,inodd(ngmx,nelem),necold(ncmx,nelem),
     &        nbd,ibdnod(nbd),nic,ic(nic+1),nbdd,ibddd(nbd),icold(nic+1)
      real*8 x(nnode),y(nnode),xold(nnode),yold(nnode)
c
      integer i,j
c
      nvrtd = nvert
      nnodd = nnode
      nelmd = nelem
c
      do i=1,nnode
        xold(i) = x(i)
        yold(i) = y(i)
      enddo
c
      do j=1,ngmx
        do i=1,nelem
          inodd(j,i) = inod(j,i)
        enddo
      enddo
c
      do j=1,ncmx
        do i=1,nelem
          necold(j,i) = nec(j,i)
        enddo
      enddo
c
      nbdd = nbd
      do i=1,nbd
        ibddd(i) = ibdnod(i)
      enddo
      do i=1,nic+1
        icold(i) = ic(i)
      enddo
c
      return
      end
