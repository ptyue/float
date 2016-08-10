c***********************************************************************
      subroutine prtmsh ( filename,ish,m, nvert,nnode,nelem,inod,x,y,
     &                    nref,nbd,ibdnod,nic,ic,nbound,nside,nec,iout )
c
c     INPUT
c       filename: name for the mesh file
c       nnode  : number of nodes
c       nvert  : number of vertices
c       nelem  : number of elements
c       inod   : element description table
c       x,y    : x and y coordinates
c       nref   : reference number for each element
c
c       nbd    : number of boundary nodes
c       ibdnod : pointer for boundary nodes
c       nic  : number of boundary sections
c       ic : index of nodes on boundary sections
c       nbound : number of closed boundaries
c       nside  : number of sections on each closed boundary
c
c       iout   : write unit number
c***********************************************************************
      integer nnode,nvert,nelem,inod(m,nelem),nref(nelem),iout
      integer nbd,ibdnod(nbd),nic,ic(nic+1),nbound,nside(nbound)
      integer nec(3,nelem)
      double precision x(nnode),y(nnode)
      character*8 filename
      integer i,k
c
      open(unit=iout,file=filename)
c
      write(iout,1010) nvert,nnode,nelem
      do i=1,nelem
        write(iout,1010) (inod(k,i),k=1,ish),nref(i),(nec(j,i),j=1,3)
      enddo
      do k=1,nnode
        write(iout,'(2(1pe13.6))') x(k),y(k)
      enddo
c
      write(iout,1010) nbd
      write(iout,1010) (ibdnod(k),k=1,nbd)
      write(iout,1010) nic
      write(iout,1010) (ic(k),k=1,nic+1)
      write(iout,1010) nbound
      write(iout,1010) (nside(k),k=1,nbound)
c
 1010 format(13i6)
c
      close(iout)
c
      return
      end
