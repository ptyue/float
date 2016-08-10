c***********************************************************************
      subroutine INIBD (nbd,ibdu,ibdt,ibds,ibdp,densu,denst,denss,densp,
     &            ncpu,ncps)
c
c     This routine initializes the boundary conditions
c***********************************************************************
      implicit none
      integer nbd,ncpu,ncps,ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),ibdp(nbd)
      real*8 densu(ncpu,nbd),denst(nbd),denss(ncps,nbd),densp(nbd)
c
      integer i,j
c
      do i=1,nbd
c
        do j=1,ncpu
           ibdu(j,i)  = 0
           densu(j,i) = 0.d0
        enddo
c
        ibdt(i)  = 0
        denst(i) = 0.d0
c
        ibds(i) = 0
        do j=1,ncps
           denss(j,i) = 0.d0
        enddo
c
      enddo
      ibdp(1:nbd)=0
      densp(1:nbd)=0.d0  
c
      return
      end
