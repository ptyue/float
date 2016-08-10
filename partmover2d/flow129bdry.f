      subroutine flow129bdry(iou,pwall,xd,yd,rd)
c***********************************************************************
c     This subroutine generates the boundary file for two staged pipe
c     input:
c        flow129data (data file)
c     output:
c        amphi.bdry (data file)
c        pwall      (geometry info)
c        yd,rd   (interface info)
c        xd      (the height of the transition, i.e., h1)
c     Dec 14, 2007, Pengtao Yue
c***********************************************************************
      implicit none
      integer iou
      real(8) pwall(5,1),xd,yd,rd
c
      character(1) title
      real(8)r1,r2,h1,h2

c
      open(iou,file='flw129data')
      read(iou,1001) title
      read(iou,1001) title
      read(iou,1001) title
      read(iou,1001) title
	read(iou,*)r1,r2,h1,h2
      read(iou,1001) title
	read(iou,*)yd,rd
      xd=h1
      close(iou)
c
      pwall(1,1)=r1
      pwall(2,1)=r2
      pwall(3,1)=h1
      pwall(4,1)=h2
c
      open(iou,file='amphi.bdry')
      write(iou,1010)6,6
      write(iou,1010)0.d0,-h1
      write(iou,1010)r1,-h1
      write(iou,1010)r1,0.d0
      write(iou,1010)r2,0.d0
      write(iou,1010)r2,h2
      write(iou,1010)0,h2

      write(iou,1020)'r',1,'b',1,2,0,1
      write(iou,1020)'r',1,'b',2,2,1,2
      write(iou,1020)'r',1,'b',3,2,2,3
      write(iou,1020)'r',1,'b',4,2,3,4
      write(iou,1020)'r',1,'b',5,2,4,5
      write(iou,1020)'r',1,'b',6,2,5,0
      close(iou)
c
1001  format(a1)
1010  format(2(1x,g12.5))
1020  format(1x,'polyline',2(1x,a4,1x,i4),5(1x,i4))
	end
