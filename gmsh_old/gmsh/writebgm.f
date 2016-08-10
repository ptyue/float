        subroutine writebgm(nVertices,nElements,ngmx,ncmx,
     &             x,y,inod,hsize,bgmfile)
        implicit none
        intent(in):: nVertices,nElements,ngmx,ncmx,x,y,inod,hsize
        integer nVertices, nElements,ngmx,ncmx
        real(8) x(nVertices),y(nVertices),hsize(nVertices)
        integer inod(ngmx,nElements) 
        character(100) bgmfile
c       local variables        
        integer iou
        data iou/21/
        integer i,j
        
        open(iou,file=bgmfile)
        write(iou,*)'View "background mesh" {' 
        do i=1,nElements
          write(iou,'("ST(")',advance='no')
          do j=1,ncmx
            if(j<ncmx)then
              write(iou,'(3(g12.5,","))',advance='no')
     &                      x(inod(j,i)),y(inod(j,i)),0.d0
            else
              write(iou,'(2(g12.5,","),g12.5,"){")',advance='no')
     &                      x(inod(j,i)),y(inod(j,i)),0.d0
            endif
          enddo
          do j=1,ncmx
            if(j<ncmx)then
              write(iou,'(g12.5,",")',advance='no')hsize(inod(j,i))
            else
              write(iou,'(g12.5,"};")')hsize(inod(j,i))
            endif
          enddo
        enddo
        write(iou,*)"};"
        close(iou)
        end
        
        
             
