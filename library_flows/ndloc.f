C***********************************************************************
      function ndloc ( ishape, iorder )
c
c     INPUT
c       ishape: shape of the element, =0 triangular/ =1 quadrilateral
c       iorder: order of interpolation, =0 / =1 / =2
c
c     OUTPUT
c       ndloc: number of local nodes 
c
      implicit none
      integer ndloc, ishape, iorder
c
      if ( iorder.eq.0 ) then
        ndloc = 1
      elseif ( iorder.eq.1 ) then
        ndloc = 3 + ishape
      elseif ( iorder.eq.2 ) then
        ndloc = 6 + 3*ishape
      endif
c
      return
      end
c
c***********************************************************************
      function ndglb ( iorder )
c
c     INPUT
c       iorder: order of interpolation, =0 / =1 / =2
c
c     OUTPUT
c       ndglb: number of the global nodes
c
      implicit none
      integer ndglb, iorder, nnode,nvert,nelem
      common /mesh/ nnode,nvert,nelem
c
      if ( iorder.eq.0 ) then
        ndglb = nelem
      elseif ( iorder.eq.1 ) then
        ndglb = nvert
      elseif ( iorder.eq.2 ) then
        ndglb = nnode
      endif
c
      return
      end

