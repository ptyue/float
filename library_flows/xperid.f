c***********************************************************************
      double precision function xcoor(x0,x)	
      double precision x0,x,xlength
      logical xperid
      common /xperid/ xlength,xperid
      xcoor = x
      if ( xperid ) then
         if ( (x-x0).lt.-0.5*xlength ) xcoor = x+xlength
         if ( (x-x0).gt. 0.5*xlength ) xcoor = x-xlength
      endif
      return
      end
c
c***********************************************************************
      double precision function xploc(x)
      double precision x,xlength
      logical xperid
      common /xperid/ xlength,xperid
      xploc = x
      if ( xperid ) then
        if ( x.lt.  0.0 ) xploc = x+xlength
        if ( x.gt.xlength ) xploc = x-xlength
      endif
      return
      end
