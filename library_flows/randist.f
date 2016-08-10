cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine to generate random particle positions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine randist(np,l1,h1,part,diaa,gapp)
c
c     INPUT:
c       np = number of particles to be placed
c       l1 = length of the domain in x-direction
c       h1 = length of the domain in y-direction
c       diaa(np) = diameter of the particles to be placed
c       gapp = gap parameter
c
c     OUTPUT:
c       part(3,np) = position of the particles placed in the domain
c      
c     This subroutine generates random particle (circular of any 
c     specified diameter) positions in a rectangular domain with walls
c     perpendicular to the y direction. The domain can be periodic
c     or non-periodic in the x direction
c     
c     Provision for walls perpendicular to the x direction can be
c     made through the input parameter l1 or by suitably appending
c     the subroutine if_touch
c
c     Modified on: August 14th, 1997
c-----------------------------------------------------------------------
      double precision l1,h1,part(3,np),diaa(np),gapp
      double precision ran1
      integer np,increment,itr
      logical touch
c
      do j=1,3
        do i=1,np
          part(j,i) = 0.
        enddo
      enddo
c
      increment = 7
c
      do 2 i=1,np
	itr = 0
 3      itr = itr + 1
	increment = increment + 1
	if (itr.gt.500) then
	   write(*,*) 'no place to put new particle ',i
	   stop
        endif
c-----------------------------------------------------------------------
c       generating x position of particle i
c-----------------------------------------------------------------------
c
	part(1,i) = ran1(increment)*l1   
	increment = increment + 1
c
c-----------------------------------------------------------------------
c       generating y position of particle i
c-----------------------------------------------------------------------
c
	part(2,i) = ran1(increment)*h1
c
	call if_touch(np,l1,h1,i,part,diaa,gapp,touch)
	if (touch) goto 3
 2    continue
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine to check if the new particle 
c     touches the existing particles and the walls
c     presently provision is made only for walls perpendicular to
c     the y direction
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine if_touch(np,l1,h1,ipart,part,diaa,gapp,touch)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer np,ipart
      double precision part(3,np),diaa(np),sx,sy,l1,h1
      double precision gapp,dsqr,xcoor,x0,delt
      logical touch
c
      touch = .false.
c
c-----------------------------------------------------------------------
c     check if the particles are colliding with the side walls 
c     similar statement can be included for the walls perpendicular
c     to the x direction
c-----------------------------------------------------------------------
c
      if (part(2,ipart).le.((gapp+0.5)*diaa(ipart)).or.
     &    part(2,ipart).ge.(h1-(gapp+0.5)*diaa(ipart))) then
	 touch = .true.
	 return
      endif
c
c-----------------------------------------------------------------------
c     check if the particles collide
c-----------------------------------------------------------------------
c
      do i = 1,ipart-1
	x0 = part(1,ipart)
	sx = xcoor(x0,part(1,i))-x0
	sy = part(2,i) - part(2,ipart)
	dsqr = sx**2+sy**2
	delt = gapp*max(diaa(i),diaa(ipart))
	if (dsqr.le.((0.5*(diaa(i)+diaa(ipart))+delt)**2)) then
	   touch = .true.
	   return
        endif
      enddo
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     function to generate a real random number
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function ran1(iseed)
c
      integer iseed
      double precision aux1
c
      aux1 = float(iseed)
      aux1 = sqrt(aux1)*546423.0
      ran1 = (1+cos(aux1))/2.0
c       ran1 = cos(aux1)**2
c
      return
      end
c
c-------------------------------------------------------------------------
