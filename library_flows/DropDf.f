c***********************************************************************
      subroutine DropDf(t,io,iout,ishm,kaxi,nvert,nelem,ngmx,ncmx,
     &                  inod,nec,nrdgeom,nrdphi,ndggeom,ndgphi,x,y,phi,
     &                  iwork,rwork,dropone)

c
c     This subroutine calculate the deformation parameter of drop and
c        write the information into file deform.dat:
c        t,D,L,B,thmax,thmin,x0,y0
c     Input:
c        t:       time
c        io:      the file is opened and closed within this subroutine
c                 (input as io5)
c        iout:    i/o unit for error messages(input as io6)
c        ishm:    0 for triangular cells and 1 for quadrilateral cells
c        kaxi:    0 for 2D planar; 1 for 2D axisymmetric
c        dropone: T, phi=1 in drop; F, phi=-1 in drop
c     Output:      
c     Working:
c        iwork,rwork
c     Written Feb 23, 2005 by Pengtao Yue
c***********************************************************************
      implicit none
      integer io,iout,ishm,kaxi,nvert,nelem,ngmx,ncmx,
     &        nrdgeom,nrdphi,ndggeom,ndgphi
      integer inod(ngmx,nelem),nec(ncmx,nelem),iwork(nelem)
      real*8  x(ndggeom),y(ndggeom),phi(ndgphi),rwork(nelem*2+nvert),t
      logical dropone
c
      real*8 D,L,B,thmax,thmin,x0,y0,f0,rr,x1,y1,pi
      integer npt,i,Imax,Imin
      integer lrintr,lrphi
c
      pi=acos(-1.d0)
c
      lrintr=1
      lrphi=lrintr+2*nelem    
c     get the phi=f0 contour
      f0=0.d0
      call getintrfc(iout,ishm,nvert,nelem,ngmx,inod,ncmx,nec,
     &               x,y,phi,f0,npt,rwork(lrintr),iwork,rwork(lrphi))
c     get the center of the enclosed area
      if(kaxi.eq.0)then
         call getCntrGeom(iout,x0,y0,npt,rwork(lrintr))
      else
         call getCntr(iout,ishm,kaxi,nelem,ngmx,inod,
     &                nrdgeom,nrdphi,ndggeom,ndgphi,
     &                x,y,phi,dropone,x0,y0)
      endif
      rr=(rwork(lrintr)-x0)**2+(rwork(lrintr+1)-y0)**2
      L=rr
      B=rr
      Imax=1
      Imin=1
      do i=2,npt
         rr=(rwork(lrintr+(i*2)-2)-x0)**2+(rwork(lrintr+(i*2)-1)-y0)**2
         if(rr.gt.L)then
            Imax=i
            L=rr
         endif
         if(rr.lt.B)then
            Imin=i
            B=rr
         endif
      enddo
      L=sqrt(L)
      B=sqrt(B)
c     polar angle of the long axis
      x1=rwork(lrintr+(Imax*2)-2)-x0
      y1=rwork(lrintr+(Imax*2)-1)-y0
      if(y1.ge.0.d0)then
         thmax=acos(x1/L)
      else
         thmax=pi-acos(x1/L)
      endif
c     polar angle of the short axis
      x1=rwork(lrintr+(Imin*2)-2)-x0
      y1=rwork(lrintr+(Imin*2)-1)-y0
      if(y1.ge.0.d0)then
         thmin=acos(x1/L)
      else
         thmin=pi-acos(x1/L)
      endif
c
      D=(L-B)/(L+B)
c
      open(io,file='deform.dat',access='append')
      write(io,100)t,D,L,B,thmax,thmin,x0,y0
100   format(8(1x,g12.5))
      close(io)
      end
c***********************************************************************
      subroutine getCntr(iout,ishm,kaxi,nelem,ngmx,inod,
     &                   nrdgeom,nrdphi,ndggeom,ndgphi,
     &                   x,y,phi,dropone,x0,y0)
c     the subroutine calculate the center of the drop
c     Input:
c        kaxi:    0, planar;1, axisymmetric
c        ishm:    0, triangular element; 1, quadrilateral element
c        dropone:    T, phi=1 in drop phase; F, phi=-1 in drop phase
c     Output:
c        x0,y0:      coordinates of the drop center
c     written Feb 22, 2005, by Pengtao Yue
c***********************************************************************
      implicit none
      integer ishm,kaxi,nelem,ngmx,ndggeom,ndgphi,nrdgeom,
     &        nrdphi,iout
      integer inod(ngmx,nelem)
      logical dropone
      real*8  x(ndggeom),y(ndggeom),phi(ndgphi),x0,y0
c
      real*8  shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      integer ndloc,ndlgeom,ndlphi
      integer ne,l,n1,n2,i
      real*8  xcoor,fn(9),xn(9),yn(9),fl,xx,yy,dxxi,dxet,dyxi,
     &        dyet,area,jcbn,aljcwe,x1
c
      if(nrdgeom.gt.2.or.nrdphi.gt.2)then
         write(iout,*)'nrdgeom or nrdphi out of range in getCntr'
         stop
      endif
      ndlgeom=ndloc(ishm,nrdgeom)
      ndlphi=ndloc(ishm,nrdphi)         
      n1=ishm*2+nrdgeom
      n2=ishm*2+nrdphi
c     To calculate the center of the drop
      call shap2d ( shg2,shg2x,shg2y,we2)
      x0=0.d0
      y0=0.d0
      area=0.d0
      do ne=1,nelem
         if(dropone)then
            do i=1,ndlphi
               fn(i)=(1.+phi(inod(i,ne)))/2.
	      enddo
         else
            do i=1,ndlphi
               fn(i)=(1.-phi(inod(i,ne)))/2.
	      enddo
         endif
         x1=x(inod(1,ne))
         do i=1,ndlgeom
            xn(i)=xcoor(x1,x(inod(i,ne)))
            yn(i)=y(inod(i,ne))
         enddo
         do l=1,ishm*2+7
c     calculate values at the quadrature points
            xx   = sum(xn(1:ndlgeom)*shg2 (1:ndlgeom,l,n1))
            yy   = sum(yn(1:ndlgeom)*shg2 (1:ndlgeom,l,n1))
            dxxi = sum(xn(1:ndlgeom)*shg2x(1:ndlgeom,l,n1))
            dxet = sum(xn(1:ndlgeom)*shg2y(1:ndlgeom,l,n1))
            dyxi = sum(yn(1:ndlgeom)*shg2x(1:ndlgeom,l,n1))
            dyet = sum(yn(1:ndlgeom)*shg2y(1:ndlgeom,l,n1))
            fl   = sum(fn(1:ndlphi) *shg2 (1:ndlphi,l,n2))
c     Jacobian
            jcbn  = dxxi*dyet - dxet*dyxi            
	      aljcwe = (1.d0 + kaxi* (xx - 1.d0))*jcbn*we2(l,n1)
c
            x0  =x0  +xx*fl*aljcwe
            y0  =y0  +yy*fl*aljcwe
            area=area+   fl*aljcwe
         enddo
      enddo
      x0=x0/area
      y0=y0/area
      if(kaxi.ne.0)x0=0.d0
      end
c***********************************************************************
      subroutine getCntrGeom(iout,x0,y0,npt,intr)
c     the subroutine calculate the center of the drop in a geometric way
c     Only the line segements which encircles the drop are used. This
c     subroutine only works for 2D planar.   
c     Input:
c        intr:    coordinates of the points along the interface
c        npt:     number of interface points
c     Output:
c        x0,y0:   coordinates of the drop center
c     Jun 02, 2005, by Pengtao Yue
c***********************************************************************
      implicit none
      intent(in)  :: iout,npt,intr
      intent(out) :: x0,y0
      integer npt,iout
      real(8) intr(2,npt),x0,y0
c
      real(8) x1,y1,x2,y2,x3,y3,s,area
      integer i
c
      if(npt.le.2)then
         write(iout,*)'interface info not correct in getCntrGeom'
         stop
      endif
c
c     planar 2D case
c     --------------
      x0=0.d0
      y0=0.d0
      area=0.d0
      x1=intr(1,npt)
      y1=intr(2,npt)
      do i=1,npt-2
         x2=intr(1,i)
         y2=intr(2,i)
         x3=intr(1,i+1)
         y3=intr(2,i+1)
         s=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
         area=area+s
         x0=x0+(x1+x2+x3)/3.d0*s
         y0=y0+(y1+y2+y3)/3.d0*s
      enddo
      x0=x0/area
      y0=y0/area
      end

c***********************************************************************
      subroutine getintrfc(iout,ishm,nvert,nelem,ngmx,inod,ncmx,nec,
     &                     x,y,phi,f0,npt,intr,iflag,f)
c     the subroutine calculate the phi=f0 level set and returns a 1d 
c     array representing it. The calculation is based on P1 elements
c     Input:
c     Output:
c        npt:        number of points on the contour
c        intr(2,npt):coordinates of the points
c     working:
c        iflag(nelem),f(nvert)
c     written Feb 22, 2005, by Pengtao Yue
c     Note: only works for the case that interface doesn't intersect the
c           boundary. 
c***********************************************************************
      implicit none
      integer iout,ishm,nvert,nelem,ngmx,ncmx,npt
      integer inod(ngmx,nelem),nec(ncmx,nelem),iflag(nelem)
      real*8  x(nvert),y(nvert),phi(nvert),f(nvert),intr(2,nelem)
      real*8  f0
c
      real*8  eps,fmax,fmin,f1,f2,d1,d2,fn(4),xn(4),yn(4)
      integer i,ne,ndl,n,i2,n1
      real*8  xcoor,x0
c
      if(ishm.eq.0)then
         ndl=3
      elseif(ishm.eq.1)then
         ndl=4
      else
         write(iout,*)'ishm wrong value in getintrfc'
         stop
      endif
c     calculating the numerical tolerance of phi
      eps=1.d-8
      fmax=0.
	do i=1,nvert
         fmax=max(fmax,abs(phi(nvert)))
      enddo
      eps=eps*max(fmax,1.d0)
c     shift f(phi) a little bit the avoid the case that the contour goes
c        through any vertices
      f(1:nvert)=phi(1:nvert)            
      do i=1,nvert
         if(abs(f(i)-f0).lt.eps)then
            if(f(i).gt.f0)then
               f(i)=f0+eps
            else
               f(i)=f0-eps
            endif
         endif
      enddo
c
      iflag(1:nelem)=0
c     only returns one connected contour
      npt=0
      do n=1,nelem
         ne=n
c     if the cell has been flagged, then look for the next cell
c     (only meaningful when calculating more than one contours)
         if(iflag(ne).eq.1) cycle
         iflag(ne)=1
         do i=1,ndl
            fn(i)=f(inod(i,ne))
         enddo
         fmin=fn(1)
         fmax=fmin
         do i=2,ndl
            fmin=min(fmin,fn(i))
            fmax=max(fmax,fn(i))
         enddo
         if(f0.gt.fmax.or.f0.lt.fmin)cycle
         n1=ne
100      continue  
         x0=x(inod(1,ne))       
         do i=1,ndl
            xn(i)=xcoor(x0,x(inod(i,ne)))
            yn(i)=y(inod(i,ne))
         enddo
c
         do i=1,ndl
c     doesn't allow the line to go back
            if(n1.eq.nec(i,ne))cycle
            i2=mod(i,ndl)+1
            f1=fn(i)
            f2=fn(i2)
            fmin=min(f1,f2)
            fmax=max(f1,f2)
            if(f0.gt.fmax.or.f0.lt.fmin)cycle
            d2=(f1-f0)/(f1-f2)
            d1=(f0-f2)/(f1-f2)
            npt=npt+1
            intr(1,npt)=d1*xn(i)+d2*xn(i2)
            intr(2,npt)=d1*yn(i)+d2*yn(i2)
            n1=ne
            ne=nec(i,n1)
            exit
         enddo
         if(ne.le.0)then
c     meet the boundary ( if this happens, this subroutine may not give 
c     the complete contour."exit" should be delted if all the contours 
c     are to be calculated)
            exit
         elseif(iflag(ne).eq.0)then
c     process next point on the contour
            iflag(ne)=1
            do i=1,ndl
               fn(i)=f(inod(i,ne))
            enddo
            goto 100
         else
c     The contour is closed
            exit
         endif
      enddo
c      open(15,file='intr.dat')
c      do i=1,npt
c      write(15,*)intr(1,i),intr(2,i)
c      enddo
c      close(15)
      return
      end
c***********************************************************************
      subroutine getintrfcbd(iout,ishm,nvert,nelem,ngmx,inod,ncmx,nec,
     &                     x,y,phi,f0,npt,intr,iflag,f)
c     the subroutine calculate the phi=f0 level set and returns a 1d 
c     array representing it. The array starts and ends at boundary.
c
c     The calculation is based on P1 elements
c
c     Input:
c     Output:
c        npt:        number of points on the contour
c        intr(2,npt):coordinates of the points
c     working:
c        iflag(nelem),f(nvert)
c     written Mar 7, 2009, Pengtao Yue ( modified from sub getintrfc)
c                
c     requires that nec(i,ne)<=0 if the i-the edge of ne-th element is
c     boundary
c
c     Note: works for level sets intersecting boundaries. For now, it 
c     only outputs one connected level curve. A small modification can 
c     output all the level curves. 
c***********************************************************************
      implicit none
      integer iout,ishm,nvert,nelem,ngmx,ncmx,npt
      integer inod(ngmx,nelem),nec(ncmx,nelem),iflag(nelem)
      real*8  x(nvert),y(nvert),phi(nvert),f(nvert),intr(2,nelem)
      real*8  f0
c
      real*8  eps,fmax,fmin,f1,f2,d1,d2,fn(4),xn(4),yn(4)
      integer i,ne,ndl,n,i2,n1
      real*8  xcoor,x0
c
      if(ishm.eq.0)then
         ndl=3
      elseif(ishm.eq.1)then
         ndl=4
      else
         write(iout,*)'ishm wrong value in getintrfc'
         stop
      endif
c     calculating the numerical tolerance of phi
      eps=1.d-8
      fmax=0.
	do i=1,nvert
         fmax=max(fmax,abs(phi(nvert)))
      enddo
      eps=eps*max(fmax,1.d0)
c     shift f(phi) a little bit the avoid the case that the contour goes
c        through any vertices
      f(1:nvert)=phi(1:nvert)            
      do i=1,nvert
         if(abs(f(i)-f0).lt.eps)then
            if(f(i).gt.f0)then
               f(i)=f0+eps
            else
               f(i)=f0-eps
            endif
         endif
      enddo
c
      iflag(1:nelem)=0
c     only returns one connected contour
      npt=0
      do n=1,nelem
         ne=n
c
c     getting the starting point of a level curve
c
c     if the cell has been flagged, then look for the next cell
c     (only meaningful when calculating more than one contours)
         if(iflag(ne).eq.1) cycle
         do i=1,ndl
            fn(i)=f(inod(i,ne))
         enddo
         fmin=fn(1)
         fmax=fmin
         do i=2,ndl
            fmin=min(fmin,fn(i))
            fmax=max(fmax,fn(i))
         enddo
         if(f0.gt.fmax.or.f0.lt.fmin)then
            iflag(ne)=1
            cycle
         endif
         if(nec(1,ne).ne.0.and.nec(2,ne).ne.0.and.nec(3,ne).ne.0)cycle
c
         x0=x(inod(1,ne))       
         do i=1,ndl
            xn(i)=xcoor(x0,x(inod(i,ne)))
            yn(i)=y(inod(i,ne))
         enddo
c
         n1=1
         do i=1,ndl
            i2=mod(i,ndl)+1
            f1=fn(i)
            f2=fn(i2)
            fmin=min(f1,f2)
            fmax=max(f1,f2)
            if(f0.gt.fmax.or.f0.lt.fmin)cycle
            n1=nec(i,ne)
            if(n1<=0)then
               npt=npt+1
               d2=(f1-f0)/(f1-f2)
               d1=(f0-f2)/(f1-f2)
               intr(1,npt)=d1*xn(i)+d2*xn(i2)
               intr(2,npt)=d1*yn(i)+d2*yn(i2)
               iflag(ne)=1
               exit
            endif
         enddo
c     the level set doesn't interset the boundary
         if(n1>0)cycle
c
         n1=ne
         do i=1,ndl
            i2=mod(i,ndl)+1
            f1=fn(i)
            f2=fn(i2)
            fmin=min(f1,f2)
            fmax=max(f1,f2)
            if(f0.gt.fmax.or.f0.lt.fmin)cycle
            ne=nec(i,n1)
            if(ne>0)then
               npt=npt+1
               d2=(f1-f0)/(f1-f2)
               d1=(f0-f2)/(f1-f2)
               intr(1,npt)=d1*xn(i)+d2*xn(i2)
               intr(2,npt)=d1*yn(i)+d2*yn(i2)
               iflag(ne)=1
               exit
            endif
         enddo
c         if(iflag(ne)==0)then
c         iflag(ne)=1
c         else
c     if this happens, there is something wrong with the phi field.
c     Probably some numerical noise on the boundary.
c            exit
c         endif

c         n1=ne
100      continue  
         do i=1,ndl
            fn(i)=f(inod(i,ne))
         enddo
         x0=x(inod(1,ne))       
         do i=1,ndl
            xn(i)=xcoor(x0,x(inod(i,ne)))
            yn(i)=y(inod(i,ne))
         enddo
c
         do i=1,ndl
c     doesn't allow the line to go back
            if(n1.eq.nec(i,ne))cycle
            i2=mod(i,ndl)+1
            f1=fn(i)
            f2=fn(i2)
            fmin=min(f1,f2)
            fmax=max(f1,f2)
            if(f0.gt.fmax.or.f0.lt.fmin)cycle
            d2=(f1-f0)/(f1-f2)
            d1=(f0-f2)/(f1-f2)
            npt=npt+1
            intr(1,npt)=d1*xn(i)+d2*xn(i2)
            intr(2,npt)=d1*yn(i)+d2*yn(i2)
            n1=ne
            ne=nec(i,n1)
            exit
         enddo
         if(ne<=0)then
c     meet the boundary ( the expected ending condition)
            exit
         elseif(iflag(ne).eq.0)then
c     process next point on the contour
            iflag(ne)=1
            goto 100
         else
c     The contour is closed
            exit
         endif
      enddo
c      open(15,file='intr.dat')
c      do i=1,npt
c      write(15,*)intr(1,i),intr(2,i)
c      enddo
c      close(15)
      return
      end
c***********************************************************************
      subroutine DropDfExt(t,io,nrdgeom,ndggeom,x,y,nrdphi,ndgphi,phi,
     &                     nic,ic,nbd,ibdnod)
c
c     calcaute drop deformation ONLY for extensional flow
c     symmetry conditions are used, the drop is centered at (0,0)
c        deform.dat:  t,D,L,B
c     Input:
c        t:       time
c        io:      the file is opened and closed within this subroutine
c                 (input as io5)
c     Output:      
c     Note:  applicable range: 2>=nrdgeom>=nrdphi
c            center of the drop is (0,0)
c
c     Written Apr 15, 2005 by Pengtao Yue
c***********************************************************************
      implicit none
      integer io,nrdgeom,nrdphi,ndggeom,ndgphi,nic,nbd
      integer ic(nic+1),ibdnod(nbd)
      real*8  x(ndggeom),y(ndggeom),phi(ndgphi),t
c
      real(8) f0,xr,yr,L,B,D,f1,f2,c1,c2
      integer i,n1,n2,k
c
      f0=0.d0
      xr=1.
      yr=1.
      if(nrdgeom.eq.2.and.nrdphi.eq.1)then
         k=2
      else
         k=1
      endif
c     along x axis
      do i=ic(1),ic(2)-k-1,k
         n1=ibdnod(i)
         n2=ibdnod(i+k)
         f1=phi(n1)
         f2=phi(n2)
         if(min(f1,f2).gt.f0.or.max(f1,f2).lt.f0)cycle
         if(f1.eq.f0)then
            xr=x(n1)
            exit
         elseif(f2.eq.f0)then
            xr=x(n2)
            exit
         else
            c1=(f2-f0)/(f2-f1)
            c2=(f0-f1)/(f2-f1)
            xr=c1*x(n1)+c2*x(n2)
            exit
         endif
      enddo
c     along y axis
      do i=ic(4),ic(5)-k-1,k
         n1=ibdnod(i)
         n2=ibdnod(i+k)
         f1=phi(n1)
         f2=phi(n2)
         if(min(f1,f2).gt.f0.or.max(f1,f2).lt.f0)cycle
         if(f1.eq.f0)then
            yr=y(n1)
            exit
         elseif(f2.eq.f0)then
            yr=y(n2)
            exit
         else
            c1=(f2-f0)/(f2-f1)
            c2=(f0-f1)/(f2-f1)
            yr=c1*y(n1)+c2*y(n2)
            exit
         endif
      enddo
      L=max(xr,yr)
      B=min(xr,yr)
      D=(L-B)/(L+B)
      open(io,file='deform.dat',access='append')
      write(io,100)t,D,L,B
100   format(8(1x,g12.5))
      close(io)
      end
c***********************************************************************
      subroutine DropLenWid(t,io,kaxi,ngmx,nelem,ndggeom,
     &                      inod,x,y,nrdphi,ndgphi,phi)
c
c     calcaute drop drop length and width
c        contents of deform.dat:  
c           t,xmax-xmin,ymax-ymin,(xmin+xmax)/2.,(ymin+ymax)/2.
c     Input:
c        t:       time
c        io:      the file is opened and closed within this subroutine
c                 (input as io5)
c     Output:      
c
c     Jun 13, 2005  Pengtao Yue
c        Note: Only works for 
c              1.non-periodic mesh, and triangular elements
c              2. 2>=nrdgeom>=nrdphi
c***********************************************************************
      implicit none
      intent(in) :: t,io,kaxi,ngmx,nelem,ndggeom,inod,x,y,nrdphi,ndgphi,
     &              phi
      integer io,ngmx,nelem,ndggeom,nrdphi,ndgphi,inod(ngmx,nelem),kaxi
      real*8  x(ndggeom),y(ndggeom),phi(ndgphi),t
c
      real(8)  xmin,xmax,ymin,ymax,f0,xx,yy,c1,c2,f1,f2
      integer  ndl,ne,i1,i2,k,edge(2,6)
c
c     phi=f0 contour is used
      f0=0.
c
      xmin=1.d10
      xmax=-xmin
      ymin=1.d10
      ymax=-ymin
      if(nrdphi.eq.1)then
         edge(1,1)=1
         edge(2,1)=2
         edge(1,2)=2
         edge(2,2)=3
         edge(1,3)=3
         edge(2,3)=1
         ndl=3
      elseif(nrdphi.eq.2)then
         edge(1,1)=1
         edge(2,1)=4
         edge(1,2)=4
         edge(2,2)=2
         edge(1,3)=2
         edge(2,3)=5
         edge(1,4)=5
         edge(2,4)=3
         edge(1,5)=3
         edge(2,5)=6
         edge(1,6)=6
         edge(2,6)=1
         ndl=6
      else
         stop 'error in nrdphi, subroutine DropLenWid'
      endif
         
      do ne=1,nelem      
      do k=1,ndl
         i1=inod(edge(1,k),ne)
         i2=inod(edge(2,k),ne)
         f1=phi(i1)
         f2=phi(i2)
         if(min(f1,f2).gt.f0.or.max(f1,f2).lt.f0)cycle
         if(f1.eq.f0)then
            xx=x(i1)
            yy=y(i1)
         elseif(f2.eq.f0)then
            xx=x(i1)
            yy=y(i2)
         else
            c1=(f2-f0)/(f2-f1)
            c2=(f0-f1)/(f2-f1)
            xx=c1*x(i1)+c2*x(i2)
            yy=c1*y(i1)+c2*y(i2)
         endif
         xmin=min(xmin,xx)
         xmax=max(xmax,xx)
         ymin=min(ymin,yy)
         ymax=max(ymax,yy)
      enddo
      enddo
      if(kaxi.eq.1)xmin=-xmax
      open(io,file='deform.dat',access='append')
      write(io,100)t,xmax-xmin,ymax-ymin,(xmin+xmax)/2.,(ymin+ymax)/2.
100   format(6(1x,g12.5))
      close(io)
      end     
 