c***********************************************************************
      subroutine timestep(itime,ishm,dt,dt1,dtmin,dtmax,dtrate,
     &                    cflflow,cflphi,
     &                    nelem,ngmx,ncpg,inod,
     &                    nrdgeom,nrdvelo,nrdphi,ndggeom,ndgvelo,ndgphi,
     &                    x,y,u,phi)
c-----------------------------------------------------------------------
c     Calcuate the time step according to 3 conditions:
c        1. time step growh rate
c        2. CFL condition due to flow speed(cflfl)
c        3. CFL condition due to interface moving speed(cflph)
c     Input:
c        itime:   the time step(1 for the 1st step)
c        ishm:    0 for triangular cells and 1 for quadrilateral cells
c        dt:      time step in previous step, return as the new value
c        dtmin,dtmax,dtrate:
c                 minimum,maximum, and growth rate of time step
c        cflflow: courant number for flow (can be large, such as 5)
c        cflphi:  cournat number for interface movement(<=1, otherwise
c                 the interface may move more than one cell in one step)
c        xxxxx
c     Output:
c        dt:      new time step
c        dt1:     previous time step
c
c-----------------------------------------------------------------------
c     Mar 18, 2005, Pengtao Yue
c-----------------------------------------------------------------------
      implicit none
      intent    (in):: itime,ishm,dtmin,dtmax,dtrate,cflflow,cflphi,
     &                 nrdgeom,nrdvelo,nrdphi,ndggeom,ndgvelo,ndgphi,
     &                 u,x,y,phi,ncpg,ngmx
      intent   (out):: dt1
      intent (inout):: dt
      integer itime,ishm,nelem,ngmx,ncpg,nrdgeom,nrdvelo,nrdphi,ndggeom,
     &        ndgvelo,ndgphi
      integer inod(ngmx,nelem)
      real(8) dt,dt1,dtmin,dtmax,dtrate,cflflow,cflphi
      real(8) x(ndggeom),y(ndggeom),u(ncpg,ndgvelo),phi(ndgphi)
c
c      logical, save :: first=.true.
      real(8), parameter :: eps=1.d-10, phic=0.5d0
      real(8) shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      real(8) ul,vl,fl,fxl,fyl,dxxi,dxet,dyxi,dyet,dxix,dxiy,detx,dety,
     &        jcbn,f1,f2,nx,ny,len,vel
      real(8) xcoor,xn(9),yn(9),x0
      integer l,n1,n2,n3,i,k,ne
      integer ndlgeom,ndlvelo,ndlphi,ndloc
      real(8) dtflow,dtphi
c   
      dt1=dt
c
c     calculate the time step according to growth rate
c     ------------------------------------------------
      if(itime.eq.1)then
c         first=.false.
         dt=dtmin
      else
         dt=min(dt*dtrate,dtmax)
      endif
c
c     calculate time step due to CFL contraints
c     -----------------------------------------
c     get the shape functions and their derivatives      
      call shap2d ( shg2,shg2x,shg2y,we2 )
c
      n1 = ishm*2 + nrdgeom
      n2 = ishm*2 + nrdvelo
      n3 = ishm*2 + nrdphi
c     the operations are done on the centroid of cell, i.e., l=1
      l=1
c
      ndlgeom = ndloc(ishm, nrdgeom)
      ndlvelo = ndloc(ishm, nrdvelo)
      ndlphi  = ndloc(ishm, nrdphi )
c
      dtflow=0.d0
      dtphi=0.d0
      do ne=1,nelem
c     mesh derivatives
         x0 = x(inod(1,ne))
         do i=1,ndlgeom
            k = inod(i,ne)
            xn(i) = xcoor(x0,x(k))
            yn(i) = y(k)
         enddo
         dxxi = sum(xn(1:ndlgeom)*shg2x(1:ndlgeom,l,n1))
         dxet = sum(xn(1:ndlgeom)*shg2y(1:ndlgeom,l,n1))
         dyxi = sum(yn(1:ndlgeom)*shg2x(1:ndlgeom,l,n1))
         dyet = sum(yn(1:ndlgeom)*shg2y(1:ndlgeom,l,n1))
         jcbn = dxxi*dyet - dxet*dyxi
         dxix =  dyet/jcbn 
         dxiy = -dxet/jcbn
         detx = -dyxi/jcbn
         dety =  dxxi/jcbn
c
c     CFL contraint due to fluid velocity
c     -----------------------------------
c        velocity
         ul = sum(u(1,inod(1:ndlvelo,ne))*shg2(1:ndlvelo,l,n2))
         vl = sum(u(2,inod(1:ndlvelo,ne))*shg2(1:ndlvelo,l,n2))
         vel= sqrt(ul**2+vl**2)+eps
         nx = ul/vel
         ny = vl/vel
c        local length along the flow directoin
         f1 = nx*dxix + ny*dxiy
         f2 = nx*detx + ny*dety
c         len = 1.d0 / (sqrt( f1*f1 + f2*f2 )+eps)
c         dtflow=max(dtflow,vel/len)
         len = sqrt( f1*f1 + f2*f2 )
         dtflow=max(dtflow,vel*len)
c        
c     CFL contraint due to interface movement
c     ---------------------------------------
         fl = sum(phi(inod(1:ndlphi,ne))* shg2(1:ndlphi,l,n3))
         if(abs(fl).lt.phic)then
c     only the interfacial cells will be considered
            f1 = sum(phi(inod(1:ndlphi,ne))*shg2x(1:ndlphi,l,n3))
            f2 = sum(phi(inod(1:ndlphi,ne))*shg2y(1:ndlphi,l,n3))
            fxl= f1*dxix+f2*detx
            fyl= f1*dxiy+f2*dety
            fl = sqrt(fxl**2+fyl**2)+eps
            nx = fxl/fl
            ny = fyl/fl
c     velocity along the interface normal
            vel= abs(ul*nx+vl*ny)
c     local length along the interface normal
            f1 = nx*dxix + ny*dxiy
            f2 = nx*detx + ny*dety
c            len = 1.d0 / (sqrt( f1*f1 + f2*f2 )+eps)
c            dtphi=max(dtphi,abs(vel)/len)
            len = sqrt( f1*f1 + f2*f2 )
            dtphi=max(dtphi,vel*len)
         endif
      enddo
      dtflow=cflflow/(dtflow+eps)
      dtphi =cflphi /(dtphi +eps)
c
      if(dt.gt.dtflow)then
         write(*,'("CFL condition(flow) enforced!, dtflow=",g12.5)')
     &         dtflow
         dt=dtflow
      endif
      if(dt.gt.dtphi)then
         write(*,'("CFL condition(interface) enforced!, dtphi=",g12.5)')
     &         dtphi
         dt=dtphi
      endif
      return
c-----------------------------------------------------------------------
      end
