c***********************************************************************
      subroutine UPDATE ( itime,itordr,dt,dt1,dtmin,dtmax,dtrate,
     &                    nnode,nelem,inod,x,y,area,aspr,ngmx,reft,
     &                    ncpg,xold,yold,umesh,lrmsh,
     &                    nbd,ibdnod,nic,ic,mshslpp,mshslpw,ifixpt,
     &                    nbrigid,nbfluid,nbelast, ncpb,xpos,xpod,
     &                     uvom,uvod1,fxym,amas,fpart,fbody,fcols,
     &                    ioptm,lfilm,maxitm,epsm,
     &                    diaa,diab,disinc,velinc,rinc,ncr,
     &                    icollision,epscollision,flowtype,pwall,bdseg,
     &                    nrdmsh,nrdgeom,nrdvelo,
     &                    iwork,maxiwk,rwork,maxrwk,iou)
c
c     This routine determines the appropriate time step and 
c     updates the position and condition of the particles.
c
c     collision models
c       icollision = 0 no model
c       icollision = 1 modify the particle position only
c       icollision = 2 modify the particle position and velocity,
c                      and the body force on the particle!
c
c     Note : fxym is changed on return!
c
c     By Howard Hu & Mingyu Zhu. July,14,1997
c***********************************************************************
      implicit none
      integer itime,itordr,ncpb,ncpg,nic,ic(nic+1),nbd,ibdnod(nbd)
      integer nnode,nelem,ngmx,inod(ngmx,nelem),icollision,flowtype
      integer nbrigid,nbfluid,nbelast,bdseg, reft(nelem),ncr
      double precision dt,dt1,dtmax,dtmin,dtrate,
     &                 x(nnode),y(nnode),xold(nnode),yold(nnode),
     &                 area(nelem),aspr(nelem),umesh(ncpg,nnode),
     &                 xpos(ncpb,nbrigid),xpod(ncpb,nbrigid),
     &                 uvom(ncpb,nbrigid),uvod1(ncpb,nbrigid),
     &                 fxym(ncpb,nbrigid),amas(ncpb,nbrigid),
     &                 fpart(ncpb,nbrigid),fbody(ncpb,nbrigid),
     &                 fcols(ncpb,nbrigid), 
     &                 diaa(nbrigid),diab(nbrigid),disinc,velinc,
     &                 rinc,pwall(5,bdseg),epscollision,epsm,smax
      logical mshslpw,mshslpp,lrmsh,ifixpt(ncpb)
      integer nrdmsh,nrdgeom,nrdvelo, ioptm,lfilm,maxitm
      integer maxiwk,maxrwk,maxmt,iou,iwork(maxiwk)
      double precision rwork(maxrwk)
c
c     local variables
      integer ndgmsh,ndggeom,ndgvelo
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      integer limp,liperm,ljperm,liwk,ljwk,lia,lja,liwork,
     &        larhs,lvec,la,lrwork,licont,leniwk,lpnew,
     &        lg,lfxym,lenrwk,lubdmsh,lvbdmsh,lshape,lamesh
      integer i,j,k, ia,ib,nd,k1,n0,n,nwall,ierr, mxiwk,mxrwk
      double precision eps,dis,dis1,ymax,ymin,csc,snc,tnc,csf,snf,r,
     &                 xcoor,xploc,tmp1,tmp2,xp,yp,xk,yk,
     &                 xkold,ykold,dceta,si,co,xmin,xmax
      logical oldmsh
      integer ncontmax,nelt, ndglb
      integer i1,i2,nel,m
      real*8 ds0,ds,ds1,xn(3),yn(3), x1,x3,y1,y3,a,b,xc,yc,r2
c
c     Save current positions and conditions
c     -------------------------------------
      do j=1,ncpb
        do i=1,nbrigid
          xpod(j,i) = xpos(j,i)
          uvod1(j,i) = uvom(j,i)
        enddo
      enddo
      do i=1,nnode
        xold(i) = x(i)
        yold(i) = y(i)
      enddo
      dt1 = dt
      if ( dtrate.le.1.d0 ) dtrate=1.5d0
      lrmsh = .false.
      nwall = nic-nbrigid-nbfluid-nbelast
c
c     Choose initial time step
c     ------------------------
      if ( itime.eq.1 ) then
         do k=1,nbrigid
           do i=1,ncpb
             fbody(i,k) = fpart(i,k)
             fcols(i,k) = 0.d0
             fxym (i,k) = fpart(i,k)
           enddo
         enddo
         goto 80
      endif
c
      dt = min(dtmax,dt1*dtrate)
      if ( ifixpt(1) .and. ifixpt(2) .and. ifixpt(3) ) return
c
      oldmsh = .false.
!      eps = 0.0001*diaa(1)

      smax = 1.d0
c
      ndgvelo = ndglb(nrdvelo)
      ndggeom = ndglb(nrdgeom)
      ndgmsh  = ndglb(nrdmsh)
c
      lshape  = 1
      limp    = lshape  + nelem
      liperm  = limp    + ndgmsh
      ljperm  = liperm  + ndgmsh
      liwk    = ljperm  + ndgmsh
      ljwk    = liwk    + ndgvelo
      liwork  = ljwk    + ndgmsh
c
      do i=1,nelem
         iwork(lshape+i-1) = 0
      enddo
      call matlen1 ( maxmt, nelem,inod,ngmx,iwork(lshape), nrdmsh,
     &              iwork(liwk),iwork(ljwk),iwork(liwork) )
      nelt  = maxmt
      mxiwk = max(3*nelem,maxmt+ndgmsh*(2*lfilm+5)+3 )
      mxrwk = maxmt + ndgmsh*(2*lfilm+7)+2
      ncontmax = 10*nbrigid
c
      lia     = ljwk    + ndgmsh
      lja     = lia     + ndgmsh+1
      liwork  = lja     + maxmt
      licont  = liwork  + mxiwk
      leniwk  = licont  + 3*ncontmax
      if ( leniwk.gt.maxiwk ) then
         write(*,*) 'iwork too short, should be larger than',leniwk
         stop
      endif
c
      lamesh  = 1
      lubdmsh = lamesh  + ncpg*ndgvelo
      lvbdmsh = lubdmsh + nbd
      larhs   = lvbdmsh + nbd
      lvec    = larhs   + ndgmsh
      la      = lvec    + ndgmsh
      lrwork  = la      + maxmt
      lpnew   = lrwork  + mxrwk
      lg      = lpnew   + 2*nbrigid
      lfxym   = lg      + 5*ncontmax
      lenrwk  = lfxym   + 3*nbrigid
      if ( lenrwk.gt.maxrwk ) then
         write(*,*) 'rwork too short, should be larger than',lenrwk
         stop
      endif
c
c     particle acceleration due to forces other than the collision force
c
      if ( icollision.le.1 ) then
        do i=1,nbrigid
           do j=1,ncpb
             fxym(j,i) = fxym(j,i)/amas(j,i)
           enddo
        enddo
      elseif ( icollision.eq.2 ) then
        do i=1,nbrigid
           do j=1,ncpb
             fxym(j,i) = (fxym(j,i)-fcols(i,j))/amas(j,i)
           enddo
        enddo
      endif
c
      call shap2d ( shg2,shg2x,shg2y,we2 )
c
c     Determine the appropriate time step
c     -----------------------------------
 40   continue
      tmp1 = 0.5d0*dt*dt
c
      if ( icollision.eq.1 .or. icollision.eq.2 )
     &   call collision (nbrigid,bdseg,epscollision,dt,diaa,amas,xpod,
     &                   uvod1,fxym,pwall,fcols,ncontmax,
     &                   rwork(lpnew),rwork(lg),iwork(licont) )
c
      tmp2 = 0.d0
      do i=1,nbrigid
         do j = 1,ncpb
           tmp2 = tmp2 + abs(fcols(j,i))
         enddo
      enddo
      if ( tmp2.gt.0 ) then
         print*,'particle collision occurs'
c         do i=1,nbrigid
c           write(*,*) 'fcols=',(fcols(j,i),j=1,2)
c         enddo
      endif
c
c     new particle accleration with collision force
c
      do i=1,nbrigid
         do j = 1,ncpb
            rwork(lfxym+(i-1)*ncpb+j-1) = fxym(j,i) + fcols(j,i)
         enddo
      enddo
c
      if ( icollision.ge.1 .or. itordr.ge.2 ) then
         do k=1,nbrigid
            do i=1,ncpb
               uvom(i,k) = uvod1(i,k) + dt*fxym(i,k)
               xpos(i,k) = xpod(i,k) +  dt*uvod1(i,k) +
     &                     tmp1*rwork(lfxym+i-1+(k-1)*ncpb)
            enddo
         enddo
      endif
      if ( icollision.eq.0 .and. itordr.eq.1 ) then
         do k=1,nbrigid
            do i=1,ncpb
               uvom(i,k) = uvod1(i,k)
               xpos(i,k) = xpod(i,k) + dt*uvod1(i,k)
            enddo
         enddo
      endif
      if ( itordr.eq.1 ) then
         do k=1,nbrigid
            do i=1,ncpb
               uvom(i,k) = uvod1(i,k)
            enddo
         enddo
      endif
c
c     particles move too much (> disinc)
c
      dis = 0.d0
      do k=1,nbrigid
         dis1=(xpos(1,k)-xpod(1,k))**2+(xpos(2,k)-xpod(2,k))**2
         dis = max(dis,dis1)
      enddo
      if(dis .gt. disinc**2) then
         dt = dt/dtrate
         go to 40
      endif
c
c     velocities change too much (>velinc)
c
      dis = 0.d0
      do k=1,nbrigid
         dis1=max(dabs(uvom(1,k)-uvod1(1,k)),
     &            dabs(uvom(2,k)-uvod1(2,k)),
     &            dabs(uvom(3,k)-uvod1(3,k))/10.)
         dis = max(dis,dis1)
      enddo
      if(dis .gt. velinc) then
         dt = dt/dtrate
         go to 40
      endif
c
c     particles rotate too much (>rinc)
c
      dis = 0.d0
      do k=1,nbrigid
         dis = max(dis,dabs(xpos(3,k)-xpod(3,k)))
      enddo
      if(dis .gt. rinc) then
         dt = dt/dtrate
         go to 40
      endif
c
c     for collision scheme No.2, modify the particle velocity
c     -------------------------------------------------------
      do k=1,nbrigid
         do i=1,ncpb
            fcols(i,k)  = amas(i,k)*fcols(i,k)
         enddo
      enddo
c
      if ( icollision.eq.2 ) then
         do k=1,nbrigid
           do i=1,ncpb
             uvom(i,k)  = uvod1(i,k) + dt*rwork(lfxym+i-1+(k-1)*ncpb)
             fbody(i,k) = fpart(i,k) + fcols(i,k)
           enddo
         enddo
      else
         do k=1,nbrigid
           do i=1,ncpb
             fbody(i,k) = fpart(i,k)
           enddo
         enddo
      endif
c
c     for fixed particles
c
      do i=1,ncpb
         if ( ifixpt(i) ) then
            do k=1,nbrigid
               uvom(i,k) = uvod1(i,k)
               xpos(i,k) = xpod(i,k)
            enddo
         endif
      enddo
c
c     update periodic position
c     ------------------------
      do k=1,nbrigid
         xpos(1,k) = xploc(xpos(1,k))
      enddo
c
c     calculate the mesh acceleration
c     -------------------------------
      call bdmacl (rwork(lubdmsh),rwork(lvbdmsh),mshslpp,mshslpw,
     &             rwork(lfxym),nbd,ncpb,nbrigid,nbfluid,nbelast,
     &             nic,ic,ibdnod,uvod1,xpod,x,y,nnode,ifixpt)
c
      do i=1,ncpg*ndgmsh
         rwork(lamesh+i-1) = 0.d0
      enddo
      call mshmov (oldmsh,nrdmsh,nrdgeom,ndgmsh,ndggeom,ndgvelo,
     &             iou,rwork(lamesh),ncpg,rwork(lubdmsh),rwork(lvbdmsh),
     &             nelem,inod,ngmx,x,y,iwork(lshape),nic,ic,nbd,ibdnod,
     &             shg2,shg2x,shg2y,we2,
     &             ioptm,lfilm,maxitm,epsm,smax,
     &             rwork(la),iwork(lja),nelt,iwork(lia),
     &             iwork(liperm),iwork(ljperm),iwork(limp),
     &             rwork(larhs),rwork(lvec),
     &             iwork(liwk),iwork(ljwk),iwork(liwork),mxiwk,
     &             rwork(lrwork),mxrwk )
c
c     update the mesh 
c     ---------------
      tmp1 = 0.5d0*dt*dt
      if ( icollision.ge.1 .or. itordr.ge.2 ) then
         do i=1,nnode
            x(i) = xploc( xold(i) + dt*umesh(1,i) + 
     &                    tmp1*rwork(lamesh+ncpg*(i-1)) )
            y(i) = yold(i) + dt*umesh(2,i) + 
     &                    tmp1*rwork(lamesh+ncpg*(i-1)+1)
         enddo
      endif
      if ( icollision.eq.0 .and. itordr.eq.1 ) then
         do i=1,nnode
            x(i) = xploc( xold(i) + dt*umesh(1,i) )
            y(i) =        yold(i) + dt*umesh(2,i)
         enddo
      endif
c
c     reset the nodes on particle surface - rigid body motion
c     -------------------------------------------------------
      if (.not.mshslpp) then
         do k=1,nbrigid
            xk = xpos(1,k)
            yk = xpos(2,k)
            xkold = xpod(1,k)
            ykold = xpod(2,k)
            dceta = xpos(3,k) - xpod(3,k)
            si = dsin(dceta)
            co = dcos(dceta)
            do i=ic(nwall+k),ic(nwall+k+1)-1
               nd = ibdnod(i)
               xp = xcoor(xkold,xold(nd)) - xkold
               yp = yold(nd) - ykold
               tmp1 = xk + xp*co - yp*si
               tmp2 = yk + yp*co + xp*si
c              write(33,'(i5,2e13.5)') i,tmp1-x(nd),tmp2-y(nd)
               x(nd) = xploc( tmp1 )
               y(nd) =  tmp2
            enddo
         enddo
      endif
c
c     Check mesh quality
c     ------------------
 80   call MSHCHK (lrmsh,nnode,nelem,inod,x,y,area,aspr,xold,yold,ierr)
      if ( ierr.eq.2 .or. (ierr.eq.1 .and. icollision.eq.0) ) then
         dt = dt/dtrate
         lrmsh = .true.
         goto 40
      endif
c
c     check the interfacial nodes
c     ---------------------------
      do k=1,nbfluid+nbelast
         ds0 = (2.d0*3.14159265358979d0/float(ncr))**2
         r2  = 0.25*diaa(k)**2
         i1 = ic(nwall+nbrigid+k)
         i2 = ic(nwall+nbrigid+k+1)
         nel = (i2-i1)/2
         do m=1,nel
            ib = i1 + 2*(m-1)
            n = ibdnod(ib)
            xn(1) = x(n)
            yn(1) = y(n)
            do i=2,3
              j = ib + i-1
              if ( j.eq.i2 ) j=i1
              n = ibdnod(j)
              xn(i) = xcoor(xn(1),x(n))
              yn(i) = y(n)
            enddo
c
            x1 = xn(1) - xn(2)
            x3 = xn(3) - xn(2)
            y1 = yn(1) - yn(2)
            y3 = yn(3) - yn(2)
            ds = 2*(x1*y3-x3*y1)
            if ( abs(ds).gt.1.d-20 ) then
               ds = 1.d0/ds
               a  = x1**2+y1**2
               b  = x3**2+y3**2
               xc = (a*y3-b*y1)*ds
               yc = (b*x1-a*x3)*ds
               a  = xc**2 + yc**2
               a = min(r2,a)
               a = max(0.1*r2,a)
            else
               a  = r2
            endif
c
            ds1 = (xn(2)-xn(1))**2 + (yn(2)-yn(1))**2
            ds  = (xn(3)-xn(1))**2 + (yn(3)-yn(1))**2
            if ( ds1 .lt. ds*(1.d0/3.d0)**2 .or.
     &           ds1 .gt. ds*(2.d0/3.d0)**2 ) lrmsh = .true.
            ds = ds/a
            if ( ds.ge.3.d0*ds0 .or. ds.le.0.25*ds0 )
     &         lrmsh = .true.
        enddo
      enddo
c
c     update the flow domain - pwall
c     ------------------------------
      if ( flowtype.eq.2 ) then
         xmin = x(1)
         xmax = x(1)
         do k=1,nwall
            do i=ic(k),ic(k+1)-1
               xmax = max(xmax,x(ibdnod(i)))
               xmin = min(xmin,x(ibdnod(i)))
            enddo
         enddo
         pwall(2,1) = xmin
         pwall(4,1) = xmin
         pwall(2,2) = xmin
         pwall(4,2) = xmax
         pwall(2,3) = xmax
         pwall(4,3) = xmax
         pwall(2,4) = xmax
         pwall(4,4) = xmin
      endif
c
      return
      end
c
c***********************************************************************
      subroutine MSHCHK (lrmsh,nnode,nelem,inod,x,y,area,aspr,xold,yold,
     &                   ierr )
c
c     This routine checks mesh distortion (element area & aspect ratio)
c
c     INPUT:
c       nelem,nnode,inod,x,y,area,aspr: information about the old mesh
c
c     OUTPUT:
c       lrmsh = .true.  remeshing is needed
c             = .false. remeshing is not needed
c***********************************************************************
      implicit none
      integer nelem,nnode, inod(6,nelem),ierr
      real*8 x(nnode),y(nnode),area(nelem),aspr(nelem),
     &       xold(nnode),yold(nnode)
      logical lrmsh
c
      integer i, n1,n2,n3
      real*8 fa,far,fam,farm,x1,x2,x3,y1,y2,y3,xcoor,
     &                 x21,x31,y21,y31,tmp,tmp1,dd
      data dd/1.3d0/
c
      ierr = 0
      fam  = 0.0
      farm = 0.0
      do i=1,nelem
        n1 = inod(1,i)
        n2 = inod(2,i)
        n3 = inod(3,i)
c
        x1 = xold(n1)
        x2 = xcoor(x1,xold(n2))
        x3 = xcoor(x1,xold(n3))
        y1 = yold(n1)
        y2 = yold(n2)
        y3 = yold(n3)
        tmp1 = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
c
        x1 = x(n1)
        x2 = xcoor(x1,x(n2))
        x3 = xcoor(x1,x(n3))
        y1 = y(n1)
        y2 = y(n2)
        y3 = y(n3)
        x21 = x2-x1
        x31 = x3-x1
        y21 = y2-y1
        y31 = y3-y1
        tmp = x21*y31 - y21*x31
c	write(21,'(i5,3e14.5)') i,tmp,area(i),tmp/area(i)
        if ( tmp.le.0.d0 ) then
          ierr = 2
          return
        endif
        if ( dabs(dlog(tmp/tmp1)).gt.dd ) then
          ierr = 1
        endif
c
        fa  = dlog(tmp/area(i))
        far = dlog(max(x21**2+y21**2,x31**2+y31**2,
     &                (x3-x2)**2+(y3-y2)**2) / (tmp*aspr(i)))
        fam  = max(fam,dabs(fa))
        farm = max(farm,dabs(far))
        if ( fam.gt.dd .or. farm.gt.dd ) then
          lrmsh = .true.
        endif
c
      enddo
c
      write(*,'(A,2f13.5)') 'mesh deformation: ', fam,farm
      return
      end
