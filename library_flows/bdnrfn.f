c***********************************************************************
      subroutine BDNRFN (refine,coor,h,bdsec,seclst,nbd,nbound,nside,
     &                  xymid,mbdlst,pwall,bdseg,nbrigid,diaa,diab,xpos,
     &                   maxnbd,link,rwork,ncb,iscb)
c
c-----------------------------------------------------------------------
c    This subroutine refines grids on bounary sections.
c
c    INPUT:
c         refine  = .true. refines the segments on boundary sections
c                 = .false. uses uniform segments on boundary sections
c         coor    = coordinates of the nodes
c         bdsec  = number of sections within the bounrary nodes
c         seclst = location of the first node in each section
c         nbd    = total number of boundary nodes
c         nbound = number of the closed boundary
c         nside  = number of sides for each closed boundary
c         nbrigid  = number of particles
c         diaa,diab = diameters for each particle
c         xpos   = postion for each particle
c         maxnbd = maximum number of boundary nodes
c         xymid = coordinates for the middle nodes
c         mbdlst = information on the boundary sections
c         pwall = information on the outer boundary sections
c         rwork(3*maxnbd), link(maxnbd), ncb(nbound),iscb(nbound+1)
c                = working arrays
c    OUTPUT:
c         coor = the refined the x,y coordinates of the nodes
c         h      = the distance between the neighbouring nodes
c
c    Mingyu Zhu & Howard Hu,  Sept. 28 , 1996/Jan. 1999
c***********************************************************************
      implicit none
      logical refine
      integer nbd,nbound,nside(nbound),bdsec,seclst(bdsec+1),maxnbd,
     &        nbrigid,mbdlst(bdsec),link(maxnbd),ncb(nbound),
     &        iscb(nbound+1),bdseg
      double precision coor(2,maxnbd),xpos(3,nbrigid),pwall(5,bdseg),
     &                 diab(nbrigid),diaa(nbrigid),h(maxnbd),
     &                 rwork(3*maxnbd),xymid(2,maxnbd)
c
      integer nwall,j1,j2,k,ngb,locx,locy,locx1,locy1
      integer i,kp,ia1,ib1,ncur1,nnxt1,nprv1,ncur2,nprv2
      double precision xcur1,xnxt1,xprv1,xcur2,xprv2,xmid1,ymid1,
     &                 ycur1,ynxt1,yprv1,ycur2,yprv2,xmid2,ymid2,
     &                 dis,ratio
      integer lcrv1,lcrv2
      double precision factor,factor1, x0,xp,xcoor,xploc
      data factor/1.5d0/,factor1/0.75d0/
c
c     build node link list
c     --------------------
      do i=1,nbd
         link(i) = i+1
      enddo
c
      nwall = nbound - nbrigid
      seclst(bdsec+1) = nbd+1
      iscb(1) = 1
      do i=1,nbound
         k = iscb(i) + nside(i)
         iscb(i+1) = k
         ia1 = seclst(iscb(i))
         ib1 = seclst(k)
         ncb(i) = ib1 - ia1
         link(ib1-1) = ia1
      enddo
c
c     initial mesh size
c     -----------------
      do i=1,nbound
         if ( i.le.nwall ) ratio = 1.d0/factor1
         ncur1 = seclst(iscb(i))
         xcur1 = coor(1,ncur1)
         ycur1 = coor(2,ncur1)
         do j1=1,ncb(i)
            nprv1 = ncur1
            xprv1 = xploc(xcur1)
            yprv1 = ycur1
            ncur1 = link(nprv1)
            xcur1 = xcoor(xprv1,coor(1,ncur1))
            ycur1 = coor(2,ncur1)
            h(ncur1) = sqrt((xcur1-xprv1)**2 + (ycur1-yprv1)**2)
         enddo
      enddo
c
c     uniform mesh
c     ------------
      if ( .not.refine ) return
c 
c     check through all nodes
c     -----------------------
      do i=1,nbound
        kp = i-nwall
c
c       loop over nodes on current boundary section
        lcrv1 = mbdlst(iscb(i))
        ncur1 = seclst(iscb(i))
        xcur1 = coor(1,ncur1)
        ycur1 = coor(2,ncur1)
        nnxt1 = link(ncur1)
        xnxt1 = xcoor(xcur1,coor(1,nnxt1))
        ynxt1 = coor(2,nnxt1)
        j1 = 1
        do 20 while ( j1.le.ncb(i) )
           j1    = j1 + 1
           x0    = xploc(xnxt1)
           nprv1 = ncur1
           xprv1 = xcoor(x0,xcur1)
           yprv1 = ycur1
           ncur1 = nnxt1
           xcur1 = x0
           ycur1 = ynxt1
           nnxt1 = link(ncur1)
           xnxt1 = xcoor(x0,coor(1,nnxt1))
           ynxt1 = coor(2,nnxt1)
c
c         search over all nodes on moving particles 
c         -----------------------------------------
          do 40 k=max(1,i-nwall+1),nbrigid
            ngb = k + nwall
            lcrv2 = mbdlst(iscb(ngb))
70          continue
            ncur2 = seclst(iscb(ngb))
            xcur2 = xcoor(x0,coor(1,ncur2))
            ycur2 = coor(2,ncur2)
            xmid1 = xcoor(x0,xymid(1,nprv1))
            ymid1 = xymid(2,nprv1)
            j2 = 1
            do 50 while( j2.le.ncb(ngb) )
               j2    = j2 + 1
               nprv2 = ncur2
               xprv2 = xcur2
               yprv2 = ycur2
               ncur2 = link(ncur2)
               xcur2 = xcoor(x0,coor(1,ncur2))
               ycur2 = coor(2,ncur2)
c
 60            continue
               xmid2 = xcoor(x0,xymid(1,nprv2))
               ymid2 = xymid(2,nprv2)
               dis = sqrt( (xmid2-xmid1)**2+(ymid2-ymid1)**2 )
c             write(*,'(7e12.5)') xcur1,ycur1,xcur2,ycur2,dis,h(ncur1)
c
c             mesh refinement on section 2
              if ( factor*dis .lt. h(ncur2) ) then
c                write(*,*) 'h2 refinement'
                link(nprv2) = nbd+1
                link(nbd+1) = ncur2
                xymid(1,nprv2) = 0.375*xprv2+0.75*xmid2-0.125*xcur2
                xymid(2,nprv2) = 0.375*yprv2+0.75*ymid2-0.125*ycur2
                coor(1,nbd+1)  = xmid2
                coor(2,nbd+1)  = ymid2
                xymid(1,nbd+1) = -0.125*xprv2+0.75*xmid2+0.375*xcur2
                xymid(2,nbd+1) = -0.125*yprv2+0.75*ymid2+0.375*ycur2
                h(ncur2) = sqrt((xcur2-xmid2)**2+(ycur2-ymid2)**2)
                h(nbd+1) = sqrt((xmid2-xprv2)**2+(ymid2-yprv2)**2)
                if ( lcrv2.lt.0 ) then
                   call movcoor(xpos(1,k),diaa(k),diab(k),
     &                          xymid(1,nprv2),xymid(2,nprv2))
                   call movcoor(xpos(1,k),diaa(k),diab(k),
     &                          xymid(1,nbd+1),xymid(2,nbd+1))
                elseif ( lcrv2.gt.0 .and. lcrv2.le.bdseg ) then
                   if ( int(pwall(1,lcrv2))/100.eq.1 ) then
                      call movcoor0 (pwall(2,lcrv2),pwall(3,lcrv2)
     &                    ,pwall(4,lcrv2),xymid(1,nprv2),xymid(2,nprv2))
                      call movcoor0 (pwall(2,lcrv2),pwall(3,lcrv2)
     &                    ,pwall(4,lcrv2),xymid(1,nbd+1),xymid(2,nbd+1))
                   else
                      stop 'pwall(1,*).ne.1** is not implemented'
                   endif
                endif
                ncur2 = nbd + 1
                xcur2 = xcoor(x0,coor(1,ncur2))
                ycur2 = coor(2,ncur2)
                nbd   = nbd + 1
                ncb(ngb) = ncb(ngb) + 1
                goto 60
              endif
c
c             mesh refinement on section 1
              if ( factor*dis .lt. h(ncur1) ) then
c               write(*,*) 'h1 refinement'
                link(nprv1) = nbd+1
                link(nbd+1) = ncur1
                xymid(1,nprv1) = 0.375*xprv1+0.75*xmid1-0.125*xcur1
                xymid(2,nprv1) = 0.375*yprv1+0.75*ymid1-0.125*ycur1
                coor(1,nbd+1) = xmid1
                coor(2,nbd+1) = ymid1
                xymid(1,nbd+1) = -0.125*xprv1+0.75*xmid1+0.375*xcur1
                xymid(2,nbd+1) = -0.125*yprv1+0.75*ymid1+0.375*ycur1
                h(ncur1) = sqrt((xcur1-xmid1)**2+(ycur1-ymid1)**2)
                h(nbd+1) = sqrt((xmid1-xprv1)**2+(ymid1-yprv1)**2)
                if ( lcrv1.lt.0 ) then
                   call movcoor(xpos(1,kp),diaa(kp),diab(kp),
     &                          xymid(1,nprv1),xymid(2,nprv1))
                   call movcoor(xpos(1,kp),diaa(kp),diab(kp),
     &                          xymid(1,nbd+1),xymid(2,nbd+1))
                elseif ( lcrv1.gt.0 .and. lcrv1.le.bdseg ) then
                   if ( int(pwall(1,lcrv1))/100.eq.1 ) then
                      call movcoor0 (pwall(2,lcrv1),pwall(3,lcrv1)
     &                    ,pwall(4,lcrv1),xymid(1,nprv1),xymid(2,nprv1))
                      call movcoor0 (pwall(2,lcrv1),pwall(3,lcrv1)
     &                    ,pwall(4,lcrv1),xymid(1,nbd+1),xymid(2,nbd+1))
                  else
                    stop 'pwall(1,*).ne.1** is not implemented'
                  endif
                endif
                nnxt1 = ncur1
                xnxt1 = xcoor(x0,xcur1)
                ynxt1 = ycur1
                ncur1 = nbd + 1
                xcur1 = xcoor(x0,coor(1,ncur1))
                ycur1 = coor(2,ncur1)
                nbd   = nbd + 1
                ncb(i)= ncb(i) + 1
                goto 70
              endif
c
  50        continue
c
 40       continue
c
c         mesh smoothing
c         --------------
80        ratio = h(nnxt1)/h(ncur1)
          if ( ratio .gt. 2.5 ) then
c            write(*,*) 'mesh smoothing 1',h(ncur1),h(nnxt1)
            link(ncur1) = nbd+1
            link(nbd+1) = nnxt1
            xmid1 = xcoor(x0,xymid(1,ncur1))
            ymid1 = xymid(2,ncur1)
            xymid(1,ncur1) = 0.375*xcur1+0.75*xmid1-0.125*xnxt1
            xymid(2,ncur1) = 0.375*ycur1+0.75*ymid1-0.125*ynxt1
            coor(1,nbd+1) = xmid1
            coor(2,nbd+1) = ymid1
            xymid(1,nbd+1) = -0.125*xcur1+0.75*xmid1+0.375*xnxt1
            xymid(2,nbd+1) = -0.125*ycur1+0.75*ymid1+0.375*ynxt1
            h(nnxt1) = sqrt((xnxt1-xmid1)**2+(ynxt1-ymid1)**2)
            h(nbd+1) = sqrt((xmid1-xcur1)**2+(ymid1-ycur1)**2)
            if ( lcrv1.lt.0 ) then
               call movcoor(xpos(1,kp),diaa(kp),diab(kp),
     &                      xymid(1,ncur1),xymid(2,ncur1))
               call movcoor(xpos(1,kp),diaa(kp),diab(kp),
     &                      xymid(1,nbd+1),xymid(2,nbd+1))
            elseif ( lcrv1.gt.0 .and. lcrv1.le.bdseg ) then
               if ( int(pwall(1,lcrv1))/100.eq.1 ) then
                  call movcoor0 (pwall(2,lcrv1),pwall(3,lcrv1)
     &                 ,pwall(4,lcrv1),xymid(1,ncur1),xymid(2,ncur1))
                  call movcoor0 (pwall(2,lcrv1),pwall(3,lcrv1)
     &                 ,pwall(4,lcrv1),xymid(1,nbd+1),xymid(2,nbd+1))
               else
                  stop 'pwall(1,*).ne.1** is not impelemented'
               endif
            endif
            nnxt1 = nbd + 1
            xnxt1 = xcoor(x0,coor(1,nnxt1))
            ynxt1 = coor(2,nnxt1)
            nbd   = nbd + 1
            ncb(i)= ncb(i) + 1
            goto 80
          elseif ( ratio .lt. 0.4 ) then
c            write(*,*) 'mesh smoothing 2',h(ncur1),h(nnxt1)
            link(nprv1) = nbd+1
            link(nbd+1) = ncur1
            xmid1 = xcoor(x0,xymid(1,nprv1))
            ymid1 = xymid(2,nprv1)
            xymid(1,nprv1) = 0.375*xprv1+0.75*xmid1-0.125*xcur1
            xymid(2,nprv1) = 0.375*yprv1+0.75*ymid1-0.125*ycur1
            coor(1,nbd+1) = xmid1
            coor(2,nbd+1) = ymid1
            xymid(1,nbd+1) = -0.125*xprv1+0.75*xmid1+0.375*xcur1
            xymid(2,nbd+1) = -0.125*yprv1+0.75*ymid1+0.375*ycur1
            h(ncur1) = sqrt((xcur1-xmid1)**2+(ycur1-ymid1)**2)
            h(nbd+1) = sqrt((xmid1-xprv1)**2+(ymid1-yprv1)**2)
            if ( lcrv1.lt.0 ) then
               call movcoor(xpos(1,kp),diaa(kp),diab(kp),
     &                      xymid(1,nprv1),xymid(2,nprv1))
               call movcoor(xpos(1,kp),diaa(kp),diab(kp),
     &                      xymid(1,nbd+1),xymid(2,nbd+1))
            elseif ( lcrv1.gt.0 .and. lcrv1.le.bdseg ) then
               if ( int(pwall(1,lcrv1))/100.eq.1 ) then
                  call movcoor0 (pwall(2,lcrv1),pwall(3,lcrv1)
     &                 ,pwall(4,lcrv1),xymid(1,nprv1),xymid(2,nprv1))
                  call movcoor0 (pwall(2,lcrv1),pwall(3,lcrv1)
     &                 ,pwall(4,lcrv1),xymid(1,nbd+1),xymid(2,nbd+1))
               else
                  stop 'pwall(1,*).ne.1** is not impelemented'
               endif
            endif
            nnxt1 = ncur1
            xnxt1 = xcoor(x0,xcur1)
            ynxt1 = ycur1
            ncur1 = nbd + 1
            xcur1 = xcoor(x0,coor(1,ncur1))
            ycur1 = coor(2,ncur1)
            nbd   = nbd + 1
            ncb(i)= ncb(i) + 1
            goto 80
          endif 
c
20      continue
c
      enddo
c
c     reordering of nodes
c     -------------------
      do i=1,nbd
         rwork(i) = h(i)
      enddo
      locx = nbd
      locy = locx + nbd
      locx1 = locy + nbd
      locy1 = locx1 + nbd
      j1 = 0
      j2 = 1
      k  = 0
      do i=1,nbound
         ncur1 = seclst(iscb(i))
         nnxt1 = link(ncur1)
         seclst(j2) = j1+1
         j2 = j2 + 1
         k = k + ncb(i)
         do while ( j1.lt.k)
            j1 = j1 + 1
            rwork(locx+j1)  = coor(1,ncur1)
            rwork(locy+j1)  = coor(2,ncur1)
            rwork(locx1+j1) = xymid(1,ncur1)
            rwork(locy1+j1) = xymid(2,ncur1)
            h(j1) = 0.5*(rwork(ncur1)+rwork(nnxt1))
            ncur1 = nnxt1
            nnxt1 = link(nnxt1)
            if ( ncur1 .eq. seclst(j2) ) then
               seclst(j2) = j1+1
               j2 = j2 + 1
            endif 
         enddo
      enddo
      seclst(bdsec+1) = nbd + 1
c
      do i=1,nbd
         coor(1,i)  = xploc(rwork(locx+i))
         coor(2,i)  =       rwork(locy+i)
         xymid(1,i) = xploc(rwork(locx1+i))
         xymid(2,i) =       rwork(locy1+i)
c         write(20,'(i4,5e13.5)') i,coor(1,i),coor(2,i),h(i)
c     &                            ,xymid(1,i),xymid(2,i)
      enddo
c
      if ( 2*nbd.gt.maxnbd ) then
         write(*,*) 'maxnbd is too small, it must >',2*nbd
         stop
      endif
c
      return
      end
c
c
      subroutine MOVCOOR ( xpos,a,b, x,y )
c
c     Move the node (x,y) onto the surface of a given particle
c
      implicit none
      double precision xpos(3),a,b, x,y
      double precision dx,dy,d2,co,si,r, xcoor,xploc
c
      dx = xcoor(xpos(1),x) - xpos(1)
      dy = y - xpos(2)
      d2 = dx*dx + dy*dy
      co = dcos(xpos(3))
      si = dsin(xpos(3))
      r  = 0.5*a/dsqrt( d2 + ((a/b)**2-1.d0)*(dx*co+dy*si)**2)
      x  = xploc(xpos(1) + dx*r)
      y  =       xpos(2) + dy*r
      return
      end   
c
      subroutine MOVCOOR0 ( x0,y0,r, x,y )
c
c     Move the node (x,y) onto the surface of a circle (x0,y0,r)
c
      implicit none
      double precision x0,y0,r, x,y
      double precision dx,dy,d2,a
c
      dx = x - x0
      dy = y - y0
      d2 = dsqrt(dx*dx + dy*dy)
      a  = abs(r)/d2
      x  = x0 + dx*a
      y  = y0 + dy*a
      return
      end
