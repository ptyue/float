c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c  Automatic mesh generator "RMESH" for partflow.f
c    with periodic boundary condition in the x-direction
c
c   by Howard H. Hu, Univ. of Pennsylvania
c
c   modified  Aug.20, 1991, for elliptic cylinders
c             Oct.25, 1991, using INCLUDE statement
c             Nov.16, 1991, for curved boundaries
c             Mar.12, 1992, move rmhini & mshptg to polylib
c             July 18, 1993, updated
c             Jan.20, 1995, periodic mesh generator!
c             Feb.9, 1995, rmhini removed, local refinement added
c             Feb.14, 1995, cleanup
c             Oct.25, 1995, working version on sgi, sun ...
c             Nov.20, 1995, element area and aspect ratio are calculated
c             March 4, 1996, mapping is applied only for boundary nodes
c             May 3, 1996, removed all implicit double precision and
c                          parameter statements, better memory mangement
c             June 24, 1996, the periodic and non-periodic version of
c                            mesh generators are merged.
c             Sept.13, 1996, rectangluar particles are included
c             Dec.30, 1996, new mesh refinement scheme is added
c             Oct.13, 1998, for multi-domain mesh generation
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c***********************************************************************
      subroutine RMESH (itime,flowtype,pwall,bdseg,bdrefn,ncr,ny,
     &                  nbrigid,nbfluid,nbelast,nbtotal,
     &                  nbshp,diaa,diab,xpos,
     &                  nvert,nnode,nelem,inod,x,y,
     &                  elmref,elmnec,elmara,elmasp,
     &                  nbd,ibdnod,bdsec,seclst,nbound,nside,
     &                  maxvrt,maxnod,maxelt,maxnbd,maxbdp, 
     &                  iwork,leniw,rwork,lenrw )
c
c   INPUT:
c       itime   = time step index
c       flowtype = index for the domain
c               =  1: periodic flow in x-direction
c               =  2: flow in an infinite channel in x-direction
c               =  3: flow in a channel closed in x-direction
c               =  4: fluidization in x-direction
c               = 11: flow inside a circular cylinder
c               = 12: flow in an annulus between two cylinders
c               = 99: flow in an arbitrary geometry
c       pwall   = parameters for the domain - refer to the manual!
c       bdseg   = number of boundary segments
c                 (include periodic boundary segments)
c       bdrefn  = .true. refines the segments on boundary sections
c               = .false. uses uniform segments on boundary sections
c       ncr     : number of segments on particle surface
c       ny      : number of segments across the channel width
c       nbrigid = number of rigid particles
c       nbfluid : number of deformable fluid particles
c       nbelast : number of deformable elastic particles
c       nbtotal = nbrigid+nbfluid+nbelast
c       nbshp   = xyz  shape index of the particle
c                 where  x = 0 rigid particle
c                        x = 1 deformable fluid particle
c                        x = 2 deformable elastic solid partile
c                        yz = 1 circular / elliptical
c                        yz = 2 rectangular
c       diaa    : diameter a of particle
c       diab    : diameter b of particle
c       xpos    : positions of particles 
c       maxvrt  = maximum number of vertices
c       maxnod  = maximum number of nodes
c       maxelt  = maximum number of elements
c       maxnbd  = maximum number of boundary nodes
c       maxbdp  = maximum number of boundary sections
c
c   OUTPUT:
c       nvert : number of vertices 
c       nnode : number of nodes
c       nelem : number of elements
c       inod : element description table
c       x,y :  x and y coordinates
c       elmref: element reference number (subdomain number)
c       elmnec: element connectivity
c       elmara: area of element
c       elmasp: aspect ratio of element
c       nbd : number of boundary nodes
c       ibdnod : list of boundary nodes
c       bdsec : number of boundary sections
c       seclst: index of the first node on each boundary section
c       nbound : number of closed boundaries
c       nside : number of sections on each closed boundary
c
c   WORKING ARRAY:
c       iwork(leniw), rwork(lenrw)
c         leniw = 12*maxelt+6*maxvrt+4*maxbdp+maxnbd+3
c         lenrw = 3*maxvrt+2*maxnbd+1
c
c***********************************************************************
      implicit none
      logical bdrefn
      integer flowtype,ncr,ny,nbrigid,maxvrt,maxnod,maxelt,maxnbd,maxbdp
     &       ,bdseg,nbfluid,nbelast,itime, nbtotal
      real*8 pwall(5,bdseg),diaa(nbtotal),diab(nbtotal),xpos(3,nbtotal)
      integer nvert,nnode,nelem,inod(6,maxelt),elmnec(3,maxelt),
     &        elmref(maxelt),nbd,ibdnod(maxnbd),bdsec,seclst(maxbdp),
     &        nbound,nside(maxbdp), nbshp(nbtotal)
      real*8 x(maxnod),y(maxnod),elmara(maxelt),elmasp(maxelt)
      integer leniw,iwork(leniw),lenrw
      real*8 rwork(lenrw)
      real*8 h0
c
c     local variables
c     ---------------
      integer itermx
      real*8 omega,eps,puis,coef
      parameter (itermx=20,omega=1.0,eps=0.005,puis=0.25,coef=0.75)
      logical regul
      integer nsd
c
      integer i,nba,ierr,k,j,n,n1,n2,n3
      real*8 scal,tmp,c1,c2,c3,y1,y2,y3,x21,x31,y21,y31,xcoor
      real*8 xmax,xmin
      integer locnu,loctri,lociendpt,locibdnod,locrmhnel,
     &        locrmhlen,locmbdlst,lociw,loccoor,loch,loctmp,locrw,locsd
     &       ,locxymid,lcoorcp
      logical xperid
c
      real*8 pi
      data pi/3.14159265358979d0/
C     Support for Petsc timing
c      include '../include/petsc.h'
c
c      call PLogEventBegin(REMESH_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     initializtaion
c     --------------
      locnu     = 1
      loctri    = locnu     + 6*maxelt
      lociendpt = loctri    + 3*maxelt+maxvrt+1
      locibdnod = lociendpt + 2*maxvrt
      locrmhnel = locibdnod + maxnbd
      locrmhlen = locrmhnel + maxelt
      locmbdlst = locrmhlen + maxelt
      locsd     = locmbdlst + maxbdp
      lociw     = locsd     + 2*(nbfluid+nbelast+1)
      if ( lociw.gt.leniw ) then
         print*, 'error in RMESH: leniw is too small, must>',lociw
         stop
      endif
c
      loccoor  = 1
      loch     = loccoor  + 2*maxvrt
      loctmp   = loch     + 2*maxvrt
      locxymid = loctmp   + 3*maxnbd
      lcoorcp  = locxymid + 2*maxnbd
      locrw    = lcoorcp  + 2*maxnbd
      if ( locrw.gt.lenrw ) then
         print*, 'error in RMESH: lenrw is too small, must>',locrw
         stop
      endif
c
      do i=1,maxbdp
         iwork(locmbdlst+i-1) = 0
         nside(i) = 1
      enddo
c
c     check pwall
c     -----------
      xperid = .false.
      xmax = pwall(2,1)
      xmin = pwall(2,1)
      do i=1,bdseg
         if ( abs(pwall(1,i)).lt.0.1 ) xperid = .true.
         xmax = max(xmax,pwall(2,i))
         xmin = min(xmin,pwall(2,i)) 
      enddo
c
c     generate points on the outer boundary
c     -------------------------------------
c     print*, '  generating outer boundary'
      call OUTBD ( flowtype,pwall,bdseg,ny,rwork(loccoor),nbd,bdsec,
     &             seclst,nbound,nside,h0,
     &             rwork(locxymid),iwork(locmbdlst),maxnbd,maxbdp)
c
c     generate points on rigid particle surface
c     -----------------------------------------
c     print*, '  generating points on particle surfaces'
      if ( itime.eq.0 ) then
         call PARTCBD ( ncr,nbtotal,nbshp,diaa,diab,xpos,
     &                  rwork(loccoor),nbd,bdsec,seclst,nbound,
     &                  rwork(locxymid),iwork(locmbdlst),maxnbd,maxbdp )
      else
         call PARTCBD ( ncr,nbrigid,nbshp,diaa,diab,xpos,
     &                  rwork(loccoor),nbd,bdsec,seclst,nbound,
     &                  rwork(locxymid),iwork(locmbdlst),maxnbd,maxbdp )
c
c     generate points on deformable particle surface
c     ----------------------------------------------
         if(nbfluid+nbelast.ne.0)then
         call DEFPART ( nbfluid+nbelast,ncr,diaa(nbrigid+1),ibdnod,x,y,
     &                  rwork(loccoor),nbd,bdsec,seclst,nbound,
     &                  rwork(locxymid),iwork(locmbdlst),maxnbd,maxbdp )
         endif
      endif
c
c      do i=1,nbd
c         n = loccoor-2+2*i
c         x(i) = rwork(n)
c         y(i) = rwork(n+1)
c         write(10,'(i5,4e13.5)') i,x(i),y(i)
c     &                   ,rwork(locxymid-2+2*i),rwork(locxymid-1+2*i)
c      enddo
c
c     refine the boundary nodes
c     -------------------------
c      print*,'  boundary node refinement'
      call BDNRFN (bdrefn,rwork(loccoor),rwork(loch),bdsec,seclst,nbd,
     &             nbound,nside,rwork(locxymid),iwork(locmbdlst),
     &             pwall,bdseg,nbtotal,
     &             diaa,diab,xpos,maxnbd,iwork(lociendpt),rwork(loctmp),
     &             iwork(loctri),iwork(loctri+maxbdp) )
c
      do i=1,nbd
         iwork(lociendpt+2*i-2) = i
         iwork(lociendpt+2*i-1) = i+1
         iwork(locibdnod-1+i) = i
      enddo
      k = 1
      do i=1,nbound
         n1 = seclst(k)
         n2 = seclst(k+nside(i))
         iwork(lociendpt+2*n2-3) = n1
         k = k+nside(i)
      enddo
c
      do i=1,nbd
         rwork(lcoorcp+2*i-2) = rwork(loccoor+2*i-2)
         rwork(lcoorcp+2*i-1) = rwork(loccoor+2*i-1)
c         write(11,'(i5,5e12.5)') i,rwork(loccoor+2*i-2),
c     &                           rwork(loccoor+2*i-1),rwork(loch-1+i),
c     &                      rwork(locxymid-2+2*i),rwork(locxymid-1+2*i)
      enddo
c
c
c     domain and meshing information
c     ------------------------------
      iwork(locsd)   = 1
      iwork(locsd+1) = 1
      nsd = 1
      do k=1,nbfluid+nbelast
         nsd = nsd + 1
         iwork(locsd+2*k)   = -seclst(bdsec-nbfluid-nbelast+k)
         iwork(locsd+2*k+1) = 1+nbrigid+k
      enddo
      nvert = nbd
      nba = nbd
c
      if ( xperid ) then
c       periodic domain, mapping the coordinates
c       ----------------------------------------
c        print*,'  mapping coordinates'
        scal = 2.*pi/(xmax-xmin)
        call MAP ( rwork(loccoor), nbd, scal )
c
c        do i=1,nbd
c          n = loccoor-2+2*i
c          write(12,'(i5,2e13.5)') i,rwork(n),rwork(n+1)
c        enddo
c
      endif
c
c     mesh generation
c     ---------------
c      print*,'  mesh generation'
      regul = .true.
      call  MESH2D(rwork(loccoor),rwork(lcoorcp),rwork(loch),
     &             iwork(locnu),nvert,
     &             maxvrt,iwork(lociendpt),nba,iwork(locsd),nsd,
     &             elmref,nelem,regul,h0,0, ierr,iwork(loctri))
      if ( ierr.ne.0 ) then
         print*, ' error in mesh generation mesh2d; ERR = ',ierr
         stop
      endif
c
c      do i=1,nvert
c         n = loccoor-2+2*i
c         x(i) = rwork(n)
c         y(i) = rwork(n+1)
c      enddo
c      do i=1,nelem
c         do j=1,3
c           inod(j,i)   = iwork(locnu+6*i+j-7)
c           elmnec(j,i) = iwork(locnu+6*i+j-4)/8
c         enddo
c      enddo
c      call PRTMSH('hh.msh1',3,6,nvert,nvert,nelem,inod,x,y,
c     &             elmref,nbd,iwork(locibdnod),bdsec,seclst,
c     &             nbound,nside, elmnec,15 )
c      pause '111'
c 
c
c     re-number the elements
c     ----------------------
c      print*, '  renumbering mesh'
      call RNBNEL ( iwork(locnu),nvert,nelem,
     &              iwork(locrmhnel),iwork(locrmhlen) )
c
c     generate mid-nodes
c     ------------------
c      print*,'  generating mid-nodes'
      call MNDNOD ( nvert,nnode,nelem,iwork(locnu),rwork(loccoor),
     &              iwork(locrmhlen),elmref,
     &              nbd,iwork(locibdnod),bdsec,seclst,nbound,nside,
     &              nbtotal,diaa,diab,xpos,inod,x,y,elmnec, 
     &              ibdnod,rwork(locxymid),iwork(locmbdlst),pwall,bdseg,
     &              maxnod,maxnbd,iwork(loctri) )
c
c     compute the area and aspect ratio of elements
c     ---------------------------------------------
      do i=1,nelem
c         do j=1,3
c            if ( elmnec(j,i).gt.nelem ) then
c               write(*,*) 'RMESH: elmnec error:',j,i,elmnec(j,i)
c               stop
c            endif
c         enddo
         n1 = inod(1,i)
         n2 = inod(2,i)
         n3 = inod(3,i)
         c1 = x(n1)
         c2 = xcoor(c1,x(n2))
         c3 = xcoor(c1,x(n3))
         y1 = y(n1)
         y2 = y(n2)
         y3 = y(n3)
         x21 = c2-c1
         x31 = c3-c1
         y21 = y2-y1
         y31 = y3-y1
         tmp = x21*y31 - y21*x31
         if ( tmp.le.0.d0 ) then
            write(*,*) 'RMESH: element with negative area!'
            stop
         endif
         elmara(i) = tmp
         elmasp(i) = max(x21**2+y21**2,x31**2+y31**2,
     &                   (c3-c2)**2+(y3-y2)**2) / tmp
      enddo
c
C     Support for Petsc timing
c      call PLogEventEnd(REMESH_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c***********************************************************************
      subroutine OUTBD ( flowtype,pwall,bdseg,ny,coor, nbd,bdsec,seclst,
     &                   nbound,nside,h0,xymid,mbdlst, maxnbd,maxbdp )
c
c     Generates points on the outer boundary - periodic channel
c	
c     INPUT: 
c       flowtype = index for the domain
c       pwall  : channel information
c       ny     : number of segments across channel
c	maxnbd : max. number of boundary nodes
c	maxbdp : max. number of boundary parts
c       h0    : background mesh size
c
c     OUTPUT:
c	coor  : x and y coordinates of the points generated
c       nbd   : number of boundary nodes
c	bdsec : number of boundary sections
c	seclst: index of the first node on each boundary section
c       nside : number of sections on each closed boundary
c               = 1 by default
c       xymid = coordinates for the mid nodes
c       mbdlst = n: index for the boundary segment in pwall
c                n=0 for straight segment (by default)
c***********************************************************************
      implicit none
      integer flowtype
      integer ny,maxnbd,maxbdp,nbd,bdsec,seclst(maxbdp),nbound
     &       ,nside(maxbdp),mbdlst(maxbdp), bdseg
      real*8 pwall(5,bdseg),coor(2,maxnbd),h0,xymid(2,maxnbd)
c
c     local variables
      integer i,n,k
      real*8 xleng,rad,ceta,tmp,dx,dy,ymax,ymin
      logical xperid
c
      if ( flowtype.eq.1 ) then
c     flow in a periodic channel
c     --------------------------
c
         xleng = pwall(4,2)-pwall(2,2)
         n     = ny*xleng/(pwall(3,1)-pwall(5,1)) + 0.5
         h0    = xleng/float(n)
c
         nbd = 0
         bdsec = 1
         seclst(bdsec) = nbd+1
         do i=1,n
            nbd = nbd+1
            coor(1,nbd)  = pwall(2,2)+(i-1)*h0
            coor(2,nbd)  = pwall(3,2)
            xymid(1,nbd) = pwall(2,2)+(i-0.5)*h0
            xymid(2,nbd) = pwall(3,2)
         enddo
c 
         bdsec = 2
         seclst(bdsec) = nbd+1
         do i=1,n
            nbd = nbd+1
            coor(1,nbd)  = pwall(2,4)-(i-1)*h0
            coor(2,nbd)  = pwall(3,4)
            xymid(1,nbd) = pwall(2,4)-(i-0.5)*h0
            xymid(2,nbd) = pwall(3,4)
         enddo
c
         nbound = 2
c
      elseif ( flowtype.eq.2 .or. flowtype.eq.3 .or. 
     &         flowtype.eq.4 ) then
c     flow in a non-periodic channel
c     ------------------------------
c
         xleng = pwall(4,2)-pwall(2,2)
         n     = ny*xleng/(pwall(3,1)-pwall(5,1)) + 0.5
         h0    = xleng/float(n)
c
         nbd = 0
         bdsec = 1
         seclst(bdsec) = nbd+1
         do i=1,ny
            nbd = nbd + 1
            coor(1,nbd)  = pwall(2,1)
            coor(2,nbd)  = pwall(3,1) + 
     &                     (i-1)*(pwall(5,1)-pwall(3,1))/float(ny)
            xymid(1,nbd) = pwall(2,1)
            xymid(2,nbd) = pwall(3,1) +
     &                     (i-0.5)*(pwall(5,1)-pwall(3,1))/float(ny)
         enddo
c
         bdsec = 2
         seclst(bdsec) = nbd+1
         do i=1,n
            nbd = nbd+1
            coor(1,nbd)  = pwall(2,2)+(i-1)*h0
            coor(2,nbd)  = pwall(3,2)
            xymid(1,nbd) = pwall(2,2)+(i-0.5)*h0
            xymid(2,nbd) = pwall(3,2)
         enddo
c
         bdsec = 3
         seclst(bdsec) = nbd+1
         do i=1,ny
            nbd = nbd + 1
            coor(1,nbd)  = pwall(2,3)
            coor(2,nbd)  = pwall(3,3) + 
     &                    (i-1)*(pwall(5,3)-pwall(3,3))/float(ny) 
            xymid(1,nbd) = pwall(2,3)
            xymid(2,nbd) = pwall(3,3) +
     &                    (i-0.5)*(pwall(5,3)-pwall(3,3))/float(ny)
         enddo
c
         bdsec = 4
         seclst(bdsec) = nbd+1
         do i=1,n
            nbd = nbd+1
            coor(1,nbd)  = pwall(2,4)-(i-1)*h0
            coor(2,nbd)  = pwall(3,4)
            xymid(1,nbd) = pwall(2,4)-(i-0.5)*h0
            xymid(2,nbd) = pwall(3,4)
         enddo
c
         nbound   = 1
         nside(1) = 4
c
      elseif ( flowtype.eq.11 ) then
c     flow in a circular cylinder
c     ---------------------------
         rad = abs(pwall(4,1))
         h0  = 6.28*rad/float(ny)
c
         nbd = 0
         bdsec = 1
         seclst(bdsec) = nbd+1
         do i=1,ny
            nbd = nbd+1
            ceta = (i-1)*6.28318530/ny
            coor(1,nbd)  = pwall(2,1) + rad*dcos(ceta)
            coor(2,nbd)  = pwall(3,1) + rad*dsin(ceta)
            ceta = (i-0.5)*6.28318530/ny
            xymid(1,nbd) = pwall(2,1) + rad*dcos(ceta)
            xymid(2,nbd) = pwall(3,1) + rad*dsin(ceta)
         enddo
         mbdlst(bdsec) = 1
c 
         nbound = 1
c
      elseif ( flowtype.eq.12 ) then
c     flow in an annulus between two circular cylinders
c     -------------------------------------------------
c        outer cylinder
         rad = abs(pwall(4,1))
         h0  = 6.28*rad/float(ny)
c
         nbd = 0
         bdsec = 1
         seclst(bdsec) = nbd+1
         do i=1,ny
            nbd = nbd+1
            ceta = (i-1)*6.28318530/ny
            coor(1,nbd)  = pwall(2,1) + rad*dcos(ceta)
            coor(2,nbd)  = pwall(3,1) + rad*dsin(ceta)
            ceta = (i-0.5)*6.28318530/ny
            xymid(1,nbd) = pwall(2,1) + rad*dcos(ceta)
            xymid(2,nbd) = pwall(3,1) + rad*dsin(ceta)
         enddo
         mbdlst(bdsec) = 1
c
c        inner cylinder
         rad = abs(pwall(4,2))
         n   = ny*abs(pwall(4,2)/pwall(4,1))+0.5
c
         bdsec = 2
         seclst(bdsec) = nbd+1
         do i=n,1,-1
            nbd = nbd+1
            ceta = (i-1)*6.28318530/n
            coor(1,nbd)  = pwall(2,2) + rad*dcos(ceta)
            coor(2,nbd)  = pwall(3,2) + rad*dsin(ceta)
            ceta = (i-1.5)*6.28318530/n
            xymid(1,nbd) = pwall(2,2) + rad*dcos(ceta)
            xymid(2,nbd) = pwall(3,2) + rad*dsin(ceta)
         enddo
         mbdlst(bdsec) = 2
c 
         nbound = 2
c
      elseif ( flowtype.eq.99 ) then
c     flow in an arbitrary geometry
c     -----------------------------
         ymax = pwall(3,1)
         ymin = pwall(3,1)
         do i=2,bdseg
            ymax = max(ymax,pwall(3,i))
            ymin = min(ymin,pwall(3,i))
         enddo
         h0 = (ymax-ymin)/float(ny)
c
         xperid = .false.
         nbd = 0
         bdsec = 0
         do k=1,bdseg
            if ( abs(pwall(1,k)).gt.0.01) then
               bdsec = bdsec+1
               seclst(bdsec) = nbd+1
               tmp = dsqrt((pwall(4,k)-pwall(2,k))**2 +
     &                     (pwall(5,k)-pwall(3,k))**2 )
               n  = int(tmp/h0 + 0.5)
               dx = (pwall(4,k)-pwall(2,k))/float(n)
               dy = (pwall(5,k)-pwall(3,k))/float(n)
               do i=1,n
                  nbd = nbd + 1
                  coor(1,nbd)  = pwall(2,k) + (i-1)*dx
                  coor(2,nbd)  = pwall(3,k) + (i-1)*dy
                  xymid(1,nbd) = pwall(2,k) + (i-0.5)*dx
                  xymid(2,nbd) = pwall(3,k) + (i-0.5)*dy
               enddo
            else
               xperid = .true.
               nside(1) = bdsec
            endif
         enddo
c
         if ( .not.xperid ) then
            nbound   = 1
            nside(1) = bdsec
         else
            nbound   = 2
            nside(2) = bdsec-nside(1)
         endif
c
      else
c
         print*,'RMESH: flowtype=',flowtype,' is not implemented!'
         stop
      endif
      return
      end
c
c***********************************************************************
      subroutine PARTCBD (ncr,nbrigid,nbshp,diaa,diab,xpos,
     &                    coor,nbd,bdsec,seclst,nbound,
     &                    xymid,mbdlst,maxnbd,maxbdp )
c
c     This routine generates points on surface of moving rigid particles
c     Note: arrangement of nodes should be in the clockwise direction
c
c   INPUT:
c       diaa(nbrigid),diab(nbrigid): parameters for particle and domain
c       xpos(3,nbrigid):  positions of particles
c       ncr     : number of segments on particle surface
c       maxnbd : max. number of boundary nodes
c       maxbdp : max. number of boundary parts
c
c   OUTPUT:
c       nbd    : number of boundary nodes
c       coor   : x and y coordinates of nodes generates
c       bdsec  : number of boundary sections
c       seclst : index of the first node on each boundary section
c       xymid = coordinates for the mid nodes
c       mbdlst : index for the curved boundary
c                = -k for particle k
c***********************************************************************
      implicit none
      integer ncr,nbrigid,maxnbd,maxbdp, nbd,
     &        bdsec,seclst(maxbdp),nbound,mbdlst(maxbdp),nbshp(nbrigid)
      real*8 diaa(nbrigid),diab(nbrigid),xpos(3,nbrigid),coor(2,maxnbd)
     &      ,xymid(2,maxnbd)
c
c     local variables
      integer i,k,ncr1
      real*8 pi2,dbt0,dbt,beta,alp,r,r1,r2,fac, xploc
      real*8 xycoor(2,200)
c
      pi2 = 2.*3.14159265358979d0
      dbt0 = pi2/float(ncr)
c
      do k=1,nbrigid
         bdsec = bdsec + 1
         seclst(bdsec) = nbd+1
c
        if ( mod(nbshp(k),100).eq.1 ) then
c
c         circular/elliptical particles
c         -----------------------------
           mbdlst(bdsec) = -k
           beta = 0.d0
           r1   = 0.5*diaa(k)
           dbt = dbt0*0.5*diab(k)/r1
           do while ( beta .gt. -pi2+0.3*dbt )
              alp = xpos(3,k) + beta
              nbd = nbd + 1
              r1 = 0.5*diaa(k)*diab(k)/dsqrt( diab(k)**2+(diaa(k)**2
     &            -diab(k)**2)*(sin(beta))**2 )
              coor(1,nbd) = xploc(xpos(1,k) - r1*sin(alp))
              coor(2,nbd) =       xpos(2,k) + r1*cos(alp)
c
              if (  beta-1.3*dbt .le. -pi2 ) dbt  = pi2 + beta
              r2 = 0.5*diaa(k)*diab(k)/dsqrt( diab(k)**2+(diaa(k)**2
     &            -diab(k)**2)*(sin(beta-0.5*dbt))**2 )
              fac = r2/(r1+r2)
              xymid(1,nbd) = xploc(xpos(1,k) - r2*sin(alp-fac*dbt))
              xymid(2,nbd) =       xpos(2,k) + r2*cos(alp-fac*dbt)
c
              dbt = dbt0*0.5*diab(k)/r2
              beta = beta - dbt
           enddo
c
        elseif ( mod(nbshp(k),100).eq.2 ) then
c
c         rectangular particles
c         ---------------------
           ncr1 = ncr
           call RECTMSH (0.5*diaa(k),0.5*diab(k),ncr1,xpos(3,k),xycoor)
           do i=1,ncr1,2
              nbd = nbd + 1
              coor(1,nbd)  = xploc(xpos(1,k) + xycoor(1,2*i-1))
              coor(2,nbd)  =       xpos(2,k) + xycoor(2,2*i-1)
              xymid(1,nbd) = xploc(xpos(1,k) + xycoor(1,2*i))
              xymid(2,nbd) =       xpos(2,k) + xycoor(2,2*i)
           enddo
c
        endif
c
      enddo
c
      nbound = nbound + nbrigid
c
      if ( 2*nbd.gt.maxnbd ) then
         write(*,*) 'maxnbd is too small, it must >',2*nbd
         stop
      endif
      if ( bdsec+1.gt.maxbdp ) then
         write(*,*) 'maxbdp is too small, it must >',bdsec
         stop
      endif
c
      return
      end
c
c***********************************************************************
      subroutine DEFPART (nbfluid,ncr,diaa,ibdnod,x,y, coor,nbd,bdsec,
     &                    seclst,nbound, xymid,mbdlst,maxnbd,maxbdp )
c
c     This routine generates points on the surface of moving deformable
c       particles
c
c   INPUT:
c       nbfluid : number of deformable particles
c       seclst : section list for the previous mesh
c       ibdnod : boundary nodes for the previous mesh
c       x,y:     coordinates for the previous mesh
c       ncr,diaa: parameters for particle
c       maxnbd : max. number of boundary nodes
c       maxbdp : max. number of boundary parts
c
c   OUTPUT:
c       nbd    : number of boundary nodes
c       coor   : x and y coordinates of nodes generates
c       bdsec  : number of boundary sections
c       seclst : index of the first node on each boundary section
c       xymid  : coordinates for the mid nodes
c       mbdlst : index for the curved boundary
c                = -k for particle k
c***********************************************************************
      implicit none
      integer ncr,nbfluid,maxnbd,maxbdp,
     &        nbd,bdsec,seclst(maxbdp),nbound,mbdlst(maxbdp),ibdnod(*)
      real*8 coor(2,maxnbd),xymid(2,maxnbd),diaa(nbfluid),x(*),y(*)
c
c     local variables
      integer k,i1,i2,nel,m,n,ib,i,j
      real*8 ds0,ds,ds1,xn(3),yn(3),xcoor, x1,x3,y1,y3,a,b,xc,yc,r2
c
      do k=1,nbfluid
         ds0 = (2.d0*3.14159265358979d0/float(ncr))**2
         r2  = 0.25*diaa(k)**2
         bdsec = bdsec + 1
         i1 = seclst(bdsec)
         i2 = seclst(bdsec+1)
         nel = (i2-i1)/2
c
         seclst(bdsec) = nbd+1
         mbdlst(bdsec) = 100+k
c
         m = 1
         do while ( m.le.nel)
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
c           radius of curvature
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
               if ( a.gt.r2 ) then
                  a = r2
               elseif ( a.lt.0.1*r2 ) then
                  a = 0.1*r2
               endif
            else
               a  = r2
            endif
c
            ds1 = (xn(2)-xn(1))**2 + (yn(2)-yn(1))**2
            ds  = (xn(3)-xn(1))**2 + (yn(3)-yn(1))**2
            if ( ds1 .lt. ds*(1.d0/3.d0)**2 ) then
c              reset the middle node at eta=0.75
               xn(2) = -0.125*xn(1)+0.75*xn(2)+0.375*xn(3)
               yn(2) = -0.125*yn(1)+0.75*yn(2)+0.375*yn(3)
            elseif ( ds1 .gt. ds*(2.d0/3.d0)**2 ) then
c              reset the middle node at eta=0.25
               xn(2) = 0.375*xn(1)+0.75*xn(2)-0.125*xn(3)
               yn(2) = 0.375*yn(1)+0.75*yn(2)-0.125*yn(3)
            endif 
c
c           refinement of the segment
c           -------------------------
            ds = ds/a
            if ( ds.ge.3.d0*ds0 ) then
               nbd = nbd + 1
               coor(1,nbd)  = xn(1)
               coor(2,nbd)  = yn(1)
               xymid(1,nbd) = 0.375*xn(1)+0.75*xn(2)-0.125*xn(3)
               xymid(2,nbd) = 0.375*yn(1)+0.75*yn(2)-0.125*yn(3)
               nbd = nbd + 1
               coor(1,nbd)  = xn(2)
               coor(2,nbd)  = yn(2)
               xymid(1,nbd) = -0.125*xn(1)+0.75*xn(2)+0.375*xn(3)
               xymid(2,nbd) = -0.125*yn(1)+0.75*yn(2)+0.375*yn(3)
               m = m + 1
c
c           merge two segments
c           ------------------
            elseif ( ds.le.0.25*ds0 ) then
               if ( m.lt.nel ) then
                  nbd = nbd + 1
                  coor(1,nbd)  = xn(1)
                  coor(2,nbd)  = yn(1)
                  xymid(1,nbd) = xn(3)
                  xymid(2,nbd) = yn(3)
               endif
               m = m + 2
c
c           do nothing
c           ----------
            else
               nbd = nbd + 1
               coor(1,nbd)  = xn(1)
               coor(2,nbd)  = yn(1)
               xymid(1,nbd) = xn(2)
               xymid(2,nbd) = yn(2)
               m = m + 1
            endif 
        enddo
c
      enddo
c
      nbound = nbound + nbfluid
c
      if ( 2*nbd.gt.maxnbd ) then
         write(*,*) 'maxnbd is too small, it must >',2*nbd
         stop
      endif
      if ( bdsec+1.gt.maxbdp ) then
         write(*,*) 'maxbdp is too small, it must >',bdsec
         stop
      endif
c
      return
      end
c
c***********************************************************************
      subroutine MAP ( coor,nbd,scal )
c
c     this routine maps the physical domain to a mathematical surface
c     for mesh generation
c
      implicit none
      integer nbd
      real*8 coor(2,nbd),scal
c
      integer i
      double precision ceta,amp
c
      do i=1,nbd
         ceta = scal*coor(1,i)
         amp = 1+coor(2,i)
         coor(1,i) =  amp*dcos(ceta)
         coor(2,i) = -amp*dsin(ceta)
      enddo
c
      return
      end
c
c***********************************************************************
      subroutine RNBNEL ( nu,nvert,nelem, rmhnel,rmhlen )
c
c     The routine renumbers the mesh using Cuthill & McKee's algorithm
c
c     OUTPUT 
c       rmhnel(new elem numb) = old elem numb
c       rmhlen(old elem numb) = new elem numb
c***********************************************************************
      implicit none
      integer nvert,nelem,nu(6,nelem),rmhnel(nelem),rmhlen(nelem)
c
      integer i,ie,nact1,nact2,mact,ne,nne
c
      do i=1,nelem
         rmhnel(i) = 0
         rmhlen(i) = 0
      enddo
      rmhnel(1) = 1
      rmhlen(1) = 1
c
      nact1 = 1
      nact2 = 1
      do while (nact2.lt.nelem )
c
         mact = 0
         do ie=nact1,nact2
            ne = rmhnel(ie)
            do i=4,6
               nne = abs(nu(i,ne))/8
               if ( nne.ne.0 ) then
                  if ( rmhlen(nne).eq.0 ) then
                     mact = mact + 1
                     rmhnel(nact2+mact) = nne
                     rmhlen(nne) = nact2+mact
                  endif
               endif
            enddo
         enddo
c
         nact1 = nact2 + 1
         nact2 = nact2 + mact
c
      enddo
c
      return
      end
c
c
c***********************************************************************
      subroutine MNDNOD ( nvert,nnode,nelem,nu,coor,rmhlen,
     &                    elmref,nbd,ibdnod,bdsec,seclst,nbound,nside,
     &                    nbrigid,diaa,diab,xpos,
     &                    rmhind,x,y,rmhnec, rmhibd,
     &                    xymid,mbdlst,pwall,bdseg,maxnod,maxnbd,nflg)
c
c     This routine generates midnods for the mesh
c***********************************************************************
      implicit none
      integer nvert,nnode,nelem,nu(6,nelem),
     &        rmhlen(nelem),maxnod,maxnbd,elmref(nelem),
     &        nbd,ibdnod(nbd),bdsec,seclst(bdsec+1),nbound,
     &        nside(nbound), mbdlst(bdsec), nflg(maxnbd), nbrigid,bdseg
      integer rmhind(6,nelem),rmhibd(2*nbd),rmhnec(3,nelem)
      real*8 coor(2,nvert),x(maxnod),y(maxnod),xymid(2,nbd),
     &       diaa(nbrigid),diab(nbrigid),xpos(3,nbrigid),pwall(5,bdseg)
c
      integer i,j,k,n,nod1,nod2,ne,ni
      double precision xcoor,xploc,x1,x2
      integer p3(3)
      data p3/2,3,1/
c
c     copy mesh data
c     --------------
      do i=1,nbd
         rmhibd(2*i-1) = ibdnod(i)
      enddo
      do i=1,nvert
         x(i) = coor(1,i)
         y(i) = coor(2,i)
      enddo
c
c     generate new element reference number
c     -------------------------------------
      do i=1,nelem
         rmhnec(1,i) = elmref(i)
      enddo
      do i=1,nelem
         elmref(rmhlen(i)) = rmhnec(1,i)
      enddo
c
c     generate mid-nodes on boundary and update rmhibd
c     ------------------------------------------------
      nnode = nvert
      do i=1,nbd
         nnode = nnode+1
         nod1  = ibdnod(i)
         x(nnode) = xploc(xymid(1,nod1))
         y(nnode) =       xymid(2,nod1)
         rmhibd(2*i) = nnode
         nflg(nod1)  = nnode
      enddo
c
c     generate other mid-nodes
c     ------------------------
      do i=1,nelem
         ne = rmhlen(i)
         do j=1,3
            rmhind(j,ne) = nu(j,i)
            nod1 = nu(j,i)
            nod2 = nu(p3(j),i)
            n = nu(j+3,i)
            ni = n/8
            if ( ni.eq.0 ) then
               rmhnec(j,ne) = 0
               rmhind(j+3,ne) = nflg(nod1)
               if ( elmref(ne).gt.1 ) then
                  rmhind(j+3,ne) = nflg(nod2)
               endif
            elseif ( ni.lt.0 ) then
               ni = -ni
               rmhnec(j,ne)   = rmhlen(ni)
               rmhind(j+3,ne) = nflg(nod1)
               if ( elmref(ne).gt.1 ) then
                  rmhind(j+3,ne) = nflg(nod2)
               endif
            else
               rmhnec(j,ne) = rmhlen(ni)
               if ( ni.gt.i ) then
                  nnode = nnode+1
                  rmhind(j+3,ne) = nnode
                  x1 = x(nod1)
                  x2 = xcoor(x1,x(nod2))
                  x(nnode) = xploc(0.5*(x1+x2))
                  y(nnode) = 0.5*(y(nod1)+y(nod2))
               else
                  k = n - ni*8
                  rmhind(j+3,ne) = rmhind(k,rmhlen(ni))
               endif
            endif
         enddo
c
        if ( nnode+4.gt.maxnod ) then
           write(*,'(A,i5)') 
     &      'Error from MNDNOD, maxnod is too small, must >',nnode+4
           stop
        endif
c
      enddo
c
c     modify the the location of mid-nodes for a smoothing mesh 
c
c      call midsmth ( x,y, nnode,nelem, rmhind,rmhnec)
c
      nbd = 2*nbd
      do i=1,bdsec
         seclst(i) = 2*seclst(i)-1
      enddo
      seclst(bdsec+1) = nbd+1
c
      return
      end
c
c
c**********************************************************************
      subroutine midsmth ( x,y, nnode,nelem, inod,nec)
c
c     Smoothing of generated mid-nodes -> the centroid of the polygon
c----------------------------------------------------------------------
      implicit none
      integer nnode,nelem,inod(6,nelem),nec(3,nelem)
      real*8 x(nnode),y(nnode)
c
      integer itermx
      real omega,eps      
      integer p3(3)
      data p3/2,3,1/
      integer iter,ne1,ne2,i,i1,n0,n
      real*8 err, x0,y0,dx,dy,xcoor,ymin,ymax,yl
c
      itermx = 10
      omega  = 1.0
      eps    = 0.002
c
      ymin = 1.d10
      ymax = -1.d10
      do i=1,nnode
         ymin = min(ymin,y(i))
         ymax = max(ymax,y(i))
      enddo
      yl = ymax-ymin
c
c     update nodes to the centroid of the polygon
c     -------------------------------------------
      do iter=1,itermx
         err = 0
         do ne1 = 1,nelem
            do i=1,3
               ne2 = nec(i,ne1)
               if ( ne2.gt.0 ) then
                  n0 = inod(i+3,ne1)
                  x0 = x(n0)
                  y0 = y(n0)
                  n  = inod(i,ne1)
                  dx = -8.d0*x0 + xcoor(x0,x(n))
                  dy = -8.d0*y0 + y(n)
                  n  = inod(p3(i),ne1)
                  dx = dx + xcoor(x0,x(n))
                  dy = dy + y(n)
                  do i1=4,6
                     n  = inod(i1,ne1)
                     dx = dx + xcoor(x0,x(n))
                     dy = dy + y(n)
                     n =  inod(i1,ne2)
                     dx = dx + xcoor(x0,x(n))
                     dy = dy + y(n)
                  enddo
                  dx = dx/6.d0
                  dy = dy/6.d0
                  x(n0) = x0 + omega*dx
                  y(n0) = y0 + omega*dy
                  err = max(dabs(dx),dabs(dy),err)
               endif
            enddo
         enddo
         if ( err.le.eps*yl ) return
      enddo
c
      return
      end
c
c
c***********************************************************************
      subroutine RECTMSH (sidea,sideb,ncr,ceta, xycoor)
c
c     generating grids on surface of a rectangular particles
c
c   INPUT:
c       sidea, sideb = a and b of the the rectangular particle
c       ncr     : number of segments on particle surface
c       ceta    = orientation of the particle
c
c   OUTPUT:
c       xycoor = x & y coordinates of the grid
c***********************************************************************
      implicit none
      integer ncr
      real*8 sidea,sideb,xycoor(2,2*(ncr+4)), ceta
c
c     local variables
      integer nx,ny,nd,i
      real*8 dx,dy,cs,si
c
      nx = 2*int(ncr*(0.5*sidea/(sidea+sideb)))
      ny = 2*int(ncr*(0.5*sideb/(sidea+sideb)))
      ncr = 2*(nx+ny)
      dx = 2.d0*sidea/nx
      dy = 2.d0*sideb/ny
      nd = 0
c
c     along side 1
      do i=1,nx
        nd = nd+1
        xycoor(1,nd) = -sidea + (i-1)*dx
        xycoor(2,nd) =  sideb
      enddo
c
c     along side 2
      do i=1,ny
        nd = nd+1
        xycoor(1,nd) = sidea
        xycoor(2,nd) = sideb - (i-1)*dy
      enddo
c
c     along side 3
      do i=1,nx
        nd = nd+1
        xycoor(1,nd) =  sidea - (i-1)*dx
        xycoor(2,nd) = -sideb
      enddo
c
c     along side 4
      do i=1,ny
        nd = nd+1
        xycoor(1,nd) = -sidea
        xycoor(2,nd) = -sideb + (i-1)*dy
      enddo
c
c     rotate the coordinates
      cs = dcos(ceta)
      si = dsin(ceta)
      do i=1,ncr
        dx          = xycoor(1,i)*cs - xycoor(2,i)*si
        xycoor(2,i) = xycoor(1,i)*si + xycoor(2,i)*cs
        xycoor(1,i) = dx
      enddo
c
      return
      end
