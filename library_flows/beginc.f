c***********************************************************************
      subroutine beginc ( icase,iax,kaxi,neq,iou, x,y, z,zs,
     &                    trans,tsc,tsc1, shg1,shg1x,we1,
     &                    SLVvelo,SLVtemp,SLVelas,SLVpres, 
     &                    ndggeom,ndgvelo,ndgtemp,ndgelas,ndgpres,
     &                    nrdgeom,indgu,indgp,indgtm,indgr,
     &                    nbd,ibdnod,nic,ic,nbound,nside, ipnode,
     &                    ibdu,ibdt,ibds,ibdp,densu,denst,denss,densp,
     &                    ncpu,ncps )
c
c     this routine prepares the boundary conditions
c
c     INPUT:
c           icase = index for the problem
c           iax, kaxi = geometry of the problem
c           neq = number of equations for the problem
c           iou = unit number for write
c           trans = logic for transient problem
c           ndb = number of bounrary nodes
c           ibdnod(nbd) = index for boundary nodes
c           nic, ic(nic) = index for boundary segements
c           nbound = number of closed boundaries
c           nside(nbound) = number of boundary segements for each 
c                           closed boundaries
c           ibdu(ncpu,nbd) = boundary condition index for velocities
c           ibdt(nbd) = boundary condition index for temperature
c           ibds(nbd) = boundary condition index for stress
c           ibdp(nbd) = boundary condition index for pressure 
c                       (when ipnode=-1)
c           SLVvelo,SLVtemp,SLVelas = logics for each variables
c           ndggeom,ndgvelo,ndgtemp,ndgelas = number of nodes for 
c                                             each variables
c           nrdgeom = orders of interpolations for the mesh
c           indgu,indgtm,indgr = location index for variables
c           x(ndggeom),y(ndggeom) = coordinates for the mesh
c           z(neq) = variables
c           zs(neq) = history of variables
c           tsc,tsc1 = time step
c           densu(ncpu,nbd) = boundary velocity
c           denst(nbd) = boundary temperature
c           denss(ncps,nbd) = boundary stress components
c           densp(ncps,nbd) = boundary pressure
c           ncpu = maximum number of volocity components
c           ncps = maximum number of stress components
c
c     10/16/1996 by Howard H. Hu
c***********************************************************************
      implicit none
      integer iou,icase,iax,kaxi,nbd,neq,nic,ic(nic),ibdnod(nbd),
     &        nbound,nside(nbound),ipnode,
     &        ncpu,ncps,ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),ibdp(nbd),
     &        ndggeom,ndgvelo,ndgtemp,ndgelas,ndgpres,nrdgeom,
     &        indgu,indgp,indgtm,indgr, maxrwk
      real*8 z(neq),zs(neq),tsc,tsc1,x(ndggeom),y(ndggeom),
     &       densu(ncpu,nbd),denst(nbd),denss(ncps,nbd),densp(nbd)
      real*8 shg1(3,3,2),shg1x(3,3,2),we1(3,2)
      logical SLVvelo,SLVtemp,SLVelas,SLVpres,trans
c
      integer i
c
c      call prebnd(icase,kaxi, nrdgeom,ndggeom,x,y, shg1,shg1x,we1,
c     &            nbd,ibdnod,nic,ic,nbound,nside,
c     &            ibdu,ibdt, densu,denst,ncpu )
c
      if (.not.trans ) then
c
c       steady scheme
c       -------------
        tsc = 0.d0
        tsc1 = 0.d0
        do i=1,neq
          zs(i) = 0.d0
        enddo
c
      else
c
c       transient scheme
c       ----------------
        call impess ( iax,z, nbd,ibdnod, ipnode,
     &                ibdu,ibdt,ibds,ibdp,densu,denst,denss,densp,
     &                ncpu,ncps,indgu,indgp,indgtm,indgr, 
     &                ndgvelo,ndgtemp,ndgelas,ndgpres,
     &                SLVvelo,SLVtemp,SLVelas,SLVpres)
c
      endif
c
      return
      end
c
c
c***********************************************************************
c      subroutine prebnd(icase,kaxi, nrdgeom,ndggeom,x,y, shg1,shg1x,we1,
c     &                  nbd,ibdnod,nic,ic,nbound,nside,
c     &                  ibdu,ibdt, densu,denst,ncpu )
c
c     this routine calculates nodal forces on boundary points
c***********************************************************************
c      implicit none
c      integer icase,kaxi,nrdgeom, nbd,ibdnod(nbd),nic,ic(nic+1),nbound,
c     &       nside(nbound),ibdu(ncpu,nbd),ibdt(nbd), ndggeom, ncpu
c      real*8 x(ndggeom),y(ndggeom), densu(ncpu,nbd),denst(nbd)
c      real*8 shg1(3,3,2),shg1x(3,3,2),we1(3,2)
c
c      integer nidx,nd,itype,nb,kk,i,j,k,iside,ica,icb,ica1,nel,m,ib,
c     &        ibm1,iaside,n
c      real*8 xn(3),yn(3),b1dx(3,3),b1dy(3,3),b1ds(3,3), xcoor
c
c      nidx = nrdgeom
c      nd = nidx+1
c
c     1. initialization
c     -----------------
c      do iside=1,nic
c        ica = ic(iside)
c        icb = ic(iside+1)-1
c        do i=ica,icb
c          itype = ibd(i)
c          if ( itype.eq.2 ) then
c            valx(i) = 0.d0
c            valy(i) = 0.d0
c          elseif ( itype.eq.3 ) then
c           mixed : un,fs
c            valy(i) = 0.d0
c          elseif ( itype.eq.4 ) then
c           mixed : fn,us
c            valx(i) = 0.d0
c          endif
c          if ( ibdu(3,i).eq.0 ) valw(i) = 0.d0
c          if ( ibdt(i).eq.0 ) valtmp(i) = 0.d0
c        enddo
c      enddo
c
c     2. calculation of nodal forces on each side
c     -------------------------------------------
c      iside = 0
c      do nb=1,nbound
c        ica1 = ic(iside+1)
c        do kk = 1,nside(nb)
c          iside = iside + 1
c          ica = ic(iside)
c          icb = ic(iside + 1)
c          nel = (icb-ica)/nidx
c
c          do 22 m=1,nel
c            ib = ica + nidx*(m-1)
c            ibm1 = ib - 1
c            iaside = ib + iside - 2
c
c            n = ibdnod(ib)
c            xn(1) = x(n)
c            yn(1) = y(n)
c            do i=2,nd
c              k = ibm1 + i
c              if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
c              n = ibdnod(k)
c              xn(i) = xcoor(xn(1),x(n))
c              yn(i) = y(n)
c            enddo
cc
c            call mat1ds (xn,yn,nidx,kaxi,b1dx,b1dy,b1ds,shg1,shg1x,we1)
cc
c            do 22 i=1,nd
c              k = ibm1 + i
c              if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
c              do 22 j=1,nd
c              n = iaside + j
cc
c              if ( ibd(k).eq.2 ) then
c                valx(k) = valx(k) + b1dy(j,i)*densn(n)
c     &                            + b1dx(j,i)*denss(n)
c                valy(k) = valy(k) - b1dx(j,i)*densn(n)
c     &                            + b1dy(j,i)*denss(n)
c              elseif ( ibd(k).eq.3 ) then
c                valy(k) = valy(k) + b1ds(j,i)*denss(n)
c              elseif ( ibd(k).eq.4 ) then
c                valx(k) = valx(k) + b1ds(j,i)*densn(n)
c              endif
c
c              if ( ibdw(k).eq.2 ) 
c     &                           valw(k) = valw(k) + b1ds(j,i)*densw(n)
c              if ( ibdt(k).eq.2 ) 
c     &                       valtmp(k) = valtmp(k) - b1ds(j,i)*denst(n)
c
c  22      continue
c
c        enddo
c      enddo
cc
c      return
c      end
c
c
c***********************************************************************
      subroutine impess ( iax,z, nbd,ibdnod, ipnode,
     &                    ibdu,ibdt,ibds,ibdp,densu,denst,denss,densp,
     &                    ncpu,ncps,indgu,indgp,indgtm,indgr, 
     &                    ndgvelo,ndgtemp,ndgelas,ndgpres,
     &                    SLVvelo,SLVtemp,SLVelas,SLVpres)
c
c     this routine is imposing the equality between the essential
c       boundary conditions and the predicted values on the boundaries
c***********************************************************************
      implicit none
      integer ndgvelo,ndgtemp,ndgelas,ndgpres,indgu,indgp,indgtm,indgr
      integer ncpu,ncps, ipnode
      integer iax,nbd,ibdnod(nbd),ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),
     &        ibdp(nbd)
      real*8 densu(ncpu,nbd),denst(nbd),denss(ncps,nbd),densp(nbd),z(*)
      logical SLVvelo,SLVtemp,SLVelas,SLVpres
c
      integer i,k
c
c     imposition of the elastic stresses
c     ----------------------------------
      if ( SLVelas ) then
        do i=1,nbd
          k = ibdnod(i)
          if ( ibds(i).eq.1 .and. k.le.ndgelas ) then
            z(indgr+k)           = denss(1,i)
            z(indgr+ndgelas+k)   = denss(2,i)
            z(indgr+2*ndgelas+k) = denss(3,i)
            if ( iax.eq.1 ) z(indgr+3*ndgelas+k) = denss(4,i)
          endif
        enddo
      endif
c
c     imposition of the velocities
c     ----------------------------
      if ( SLVvelo ) then
        do i=1,nbd
          k = ibdnod(i)
          if ( (ibdu(1,i).eq.1.or.ibdu(1,i).lt.0) .and. k.le.ndgvelo) 
     &       z(indgu+k) = densu(1,i)
          if ( (ibdu(2,i).eq.1.or.ibdu(2,i).lt.0) .and. k.le.ndgvelo) 
     &       z(indgu+ndgvelo+k) = densu(2,i)
          if ( iax.eq.2 .and. (ibdu(3,i).eq.1.and.k.le.ndgvelo) )
     &       z(indgu+2*ndgvelo+k) = densu(3,i)
        enddo
      endif
c
c     imposition of the pressure
c     --------------------------
      if(SLVpres)then
      if ( ipnode.gt.0 ) then
        k = ibdnod(ipnode)
        if(k<=ndgpres)then
          z(indgp+k) = 0.d0
        else
          write(*,*)'subroutine impess: invalid ipnode=',ipnode
          stop
        endif
      elseif(ipnode.eq.-1)then
        do i=1,nbd
          k = ibdnod(i)
          if(ibdp(i).eq.1.and.k.le.ndgpres) z(indgp+k)=densp(i)
        enddo
      endif
      endif
c
c     imposition of the temperature
c     -----------------------------
      if ( SLVtemp ) then
        do i=1,nbd
          k = ibdnod(i) 
          if ( ibdt(i).eq.1 .and. k.le.ndgtemp ) z(indgtm+k)=denst(i)
        enddo
      endif
c
      return
      end
