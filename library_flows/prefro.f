c***********************************************************************
      subroutine prefro(neq,neqA,neqB, iout, nelem,inod,ngmx,ishape,
     &                  nbd,ibdnod,nic,ic, ibdu,ibdt,ibds,ibdp,ncpu,
     &                  ipnode,
     &                  nbrigid,ncpb,nveli,ia,ja,ib,jb,neltA,neltB,
     &                  isx,jsx,neltS,ioptp, imp,iperm,jperm,matglb,
     &                  iwork,iwk,jwk, iglo,ntot,
     &                  indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                  ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                  nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &                  ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi)
c
c     this subroutine is a preprocessor preparing the matrix
c
c     ncpxxxx - number of components for each group of variable.
c     nrdxxxx - order of interpolation for each group of variable.
c     ndgxxxx - number of gobal nodes for each variable.
c        xxxx = geom, velo, pres, temp, elas, strm
c
c     NOTE:
c       saddle-point problem A & B:  | A  B^t |
c                                    | B  0   |
c       for global matrix system B=0
c
c     called by front
c
c     rewritten 6/23/99  by H.H. Hu
c***********************************************************************
      implicit none
      integer neq,neqA,neqB,nbd,nelem,ngmx,ntot,ncpu,ncpb,neltA,neltB,
     &        nbrigid,inod(ngmx,nelem),ibdnod(nbd),ipnode,nveli(2),iout,
     &        ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),ibdp(nbd),ishape(nelem)
     &        ,ia(neqA+1),ja(neltA),imp(neq),iperm(neq),jperm(neq),
     &        nic,ic(nic+1), 
     &        iwk(neq+1),jwk(neq),iwork(ntot*nelem),
     &        iglo(ntot),ib(neqB+1),jb(neltB)
      integer neltS,isx(neqB+1),jsx(neltS),ioptp
      integer ncpvelo,ncppres,ncptemp,ncpelas,ncpphi
      integer nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi
      integer ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi
      integer indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi
      logical matglb
c
      integer ndloc,ishm,ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,nvel
      integer ndlgeom,ne,i,l,k,nd,n,np,j,m,nv,nwall,
     &        m0,mi,n1,n2,neqt,nvel1, mb,ms,mb0,ms0
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(PRE_MAT_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     calculation of nveli
c     --------------------
      nveli(1) = ndloc(0,nrdelas)*ncpelas + ndloc(0,nrdvelo)*ncpvelo +
     &           ndloc(0,nrdpres)*ncppres + ndloc(0,nrdtemp)*ncptemp +
     &           ndloc(0,nrdphi) *ncpphi*2
      nveli(2) = ndloc(1,nrdelas)*ncpelas + ndloc(1,nrdvelo)*ncpvelo +
     &           ndloc(1,nrdpres)*ncppres + ndloc(1,nrdtemp)*ncptemp +
     &           ndloc(1,nrdphi) *ncpphi*2
c
c     computation of iglo(local var#,element#) = global var#
c     ------------------------------------------------------
      neqt = neqA+neqB
      ishm = ishape(1)
      ndlgeom = ndloc(ishm, nrdgeom)
      ndlelas = ndloc(ishm, nrdelas)
      ndlvelo = ndloc(ishm, nrdvelo)
      ndlpres = ndloc(ishm, nrdpres)
      ndltemp = ndloc(ishm, nrdtemp)
      ndlphi  = ndloc(ishm, nrdphi )
      nvel = nveli(ishm+1)
c
c     generate boundary condition index: imp
c     --------------------------------------
      call bcimp ( imp,neq,nbd,ibdnod,ibdu,ibdt,ibds,ibdp,ipnode,ncpu,
     &             indgu,indgp,indgtm,indgr,indgnb,
     &             ncpvelo,ncppres,ncptemp,ncpelas,
     &             ndgvelo,ndgpres,ndgtemp,ndgelas )
c
c     computation of the order of elimination ( permutation idex)
c     -----------------------------------------------------------
c     calculate iwk(global var#) = last element# related
      do i=1,neq
         jwk(i) = 0
      enddo
c     debug
c      write(*,*)'elas',ndgelas,ncpelas
c      write(*,*)'velo',ndgvelo,ncpvelo
c      write(*,*)'temp',ndgtemp,ncptemp
c      write(*,*)'pres',ndgpres,ncppres
c      write(*,*)'phi',ndgphi,ncpphi
c      write(*,*)'neq',neq
c
      do ne=1,nelem

         call FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)

         do j=1,nvel
            l = iglo(j)
            iwk(l) = ne
c     debug
c      write(*,*)'***',ne,l,neq
c
            jwk(l) = jwk(l) + 1
         enddo
      enddo
c
c     iperm(global var#) = order of elimination
c     jperm(order of elimination) = global var#
c     -----------------------------------------
      do i=1,neq
         iperm(i) = i
         jperm(i) = i
      enddo
      n = ncpb*nbrigid
      m = neqt
c
      if ( matglb ) then
c
c        global matrix
         nvel1 = nvel
         do ne=1,nelem
            call FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)
            do j=1,nvel
               l = iglo(j)
               if ( iwk(l).eq.ne ) then
                  if ( imp(l).eq.0 ) then
                     n = n + 1
                     iperm(l) = n
                     jperm(n) = l
                  else
                     m = m + 1
                     iperm(l) = m
                     jperm(m) = l
                  endif
               endif
            enddo
         enddo
c
      else
c
c        separate matrices A & B
         np = neqA
c        n1 = neq-ndgpres
         n1 = indgp+ndgpres*ncppres
cdebug
c      write(*,*)"indgp,ndgpres,ncppres",indgp,ndgpres,ncppres
c      
         nvel1 = nvel-ncppres*ndlpres
         do ne=1,nelem
            call FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)

            do j=1,nvel
               l = iglo(j)
               if ( iwk(l).eq.ne) then
                  if ( imp(l).eq.0 ) then
                     if ( l.gt.indgp.and.l.le.n1) then
                        np = np + 1
                        iperm(l)  = np
                        jperm(np) = l
                     else
                        n = n + 1
                        iperm(l) = n
                        jperm(n) = l
                     endif
                  else
                     m = m + 1
                     iperm(l) = m
                     jperm(m) = l
                  endif
               endif
            enddo
         enddo
c
      endif
c
c     original ordering
c     -----------------
c      do i=1,neq
c         iperm(i) = i
c         jperm(i) = i
c      enddo
c
c      do i=1,neq
c         write(15,*) i,iperm(i),jperm(i)
c      enddo
cdab
c
c     preparation of storage of the matrix in the CSR format
c     ------------------------------------------------------
      iwk(1) = 0
      do i=1,neq
         iwk(i+1)= iwk(i) + jwk(i)
      enddo

      do ne=1,nelem
            call FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)

         do j=1,nvel
            i = iglo(j)
            iwk(i) = iwk(i) + 1
            iwork(iwk(i)) = ne
         enddo
      enddo

      do i=neq+1,2,-1
         iwk(i) = iwk(i-1)
      enddo
      iwk(1) = 0
c
c      write(25,'(10i6)') (iwk(i),i=1,neq+1)
c      write(25,'(12i6)') (iwork(i),i=1,ntot*nelem)
c
c     equations for velocities of moving particles
c     note: diagonal element is at the first
c     --------------------------------------------
      do i=1,neq
         jwk(i) = 0
      enddo
      m = 0
      ia(1) = 1
      nwall = nic-nbrigid
      do np=1,nbrigid
         do k=1,ncpb
            nv = indgnb + ncpb*(np-1)+k
            ja(m+1) = nv
            m = m + 1
         enddo
         nv = indgnb + ncpb*(np-1)+1
         jwk(nv) = nv
         m0 = m+1
         do i=ic(nwall+np),ic(nwall+np+1)-1
            nd = indgu+ibdnod(i)
            do j=iwk(nd)+1,iwk(nd+1)
               ne = iwork(j)
            call FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)
               do k=1,nvel1
                  n1 = iglo(k)
                  n = iperm(n1)
                  if ( jwk(n).ne.nv ) then
                     jwk(n) = nv
                     if ( imp(n1).eq.0 ) then
                        mi = m
                        do while ( ja(mi).gt.n .and. mi.ge.m0 )
                           ja(mi+1) = ja(mi)
                           mi = mi - 1
                        enddo
                        ja(mi+1) = n
                        m = m + 1
c
                     elseif ( imp(n1).lt.0 ) then
                        n2 = indgnb + ncpb*(abs(imp(n1))-1)+1
                        if ( jwk(n2).ne.nv ) then
                           jwk(n2) = nv
                           mi = m
                           do while ( ja(mi).gt.n2 .and. mi.ge.m0 )
                              ja(mi+ncpb) = ja(mi)
                              mi = mi - 1
                           enddo
                           do l=1,ncpb
                              ja(mi+l) = n2+l-1
                              m = m + 1
                           enddo
                        endif
c
                     endif
                  endif
               enddo
             enddo
         enddo
         ia(nv+1) = m+1
         n = ia(nv+1) - ia(nv)
c
         do k=2,ncpb
            do j=1,n
               ja(m+j) = ja(m-n+j)
            enddo
            n1 = ja(m+1)
            ja(m+1) = nv + k-1
            ja(m+k) = n1
            m = m+n
            ia(nv+k) = m+1
         enddo
c
      enddo
c
      do i=ncpb*nbrigid+1,neqA
         nd = jperm(i)
         m = m + 1
         ja(m) = i
         jwk(i) = i
         m0 = m+1
         do j=iwk(nd)+1,iwk(nd+1)
            ne = iwork(j)
            call FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)
            do k=1,nvel1
               n1 = iglo(k)
               n = iperm(n1)
               if ( jwk(n).ne.i ) then
                  jwk(n) = i
                  if ( imp(n1).eq.0 ) then
                     mi = m
                     do while ( ja(mi).gt.n .and. mi.ge.m0 )
                        ja(mi+1) = ja(mi)
                        mi = mi - 1
                     enddo
                     ja(mi+1) = n
                     m = m + 1
c
                  elseif ( imp(n1).lt.0 ) then
                     n2 = indgnb + ncpb*(abs(imp(n1))-1)+1
                     if ( jwk(n2).ne.i ) then
                        jwk(n2) = i
                        mi = m
                        do while ( ja(mi).gt.n2 .and. mi.ge.m0 )
                           ja(mi+ncpb) = ja(mi)
                           mi = mi - 1
                        enddo
                        do l=1,ncpb
                           ja(mi+l) = n2+l-1
                           m = m + 1
                        enddo
                     endif
c
                  endif
               endif
            enddo
         enddo
         ia(i+1) = m+1
c
      enddo
c
      mb = 0
      ib(1) = 1
      ms = 0
      isx(1) = 1

      do i=neqA+1,neqt
         nd = jperm(i)
         mb0 = mb+1

         if ( ioptp/10.eq.2 ) then
            ms = ms + 1
            jsx(ms) = i-neqA
            jwk(i) = i
            ms0 = ms+1
         endif

         do j=iwk(nd)+1,iwk(nd+1)
            ne = iwork(j)
            call FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)

            do k=1,nvel1
               n1 = iglo(k)
               n  = iperm(n1)
               if ( jwk(n).ne.i ) then
                  jwk(n) = i
                  if ( imp(n1).eq.0 ) then
                     mi = mb
                     do while ( mi.ge.mb0 .and. jb(mi).gt.n )
                        jb(mi+1) = jb(mi)
                        mi = mi - 1
                     enddo
                     jb(mi+1) = n
                     mb = mb + 1
                  elseif ( imp(n1).lt.0 ) then
                     n2 = indgnb + ncpb*(abs(imp(n1))-1)+1
                     if ( jwk(n2).ne.i ) then
                        jwk(n2) = i
                        mi = mb
                        do while ( mi.ge.mb0 .and. jb(mi).gt.n2 )
                           jb(mi+ncpb) = jb(mi)
                           mi = mi - 1
                        enddo
                        do l=1,ncpb
                           jb(mi+l) = n2+l-1
                           mb = mb + 1
                        enddo
                     endif
                  endif
               endif
            enddo
c
            if ( ioptp/10.eq.2 ) then
               do k=nvel1+1,nvel
                  n1 = iglo(k)
                  n  = iperm(n1)
                  if ( jwk(n).ne.i ) then
                     jwk(n) = i
                     if ( imp(n1).eq.0 ) then
                        mi = ms
                        do while ( jsx(mi).gt.n-neqA .and. mi.ge.ms0 )
                           jsx(mi+1) = jsx(mi)
                           mi = mi - 1
                        enddo
                        jsx(mi+1) = n-neqA
                        ms = ms + 1
                     endif
                  endif
               enddo
            endif

         enddo

         ib(i-neqA+1)  = mb+1
         isx(i-neqA+1) = ms+1

      enddo
c
      if ( ioptp/10.eq.3 ) then
         do i=1,neqA
            jwk(i) = 0
         enddo
         do i=1,neqB
           ms = ms + 1
           jsx(ms) = i
           do k=ib(i),ib(i+1)-1
              jwk(jb(k)) = i
           enddo

           do 100 j=1,neqB
              if ( j.ne.i ) then
                 do k=ib(j),ib(j+1)-1
                    if ( jwk(jb(k)).eq.i ) then
                       ms = ms + 1
                       jsx(ms) = j
                       go to 100
                    endif
                 enddo
              endif
100        continue

           isx(i+1) = ms+1
         enddo
      endif
c
c     check the length of the matrices
c     --------------------------------
      m = ia(neqA+1) - 1
      if ( m.ne.neltA ) then
         write(iout,'(/A,i9)') 'The length of the A matrix is',m
         if ( m.gt.neltA ) stop
      endif
      m = ib(neqB+1) - 1
      if ( m.ne.neltB ) then
         write(iout,'(/A,i9)') 'The length of the B matrix is',m
         if ( m.gt.neltB ) stop
      endif
      m = isx(neqB+1) - 1
      if ( m.ne.neltS ) then
         write(iout,'(/A,i9)') 'The length of the S matrix is',m
c         neltS = m
         if ( m.gt.neltS ) stop
      endif
c
c     permute the index for boundary conditions
c     -----------------------------------------
      do i=1,neq
         iwk(i) = imp(i)
      enddo
      do i=1,neq
         j=iperm(i)
         imp(j) = iwk(i)
      enddo
c
cC     Support for Petsc timing
c      call PLogEventEnd(PRE_MAT_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c***********************************************************************
      subroutine bcimp( imp,neq,nbd,ibdnod,ibdu,ibdt,ibds,ibdp,ipnode,
     &                  ncpu,indgu,indgp,indgtm,indgr,indgnb,
     &                  ncpvelo,ncppres,ncptemp,ncpelas,
     &                  ndgvelo,ndgpres,ndgtemp,ndgelas)
c
c     generate boundary condition index: imp
c
c     ncpxxxx - number of components for each group of variable.
c     ndgxxxx - number of gobal nodes for each variable.
c        xxxx = velo, pres, temp, elas, strm
c***********************************************************************
      implicit none
      integer neq,imp(neq),nbd,ibdnod(nbd),ipnode,
     &        ncpu,ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),ibdp(nbd)
      integer indgu,indgp,indgtm,indgr,indgnb
      integer ncpvelo,ncppres,ncptemp,ncpelas
      integer ndgvelo,ndgpres,ndgtemp,ndgelas
c
      integer i,l,nd
c
      do i=1,neq
         imp(i) = 0
      enddo
c
      do l=1,ncpelas
         do i=1,nbd
            nd = ibdnod(i)
            if ( nd.le.ndgelas ) imp(indgr+(l-1)*ndgelas+nd)=ibds(i)
         enddo
      enddo
c
      do l=1,ncpvelo
         do i=1,nbd
            nd = ibdnod(i)
            if ( nd.le.ndgvelo ) imp(indgu+(l-1)*ndgvelo+nd)=ibdu(l,i)
         enddo
      enddo
c
      do l=1,ncppres

         if ( ipnode.gt.0 ) then
            imp(indgp+ibdnod(ipnode)) = 1
         elseif(ipnode.eq.-1)then
            do i=1,nbd
               nd=ibdnod(i)
               if(nd.le.ndgpres)imp(indgp+nd) = ibdp(i)
            enddo
         endif
      enddo
c
      do l=1,ncptemp
         do i=1,nbd
            nd = ibdnod(i)
            if ( nd.le.ndgtemp ) imp(indgtm+nd) = ibdt(i)
         enddo
      enddo
c
      return
      end
c
c
c***********************************************************************
      subroutine FORMiglo ( iglo,nvel, ne,inod,ngmx,nrdpres,
     &                      ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                      ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                      ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                      indgr,indgu,indgp,indgtm,indgphi,indgpsi)
c
      implicit none
      integer nvel,iglo(nvel),ne,ngmx,inod(ngmx,*)
      integer ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &        ndgelas,ndgvelo,ndgpres,ndgtemp,
     &        ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,nrdpres,
     &        indgr,indgu,indgp,indgtm,indgphi,indgpsi
      
      integer nnl,ike,l,k
c
c       elastic-stress components
        nnl = 0
        do l=1,ncpelas
           do k=1,ndlelas
              ike = nnl + (l-1)*ndlelas + k
              iglo(ike) = indgr + (l-1)*ndgelas + inod(k,ne)
           enddo
        enddo
c
c       velocity components
        nnl = nnl + ncpelas*ndlelas
        do l=1,ncpvelo
           do k=1,ndlvelo
              ike = nnl + (l-1)*ndlvelo + k
              iglo(ike) = indgu + (l-1)*ndgvelo + inod(k,ne)
           enddo
        enddo
c
c       temperature component
        nnl = nnl + ncpvelo*ndlvelo
        do l=1,ncptemp
           do k=1,ndltemp
              ike = nnl + (l-1)*ndltemp + k
              iglo(ike) = indgtm + (l-1)*ndgtemp + inod(k,ne)
           enddo
        enddo

c
c       phase-field variables
        nnl = nnl + ncptemp*ndltemp
        if(ncpphi.eq.1)then
           do k=1,ndlphi
              iglo(nnl+k) = indgphi + inod(k,ne)
              iglo(nnl+ndlphi+k) = indgpsi + inod(k,ne)
           enddo
        endif
c
c       pressure component
        nnl = nnl + 2*ncpphi*ndlphi
        do l=1,ncppres
           do k=1,ndlpres
              ike = nnl + (l-1)*ndlpres + k
              if ( nrdpres.eq.0 ) then
                 iglo(ike) = indgp + (l-1)*ndgpres + ne
              else
                 iglo(ike) = indgp + (l-1)*ndgpres + inod(k,ne)
              endif
           enddo
        enddo
c
      return
      end
