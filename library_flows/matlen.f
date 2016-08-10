c***********************************************************************
      subroutine matlen1 ( nelt, nelem,inod,ngmx,ishape, norder,
     &                     iwk,jwk,iwork )
c
c     this subrourine calculates matrix size for node based variables
c     output: nelt (number of non-zero elements in the sparse matrix)
c     the length of iwk should >= ndggeom (or nnode) +1
c     the length of jwk should >= ndggeom (or nnode)
c***********************************************************************
      implicit none
      integer nelt,nelem,ishape(nelem),ngmx,inod(ngmx,nelem),norder,
     &        iwk(*),jwk(*),iwork(ngmx*nelem)
c
      integer i,j,k,n,ne,ndloc,ndglb,ndl,ndg,ishm, is
c
      ndg = ndglb(norder)
      ishm = ishape(1)
      ndl = ndloc(ishm,norder)
c
      do i=1,ndg+1
         iwk(i) = 0
      enddo
      do i=1,nelem
         do j=1,ndl
            is = inod(j,i)+1
            iwk(is) = iwk(is) + 1
         enddo
      enddo
      do i=2,ndg+1
         iwk(i)= iwk(i-1) + iwk(i)
      enddo

      do i=1,nelem
         do j=1,ndl
            is = inod(j,i)
            iwk(is) = iwk(is) + 1
c            write(*,'(a,5(1x,i5))')'&&',i,is,iwk(is),ngmx,nelem
            iwork(iwk(is)) = i
         enddo
      enddo
      do i=ndg+1,2,-1
         iwk(i) = iwk(i-1)
      enddo
      iwk(1) = 0
c
      do i=1,ndg
         jwk(i) = 0
      enddo
      nelt = 0
      do i=1,ndg
         do j=iwk(i)+1,iwk(i+1)
            ne = iwork(j)
            do k=1,ndl
               n = inod(k,ne)
               if ( jwk(n).ne.i ) then
                  nelt = nelt + 1
                  jwk(n) = i
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
      subroutine matlen2 ( neltA,neltB,neltS,neq,neqA,neqB,nelem,inod,
     &                 ishape,nbd,ibdnod,nic,ic, ibdu,ibdt,ibds,ibdp,
     &                 ipnode,nbrigid,imp,iwork,iwk,jwk,jwk1,
     &                 ngmx,ncpu,ncpb, matglb,
     &                 indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                 ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                 nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &                 ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi)
c
c     this subroutine calculates storage space for the matrices
c
c     OUTPUT:
c     neltA: matrix size for A
c     neltB: matrix size for B - saddle point problem
c     neltS: matrix size for S - approximate Schur complement
c     neqA: number of equations in A
c     neqB: number of equations in B
c     imp: boundary condition indx for the global variables
c     
c     ncpxxxx - number of components for each group of variable.
c     nrdxxxx - order of interpolation for each group of variable.
c     ndgxxxx - number of gobal nodes for each variable.
c        xxxx = geom, velo, pres, temp, elas, strm
c    
c     NOTE:
c       saddle-point problem :  | A  B^t |
c                               | B  0   |
c       for global matrix system B=0
c     
c     Modified by Pengtao Yue
c        Feb 25, 2005,  modified some array bounds, the size of iwk 
c           should be nnode+1
c***********************************************************************
      implicit none
      integer nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi
      integer ncpvelo,ncppres,ncptemp,ncpelas,ncpphi
      integer ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi
      integer indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi
      integer neltA,neltB,neltS,neq,neqA,neqB,nbd,nelem,ishape(nelem),
     &        nbrigid,ngmx,ncpu,inod(ngmx,nelem),ibdnod(nbd),ipnode,
     &        ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),ibdp(nbd),imp(neq),
     &        nic,ic(nic+1),ncpb,iwk(ndggeom+1),
     &        jwk(ndggeom+nbrigid),jwk1(nelem),iwork(ngmx*nelem)
      logical matglb
c
      integer ndloc,ndglb,ishm,nvelo,npres,nelas,ntemp,nphi,npsi
      integer ne,i,j,l,k,nd,n,n1,np,nv,nwall,ndl,ndg,is
c
c     generate boundary condition index: imp
c     --------------------------------------
      call bcimp ( imp,neq,nbd,ibdnod,ibdu,ibdt,ibds,ibdp,ipnode,ncpu,
     &             indgu,indgp,indgtm,indgr,indgnb,
     &             ncpvelo,ncppres,ncptemp,ncpelas,
     &             ndgvelo,ndgpres,ndgtemp,ndgelas )
c
c     compute the linked elements for each node
c     ----------------------------------------- 
      ishm = ishape(1)
      ndg  = ndglb(2)
      ndl  = ndloc(ishm,2)
c
      do i=1,ndg+1
         iwk(i) = 0
      enddo
      do i=1,nelem
         do j=1,ndl
            is = inod(j,i)+1
            iwk(is) = iwk(is) + 1
         enddo
      enddo
      do i=2,ndg+1
         iwk(i)= iwk(i-1) + iwk(i)
      enddo

      do i=1,nelem
         do j=1,ndl
            is = inod(j,i)
            iwk(is) = iwk(is) + 1
            iwork(iwk(is)) = i
         enddo
      enddo
      do i=ndg+1,2,-1
         iwk(i) = iwk(i-1)
      enddo
      iwk(1) = 0
c
c     calculate the length for the global matrix
c     ------------------------------------------
      do i=1,ndg+nbrigid
         jwk(i) = 0
      enddo
      do i=1,nelem
         jwk1(i) = 0
      enddo
c      do i=1,iwk(ndg+1)
c         if(iwork(i).eq.0)then
c         write(*,*)'i,iwork(i)',i
c         endif
c      enddo

      neltA = 0
      neltB = 0
      neltS = 0
      neqA  = 0
      neqB  = 0
      do i=1,ndg
c
         nvelo = 0
         npres = 0
         nelas = 0
         ntemp = 0
         nphi  = 0
         npsi  = 0
         do j=iwk(i)+1,iwk(i+1)
            ne = iwork(j)
c
            if ( jwk1(ne).ne.i ) then
               jwk1(ne) = i
               if (  nrdpres.eq.0 ) then
                  do l=1,ncppres
                     n1 = indgp + (l-1)*ndgpres + ne
                     if( imp(n1).eq.0 ) npres = npres + 1
                  enddo
               endif
            endif
c
            do k=1,ndl
               n = inod(k,ne)
               if ( jwk(n).ne.i ) then
                  jwk(n) = i
c
                  if ( nrdpres.gt.0 .and. n.le.ndgpres ) then
                     do l=1,ncppres
                        n1 = indgp + (l-1)*ndgpres + n
                        if( imp(n1).eq.0 ) npres = npres + 1
                     enddo
                  endif
c
                  if ( nrdelas.gt.0 .and. n.le.ndgelas ) then
                     do l=1,ncpelas
                        n1 = indgr + (l-1)*ndgelas + n
                        if ( imp(n1).eq.0 ) nelas = nelas + 1
                     enddo
                  endif
c
	            if ( nrdphi.gt.0 .and. n.le.ndgphi )then
                     do l=1,ncpphi
                        if(imp(indgphi+n).eq.0) nphi = nphi + 1
                        if(imp(indgpsi+n).eq.0) npsi = npsi + 1
                     enddo
                  endif
c
                  if ( nrdvelo.gt.0 .and. n.le.ndgvelo ) then
                     do l=1,ncpvelo
                        n1 = indgu + (l-1)*ndgvelo + n
                        if ( imp(n1).eq.0 ) then
                           nvelo = nvelo + 1
                        elseif ( imp(n1).lt.0 ) then
                           nv = ndg + abs(imp(n1))
                           if ( jwk(nv).ne.i ) then
                              jwk(nv) = i
                              nvelo = nvelo + ncpb
                           endif
                        endif
                     enddo
                  endif
c
               endif
            enddo
         enddo
c
         do l=1,ncpvelo
            if ( i.le.ndgvelo ) then
               n1 = indgu + (l-1)*ndgvelo + i
               if ( imp(n1) .eq. 0 ) then
                  neqA = neqA + 1
                  neltA = neltA + nvelo+nelas+nphi+npsi+ntemp
                  if ( matglb ) neltA = neltA + npres
               endif
            endif
         enddo
c
         do l=1,ncpelas
            if ( i.le.ndgelas ) then
               n1 = indgr + (l-1)*ndgelas + i
               if ( imp(n1) .eq. 0 ) then
                  neqA = neqA + 1
                  neltA = neltA + nvelo+nelas+nphi+npsi+ntemp
                  if ( matglb ) neltA = neltA + npres
               endif
            endif
         enddo
c
         do l=1,ncppres
            if ( i.le.ndgpres ) then
               n1 = indgp + (l-1)*ndgpres + i
               if ( imp(n1) .eq. 0 ) then
                  neqB = neqB + 1
                  if ( matglb ) then
                    if ( nrdpres.eq.0 ) then
                       neltA = neltA + ncpelas*(3+(nrdelas-1)*3 )
     &                               + ncpvelo*(3+(nrdvelo-1)*3 )
     &                               + ncptemp*(3+(nrdtemp-1)*3 )
     &                               + ncpphi *(3+(nrdphi -1)*3 )*2+1
                    else
                     neltA = neltA + (nvelo+nelas+npres+nphi+npsi+ntemp)
                    endif
                  else
                    if ( nrdpres.eq.0 ) then
                       neltB = neltB + ncpelas*(3+(nrdelas-1)*3 )
     &                               + ncpvelo*(3+(nrdvelo-1)*3 )
     &                               + ncptemp*(3+(nrdtemp-1)*3 )
     &                               + ncpphi *(3+(nrdphi-1)*3 )*2
                       neltS = neltS + 1
                    else
                       neltB = neltB + (nvelo+nelas+nphi+npsi+ntemp)
                       neltS = neltS + npres
                    endif
                  endif
               endif
            endif
         enddo
c
         do l=1,ncptemp
            if ( i.le.ndgtemp ) then
               n1 = indgtm + (l-1)*ndgtemp + i
               if ( imp(n1) .eq. 0 ) then
                  neqA = neqA + 1
                  neltA = neltA + nvelo+nelas+nphi+npsi
                  if ( matglb ) neltA = neltA + npres
               endif
            endif
         enddo
c
         do l=1,ncpphi
            if(i.le.ndgphi)then
               if(imp(indgphi+i).eq.0)then
                  neqA=neqA+1
                  neltA=neltA + nvelo + nelas + nphi +npsi + ntemp
               endif
               if(imp(indgpsi+i).eq.0)then
                  neqA=neqA+1
                  neltA=neltA + nvelo + nelas + nphi +npsi + ntemp
               endif
               if(matglb) neltA=neltA+npres*2
            endif
         enddo
                
c
      enddo
c
c     for particle equations
c
      neqA = neqA + ncpb*nbrigid
      nwall = nic-nbrigid
      do np=1,nbrigid
        nvelo = 0
        nelas = 0
        npres = 0
        ntemp = 0
        nphi  = 0
        npsi  = 0
        nv = ndg + np
        do nd=ic(nwall+np),ic(nwall+np+1)-1
          i = ibdnod(nd)
          do j=iwk(i)+1,iwk(i+1)
            ne = iwork(j)
c
            if ( jwk1(ne).ne.nv ) then
               jwk1(ne) = nv
               if (  nrdpres.eq.0 ) then
                  do l=1,ncppres
                     n1 = indgp + (l-1)*ndgpres + ne
                     if( imp(n1).eq.0 ) npres = npres + 1
                  enddo
               endif
            else
               cycle
            endif
c
            do k=1,ndl
               n = inod(k,ne)
               if ( jwk(n).ne.nv ) then
                  jwk(n) = nv
c
                  if ( nrdpres.gt.0 .and. n.le.ndgpres ) then
                     do l=1,ncppres
                        n1 = indgp + (l-1)*ndgpres + n
                        if( imp(n1).eq.0 ) npres = npres + 1
                     enddo
                  endif
c
                  if ( nrdelas.gt.0 .and. n.le.ndgelas ) then
                     do l=1,ncpelas
                        n1 = indgr + (l-1)*ndgelas + n
                        if ( imp(n1).eq.0 ) nelas = nelas + 1
                     enddo
                  endif
c
                  if ( nrdvelo.gt.0 .and. n.le.ndgvelo ) then
                     do l=1,ncpvelo
                        n1 = indgu + (l-1)*ndgvelo + n
                        if ( imp(n1).eq.0 ) then
                           nvelo = nvelo + 1
                        elseif ( imp(n1).lt.0 ) then
                           n1 = ndg + abs(imp(n1))
                           if ( jwk(n1).ne.nv ) then
                              jwk(n1) = nv
                              nvelo = nvelo + ncpb
                           endif
                        endif
                     enddo
                  endif
c
                  if(ncpphi.eq.1)then
                  if ( nrdphi.gt.0 .and. n.le.ndgphi ) then
                     if( imp(indgphi+n).eq.0 ) nphi = nphi + 1
                     if( imp(indgpsi+n).eq.0 ) npsi = npsi + 1
                  endif
                  endif
c
               endif
             enddo
           enddo
c
        enddo
        neltA = neltA + ncpb*(nvelo+nelas+nphi+npsi+ntemp)
        if ( matglb ) neltA = neltA + ncpb*npres
c
      enddo
c
      if ( matglb ) then
         neqA = neqA+neqB
         neqB = 0
      endif
c
      return
      end
