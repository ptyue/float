c***********************************************************************
      subroutine forces ( fxnbd,fynbd,fcspri, iou,iax,kaxi, nelt,neq, 
     &                    neqA,neqB,neltA,neltB,matglb,
     &                    upwind,inert,newton,trans,tsc,
     &                    nnode,nelem,inod,ngmx,x,y,ishape,reft, 
     &                    z,zs,umesh,ncpg,
     &                    nbd,ibdnod,nveli,iglo,ntot,nbrigid,ncpb,xpos,
     &                    ivisc,itrel, ro,grav,pvis,pelas,pterm,icoe,
     &                    SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &                    ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                  ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &                  nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &                  indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                    shg2,shg2x,shg2y,we2,
     &                    a,ia,ja,b,ib,jb,imp,iperm,arhs, jwk,st,rh,
     &                    phgam,phlam,pheps,phsft)

c
c     this subroutine computes the forces along the boundaries
c       ( in contact with the continuous phase only)
c
c     July 20, 1994, by Howard H. Hu
c     cleaned-up 5/16/1996
c     Feb 18, 2005, by Pengtao Yue
c        some modifications due to phase-field equation
c***********************************************************************
      implicit none
      integer nbd,nnode,nelem,nelt,neq,iax,kaxi, ntot,
     &        iperm(neq),ishape(nelem),iou,iglo(ntot),nbrigid,
     &        ngmx,inod(ngmx,nelem),reft(nelem),ibdnod(nbd),nveli(2),
     &        ivisc,itrel,neltA,neltB,neqA,neqB,ia(neqA+1),ja(neltA),
     &        ib(neqB+1),jb(neltB),imp(neq),jwk(neq),ncpg,ncpb, icoe
      integer ncpvelo,ncppres,ncptemp,ncpelas,ndggeom,ndgvelo,
     &        ndgelas,nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &        ndgpres,ndgtemp,indgu,indgp,indgtm,indgr,indgnb,
     &        indgphi,indgpsi,ncpphi,ndgphi
      real*8 fxnbd(nbd),fynbd(nbd),z(neq),zs(neq),umesh(ncpg,ndgvelo),
     &       x(ndggeom),y(ndggeom),arhs(neq),xpos(ncpb,nbrigid),tsc,
     &       a(neltA),b(neltB),st(ntot,ntot),rh(ntot)
      double precision ro(icoe),grav(2,icoe),pvis(10,icoe),
     &                 pterm(10,icoe),pelas(10,icoe)
      real*8 phgam,phlam,pheps,phsft
c
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      logical SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,fcspri,
     &        upwind,inert,newton,trans,matglb
c
      integer nbd0,i,m,ishm,nvel,l,ndlgeom,k,nd,nk,ik,imat, ndloc
      logical jacoct
c
c     1. initialization
c     -----------------
      nbd0 = nbd/2
      do i=1,neq
         arhs(i) = 0.d0
      enddo
      jacoct = .true.
c
c     2. loop over element
c     --------------------
      do 20 m=1,nelem
c
c     if(mod(m,100).eq.0) print*, 'For m=',m,' reft=',reft(m)
c     if(reft(m).ne.1) print*, 'ref for m=',m,' is',reft(m)
c     goto 20
c
        ishm = ishape(m)
        nvel = nveli(ishm+1)
c
        ndlgeom = ndloc(ishm, nrdgeom)
        do l=1,ndlgeom
          k = inod(l,m)
          if ( k.le.nbd0 ) goto 25
        enddo
        goto 20
c
c       3. calculation of forces on boundary nodes
c       ------------------------------------------
25      imat = reft(m)
        if ( imat.ne.1 ) then
           print *, 'FORCES: bounary section not in contact with',
     &              ' the continuous phase'
           stop
        endif
c
c     print*, 'Matnew for element=',m,' at (x,y)='
c     write(*,*) (x(inod(i,m)),y(inod(i,m)),i=1,3)
c
        call matnew(m,ishm,iax,kaxi,nelem,inod,ngmx,x,y,
     &              nbrigid,ncpb,xpos,umesh,ncpg,z,zs,
     &              ivisc,itrel, ro(imat),grav(1,imat),
     &              pvis(1,imat),pelas(1,imat),pterm(1,imat),
     &              trans,tsc,jacoct,inert,newton,upwind,
     &              SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &              ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &              nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &              ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &              indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &              matglb,neq,neqA,neqB,nvel,neltA,neltB,iglo,ntot, 
     &              a,ja,ia,b,jb,ib,arhs,imp,iperm,jwk,rh,st,
     &              shg2,shg2x,shg2y,we2,
     &              phgam,phlam,pheps,phsft)
20    continue
c
      do i=1,nbd
         nd = ibdnod(i)
         nk = indgu+nd
         ik = iperm(nk)
         fxnbd(i) =  -arhs(ik)
         nk = indgu+ndgvelo+nd
         ik = iperm(nk)
         fynbd(i) =  -arhs(ik)
      enddo
c
c     4. printing of results
c        -------------------
      if ( fcspri ) then
c
        write(iou,*) '**************************'
        write(iou,*) '* forces   postprocessor *'
        write(iou,*) '**************************'
        write(iou,1200)
1200    format(/72('-')/25x,'forces on boundary nodes'/6x,
     &         'node',22x,'fx or fr',5x,'fy or fz')
        do i=1,nbd
          write(iou,'(21x,i3,2(3x,1pe14.7))') i,fxnbd(i),fynbd(i)
        enddo
c
      endif
c
      return
      end
