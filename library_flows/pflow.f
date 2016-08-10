c***********************************************************************
      subroutine pflow (iou, tsc,tsc1,dtin,endtim,flowtype,
     &                  neq,nvar, mshslpw,mshslpp, nnew1,nnew2,
     &                  nvert,nnode,nelem,inod,x,y,reft,
     &                  nbd,ibdnod,nic,ic,nbound,nside,
     &                   nbd0,nic0,nbound0,
     &                  maxnod,maxelt,maxnbd,maxbdp,maxcoe,maxvar,maxnz,
     &                  nbrigid,nbfluid,nbelast,
     &                   xpos,fxym,amas,fpart,ifixpt,
     &                  ioptu,imu,lfilu,maxitu,epsu,lfilter,
c     &                 ioptphi,imphi,lfilphi,maxitphi,epsphi,lfilterphi, 
     &                  iopts,lfils,maxits,epss,
     &                  ioptm,lfilm,maxitm,epsm,
     &                  ioptp,maxitp,epsp,
     &                  nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,
     &                   nrdstrm,nrdmsh,nrdphi,
     &                  ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,
     &                   ndgstrm,ndgmsh,ndgphi,
     &                  SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &                  ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                  indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                   indgstrm,
     &                  ngmx,ncpu,ncpg,ncps,ncpb,
     &                  fxnbd,fynbd, umesh,z,zs, stress,dt,
     &                  ibdu,ibdt,ibds,ibdp,ibdn, 
     &                  densu,denst,denss,densp,densn,
     &                  noexec,inert,ivheat,newton,upwind,
     &                  bdpr,resold,rstart,trans,
     &                  icase,iax,ifluid,ivisc,itemp,itrel,icoe,
     &                  iforce,istres,ipnode,isnode,itmax,iorder,
     &                  ro,grav,pvis,pterm,pelas,zconv,
     &                  gamma,plambda,epsilon,shift,
     &                  wallrela,wallener,
     &                  iwork,maxiwk,rwork,maxrwk,
     &                  refr,refu,refp,dropone)
c
c     this routine calls the main subroutines
c
c     modified version for particle simulation - april, 1991
c     July 20, 1994, cleaned up for viscoelastic models, Howard Hu
c     Modified  by Pengtao Yue
c        Feb 18,2005,   add the Cahn-Hillard equation
c        Feb 25,2005,   changed ljwk and liwork
c     Notes for amphi:
c        nbd0 and nbd are the same.
c***********************************************************************
      implicit none
      integer ngmx,ncpu,ncpg,ncps,ncpb,flowtype
      integer iou,nvert,nnode,nelem,inod(ngmx,nelem),reft(nelem),
     &        nbd,ibdnod(nbd),nic,ic(nic+1),nbound,nside(nbound),
     &        nbd0,nic0,nbound0,
     &        nbrigid,nbfluid,nbelast,nvar,neq,nnew1,nnew2, 
     &        ioptu,imu,lfilu,maxitu,
     &        iopts,lfils,maxits, ioptm,lfilm,maxitm,
     &        ioptp,maxitp,
     &       ibdu(ncpu,nbd0),ibdt(nbd0),ibds(nbd0),ibdp(nbd0),ibdn(nbd0)
      integer nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdstrm,nrdmsh
      integer ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgstrm,ndgmsh
      integer indgu,indgp,indgtm,indgr,indgnb,indgstrm
      integer maxnod,maxelt,maxnbd,maxbdp,maxcoe,maxvar,maxnz
      integer ncpvelo,ncppres,ncptemp,ncpelas
      logical SLVvelo,SLVpres,SLVtemp,SLVelas
      double precision x(nnode),y(nnode),fxnbd(nbd0),fynbd(nbd0),
     &                 xpos(ncpb,nbrigid),fxym(ncpb,nbrigid),
     &                 amas(ncpb,nbrigid),fpart(ncpb,nbrigid),
     &                 tsc,tsc1,umesh(ncpg,nnode),
     &                 z(maxnz),zs(maxnz),dtin,endtim,epsu,epss,epsm,
     &                 epsp,densu(ncpu,nbd0),denst(nbd0),
     &                 denss(ncps,nbd0),densp(nbd0),densn(nbd0),
     &                 stress(ncps,ndgelas),dt
      logical mshslpw,mshslpp,ifixpt(ncpb)
      integer maxmt1,maxmt2,maxiwk,maxrwk
      integer iwork(maxiwk)
      double precision rwork(maxrwk)
      logical noexec,inert,ivheat,newton,upwind
      logical bdpr,resold,rstart,trans, lfilter,matglb
      integer icase,iax,ifluid,icoe,ivisc(icoe),itemp,itrel
      integer iforce,istres,ipnode,isnode,itmax,iorder
      double precision ro(icoe),grav(2,icoe),pvis(10,icoe),
     &                 pterm(10,icoe),pelas(10,icoe),zconv
      real(8) refr,refu,refp
c*CH start
      integer ncpphi,nrdphi,ndgphi
      integer indgphi, indgpsi
      logical SLVphi
	real*8 gamma,epsilon,plambda,shift
c     moving contact-line
      real(8) wallrela(nic),wallener(nic)
c
      logical dropone
c      integer ioptphi,imphi,lfilphi,maxitphi
c      real*8 epsphi
c      logical lfilterphi
c*CH end
c
c     local variables
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
     &      ,shg1(3,3,2),shg1x(3,3,2),we1(3,2), smax,dsmax
      logical jacoct,oldmsh,fcspri
      integer kaxi,iter,nveli(2), ntot,ndloc
c
      integer itte,ittg, i,j,k,k1,n, neqA,neqB,neqt,neltA,neltB,neltS
      integer nelt,nelt2
      integer lshape,liglo,limp,liperm,ljperm,
     &        liwk,ljwk,lia,lja,limp2,liperm2,ljperm2,lia2,lja2,
     &        lb,lib,ljb,
     &        liwork2,liwork,leniwk, nn,mxiwk2,mxrwk2,mxiwk1,mxrwk1,
     &        lubdmsh,lvbdmsh,larhs,lvec,lvec2,la,la2,lrwork2,
     &        lrwork,lenrwk, lrh,lst,ierr
     &       ,lsx,lisx,ljsx, mxrwk3,lengS,lengw
      logical conv,linesearch
      integer lrscli,lrsclj
c      real(8),save:: refr=1.,refu=1.,refp=1.
c
      logical, save :: init=.true.
      integer njacot
c     

C
C     Inexact-Newton variables
C     res            -- linear system residual
C     eta            -- forcing term (tolerance)
C     oldRhs, newRhs -- norms of function evaluations
C     globits        -- global iteration count for nonlinear solve
      double precision res, eta, oldRhs, newRhs,intRhs, tauline
      integer globits
c
      if(init)then
         njacot=1
         init=.false.
      else
         njacot=nnew1
      endif
c
      intRhs = 1.d0
      oldRhs = 1.d20
      newRhs = 0.d0
      globits = 1
c
      ierr = 0
      fcspri = .true.
      oldmsh = .false.
      jacoct = .false.
      conv   = .false.
c
      nvar = ncpelas*ndgelas + ncpvelo*ndgvelo +
     &       ncppres*ndgpres + ncptemp*ndgtemp + ncpphi*ndgphi*2
      neq = nvar + ncpb*nbrigid
      ntot = ncpvelo*ndloc(0,nrdvelo) + ncppres*ndloc(0,nrdpres) +
     &       ncptemp*ndloc(0,nrdtemp) + ncpelas*ndloc(0,nrdelas) +
     &       ncpphi *ndloc(0,nrdphi )*2
c
      matglb = .true.
      if ( ioptu/100 .ne. 0 ) matglb = .false.
c
c     check data
c     ----------
      lshape = 1
      call chkdat (iou, icase,iax,kaxi, neq,nvar,
     &             resold,rstart,upwind,newton,inert,
     &             trans,dtin,endtim,bdpr,iforce,itmax,zconv,
     &             nvert,nnode,nelem,inod,ngmx,x,y,iwork(lshape),
     &             nbd0,ibdnod,nic0,ic, ipnode,isnode,
     &             ibdu,ibdt,ibds,densu,denst,denss,ncpu,ncps,
     &             ifluid,ivisc(1),itemp,itrel,ivheat,istres,icoe,
     &             ro,grav,pvis,pterm,pelas, iorder, ndgstrm,
     &             maxnod,maxelt,maxnbd,maxbdp,maxcoe,maxvar,maxnz)
c
c     preprocessors
c     -------------
      call shap2d ( shg2,shg2x,shg2y,we2 )
      call shap1d ( shg1,shg1x,we1 )
c
      call beginc ( icase,iax,kaxi,neq,iou, x,y, z,zs,
     &              trans,tsc,tsc1, shg1,shg1x,we1,
     &              SLVvelo,SLVtemp,SLVelas,SLVpres,
     &              ndggeom,ndgvelo,ndgtemp,ndgelas,ndgpres,
     &              nrdgeom,indgu,indgp,indgtm,indgr,
     &              nbd0,ibdnod,nic0,ic,nbound0,nside, ipnode,
     &              ibdu,ibdt,ibds,ibdp,densu,denst,denss,densp,
     &              ncpu,ncps )
c
      liglo   = lshape  + nelem
      limp    = liglo   + ntot
      liwk    = limp    + neq
      ljwk    = liwk    + max(nnode,neq)+1
      liwork  = ljwk    + nnode+nbrigid+nelem
c
      call matlen1 ( nelt2, nelem,inod,ngmx,iwork(lshape), nrdmsh,
     &               iwork(liwk),iwork(ljwk),iwork(liwork) )
c
      call matlen2 ( neltA,neltB,neltS,neq,neqA,neqB,nelem,inod,
     &              iwork(lshape),
     &              nbd0,ibdnod,nic0,ic, ibdu,ibdt,ibds,ibdp,ipnode,
     &              nbrigid, iwork(limp),iwork(liwork),iwork(liwk),
     &              iwork(ljwk),iwork(ljwk+ndggeom+nbrigid),
     &              ngmx,ncpu,ncpb, matglb,
     &              indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &              ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &              nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &              ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi)
c
      dsmax = 0.d0
      smax  = 1.d0
c      smax  = dnrm2(neqA+neqB,z,1)      
c
      if ( ioptu/100.ge.2 .and. ioptp/10.eq.3 ) neltS = neltB

      neqt = neqA+neqB
      write(iou,*) 'number of equations,       neqA  = ',neqA
      write(iou,*) 'number of equations,       neqB  = ',neqB
      write(iou,'(/A,i9)') 'Length of fluid-particle matrix A =',neltA
      write(iou,'(A,i9)')  'Length of fluid-particle matrix B =',neltB
      write(iou,'(A,i9)')  'Length of app. Schur complement S =',neltS
      write(iou,'(A,i9/)') 'Length of matrix for mesh         =',nelt2
c
c     storage management - integer
c     ----------------------------
      nn = max(ndgstrm,ndgmsh)
      if ( mod(ioptm,10).eq.1 ) then
         maxmt2 = nn
         mxiwk2 = maxmt2
      elseif ( mod(ioptm,10).eq.2 ) then
         maxmt2 = nelt2  + 1
         mxiwk2 = maxmt2 + 2*nn
      elseif ( mod(ioptm,10).eq.3 ) then
         maxmt2 = nelt2  + 2*lfilm*nn + 1
         mxiwk2 = maxmt2 + 3*nn
      endif
c
      k = mod(ioptu,10)
      if ( k.eq.1 ) then
         maxmt1 = neqA
         lengw  = 0
         mxiwk1 = maxmt1
      elseif ( k.eq.2) then
         maxmt1 = neltA+1
         lengw  = neqA
         mxiwk1 = maxmt1 + neqA
      elseif ( k.eq.3 ) then
         maxmt1 = neltA+1 + 2*lfilu*neqA
         lengw  = 2*neqA
         mxiwk1 = maxmt1 + neqA
      endif

      if ( ioptu/100.ge.2 ) then
         if (ioptp/10.eq.1 ) then
            lengS = neqB
            mxiwk1 = mxiwk1 + lengS
         elseif ( ioptp/10.ge.2 ) then
            lengS = neltS+1
            lengw = max(lengw,neqB)
            mxiwk1 = mxiwk1 + lengS + neqB
         endif
      endif
      mxiwk1 = mxiwk1 + lengw
      mxiwk1 = max(mxiwk1, ntot*nelem)
c
      liperm  = ljwk    + neq
      ljperm  = liperm  + neq
      limp2   = ljperm  + neq
      liperm2 = limp2   + nn
      ljperm2 = liperm2 + nn
      lia2    = ljperm2 + nn
      lia     = lia2    + nn+1
      lib     = lia     + neqA+1
      lisx    = lib     + neqB+1
      lja2    = lisx    + neqB+1
      liwork2 = lja2    + nelt2
      lja     = liwork2 + mxiwk2
      ljb     = lja     + neltA
      ljsx    = ljb     + neltB
      liwork  = ljsx    + neltS
      leniwk  = liwork  + mxiwk1

      write(iou,'(A,i10,A,i10)') 'iwork size: assigned=',maxiwk,
     &                           ', used=',leniwk
      if ( leniwk.gt.maxiwk ) then
         write(*,*) 'iwork too short, should be larger than',leniwk
         write(*,'(80("-"))')
         write(*,*) 'maxnod',maxnod,'maxelt',maxelt,'maxvar',maxvar
         write(*,*)'nnode',nnode,'nelem',nelem,'neq',neq
         write(*,*)'leniwk/neq',real(leniwk)/neq
         write(*,'(80("-"))')
         ierr = 1
      endif
c
c     storage management - real
c     -------------------------
      mxrwk2 = maxmt2 + 3*nn

      if ( mod(ioptu/10,10).eq.1 ) then
         mxrwk1 = maxmt1 + (imu+2)*neqt
      elseif ( mod(ioptu/10,10).eq.2 ) then
         mxrwk1 = maxmt1 + 7*neqt
      endif

      mxrwk3 = neqA
      if ( (ioptu/100).ge.2 ) then
        k = mod(ioptp,10)
        if ( k.eq.2 .or. k.eq.3 ) then
           mxrwk3  = mxrwk3 + 3*neqB + neqA
        elseif ( k.eq.5 .or. k.eq.6 ) then
           mxrwk3  = mxrwk3 + 5*neqB + neqA
        endif
        mxrwk1 = mxrwk1 + lengS
      endif
      mxrwk1 = mxrwk1 + mxrwk3
c
      lubdmsh = 1
      lvbdmsh = lubdmsh + nbd
      lrh     = lvbdmsh + nbd
      lst     = lrh     + ntot
      larhs   = lst     + ntot*ntot
      lvec    = larhs   + neqt
      lvec2   = lvec    + neqt
      la      = lvec2   + nn
      lb      = la      + neltA
      la2     = lb      + neltB
      lsx     = la2     + nelt2
      lrwork2 = lsx     + neltS
      lrwork  = lrwork2 + mxrwk2
      lenrwk  = lrwork  + mxrwk1
c     array for rescaling
      lrscli  = lenrwk + 1 
      lrsclj  = lrscli + neq
      lenrwk  = lrsclj + neq     
c
      write(iou,'(A,i10,A,i10/)') 'rwork size: assigned=',maxrwk,
     &                           ', used=',lenrwk
      if ( lenrwk.gt.maxrwk ) then
         write(*,*) 'rwork too short, should be larger than',lenrwk
         write(*,'(80("-"))')
         write(*,*) 'maxnod',maxnod,'maxelt',maxelt,'maxvar',maxvar
         write(*,*)'nnode',nnode,'nelem',nelem,'neq',neq
         write(*,*)'lenrwk/neq',real(lenrwk)/neq
         write(*,'(80("-"))')
         ierr = 1
      endif
      if ( ierr.ne.0 ) stop
c
      if ( noexec ) then
         write(iou,*) 'noexec = true: stop after testing data'
         stop
      endif
c
      if ( (resold.or.rstart) .and. trans ) goto 41
c
c     steady loop
c     -----------
      resold = .true.
      rstart = .true.
c
c     calculation of velocities and temperature
c
      do iter=1,itmax
        call front(iax,kaxi, iou, neq,z,zs,umesh,ncpg, smax,dsmax,
     &             nvert,nnode,nelem,inod,ngmx,x,y,iwork(lshape),
     &             reft,nbd0,ibdnod,nic0,ic, ipnode, neqA,neqB,
     &             nbound,nside,
     &             nbrigid,nbfluid,nbelast,
     &             ncpb,xpos,fxym,amas,fpart,ifixpt,
     &             nveli,iwork(liglo),ntot,
     &             ivisc,itrel,ro,grav,pvis,pelas,pterm,icoe,
     &             .false.,tsc, inert,newton,jacoct,upwind,
     &             ibdu,ibdt,ibds,ibdp,ncpu,
     &             SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &             ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &             ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi,
     &             nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &             indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi, 
     &             shg2,shg2x,shg2y,we2,shg1,shg1x,we1,
     &             ioptu,imu,lfilu,maxitu,epsu,lfilter,
     &             rwork(la),iwork(lia),iwork(lja),iwork(limp),
     &             iwork(liperm),iwork(ljperm),rwork(larhs),
     &             rwork(lb),iwork(lib),iwork(ljb),neltA,neltB,matglb,
     &             rwork(lsx),iwork(lisx),iwork(ljsx),neltS,
     &             rwork(lvec),iwork(liwk),iwork(ljwk),
     &             iwork(liwork),mxiwk1, rwork(lrwork),mxrwk1,
     &             rwork(lrh),rwork(lst),
     &             res, eta, oldRhs,newRhs,intRhs, globits, conv
     &            ,linesearch,tauline, stress,ncps,dt,
     &             ioptp,maxitp,epsp,smax, 
     &             gamma,plambda,epsilon,shift,
     &             wallrela,wallener,ibdn,densn,
     &             refr,refu,refp,rwork(lrscli),rwork(lrsclj),
     &             flowtype,dropone)

c
        if ( conv. or. dabs(dsmax).le.zconv ) goto 32
      enddo
32    if ( .not.trans ) goto 51
c
c     transient loop
c     --------------
41    itte = 0
      ittg = 0
 45   continue
c      tauline = 1.d0
c
c 46   continue
c
c     mesh velocity
c     -------------
      if ( nbrigid+nbfluid+nbelast.gt.0 ) then
         call bdmvel ( rwork(lubdmsh),rwork(lvbdmsh),mshslpp,mshslpw,
     &                 nbd,ncpb,nbrigid,nbfluid,nbelast,
     &                 nic,ic,nnode,x,y,ibdnod,
     &                 z(indgu+1),z(indgu+ndgvelo+1),z(indgnb+1),xpos )
c
         if(flowtype/100.eq.0)
     &   call mshmov (oldmsh,nrdmsh,nrdgeom,ndgmsh,ndggeom,ndgvelo,
     &                iou, umesh,ncpg,rwork(lubdmsh),rwork(lvbdmsh),
     &                nelem,inod,ngmx,x,y,iwork(lshape),
     &                nic,ic,nbd,ibdnod,
     &                shg2,shg2x,shg2y,we2,
     &                ioptm,lfilm,maxitm,epsm,smax,
     &                rwork(la2),iwork(lja2),nelt2,iwork(lia2),
     &                iwork(liperm2),iwork(ljperm2),iwork(limp2),
     &                rwork(larhs),rwork(lvec2),
     &                iwork(liwk),iwork(ljwk),iwork(liwork2),mxiwk2,
     &                rwork(lrwork2),mxrwk2 )
c
      endif
c
c     solving the flow field
c     ----------------------
c      if(mod(itte,5).eq.0)jacoct=.false.
        call front(iax,kaxi, iou, neq,z,zs,umesh,ncpg, smax,dsmax,
     &             nvert,nnode,nelem,inod,ngmx,x,y,iwork(lshape),
     &             reft,nbd0,ibdnod,nic0,ic, ipnode, neqA,neqB,
     &             nbound,nside,
     &             nbrigid,nbfluid,nbelast,
     &             ncpb,xpos,fxym,amas,fpart,ifixpt,
     &             nveli,iwork(liglo),ntot,
     &             ivisc,itrel,ro,grav,pvis,pelas,pterm,icoe,
     &             trans,tsc, inert,newton,jacoct,upwind,
     &             ibdu,ibdt,ibds,ibdp,ncpu,
     &             SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &             ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &             ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi,
     &             nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &             indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi, 
     &             shg2,shg2x,shg2y,we2,shg1,shg1x,we1,
     &             ioptu,imu,lfilu,maxitu,epsu,lfilter,
     &             rwork(la),iwork(lia),iwork(lja),iwork(limp),
     &             iwork(liperm),iwork(ljperm),rwork(larhs),
     &             rwork(lb),iwork(lib),iwork(ljb),neltA,neltB,matglb,
     &             rwork(lsx),iwork(lisx),iwork(ljsx),neltS,
     &             rwork(lvec),iwork(liwk),iwork(ljwk),
     &             iwork(liwork),mxiwk1, rwork(lrwork),mxrwk1,
     &             rwork(lrh),rwork(lst),
     &             res, eta, oldRhs,newRhs,intRhs, globits, conv
     &            ,linesearch,tauline, stress,ncps,dt,
     &             ioptp,maxitp,epsp,smax, 
     &             gamma,plambda,epsilon,shift,
     &             wallrela,wallener,ibdn,densn,
     &             refr,refu,refp,rwork(lrscli),rwork(lrsclj),
     &             flowtype,dropone)
c
      if ( itte.eq.0 ) intRhs = newRhs
c
c      if ( linesearch ) goto 46
      if ( conv ) goto 51
c
      itte = itte+1
      if (mod(itte,njacot).eq.0 ) then
         jacoct = .false.
      else
         jacoct = .true.
      endif
c
c      write(*,'(1pe13.4)') dsmax
c      if (ittg.lt.1 .and. itte.ge.nnew1 .and. abs(dsmax).gt.zconv) then
      if ( itte.ge.nnew2 ) then
         write(*,*)'itte: ',itte,' >= nnew2: ',nnew2
c         stop
         goto 51
      endif
      if ( dabs(dsmax).lt.zconv ) then
         conv = .true.
         goto 51
      endif
      go to 45
c
51    continue
c
c     computation of stream function
c     ------------------------------
      if ( bdpr ) 
     &  call masbal ( iou, kaxi, nrdgeom,ndggeom,x,y, 
     &                ndgvelo,z(indgu+1),z(indgu+ndgvelo+1),
     &                nbd0,ibdnod,nic0,ic,nbound0,nside, shg1,shg1x,we1)
c
      if ( icase.eq.1 .or. icase.eq.2 .or. icase.eq.4 )
     &  call streamf( z(neq+1), iou,kaxi,
     &                nrdstrm,nrdgeom,nrdvelo, ndgstrm,ndggeom,ndgvelo,
     &                nelem,inod,ngmx,x,y,iwork(lshape),
     &                 z(indgu+1),z(indgu+ndgvelo+1),
     &                nbd0,ibdnod,nic0,ic,nbound0,nside, isnode,
     &                shg2,shg2x,shg2y,we2, shg1,shg1x,we1,
     &                iopts,lfils,maxits,epss,
     &                rwork(la2),iwork(lja2),maxmt2,iwork(lia2),
     &                 iwork(liperm2),iwork(ljperm2),rwork(larhs),
     &                 iwork(liwk),iwork(ljwk),
     &                 iwork(liwork2),mxiwk2,rwork(lrwork2),mxrwk2 )
c
c     postprocesor
c     ------------
      if ( iforce.eq.1 ) 
     &  call  forces ( fxnbd,fynbd,fcspri, iou,iax,kaxi, nelt,neq, 
     &                 neqA,neqB,neltA,neltB,matglb,
     &                 upwind,inert,newton,trans,tsc,
     &                 nnode,nelem,inod,ngmx,x,y,iwork(lshape),reft,
     &                 z,zs,umesh,ncpg,
     &                 nbd,ibdnod, nveli,iwork(liglo), ntot,
     &                 nbrigid,ncpb,xpos,
     &                 ivisc(1),itrel, ro,grav, pvis,pelas,pterm,icoe,
     &                 SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &                 ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                 ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &                 nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &                 indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                 shg2,shg2x,shg2y,we2,
     &                 rwork(la),iwork(lia),iwork(lja),rwork(lb),
     &                 rwork(lib),rwork(ljb),iwork(limp),iwork(liperm),
     &                 rwork(larhs), iwork(ljwk),rwork(lst),rwork(lrh),
     &                 gamma,plambda,epsilon,shift)
c
c     reset mesh velocity inside elastic particles
c     --------------------------------------------
      if ( nbelast.gt.0 ) then
         k1 = 1+nbrigid+nbfluid
         do i=1,nelem
            k = reft(i)
            if ( k.gt.k1 ) then
               do j=1,6
                  n = inod(j,i)
                  umesh(1,n) = z(indgu+n)
                  umesh(2,n) = z(indgu+ndgvelo+n)
               enddo
            endif
         enddo
      endif
c
      return
      end
