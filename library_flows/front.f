c***********************************************************************
      subroutine front (iax,kaxi, iout, neq,z,zs,umesh,ncpg,smax,dsmax,
     &                  nvert,nnode,nelem,inod,ngmx,x,y,ishape,reft,
     &                  nbd,ibdnod,nic,ic, ipnode,neqA,neqB,
     &                  nbound,nside,
     &                  nbrigid,nbfluid,nbelast,
     &                  ncpb,xpos,fxym,amas,fpart,ifixpt, 
     &                  nveli,iglo,ntot,
     &                  ivisc,itrel,ro,grav,pvis,pelas,pterm,icoe,
     &                  trans,tsc, inert,newton,jacot,upwind,
     &                  ibdu,ibdt,ibds,ibdp,ncpu,
     &                  SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &                  ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                  ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi,
     &                  nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &                  indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                  shg2,shg2x,shg2y,we2,shg1,shg1x,we1,
     &                  ioptu,im,lfil,maxitu,epsu,lfilter,
     &                  a,ia,ja,imp,iperm,jperm,arhs,
     &                  b,ib,jb,neltA,neltB,matglb, sx,isx,jsx,neltS,
     &                  vec,iwk,jwk,iwork,maxiwk,rwork,maxrwk,rh,st,
     &                  res,eta,oldRhs,newRhs,intRhs, globits, conv
     &                 ,linesearch,tauline, stress,ncps,dt,
     &                  ioptp,maxitp,epsp,smaxp,
     &                  phgam,phlam,pheps,phsft,
     &                  wallrela,wallener,ibdn,densn,
     &                  refr,refu,refp,scali,scalj,
     &                  flowtype,dropone)
c
c     input:
c        jacot:     T. use the old matrix(constant-slope Newont's 
c                       method)
c                    F. calculate a new matrix
c        refr,refu,refp:
c                    Characteristic valcues for extra stress, velocity
c                    and pressure     
c        working:
c           arhs(destroyed at output)
c     this routine computes the flow field 
c     Modified by Pengtao Yue, Feb 17,2005
c        add the Cahn-Hillard eqn
c     Modified by Pengtao Yue, Apr 08, 2005
c        matamphi is added to account for different viscosities and 
c        densities in the interfacial flow.
c     Modified by Pengtao Yue, May 12, 2005
c        An Inexact Newton Backtracking method is used. To enhance
c        robustness furthermore, the system is preprocessed by a left
c        scaling matrix before Krylov solvers(this scaling can only be
c        used for global matrix).
c        [1] R.P. Pawlowski, J.N. Shadid, J.P. Simonis, 
c            and H.F. Walker, SAND 2004-1777
c        [2] S.C. Eisenstat and H.F. Walker, SISC, 17(1):16-32,1996
c
c     Notes for amphi:
c        neqA: total number of unknowns (neq in some subs)
c        neq:  total number of variables, including the boudnary ones
c              (nvar in some subs)
c
c***********************************************************************
      implicit none
      integer ncpvelo,ncppres,ncptemp,ncpelas,
     &        ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi,
     &        nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &        indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &        ncpg,ncpb,ncpu,ncps,ncpphi
      integer nnode,nvert,nbd,nelem,neq,neqA,neqB,iax,kaxi,
     &        neltA,neltB,ia(neqA+1),ja(neltA),imp(neq),
     &        iperm(neq),jperm(neq),ib(neqB+1),jb(neltB),
     &        neltS,isx(neqB+1),jsx(neltS),
     &        ngmx,inod(ngmx,nelem),reft(nelem),ibdnod(nbd),
     &        nbound, nside(nbound),nbrigid,nbfluid,nbelast,
     &        nic,ic(nic+1),ishape(nelem),iout,ntot,iglo(ntot),
     &        nveli(2),ioptu,im,lfil,maxitu,ipnode,
     &        ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),ibdp(nbd),
     &        iwk(neq),jwk(neq),maxiwk,maxrwk,iwork(maxiwk),
     &        icoe,ivisc(icoe),itrel
      real*8 z(neq),zs(neq),umesh(ncpg,ndgvelo),a(neltA),
     &       x(ndggeom),y(ndggeom),b(neltB),arhs(neq),vec(neq),
     &       rwork(maxrwk),xpos(ncpb,nbrigid),fxym(ncpb,nbrigid),
     &       amas(ncpb,nbrigid),fpart(ncpb,nbrigid),stress(ncps,ndgelas)
     &      ,dt,tsc,smax,dsmax,epsu, rh(ntot),st(ntot,ntot), sx(neltS)
      double precision ro(icoe),grav(2,icoe),pvis(10,icoe),
     &                 pterm(10,icoe),pelas(10,icoe)
      logical SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &        upwind,inert,newton,trans,jacot, ifixpt(ncpb), matglb
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      real(8) shg1(3,3,2),shg1x(3,3,2),we1(3,2)
      real*8 res, eta, oldRhs, newRhs,intRhs, tauline
      logical conv, linesearch, lfilter
      real(8) refr,refu,refp,scali(neq),scalj(neq)
c     moving contact line
      real(8) wallrela(nic),wallener(nic)
c     stress bc
      integer ibdn(nbd)
      real(8) densn(nbd)
c
      integer flowtype
      logical dropone
c
      integer ioptp,maxitp
      real*8 epsp,smaxp
c
      real*8 phgam,phlam,pheps,phsft
c     local variables
      integer i,neltv,m,ishm,nvel,ierr,j, imat
      integer neqt
      integer nwall,nv,nd,np
      real*8 xp,yp,dx,dy,xcoor,tmp
      logical jacoct,jacoct1
      integer searchiter

c
C     variables for Inexact-Newton Backtracking (INB)
C     res            -- linear system residual
C     eta            -- forcing term (tolerance)
C     eta_max        -- 0.9 (largest allowed eta)
C     eta_0          -- 0.01 (initial eta)
c     theta_min      -- 0.1 (minimum damping factor)
c     theta_max      -- 0.5 (maximum ....)     
c     t              -- 10^-4
c     p0             -- p(0)
c     p1             -- p(1)
c     dp             -- p'(0)
c     dtheta         -- the step size for the FD approximation of p'(0)
C     safegrd        -- safeguarded value of eta
C     eps            -- linear system absolute tolerance
C     oldRhs, newRhs -- norms of function evaluations
C     globits        -- global iteration count for nonlinear solve
c     backtrack      -- 1, use p2 to get the optimal theta
c                       2, use subdivision by 0.5 to find theta
c                       else, no backtracking used
      real*8 safegrd, eps, phi, dnrm2, theta, p0, p1, dp
      real(8),parameter :: eta_0=0.01, eta_max=0.2, eta_min=1.d-4
      real(8),parameter :: theta_min=0.1,theta_max=0.5,t=1.d-4,
     &                     dtheta=1.d-4
      integer,parameter :: backtrack = 1 
      integer globits
      phi = 0.5d0*(1.d0 + sqrt(5.0d0))
      linesearch = .false.
      tauline = 1.d0
      searchiter = 0
      jacoct1 = jacot
10    continue
      jacoct=jacoct1
C
c     initialization
c     --------------
      neqt = neqA + neqB
      do i=1,neqt
         arhs(i) = 0.d0
      enddo
c
      if ( .not.jacoct ) then
         call prefro ( neq,neqA,neqB,iout, nelem,inod,ngmx,ishape,
     &                nbd,ibdnod,nic,ic,ibdu,ibdt,ibds,ibdp,ncpu,ipnode,
     &                nbrigid,ncpb, nveli,
     &                ia,ja,ib,jb,neltA,neltB,isx,jsx,neltS,ioptp,
     &                imp,iperm,jperm,matglb,
     &                iwork,iwk,jwk, iglo,ntot,
     &                indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &                ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi)
         do i=1,neltA
            a(i) = 0.d0
         enddo
         do i=1,neltB
            b(i) = 0.d0
         enddo
         do i=1,neq
            jwk(i) = 0
         enddo
      endif
c
c     construction of the matrix and vector
c     -------------------------------------
      do m=1,nelem
        ishm = ishape(m)
        nvel = nveli(ishm+1)
c
        imat = reft(m)
        if ( imat.le.1+nbrigid+nbfluid ) then
c
c       fluid
c       -----
         if(SLVphi.and.(icoe.eq.2))then
c     two phase systems with different densities and viscosities.
         imat=1
        call matamphi(m,ishm,iax,kaxi,nelem,inod,ngmx,x,y,
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
     &              phgam,phlam,pheps,phsft,
     &              flowtype,dropone)
         else
        imat = min(2,imat)
        call matnew(m,ishm,iax,kaxi,nelem,inod,ngmx,x,y,
     &              nbrigid,ncpb,xpos,umesh,ncpg,z,zs,
     &              ivisc(1),itrel, ro(imat),grav(1,imat),
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
         endif
c
        elseif ( imat.le.1+nbrigid+nbfluid+nbelast ) then
c
c       elastic solid
c       -------------
c       change for A&B matrices
        imat = icoe
        call matsld(m,ishm,nelem,inod,ngmx,x,y,ncpg,z,zs,
     &              ro(imat),grav(1,imat),pelas(1,imat), 
     &              trans,tsc,jacoct,
     &              SLVvelo,SLVpres, 
     &              ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &              nrdgeom,nrdvelo,nrdpres,nrdelas,nrdtemp,nrdphi,
     &              ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &              indgr,indgu,indgp,indgtm,indgphi,indgpsi,matglb,
     &              neq,neqA,neqB,nvel,neltA,neltB,iglo,ntot,
     &              a,ja,ia,b,jb,ib,arhs,imp,iperm,jwk,rh,st,
     &              shg2,shg2x,shg2y,we2, stress,ncps,dt )
        endif
c
      enddo

c
      call partbc ( indgnb,neq,neltA,nbrigid,ncpb, z,zs,tsc,
     &              fxym,fpart,amas, a,ia,arhs,jacoct,ifixpt )
c
c       moving contact settings
c       -----------------------
      if(flowtype/10==12)
     &   call contact(kaxi,nic,nbd,ic,nbound,ibdnod,nside,
     &                nrdgeom,nrdvelo,nrdphi,ndggeom,ndgvelo,ndgphi,
     &                indgu,indgphi,x,y,
     &                neqA,neq,neltA,z,zs,a,ia,ja,iperm,arhs,
     &                wallrela,wallener,tsc,phlam,pheps,
     &                shg1,shg1x,we1,jwk,jacoct)
c
c     stress boundary conditions
c     --------------------------
      if(ipnode==0)
     &   call boundpre(kaxi,nic,nbd,ic,ibdnod,nbound,nside,
     &                 nrdgeom,ndggeom,x,y,nrdvelo,ndgvelo,indgu,
     &                 neq,iperm,imp,neqA,arhs,
     &                 ibdn,densn,shg1,shg1x,we1)
c
      if ( .not.jacoct.and.matglb ) then
c
c        (right)rescale matrix according to reference parameters
c        -------------------------------------------------------
         scalj(1:neqt)=1.0
         do i=1,neqt
           j=jperm(i)  !j is the global variable index
           if(j>=indgp+1.and.j<=indgp+ndgpres*ncppres)then
             scalj(i)=refp
           elseif(j>=indgu+1.and.j<=indgu+ndgvelo*ncpvelo)then
             scalj(i)=refu
           elseif(j>=indgr+1.and.j<=indgr+ndgelas*ncpelas)then
             scalj(i)=refr
           endif
         enddo
c         do i=indgp+1,indgp+ndgpres*ncppres
c            scalj(iperm(i))=refp
c         enddo
c         do i=indgu+1,indgu+ndgvelo*ncpvelo
c            scalj(iperm(i))=refu
c         enddo
c         do i=indgr+1,indgr+ndgelas*ncpelas
c            scalj(iperm(i))=refr
c         enddo
c
         do i=1,neltA
            a(i)=a(i)*scalj(ja(i))
         enddo
c
c        (left)rescale matrix according to the summation of each row
c        -----------------------------------------------------------
         do i=1,neqA
            tmp=0.d0
            do j=ia(i),ia(i+1)-1
               tmp=tmp+abs(a(j))
c              tmp=max(tmp,abs(a(j)))
            enddo
            scali(i)=tmp
         enddo
c
         do i=1,neqA
            a(ia(i):ia(i+1)-1)=a(ia(i):ia(i+1)-1)/scali(i)
         enddo

      endif
c
c     rescale the RHS vector
c     ----------------------
      if(matglb)then
        arhs(1:neqA)=arhs(1:neqA)/scali(1:neqA)
      endif
      newRhs = dnrm2(neqt,arhs,1)
cdebug
c      tmp=scali(1)
c      tmp1=tmp
c      do i=2,neqA
c         tmp=max(tmp,scali(i))
c         tmp1=min(tmp1,scali(i))
c      enddo
c      write(*,*)'scaling: max & min:',tmp,tmp1
c      write(*,*)'neqA',neqA,neqt, ndgphi
c      open(1,file='scaliphi.dat')
c      do i=1,ndgphi
c      j=indgphi+i
c      write(1,'(1x,i5, 3(1x,g12.5))')i,x(i),y(i),scali(iperm(j))
c      enddo
c      close(1)
c      open(1,file='scalipsi.dat')
c      do i=1,ndgphi
c      j=indgpsi+i
c      write(1,'(1x,i5, 3(1x,g12.5))')i,x(i),y(i),scali(iperm(j))
c      enddo
c      close(1)
c      pause   
c      
cdebugend
c
c     Inexact line search scheme
c     --------------------------
      if(backtrack.eq.1)then
      if (searchiter.eq.0.and.newRhs.ge.oldRhs)then
c      if (searchiter.eq.0.and.newRhs.gt.(1.-t*(1.-eta))*oldRhs)then
         linesearch = .true.
         jacoct1 = .true.
         searchiter = searchiter+1
         p0=oldRHS**2
         p1=newRHS**2
         tmp=1.-dtheta
         do i=1,neqt
            j = jperm(i)
            z(j) = z(j) - tmp*vec(i)
         enddo
         goto 10
      endif
      if(searchiter.eq.1)then
         searchiter=searchiter+1
c     calcuate p'(0)
         dp=(newRHS**2-p0)/dtheta
c     calcuate theta wich minimizes the 2nd order polymial p(x)
         if(dp.ge.0.d0)then
            theta=theta_min
         else
            theta=max(theta_min,min(theta_max,-dp/(2.*(p1-p0-dp))))
         endif
         tmp=theta-dtheta
         do i=1,neqt
            j = jperm(i)
            z(j) = z(j) + tmp*vec(i)
         enddo
         linesearch=.false.
c         jacoct1 = jacot
         jacoct1 = .false.
         eta=1.-theta*(1.-eta)
         write(*,'(A,g12.5)')"backtracking invoked, theta=", theta
         goto 10
      endif
      elseif(backtrack.eq.2)then
      if ( newRhs.gt.oldRhs .and. searchiter.lt.5) then
         linesearch = .true.
         jacoct1=.true.
         searchiter = searchiter+1
         tauline    = 0.5d0*tauline
!         write(*,*) tauline
         do i=1,neqt
            j = jperm(i)
            z(j) = z(j) - tauline*vec(i)
         enddo
         goto 10
      endif
      if (linesearch) then
         write(*,'(A,I2)')"line search iter=", searchiter
         if(.not.jacot)then
c            jacoct1=jacot
            jacoct1 = .false.
            linesearch=.false.
            goto 10
         endif
      endif
      endif
c
      write(*,'(A,2e12.5)') 'nonlinear residual= ',newRhs,smax
c      if ( newRhs.lt.epsu*smax .or. newRhs.lt.epsu*intRhs ) then
      if ( newRhs.lt.epsu*max(intRhs,min(smax,sqrt(real(neqt))))) then
c      if ( newRhs.lt.epsu*max(intRhs,sqrt(real(neqt)))) then
         conv = .true.
         write(*,*) 'Newton iteration converged'
      endif
      if ( conv ) return
c
c     Inexact-Newton forcing term control
c     -- from Eisenstat and Walker, SISC, 17, 1 Jan. 1996
c     ---------------------------------------------------
      if ( oldRhs .gt. 1.d15 ) then
c       in the first iteration, oldRHS=1.d20      
         eta = eta_0
      else
         safegrd = eta**phi
         eta = dabs(newRhs-res)/oldRhs
         if ( safegrd.gt.0.1d0) then
            eta=max(eta,safegrd)
c            eta=min(eta_max, max(eta,safegrd))
c            write(6,*) 'Using safeguard eta'
         endif
         eta=min(eta,eta_max)
c         write(6,*) 'New eta: ',eta, safegrd
      endif
      
      oldRhs = newRhs
      eps=max(eta, eta_min)*newRhs
c      write(6,*) 'Linear system tolerance: ',eps
c
c     solver
c     ------

      do i=1,neqt
         vec(i)  = 0.d0
      enddo
      neltv = maxrwk-2*neqt-1
      call solveru (neqA,a,ja,ia,neqB,b,jb,ib,arhs,vec,neltA,neltB,neqt,
     &              ioptu,maxitu,iout,eps,im,lfil,jacoct,ierr,lfilter,
     &              rwork,maxrwk,iwork,maxiwk,neltv,res,globits, 
     &              ioptp,maxitp,epsp,smaxp, sx,isx,jsx,neltS )
      if (ierr.gt.1 ) then
         write(*,*) 'error in solveru, ierr=', ierr
         stop
      endif
C
      if ( ioptu.eq.0 .and. .not.jacoct ) write(iout, '(A,i8)')
     &        'Length of LDU matrix for global variables is',neltv
c
c     check of convergence
c     --------------------
      dsmax=dnrm2(neqt,vec,1)
      if(matglb) vec(1:neqt)=vec(1:neqt)*scalj(1:neqt)
      do i=1,neqt
         j = jperm(i)
         z(j) = z(j) + vec(i)
      enddo
      if(matglb)then
c        calculate the L2 norm of z/scal
         do i=1,neqt
            arhs(i)=z(jperm(i))/scalj(i)
         enddo
      else
         do i=1,neqt
            arhs(i)=z(jperm(i))
         enddo
      endif
      smax  = dnrm2(neqt,arhs,1)
      dsmax = dsmax/(smax+1.d-20)
c
c     update the velocity on particle surface
c     ---------------------------------------
      nwall = nic-nbrigid
      do np=1,nbrigid
         nv = indgnb + ncpb*(np-1) + 1
         xp = xpos(1,np)
         yp = xpos(2,np)
         do i=ic(np+nwall),ic(np+nwall+1)-1
            nd = ibdnod(i)
            dx = xcoor(xp,x(nd)) - xp
            dy = y(nd) - yp
            z(indgu+nd)         = z(nv)-z(nv+2)*dy
            z(indgu+ndgvelo+nd) = z(nv+1)+z(nv+2)*dx
         enddo
      enddo
c
      return
      end
