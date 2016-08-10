c***********************************************************************
      subroutine matamphi(m,ishm,iax,kaxi,nelem,inod,ngmx,x,y,
     &                  nbrigid,ncpb,xpos,umesh,ncpg,z,zs,
     &                  ivisc,itrel,ro,grav,pvis,pelas,pterm, 
     &                  trans,tsc,jacoct,inert,newton,upwind,
     &                  SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &                  ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                  nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &                  ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &                  indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                  matglb,neq,neqA,neqB,nvel,neltA,neltB,iglo,ntot, 
     &                  a,ja,ia,b,jb,ib,arhs,imp,iperm,jwk,rh,st, 
     &                  shg2,shg2x,shg2y,we2,
     &                  phgam,phlam,pheps,phsft,
     &                  flowtype,dropone)
c
c     by Howard Hu, 3/20/94
c     last modified, 8/18/97, iglo is calculated locally
c                    6/23/99, separate A and B matrices
c          for saddle-point problem A & B:  | A  B^t |
c                                           | B  0   |
c     Modified by Pengtao Yue, Feb 9, 2005
c        add the interfacial tension terms
c     Modified by Pengtao Yue, Feb 17, 2005
c        Cahn-Hillard Eq is added in a fully coupled manner
c     Modified by Pengtao Yue, Mar 22-Apr 8, 2005
c        subroutine name is changed from matnew to matamphi. It's used 
c        to calculate interfacial flows with different viscosities and
c        densities. (The derivatives of viscosities w.r.t. velocities 
c        etc. are not tested.)
c     Modified by Pentao Yue, Jul 28, 2005
c        modifications are made on the body force term for drop/interface
c        impact calculations.
c     Modified by P. Yue, 1/14/2015
c        different forms of surface tension are incorporated
c        surftens=1, body force G*\nabla \phi
c                 2, tensor form fI-\lambda*(nabla\phi otimes \nabla\phi)
c                 3, tensor form -\lambda*(nabla\phi otimes \nabla\phi)
c***********************************************************************
      implicit none 
      integer, parameter :: surftens=1
      integer m,ishm,nelem, ngmx,inod(ngmx,nelem),neq, ivisc(2),itrel
      logical trans,jacoct,inert,newton
      logical SLVvelo, SLVpres, SLVtemp, SLVelas, upwind,matglb,SLVphi
      integer ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &        nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdphi,
     &        ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &        indgu,indgp,indgtm,indgr,indgnb, indgphi,indgpsi,iax,kaxi
      integer nvel,nbrigid,ntot,neqA,neqB,neltA,neltB,ncpb,ncpg,
     &        iglo(ntot),ia(neqA+1),ja(neltA),imp(neq),iperm(neq),
     &        jwk(neq),ib(neqB+1),jb(neltB)
      real*8 x(ndggeom),y(ndggeom),umesh(ncpg,ndgvelo), z(*),
     &       zs(*),tsc,ro(2),grav(2),pvis(10,2),pelas(10),pterm(10),
     &       xpos(ncpb,nbrigid),a(neltA),b(neltB),arhs(neq)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      real*8 rh(ntot),st(ntot,ntot)
c     variables for drop/interface impact calculations
      integer flowtype
      logical dropone
c
      real*8 phgam,phlam,pheps,phsft
c
c=======================================================================
c         0.   Declarations 
c=======================================================================
c     aUPLWD - =1/=0: upper/lower convected derivatives.
c     epsPTT - eps in Phan-Thien & Tanner model.
c     alfaGL - alfa in Giesekus model.
      real*8 aUPLWD, epsPTT, alfaGL
c
      real*8 ps(9),dpsdx(9),dpsdy(9),wps(9)
      real*8 pu(9),dpudx(9),dpudy(9),
     &       dvscdu(9),dvscdv(9),dlmbdu(9),dlmbdv(9) 
      real*8 pf(9),dpfdx(9),dpfdy(9)
c-----------------------------------------------------------------------
c     working variables
      integer  i, j, l, k,n,ndloc, n1,n2,n3,n4, n5,i1,i2,j1,j2
      real*8 x0,xcoor
      real*8 tmp,tmp1,tmp2,tmp3,tmp4
      real*8 xx,yy,dxxi,dxet,dyxi,dyet,dxix,dxiy,detx,dety,
     &       alwejc,jcbn,rjcbn,alpha,xsmfct 
      real*8 rl,sl,tl,ql, ul,vl, pl,
     &       rsl,ssl,tsl,qsl,usl,vsl,
     &       drdxl,drdyl,dsdxl,dsdyl,dtdxl,dtdyl,dqdxl,dqdyl,
     &       dudxl,dudyl,dvdxl,dvdyl
      real*8 resr,ress,rest,resq,
     &       resu, resu1, resu2, resv, resv1, resv2, resp
c     Coef. for material properties
      real*8 dvscds,dvscdt,ratio,cratio,phitm,lambda,lambdt,
     &       dlmbds,aplow,apupp,pttelv,ptty,argexp,gksalv,dsn(3,3),
     &       visct1,visct2, appml,apupp2,aplow2,appml2,ulxsm,roalwe,
     &       cfrdr,cfrds,cfrdt,cfrdq,cfsdr,cfsds,cfsdt,cfsdq,cftdr,
     &       cftds,cftdt,cftdq,cfqdr,cfqds,cfqdq,rlalwe,vsalwe1,vsalwe2,
     &       cfrdu,cfrdv,cfsdu,cfsdv,cftdu,cftdv,cfqdu,cfqdv,psj,psk,
     &       cfudu1,cfudu2,cfudu3,cfudv,cfvdu,cfvdv1,cfvdv2,cfvdv3,
     &       cfudv3,cfvdu3,visc1p,visc1s,visc1,visc2
      real*8 um,vm,ulm,vlm,unorm,ubar,vbar,tau,txu,tyu,
     &       dtxdu,dtydu,dtxdv,dtydv
c
      integer ndlgeom,ndlvelo,ndlpres,ndltemp,ndlelas
      integer indu,indv,indp,indtm,indr,inds,indt,indq,indphi,indpsi,nnl
      real*8 xn(9),yn(9),u(9),v(9),p(9),temp(9),
     &       r(9),s(9),t(9),q(9),us(9),vs(9),temps(9),
     &       rs(9),ss(9),ts(9),qs(9), umsh(9),vmsh(9),
     &       f(9),g(9),fs(9),gs(9)
c
      integer in,ik,im,jn,jk,jm,nd,np,np1,irowst,ilast,
     &        iv,iv1,iv2,jv,nvel1
      real*8 xp,yp,ulmdu,vlmdv
c
      real*8 fl,gl,fsl,gsl,dfdxl,dfdyl,dgdxl,dgdyl,eps2,ga1,ga2
      integer ndlphi
c
      real(8)  a1,a2,c1,c2,tmpu1,tmpu2,tmpu3,tmpv2
C     Support for Petsc timing
c      include '../include/petsc.h'
c
c      call PLogEventBegin(ELEM_MAT_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c=======================================================================
c     0. Global parameters
c=======================================================================
      epsPTT = 0.d0
c     Gieskus model
      alfaGL = pelas(4)
c
      aUPLWD = 1.d0
c
      ratio  = pelas(6)
      apupp  = aUPLWD
      aplow  = 1.d0 - apupp
      appml  = apupp - aplow
      apupp2 = 2.d0 * apupp
      aplow2 = 2.d0 * aplow
      appml2 = 2.d0 * appml
c
      eps2=pheps**2
      ga1=phlam/eps2
      ga2=phgam*ga1
      a1=(ro(1)-ro(2))/2.
c
c     1. load information for the element
c     -----------------------------------
c     number of local nodes for each group of variables
      ndlgeom = ndloc(ishm, nrdgeom)
      ndlvelo = ndloc(ishm, nrdvelo)
      ndlpres = ndloc(ishm, nrdpres)
      ndltemp = ndloc(ishm, nrdtemp)
      ndlelas = ndloc(ishm, nrdelas)
      ndlphi  = ndloc(ishm, nrdphi)
c
c     shifts for numbering in local elements
      nnl=0
      if ( SLVelas ) then
         indr = nnl
         inds = indr + ndlelas
         indt = inds + ndlelas
         if ( iax.eq.1 ) indq = indt + ndlelas
      endif
      nnl=nnl+ncpelas*ndlelas
      if ( SLVvelo ) then
         indu = nnl
         indv = indu + ndlvelo
      endif
      nnl=nnl+ncpvelo*ndlvelo
      if ( SLVtemp ) indtm = nnl
      nnl=nnl+ncptemp*ndltemp
      if ( SLVphi )then
         indphi = nnl
         indpsi = indphi + ndlphi
      endif
      nnl=nnl+ncpphi*ndlphi*2
      if ( SLVpres ) indp=nnl
c
c     local values of varaibles (velo., pres, temp, stresses)
      x0 = x(inod(1,m))
      do i=1,ndlgeom
         k = inod(i,m)
         xn(i) = xcoor(x0,x(k))
         yn(i) = y(k)
      enddo
c
      call FORMiglo ( iglo,nvel, m,inod,ngmx,nrdpres,
     &                ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                ndgelas,ndgvelo,ndgpres,ndgtemp,
     &                ndlelas,ndlvelo,ndlpres,ndltemp,ndlphi,
     &                indgr,indgu,indgp,indgtm,indgphi,indgpsi)
c
c     elastic-stress
      if ( ncpelas.ne.0) then
        n = indgr
        do i=1,ndlelas
          k = n + inod(i,m)
          r(i) = z(k)
          rs(i) = zs(k)
          k = k + ndgelas
          s(i) = z(k)
          ss(i) = zs(k)
          k = k + ndgelas
          t(i) = z(k)
          ts(i) = zs(k)
          if ( iax.eq.1 ) then
             k = k + ndgelas
             q(i) = z(k)
             qs(i) = zs(k)
          endif
        enddo
      endif
c
c     velocity
      if ( ncpvelo.ne.0 ) then
        n = indgu
        do i=1,ndlvelo
          k = inod(i,m)
          umsh(i) = umesh(1,k)
          vmsh(i) = umesh(2,k)
          k = n + inod(i,m)
          u(i) = z(k)
          us(i) = zs(k)
          k = k + ndgvelo
          v(i) = z(k)
          vs(i) = zs(k)
c          if ( iax.eq.2 ) then
c            k = k + ndgvelo
c            w(i) = z(k)
c            ws(i) = zs(k)
c          endif
        enddo
      endif
c
c     temperature
      if ( ncptemp.ne.0) then
        n = indgtm
        do i=1,ndltemp
           k = n + inod(i,m)
           temp(i) = z(k)
           temps(i) = zs(k)
        enddo
      endif
c
c     phase-field variables
      if (ncpphi.eq.1)then
         do i=1,ndlphi
            k=iglo(indphi+i)
            f(i)=z(k)
            fs(i)=zs(k)
            k=iglo(indpsi+i)
            g(i)=z(k)
            gs(i)=zs(k)
         enddo
      endif
c
c     pressure
      if ( ncppres.ne.0) then
        n = indgp
        do i=1,ndlpres
          if ( nrdpres.eq.0 ) then
             k = n + m
          else
             k = n + inod(i,m)
          endif
          p(i) = z(k)
        enddo
      endif
c=======================================================================
c     2. Initialization.
c=======================================================================
      n1 = ishm*2 + nrdgeom
      n2 = ishm*2 + nrdvelo
      n4 = ishm*2 + nrdelas
      n5 = ishm*2 + nrdphi
c
      if ( nrdpres.eq.0 ) then
         n3 = 6
      else
         n3 = ishm*2 + nrdpres
      endif

      rh(1:ntot) = 0.d0
      if (.not.jacoct) st(1:ntot,1:ntot)=0.d0
c
c=======================================================================
c     3. loop on integration points
c=======================================================================
      do 30 l=1,7+2*ishm
c
c       geometry
c       --------
        xx   = 0.d0
        yy   = 0.d0
        dxxi = 0.d0
        dxet = 0.d0
        dyxi = 0.d0
        dyet = 0.d0
        do i=1,ndlgeom
           xx   = xx   + xn(i)*shg2(i,l,n1)
           yy   = yy   + yn(i)*shg2(i,l,n1)
           dxxi = dxxi + xn(i)*shg2x(i,l,n1)
           dxet = dxet + xn(i)*shg2y(i,l,n1)
           dyxi = dyxi + yn(i)*shg2x(i,l,n1)
           dyet = dyet + yn(i)*shg2y(i,l,n1)
        enddo
        jcbn  = dxxi*dyet - dxet*dyxi
        rjcbn = 1.0/jcbn
        alpha = 1.d0 + kaxi* (xx - 1.d0)
        xsmfct = kaxi / alpha
        dxix =  dyet*rjcbn 
        dxiy = -dxet*rjcbn
        detx = -dyxi*rjcbn
        dety =  dxxi*rjcbn
c
c       Tensor Sij
c       ----------
        if ( SLVelas ) then
          rl    = 0.d0
          drdxl = 0.d0
          drdyl = 0.d0
          sl    = 0.d0
          dsdxl = 0.d0
          dsdyl = 0.d0
          tl    = 0.d0
          dtdxl = 0.d0
          dtdyl = 0.d0
          do i=1,ndlelas
            ps(i) = shg2(i,l,n4)
            dpsdx(i) = shg2x(i,l,n4)*dxix + shg2y(i,l,n4)*detx
            dpsdy(i) = shg2x(i,l,n4)*dxiy + shg2y(i,l,n4)*dety
            rl    = rl    + r(i)*ps(i)
            drdxl = drdxl + r(i)*dpsdx(i)
            drdyl = drdyl + r(i)*dpsdy(i)
            sl    = sl    + s(i)*ps(i)
            dsdxl = dsdxl + s(i)*dpsdx(i)
            dsdyl = dsdyl + s(i)*dpsdy(i)
            tl    = tl    + t(i)*ps(i)
            dtdxl = dtdxl + t(i)*dpsdx(i)
            dtdyl = dtdyl + t(i)*dpsdy(i)
          enddo
c
          if ( trans ) then
            rsl = 0.d0
            ssl = 0.d0
            tsl = 0.d0
            do i=1,ndlelas
              rsl = rsl + rs(i)*ps(i)
              ssl = ssl + ss(i)*ps(i)
              tsl = tsl + ts(i)*ps(i)
            enddo
          endif
c
            ql    = 0.d0
            dqdxl = 0.d0
            dqdyl = 0.d0
          if( iax.eq.1 ) then
            do i=1,ndlelas
              ql    = ql    + q(i)*ps(i)
              dqdxl = dqdxl + q(i)*dpsdx(i)
              dqdyl = dqdyl + q(i)*dpsdy(i)
            enddo
            if ( trans ) then
              qsl = 0.d0
              do i=1,ndlelas
                qsl = qsl + qs(i)*ps(i)
              enddo
            endif
          endif
c
        endif
c
c       velocities
c       ----------
        if ( SLVvelo ) then
          ul    = 0.d0
          dudxl = 0.d0
          dudyl = 0.d0
          vl    = 0.d0
          dvdxl = 0.d0
          dvdyl = 0.d0
          um = 0.d0
          vm = 0.d0
          do i=1,ndlvelo
            pu(i) = shg2(i,l,n2)
            dpudx(i) = shg2x(i,l,n2)*dxix + shg2y(i,l,n2)*detx
            dpudy(i) = shg2x(i,l,n2)*dxiy + shg2y(i,l,n2)*dety
            ul    = ul    + u(i)*pu(i)
            dudxl = dudxl + u(i)*dpudx(i)
            dudyl = dudyl + u(i)*dpudy(i)
            vl    = vl    + v(i)*pu(i)
            dvdxl = dvdxl + v(i)*dpudx(i)
            dvdyl = dvdyl + v(i)*dpudy(i)
            um    = um    + umsh(i)*pu(i)
            vm    = vm    + vmsh(i)*pu(i)
          enddo
          ulxsm = ul * xsmfct
          ulm = ul - um
          vlm = vl - vm
c
c          ulmdu = dabs( ulm/(ul+1.d-40) )
c          vlmdv = dabs( vlm/(vl+1.d-40) )
c          if ( ulmdu.gt.1.d0 ) ulmdu = 1.d0
c          if ( vlmdv.gt.1.d0 ) vlmdv = 1.d0
          ulmdu = 1.d0
          vlmdv = 1.d0
c
          if ( trans ) then
            usl = 0.d0
            vsl = 0.d0
            do i=1,ndlvelo
              usl = usl + us(i)*pu(i)
              vsl = vsl + vs(i)*pu(i)
            enddo
          endif
c
c         rate of strain tensor
c         ---------------------
          dsn(1,1) = 2.d0*dudxl
          dsn(1,2) = dudyl + dvdxl
          dsn(1,3) = 0.d0
          dsn(2,2) = 2.d0*dvdyl
          dsn(2,3) = 0.d0
          dsn(3,3) = 2.d0*ulxsm
c
        else
          ulm=0.d0
          vlm=0.d0
        endif
c
c        Phase-field variables
c        ---------------------
        if (SLVphi) then

         do i=1,ndlphi
            pf(i) = shg2(i,l,n5)
            dpfdx(i) = shg2x(i,l,n5)*dxix + shg2y(i,l,n5)*detx
            dpfdy(i) = shg2x(i,l,n5)*dxiy + shg2y(i,l,n5)*dety
         enddo
         fl    = sum(f(1:ndlphi)* pf  (1:ndlphi))
         dfdxl = sum(f(1:ndlphi)*dpfdx(1:ndlphi))
         dfdyl = sum(f(1:ndlphi)*dpfdy(1:ndlphi))
         gl    = sum(g(1:ndlphi)* pf  (1:ndlphi))
         dgdxl = sum(g(1:ndlphi)*dpfdx(1:ndlphi))
         dgdyl = sum(g(1:ndlphi)*dpfdy(1:ndlphi))
         fsl   = 0.d0
         gsl   = 0.d0
         if(trans)then
            fsl=sum(fs(1:ndlphi)*pf(1:ndlphi))
         endif  
         c1=min(1.d0,max(0.d0,(1.d0+fl)/2.d0))
         c2=1.-c1
        else
	   c1=1.
         c2=0.
        endif  
c
c       pressure
c       --------
        if ( SLVpres ) then
          pl = 0.
          do i=1,ndlpres
             pl = pl + p(i)*shg2(i,l,n3)
          enddo
        endif
c
c       temperature
c       -----------
c       NOT IMPLEMENTED YET!

c=======================================================================
c       4. evaluation of viscosity and its derivatives
c=======================================================================
c     fluid one(phi=1, may be non-newtonian fluids)
        call viscos( ivisc(1),pvis(1,1),dsn,1,pterm,0.d0,visc1,
     &               dvscds,dvscdt )
        cratio = 1.d0 - ratio
        visc1p  = visc1 * cratio
        visc1s  = visc1 * ratio
c     fluid two(phi=-1, can only be (generalized) newtonian fluids)
        call viscos( ivisc(2),pvis(1,2),dsn,1,pterm,0.d0,visc2,
     &               dvscds,dvscdt )
c
c       evalution of relaxation time and its derivatives
c       ------------------------------------------------
        call relx( itrel, pelas, dsn, lambda, dlmbds )
c
c       Influence of temperature, not used!!
c       ------------------------------------
        phitm = 1.d0
        lambdt = lambda * phitm
        visct1 = visc1p * phitm
        visct2 = (visc1s*c1+visc2*c2) * phitm
        a2=(visc1s-visc2)/2.
c
c       Jacobian and other expressions
c       ------------------------------
        alwejc  = alpha * jcbn * we2(l,n1)
        roalwe  = (c1*ro(1)+c2*ro(2))* alwejc
        vsalwe1 = visct1 * alwejc
        vsalwe2 = visct2 * alwejc
        rlalwe  = lambdt * alwejc
c
c       Upwinding weighting functions
c       -----------------------------
        do i=1,ndlelas
          wps(i) = ps(i)
        enddo
        if ( upwind ) then
          unorm = 1.d0 / (dsqrt( ulm*ulm + vlm*vlm ) + 1.d-20)
          ubar = ulm*unorm
          vbar = vlm*unorm
          tmp1 = ubar*dxix + vbar*dxiy
          tmp2 = ubar*detx + vbar*dety
          tau = 1.d0 / (dsqrt( tmp1*tmp1 + tmp2*tmp2 )+1.d-20)
          do i=1,ndlelas
            wps(i) = wps(i) + tau * ( ubar*dpsdx(i) + vbar*dpsdy(i) )
          enddo
          tau = tau*tau*tau
          txu = ( dxiy*tmp1 + dety*tmp2 )*tau
          tyu = ( dxix*tmp1 + detx*tmp2 )*tau
          dtxdu =  txu*vbar*unorm
          dtydu = -tyu*vbar*unorm
          dtxdv = -txu*ubar*unorm
          dtydv =  tyu*ubar*unorm
        endif
c
c=======================================================================
c       5. elemental vectors (residuals of the equations)
c=======================================================================
c
        if ( SLVelas ) then
c         constitutive equations
c         ----------------------
          if ( visc1p.ne.0.d0 ) then
            pttelv = epsPTT * lambda / visc1p
            argexp = pttelv * (rl+sl)
            if ( iax.eq.1 ) argexp = pttelv * (rl+sl+ql)
            gksalv = alfaGL * lambda / visc1p
            if ( argexp.eq.0 ) then
              pttelv = 0.d0
              ptty   = 1.d0
            elseif ( argexp.gt.50.d0 ) then
              stop 'Overflow in PTT fluid!'
            else
              ptty   = dexp ( argexp )
              pttelv = pttelv * ptty
            endif
          else
            pttelv = 0.d0
            ptty   = 1.d0
            gksalv = 0.d0
          endif
c
          resr = ( ptty * rl
     &           + gksalv * ( rl*rl + tl*tl )
     &           - lambdt * apupp2 * ( dudxl * rl + dudyl * tl )
     &           + lambdt * aplow2 * ( dudxl * rl + dvdxl * tl ) 
     &           - visct1*dsn(1,1)
     &           + lambdt * ( ulm * drdxl + vlm * drdyl ) ) * alwejc
          ress = ( ptty * sl
     &           + gksalv * ( sl*sl + tl*tl )
     &           - lambdt * apupp2 * ( dvdxl * tl + dvdyl * sl )
     &           + lambdt * aplow2 * ( dudyl * tl + dvdyl * sl ) 
     &           - visct1*dsn(2,2)
     &           + lambdt * ( ulm * dsdxl + vlm * dsdyl ) ) * alwejc
          rest = ( ptty * tl
     &           + gksalv * ( rl + sl ) * tl
     &           - lambdt * apupp * ( (dudxl + dvdyl) * tl
     &                               + dudyl * sl + dvdxl * rl )
     &           + lambdt * aplow * ( (dudxl + dvdyl) * tl
     &                               + dudyl * rl + dvdxl * sl )
     &           - visct1*dsn(1,2)
     &           + lambdt * ( ulm * dtdxl + vlm * dtdyl ) ) * alwejc
c
          if ( trans ) then
            resr = resr + (rl*tsc + rsl) * rlalwe
            ress = ress + (sl*tsc + ssl) * rlalwe
            rest = rest + (tl*tsc + tsl) * rlalwe
          endif
c
          do i = 1, ndlelas
            rh(indr+i) = rh(indr+i) - wps(i)*resr
            rh(inds+i) = rh(inds+i) - wps(i)*ress
            rh(indt+i) = rh(indt+i) - wps(i)*rest
          enddo
c
          if( iax.eq.1 ) then
            resq = ( ptty * ql
     &           + gksalv * ql*ql
     &           - lambdt * appml2 * ulxsm * ql
     &           - visct1*dsn(3,3)
     &           + lambdt * ( ulm * dqdxl + vlm * dqdyl ) ) * alwejc
            if ( trans ) resq = resq + (ql*tsc + qsl) * rlalwe
            do i=1,ndlelas
              rh(indq+i) = rh(indq+i) - wps(i)*resq
            enddo
          endif
c
        endif
c
        if ( SLVvelo ) then
c         momentum equations
c         ------------------
          resu  = visct2*dsn(3,3) - pl
          resu1 = visct2*dsn(1,1) - pl
          resu2 = visct2*dsn(1,2)
          resv2 = visct2*dsn(2,2) - pl
          if ( SLVelas ) then
            if (iax.eq.1) resu  = resu + ql*c1
            resu1 = resu1 + rl*c1
            resu2 = resu2 + tl*c1
            resv2 = resv2 + sl*c1
          endif
          resu  = resu * xsmfct * alwejc
          resu1 = resu1 * alwejc
          resu2 = resu2 * alwejc
          resv  = 0.d0
          resv1 = resu2
          resv2 = resv2 * alwejc
c
         if(flowtype.eq.106.or.flowtype.eq.116.or.flowtype.eq.107.or.
     &      flowtype.eq.117.or.flowtype.eq.124)then
c        rho_matrix*g*h is taken into pressure
            if(dropone)then
c              phase 2 (phi=-1) is the matrix phase
               tmp=c1*(ro(1)-ro(2))*alwejc
               resu = resu - grav(1)*tmp
               resv = resv - grav(2)*tmp
            else
c              phase 1 (phi=+1) is the matrix phase
               tmp=c2*(ro(2)-ro(1))*alwejc
               resu = resu - grav(1)*tmp
               resv = resv - grav(2)*tmp
            endif      
         else
          resu = resu - grav(1)*roalwe
          resv = resv - grav(2)*roalwe
         endif
c
          if ( inert ) then
            resu = resu + ( ulm*dudxl + vlm*dudyl) * roalwe
            resv = resv + ( ulm*dvdxl + vlm*dvdyl) * roalwe
          endif
          if ( trans ) then
            resu = resu + ( ul*tsc + usl ) * roalwe
            resv = resv + ( vl*tsc + vsl ) * roalwe
          endif
c     interfacial tension
          if(SLVphi)then
            select case (surftens)
              case(1)
c surface tension form 1: G*nabla phi              
                tmp=ga1*(gl+phsft*fl)*alwejc
                resu = resu - tmp * dfdxl
                resv = resv - tmp * dfdyl
              case(2)
c   surface tension form 2: tensor form fI-lambda*(nabla phi nabla phi)
c   P. Yue, 1/13/2016
                tmp=(fl**2-1.d0)**2/(4.d0*pheps**2)
     &             +0.5*(dfdxl**2+dfdyl**2)
                tmp1=phlam*(dfdxl*dfdyl)*alwejc
                resu1=resu1+phlam*(tmp-dfdxl**2)*alwejc
                resu2=resu2-tmp1
                resv1=resv1-tmp1
                resv2=resv2+phlam*(tmp-dfdyl**2)*alwejc
                if(iax==1)resu=resu+tmp*phlam*xsmfct*alwejc
              case(3)
c   surface tension form 3: tensor form -lambda*(nabla phi nabla phi)
c   in this form the pressure term also contains contribution from 
c   mixing energy. P. Yue, 1/13/2016
                tmp=phlam*alwejc
                resu1=resu1-dfdxl*dfdxl*tmp
                resu2=resu2-dfdxl*dfdyl*tmp
                resv1=resv1-dfdxl*dfdyl*tmp
                resv2=resv2-dfdyl*dfdyl*tmp
            end select
          endif 
c
          do i=1,ndlvelo
            tmp = pu(i)
            tmp1 = dpudx(i)
            tmp2 = dpudy(i)
            rh(indu+i) = rh(indu+i) - tmp1*resu1-tmp2*resu2-tmp*resu 
            rh(indv+i) = rh(indv+i) - tmp1*resv1-tmp2*resv2-tmp*resv
          enddo
        endif
        if (SLVphi)then
c
c        Cahn-Hillard Equation
c        ---------------------
         if(SLVvelo)then
            tmp1=ulm*dfdxl+vlm*dfdyl
         else
            tmp1=0.d0
         endif
         tmp2=(fl**2-1.d0-phsft)*fl-gl
         do i=1,ndlphi
            tmp3=dfdxl*dpfdx(i)+dfdyl*dpfdy(i)
            j=i+indphi
            rh(j)=rh(j)-(tmp1*pf(i)+ga2*(dgdxl*dpfdx(i)
     &           +dgdyl*dpfdy(i)+phsft*tmp3))*alwejc
            j=i+indpsi
            rh(j)=rh(j)+(tmp2*pf(i)+eps2*tmp3)*alwejc
         enddo
         i1=indphi+1
         i2=indphi+ndlphi
         if(trans)then
            tmp1=(tsc*fl+fsl)*alwejc
            rh(i1:i2)=rh(i1:i2)-tmp1*pf(1:ndlphi)
         endif
        endif

        if ( SLVpres ) then
c         continuity equation
c         -------------------
          resp = -( dudxl + dvdyl + ulxsm ) * alwejc
          do i=1,ndlpres
            rh(indp+i) = rh(indp+i) - shg2(i,l,n3)*resp
          enddo
        endif
c
c=======================================================================
c       6. elemental matrices (derivatives of the equations)
c=======================================================================
        if( jacoct ) go to 30
c
        if ( SLVelas ) then
c
c         Constitutive equations for elastic stresses Sij
c         -----------------------------------------------
          cfrdq = pttelv * rl * alwejc
          cfrdr = cfrdq + ( ptty + gksalv * 2.d0 * rl
     &                    - lambdt * appml2 * dudxl ) * alwejc
          cfrds = cfrdq
          cfrdt = ( gksalv * 2.d0 * tl - lambdt * apupp2 * dudyl
     &            + lambdt * aplow2 * dvdxl ) * alwejc
          cfsdq = pttelv * sl * alwejc
          cfsdr = cfsdq
          cfsds = cfsdq + ( ptty + gksalv * 2.d0 * sl 
     &                    - lambdt * appml2 * dvdyl ) * alwejc
          cfsdt = ( gksalv * 2.d0 * tl - lambdt * apupp2 * dvdxl
     &            + lambdt * aplow2 * dudyl ) * alwejc
          cftdq = pttelv * tl * alwejc
          cftdr = cftdq + ( gksalv * tl - lambdt * apupp * dvdxl
     &                    + lambdt * aplow * dudyl ) * alwejc
          cftds = cftdq + ( gksalv * tl - lambdt * apupp * dudyl
     &                    + lambdt * aplow * dvdxl ) * alwejc
          cftdt = ( ptty  + gksalv * (rl + sl)
     &            - lambdt * appml * (dudxl + dvdyl) ) * alwejc
c
c         derivatives w.r.t. stresses Sij
          do j=1,ndlelas
            psj = ps(j)
            psk = ( ulm*dpsdx(j) + vlm*dpsdy(j) ) * rlalwe
            do i=1,ndlelas
              tmp  = wps(i)*psj
              tmp1 = wps(i)*psk
              st(indr+i,indr+j) = st(indr+i,indr+j) + tmp*cfrdr + tmp1
              st(indr+i,inds+j) = st(indr+i,inds+j) + tmp*cfrds
              st(indr+i,indt+j) = st(indr+i,indt+j) + tmp*cfrdt
              st(inds+i,indr+j) = st(inds+i,indr+j) + tmp*cfsdr
              st(inds+i,inds+j) = st(inds+i,inds+j) + tmp*cfsds + tmp1
              st(inds+i,indt+j) = st(inds+i,indt+j) + tmp*cfsdt
              st(indt+i,indr+j) = st(indt+i,indr+j) + tmp*cftdr
              st(indt+i,inds+j) = st(indt+i,inds+j) + tmp*cftds
              st(indt+i,indt+j) = st(indt+i,indt+j) + tmp*cftdt + tmp1
            enddo
          enddo
c
c         derivatives w.r.t. velocities U & V
          do j=1,ndlvelo
            cfrdu = -(appml2*dpudx(j)*rl+apupp2*dpudy(j)*tl)*rlalwe
     &              -2.d0*dpudx(j)*vsalwe1
            cfrdv = aplow2*dpudx(j)*tl * rlalwe
            cfsdu = aplow2*dpudy(j)*tl * rlalwe
            cfsdv = -(apupp2*dpudx(j)*tl+appml2*dpudy(j)*sl)*rlalwe
     &              -2.d0*dpudy(j)*vsalwe1
            cftdu = ( -appml*dpudx(j)*tl
     &             + ( -apupp*sl + aplow*rl )*dpudy(j) ) * rlalwe
     &             -dpudy(j)*vsalwe1
            cftdv = ( -appml*dpudy(j)*tl
     &             + ( -apupp*rl + aplow*sl )*dpudx(j) ) * rlalwe
     &             -dpudx(j)*vsalwe1
            cfrdu = cfrdu + pu(j)*drdxl * rlalwe * ulmdu
            cfrdv = cfrdv + pu(j)*drdyl * rlalwe * vlmdv
            cfsdu = cfsdu + pu(j)*dsdxl * rlalwe * ulmdu
            cfsdv = cfsdv + pu(j)*dsdyl * rlalwe * vlmdv
            cftdu = cftdu + pu(j)*dtdxl * rlalwe * ulmdu
            cftdv = cftdv + pu(j)*dtdyl * rlalwe * vlmdv
            do i=1,ndlelas
              st(indr+i,indu+j) = st(indr+i,indu+j) + wps(i)*cfrdu
              st(indr+i,indv+j) = st(indr+i,indv+j) + wps(i)*cfrdv
              st(inds+i,indu+j) = st(inds+i,indu+j) + wps(i)*cfsdu
              st(inds+i,indv+j) = st(inds+i,indv+j) + wps(i)*cfsdv
              st(indt+i,indu+j) = st(indt+i,indu+j) + wps(i)*cftdu
              st(indt+i,indv+j) = st(indt+i,indv+j) + wps(i)*cftdv
            enddo
          enddo
c
          if ( iax.eq.1 ) then
            cfqdr = pttelv * ql * alwejc
            cfqds = cfqdr
            cfqdq = cfqdr + ( ptty + gksalv * 2.d0 * ql
     &                    - lambdt * appml2 * ulxsm ) * alwejc
            do j=1,ndlelas
              psj = ps(j)
              psk = ( ulm*dpsdx(j) + vlm*dpsdy(j) ) * rlalwe
              do i=1,ndlelas
                tmp  = wps(i)*psj
                tmp1 = wps(i)*psk
                st(indr+i,indq+j) = st(indr+i,indq+j) + tmp*cfrdq
                st(inds+i,indq+j) = st(inds+i,indq+j) + tmp*cfsdq
                st(indt+i,indq+j) = st(indt+i,indq+j) + tmp*cftdq
                st(indq+i,indr+j) = st(indq+i,indr+j) + tmp*cfqdr
                st(indq+i,inds+j) = st(indq+i,inds+j) + tmp*cfqds
                st(indq+i,indq+j) = st(indq+i,indq+j) + tmp*cfqdq+tmp1
              enddo
            enddo
            do j=1,ndlvelo
              cfqdu = pu(j)*appml2*xsmfct*ql * rlalwe
     &               -2.d0*pu(j)*xsmfct*vsalwe1
              cfqdv = 0.
              cfqdu = cfqdu + pu(j) * dqdxl * rlalwe * ulmdu
              cfqdv = cfqdv + pu(j) * dqdyl * rlalwe * vlmdv
              do i=1,ndlelas
                st(indq+i,indu+j) = st(indq+i,indu+j) + wps(i)*cfqdu
                st(indq+i,indv+j) = st(indq+i,indv+j) + wps(i)*cfqdv
              enddo
            enddo
          endif
c
        endif
c            SLVelas
c
c
        if ( SLVvelo ) then
c
c         Momentum equations for velocities
c         ---------------------------------
c         derivatives w.r.t. pressure P
          do i=1,ndlvelo
            psj = -( dpudx(i) + xsmfct * pu(i) ) * alwejc
            psk = - dpudy(i) * alwejc
            do j=1,ndlpres
              tmp1 = shg2(j,l,n3)*psj
              tmp2 = shg2(j,l,n3)*psk
              st(indu+i,indp+j) = st(indu+i,indp+j) + tmp1
              st(indp+j,indu+i) = st(indp+j,indu+i) + tmp1
              st(indv+i,indp+j) = st(indv+i,indp+j) + tmp2
              st(indp+j,indv+i) = st(indp+j,indv+i) + tmp2
            enddo
          enddo
c
          if ( SLVelas ) then
c           derivatives w.r.t. extra-stresses Sij
            do j=1,ndlelas
              psj = ps(j)*alwejc*c1
              do i=1,ndlvelo
                tmp1 = dpudx(i)*psj
                tmp2 = dpudy(i)*psj
                st(indu+i,indr+j) = st(indu+i,indr+j) + tmp1
                st(indv+i,indt+j) = st(indv+i,indt+j) + tmp1
                st(indu+i,indt+j) = st(indu+i,indt+j) + tmp2
                st(indv+i,inds+j) = st(indv+i,inds+j) + tmp2
              enddo
            enddo
            if ( iax.eq.1 ) then
              do j=1,ndlelas
                psj = ps(j)*alwejc*xsmfct*c1
                do i=1,ndlvelo
                  st(indu+i,indq+j) = st(indu+i,indq+j) + pu(i)*psj
                enddo
              enddo
            endif
          endif
c
c         derivatives w.r.t. velocities U & V
          do j=1,ndlvelo
            cfudv = dpudx(j)*vsalwe2
            cfvdu = dpudy(j)*vsalwe2
            cfudu1 = 2.d0*cfudv
            cfudu2 = cfvdu
            cfudu3 = 2.d0*pu(j)*xsmfct*xsmfct*vsalwe2
            cfvdv1 = cfudv
            cfvdv2 = 2.d0*cfvdu
            cfvdv3 = 0.d0
            cfudv3 = 0.d0
            cfvdu3 = 0.d0
            if ( inert ) then
              tmp = ( ulm*dpudx(j)+vlm*dpudy(j) ) * roalwe
              cfudu3 = cfudu3 + tmp + pu(j)*dudxl * roalwe * ulmdu
              cfvdv3 =          tmp + pu(j)*dvdyl * roalwe * vlmdv
              cfudv3 = pu(j)*dudyl * roalwe * vlmdv
              cfvdu3 = pu(j)*dvdxl * roalwe * ulmdu
            endif
            do i=1,ndlvelo
              st(indu+i,indu+j) = st(indu+i,indu+j) + dpudx(i)*cfudu1
     &                               + dpudy(i)*cfudu2 + pu(i)*cfudu3
              st(indu+i,indv+j) = st(indu+i,indv+j) + dpudy(i)*cfudv
     &                                              + pu(i)*cfudv3
              st(indv+i,indu+j) = st(indv+i,indu+j) + dpudx(i)*cfvdu
     &                                              + pu(i)*cfvdu3
              st(indv+i,indv+j) = st(indv+i,indv+j) + dpudx(i)*cfvdv1
     &                               + dpudy(i)*cfvdv2 + pu(i)*cfvdv3
            enddo
          enddo
c
         if(SLVphi)then 
c
c         derivatives w.r.t. phase-field variables
c        -----------------------------------------
c     1.  derivatives due to interfacial tension
            do j=1,ndlphi
              select case (surftens)            
                case(1)
c surface tension form 1: G*nabla phi            
              tmp1=ga1*(phsft*dfdxl*pf(j)+(gl+phsft*fl)*dpfdx(j))*alwejc
              tmp2=ga1*(phsft*dfdyl*pf(j)+(gl+phsft*fl)*dpfdy(j))*alwejc
              tmp3=ga1*dfdxl*pf(j)*alwejc
              tmp4=ga1*dfdyl*pf(j)*alwejc
                  do i=1,ndlvelo
                    st(indu+i,indphi+j)=st(indu+i,indphi+j)-tmp1*pu(i)
                    st(indv+i,indphi+j)=st(indv+i,indphi+j)-tmp2*pu(i)
                    st(indu+i,indpsi+j)=st(indu+i,indpsi+j)-tmp3*pu(i)
                    st(indv+i,indpsi+j)=st(indv+i,indpsi+j)-tmp4*pu(i)
                  enddo
                case(2)
c   surface tension form 2: tensor form fI-lambda*(nabla phi nabla phi)
                  tmp1=(fl**2-1)*fl/(pheps**2)*pf(j)
     &                +(dfdxl*dpfdx(j)+dfdyl*dpfdy(j))
                  tmp2=(tmp1-2*dfdxl*dpfdx(j))*phlam*alwejc
                  tmp3=(tmp1-2*dfdyl*dpfdy(j))*phlam*alwejc
                  tmp4=(dfdxl*dpfdy(j)+dfdyl*dpfdx(j))*phlam*alwejc
                  tmp1=tmp1*phlam*alwejc*xsmfct
                  do i=1,ndlvelo
                  st(indu+i,indphi+j)=st(indu+i,indphi+j)+tmp2*dpudx(i)
     &                                        -tmp4*dpudy(i)+tmp1*pu(i)
                  st(indv+i,indphi+j)=st(indv+i,indphi+j)+tmp3*dpudy(i)
     &                                                   -tmp4*dpudx(i)
                  enddo
                case(3)
c   surface tension form 3: tensor form -lambda*(nabla phi nabla phi)
c   in this form the pressure term also contains contribution from 
c   mixing energy. P. Yue, 1/13/2016
                  tmp=phlam*alwejc
                  tmp1=2.d0*tmp*dfdxl*dpfdx(j)
                  tmp2=2.d0*tmp*dfdyl*dpfdy(j)
                  tmp3=tmp*(dfdxl*dpfdy(j)+dfdyl*dpfdx(j))
                  do i=1,ndlvelo
                  st(indu+i,indphi+j)=st(indu+i,indphi+j)-tmp1*dpudx(i)
     &                                                   -tmp3*dpudy(i)
                  st(indv+i,indphi+j)=st(indv+i,indphi+j)-tmp3*dpudx(i)
     &                                                   -tmp2*dpudy(i)
                  enddo
              end select
            enddo
c     2. derivatives due to density and viscosity
c     inertia and body force terms
            resu=-grav(1)
            resv=-grav(2)
            if ( inert ) then
               resu = resu + ( ulm*dudxl + vlm*dudyl)
               resv = resv + ( ulm*dvdxl + vlm*dvdyl)
            endif
            if ( trans ) then
               resu = resu + ( ul*tsc + usl )
               resv = resv + ( vl*tsc + vsl )
            endif
c     body force
            resu=resu*a1*alwejc
            resv=resv*a1*alwejc
c     stress terms
            if(SLVelas)then
               tmpu1=(a2*dsn(1,1)+rl/2.)*alwejc
               tmpu2=(a2*dsn(1,2)+tl/2.)*alwejc
               tmpu3=(a2*dsn(3,3)+ql/2.)*xsmfct*alwejc
               tmpv2=(a2*dsn(2,2)+sl/2.)*alwejc
            else
               tmpu1=a2*dsn(1,1)*alwejc
               tmpu2=a2*dsn(1,2)*alwejc
               tmpu3=a2*dsn(3,3)*xsmfct*alwejc
               tmpv2=a2*dsn(2,2)*alwejc
            endif
c
            do i=1,ndlvelo
               resu=tmpu1*dpudx(i)+tmpu2*dpudy(i)+(resu+tmpu3)*pu(i)
               resv=tmpu2*dpudx(i)+tmpv2*dpudy(i)+resv*pu(i)
               do j=1,ndlphi
                  st(indu+i,indphi+j)=st(indu+i,indphi+j)+resu*pf(j)
                  st(indv+i,indphi+j)=st(indv+i,indphi+j)+resv*pf(j)
               enddo
            enddo
         endif   
c
        endif
c
        if(SLVphi)then
c
c        Cahn-Hillard Equation
c        ---------------------
c     derivatives w.r.t. phi & psi
         tmp1=(3.d0*fl**2-1.d0-phsft)*alwejc
         do j=1,ndlphi
          tmp2=(ulm*dpfdx(j)+vlm*dpfdy(j))*alwejc
          tmp3=pf(j)*alwejc
          j1=j+indphi
          j2=j+indpsi
         do i=1,ndlphi
            i1=i+indphi
            i2=i+indpsi
            tmp4=(dpfdx(i)*dpfdx(j)+dpfdy(i)*dpfdy(j))*alwejc
            st(i1,j1)=st(i1,j1)+(pf(i)*tmp2+phsft*ga2*tmp4)
            st(i1,j2)=st(i1,j2)+ga2*tmp4
            st(i2,j1)=st(i2,j1)-(tmp1*pf(i)*pf(j)+eps2*tmp4)
            st(i2,j2)=st(i2,j2)+pf(i)*tmp3
         enddo
         enddo
         if(SLVvelo)then
c     derivatives w.r.t. u & v
            i1=indphi+1
            i2=indphi+ndlphi
            do j=1,ndlvelo
               tmp1=dfdxl*ulmdu*pu(j)*alwejc
               tmp2=dfdyl*vlmdv*pu(j)*alwejc
               j1=indu+j
               j2=indv+j
               st(i1:i2,j1)=st(i1:i2,j1)+pf(1:ndlphi)*tmp1
               st(i1:i2,j2)=st(i1:i2,j2)+pf(1:ndlphi)*tmp2
            enddo
         endif
        endif 
c            SLVvelo
c
c
c       transient terms
c       ---------------
        if ( trans ) then
          if ( SLVelas ) then
            do j=1,ndlelas
              psj = ps(j)*rlalwe*tsc
              do i=1,ndlelas
                tmp = wps(i)*psj
                st(indr+i,indr+j) = st(indr+i,indr+j) + tmp
                st(inds+i,inds+j) = st(inds+i,inds+j) + tmp
                st(indt+i,indt+j) = st(indt+i,indt+j) + tmp
              enddo
            enddo
c
            if ( iax.eq.1 ) then
              do j=1,ndlelas
                psj = ps(j)*rlalwe*tsc
                do i=1,ndlelas
                  st(indq+i,indq+j) = st(indq+i,indq+j) + wps(i)*psj
                enddo
              enddo
            endif
          endif
c
          if ( SLVvelo ) then
            do j=1,ndlvelo
              psj = pu(j)*roalwe*tsc
              do i=1,ndlvelo
                tmp = pu(i)*psj
                st(indu+i,indu+j) = st(indu+i,indu+j) + tmp
                st(indv+i,indv+j) = st(indv+i,indv+j) + tmp
              enddo
            enddo
          endif
c
         if ( SLVphi)then
            i1=indphi+1
            i2=indphi+ndlphi
            tmp1=tsc*alwejc
            do j=1,ndlphi
               tmp2=tmp1*pf(j)
               j1=indphi+j
               st(i1:i2,j1)=st(i1:i2,j1)+tmp2*pf(1:ndlphi)
            enddo
         endif            
        endif
c            trans
c
c
c       Newton-Raphson for White-Metzner fluids
c       ---------------------------------------
c     This part is not tested in the diffuse-interface calculation
        if ( newton ) then
c
c         derivatives of viscosity
          if ( ivisc(1).ne.1 ) then
            do i=1,ndlvelo
c             derivatives of viscosity w.r.t. velocities U & V
              dvscdu(i) = dvscds*phitm*0.5*( dsn(1,1)*dpudx(i)
     &                                   + dsn(1,2)*dpudy(i)
     &                                   + dsn(3,3)*pu(i)*xsmfct )
              dvscdv(i) = dvscds*phitm*0.5*( dsn(1,2)*dpudx(i)
     &                                   + dsn(2,2)*dpudy(i) )
            enddo
c
            if ( SLVelas ) then
c             Equations for elastic stress Sij
c             --------------------------------
              do i=1,ndlelas
                tmp1 = ps(i) * dsn(1,1) * cratio * alwejc
                tmp2 = ps(i) * dsn(2,2) * cratio * alwejc
                tmp3 = ps(i) * dsn(1,2) * cratio * alwejc
                do j=1,ndlvelo
                  st(indr+i,indu+j)=st(indr+i,indu+j) - tmp1*dvscdu(j)
                  st(indr+i,indv+j)=st(indr+i,indv+j) - tmp1*dvscdv(j)
                  st(inds+i,indu+j)=st(inds+i,indu+j) - tmp2*dvscdu(j)
                  st(inds+i,indv+j)=st(inds+i,indv+j) - tmp2*dvscdv(j)
                  st(indt+i,indu+j)=st(indt+i,indu+j) - tmp3*dvscdu(j)
                  st(indt+i,indv+j)=st(indt+i,indv+j) - tmp3*dvscdv(j)
                enddo
              enddo
              if ( iax.eq.1 ) then
                do i=1,ndlelas
                  tmp4 = ps(i) * dsn(3,3) * cratio * alwejc
                  do j=1,ndlvelo
                    st(indq+i,indu+j)=st(indq+i,indu+j)-tmp4*dvscdu(j)
                    st(indq+i,indv+j)=st(indq+i,indv+j)-tmp4*dvscdv(j)
                  enddo
                enddo
              endif
            endif
c
            if ( SLVvelo ) then
c             Momentum equations 
c             ------------------
              do i=1,ndlvelo
                tmp1 = ratio*alwejc * ( dpudx(i)*dsn(1,1)
     &                                + dpudy(i)*dsn(1,2)
     &                                + pu(i)*dsn(3,3)*xsmfct )
                tmp2 = ratio*alwejc * ( dpudx(i)*dsn(1,2)
     &                                + dpudy(i)*dsn(2,2) )
                do j=1,ndlvelo
                  st(indu+i,indu+j) = st(indu+i,indu+j) + tmp1*dvscdu(j)
                  st(indu+i,indv+j) = st(indu+i,indv+j) + tmp1*dvscdv(j)
                  st(indv+i,indu+j) = st(indv+i,indu+j) + tmp2*dvscdu(j)
                  st(indv+i,indv+j) = st(indv+i,indv+j) + tmp2*dvscdv(j)
                enddo
              enddo
            endif
c
          endif
c
c
c         derivatives of relaxation time
c         ------------------------------
          if ( itrel.ne.1 .and. SLVelas ) then
            do i=1,ndlvelo
c             derivatives of ralaxtion time w.r.t. velocities
              dlmbdu(i) = dlmbds*0.5*( dsn(1,1)*dpudx(i)
     &                               + dsn(1,2)*dpudy(i)
     &                               + dsn(3,3)*pu(i)*xsmfct )
              dlmbdv(i) = dlmbds*0.5*( dsn(1,2)*dpudx(i)
     &                               + dsn(2,2)*dpudy(i) )
            enddo
c
c           Equations for extra-stress Sij
            do i=1,ndlelas
              tmp1 = wps(i) * ( ( ulm*drdxl + vlm*drdyl ) 
     &                - apupp2 * ( dudxl*rl + dudyl*tl )
     &                + aplow2 * ( dudxl*rl + dvdxl*tl ) ) * alwejc
              tmp2 = wps(i) * ( ( ulm*dsdxl + vlm*dsdyl )
     &                - apupp2 * ( dvdxl*tl + dvdyl*sl )
     &                + aplow2 * ( dudyl*tl + dvdyl*sl ) ) * alwejc
              tmp3 = wps(i) * ( ( ulm*dtdxl + vlm*dtdyl )
     &           - apupp * ( (dudxl+dvdyl)*tl + dudyl*sl + dvdxl*rl )
     &           + aplow * ( (dudxl+dvdyl)*tl + dudyl*rl + dvdxl*sl )
     &                       ) * alwejc
              tmp4 = wps(i) * ( ( ulm*dqdxl + vlm*dqdyl )
     &                        - appml2*ulxsm*ql ) * alwejc
              if ( trans ) then         
                tmp1 = tmp1 + wps(i) * ( rl*tsc + rsl ) * alwejc
                tmp2 = tmp2 + wps(i) * ( sl*tsc + ssl ) * alwejc
                tmp3 = tmp3 + wps(i) * ( tl*tsc + tsl ) * alwejc
                tmp4 = tmp4 + wps(i) * ( ql*tsc + qsl ) * alwejc
              endif
              do j=1,ndlvelo
                st(indr+i,indu+j) = st(indr+i,indu+j) + tmp1*dlmbdu(j)
                st(indr+i,indv+j) = st(indr+i,indv+j) + tmp1*dlmbdv(j)
                st(inds+i,indu+j) = st(inds+i,indu+j) + tmp2*dlmbdu(j)
                st(inds+i,indv+j) = st(inds+i,indv+j) + tmp2*dlmbdv(j)
                st(indt+i,indu+j) = st(indt+i,indu+j) + tmp3*dlmbdu(j)
                st(indt+i,indv+j) = st(indt+i,indv+j) + tmp3*dlmbdv(j)
                if ( iax.eq.1 ) then
                  st(indq+i,indu+j) = st(indq+i,indu+j) + tmp4*dlmbdu(j)
                  st(indq+i,indv+j) = st(indq+i,indv+j) + tmp4*dlmbdv(j)
                endif
              enddo
            enddo
c
          endif
c
          if ( upwind ) then
            do i=1,ndlelas
              tmp1 = ( dtxdu*dpsdx(i) + dtydu*dpsdy(i) )*ulmdu
              tmp2 = ( dtxdv*dpsdx(i) + dtydv*dpsdy(i) )*vlmdv
              do j=1,ndlvelo
                st(indr+i,indu+j) = st(indr+i,indu+j) + resr*tmp1*pu(j)
                st(indr+i,indv+j) = st(indr+i,indv+j) + resr*tmp2*pu(j)
                st(inds+i,indu+j) = st(inds+i,indu+j) + ress*tmp1*pu(j)
                st(inds+i,indv+j) = st(inds+i,indv+j) + ress*tmp2*pu(j)
                st(indt+i,indu+j) = st(indt+i,indu+j) + rest*tmp1*pu(j)
              enddo
            enddo
            if ( iax.eq.1 ) then
              do i=1,ndlelas
                tmp1 = ( dtxdu*dpsdx(i) + dtydu*dpsdy(i) )*ulmdu
                tmp2 = ( dtxdv*dpsdx(i) + dtydv*dpsdy(i) )*vlmdv
                do j=1,ndlvelo
                  st(indq+i,indu+j)=st(indq+i,indu+j) + resq*tmp1*pu(j)
                  st(indq+i,indv+j)=st(indq+i,indv+j) + resq*tmp2*pu(j)
                enddo
              enddo
            endif
          endif
c
        endif
c            newton
c
30    continue
c
c     8. assembly
c     -----------
      do i=1,nvel
        in = iglo(i)
        ik = iperm(in)
        im = imp(ik)
        if ( im.eq.0 ) then
           arhs(ik) = arhs(ik) + rh(i)
        elseif ( im.lt.0 ) then
           nd = mod(in-indgu,ndgvelo)
           np = -im
           iv = indgnb + ncpb*(np-1) + 1
           if ( i.le.indv ) then
              yp = xpos(2,np)
              arhs(iv)   = arhs(iv)   + rh(i)
              arhs(iv+2) = arhs(iv+2) - rh(i)*(y(nd)-yp)
           else
              xp = xpos(1,np)
              tmp = xcoor(xp,x(nd))-xp
              arhs(iv+1) = arhs(iv+1) + rh(i)
              arhs(iv+2) = arhs(iv+2) + rh(i)*tmp
           endif
        endif
      enddo
c
      if ( jacoct ) goto 100
c
      if ( matglb ) then
         nvel1 = nvel
      else
         nvel1 = nvel - ncppres*ndlpres
      endif
c
      do i=1,nvel1
        in = iglo(i)
        ik = iperm(in)
        im = imp(ik)
        if ( im.eq.0 ) then
c
c          for normal interior variables
c          -----------------------------
           irowst = ia(ik)
           ilast  = ia(ik+1)-1
           do k=irowst,ilast
              jwk(ja(k)) = k
           enddo
c
           do j=1,nvel1
              jn = iglo(j)
              jk = iperm(jn)
              jm = imp(jk)
              if ( jm.eq.0 ) then
                 k = jwk(jk)
                 a(k) = a(k) + st(i,j)
              elseif ( jm.lt.0 ) then
                 nd = mod(jn-indgu,ndgvelo)
                 np1 = -jm
                 jv  = indgnb + ncpb*(np1-1) + 1
                 if ( j.le.indv ) then
                    yp  = xpos(2,np1)
                    k = jwk(jv)
                    a(k) = a(k) + st(i,j)
                    k = jwk(jv+2)
                    a(k) = a(k) - st(i,j)*(y(nd)-yp)
                 else
                    xp  = xpos(1,np1)
                    k = jwk(jv+1)
                    a(k) = a(k) + st(i,j)
                    k = jwk(jv+2)
                    tmp =  (xcoor(xp,x(nd))-xp)
                    a(k) = a(k) + st(i,j)*tmp
                 endif
              endif
           enddo
c
           do k=irowst,ilast
              jwk(ja(k)) = 0
           enddo
c
        elseif ( im.lt.0 ) then
c
c          velocities on moving particles
c          ------------------------------
           nd = mod(in-indgu,ndgvelo)
           np = -im
           iv = indgnb + ncpb*(np-1) + 1
           iv2 = iv+2
           if ( i.le.indv ) then
              tmp = -(y(nd)-xpos(2,np))
              iv1 = iv
           else
              xp = xpos(1,np)
              tmp =  (xcoor(xp,x(nd))-xp)
              iv1 = iv+1
           endif
c
           irowst = ia(iv1)
           ilast  = ia(iv1+1)-1
           do k=irowst,ilast
              jwk(ja(k)) = k
           enddo
c
           do j=1,nvel1
              jn = iglo(j)
              jk = iperm(jn)
              jm = imp(jk)
              if ( jm.eq.0 ) then
                 k = jwk(jk)
                 a(k) = a(k) + st(i,j)
              elseif ( jm.lt.0 ) then
                 nd = mod(jn-indgu,ndgvelo)
                 np1 = -jm
                 jv  = indgnb + ncpb*(np1-1) + 1
                 if ( j.le.indv ) then
                    yp  = xpos(2,np1)
                    k = jwk(jv)
                    a(k) = a(k) + st(i,j)
                    k = jwk(jv+2)
                    a(k) = a(k) - st(i,j)*(y(nd)-yp)
                 else
                    xp  = xpos(1,np1)
                    k = jwk(jv+1)
                    a(k) = a(k) + st(i,j)
                    k = jwk(jv+2)
                    tmp1 = (xcoor(xp,x(nd))-xp)
                    a(k) = a(k) + st(i,j)*tmp1
                 endif
              endif
           enddo
c
           do k=irowst,ilast
              jwk(ja(k)) = 0
           enddo
c
           irowst = ia(iv2)
           ilast  = ia(iv2+1)-1
           do k=irowst,ilast
              jwk(ja(k)) = k
           enddo
c
           do j=1,nvel1
              jn = iglo(j)
              jk = iperm(jn)
              jm = imp(jk)
              if ( jm.eq.0 ) then
                 k = jwk(jk)
                 a(k) = a(k) + st(i,j)*tmp
              elseif ( jm.lt.0 ) then
                 nd = mod(jn-indgu,ndgvelo)
                 np1 = -jm
                 jv  = indgnb + ncpb*(np1-1) + 1
                 if ( j.le.indv ) then
                    yp  = xpos(2,np1)
                    k = jwk(jv)
                    a(k) = a(k) + st(i,j)*tmp
                    k = jwk(jv+2)
                    a(k) = a(k) - st(i,j)*(y(nd)-yp)*tmp
                 else
                    xp  = xpos(1,np1)
                    k = jwk(jv+1)
                    a(k) = a(k) + st(i,j)*tmp
                    k = jwk(jv+2)
                    tmp1 = (xcoor(xp,x(nd))-xp)
                    a(k) = a(k) + st(i,j)*tmp1*tmp
                 endif
              endif
           enddo
c
           do k=irowst,ilast
              jwk(ja(k)) = 0
           enddo
c
        endif
      enddo
c
      do i=nvel1+1,nvel
        in = iglo(i)
        ik = iperm(in)
        im = imp(ik)
        if ( im.eq.0 ) then
c
           irowst = ib(ik-neqA)
           ilast  = ib(ik-neqA+1)-1
           do k=irowst,ilast
              jwk(jb(k)) = k
           enddo
c
           do j=1,nvel1
              jn = iglo(j)
              jk = iperm(jn)
              jm = imp(jk)
              if ( jm.eq.0 ) then
                 k = jwk(jk)
                 b(k) = b(k) + st(i,j)
              elseif ( jm.lt.0 ) then
                 nd = mod(jn-indgu,ndgvelo)
                 np1 = -jm
                 jv  = indgnb + ncpb*(np1-1) + 1
                 if ( j.le.indv ) then
                    yp  = xpos(2,np1)
                    k = jwk(jv)
                    b(k) = b(k) + st(i,j)
                    k = jwk(jv+2)
                    b(k) = b(k) - st(i,j)*(y(nd)-yp)
                 else
                    xp  = xpos(1,np1)
                    k = jwk(jv+1)
                    b(k) = b(k) + st(i,j)
                    k = jwk(jv+2)
                    tmp =  (xcoor(xp,x(nd))-xp)
                    b(k) = b(k) + st(i,j)*tmp
                 endif
              endif
           enddo
c
           do k=irowst,ilast
              jwk(jb(k)) = 0
           enddo
        endif
c
      enddo
c
C     Support for Petsc timing
100   continue
c      call PLogEventEnd(ELEM_MAT_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
