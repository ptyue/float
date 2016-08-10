c***********************************************************************
      subroutine matsld(m,ishm,nelem,inod,ngmx,x,y,ncpg,z,zs,
     &                   ro,grav,pelas, trans,tsc,jacoct,
     &                   SLVvelo,SLVpres, 
     &                   ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &                   nrdgeom,nrdvelo,nrdpres,nrdelas,nrdtemp,nrdphi,
     &                   ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &                   indgr,indgu,indgp,indgtm,indgphi,indgpsi, 
     &                   matglb,
     &                   neq,neqA,neqB,nvel,neltA,neltB,iglo,ntot,
     &                   a,ja,ia,b,jb,ib,arhs,imp,iperm,jwk,rh,st,
     &                   shg2,shg2x,shg2y,we2, stress,ncps,dt )
c
c     this routine generates element matrix for elastic solid material
c       using Lagrangian frame (Umesh=U)
c 
c     by Hanping Xu and Howard Hu,  Jan. 6, 1999
c
c     iax =2/3  may not work!
c     Modified Feb 18, 2005, by Pengtao Yue
c        to match the arguments of subprograms(FORMiglo), not tested
c***********************************************************************
      implicit none 
      integer m,ishm,nelem,ngmx,neq,inod(ngmx,nelem)
      logical trans,jacoct, SLVvelo, SLVpres, matglb
      integer ncpelas,ncpvelo,ncppres,ncptemp,ncpphi,
     &        nrdgeom,nrdvelo,nrdpres,nrdelas,nrdtemp,nrdphi,
     &        ndggeom,ndgvelo,ndgelas,ndgpres,ndgtemp,ndgphi,
     &        indgu,indgp,indgtm,indgr,indgphi,indgpsi
      integer nvel,ntot,ncpg,ncps,neqA,neqB,neltA,neltB,iglo(ntot),
     &        ia(neqA+1),ja(neltA),imp(neq),iperm(neq),jwk(neq), 
     &        ib(neqB+1),jb(neltB)
      real*8 x(ndggeom),y(ndggeom), z(neq),zs(neq),
     &       tsc,ro,grav(2),pelas(10), a(neltA),b(neltB),arhs(neq)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      real*8 rh(ntot),st(ntot,ntot)
      real*8 stress(ncps,ndgelas),dt
c
c=======================================================================
c         0.   Declarations 
c=======================================================================
      real*8 pu(9),dpudx(9),dpudy(9)
c-----------------------------------------------------------------------
c     working variables
      integer  i, j, l, k,n,ndloc, n1,n2,n3,n4
      real*8 x0,xcoor
      real*8 tmp,tmp1,tmp2
      real*8 dxxi,dxet,dyxi,dyet,dxix,dxiy,detx,dety,
     &       alwejc,jcbn,rjcbn
      real*8 ul,vl,dudxl,dudyl,dvdxl,dvdyl,txxl,txyl
     &      ,usl,vsl,pl, psk,psj
      real*8 resu, resu1, resu2, resv, resv1, resv2, resp
      real*8 roalwe,vsalwe,vsalwt, cfx,cfy
c
      integer ndlgeom,ndlvelo,ndlpres,ndlelas,ndltemp,ndlphi
      integer indu,indv,indp, nnl
      real*8 xn(9),yn(9),u(9),v(9),p(9),us(9),vs(9),txx(9),txy(9)
c
      integer nk,ik,im,irowst,ilast,jk,ir, nvel1
c
c=======================================================================
c     1. load information for the element
c     -----------------------------------
c     number of local nodes for each group of variables
      ndlgeom = ndloc(ishm, nrdgeom)
      ndlvelo = ndloc(ishm, nrdvelo)
      ndlpres = ndloc(ishm, nrdpres)
      ndlelas = ndloc(ishm, nrdelas)
      ndltemp = ndloc(ishm, nrdtemp)
      ndlphi  = ndloc(ishm, nrdphi )
c
c     shifts for numbering in local elements
      if ( SLVvelo ) then
         indu = 0
         indv = indu + ndlvelo
      endif
      if ( SLVpres ) indp = ncpvelo*ndlvelo
c
c     local values of varaibles
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
c     velocity
      nnl = 0
      if ( ncpvelo.ne.0 ) then
         n = indgu
         do i=1,ndlvelo
            k = n + inod(i,m)
            u(i) = z(k)
            us(i) = zs(k)
            k = k + ndgvelo
            v(i) = z(k)
            vs(i) = zs(k)
        enddo
c
c       stress field
        do i=1,ndlelas
           k = inod(i,m)
           txx(i) = stress(1,k)
           txy(i) = stress(2,k)
        enddo
      endif
c
c     pressure
      nnl = nnl + ncpvelo*ndlvelo
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
c
c=======================================================================
c     2. Initialization.
c=======================================================================
      n1 = ishm*2 + nrdgeom
      n2 = ishm*2 + nrdvelo
      n4 = ishm*2 + nrdelas
      if ( nrdpres.eq.0 ) then
         n3 = 6
      else
         n3 = ishm*2 + nrdpres
      endif
c
      do i=1,ntot
         rh(i) = 0.d0
      enddo
      if (.not.jacoct) then
        do i=1,ntot
           do j=1,ntot
              st(i,j) = 0.d0
           enddo
        enddo
      endif
c
c=======================================================================
c     3. loop on integration points
c=======================================================================
      do 30 l=1,7+2*ishm
c
c       geometry
c       --------
        dxxi = 0.d0
        dxet = 0.d0
        dyxi = 0.d0
        dyet = 0.d0
        do i=1,ndlgeom
           dxxi = dxxi + xn(i)*shg2x(i,l,n1)
           dxet = dxet + xn(i)*shg2y(i,l,n1)
           dyxi = dyxi + yn(i)*shg2x(i,l,n1)
           dyet = dyet + yn(i)*shg2y(i,l,n1)
        enddo
        jcbn  = dxxi*dyet - dxet*dyxi
        rjcbn = 1.0/jcbn
        dxix =  dyet*rjcbn 
        dxiy = -dxet*rjcbn
        detx = -dyxi*rjcbn
        dety =  dxxi*rjcbn
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
          enddo
c
          txxl  = 0.d0
          txyl  = 0.d0
          do i=1,ndlelas
             txxl  = txxl  + txx(i)*shg2(i,l,n4)
             txyl  = txyl  + txy(i)*shg2(i,l,n4)
          enddo
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
c       Jacobian and other expressions
c       ------------------------------
        alwejc = jcbn * we2(l,n1)
        roalwe = ro * alwejc
        vsalwe = 0.5*pelas(1)/(1+pelas(2)) * alwejc
        vsalwt = vsalwe*dt
c
c=======================================================================
c       elemental vectors (residuals of the equations)
c=======================================================================
c
        if ( SLVvelo ) then
c         momentum equations
c         ------------------
          resu  = -grav(1)*roalwe
          resu1 = (txxl-pl)*alwejc + 2.d0*dudxl*vsalwt
          resu2 = txyl*alwejc + (dudyl+dvdxl)*vsalwt
          resv  = -grav(2)*roalwe
          resv1 = resu2
          resv2 = (-txxl-pl)*alwejc + 2.d0*dvdyl*vsalwt
          if ( trans ) then
             resu = resu + ( ul*tsc + usl ) * roalwe
             resv = resv + ( vl*tsc + vsl ) * roalwe
          endif
c
          do i=1,ndlvelo
             tmp  = pu(i)
             tmp1 = dpudx(i)
             tmp2 = dpudy(i)
             rh(indu+i) = rh(indu+i) - tmp1*resu1-tmp2*resu2-tmp*resu 
             rh(indv+i) = rh(indv+i) - tmp1*resv1-tmp2*resv2-tmp*resv
          enddo
        endif
c
        if ( SLVpres ) then
c         continuity equation
c         -------------------
          resp = -(dudxl+dvdyl) * alwejc
          do i=1,ndlpres
             rh(indp+i) = rh(indp+i) - shg2(i,l,n3)*resp
          enddo
        endif
c
c=======================================================================
c       elemental matrices (derivatives of the equations)
c=======================================================================
        if( jacoct ) go to 30
c
        if ( SLVvelo ) then
c
c         Momentum equations for velocities
c         ---------------------------------
c         derivatives w.r.t. pressure P
          do i=1,ndlvelo
            psj = -dpudx(i) * alwejc
            psk = -dpudy(i) * alwejc
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
c         derivatives w.r.t. velocities U & V
          do j=1,ndlvelo
            cfx = dpudx(j)*vsalwt
            cfy = dpudy(j)*vsalwt
            do i=1,ndlvelo
              st(indu+i,indu+j) = st(indu+i,indu+j) + dpudx(i)*2.0*cfx
     &                                              + dpudy(i)*cfy
              st(indu+i,indv+j) = st(indu+i,indv+j) + dpudy(i)*cfx
              st(indv+i,indu+j) = st(indv+i,indu+j) + dpudx(i)*cfy
              st(indv+i,indv+j) = st(indv+i,indv+j) + dpudx(i)*cfx
     &                                              + dpudy(i)*2.0*cfy
            enddo
          enddo
c
          if ( trans ) then
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
        endif
c            SLVvelo
c
30    continue
c
c     ----------
c     assembly
c     -----------
      do i=1,nvel
         nk = iglo(i) 
         ik = iperm(nk)
         arhs(ik) = arhs(ik) + rh(i)
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
        nk = iglo(i)
        ik = iperm(nk)
        im = imp(ik)
        if ( im.eq.0 ) then
c
           irowst = ia(ik)
           ilast  = ia(ik+1)-1
           do k=irowst,ilast
              jwk(ja(k)) = k
           enddo

           do j=1,nvel1
              nk = iglo(j)
              jk = iperm(nk)
              k = jwk(jk)
c              if ( k.eq.0 ) then
c                write(*,'(A,3i5)') 'error in assmbl',m,i,ik
c                stop
c              endif
              a(k) = a(k) + st(i,j)
           enddo

           do k=irowst,ilast
             jwk(ja(k)) = 0
           enddo
c
        else
c
           ir = ia(ik)
           a(ir) = 1.0d0
        endif
c
      enddo
c
      do i=nvel1+1,nvel
        nk = iglo(i)
        ik = iperm(nk)
        im = imp(ik)
        if ( im.eq.0 ) then
c
           irowst = ib(ik-neqA)
           ilast  = ib(ik-neqA+1)-1
           do k=irowst,ilast
              jwk(jb(k)) = k
           enddo

           do j=1,nvel1
              nk = iglo(j)
              jk = iperm(nk)
              k = jwk(jk)
              b(k) = b(k) + st(i,j)
           enddo

           do k=irowst,ilast
              jwk(jb(k)) = 0
           enddo
        endif
c
      enddo
c
100   continue
      return
      end
