c***********************************************************************
      subroutine streamf( strm, iout,kaxi,
     &                    nrdstrm,nrdgeom,nrdvelo,
     &                    ndgstrm,ndggeom,ndgvelo,
     &                    nelem,inod,ngmx,x,y,ishape, u,v,
     &                    nbd,ibdnod,nic,ic,nbound,nside, isnode,
     &                    shg2,shg2x,shg2y,we2, shg1,shg1x,we1,
     &                    iopt,lfil,maxit,eps,
     &                    a,ja,mxnelt,ia,iperm,jperm,arhs,
     &                    iwk,jwk,iwork,maxiwk,rwork,maxrwk )
c
c    package for streamfunction calculation
c      (neuman boundary conditions are used d(fee)/dn = -v*nx+u*ny )
c
c     INPUT
c       kaxi = index for geometry, 0: planar; 1: axisymmetric
c       iout = output unit number
c       nrdstrm = interpolation order for stream function
c       nrdgeom = order of interpolation for coordinates
c       nrdvelo = order of interpolation for velocities
c       ndgstrm = number of nodes for stream function
c       ndggeom = number of nodes in the mesh
c       ndgvelo = number of nodes for the velocities
c       nelem = number of elements
c       inod(ngmx,nelem) = element description table
c       x(ndggeom),y(ndggeom) = mesh coordinates
c       ishape(nelem) = element shape, 0:triangular; 1=quadrilateral
c       nbd = number of boundary nodes
c       ibdnod(nbd) = index for boundary nodes
c       nic   : number of boundary sections
c       ic    : pointer for the boundary sections
c       nbound,nside(nbound),nic,ic = boundary information
c       isnode = node number where stram function is zero
c       iopt    = option for the solver
c         lfil  = number of fill-ins in the ilut preconditioner
c         maxit = maximum number of iterations allowed in GMRES
c         eps   = tolerance for stopping criterion
c       mxnelt = maximum length of matrix
c
c     OUTPUT 
c       strm(ndgstrm) = stream function
c
c     WORKING VARIABLES
c       a(mxnelt),ja(mxnelt),ia(ndgstrm) = matrix in CSR format
c       iperm(ndgstrm) : iperm(global var#) = order of elimination
c       jperm(ndgstrm) : see below
c       arhs(ndgstrm) = right-hand-side vector
c       iwk(ndgstrm),jwk(ndgstrm),iwork(maxiwk),rwork(maxrwk)
c
c     Howard Hu, Jan 27, 1996
c     
c     Jun 10, 2005, Pengtao Yue
c        extended to axisymmetric geometry. Note that r*psi instead of 
c        the stream function psi itself is output. So that the contours
c        of strm coincide with streamlines in axisymmetric cases.
c***********************************************************************
      implicit none
      integer nrdgeom,ndggeom,nrdvelo,ndgvelo,nrdstrm,ndgstrm,mxnelt,
     &        iout,nbd,nelem,nbound,kaxi, isnode,ishape(nelem),
     &        ngmx,inod(ngmx,nelem),ibdnod(nbd),
     &        nside(nbound),ia(ndgstrm+1),ja(mxnelt),
     &        iperm(ndgstrm),jperm(ndgstrm),nic,ic(nic+1),
     &        iopt,lfil,maxit, maxiwk,maxrwk
      real*8 eps,x(ndggeom),y(ndggeom),u(ndgvelo),v(ndgvelo),a(mxnelt),
     &       arhs(ndgstrm), strm(ndgstrm),rwork(maxrwk)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      real*8 shg1(3,3,2),shg1x(3,3,2),we1(3,2)
      integer iwk(ndgstrm+1),jwk(ndgstrm),iwork(maxiwk)
c
c     local variables
      integer i,m,ishm,ndloc,ndlgeom,ndlvelo,ndlstrm,ierr,
     &        nd,nelt,lenrw,leniw,neltv
      real*8 fact
c
c     1. calculation of the boundary integral for stream function
c     -----------------------------------------------------------
      call prstrm( kaxi, nrdstrm,ndgstrm,nrdgeom,ndggeom,ndgvelo,
     &             nelem,inod,x,y,ishape,ngmx,
     &             nbd,ibdnod,nic,ic,nbound,nside, isnode,
     &             u,v,strm,rwork, ia,ja,iperm,jperm,nelt,mxnelt,
     &             iwk,jwk,iwork, iout, shg1,shg1x,we1 )
c
      do i=1,ndgstrm
         jwk(i) = 0
         arhs(i) = 0.d0
      enddo
      do i=1,nelt
         a(i) = 0.d0
      enddo
c
c     2. construction of the stiffness matrix and triangulation
c     ---------------------------------------------------------
      ishm = ishape(1)
      ndlgeom = ndloc(ishm,nrdgeom)
      ndlstrm = ndloc(ishm,nrdstrm)
      ndlvelo = ndloc(ishm,nrdvelo)
      do m=1,nelem
        call matstr(m,ishm,kaxi,nrdgeom,ndggeom,ndlgeom,nrdstrm,ndgstrm,
     &             ndlstrm,nrdvelo,ndgvelo,ndlvelo,nelem,inod, ngmx,x,y,
     &          u,v,nelt,a,ja,ia, strm,iperm,jwk, shg2,shg2x,shg2y,we2 )
      enddo
      if(kaxi.eq.1)then
        call matstr_axi(nrdstrm,ndgstrm,nrdgeom,ndggeom,x,y,
     &                  nbd,ibdnod,nic,ic,nbound,nside,
     &                  nelt,a,ia,ja,iperm,jwk,shg1,shg1x,we1)
      endif
c
c     scaling
      fact = 1.d0 + kaxi*(2.d0*3.1415926536d0-1.d0)
      do i=1,ndgstrm
         strm(i) = strm(i)*fact
      enddo
c
c     boundary condition
      nd = iperm(isnode)
      a(ia(nd)) = 1.d50
      strm(nd)  = 0.d0
c
c     3. solver
c     ---------
      neltv = maxrwk-2*ndgstrm-1
      call solver (arhs, ndgstrm,strm,a,ja,ia,nelt,
     &             iopt,lfil,maxit,eps,iout,.false.,ierr,
     &             rwork,maxrwk,iwork,maxiwk,neltv )
      if ( iopt.eq.0 ) write(iout,'(A,i8)')
     &              'Length of LDU matrix for stream function is ',neltv
c
c     4. copy the result
c     ------------------
      do i=1,ndgstrm
         strm(i) = arhs(iperm(i))
      enddo
      if(kaxi.eq.1)strm(1:ndgstrm)=strm(1:ndgstrm)*x(1:ndgstrm)
c
c      do i=1,nbd
c         m = ibdnod(i)
c
c Check if strm is correct on boundary 
c
c     if ( m.le.ndgstrm ) write(*,222) i,m,x(m),y(m),strm(m)
c 222    format('Bound. node=',i3,' Node=',i6,' (x,y)=',2F9.4,
c    +          ' Stream Func. =',F10.6)
c
c      enddo
      return
      end
c
c***********************************************************************
      subroutine prstrm( kaxi, nrdstrm,ndgstrm,nrdgeom,ndggeom,ndgvelo,
     &                   nelem,inod,x,y,ishape,ngmx,
     &                   nbd,ibdnod,nic,ic,nbound,nside, isnode,
     &                   u,v,rhloc,valstr, 
     &                   ia,ja,iperm,jperm,nelt,mxnelt,
     &                   iwk,jwk,iwork, iout, shg1,shg1x,we1 )
c
c     this routine calculates the boundary integral for stream function.
c***********************************************************************
      implicit none
      integer nbd,nelem,ngmx,ishape(nelem),inod(ngmx,nelem),ibdnod(nbd),
     &        isnode,nbound,kaxi,nic,ic(nic+1),
     &        nside(nbound),nrdstrm,ndgstrm,nrdgeom,ndggeom,ndgvelo,
     &        mxnelt,ia(ndgstrm+1),ja(mxnelt),nelt,iperm(ndgstrm),
     &        jperm(ndgstrm),iwk(ndgstrm+1),jwk(ndgstrm),
     &        iwork(3*nelem), iout
      real*8 u(ndgvelo),v(ndgvelo),x(ndggeom),y(ndggeom),rhloc(ndgstrm),
     &       valstr(ndgstrm), shg1(3,3,2),shg1x(3,3,2),we1(3,2)
c
      integer nidx,nd,ndlstrm,ndloc,i,j,k,kk,ne,n,ishm,nb,
     &        m,iside,ica1,ica,icb,nel,ib,ibm1,m0,mi, is
      real*8 xn(3),yn(3),b1dx(3,3),b1dy(3,3),b1ds(3,3),c(3), xcoor
c
c     1. initialization
c     -----------------
      nidx = nrdgeom
      nd = nidx+1
      do i=1,nbd+1
        valstr(i) = 0.d0
      enddo
c
c     2. loop over each boundary section
c     ----------------------------------
      iside = 0
      do nb=1,nbound
        ica1 = ic(iside+1)
        do kk=1,nside(nb)
          iside = iside + 1
          ica = ic(iside)
          icb = ic(iside+1)
          nel = (icb-ica)/nidx
c
c         3. loop over each boundary segment
c         ----------------------------------
          do m=1,nel
            ib = ica + nidx*(m-1)
            ibm1 = ib - 1
c
            n = ibdnod(ib)
            xn(1) = x(n)
            yn(1) = y(n)
            do i=2,nd
              k = ibm1 + i
              if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
              n = ibdnod(k)
              xn(i) = xcoor(xn(1),x(n))
              yn(i) = y(n)
            enddo
c
            call mat1ds (xn,yn,nidx,kaxi,b1dx,b1dy,b1ds,shg1,shg1x,we1)
c
c           3.1 calculating boundary integrals
            do i=1,nd
              c(i) = 0.
            enddo
            do i=1,nd
              do j=1,nd
                k = ibm1 + j
                if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
                n = ibdnod(k)
                c(i) = c(i) - b1dx(j,i)*u(n) - b1dy(j,i)*v(n)
              enddo
            enddo
c
c           3.2 distributing the boundary integrals
            if ( nrdstrm.eq.nrdgeom ) then
              do i=1,nd
                k = ibm1 + i
                if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
                valstr(k) = valstr(k) + c(i)
              enddo
            endif
            if ( nrdstrm.eq.1 .and. nrdgeom.eq.2 ) then 
              do i=1,3,2
                k = ibm1 + i
                if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
                valstr(k) = valstr(k) + c(i) + 0.5*c(2)
              enddo
            endif
c
          enddo
c
        enddo
      enddo
c
c     5. computation of the order of elimination (permutation index)
c     --------------------------------------------------------------
c     calculate iwk(global var#) = last element# related
      ishm = ishape(1)
      ndlstrm = ndloc(ishm,nrdstrm)
      do ne=1,nelem
         do j=1,ndlstrm
           i= inod(j,ne)
           iwk(i) = ne
         enddo
      enddo
c     iperm(global var#) = order of elimination
c     jperm(order of elimination) = global var#
      n = 0
      do ne=1,nelem
         do j=1,ndlstrm
            i = inod(j,ne)
            if ( iwk(i).eq.ne ) then
               n = n + 1
               iperm(i) = n
               jperm(n) = i
            endif
         enddo
      enddo
c
c     load rhloc
      do i=1,ndgstrm
        rhloc(i) = 0.d0
      enddo
      do i=1,nbd
         n = ibdnod(i)
         if ( n.le.ndgstrm ) then
           k = iperm(n)
           rhloc(k) = rhloc(k) + valstr(i)
         endif
c
c Check if the Neumann B.C. (surf. integn.) has been done correctly
c      write(*,222) i,n,x(n),y(n),valstr(i)
c  222 format('i=',i3,' n=',i3,' (x,y)=',2F7.3,' valstr=',F10.4)
c
      enddo
c
c     6. preparation of storage of the matrix in the CSR format
c     --------------------------------------------------------
      do i=1,ndgstrm+1
         iwk(i) = 0
      enddo
      do i=1,nelem
         do j=1,ndlstrm
            is = inod(j,i)+1
            iwk(is) = iwk(is) + 1
         enddo
      enddo
      do i=2,ndgstrm+1
         iwk(i)= iwk(i-1) + iwk(i)
      enddo

      do i=1,nelem
         do j=1,ndlstrm
            is = inod(j,i)
            iwk(is) = iwk(is) + 1
            iwork(iwk(is)) = i
         enddo
      enddo
      do i=ndgstrm+1,2,-1
         iwk(i) = iwk(i-1)
      enddo
      iwk(1) = 0
c
      do i=1,ndgstrm
        jwk(i) = 0
      enddo
      m = 0
      ia(1) = 1
      do i=1,ndgstrm
         nd = jperm(i)
c
         m = m + 1
         ja(m) = i
         jwk(i) = i 
         m0 = m+1
         do j=iwk(nd)+1,iwk(nd+1)
            ne = iwork(j)
            do k=1,ndlstrm
               n = iperm(inod(k,ne))
               if ( jwk(n).ne.i ) then
                  mi = m
                  do while ( ja(mi).gt.n .and. mi.ge.m0 )
                     ja(mi+1) = ja(mi)
                     mi = mi - 1
                  enddo
                  ja(mi+1) = n
                  m = m + 1
                  jwk(n) = i
               endif
            enddo
         enddo
c
         ia(i+1) = m+1
c
      enddo
c
      nelt = ia(ndgstrm+1) - 1
      write(iout,'(/A,i8)')
     &                   'Length of matrix for stream function is ',nelt
      if ( nelt.gt. mxnelt) then
        write(iout,'(A,i7,A,i7,A)')
     &    'current maxln=',mxnelt,' is too small, must > ',nelt,' STOP!'
        stop
      endif
c
      return
      end
c
c
c***********************************************************************
      subroutine matstr ( m,ishm,kaxi, nrdgeom,ndggeom,ndlgeom, nrdstrm,
     &                    ndgstrm,ndlstrm,nrdvelo,ndgvelo,ndlvelo,nelem,
     &                    inod,ngmx,x,y,u,v,nelt,a,ja,ia,arhs,iperm,jwk,
     &                    shg2,shg2x,shg2y,we2 )
c
c     this subroutine calculates the local stiffness matrix for the
c       stream function AND assmbles into the global matrix
c***********************************************************************
      implicit none
      integer m,ishm,nrdgeom,ndlgeom,nrdstrm,ndlstrm,nrdvelo,ndlvelo,
     &        kaxi,nelem,ndggeom,ndgvelo,ngmx,inod(ngmx,nelem), ndgstrm,
     &        nelt,ia(ndgstrm+1),ja(nelt),iperm(ndgstrm),jwk(ndgstrm)
     &        
      real*8 x(ndggeom),y(ndggeom),u(ndgvelo),v(ndgvelo),
     &       a(nelt),arhs(ndgstrm)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
c
c     local variables
      real*8 xn(9),yn(9),un(9),vn(9),at(9,9),hs(9),psix(9),psiy(9)
      real*8 wa, xx,dxxi,dxet,dyxi,dyet,ajac,alpha,dudyl,dvdxl,
     &       ajx,ajy,x0,xcoor,rajac,dxix,dxiy,detx,dety
      integer i,j,l,k, nk,ik,jk,irowst,ilast,n1,n2,n3
c
c     0. initialization of matrices
c     -----------------------------
      n1 = ishm*2 + nrdgeom
      n2 = ishm*2 + nrdstrm
      n3 = ishm*2 + nrdvelo
      do i=1,ndlstrm
        do j=1,ndlstrm
          at(i,j) = 0.d0
        enddo
      enddo
      do i=1,ndlstrm
        hs(i) = 0.d0
      enddo
c
c     1. coordinates of nodes
c     -----------------------
      x0 = x(inod(1,m))
      do i=1,ndlgeom
        k = inod(i,m)
        xn(i) = xcoor(x0,x(k))
        yn(i) = y(k)
      enddo
c
      do i=1,ndlvelo
        k = inod(i,m)
        un(i) = u(k)
        vn(i) = v(k)
      enddo
c
c     2. loop on integration points
c     -----------------------------
      do l=1,7+2*ishm
c
        xx   = 0.d0
        dxxi = 0.d0
        dxet = 0.d0
        dyxi = 0.d0
        dyet = 0.d0
        do i=1,ndlgeom
          xx   = xx   + xn(i)*shg2(i,l,n1)
          dxxi = dxxi + xn(i)*shg2x(i,l,n1)
          dxet = dxet + xn(i)*shg2y(i,l,n1)
          dyxi = dyxi + yn(i)*shg2x(i,l,n1)
          dyet = dyet + yn(i)*shg2y(i,l,n1)
        enddo
        ajac = dxxi * dyet - dxet * dyxi
        rajac= 1.d0/ajac
        dxix =  dyet*rajac
        dxiy = -dxet*rajac
        detx = -dyxi*rajac
        dety =  dxxi*rajac
        alpha = 1.d0 + kaxi*(xx-1.d0)
c
        dudyl = 0.d0
        dvdxl = 0.d0
        do i=1,ndlvelo
          dudyl = dudyl + un(i)*(shg2x(i,l,n3)*dxiy+shg2y(i,l,n3)*dety)
          dvdxl = dvdxl + vn(i)*(shg2x(i,l,n3)*dxix+shg2y(i,l,n3)*detx)
        enddo
c
c       3. stream function
c       ------------------
        do i=1,ndlstrm
          psix(i) = shg2x(i,l,n2)*dxix + shg2y(i,l,n2)*detx
          psiy(i) = shg2x(i,l,n2)*dxiy + shg2y(i,l,n2)*dety
        enddo
c
c       4. computes coefficients in the local matrix
c       --------------------------------------------
c        wa = we2(l,n2)*ajac
c     extra terms due to axisymmetry
        if(kaxi.eq.1)then
        wa=we2(l,n2)*ajac/xx
        do j=1,ndlstrm
          do i=1,ndlstrm
            at(i,j)=at(i,j)+shg2(i,l,n2)*shg2(j,l,n2)*wa
          enddo
        enddo
        endif
c
        wa = we2(l,n2)*ajac*alpha
        do j=1,ndlstrm
          ajx = wa*psix(j)
          ajy = wa*psiy(j)
          do i=1,ndlstrm
            at(i,j) = at(i,j) + ajx*psix(i) + ajy*psiy(i)
          enddo
        enddo
         

c
c        wa = wa * alpha * (dvdxl - dudyl)
        wa = wa *(dvdxl - dudyl)
        do i=1,ndlstrm
          hs(i) = hs(i) + wa*shg2(i,l,n2)
        enddo
c
      enddo
c
c     5. assembles the matrix and vector
c     ----------------------------------
      do i=1,ndlstrm
        nk = inod(i,m)
        ik = iperm(nk)
        arhs(ik) = arhs(ik) + hs(i)
      enddo
c
      do i=1,ndlstrm
        nk = inod(i,m)
        ik = iperm(nk)
c
        irowst = ia(ik)
        ilast  = ia(ik+1)-1
        do k=irowst,ilast
          jwk(ja(k)) = k
        enddo
c
        do j=1,ndlstrm
           nk = inod(j,m) 
           jk = iperm(nk)
           k = jwk(jk)
           a(k) = a(k) + at(i,j)
c           if ( k.eq.0 ) then
c             write(*,'(A,3i5)') 'error in assmbls',m,i,ik
c             stop
c           endif
        enddo
c
        do k=irowst,ilast
          jwk(ja(k)) = 0
        enddo
c
      enddo
c
      return
      end
c
c***********************************************************************
      subroutine matstr_axi(nrdstrm,ndgstrm,nrdgeom,ndggeom,x,y,
     &                      nbd,ibdnod,nic,ic,nbound,nside,
     &                      nelt,a,ia,ja,iperm,jwk,shg1,shg1x,we1)
c
c     This subroutine calcuated the extra terms to the matrix due to 
c     axisymmetry, and modify matrix a at output      
c
c     Note:
c        this subroutine assumes that nrdstrm<=nrdgeom<=2
c
c***********************************************************************
c     Jun 09, 2005 Pengtao Yue
c***********************************************************************
      implicit none
      integer nrdstrm,ndgstrm,nrdgeom,ndggeom,nbd,nic,nbound,nelt
      integer ibdnod(nbd),ic(nic+1),nside(nbound),ia(ndgstrm+1),
     &        ja(nelt),iperm(ndgstrm),jwk(ndgstrm)
      real(8) x(ndggeom),y(ndggeom),a(nelt),
     &        shg1(3,3,2),shg1x(3,3,2),we1(3,2)
c
      integer ndlstrm,ndlgeom,i,j,k,l,kk,n,nb,
     &        m,iside,ica1,ica,icb,nel,ib,ibm1,intg

c
      integer nk,ik,jk,irowst,ilast
      real*8 xn(3),yn(3), xcoor, at(3,3),yxil,wl
      integer iglob(3)
c
c     2. loop over each boundary section
c     ----------------------------------
      iside = 0
      ndlstrm=nrdstrm+1
      ndlgeom=nrdgeom+1
      intg=3
      do nb=1,nbound
c     loop over each closed bounary
        ica1 = ic(iside+1)
        do kk=1,nside(nb)
          iside = iside + 1
          ica = ic(iside)
          icb = ic(iside+1)
          nel = (icb-ica)/nrdgeom
c
c         loop over each boundary segment
c         ----------------------------------
          do m=1,nel
            ib = ica + nrdgeom*(m-1)
            ibm1 = ib - 1
c     get the nodes coordinates
            n = ibdnod(ib)
            xn(1) = x(n)
            yn(1) = y(n)
            do i=2,ndlgeom
              k = ibm1 + i
              if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
              n = ibdnod(k)
              xn(i) = xcoor(xn(1),x(n))
              yn(i) = y(n)
            enddo
c
c           loop over the integration points
c           --------------------------------
            at(1:ndlstrm,1:ndlstrm)=0.d0   
            do l=1,intg
               yxil=sum(yn(1:ndlgeom)*shg1x(1:ndlgeom,l,nrdgeom))
               wl=we1(l,nrdstrm)*yxil
               do i=1,ndlstrm
               do j=1,ndlstrm
                  at(i,j)=at(i,j)+shg1(i,l,nrdstrm)*shg1(j,l,nrdstrm)*wl
               enddo
               enddo
            enddo     
c     get the relation between local node index and global index
            if(nrdstrm.eq.1)then
               k=ibm1+ndlgeom
               if ( k.eq.icb .and. kk.eq.nside(nb) ) k=ica1
               iglob(1)=ibdnod(ib)
               iglob(2)=ibdnod(k)            
            else
               k=ibm1+ndlgeom
               if ( k.eq.icb .and. kk.eq.nside(nb) ) then
                  iglob(1:ndlgeom-1)=ibdnod(ib:k-1)
                  iglob(ndlgeom)=ibdnod(ica1)
               else
                  iglob(1:ndlgeom)=ibdnod(ib:k)
               endif
            endif
c
c           modify the global matrix
c           ------------------------
            do i=1,ndlstrm
               nk = iglob(i)
               ik = iperm(nk)
c
               irowst = ia(ik)
               ilast  = ia(ik+1)-1
               do k=irowst,ilast
                  jwk(ja(k)) = k
               enddo
c
               do j=1,ndlstrm
                  nk = iglob(j)
                  jk = iperm(nk)
                  k = jwk(jk)
                  a(k) = a(k) + at(i,j)
               enddo
c
c            do k=irowst,ilast
c               jwk(ja(k)) = 0
c            enddo
            enddo
          enddo
        enddo
      enddo
      return
      end
