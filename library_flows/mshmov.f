c***********************************************************************
      subroutine mshmov ( oldmsh, nrdmsh,nrdgeom,ndgmsh,ndggeom,ndgvelo,
     &                    iout, uamesh,ncpg,ubdmsh,vbdmsh,
     &                    nelem,inod,ngmx,x,y,ishape, nic,ic,nbd,ibdnod,
     &                    shg2,shg2x,shg2y,we2,
     &                    iopt,lfil,maxit,epsm, smax,
     &                    a,ja,nelt,ia,iperm,jperm,imp,arhs,
     &                    vec,iwk,jwk,iwork,maxiwk,rwork,maxrwk )
c
c     this routine calculates the mesh velocity & mesh accelation
c
c     INPUT
c       oldmsh = F/T: compute the matrix / use an old one
c       nrdmsh = order of interpolation for mesh velocities
c       ndgmsh = number of nodes for mesh movement
c       ndggeom = number of nodes in mesh
c       ndgvelo = number of nodes for velocities
c       x(ndggeom),y(ndggeom) = coordinates in mesh
c       nelem = number of elements
c       inod(ngmx,nelem) = element description table
c       nic, ic(nic+1) = boundary section index
c       nbd = number of boundary nodes
c       ibdnod(nbd) = index for boundary nodes
c       ishape(nelem) = element shape, 0:triangular; 1=quadrilateral
c       nrdgeom = order of interpolation for coordinates
c       iopt    = option for the solver
c         lfil  = number of fill-ins in the ilut preconditioner
c         maxit = maximum number of iterations allowed in GMRES
c         epsm  = tolerance for stopping criterion
c         smax   = magnitudes of the variables
c       iout = write unit number 
c       nelt = length of the matrix
c       when oldmsh=.true., the following data are need; otherwise
c         they are computed in this routine       
c         a(nelt),ja(nelt),ia(ndgmsh+1) = matrix in CSR format
c         iperm(ndgmsh) : iperm(global var#) = order of elimination
c         jperm(ndgmsh) : jperm(order of elimination) = global var#
c         imp(ndgmsh) = index for boundary condition 
c
c     OUTPUT
c       uamesh(ncpg,ndgvelo) = mesh velocities or mesh accelation
c
c     WORKING VARIABLES
c       ubdmsh(nbd), vbdmsh(nbd) = mesh velocity or accelation on boundary
c       arhs(ndgmsh) = right hand side vector
c       vec(ndgmsh)  = solution vector
c       iwk(ndgvelo),jwk(ndgmsh),iwork(maxiwk),rwork(maxrwk) 
c
c     By Howard Hu & Mingyu Zhu, July 14, 1997
c***********************************************************************
      implicit none
      integer nelem,nbd,ngmx,ibdnod(nbd),inod(ngmx,nelem),
     &        ishape(nelem),nrdmsh,ndgmsh,nrdgeom,ndggeom,ndgvelo,
     &        iopt,lfil,maxit,nelt,maxiwk,maxrwk,iout,
     &        ia(ndgmsh+1),ja(nelt),imp(ndgmsh),iperm(ndgmsh),
     &        jperm(ndgmsh),iwk(ndgvelo),jwk(ndgmsh),iwork(maxiwk)
      integer ncpg,nic,ic(nic+1)
      double precision ubdmsh(nbd),vbdmsh(nbd),uamesh(ncpg,ndgvelo),
     &                 x(ndggeom),y(ndggeom),
     &                 a(nelt),arhs(ndgmsh),vec(ndgmsh),
     &                 rwork(maxrwk)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
      logical oldmsh
      double precision epsm,smax
c
      integer i,m,ishm,ndlgeom,ndlmsh,ndloc,nd,lenrw,leniw,neltv,
     &        ierr,ne,j,n1,n2, k,n
      double precision eps
      integer p3(3)
      data p3/2,3,1/
c
C     Support for Petsc timing
c      include '../include/petsc.h'
c
c      call PLogEventBegin(MESH_MOVE_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     1. preparations
c     ---------------
      if (.not.oldmsh ) then
         call premsh (nrdmsh,ndgmsh, nelem,inod,ishape,ngmx, nbd,ibdnod,
     &                ia,ja,iperm,jperm,imp, nelt,
     &                iwk,jwk,iwork,iout )
c
         do i=1,ndgmsh
            jwk(i) = 0
         enddo
         do i=1,nelt
            a(i) = 0.d0
         enddo
c     
c     2. construction of the stiffness matrix and triangulation
c     ---------------------------------------------------------
         ishm = ishape(1)
         ndlgeom = ndloc(ishm,nrdgeom)
         ndlmsh  = ndloc(ishm,nrdmsh)
         do m=1,nelem
            call matmsh (m,ishm,nrdgeom,ndlgeom,ndggeom,nrdmsh,ndlmsh,
     &                   ndgmsh,nelem,inod,ngmx,x,y,nelt,
     &                   a,ja,ia,iperm,imp,jwk,shg2,shg2x,shg2y,we2 )
         enddo
         do i=1,nbd
            nd = ibdnod(i)
            if ( nd.le.ndgmsh ) then
               j = iperm(nd)
               a(ia(j)) = 1.d0
            endif
         enddo
      endif
c
c     3. imposing boundary conditions for u-mesh
c     ------------------------------------------
      do i=1,ndgmsh
         vec(i) = uamesh(1,jperm(i))
         arhs(i) = 0.d0
      enddo
      do i=1,nbd
         nd = ibdnod(i)
         if ( nd.le.ndgmsh ) then
            j = iperm(nd)
            vec(j) = ubdmsh(i)
            arhs(j) = ubdmsh(i)
         endif
      enddo
c
c     4. solver for u-mesh
c     --------------------
      eps = epsm*smax
      neltv = maxrwk-2*ndgmsh-1
      call solver (vec, ndgmsh,arhs,a,ja,ia,nelt,
     &             iopt,lfil,maxit,eps,iout,oldmsh,ierr,
     &             rwork,maxrwk,iwork,maxiwk,neltv )
      if (ierr.gt.1 ) then
         write(*,*) 'error in solver, ierr=', ierr
         stop 
      endif
      oldmsh = .true.
c
      do i=1,ndgmsh
         uamesh(1,i) = vec(iperm(i))
      enddo
      if ( ndgmsh.lt.ndgvelo ) then
         do i=1,ndgvelo
            iwk(i) = 0
         enddo
         do ne=1,nelem
            do j=1,3
               nd = inod(3+j,ne)
               if ( iwk(nd).eq.0 ) then
                  iwk(nd) = 1
                  n1 = inod(j,ne)
                  n2 = inod(p3(j),ne)
                  uamesh(1,nd) = 0.5d0*(uamesh(1,n1)+uamesh(1,n2))
               endif
            enddo
         enddo
      endif
c
      do i=1,nbd
         nd = ibdnod(i)
         uamesh(1,nd) = ubdmsh(i)
      enddo
c
c     5. imposing boundary conditions v-mesh
c     --------------------------------------
      do i=1,ndgmsh
         vec(i) = uamesh(2,jperm(i))
         arhs(i) = 0.d0
      enddo
      do i=1,nbd
         nd = ibdnod(i)
         if ( nd.le.ndgmsh ) then
            j = iperm(nd)
            vec(j) = vbdmsh(i)
            arhs(j) = vbdmsh(i)
         endif
      enddo
c
c     6. solver for v-mesh
c     --------------------
      call solver (vec, ndgmsh,arhs,a,ja,ia,nelt,
     &             iopt,lfil,maxit,eps,iout,oldmsh,ierr,
     &             rwork,maxrwk,iwork,maxiwk,neltv )
      if (ierr.gt.1 ) then
         write(*,*) 'error in solver, ierr=', ierr
         stop
      endif
c     
      do i=1,ndgmsh
         uamesh(2,i) = vec(iperm(i))
      enddo
      if ( ndgmsh.lt.ndgvelo ) then
         do i=1,ndgvelo
            iwk(i) = 0
         enddo
         do ne=1,nelem
            do j=1,3
               nd = inod(3+j,ne)
               if ( iwk(nd).eq.0 ) then
                  iwk(nd) = 1
                  n1 = inod(j,ne)
                  n2 = inod(p3(j),ne)
                  uamesh(2,nd) = 0.5d0*(uamesh(2,n1)+uamesh(2,n2))
               endif
            enddo
         enddo
      endif
c
      do i=1,nbd
         nd = ibdnod(i)
         uamesh(2,nd) = vbdmsh(i)
      enddo
c
C     Support for Petsc timing
c      call PLogEventEnd(MESH_MOVE_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c
c***********************************************************************
      subroutine premsh ( nrdmsh,ndgmsh, nelem,inod,ishape,ngmx,
     &                    nbd,ibdnod,ia,ja,iperm,jperm,imp, nelt,
     &                    iwk,jwk,iwork,iout )
c
c     this routine prepares for mesh velocity and accelation
c***********************************************************************
      implicit none
      integer nbd,nelem,ngmx,nrdmsh,ndgmsh,ishape(nelem),
     &        inod(ngmx,nelem),ibdnod(nbd), ia(ndgmsh+1),iperm(ndgmsh),
     &        jperm(ndgmsh),imp(ndgmsh),nelt,iwk(ndgmsh+1),jwk(ndgmsh),
     &        ja(nelt),iwork(3*nelem), iout
c
      integer ndlmsh,ndloc,i,j,nd,ne,n,ishm,k,m,m0,mi,is
c
c     search of boundary conditions: imp
c     ----------------------------------
      do i=1,ndgmsh
        imp(i) = 0
      enddo
      do i=1,nbd
         nd = ibdnod(i)
         if ( nd.le.ndgmsh ) imp(nd) = 1
      enddo
c
c     computation of the order of elimination ( permutation idex)
c     -----------------------------------------------------------
c     calculate iwk(global var#) = last element# related
      ishm = ishape(1)
      ndlmsh = ndloc(ishm,nrdmsh)
      do ne=1,nelem
         do j=1,ndlmsh
           i = inod(j,ne)
           iwk(i) = ne
         enddo
      enddo
c     iperm(global var#) = order of elimination
c     jperm(order of elimination) = global var#
      n = 0
      do ne=1,nelem
         do j=1,ndlmsh
            i = inod(j,ne)
            if ( iwk(i).eq.ne ) then
               n = n + 1
               iperm(i) = n
               jperm(n) = i
            endif
         enddo
      enddo
c
c     preparation of storage of the matrix in the CSR format
c     ------------------------------------------------------
      do i=1,ndgmsh+1
         iwk(i) = 0
      enddo
      do i=1,nelem
         do j=1,ndlmsh
            is = inod(j,i)+1
            iwk(is) = iwk(is) + 1
         enddo
      enddo
      do i=2,ndgmsh+1
         iwk(i)= iwk(i-1) + iwk(i)
      enddo

      do i=1,nelem
         do j=1,ndlmsh
            is = inod(j,i)
            iwk(is) = iwk(is) + 1
            iwork(iwk(is)) = i
         enddo
      enddo
      do i=ndgmsh+1,2,-1
         iwk(i) = iwk(i-1)
      enddo
      iwk(1) = 0
c      write(25,'(10i6)') (iwk(i),i=1,ndgmsh+1)
c      write(25,'(12i6)') (iwork(i),i=1,3*nelem)
c
      do i=1,ndgmsh
         jwk(i) = 0
      enddo
      m = 0
      ia(1) = 1
      do i=1,ndgmsh
         nd = jperm(i)
c
         m = m + 1
         ja(m) = i
         if ( imp(nd).eq.0 ) then
c
           jwk(i) = i
           m0 = m+1
           do j=iwk(nd)+1,iwk(nd+1)
             ne = iwork(j)
             do k=1,ndlmsh
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
         endif

         ia(i+1) = m+1
c
      enddo
c
      m = ia(ndgmsh+1) - 1
      if ( m.gt.nelt ) then
         write(iout,'(A,i7,A,i7,A)')
     &        'MSHMOV: nelt=',nelt,' is too small, must > ',m,' STOP!'
         stop
      endif
      nelt = m
c
      return
      end
c
c
c***********************************************************************
      subroutine matmsh (m,ishm, nrdgeom,ndlgeom,ndggeom,nrdmsh,ndlmsh,
     &                   ndgmsh,nelem,inod,ngmx,x,y, nelt,a,ja,ia,
     &                   iperm,imp,jwk,shg2,shg2x,shg2y,we2 )
c
c     this subroutine calculates the element matrix for mesh velocity 
c     and mesh accelation, and assembles the matrix and vector
c***********************************************************************
      implicit none
      integer m,ishm,nrdgeom,ndlgeom,ndggeom,nelem, ngmx,nrdmsh,
     &        ndlmsh,ndgmsh,nelt,inod(ngmx,nelem),ia(ndgmsh+1),ja(nelt),
     &        iperm(ndgmsh),jwk(ndgmsh),imp(ndgmsh)
      double precision x(ndggeom),y(ndggeom),a(nelt)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
c
c     local variables
      double precision psix(9),psiy(9), xn(9),yn(9),at(9,9)
      double precision wa,dxxi,dxet,dyxi,dyet,ajac,ajx,ajy,x0,
     &                 xcoor,rajac,dxix,dxiy,detx,dety
      integer i,j,k,l, nk,ik,jk,irowst,ilast, n1,n2
c
c     initialization of matrices
c     --------------------------
      n1 = ishm*2 + nrdgeom
      n2 = ishm*2 + nrdmsh
      do i=1,ngmx
         do j=1,ngmx
            at(i,j) = 0.d0
         enddo
      enddo
c
c     coordinates of nodes
c     --------------------
      x0 = x(inod(1,m))
      do i=1,ndlgeom
         k = inod(i,m)
         xn(i) = xcoor(x0,x(k))
         yn(i) = y(k)
      enddo
c
c     ----------------------------
c     loop over integration points
c     ----------------------------
      do l=1,7+2*ishm
c
c        Jacobian
c        --------
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
        ajac  = dxxi*dyet - dxet*dyxi
        rajac = 1.0/ajac
        dxix =  dyet*rajac
        dxiy = -dxet*rajac
        detx = -dyxi*rajac
        dety =  dxxi*rajac
c
c       mesh velcoity
c       -------------
        do i=1,ndlmsh
           psix(i) = shg2x(i,l,n2)*dxix + shg2y(i,l,n2)*detx
           psiy(i) = shg2x(i,l,n2)*dxiy + shg2y(i,l,n2)*dety
        enddo
c
c       computes coefficients in the local matrix
c       -----------------------------------------
c        wa = we2(l,n2)*ajac
        wa = we2(l,n2)
        do j=1,ndlmsh
           ajx = wa*psix(j)
           ajy = wa*psiy(j)
           do i=1,ndlmsh
              at(i,j) = at(i,j) + ajx*psix(i) + ajy*psiy(i)
           enddo
        enddo
c
      enddo
c
c     assembles the matrix
c     --------------------
      do i=1,ndlmsh
         nk = inod(i,m)
         ik = iperm(nk)
c
         if ( imp(nk).eq.0 ) then
            irowst = ia(ik)
            ilast  = ia(ik+1)-1
            do k=irowst,ilast
               jwk(ja(k)) = k
            enddo
c
            do j=1,ndlmsh
               nk = inod(j,m)
               jk = iperm(nk)
               k = jwk(jk)
               a(k) = a(k) + at(i,j)
c              if ( k.eq.0 ) then
c                write(*,'(A,3i5)') 'error in assmbly mshvel',m,i,ik
c                stop
c              endif
            enddo
c
            do k=irowst,ilast
               jwk(ja(k)) = 0
            enddo
         endif
c
      enddo
c
      return
      end
