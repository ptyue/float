c***********************************************************************
      subroutine INTPLT ( nrdgeom,nrdvelo,nrdpres,nrdelas,nrdphi,
     &                    ndggeom,ndgvelo,ndgpres,ndgelas,ndgphi,
     &                    maxnod,maxvrt,
     &                    nnode,nelem,inod,ngmx,x,y,
     &                    nbd,ibdnod,nbound,nic,ic,nside, nbelast,
     &                    nnodd,nvrtd,nelmd,inodd,xold,yold,necold,ncmx,
     &                    nbdd,ibddd,icold,
     &                    u,p,dudt,ncpelas,selas,desdt,ncpg,ncps,
     &                    ncpphi,phi,psi,dphidt,
     &                    iopt,lfil,maxit,eps, local,beta,
     &                    iwork,maxiwk,rwork,maxrwk,iout )
c
c     This routine projects flow field from an old mesh onto a new mesh
c
c       nrdgeom = order of interpolation for coordinates
c       ndggeom = number of nodes in the new mesh
c       nelem = number of elements in the new mesh
c       x(ndggeom),y(ndggeom) = coordinates in the new mesh
c       inod(ngmx,nelem) = element description table, new mesh
c       nbd = number of boundary nodes
c       ibdnod(nbd) = index for boundary nodes
c       nbound,ic,nic,nside,
c       ishape(nelem) = element shape, 0:triangular; 1=quadrilateral
c       nnodd = number of nodes in the old mesh
c       nelmd = number of elements in the old mesh
c       xold(nnodd), yold(nnodd) = coordinates in the old mesh
c       inodd(ngmx,nelmd) = element description table, old mesh
c       necold(ncmx,nelmd) = neighboring element number, old mesh
c       iopt,lfil,maxit,eps = parameters for the solver
c       beta = the parameter beta in the boundary integration
c       iout = write unit number
c       maxiwk,maxrwk = size of working arrays
c       nbelast = number of deformable elastic particles
c
c     OUTPUT
c       var(nf,ndggeom) = projected values for all variables
c
c     WORKING ARRAYS

c     INPUT:
c       nnode,nelem,nbd,x,y,inod,ibdnod,
c                     : new mesh information
c       nnodd,nelmd,nvrtd,inodd,necold,xold,yold,
c                     : old mesh infomation
c       nrdvelo,nrdpres,nrdelas = order of interpolations
c       ndgvelo,ndgpres,ndgelas = number of nodes for ...
c       u,p,dudt,selas,desdt = variables
c
c     OUTPUT:
c       u,p,dudt,phi,psi,dphidt,selas,desdt = interpolated variables
c
c     WORKING ARRAYS
c       iwork(maxiwk),rwork(maxrwk)
c     
c     Modified Feb 21, 2005 by Pengtao Yue
c        Phase-field variables added
c***********************************************************************
      implicit none
      integer nnode,nelem,ngmx,inod(ngmx,nelem),ncmx,
     &        nbd,ibdnod(nbd),nbound,nic,ic(nic+1),nside(nbound),
     &        nnodd,nvrtd,nelmd,inodd(ngmx,nelmd),necold(ncmx,nelmd),
     &        nbdd,ibddd(nbdd),icold(nic+1),
     &        nrdgeom,nrdvelo,nrdpres,nrdelas,nrdphi,
     &        ndggeom,ndgvelo,ndgpres,ndgelas,ndgphi,ncpelas,ncpphi,
     &        iout,maxmt,maxiwk,maxrwk, iwork(maxiwk),
     &        iopt,lfil,maxit, ncpg,ncps, nbelast,
     &        maxnod,maxvrt
      double precision  x(nnode),y(nnode),xold(nnodd),yold(nnodd),
     &                  rwork(maxrwk), eps, beta
      real(8) u(ncpg,maxnod),dudt(ncpg,maxnod),p(maxvrt),
     &        selas(ncps,maxvrt),desdt(ncps,maxvrt),phi(maxnod),
     &        dphidt(maxnod),psi(maxnod)
      logical local
c
      integer nrd(19),larhs,lvec,lrhs,lrwk,lvar,lshape,lia,
     &        liperm,ljperm,liwk,ljwk,liwork,lja,la,mxiwk,mxrwk
      integer nf,i,n,k,k1,k2,kk
c
C     Support for Petsc timing
c      include '../include/petsc.h'
c
c      call PLogEventBegin(PROJECTION_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     number of variables to be projected
c     -----------------------------------
      kk = 5 + 3*ncpphi
      nf = kk + 2*ncpelas
      if ( nbelast.gt.0 ) nf = nf + 3
c
c     working arrays
c     --------------
      lshape = 1
      liwk   = lshape + nelem
      ljwk   = liwk   + nnode+1
      liperm = ljwk   + nnode
      ljperm = liperm + nnode
      lia    = ljperm + nnode
      liwork = lia    + nnode+1
c
      do i=1,nelem
         iwork(lshape-1+i) = 0
      enddo
      call matlen1 ( maxmt, nelem,inod,ngmx,iwork(lshape), nrdgeom,
     &               iwork(liwk),iwork(ljwk),iwork(liwork) )
c
      lja    = lia    + nnode+1
      liwork = lja    + maxmt
c
      la     = 1
      lvar   = la    + maxmt
      larhs  = lvar  + nf*max0(nnode,nnodd)
      lvec   = larhs + nnode
      lrhs   = lvec  + nnode
      lrwk   = lrhs  + nf*nnode
c
c     loading variables
c     -----------------
      do i=1,4
         nrd(i) = nrdvelo
      enddo
      do i=1,nnodd
         n = lvar-1+(i-1)*nf
c      write(*,*)'nnnn',n,i,maxrwk,ndgvelo
         rwork(n+1) = dudt(1,i)
         rwork(n+2) = dudt(2,i)
         rwork(n+3) = u(1,i)
         rwork(n+4) = u(2,i)
      enddo
c
      nrd(5) = nrdpres
      if ( nrdpres.eq.0 ) then
         k = nelmd
      elseif ( nrdpres.eq.1 ) then
         k = nvrtd
      endif
      do i=1,k
         n = lvar-1+(i-1)*nf
         rwork(n+5) = p(i)
      enddo
c
      if(ncpphi.eq.1)then
         nrd(6:kk)=nrdphi
         if(nrdphi.eq.1)then
            k=nvrtd
         elseif(nrdphi.eq.2)then
            k=nnodd
         else
            write(IOUT,*)'Error in nrdphi(subroutine intplt)'
            stop
         endif
         do i=1,k
            n = lvar-1+(i-1)*nf
            rwork(n+6)=phi(i)
            rwork(n+7)=psi(i)
            rwork(n+8)=dphidt(i)
         enddo
      endif
c
      if(nrdelas.eq.1)then
      do k=1,ncpelas
         k1 = kk+k
         k2 = k1+ncpelas
         nrd(k1) = nrdelas
         nrd(k2) = nrdelas
         do i=1,nvrtd
            n = lvar-1+(i-1)*nf
            rwork(n+k1) = selas(k,i)
            rwork(n+k2) = desdt(k,i)
         enddo
      enddo
      elseif(nrdelas.eq.2)then
      do k=1,ncpelas
         k1 = kk+k
         k2 = k1+ncpelas
         nrd(k1) = nrdelas
         nrd(k2) = nrdelas
         do i=1,nnodd
            n = lvar-1+(i-1)*nf
            rwork(n+k1) = selas(k,i)
            rwork(n+k2) = desdt(k,i)
         enddo
      enddo
      else
         write(IOUT,*)'Error in nrdelas(subroutine intplt)'
         stop
      endif
c
      if ( nbelast.gt.0 ) then
      if(nrdelas.eq.1)then
        do k=1,3
           k1 = kk + 2*ncpelas + k
           nrd(k1) = nrdelas
           do i=1,nvrtd
              n = lvar-1+(i-1)*nf
              rwork(n+k1) = selas(k,i)
           enddo
        enddo
      elseif(nrdelas.eq.2)then
        do k=1,3
           k1 = kk + 2*ncpelas + k
           nrd(k1) = nrdelas
           do i=1,nnodd
              n = lvar-1+(i-1)*nf
              rwork(n+k1) = selas(k,i)
           enddo
        enddo
      endif
      endif
c
      mxiwk = maxiwk - liwork
      mxrwk = maxrwk - lrwk
c
c     least-square projection
c     -----------------------
      if ( .not.local ) then
        call maplsq ( rwork(lvar),nrd,nf, iout, nrdgeom,ndggeom,
     &                nelem,inod,ngmx,x,y,iwork(lshape),
     &                nbd,ibdnod,nbound,nic,ic,nside,
     &                nnodd,nelmd,inodd,xold,yold,necold,ncmx, 
     &                nbdd,ibddd,icold,
     &                iopt,lfil,maxit,eps, beta,
     &                rwork(la),iwork(lja),maxmt,iwork(lia),
     &                iwork(liperm),iwork(ljperm),rwork(larhs),
     &                rwork(lvec),rwork(lrhs),
     &                iwork(liwk),iwork(ljwk),
     &                iwork(liwork),mxiwk,rwork(lrwk),mxrwk )
c
c     point-wise local interpolation
c     ------------------------------
      else
        call maploc ( rwork(lvar),nrd,nf, nrdgeom,ndggeom,
     &                nelem,inod,ngmx,x,y,iwork(lshape),
     &                nbd,ibdnod,nbound,nic,ic,nside,
     &                nnodd,nelmd,inodd,xold,yold,necold,ncmx,
     &                nbdd,ibddd,icold,
     &                rwork(lrhs),iwork(liwk),iwork(liwork),mxiwk )
      endif
c
c     update interpolated variables
c     -----------------------------
      do i=1,ndgvelo
         n = lvar-1+(i-1)*nf
         dudt(1,i) = rwork(n+1)
         dudt(2,i) = rwork(n+2)
         u(1,i) = rwork(n+3)
         u(2,i) = rwork(n+4)
      enddo
c
      do i=1,ndgpres
         n = lvar-1+(i-1)*nf
         p(i) = rwork(n+5)
      enddo
c
      if(ncpphi.eq.1)then
         do i=1,ndgphi
         n = lvar-1+(i-1)*nf
               phi(i)=rwork(n+6)
               psi(i)=rwork(n+7)
            dphidt(i)=rwork(n+8)
         enddo
      endif
c
      do k=1,ncpelas
         k1 = kk+k
         k2 = k1+ncpelas
         do i=1,ndgelas
            n = lvar-1+(i-1)*nf
            selas(k,i) = rwork(n+k1)
            desdt(k,i) = rwork(n+k2)
         enddo
      enddo
c
      if ( nbelast.gt.0 ) then
        do k=1,3
           k1 = kk+2*ncpelas+k
           do i=1,ndgelas
              n = lvar-1+(i-1)*nf
              selas(k,i) = rwork(n+k1)
           enddo
        enddo
      endif
c
C     Support for Petsc timing
c      call PLogEventEnd(PROJECTION_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c
c***********************************************************************
      subroutine maplsq ( var,nrd,nf, iout, nrdgeom,ndggeom,
     &                    nelem,inod,ngmx,x,y,ishape,
     &                    nbd,ibdnod,nbound,nic,ic,nside,
     &                    nnodd,nelmd,inodd,xold,yold,necold,ncmx,
     &                    nbdd,ibddd,icold,
     &                    iopt,lfil,maxit,eps, beta,
     &                    a,ja,maxmt,ia,iperm,jperm,arhs,vec,rhs,
     &                    iwk,jwk,iwork,maxiwk,rwork,maxrwk )
c
c     this routine performs least-square mapping for the variables
c
c     INPUT
c       nf     = total number of variables
c       var(nf,nnodd) = all variables ( values change on return)
c       nrd(nf) = order of interpolation for each variable
c       nrdgeom = order of interpolation for coordinates
c       ndggeom = number of nodes in the new mesh
c       nelem = number of elements in the new mesh
c       x(ndggeom),y(ndggeom) = coordinates in the new mesh
c       inod(ngmx,nelem) = element description table, new mesh
c       nbd = number of boundary nodes
c       ibdnod(nbd) = index for boundary nodes
c       nbound,ic,nic,nside,
c       ishape(nelem) = element shape, 0:triangular; 1=quadrilateral
c       nnodd = number of nodes in the old mesh
c       nelmd = number of elements in the old mesh
c       xold(nnodd), yold(nnodd) = coordinates in the old mesh
c       inodd(ngmx,nelmd) = element description table, old mesh
c       necold(ncmx,nelmd) = neighboring element number, old mesh
c       iopt,lfil,maxit,eps = parameters for the solver
c       beta = the parameter beta in the boundary integration
c       iout = write unit number
c       maxmt = maximum length of matrix
c       maxiwk,maxrwk = size of working arrays
c
c     OUTPUT
c       var(nf,ndggeom) = projected values for all variables
c
c     WORKING ARRAYS
c       a(maxmt),ja(maxmt),ia(ndggeom+1) = matrix in CSR format
c       iperm(ndggeom) : iperm(global var#) = order of elimination
c       jperm(ndggeom) : jperm(order of elimination) = global var#
c       arhs(ndggeom) = right-hand-side vector
c       vec(ndggeom) = solution vector
c       iwk(ndggeom+1),jwk(ndggeom),iwork(maxiwk),rwork(maxrwk)
c       rhs(nf,ndggeom)
c
c     Howard Hu, February 18, 1996
c***********************************************************************
      implicit none
      integer nf,ngmx,nrd(nf),nelem,inod(ngmx,nelem),ishape(nelem),
     &        nrdgeom,ndggeom,nnodd,nelmd,inodd(ngmx,nelmd),
     &        ncmx,maxmt,maxiwk,maxrwk,necold(ncmx,nelmd),
     &        ia(ndggeom+1),ja(maxmt),iperm(ndggeom),jperm(ndggeom),
     &        nbd,ibdnod(nbd),nbound,nic,ic(nic+1),nside(nbound),
     &        nbdd,ibddd(nbdd),icold(nic+1),
     &        iout,iopt,lfil,maxit,
     &        iwk(ndggeom),jwk(ndggeom),
     &        iwork(maxiwk)
      real*8 var(nf,*), x(ndggeom),y(ndggeom),xold(nnodd),yold(nnodd),
     &       a(maxmt),arhs(ndggeom),vec(ndggeom),
     &       rwork(maxrwk),rhs(nf,ndggeom), eps, beta
      logical oldmsh
c
c     local variables
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
     &      ,shg1(3,3,2),shg1x(3,3,2),we1(3,2)
      integer i,m,ishm,nelt,neltv,
     &        ierr,j,n,ndglb, nod,nelm
      logical lsol
c  
C     Support for Petsc timing
c      include '../include/petsc.h'
c
      call shap1d ( shg1,shg1x,we1 )
      call shap2d ( shg2,shg2x,shg2y,we2 )
c
      ishm = ishape(1)
c
c     initialization
c     --------------
      do j=1,nf
        do i=1,ndggeom
           rhs(j,i) = 0.d0
        enddo
      enddo
c
c     2nd order interpolations
c     ------------------------
      lsol = .false.
      do j=1,nf
        if ( nrd(j).eq.2 ) lsol = .true.
      enddo
      if ( .not.lsol ) goto 100
c
C     Support for Petsc timing
c      call PLogEventBegin(MAT_PROJ_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      n = ndglb(2)
      call premap ( 2, n, nelem,ishape,inod,ngmx,ia,ja,
     &              maxmt,nelt, iperm,jperm,iwk,jwk,iwork,iout)
      do i=1,ndggeom
         jwk(i) = 0
      enddo
      do i=1,nelt
         a(i) = 0.d0
      enddo
c
c     construction of the stiffness matrix
c
      nod = 0
      nelm = 1
      do i=1,nelmd
         iwork(i) = 0
      enddo
      do m=1,nelem
         call matmap ( m,ishm, nrdgeom, ndggeom, nelem, inod,
     &                x,y, nf,var,nrd, xold,yold,inodd,necold,nnodd,
     &                nelmd, ngmx,ncmx, nelm,nod,iwork, nelt,a,ja,ia,
     &                rhs,iperm,iwork(nelmd+1),jwk,
     &                shg2,shg2x,shg2y,we2 )
      enddo
c
c     additional terms contributed by the boundary integral
c      if ( beta.ne.0.d0 )
c     &  call bdterm ( nrdgeom,ndggeom, nf,var,nrd,
c     &                x,y, xold,yold,nnodd,
c     &                nbd,ibdnod,nic,ic,nbound,nside,
c     &                nbdd,ibddd,icold,
c     &                a,rhs,ia,ja,iperm,nelt,jwk,
c     &                shg1,shg1x,we1,beta )
c
C     Support for Petsc timing
c      call PLogEventEnd(MAT_PROJ_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     solver 1
C     Support for Petsc timing
c      call PLogEventBegin(SOLV_PROJ_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      neltv = maxrwk-2*ndggeom-1
      oldmsh = .false.
      do j=1,nf
        if ( nrd(j).eq.2 ) then
          do i=1,n
             vec(i) = 0.d0
             arhs(iperm(i)) = rhs(j,i)
          enddo
          call solver ( vec, n,arhs,a,ja,ia,nelt,
     &                  iopt,lfil,maxit,eps,iout,oldmsh,ierr,
     &                  rwork,maxrwk,iwork,maxiwk,neltv )
          oldmsh = .true.
          do i=1,n
             var(j,i) = vec(iperm(i))
          enddo
        endif
      enddo
c
C     Support for Petsc timing
c      call PLogEventEnd(SOLV_PROJ_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c

100   continue
      lsol = .false.
      do j=1,nf
        if ( nrd(j).eq.1 ) lsol = .true.
      enddo
      if ( .not.lsol ) goto 200

c     first order interpolations
c     --------------------------     
C     Support for Petsc timing
c      call PLogEventBegin(MAT_PROJ_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      n = ndglb(1)
      call premap ( 1, n, nelem,ishape,inod,ngmx,ia,ja,
     &              maxmt,nelt, iperm,jperm,iwk,jwk,iwork,iout)
      do i=1,n
         jwk(i) = 0
      enddo
      do i=1,nelt
         a(i) = 0.d0
      enddo
c
      ishm = ishape(1)
      do m=1,nelem
         call matmap1 ( m,ishm,nrdgeom,ndggeom, 1,n,nelem,inod,
     &            ngmx,x,y, nelt,a,ja,ia,iperm,jwk,shg2,shg2x,shg2y,we2)
      enddo
c
c      if ( beta.ne.0.d0 ) 
c     &call bdterm1( nrdgeom,ndggeom, x,y,nbd,ibdnod,nic,ic,nbound,nside,
c     &              a,ia,ja,iperm,nelt,jwk, shg1,shg1x,we1, beta)
c
C     Support for Petsc timing
c      call PLogEventEnd(MAT_PROJ_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     solver 2
c     --------
C     Support for Petsc timing
c      call PLogEventBegin(SOLV_PROJ_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      oldmsh = .false.
      do j=1,nf
        if ( nrd(j).eq.1 ) then
          do i=1,n
             vec(i) = 0.d0
             arhs(iperm(i)) = rhs(j,i)
          enddo
          call solver ( vec, n,arhs,a,ja,ia,nelt,
     &                  iopt,lfil,maxit,eps,iout,oldmsh,ierr,
     &                  rwork,maxrwk,iwork,maxiwk,neltv )
          oldmsh = .true.
          do i=1,n
             var(j,i) = vec(iperm(i))
          enddo
        endif
      enddo
c
C     Support for Petsc timing
c      call PLogEventEnd(SOLV_PROJ_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c

200   continue
      lsol = .false.
      do j=1,nf
        if ( nrd(j).eq.0 ) lsol = .true.
      enddo
      if ( .not.lsol ) goto 300

c     zeroth order interpolations
c     ---------------------------
C     Support for Petsc timing
c      call PLogEventBegin(MAT_PROJ_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      n = ndglb(0)
      call premap ( 0, n, nelem,ishape,inod,ngmx,ia,ja,
     &              maxmt,nelt, iperm,jperm,iwk,jwk,iwork,iout)
      do i=1,n
         jwk(i) = 0
      enddo
      do i=1,nelt
         a(i) = 0.d0
      enddo
c
      ishm = ishape(1)
      do m=1,nelem
         call matmap1 ( m,ishm,nrdgeom,ndggeom, 0,n,nelem,inod,
     &            ngmx,x,y, nelt,a,ja,ia,iperm,jwk,shg2,shg2x,shg2y,we2)
      enddo
c
C     Support for Petsc timing
c      call PLogEventEnd(MAT_PROJ_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     solver 3
c     --------
C     Support for Petsc timing
c      call PLogEventBegin(SOLV_PROJ_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      oldmsh = .false.
      do j=1,nf
        if ( nrd(j).eq.0 ) then
          do i=1,n
             vec(i) = 0.d0
             arhs(iperm(i)) = rhs(j,i)
          enddo
          call solver ( vec, n,arhs,a,ja,ia,nelt,
     &                  iopt,lfil,maxit,eps,iout,oldmsh,ierr,
     &                  rwork,maxrwk,iwork,maxiwk,neltv )
          oldmsh = .true.
          do i=1,n
             var(j,i) = vec(iperm(i))
          enddo
        endif
      enddo
c
C     Support for Petsc timing
c      call PLogEventEnd(SOLV_PROJ_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
300   continue
      return
      end
c
c
c***********************************************************************
      subroutine premap ( nrd,ndg, nelem, ishape,inod,ngmx,ia,ja,
     &                    maxmt,nelt,iperm,jperm,iwk,jwk,iwork, iout)
c
c     this routine prepares for least-square projection
c***********************************************************************
      implicit none
      integer nrd,ndg,nelem,ngmx,maxmt,ishape(nelem),inod(ngmx,nelem),
     &        iperm(ndg),jperm(ndg),
     &        nelt,ia(ndg+1),ja(maxmt), iout,
     &        iwk(ndg+1),jwk(ndg),iwork(ngmx*nelem)
c
      integer ndl,ndloc,i,j,nd,ne,n,ishm,k,m,m0,mi, i1,is
c
      ishm = ishape(1)
      ndl  = ndloc(ishm,nrd)
c
c     calculate iwk(global var#) = last element# related
      do ne=1,nelem
         do j=1,ndl
            if ( nrd.eq.0 ) then
              i = ne
            else
              i = inod(j,ne)
            endif
            iwk(i) = ne
         enddo
      enddo
c
c     iperm(global var#) = order of elimination
c     jperm(order of elimination) = global var#
      n = 0
      do ne=1,nelem
         do j=1,ndl
            if ( nrd.eq.0 ) then
              i = ne
            else
              i = inod(j,ne)
            endif
            if ( iwk(i).eq.ne ) then
               n = n + 1
               iperm(i) = n
               jperm(n) = i
            endif
         enddo
      enddo
c
c     prepare storage of the matrix in the CSR format
c     -----------------------------------------------
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

      do i=1,ndg
        jwk(i) = 0
      enddo
      m = 0
      ia(1) = 1
      do i=1,ndg
         nd = jperm(i)
c
         m = m + 1
         ja(m) = i
         jwk(i) = i
         m0 = m+1
         do j=iwk(nd)+1,iwk(nd+1)
            ne = iwork(j)
            do k=1,ndl
               if ( nrd.eq.0 ) then
                 i1 = ne
               else
                 i1 = inod(k,ne)
               endif
               n = iperm(i1)
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
         ia(i+1) = m+1
c
      enddo
c
      nelt = ia(ndg+1) - 1
      write(iout,'(A,i8)') 
     &            'Length of matrix for least-square projection is',nelt
      if ( nelt.gt.maxmt ) then
        write(iout,'(A,i7,A,i7,A)')
     &     'current maxmt=',maxmt,' is too small, must > ',nelt,' STOP!'
        stop
      endif
c
      return
      end
c
c
c***********************************************************************
      subroutine matmap ( m,ishm,nrdgeom,ndggeom, nelem,inod,x,y,
     &                    nf,var,nrd,xold,yold,inodd,necold,nnodd,nelmd,
     &                    ngmx,ncmx, nelm,nod,nfl, nelt,a,ja,ia,rhs,
     &                    iperm,iwk,jwk, shg2,shg2x,shg2y,we2)
c
c     this routine calculates the element matrix and vectors for 
c       least-square projection AND assmbly
c
c       nf     = total number of variables
c       var(nf,nnodd) = all variables
c       nrd(nf) = order of interpolation for each variable
c***********************************************************************
      implicit none
      integer m,ishm,nf,nrdgeom,ndggeom,nelem,ngmx,ncmx,
     &        inod(ngmx,nelem),nelmd,nnodd,inodd(ngmx,nelmd),
     &        necold(ncmx,nelmd),nrd(nf),nelt,ia(ndggeom+1),ja(nelt),
     &        iperm(ndggeom),jwk(ndggeom),nod,nelm,nfl(nelmd),iwk(nelmd)
      real*8 x(ndggeom),y(ndggeom),a(nelt),rhs(nf,ndggeom),
     &       var(nf,nnodd), xold(nnodd),yold(nnodd)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
c
c     local variables
      real*8 xn(9),yn(9),at(9,9),ht(17,9),fl(17),
     &       we,wa, xx,yy,dxxi,dxet,dyxi,dyet,x0,xcoor,ajac
      integer i,j,l,ndl,ndloc,k,ierr, nk,ik,jk,irowst,ilast,n1,n2,n3
      integer ndlgeom,iout
c
c     initialization
c     --------------
      n1 = ishm*2 + nrdgeom
      n2 = ishm*2 + 2
      n3 = ishm*2 + 1
      ndl = ndloc(ishm,1)
      ndlgeom = ndloc(ishm,nrdgeom)
      do i=1,ndlgeom
         do j=1,ndlgeom
           at(i,j) = 0.d0
         enddo
      enddo
      do j=1,nf
         do i=1,ndlgeom
           ht(j,i) = 0.d0
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
c     loop on integration points
c     --------------------------
      do l=1,7+2*ishm
c
        xx = 0.d0
        yy = 0.d0
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
        ajac = dxxi*dyet - dxet*dyxi
c
        iout = -31
c        if ( m.eq.40934 ) iout = 31
        call INTPL ( xx,yy, fl,nf,var,nrd, xold,yold,inodd,necold,
     &               nnodd,nelmd,ngmx,ncmx,nelm,nod,nfl,iwk,ierr,iout)
	if ( ierr.ne.0 ) then
          write(*,'(A,i6/2(6(1pe12.5)/))') 'in element ',m,
     &                                    (xn(i),yn(i),i=1,6)
          stop
        endif
c
c       computes coefficients in the local matrix, 2nd order
c       ----------------------------------------------------
        we = we2(l,n2) * ajac
        do j=1,ndlgeom
           wa = we*shg2(j,l,n2)
           do i=1,ndlgeom
              at(i,j) = at(i,j) + wa*shg2(i,l,n2)
           enddo
        enddo
c
        do j=1,nf
           wa = we*fl(j)
           if ( nrd(j).eq.2 ) then
              do i=1,ndlgeom
                 ht(j,i) = ht(j,i) + wa*shg2(i,l,n2)
              enddo
           elseif ( nrd(j).eq.1 ) then
              do i=1,ndl
                 ht(j,i) = ht(j,i) + wa*shg2(i,l,n3)
              enddo
           elseif ( nrd(j).eq.0 ) then
              do i=1,1
                 ht(j,i) = ht(j,i) + wa*shg2(i,l,6)
              enddo
           endif
        enddo
c
      enddo
c
c     assembly
c     --------
      do j=1,nf
         if ( nrd(j).eq.2 ) then
            do i=1,ndlgeom
               nk = inod(i,m)
               rhs(j,nk) = rhs(j,nk) + ht(j,i)
            enddo
         elseif ( nrd(j).eq.1 ) then
            do i=1,ndl
               nk = inod(i,m)
               rhs(j,nk) = rhs(j,nk) + ht(j,i)
            enddo
         elseif ( nrd(j).eq.0 ) then
            do i=1,1
               nk = m
               rhs(j,nk) = rhs(j,nk) + ht(j,i)
            enddo
         endif
      enddo
c
      do i=1,ndlgeom
        nk = inod(i,m)
        ik = iperm(nk)
c
        irowst = ia(ik)
        ilast  = ia(ik+1)-1
        do k=irowst,ilast
          jwk(ja(k)) = k
        enddo
c
        do j=1,ndlgeom
          nk = inod(j,m)
          jk = iperm(nk)
          k = jwk(jk)
c          if ( k.eq.0 ) then
c            write(*,*) 'error in assmblp',m,i,ik
c            stop
c          endif
          a(k) = a(k) + at(i,j)
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
c
c***********************************************************************
      subroutine INTPL ( xp,yp,fl,nf,var,nrd, x,y,inod,nec,nnode,nelem,
     &                          ngmx,ncmx, nelm,nod,nfl,iwk,ierr,iout )
c
c     This routine interpolates the variables at a given node
c
c     INPUT
c       xp,yp  = coordinates for the given point
c       nf     = total number of variables
c       var(nf,nnode) = all variables
c       nrd(nf) = order of interpolation for each variable
c       nnode = number of nodes in the mesh
c       nelem = number of elements
c       x(nnode),y(nnode) = coordinates of the mesh nodes
c       inod(ngmx,nelem) = element description table
c       nec(ncmx,nelem) = neighboring element number
c       nelm = starting element number for the search
c
c     OUTPUT
c       fl(nf) = interpolated variables
c       nelm   = element number to which the point belongs
c
c     WORKING VARIABLES
c       nod       = index for the calling of the routine
c                  (takes different value for different calls)
c       nfl(nelem) = flag in the routine
c       iwk(nelem) = working array
c***********************************************************************
      implicit none
      integer nf,nelem,nnode, ngmx,ncmx,inod(ngmx,nelem),
     &        nec(ncmx,nelem),nrd(nf),nelm,nod,nfl(nelem),iwk(nelem)
      real*8 xp,yp, fl(nf), x(nnode),y(nnode),var(nf,nnode)
      integer iout
c
      integer ierr,n,i,j
      real*8 f(6),xi,et,t
c
c     1. locate the element (in old mesh) the point belongs
c     -----------------------------------------------------
      nod = nod + 1
      call LOCATC ( xi,et, xp,yp, nelm, x,y,inod,nec, nnode,nelem,
     &                          ngmx,ncmx,nod,nfl,iwk, ierr,iout )
      if (ierr.ne.0) then
         write(*,'(A,2(1pe15.7),i5)') 'error in LOCATC',xp,yp,ierr
         ierr = 10
         return	
      endif
c
      do j=1,nf
        fl(j) = 0.d0
      enddo
c
c     interpolations, 2nd order
c     -------------------------
      t = 1.d0-xi-et
      f(1) = 2.d0*t*t-t
      f(2) = 2.d0*xi*xi-xi
      f(3) = 2.d0*et*et-et
      f(4) = 4.d0*xi*t
      f(5) = 4.d0*xi*et
      f(6) = 4.d0*et*t
      do i=1,6
        n = inod(i,nelm)
        do j=1,nf
          if ( nrd(j).eq.2 ) fl(j) = fl(j) + var(j,n)*f(i)
        enddo
      enddo
c
c     linear interpolations
c     ---------------------
      f(1) = t
      f(2) = xi
      f(3) = et
      do i=1,3
        n = inod(i,nelm)
        do j=1,nf
          if ( nrd(j).eq.1 ) fl(j) = fl(j) + var(j,n)*f(i)
        enddo
      enddo
c
c     constant interpolations
c     -----------------------
      do j=1,nf
        if ( nrd(j).eq.0 ) fl(j) = fl(j) + var(j,nelm)
      enddo
c
      return
      end
c
c
c***********************************************************************
      subroutine LOCATC ( xi,et, xp,yp, nelm, x,y,inod,nec, nnode,nelem,
     &                               ngmx,ncmx,nod,nfl,neli, ierr,iout )
c
c    This routine locates the element and local coordinates for a given 
c         point
c
c    INPUT:
c      xp,yp     = coordinates of the point 
c      nnode = number of nodes in the mesh
c      nelem = number of elements
c      x(nnode),y(nnode) = coordinates of the mesh
c      inod(ngmx,nelem) = element description table
c      nec(ncmx,nelem)  = neighboring element number
c      nod       = index for the calling of the routine
c                  (takes different value for different calls)
c      nfl(nelem) = nod, (flag in the routine)
c
c    OUTPUT:
c      xi,et = local coordinates
c      nelm  = the element number to which the point belongs
c      ierr  = error index
c
c    WORKING ARRAYS:
c      neli(nelem) = working array
c    "fool proof" searching scheme, Mingyu Zhu June 30,1997
c***********************************************************************
      implicit none
      integer nelm,nnode,nelem, ngmx,ncmx,inod(ngmx,nelem)
     &        ,nec(ncmx,nelem), nod,nfl(nelem),neli(nelem),ierr
      real*8 xi,et,xp,yp, x(nnode),y(nnode)
      integer iout
c
      real*8 xn(9),yn(9), 
     &       xcoor,xpm,det,eps,eps1,eps2,
     &       tmp(3),eps0
      integer ns(3),noloc,nd,i,j,il,k
      integer ipbm,nmd
      real*8 tmp0,tmp1,xiod,etod
      data eps/-0.3/,eps1/-0.002/,eps2/-0.04/,eps0/-1.0/
c
      il = 0
      ipbm = 0
      tmp0 = 10.d0
      ierr = 0
      nmd = 0
      do while ( ipbm .le. 1 )
         nd = inod(1,nelm)
         xn(1) = x(nd)
         yn(1) = y(nd)
         do i=2,3
            nd = inod(i,nelm)
            xn(i) = xcoor(xn(1),x(nd)) - xn(1)
            yn(i) = y(nd) - yn(1)
         enddo
c 
         xpm = xcoor(xn(1),xp)
         det = xn(2)*yn(3)-xn(3)*yn(2)
         xi = ((xpm-xn(1))*yn(3)-xn(3)*(yp-yn(1)))/det
         et = (xn(2)*(yp-yn(1))-(xpm-xn(1))*yn(2))/det
c
c      if(iout.gt.0) write(iout,'(i8,5e13.4)') nelm,xi,et,1-xi-et,xp,yp
c
c     for more accurate local coordinates for this point
c
         if ((1.d0-xi-et.ge.eps .and. xi.ge.eps .and. et.ge.eps) .or.
     &       (ipbm.eq.1 .and. 1.d0-xi-et.ge.eps0 .and. xi.ge.eps0 
     &        .and. et.ge.eps0) ) then
            call RST2D(xi,et,xp,yp,nelm,inod,ngmx,x,y,iout)
c     
            if ( ipbm.ge.1 ) then
               tmp1 = max(-xi,0.d0) + max(-et,0.d0)
     &              + max(-1.d0+xi+et,0.d0)
               if (tmp1.lt.tmp0) then
                  xiod = xi
                  etod = et
                  nmd = nelm
                  tmp0 = tmp1
               endif
            endif
         endif
c
c     element found
c     -------------
         if (1.d0-xi-et.ge.eps1.and.xi.ge.eps1.and.et.ge.eps1) return
c
         tmp(1) = xi
         tmp(2) = et
         tmp(3) = 1.d0-xi-et
         j=0
         do i=1,3
            ns(i) = nec(i,nelm)
            if ( ns(i).gt.0 ) then
               if ( nfl(ns(i)).eq.nod ) ns(i) = 0   
            endif
            if ( ns(i).le.0 ) then
               j=j+1
               if (tmp(i).ge.eps2 .and. tmp(mod(i,3)+1).ge.eps1 .and.
     &              tmp(mod(i+1,3)+1).ge.eps1 ) return
            else
               noloc = ns(i)
            endif
         enddo
c
         nfl(nelm) = nod
         if (j.eq.3) then
            if (il.le.0) then
               if (ipbm.eq.0) then
                  ipbm = ipbm + 1
                  nod  = nod  + 1
               else
                  xi = xiod
                  et = etod
                  nelm = nmd
                  if ( nmd .eq. 0 ) then
                     write(*,*) 'LOCALC: please check intplt.f,',
     &                          ' you need to reduce eps'
                     stop
                  endif
                  write(*,*) '**********2nd turn**********'
                  write(*,'(A,2e16.6)') 'xp ,yp = ',xp,yp
                  write(*,'(A,i8,A,3i8)') 'nelm =',nelm,'  nec(1:3) =',
     &                                    (nec(i,nelm),i=1,3)
                  do k = 1,6
                     write(*,'(2e20.6)') x(inod(k,nelm)),y(inod(k,nelm))
                  enddo
                  write(*,'(3(A,f10.5))') 'xi =',xi,' et =',et,
     &                                    ' 1-xi-et =',1-xi-et
                  return
               endif
            else
               nelm=neli(il)
               il=il-1
            endif
         elseif (j.eq.2) then
            il=il+1
            neli(il)=nelm
            nelm=noloc
         else
            il=il+1
            neli(il)=nelm
            if((xi+et-1).ge.0 .and. ns(2).gt.0 ) then
               nelm = ns(2)
            elseif (xi.le.0 .and. ns(3).gt.0 ) then
               nelm = ns(3)
            elseif (et.le.0 .and. ns(1).gt.0 ) then
               nelm = ns(1)
            elseif ((xi+et-1).ge.0) then
               nelm = ns(3)
            elseif (xi.le.0) then
               nelm = ns(1)
            else 
               nelm = ns(2)
            endif
         endif
      enddo
c
c      if (iout.gt.0) write(iout,'(i7,4e13.5)') nelm,xi,et,xp,yp
c      if (iout.gt.0) write(iout,*) xp,yp,nelm,ntime
c
      return
      end
c
c**********************************************************************
      subroutine rst2d ( xi,et, xp,yp,nelm,inod,ngmx,x,y,iout )
c
c     compute the local coordinates in a curved 2d element
c**********************************************************************
      integer nelm,inod(ngmx,*),ngmx,iout
      real*8 xi,et,xp,yp, x(*),y(*)
c
      real*8 xn(9),yn(9), xcoor,psi(9),psix(9),psiy(9),
     &       det,f,g,fx,fy,gx,gy,dx,dy
      integer i,ndl,nd,iter,ndloc
c
      ndl = ndloc(0,2)
      do i=1,ndl
        nd = inod(i,nelm)
        xn(i) = xcoor(xp,x(nd))
        yn(i) = y(nd)
      enddo
c
      iter = 0
 25   call SHAPE2(ndl,xi,et,psi,psix,psiy,1)
      f = -xp
      g = -yp
      fx = 0.d0
      fy = 0.d0
      gx = 0.d0
      gy = 0.d0
      do i=1,ndl
         f = f + xn(i)*psi(i)
         g = g + yn(i)*psi(i)
         fx = fx + xn(i)*psix(i)
         fy = fy + xn(i)*psiy(i)
         gx = gx + yn(i)*psix(i)
         gy = gy + yn(i)*psiy(i)
      enddo
      det = 1.0/(fx*gy-fy*gx+1.d-90)
      dx = -(f*gy-g*fy)*det
      dy = -(g*fx-f*gx)*det
      xi = xi + dx
      et = et + dy
      iter = iter + 1
      if ( iter.gt.25 ) then
c         write(*,'(A,f8.3,A,f8.3,A/A,i5,A,2e11.4)') 
c     &        'RST2D warning: xi=',xi,' et=',et,' does not converge',
c     &        'in element',nelm,' for point ',xp,yp
         return
      endif
      if ( dabs(dx)+dabs(dy) .gt. 1.d-10 ) go to 25
c
c      if (iout.gt.0 ) write(iout,'(A10,i7,2e13.5)') 'newton',nelm,xi,et
      return
      end
c
c
c***********************************************************************
      subroutine matmap1 ( m,ishm, nrdgeom,ndggeom, nrd,n,
     &                     nelem,inod,ngmx,x,y,nelt,a,ja,ia,iperm,jwk,
     &                     shg2,shg2x,shg2y,we2)
c
c     this routine calculates the element matrix for least-square 
c       projection AND assembles the global matrix
c       for linear and piece-wise constant variables
c***********************************************************************
      implicit none
      integer m,ishm,nrdgeom,ndggeom,nelem,ngmx,inod(ngmx,nelem),
     &        nrd,n,nelt,ia(n+1),ja(nelt),iperm(n),jwk(n)
      real*8 x(ndggeom),y(ndggeom), a(nelt)
      real*8 shg2(9,9,6),shg2x(9,9,6),shg2y(9,9,6),we2(9,6)
c
c     local variables
      real*8 xn(9),yn(9),at(4,4),we,wa, dxxi,dxet,
     &       dyxi,dyet,x0,xcoor
      integer ndloc,ndl,ndlgeom,i,j,l,k, nk,ik,jk,irowst,ilast,n1,n2
c
c     initialization
c     --------------
      ndl = ndloc(ishm,nrd)
      ndlgeom = ndloc(ishm,nrdgeom)
      n1 = ishm*2 + nrdgeom
      if ( nrd.eq.0 ) then
        n2 = 6
      else
        n2 = ishm*2 + 1
      endif
      do i=1,ndl
        do j=1,ndl
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
      do l=1,7+2*ishm
c
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
c
c       computes the element matrix
c       ---------------------------
        we = we2(l,n2) * ( dxxi*dyet - dxet*dyxi )
        do j=1,ndl
          wa = we*shg2(j,l,n2)
          do i=1,ndl
            at(i,j) = at(i,j) + wa*shg2(i,l,n2)
          enddo
        enddo
c
      enddo
c
c     assembly
c     --------
      do i=1,ndl
        if ( nrd.eq.0 ) then
          nk = m
        else
          nk = inod(i,m)
        endif
        ik = iperm(nk)
c
        irowst = ia(ik)
        ilast  = ia(ik+1)-1
        do k=irowst,ilast
          jwk(ja(k)) = k
        enddo
c
        do j=1,ndl
          if ( nrd.eq.0 ) then
            nk = m
          else
            nk = inod(j,m)
          endif
           jk = iperm(nk)
           k = jwk(jk)
c           if ( k.eq.0 ) then
c             write(*,*) 'error in assmblp1',m,i,ik
c             stop
c           endif
           a(k) = a(k) + at(i,j)
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
c
c***********************************************************************
      subroutine bdterm( nrdgeom,ndggeom, nf,var,nrd,
     &                   x,y, xold,yold,nnodd,
     &                   nbd,ibdnod,nic,ic,nbound,nside,
     &                   nbdd,ibddd,icold,
     &                   a,rhs,ia,ja,iperm,nelt,jwk,
     &                   shg1,shg1x,we1, beta)
c
c     this routine calculates the boundary integral for projection.
c***********************************************************************
      implicit none
      integer nbd,ibdnod(nbd),nbound,nic,ic(nic+1),nside(nbound),
     &        nrdgeom,ndggeom,
     &        nelt,ia(ndggeom+1),ja(nelt),iperm(ndggeom),jwk(ndggeom),
     &        nf,nrd(nf), nnodd,nbdd,icold(nic+1),ibddd(nbdd)
      real*8 x(ndggeom),y(ndggeom), xold(nnodd),yold(nnodd),
     &       var(nf,nnodd), a(nelt),rhs(nf,ndggeom),
     &       shg1(3,3,2),shg1x(3,3,2),we1(3,2), beta
c
      integer nd,nb,n,m,iside,ica1,ica,icb,nel,ib,ibm1,ierr,
     &        nk,irowst,ilast, i,j,k,kk,l,ik,jk, icad1,m0
      real*8 xn(3),yn(3),at(3,3),ht(17,3), xcoor
      real*8 xx,yy,dxxi,dyxi,dsxi,wa,waj,fl(17)
      logical lcir
c
      nd = nrdgeom+1
c
c     1. loop over each boundary section
c     ----------------------------------
      iside = 0
      do 10 nb=1,nbound
        ica1 = ic(iside+1)
        icad1 = icold(iside+1)
        do 10 kk=1,nside(nb)
          iside = iside + 1
          ica = ic(iside)
          icb = ic(iside+1)
          nel = (icb-ica)/nrdgeom
c
          lcir = .false.
          if ( kk.eq.nside(nb) ) lcir = .true.
          m0 = 1
c
c         2. loop over each boundary segment
c         ----------------------------------
          do 20 m=1,nel
            ib = ica + nrdgeom*(m-1)
            ibm1 = ib - 1
c
            n = ibdnod(ib)
            xn(1) = x(n)
            yn(1) = y(n)
            do i=2,nd
              k = ibm1 + i
              if ( k.eq.icb .and. lcir ) k=ica1
              n = ibdnod(k)
              xn(i) = xcoor(xn(1),x(n))
              yn(i) = y(n)
            enddo
c
            do i=1,nd
              do j=1,nd
                at(i,j) = 0.d0
              enddo
            enddo
            do j=1,nf
              do i=1,nd
                ht(j,i) = 0.d0
              enddo
            enddo
c
c           3. loop over integration points
c           -------------------------------
            do l=1,3
c
c             calculating derivatives
              xx = 0.d0
              yy = 0.d0
              dxxi = 0.d0
              dyxi = 0.d0
              do i=1,nd
                xx = xx + xn(i)*shg1(i,l,nrdgeom)
                yy = yy + yn(i)*shg1(i,l,nrdgeom)
                dxxi = dxxi + xn(i)*shg1x(i,l,nrdgeom)
                dyxi = dyxi + yn(i)*shg1x(i,l,nrdgeom)
              enddo
              dsxi = dsqrt(dxxi*dxxi + dyxi*dyxi)
              wa = dsxi*we1(l,nrdgeom)
c
c             compute coefficients in the local matrix
              do j=1,nd
                waj = wa*shg1(j,l,2)
                do i=1,nd
                  at(i,j) = at(i,j) + waj*shg1(i,l,2)
                enddo
              enddo
c
c             local interploation
              call INTPL1D ( xx,yy, fl,nf,var,nrd, xold,yold,
     &                     ibddd,icold,iside,m0,icad1,lcir, ierr)
              if ( ierr.ne.0 ) then
                write(*,'(A,2e12.5/A,i3,A,i3/6(1pe12.5))') 'for point',
     &                  xx,yy,'on side ',iside,' segment ',m,
     &                  (xn(i),yn(i),i=1,nd)
                stop
              endif
c
c             3.4 calculating boundary integrals
              do j=1,nf
                waj = wa*fl(j)
                if ( nrd(j).eq.2 ) then
                  do i=1,3
                    ht(j,i) = ht(j,i) + waj*shg1(i,l,2)
                  enddo
                elseif ( nrd(j).eq.1 ) then
                  do i=1,2
                    ht(j,i) = ht(j,i) + waj*shg1(i,l,1)
                  enddo
                endif
              enddo
c
            enddo
c
c           4. assembly
c           -----------
            do j=1,nf
              if ( nrd(j).eq.2 ) then
                do i=1,3
                  k = ibm1 + i
                  if ( k.eq.icb .and. lcir ) k=ica1
                  nk = ibdnod(k)
                  rhs(j,nk) = rhs(j,nk) + beta*ht(j,i)
                enddo
              elseif ( nrd(j).eq.1 ) then
                do i=1,2
                  k = ibm1 + 2*i-1
                  if ( k.eq.icb .and. lcir ) k=ica1
                  nk = ibdnod(k)
                  rhs(j,nk) = rhs(j,nk) + beta*ht(j,i)
                enddo
              endif
            enddo
c
            do i=1,nd
              k = ibm1 + i
              if ( k.eq.icb .and. lcir ) k=ica1
              nk = ibdnod(k)
              ik = iperm(nk)
c
              irowst = ia(ik)
              ilast  = ia(ik+1)-1
              do k=irowst,ilast
                jwk(ja(k)) = k
              enddo
c
              do j=1,nd
                k = ibm1 + j
                if ( k.eq.icb .and. lcir ) k=ica1
                nk = ibdnod(k)
                jk = iperm(nk)
                k = jwk(jk)
c                if ( k.eq.0 ) then
c                  write(*,'(A,3i5)') 'error in assmblp',m,i,ik
c                  stop
c                endif
                a(k) = a(k) + beta*at(i,j)
              enddo
c
              do k=irowst,ilast
                jwk(ja(k)) = 0
              enddo
c
            enddo
c
 20       continue
c
10    continue
c
      return
      end
c
c***********************************************************************
      subroutine bdterm1( nrdgeom,ndggeom, x,y,
     &                    nbd,ibdnod,nic,ic,nbound,nside,
     &                    a,ia,ja,iperm,nelt,jwk,
     &                    shg1,shg1x,we1, beta )
c
c     this routine calculates the boundary integral for projection.
c***********************************************************************
      implicit none
      integer nbd,ibdnod(nbd),nbound,nic,ic(nic+1),nside(nbound),
     &        nrdgeom,ndggeom,
     &        nelt,ia(ndggeom+1),ja(nelt),iperm(ndggeom),jwk(ndggeom)
      real*8 x(ndggeom),y(ndggeom), a(nelt),
     &       shg1(3,3,2),shg1x(3,3,2),we1(3,2), beta
c
      integer nd,nb,n,m,iside,ica1,ica,icb,nel,ib,ibm1,
     &        nk,irowst,ilast, i,j,k,kk,l,ik,jk
      real*8 xn(3),yn(3),at(3,3), xcoor
      real*8 dxxi,dyxi,dsxi,wa,waj
      logical lcir
c
      nd = nrdgeom+1
c
c     1. loop over each boundary section
c     ----------------------------------
      iside = 0
      do 10 nb=1,nbound
        ica1 = ic(iside+1)
        do 10 kk=1,nside(nb)
          iside = iside + 1
          ica = ic(iside)
          icb = ic(iside+1)
          nel = (icb-ica)/nrdgeom
c
          lcir = .false.
          if ( kk.eq.nside(nb) ) lcir = .true.
c
c         2. loop over each boundary segment
c         ----------------------------------
          do 20 m=1,nel
            ib = ica + nrdgeom*(m-1)
            ibm1 = ib - 1
c
            n = ibdnod(ib)
            xn(1) = x(n)
            yn(1) = y(n)
            do i=2,nd
              k = ibm1 + i
              if ( k.eq.icb .and. lcir ) k=ica1
              n = ibdnod(k)
              xn(i) = xcoor(xn(1),x(n))
              yn(i) = y(n)
            enddo
c
            do i=1,2
              do j=1,2
                at(i,j) = 0.d0
              enddo
            enddo
c
c           3. loop over integration points
c           -------------------------------
            do l=1,3
c
c             3.1 calculating derivatives
              dxxi = 0.d0
              dyxi = 0.d0
              do i=1,nd
                dxxi = dxxi + xn(i)*shg1x(i,l,nrdgeom)
                dyxi = dyxi + yn(i)*shg1x(i,l,nrdgeom)
              enddo
              dsxi = dsqrt(dxxi*dxxi + dyxi*dyxi)
              wa = dsxi*we1(l,nrdgeom)
c
c             3.2 compute coefficients in the local matrix
              do j=1,2
                waj = wa*shg1(j,l,1)
                do i=1,2
                  at(i,j) = at(i,j) + waj*shg1(i,l,1)
                enddo
              enddo
c
            enddo
c
c           4. assembly
c           -----------
            do i=1,2
              k = ibm1 + 2*i-1
              if ( k.eq.icb .and. lcir ) k=ica1
              nk = ibdnod(k)
              ik = iperm(nk)
c
              irowst = ia(ik)
              ilast  = ia(ik+1)-1
              do k=irowst,ilast
                jwk(ja(k)) = k
              enddo
c
              do j=1,2
                k = ibm1 + 2*j-1
                if ( k.eq.icb .and. lcir ) k=ica1
                nk = ibdnod(k)
                jk = iperm(nk)
                k = jwk(jk)
c                if ( k.eq.0 ) then
c                  write(*,'(A,3i5)') 'error in assmblp',m,i,ik
c                  stop
c                endif
                a(k) = a(k) + beta*at(i,j)
              enddo
c
              do k=irowst,ilast
                jwk(ja(k)) = 0
              enddo
c
            enddo
c
 20       continue
c
10    continue
c
      return
      end
c
c***********************************************************************
      subroutine INTPL1D ( xx,yy, fl,nf,var,nrd, xold,yold,ibddd,icold,
     &                     iside,m0,ica1,lcir, ierr)
c
c     this routine interpolates functions along a curved line
c***********************************************************************
      implicit none
      integer nf,nrd(nf),ibddd(*),icold(*),iside,ica1,m0,ierr
      real*8 xx,yy,fl(nf),var(nf,*),xold(*),yold(*)
      logical lcir
c
      integer nrdgeom,ica,icb,nel,m,ib,n1,n2,n3,k,j
      real*8 x1,y1,x2,y2,x3,y3,c,ceta,xcoor,ds
c
      ierr = 0
      nrdgeom = 2
c
      ica = icold(iside)
      icb = icold(iside+1)
      nel = (icb-ica)/nrdgeom
c
21    do m=m0,nel
        ib = ica + nrdgeom*(m-1)
        n1 = ibddd(ib)
        x1 = xcoor(xx,xold(n1))
        y1 = yold(n1)
        k = ib + 2
        if ( k.eq.icb .and. lcir ) k=ica1
        n2 = ibddd(k)
        x2 = xcoor(xx,xold(n2))
        y2 = yold(n2)
        n3 = ibddd(ib+1)
        x3 = xcoor(xx,xold(n3))
        y3 = yold(n3)
c
        if ( dabs(x2-x1).gt.dabs(y2-y1) ) then
          ds = dabs(x2-x1)
          c = ceta(xx,x1,x2,x3)
        else
          ds = dabs(y2-y1)
          c = ceta(yy,y1,y2,y3)
        endif
c
        if ( dabs(xx-x3).le.ds .and. dabs(yy-y3).le.ds
     &       .and. (c.ge.0.d0).and.(c.le.1.d0) ) then
c
c          write(32,'(i3,7(1pe12.4))') m,xx,yy,x1,y1,x2,y2,c
c
          y1 = 1.d0-c
          y2 = c
          x1 = 1.d0-3.d0*c+2.d0*c*c
          x2 = -c+2.d0*c*c
          x3 = 4.d0*c-4.d0*c*c
          do j=1,nf
            if ( nrd(j).le.1 ) then
              fl(j) = var(j,n1)*y1 + var(j,n2)*y2
            elseif ( nrd(j).eq.2 ) then
              fl(j) = var(j,n1)*x1 + var(j,n2)*x2 + var(j,n3)*x3
            endif
          enddo
          m0 = m
          return
        endif
c
      enddo
c
      ierr = ierr+1 
      if (ierr.eq.1) then
        m0 = 1
        goto 21 
      endif
c
      return
      end
c
      function ceta(x,x1,x2,x3)
c
c     find the local coordinate ceta for a given global node x
c
      implicit none
      real*8 ceta,x,x1,x2,x3
      real*8 y1,y2,f1,f2,f3,dx,dc,ceta2
      integer iter
c
      ceta = (x-x1)/(x2-x1)
      iter = 0
      y1 = -3.d0*x1-x2+4.d0*x3
      y2 = 4.d0*(x1+x2-2.d0*x3)
c
10    iter = iter + 1
      ceta2 = ceta*ceta
      f1 = 1.d0-3.d0*ceta+2.d0*ceta2
      f2 = -ceta+2.d0*ceta2
      f3 = 4.d0*ceta-4.d0*ceta2
      dx = x - x1*f1 - x2*f2 - x3*f3
      dc = dx/(y1+y2*ceta)
      ceta = ceta + dc
      if ( dabs(dx).lt.1.d-10 ) return
      if ( iter.gt.15 ) then
        write(*,*) 'warning in INTPL1D, ceta does not converge'
        stop
      endif
      goto 10
      end
c
c
c***********************************************************************
      subroutine maploc ( var,nrd,nf, nrdgeom,ndggeom,
     &                    nelem,inod,ngmx,x,y,ishape,
     &                    nbd,ibdnod,nbound,nic,ic,nside,
     &                    nnodd,nelmd,inodd,xold,yold,necold,ncmx,
     &                    nbdd,ibddd,icold,
     &                    rhs,iwk,iwork,maxiwk )
c
c     the routine performs direct local interpolation for all variables
c
c     INPUT
c       nf     = total number of variables
c       var(nf,nnodd) = all variables ( values change on return)
c       nrd(nf) = order of interpolation for each variable
c       nrdgeom = order of interpolation for coordinates
c       ndggeom = number of nodes in the new mesh
c       nelem = number of elements in the new mesh
c       x(ndggeom),y(ndggeom) = coordinates in the new mesh
c       ngmx = maximum number of nodes in the element
c       inod(ngmx,nelem) = element description table, new mesh
c       nbd = number of boundary nodes
c       ibdnod(nbd) = index for boundary nodes
c       nbound,ic,nic,nside,
c       ishape(nelem) = element shape, 0:triangular; 1=quadrilateral
c       nnodd = number of nodes in the old mesh
c       nelmd = number of elements in the old mesh
c       xold(nnodd), yold(nnodd) = coordinates in the old mesh
c       inodd(ngmx,nelmd) = element description table, old mesh
c       ncmx = maximum number of sides in the element
c       necold(ncmx,nelmd) = neighboring element number, old mesh
c       maxiwk = size of working arrays
c
c     OUTPUT
c       var(nf,ndggeom) = projected values for all variables
c
c     WORKING ARRAYS
c       iwk(ndggeom),iwork(maxiwk), rhs(nf,ndggeom)
c
c     Howard Hu, March 17, 1997
c***********************************************************************
      implicit none
      integer nf,ngmx,nrd(nf),nelem,inod(ngmx,nelem),ishape(nelem),ncmx,
     &        nrdgeom,ndggeom,nnodd,nelmd,inodd(ngmx,nelmd),
     &        necold(ncmx,nelmd),
     &        nbd,ibdnod(nbd),nbound,nic,ic(nic+1),nside(nbound),
     &        nbdd,ibddd(nbdd),icold(nic+1),
     &        iwk(ndggeom),maxiwk,iwork(maxiwk)
      real*8 var(nf,*), x(ndggeom),y(ndggeom),xold(nnodd),yold(nnodd),
     &       rhs(nf,ndggeom)
c
c     local variables
      integer i,m,ishm,ndlgeom,ndloc,ierr,j,k, nod,nelm,iout
      real*8 fl(17),xx,yy
c  
c     ---------------------------------------------------------
      do i=1,ndggeom
        iwk(i) = 0
      enddo
      nod = 0
      nelm = 1
      do i=1,nelmd
        iwork(i) = 0
      enddo
c
      ishm = ishape(1)
      ndlgeom = ndloc(ishm,nrdgeom)
      do m=1,nelem
        do i=1,ndlgeom
          k = inod(i,m)
          if ( iwk(k) .eq. 0 ) then
c
          iwk(k) = 1
          xx = x(k)
          yy = y(k)
          iout = -31
c        if ( m.eq.573 ) iout = 31
          call INTPL ( xx,yy, fl,nf,var,nrd, xold,yold,inodd,necold,
     &                 nnodd,nelmd,ngmx,ncmx,nelm,nod,iwork,
     &                 iwork(nelmd+1),ierr,iout)
          if ( ierr.ne.0 ) then
            write(*,'(A,i6/2(1pe12.5))') 'in element ',m,xx,yy
            stop
          endif
c
          do j=1,nf
            rhs(j,k) = fl(j)
          enddo
c
          endif
c
        enddo
      enddo
c
      do j=1,nf
        do i=1,ndggeom
          var(j,i) = rhs(j,i)
        enddo
      enddo
c
      return
      end
