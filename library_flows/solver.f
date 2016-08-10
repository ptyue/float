c***********************************************************************
      subroutine solveru (nA,a,ja,ia, nB,b,jb,ib,rhs,sol,nelt,neltB,n,
     &                    iopt,maxits,iout,eps,im,lfil,old,ierr,lfilter,
     &                    rwork,lenrw,iwork,leniw,neltv,res,globits,
     &                    ioptp,maxitp,epsp,smaxp,sx,isx,jsx,neltS )
c
c     solver for the fluid-partcile equations.
c     We use block preconditioners for the system
c
c       saddle-point problem A & B:  | A  B^t | | u |   | f |
c                                    | B  0   | | p | = | g |
c                  for global matrix system B=0
c
c     INPUT:
c       nA,a(nelt),ja(nelt),ia(nA+1),nelt = matrix A in CSR format
c       nB,b(neltB),jb(neltB),ib(neqB+1),neltB = B matrix in CSR
c       n = nA+nB
c       rhs(n) = right-hand-side vector of the equation
c       iopt = xyz = option for the solver
c              x = 0    -- global solver (B=0)
c              x = 1    -- multilever solver for saddle point problem
c              x = 2    -- block LU for saddle point problem
c              x = 3    -- BfB^t preconditioned solver
c              x = 4    -- B(F^-1)B^t preconditioned solver
c              yz = 11  -- GMRES/diagonal preconditioning
c              yz = 12  -- GMRES/ilu0 preconditioning
c              yz = 13  -- GMRES/ilut preconditioning
c              yz = 21  -- BCGStab/diagonal preconditioning
c              yz = 22  -- BCGStab/ilu0 preconditioning
c              yz = 23  -- BCGStab/ilut preconditioning
c       maxits = maximum number of iterations allowed in the solver
c       im   = size of the krylov subspace for GMRES
c       lfil = number of fill-ins in the ilut preconditioner
c       eps  = tolerance for stopping criterion - aboslute
c       iout = write unit number
c       old  = .false. -> new matrix, preconditioning is needed
c            = .true. -> old matrix, use existing preconditioner
c       rwork(lenrw), iwork(leniw) = working arrays
c            check below for their dimensions
c
c     OUTPUT:
c       sol(n)  = solution vector
c       neltv   = length of the matrix for direct solver
c       res     = the residual
c       globits = global iteration count for nonlinear solve
c       ierr    = error index
c
c     Note:
c       if filter is called, matrix a,ia,ja is modified on return
c
      implicit none
      integer nA,nB,im,lfil,maxits,nelt,iout,ierr,n,neltB
      integer ja(nelt),ia(nA+1),jb(neltB),ib(nB+1)
      integer iopt,lenrw,leniw,neltv,iwork(leniw),globits
      real*8  rhs(n),sol(n),a(nelt),b(neltB),eps,rwork(lenrw),res
      logical old, lfilter
c
      integer ioptp,maxitp,neltS, isx(nB+1),jsx(neltS)
      real*8 epsp,smaxp, sx(neltS)
c
c     local variables
      integer lenalu,lenslu,lenju,lensju,lenw,lenwss 
      integer localu,locw,locrw,locjlu,locju
      integer locs,locy,locv,locb, locslu,locjslu,locjsu
      integer locjw,lociw, i,k1,k2,k, locr,loci,locz,locp, locws
      double precision  tol,tnorm,drptol, tmp
      integer iopt1,iopt2
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(SOLVER_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      ierr = 0
c
      iopt1 = mod(iopt,10)
      iopt2 = mod(iopt/10,10)
c
      if ( iopt1.eq.1 ) then
         lenalu = nA
         lenju  = 0
         lenw   = 0
      elseif ( iopt1.eq.2 ) then
         lenalu = nelt + 1
         lenju  = nA
         lenw   = nA
      elseif ( iopt1.eq.3 ) then
         lenalu = nelt + 2*lfil*nA + 1
         lenju  = nA
         lenw   = 2*nA
      endif

      lenslu = 0
      lensju = 0
      if ( iopt/100.ge.2 ) then
         if ( ioptp/10.eq.1 ) then
            lenslu = nB
         elseif ( ioptp/10.ge.2 ) then
            lenslu = neltS + 1
            lensju = nB
            lenw   = max(lenw,nB)
         endif
      endif
c
c     storage management - integer
c     ----------------------------
      locjlu   = 1
      locju    = locjlu + lenalu
      locjslu  = locju  + lenju
      locjsu   = locjslu+ lenslu
      locjw    = locjsu + lensju
      lociw    = locjw +  lenw

      if ( lociw-1.gt.leniw ) then
         write(iout,*) 'solveru: iwork too short, current size=',
     &                 leniw,lociw-1,' is needed'
         if ( lociw-1.gt.leniw ) ierr =10
      endif
c
c     storage management - real
c     -------------------------
c     rwork for S-solver
      lenwss = nA
      if ( iopt/100.ge.2 ) then
         k = mod(ioptp,10)
         if ( k.eq.2 .or. k.eq.3 ) then
            lenwss = lenwss + 3*nB + nA
         elseif ( k.eq.5 .or. k.eq.6 ) then
            lenwss = lenwss + 5*nB + nA
         endif
      endif

      localu   = 1
      locslu   = localu + lenalu
      locws    = locslu + lenslu
      locw     = locws  + lenwss
c
      if ( iopt2.eq.1 ) then
        locrw  = locw   + n*(im+2)

      elseif ( iopt2.eq.2) then
        locz   = locw   + n
        locp   = locz   + n
        locs   = locp   + n
        locy   = locs   + n
        locv   = locy   + n
        locb   = locv   + n
        locrw  = locb   + n

      endif

      if ( locrw-1.ne.lenrw ) then
         write(iout,*) 'solveru: rwork too short, current size=',
     &                 lenrw,locrw-1,' is needed'
         if ( locrw-1.gt.lenrw ) ierr =10
      endif
      if( ierr.ne.0 ) return
c
      if ( old ) goto 10
c
c     preconditioning for A: diagonal, ilu0 & ilut
c     --------------------------------------------
      if ( iopt1.eq.1 ) then
         do i=1,nA
            tmp = a(ia(i))
            if ( tmp.eq.0.d0 ) then
               do k=ia(i),ia(i+1)-1
                  tmp = tmp + dabs(a(k))
               enddo
            endif
            rwork(localu-1+i) = 1.d0/tmp
         enddo
c
      elseif ( iopt1.eq.2) then
         call ilu0 (nA,a,ja,ia, rwork(localu), iwork(locjlu), 
     &              iwork(locju), iwork(locjw), nelt,ierr)
c
      elseif ( iopt1.eq.3 ) then
        tol = 1.d-3
        call ilut (nA,a,ja,ia, lfil,tol,rwork(localu),iwork(locjlu),
     &             iwork(locju),lenalu,rwork(locw),iwork(locjw),ierr)
c
      endif
      if ( ierr.ne.0 ) then
         write(iout,*) 'solver: preconditioning error A, ierr=',ierr
         stop
      endif
c
c     filter out the zeros in the matrix
c     ----------------------------------
      if ( lfilter .and. iopt2.ne.0 ) then
         drptol = 1.d-12
         lenalu = nelt
         call filter (nA, drptol, a,ja,ia, lenalu, ierr)
         write(iout,*) 'after zero-out filter length of A matrix is',
     &                 lenalu
c         write(*,*) 'after zero-out filter length of A matrix is',
c     &                 lenalu, nelt, real(lenalu)/nelt
      endif
c
c     for approximate Schur complement S
c     ----------------------------------
      if ( iopt/100.ge.2 ) then
        call FormS (nA,a,ja,ia, nB,b,jb,ib,sx,jsx,isx,ioptp,rwork(locw))
c
        if ( ioptp/10.eq.1 ) then
           do i=1,nB
              rwork(locslu+i-1) = 1.d0/sx(isx(i))
           enddo

        elseif ( ioptp/10.ge.2 ) then
           call ilu0 (nB,sx,jsx,isx, rwork(locslu), iwork(locjslu),
     &                iwork(locjsu),iwork(locjw),neltS,ierr)

        endif
        if ( ierr.ne.0 ) then
           write(iout,*) 'solver: preconditioning error S, ierr=',ierr
           stop
        endif
      endif
c
c     solvers
c     -------
10    continue
c
      if ( iopt2.eq.1 ) then
        call gmresu (sol, n,rhs,nA,a,ja,ia,nB,b,jb,ib,rwork(localu),
     &               iwork(locjlu),iwork(locju), im,rwork(locw),
     &               eps, maxits, iout,ierr,globits, res,iopt,
     &               rwork(locws),lenwss,ioptp,maxitp,epsp,smaxp,
     &               sx,isx,jsx,
     &               rwork(locslu),iwork(locjslu),iwork(locjsu))

      elseif ( iopt2.eq.2 ) then
        call bcgstabu(sol, n,rhs,nA,a,ja,ia,nB,b,jb,ib,rwork(localu),
     &                iwork(locjlu),iwork(locju),eps,maxits,iout,
     &                ierr,res,iopt,rwork(locw),rwork(locz),rwork(locp),
     &                rwork(locb),rwork(locs),rwork(locy),rwork(locv),
     &                rwork(locws),lenwss,ioptp,maxitp,epsp,smaxp,
     &                sx,isx,jsx,
     &                rwork(locslu),iwork(locjslu),iwork(locjsu))

      endif
c
cC     Support for Petsc timing
c      call PLogEventEnd(SOLVER_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c***********************************************************************
      subroutine gmresu (sol, n,rhs,nA,a,ja,ia, nB,b,jb,ib, alu,jlu,ju,
     &                   im,vv,eps,maxits,iout,ierr,globits,res,iopt,
     &                   rwork,maxrwk, ioptp,maxitp,epsp,smaxp,
     &                   sx,isx,jsx, slu,jslu,jsu)
      implicit none
      integer n, im, maxits, iout, ierr, jlu(*),ju(*)
      integer nA,ja(*),ia(nA+1), nB,jb(*),ib(nB+1)
      integer globits,iopt
      real*8 vv(n,*), rhs(n), sol(n), a(*), b(*), alu(*), eps,res
      integer isx(*),jsx(*), jsu(*),jslu(*)
      real*8 sx(*),slu(*)
c----------------------------------------------------------------------*
c                                                                      *
c                 ***  Right-Preconditioned GMRES ***                  *
c                                                                      *
c----------------------------------------------------------------------*
c This is a simple version a preconditioned GMRES algorithm.           *
c The stopping criterion utilized is based simply on reducing the      * 
c residual norm by epsilon.  We recommend using a nonzero tol          *
c (tol=.005 or .001 usually give good results) in ILUT.                *
c                                                                      *
c on entry:                                                            *
c n      = integer. The dimension of the system                        *
c im     = size of krylov subspace:  should not exceed 50 here         *
c rhs    = real vector of length n containing the right hand side.     *
c sol    = real vector of length n containing an initial guess to the  *
c          solution on input. approximate solution on output           *
c          = 0 initially                                               *
c eps    = tolerance for stopping criterion.                           *
c maxits = maximum number of iterations allowed                        *
c iout   = output unit number number for printing intermediate results *
c          if (iout .le. 0) nothing is printed out.                    *
c nA,a,ja,ia                                                           *
c        = the input matrix in compressed sparse row format:           *
c          a(1:nnz)  = nonzero elements of A stored row-wise in order  *
c          ja(1:nnz) = corresponding column indices.                   *
c          ia(1:n+1) = pointer to beginning of each row in a and ja.   *
c          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
c nB,b,jb,ib                                                           *
c  similar to A matrix
c
c alu,jlu = A matrix stored in Modified Sparse Row format containing   *
c           the L and U factors, as computed by subroutine ilut.       *
c ju      = integer array of length n containing the pointers to       *
c           the beginning of each row of U in alu, jlu as computed     *
c           by subroutine ILUT.                                        *
c iopt    = option index for the solver
c                                                                      *
c on return:                                                           *
c sol    = contains an approximate solution (upon successful return).  *
c ierr   = integer. Error message with the following meaning.          *
c          ierr = 0 --> successful return.                             *
c          ierr = 1 --> convergence not achieved in itmax iterations.  *
c          ierr =-1 --> the initial guess seems to be the exact        *
c                       solution (initial residual computed was zero)  *
c globits= Global counter for iterations in nonlinear solve            *
c res    = residual of the linear system                               *
c                                                                      *
c work arrays:                                                         *
c vv     = work array of length  n x (im+2) (store the Arnoli basis)   *
c----------------------------------------------------------------------*
      integer kmax,nf
      parameter (kmax=40)
      real*8 hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
      real*8 epsmac
      data epsmac/1.d-16/
      integer i,j,i1,k,k1,its,ii,jj
      real*8 ro,ro0,eps1,gam
      double precision dnrm2,ddot

      integer maxrwk,ioptp,maxitp
      real*8 rwork(*),epsp,smaxp
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(GMRES_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      its = 0
c-----compute initial residual vector --------------
c       call amuxu (vv, n,sol,nA,a,ja,ia, nB,b,jb,ib )
c      do j=1,n
c         vv(j,1) = rhs(j) - vv(j,1)
c      enddo
c-----for zero initial guess
      do j=1,n
         vv(j,1) = rhs(j)
      enddo
c
10    continue
      ro = dnrm2( n, vv, 1)
      ro0 = ro
      if (ro .eq. 0.0d0) goto 999
      t = 1.0d0/ro
      call dscal(n, t, vv(1,1), 1)
      if (its .eq. 0) eps1=eps
c      if (its .eq. 0) eps1=eps*ro
c     ** initialize 1-st term of rhs of hessenberg system
      rs(1) = ro
      i = 0
20    i = i+1
      its = its + 1
      i1 = i + 1
      call lusolu (rhs, n,vv(1,i),alu,jlu,ju,
     &             nA,a,ja,ia, nB,b,jb,ib, iopt,
     &             rwork(nA+1),maxrwk-nA,rwork,
     &             sx,isx,jsx,slu,jslu,jsu,
     &             ioptp,maxitp,epsp,smaxp,iout )
      call amuxu (vv(1,i1), n,rhs,nA,a,ja,ia, nB,b,jb,ib )
c-----------------------------------------
c     modified gram - schmidt...
c-----------------------------------------
      do j=1,i
         t = ddot(n, vv(1,j),1,vv(1,i1),1)
         hh(j,i) = t
         call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
      enddo
      t = dnrm2( n, vv(1,i1), 1)
      hh(i1,i) = t
      if ( t .ne. 0.0d0) call dscal(n, 1.d0/t, vv(1,i1), 1)
c
c     done with modified gram schimd and arnoldi step..
c     now  update factorization of hh
c
      do k=2,i
c--------perfrom previous transformations  on i-th column of h
         k1 = k-1
         t = hh(k1,i)
         hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
         hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
      enddo
      gam = dsqrt(hh(i,i)*hh(i,i) + hh(i1,i)*hh(i1,i))
c
c     if gamma is zero then any small value will do...
c     will affect only residual estimate
c
      if (gam .eq. 0.0d0) gam = epsmac
c
c     get next plane rotation
c
      c(i) = hh(i,i)/gam
      s(i) = hh(i1,i)/gam
      rs(i1) = -s(i)*rs(i)
      rs(i) =  c(i)*rs(i)
c
c     determine residual norm and test for convergence-
c
      hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
      ro = dabs(rs(i1))
c      write(*,*) i, ro
c     write(9,*) globits, ro
      globits = globits + 1
      if (dabs(ro).gt.1.d20) then
         write(iout,*) 'blow-up in the solver, res=',ro
         stop
      endif
      if ( (i.lt.im .and. ro.gt.eps1) .or. its.lt.maxits/10)  goto 20
      if (iout .gt. 0) write(iout, 199) its, ro0, ro, eps1
c
c     now compute solution. first solve upper triangular system.
c
      rs(i) = rs(i)/hh(i,i)
      do ii=2,i
         k=i-ii+1
         k1 = k+1
         t=rs(k)
         do j=k1,i
            t = t-hh(k,j)*rs(j)
         enddo
         rs(k) = t/hh(k,k)
      enddo
c
c     form linear combination of v(*,i)'s to get solution
c
      t = rs(1)
      do k=1,n
         rhs(k) = vv(k,1)*t
      enddo
      do j=2,i
         t = rs(j)
         call daxpy(n,t,vv(1,j),1,rhs,1)
      enddo
c:q
c     call preconditioner.
c
      call lusolu (vv(1,im+2), n,rhs,alu,jlu,ju,
     &             nA,a,ja,ia, nB,b,jb,ib, iopt,
     &             rwork(nA+1),maxrwk-nA,rwork,
     &             sx,isx,jsx,slu,jslu,jsu,
     &             ioptp,maxitp,epsp,smaxp,iout )
      do k=1,n
        sol(k) = sol(k) + vv(k,im+2)
      enddo
c
c     restart outer loop when necessary
c
      if ( ro.le.eps1 ) goto 990
      if ( its.ge.maxits ) goto 991
c
c     else compute residual vector and continue..
c
      do j=1,i
         jj = i1-j+1
         rs(jj-1) = -s(jj-1)*rs(jj)
         rs(jj) = c(jj-1)*rs(jj)
      enddo
      do j=1,i1
         t = rs(j)
         if (j .eq. 1)  t = t-1.0d0
         call daxpy(n, t, vv(1,1), 1, vv(1,j), 1)
      enddo
c
199   format('   iters =',i3,' norms=',1pe12.4,' ->',1pe12.4
     &       ,' eps=',1pe12.4)
c     restart outer loop.
      goto 10
990   ierr = 0
C     Communicate residual
      res = ro
      go to 100
991   ierr = 1
      res = ro
      goto 100
999   ierr = -1
      res = ro
c
100   continue
c
cC     Support for Petsc timing
c      call PLogEventEnd(GMRES_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
      return
      end
c
c***********************************************************************
      subroutine bcgstabu(sol, n,rhs,nA,a,ja,ia, nB,b,jb,ib, alu,jlu,ju,
     &                    eps,maxits,iout,ierr,res,iopt,r0,z0,p0,rb,s0,
     &                    y0,v0,rwork,maxrwk, ioptp,maxitp,epsp,smaxp,
     &                    sx,isx,jsx,slu,jslu,jsu )
c
c     right-preconditioned Bi conjugate gradient stabilized
c
c     n      = in       order of the matrix.
c     rhs    = in       right-hand side vector.
c     sol    = inout    on input is initial guess
c                       on output is final solution
c     eps    = in       convergence criterion
c     maxits = in       max. number of iterations
c     iter   = out      number of iterations used
c     err    = out      error estimate of solution
c     iout   = in       unit number
c     iopt   = option index for the solver
c     res    = residual of the linear system
c
      implicit none
      integer n, maxits, iout, ierr, jlu(*),ju(*), iopt
      integer nA,ja(*),ia(nA+1), nB,jb(*),ib(nB+1)
      real*8 rhs(n), sol(n), a(*), b(*), alu(*), eps,res
      real*8 r0(n), z0(n), p0(n),rb(n),s0(n),y0(n),v0(n)
      integer isx(*),jsx(*), jsu(*),jslu(*)
      real*8 sx(*),slu(*)
c
      integer maxrwk,ioptp,maxitp
      real*8 rwork(*),epsp,smaxp
c----------------------------------------------------------------------*
      integer i,j,k,k1,i1,its,ii,jj,iter
      real*8 ro,ro0,gam,beta,alpha,omega,omegan,omegad
      real*8 rho,rhold,alphad
      double precision ddot,dnrm2
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(BICG_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      its = 0
c
c-------------- compute initial residual vector --------------
      beta = dnrm2(n,rhs,1)
      ro0 = beta
      do i=1,n
         r0(i) = rhs(i)
         rb(i) = r0(i)
         p0(i) = 0.0
         v0(i) = 0.0
      enddo
      rho   = 1.0
      alpha = 1.0
      omega = 1.0

10    its=its+1

      rhold = rho
      rho = ddot(n, rb,1, r0,1)
      beta = rho/rhold*alpha/omega

      do i=1,n
         p0(i) = r0(i)+beta*(p0(i)-omega*v0(i))
      enddo

      call lusolu(y0, n,p0,alu,jlu,ju,nA,a,ja,ia,nB,b,jb,ib,iopt,
     &             rwork(nA+1),maxrwk-nA,rwork,
     &             sx,isx,jsx,slu,jslu,jsu,
     &             ioptp,maxitp,epsp,smaxp,iout )

      call amuxu (v0, n,y0,nA,a,ja,ia, nB,b,jb,ib )
 
      alphad = ddot(n, rb,1, v0,1)
      alpha  = rho/alphad

      do i=1,n
         s0(i) = r0(i)-alpha*v0(i)
      enddo

      call lusolu(z0, n,s0,alu,jlu,ju,nA,a,ja,ia,nB,b,jb,ib,iopt,
     &             rwork(nA+1),maxrwk-nA,rwork,
     &             sx,isx,jsx,slu,jslu,jsu,
     &             ioptp,maxitp,epsp,smaxp,iout )

      call amuxu (r0, n,z0,nA,a,ja,ia, nB,b,jb,ib )

      omegan = ddot(n, r0,1, s0,1)
      omegad = ddot(n, r0,1, r0,1)
      omega  = omegan/omegad

      do i=1,n
         sol(i) = sol(i)+alpha*y0(i)+omega*z0(i)
         r0(i)  = s0(i)-omega*r0(i)
      enddo
      ro = dnrm2(n,r0,1)
c      print*,'ro=',ro

      if ( ro.le.eps ) goto 990
      if ( its.ge.maxits ) goto 991

199   format('   iters =',i3,' norms=',1pe12.4,' ->',1pe12.4
     &       ,' eps=',1pe12.4)

      goto 10
990   ierr = 0
      if (iout .gt. 0) write(iout, 199) its, ro0, ro,eps
      res = ro
      goto 100
c
991   ierr = 1
      if (iout .gt. 0) write(iout, 199) its, ro0, ro,eps
      res = ro
      goto 100
c
999   ierr = -1
      res = ro
c
100   continue
c
cC     Support for Petsc timing
c      call PLogEventEnd(BICG_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
      return
      end
c
c
c***********************************************************************
      subroutine amuxu (y, n,x,nA,a,ja,ia, nB,b,jb,ib )
c-----------------------------------------------------------------------
c     matrix and vector product for saddle-point system
c        y = | A B^t | * x
c            | B 0   |
c
c     input:
c       n          = full dimension of the eqiuations
c       x          = array of length n
c       nA,a,ja,ia = input matrix A in compressed sparse row format.
c       nB,b,jb,ib = input matrix B for saddle point problem
c
c     output:
c       y          = array of length n
c-----------------------------------------------------------------------
      integer n, nA,ja(*),ia(nA+1), nB,jb(*),ib(nB+1)
      real*8  x(n), y(n), a(*), b(*), tmp,t,c
      integer i,j,k
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(MATMULT_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c     y = A*x 
c     -------
      do i= 1,nA
         tmp = 0.d0
         do k=ia(i),ia(i+1)-1
            tmp = tmp + a(k)*x(ja(k))
         enddo
         y(i) = tmp
      enddo
c
c     y = B*x and y = y + (B^t)*x
c     -----------------------
      do i=1,nB
         tmp = 0.d0
         c   = x(nA+i)
         do k=ib(i),ib(i+1)-1
            j = jb(k)
            t = b(k)
            tmp  = tmp  + t*x(j)
            y(j) = y(j) + t*c
         enddo
         y(nA+i) = tmp
      enddo
c
cC     Support for Petsc timing
c      call PLogFlops(ia(n+1)*2)
c      call PLogEventEnd(MATMULT_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
      return
      end
c
c
c***********************************************************************
      subroutine lusolu (x, n,y,alu,jlu,ju, nA,a,ja,ia,nB,b,jb,ib, iopt,
     &                   rwork,maxrwk,x1,
     &                   sx,isx,jsx,slu,jslu,jsu,
     &                   ioptp,maxitp,epsp,smaxp,iout)
c-----------------------------------------------------------------------
c     preconditioner solver for saddle-point problem
c              | F B^t | * x = y
c              | 0 -S  |
c
c     iopt/100 = 0 : global system, B=0
c     iopt/100 = 2 : block LU precoditioner, see notes
c     iopt/100 = 3 : BfB^t preconditioner S=(BB^t)(BAB^t)(^-1)(BB^t)
c     iopt/100 = 4 : approximate Schur complement S=B(F^-1)B^t
c
c     return with solution x(n)
c
c     NOTE: y should return unchanged!
c-----------------------------------------------------------------------
      implicit none
      integer n,jlu(*),ju(*)
      integer nA,ja(*),ia(nA+1), nB,jb(*),ib(nB+1)
      integer iopt,maxrwk,ioptp,maxitp,iout
      real*8 x(n),y(n),alu(*),a(*),b(*), x1(nA),rwork(*),epsp,smaxp
      integer isx(*),jsx(*),jslu(*),jsu(*)
      real*8 sx(*),slu(*)
c
c     local variables
      double precision tmp
      integer i,k,j
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(LUSOL_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      if ( iopt/100.eq.2 ) then
c
c        Block LU preconditioner
c        solve for | F  0 |   | I (F^-1)(B^t) |       | y1 |
c                  | B -S | * | 0    I        | * x = | y2 |
c        ---------------------------------------------------
c        step 1, F solver -> F * x1 = y1
         call Fsolver (x, nA,y,alu,jlu,ju,iopt )
c
c        step 2, S solver -> S * x2 = -y2 + B * x1
         do i=1,nB
            tmp = -y(nA+i)
            do k=ib(i),ib(i+1)-1
               j = jb(k)
               tmp = tmp + b(k)*x(j)
            enddo
            x1(i) = tmp
         enddo
         call Ssolver (x(nA+1), nB,x1,ioptp,maxitp,epsp,smaxp,
     &                 sx,isx,jsx,slu,jslu,jsu,nA,a,ja,ia,nB,b,jb,ib,
     &                 rwork,maxrwk,iout )
c
c        step 3, F solver -> F * x1 = y1 - B^t * x2
         do i=1,nA
            x1(i) = y(i)
         enddo
         do i=1,nB
            tmp = x(nA+i)
            do k=ib(i),ib(i+1)-1
               j = jb(k)
               x1(j) = x1(j) - b(k)*tmp
            enddo
         enddo
         call Fsolver (x, nA,x1,alu,jlu,ju,iopt )
c
         goto 100
      endif

c     -------------------------------------
c     Step 1: solver for S -> x1 = (S^-1)*y
c     -------------------------------------
      if ( iopt/100.eq.3 ) then
c
c        BfBt preconditioner
c        -------------------
         call Ssolver (x1, nB,y(nA+1),ioptp,maxitp,epsp,smaxp,
     &                 sx,isx,jsx,slu,jslu,jsu,nA,a,ja,ia,nB,b,jb,ib,
     &                 rwork,maxrwk,iout )
c
c        x=(B^t)*x1
         do i=1,nA
            x(i) = 0.d0
         enddo
         do i=1,nB
            tmp = x1(i)
            do k=ib(i),ib(i+1)-1 
               j = jb(k)
               x(j) = x(j) + b(k)*tmp
            enddo
         enddo
c
c        x1=A*x
         do i=1,nA
            tmp = 0.d0
            do k=ia(i),ia(i+1)-1
               j = ja(k)
               tmp = tmp + a(k)*x(j)
            enddo
            x1(i) = tmp
         enddo
c
c        x=B*x1
         do i=1,nB
            tmp = 0.d0
            do k=ib(i),ib(i+1)-1
               tmp = tmp + b(k)*x1(jb(k))
            enddo
            x(i) = tmp
         enddo
c
         call Ssolver (x1, nB,x,ioptp,maxitp,epsp,smaxp,
     &                 sx,isx,jsx,slu,jslu,jsu,nA,a,ja,ia,nB,b,jb,ib,
     &                 rwork,maxrwk,iout )
c
      elseif ( iopt/100.eq.4 ) then
c
c        solve B*(F^-1)*(B^t)
c        --------------------
         call Ssolver (x1, nB,y(nA+1),ioptp,maxitp,epsp,smaxp,
     &                 sx,isx,jsx,slu,jslu,jsu,nA,a,ja,ia,nB,b,jb,ib,
     &                 rwork,maxrwk,iout )
c
      endif
c
c     ------------------
c     Step 2a: set value 
c     ------------------ 
      do i=1,nB
         x(nA+i) = -x1(i)
      enddo
c
c     ------------------------------------------
c     Step 2b: update vector -> x1 = y - (B^t)*x
c     ------------------------------------------
      do i=1,nA
         x1(i) = y(i)
      enddo
      do i=1,nB
         tmp = x(nA+i)
         do k=ib(i),ib(i+1)-1
            j = jb(k)
            x1(j) = x1(j) - b(k)*tmp
         enddo
      enddo
c
c     ------------------------------------------
c     Step 3: solver for F step -> x = (F^-1)*x1
c     ------------------------------------------
      call Fsolver(x, nA,x1,alu,jlu,ju,iopt )
c
100   continue
c
cC     Support for Petsc timing
c      call PLogFlops((jlu(nf+1)+ia(n+1)-ia(nf+1))*2)
c      call PLogEventEnd(LUSOL_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
cc
      return
      end
c
c
c***********************************************************************
      subroutine Fsolver(sol, n,rhs,alu,jlu,ju,iopt )
c
c     general LU solver for systems like: F * sol = rhs
c
c     where : alu,jlu = LU matrix in MSR format
c             ju      = pointer to the diagonal in jlu
c-----------------------------------------------------------------------
      implicit none
      integer n,jlu(*),ju(*),iopt
      double precision sol(*), rhs(*),alu(*)
c
      integer i,j,k
      double precision tmp
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(LUSOL_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c
      if ( mod(iopt,10).eq.1 ) then
c
c        diagonal preconditioning
         do i=1,n
            sol(i) = rhs(i)*alu(i)
         enddo

      else
c
c        LU solver - forward step
         sol(1) = rhs(1)
         do i=2,n
            tmp = rhs(i)
            do k=jlu(i),ju(i)-1
               tmp = tmp - alu(k)*sol(jlu(k))
            enddo
            sol(i) = tmp
         enddo
c
c        LU solver - backward step
         sol(n) = sol(n)*alu(n)
         do i=n-1,1,-1
            tmp = sol(i)
            do k=ju(i),jlu(i+1)-1
               tmp = tmp - alu(k)*sol(jlu(k))
            enddo
            sol(i) = tmp*alu(i)
         enddo
      endif
c
cC     Support for Petsc timing
c      call PLogFlops(jlu(n+1)*2)
c      call PLogEventEnd(LUSOL_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c
c***********************************************************************
      subroutine solver (sol, n,rhs, a,ja,ia,nelt,
     &                   iopt,lfil,maxits,eps,iout,old,ierr,
     &                   rwork,lenrw,iwork,leniw,neltv)
c
c     genral purpose CG solver
c
c     INPUT:
c       n = dimension of the equation
c       rhs(n) = right-hand-side vector of the equation
c       a(nelt),ja(nelt),ia(n+1) = matrix in CSR format
c       nelt = length of the matrix
c       iopt = option for the solver
c            =  0  -- LDU direct solver
c            =  1  -- CG/diagonal preconditioning
c            =  2  -- CG/ilu0 preconditioning
c            =  3  -- CG/ilut preconditioning
c       lfil = number of fill-ins in the ilut preconditioner
c       maxits = maximum number of iterations allowed in the solver
c       eps  = tolerance for stopping criterion - aboslute
c       iout = write unit number
c       old  = .false. -> new matrix, preconditioning is needed
c            = .true. -> old matrix, use existing preconditioner
c       rwork(lenrw), iwork(leniw) = working arrays
c
c     OUTPUT:
c       sol(n)  = solution vector
c       neltv   = length of the matrix for direct solver
c       ierr    = error index
c
c     Note:
c       if filter is called, matrix a,ia,ja is modified on return
c
      implicit none
      integer ii
      integer n,lfil,maxits,nelt,iout,ierr,ja(nelt),ia(n+1)
      integer iopt,lenrw,leniw,neltv,iwork(leniw)
      real*8  rhs(n),sol(n),a(nelt),eps,rwork(lenrw)
      logical old, lfilter
c
c     local variables
      integer lenalu,lenju,lenw
      integer localu,locw,locrw,locjlu,locju,locz,locp
      integer locjw,lociw, i,k1,k2,k
      double precision  tol,tnorm
      double precision  dnrm2
      integer iopt1,iopt2
c
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(SOLVER_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      ierr = 0
      iopt1 = mod(iopt,10)
c
c     storage management
c     ------------------
      if ( iopt1.eq.0 ) then
         lenalu = neltv
         lenju  = n
         lenw   = n 
      elseif ( iopt1.eq.1 ) then
         lenalu = n
         lenju  = 0
         lenw   = 0
      elseif ( iopt1.eq.2 ) then
         lenalu = nelt + 1
         lenju  = n
         lenw   = n
      elseif ( iopt1.eq.3 ) then
         lenalu = nelt + 2*lfil*n + 1
         lenju  = n
         lenw   = 2*n
      endif

      locjlu = 1
      locju  = locjlu + lenalu
      locjw  = locju  + lenju
      lociw  = locjw  + lenw

      if ( lociw-1.gt.leniw ) then
         ierr = 10
         write(iout,*) 'solver: iwork too short, current size=',
     &                 leniw,lociw-1,' is needed'
      endif
c
      localu = 1
      locw   = localu + lenalu
      locz   = locw   + n
      locp   = locz   + n
      locrw  = locp   + n

      if ( locrw-1.gt.lenrw ) then
         ierr = 10
         write(iout,*) 'solver: rwork too short, current size=',
     &                 lenrw,locrw-1,' is needed'
      endif
      if( ierr.ne.0 ) return
c
      if ( old ) goto 10
c
c     preconditioning
c     ---------------
      if ( iopt.eq.0 ) then
         call lducom (n, a, ja, ia, rwork(localu), iwork(locjlu),
     &                iwork(locju),rwork(locw),iwork(locjw),neltv,ierr)
c
      elseif ( iopt1.eq.1 ) then
         do i=1,n
           if ( a(ia(i)).eq.0.d0 ) then
             rwork(localu-1+i) = 1.d0
           else
             rwork(localu-1+i) = 1.d0/dabs(a(ia(i)))
           endif
         enddo

      elseif ( iopt1.eq.2) then
         call ilu0 (n, a, ja, ia, rwork(localu), iwork(locjlu), 
     &                          iwork(locju), iwork(locjw), nelt,ierr)

      elseif ( iopt1.eq.3 ) then
        tol = 5.d-3
        call ilut (n, a, ja, ia, lfil,tol, rwork(localu),iwork(locjlu),
     &             iwork(locju),lenalu, rwork(locw),iwork(locjw), ierr)
c
      endif
c
c     solving
c     -------
10    continue
c
      if ( iopt.eq.0 ) then
         call Fsolver (sol, n,rhs,rwork(localu),iwork(locjlu),
     &                 iwork(locju), iopt)
c
      elseif ( iopt/10.eq.0 ) then
         call dcg(sol, n,rhs,a,ja,ia,rwork(localu),iwork(locjlu),
     &            iwork(locju), iopt,eps,maxits,iout, 
     &            rwork(locw),rwork(locz),rwork(locp),ierr )
c
      endif
c
c
cC     Support for Petsc timing
c      call PLogEventEnd(SOLVER_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end

c
c***********************************************************************
      subroutine dcg (sol, n,rhs, aa,ja,ia,alu,jlu,ju,
     &                          iopt,eps,maxits, iout, r,z,p, ierr )
c
c     preconditioned conjugate gradient method for symmetric systems
c
c INPUT
c    n      = order of the matrix.
c    rhs    = right-hand side vector.
c    iopt   = index for the solver
c    eps    = convergence criterion
c    maxits = max. number of iterations
c    iout   = unit number
c    res    = residual of the linear system
c    r(n),z(n),p(n) = work arrays
c
c OUTPUT
c    sol    = final solution (initial guess on input)
c    ierr   = error index
c --------------------------------------------------------------------
      implicit none
      integer n, iopt,maxits,iout,ierr, ja(*),ia(n+1),jlu(*),ju(n)
      real*8 rhs(n), sol(n), aa(*), alu(*), eps
      real*8 r(n), z(n), p(n)
      integer iter, i
      real*8 dnrm2,ddot,err0, err, rnew,rold, beta, rq,alpha
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(CG_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      iter = 0
      call amux (r, n,sol,aa,ja,ia)
      do i=1,n
         r(i) = rhs(i) - r(i)
      enddo
      err0 = dnrm2(n,r,1)
      if ( err0.eq.0.d0 ) then
         err=0.d0
         goto 200
      endif
c
      do iter=1,maxits

         call Fsolver (z, n,r,alu,jlu,ju,iopt)

c        calculate coefficient bk and direction vector p

         rnew = ddot(n,z,1,r,1)
         if ( iter.eq.1) then
            do i=1,n
               p(i) = z(i)
            enddo
         else
            beta = rnew/rold
            do i=1,n
               p(i) = z(i) + beta*p(i)
            enddo
         endif
         rold = rnew
c
c        calculate coefficient ak, new iterate x, new residual r
c        -------------------------------------------------------
         call amux (z, n,p,aa,ja,ia)
         rq    = ddot(n,p,1,z,1)
         alpha = rnew/rq
         do i=1,n
            sol(i) = sol(i) + alpha*p(i)
            r(i)   = r(i)   - alpha*z(i)
         enddo
c
c        check stopping criterion
c        ------------------------
         err = dnrm2(n,r,1)
         if (dabs(err).gt.1.d20) then
            write(iout,*) 'blow-up in the solver, res=',err
            stop
         endif
         if ( err.le.eps ) go to 200
c
      enddo
c
200   continue
      if ( iout.gt.0 ) write(iout, 199) iter, err0, err
199   format('    iters=',i3,' norms=',1pe12.4,' ->',1pe12.4)
c
cC     Support for Petsc timing
c      call PLogEventEnd(CG_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c
c***********************************************************************
      subroutine lducom (n, a,ja, ia, alu,jlu,ju, rwk, iwk, neltv, ierr)
c
      implicit none
      integer n, ja(*),ia(n+1),jlu(*),ju(n), iwk(n), neltv,ierr
      real*8 a(*), alu(*),rwk(n)
      integer ju0,i,j,k, jcol, j1,j2
      real*8 fact
c
      ju0 = n+2
      jlu(1) = ju0
      do i=1,n
        iwk(i) = 0
        rwk(i) = 0.d0
      enddo
c
      do i=1,n
c
        j1 = n
        j2 = 1
        do j=ia(i),ia(i+1)-1
          jcol = ja(j)
          rwk(jcol) = a(j)
          iwk(jcol) = jcol
          j1 = min0(j1,jcol)
          j2 = max0(j2,jcol)
        enddo
c
        do j=j1,i-1
           if ( iwk(j).ne.0 ) then
              fact = rwk(j)*alu(j)
              do k=ju(j),jlu(j+1)-1
                jcol = jlu(k)
                rwk(jcol) = rwk(jcol) - fact*alu(k)
                iwk(jcol) = jcol
              enddo
              rwk(j) = fact
              j2 = max0(j2,jcol)
           endif
        enddo
c
        do j=j1,i-1
          if ( iwk(j).ne.0 ) then
            alu(ju0) = rwk(j)
            jlu(ju0) = j
            ju0 = ju0 + 1
            iwk(j) = 0
            rwk(j) = 0.d0
          endif
        enddo
c
        if ( rwk(i).eq.0.d0 ) go to 900
        alu(i) = 1.d0/rwk(i)
        iwk(i) = 0
        rwk(i) = 0.d0
c
        ju(i) = ju0
        do j=i+1,j2
          if ( iwk(j).ne.0 ) then
            alu(ju0) = rwk(j)
            jlu(ju0) = j
            ju0 = ju0 + 1
            iwk(j) = 0
            rwk(j) = 0.d0
          endif
        enddo
        if ( ju0.gt.neltv ) go to 910
        jlu(i+1) = ju0
c
      enddo
c
      do i=1,n
        fact = alu(i)
        do j=ju(i),jlu(i+1)-1
          alu(j) = alu(j)*fact
        enddo
      enddo
      neltv = jlu(n+1)
      ierr = 0
      return
c
c     zero pivot
900   ierr = i
      return
c
c     neltv not long enough
910   ierr = -ju0
c
      return
      end
c
c***********************************************************************
      subroutine ilu0 (n, a, ja, ia, alu, jlu, ju, iw, nelt,ierr)
      implicit none
      integer n, nelt,ja(nelt), ia(n+1), ju(n), jlu(nelt+1), iw(n), ierr
      real*8 a(nelt), alu(nelt+1)
c-----------------------------------------------------------------------
c                    ***   ilu(0) preconditioner.   ***                *
c
c on entry:
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage, ordered
c
c on return:
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. 
c           The diagonal (stored in alu(1:n) ) is inverted. 
c           Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c ju      = pointer to the diagonal elements in alu, jlu.
c ierr    = integer indicating error code on return
cc         ierr = 0 --> normal return
c          ierr = k --> code encountered a zero pivot at step k.
c
c work arrays:
c iw      = integer work array of length n.
c-----------------------------------------------------------------------
      integer i,j, ii, ju0, jcol, jrow, js,jf,jm, jw, jj, j1,j2
      real*8 tl,tnorm
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(ILU0_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      ierr = 0
c
      ju0 = n+2
      jlu(1) = ju0
      do i=1,n
        iw(i) = 0
      enddo
c
c     main loop
c
      do ii=1,n
        js = ju0
c
c       copy row ii of a, ja, ia into row ii of alu, jlu matrix
c
        j1 = ia(ii)
        j2 = ia(ii+1) - 1
        tl = a(j1)
        alu(ii) = tl
        iw(ii)  = ii
        ju(ii)  = ju0
        tnorm   = dabs(tl)
c
        do j=j1+1,j2
          jcol = ja(j)
          tl   = a(j)
          if ( jcol.gt.n ) then
c           write(*,*) ii,jcol,n
          else
            alu(ju0) = tl
            jlu(ju0) = jcol
            iw(jcol) = ju0
            ju0 = ju0+1
            tnorm = tnorm + dabs(tl)
            if ( jcol.lt.ii ) ju(ii) = ju0
          endif
        enddo
        tnorm = tnorm/real(j2-j1+1)
        jlu(ii+1) = ju0
        jf = ju0-1
        jm = ju(ii)-1
c
        do j=js,jm
           jrow = jlu(j)
           tl = alu(j)*alu(jrow)
           alu(j) = tl
           do jj=ju(jrow),jlu(jrow+1)-1
             jw = iw(jlu(jj))
             if ( jw.ne.0 ) alu(jw) = alu(jw) - tl*alu(jj)
           enddo
        enddo
c
c       invert and store diagonal element.
c        if (alu(ii) .eq. 0.0d0) then
c          ierr = ii
c          goto 100
c        endif
c
c        if ( dabs(alu(ii)).lt.1.d-5*tnorm ) alu(ii)=1.d-5*tnorm
        if ( alu(ii).eq.0.d0 ) alu(ii)=tnorm
        alu(ii) = 1.0d0/alu(ii)
c
c       reset pointer iw
c
        iw(ii) = 0
        do i=js,jf
          iw(jlu(i)) = 0
        enddo
c
      enddo
c
 100   continue
c
cC     Support for Petsc timing
c      call PLogEventEnd(ILU0_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c***********************************************************************
      subroutine ilut (n, a,ja,ia, lfil,tol, alu,jlu,ju, iwk, w,jw,ierr)
      implicit none
      integer  n,ja(*),ia(n+1),jlu(*),ju(n),jw(2*n), lfil, iwk, ierr
      real*8 a(*), alu(*), w(2*n), tol
c----------------------------------------------------------------------*
c                      *** ILUT preconditioner ***                     *
c      incomplete LU factorization with dual truncation mechanism      *
c                                                                      *
c     1) Theresholding in L and U as set by tol. Any element whose size*
c        is less than some tolerance (relative to the norm of current  *
c        row in u) is dropped.                                         *
c     2) Keeping only the largest lfil+il(i) elements in the i-th row  *
c        of L and the largest lfil+iu(i) elements in the i-th row of   *
c        U where il(i), iu(i) are the original number of nonzero       *
c        elements of the L-part and the U-part of the i-th row of A    *
c                                                                      *
c input:
c n       = integer. The dimension of the matrix A.
c a,ja,ia = matrix stored in Compressed Sparse Row format.
c lfil    = integer. The fill-in parameter. 
c iwk     = integer. The minimum length of arrays alu and jlu
c tol     = real*8. Sets the threshold for dropping small terms in the
c           factorization. See below for details on dropping strategy.
c
c output:
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr  = -5   --> zero row encountered.
c
c work arrays:
c w, jw (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
c jw(n+1:2n)  stores the nonzero indicator.
c-----------------------------------------------------------------------
      integer ju0,i,ii,lenu,lenl,lenu0,lenl0,j,k,jj,nl,jrow,len,jpos,kn,
     &        len1,j1,j2
      real*8 s,fact,t,tnorm
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(ILUT_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      ju0 = n+2
      jlu(1) = ju0
      do j=1,n
        jw(n+j) = 0
      enddo
c
c     scaling factor
c     save norm in w (starting in w(n+1)). Norm= average abs value.
c
      do ii=1,n
        j1 = ia(ii)
        j2 = ia(ii+1) - 1
        tnorm = 0.0d0
        do k=j1,j2
          tnorm = tnorm+dabs(a(k))
        enddo
        if (tnorm .eq. 0.0) then
          ierr = -5
          return
        endif
        w(n+ii) = tnorm/real(j2-j1+1)
      enddo
c
c     main loop.
c
      do 500 ii=1,n
c
        j1 = ia(ii)
        j2 = ia(ii+1) - 1
        tnorm = w(n+ii)
c
c       unpack L-part and U-part of row of A in arrays w
        lenu = 1
        lenl = 0
        jw(ii) = ii
        w(ii) = 0.0
        jw(n+ii) = ii
        do j=j1,j2
          k = ja(j)
          if ( k.gt.n ) then
c           write(*,*) ii,k,n
          else
            t = a(j)
            if ( k.lt.ii ) then
              lenl = lenl+1
              jw(lenl) = k
              w(lenl) = t
              jw(n+k) = lenl
            elseif ( k.eq.ii ) then
              w(ii) = t
            else
              lenu = lenu+1
              i = ii+lenu-1
              jw(i) = k
              w(i) = t
              jw(n+k) = i
            endif
          endif
        enddo
        lenl0 = lenl
        lenu0 = lenu
c
c       eliminate previous rows 
c
        jj = 0
        nl = 0
 150    jj = jj+1
        if ( jj.gt.lenl ) goto 160
c
c       determine smallest column index & exchange
c
        jrow = jw(jj)
        k = jj
        do j=jj+1,lenl
          if (jw(j) .lt. jrow) then
            jrow = jw(j)
            k = j
          endif
        enddo
c
        if (k .ne. jj) then
          j = jw(jj)
          jw(jj) = jw(k)
          jw(k) = j
c
          jw(n+jrow) = jj
          jw(n+j) = k
c
          s = w(jj)
          w(jj) = w(k)
          w(k) = s
        endif
c
c       zero out element in row
        jw(n+jrow) = 0
c
c       get the multiplier for row to be eliminated: jrow
        fact = w(jj)*alu(jrow)
        if ( dabs(fact)*w(n+jrow).le.tol*tnorm ) goto 150
c
c       combine current row and row jrow 
c
        do k=ju(jrow),jlu(jrow+1)-1
          s = fact*alu(k)
          j = jlu(k)
          jpos = jw(n+j)
          if ( jpos.ne.0 .or. dabs(s).gt.tol*tnorm ) then
            if ( j.ge.ii ) then
              if ( jpos.eq.0 ) then
                lenu = lenu+1
                i = ii+lenu-1
                jw(i) = j
                jw(n+j) = i
                w(i) = -s
              else
                w(jpos) = w(jpos) - s
              endif
            else
              if ( jpos.eq.0 ) then
                lenl = lenl+1
                jw(lenl) = j
                jw(n+j) = lenl
                w(lenl) = -s
              else
                w(jpos) = w(jpos) - s
              endif
            endif
          endif
        enddo
c
        nl = nl+1
        w(nl) = fact
        jw(nl) = jrow
        goto 150
 160    continue
c
c       reset double-pointer to zero (U-part)
c
        do k=1,lenu
          jw(n+jw(ii+k-1)) = 0
        enddo
c
c       update l-matrix
c
        len = min0(nl,lenl0+lfil)
        call qsplit (w,jw,nl,len)
        do k=1,len
          alu(ju0) =  w(k)
          jlu(ju0) =  jw(k)
          ju0 = ju0+1
        enddo
c
c       update u-matrix 
c
        ju(ii) = ju0
        len = min0(lenu,lenu0+lfil)
        len1 = min0(lenu-lenu0,lfil)
        call qsplit(w(ii+lenu0+1),jw(ii+lenu0+1),lenu-1-lenu0,len1)
c        call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
        do k=ii+1,ii+len-1
          alu(ju0) = w(k)
          jlu(ju0) = jw(k)
          ju0 = ju0+1
        enddo
c
c       store inverse of diagonal element of u
c        if ( w(ii).eq.0.0 ) w(ii) = (0.0001+tol)*tnorm
        if ( dabs(w(ii)).lt.1.d-5*tnorm ) w(ii)=1.d-5*tnorm
        alu(ii) = 1.0d0/w(ii)
c
c       update pointer to beginning of next row of U.
        jlu(ii+1) = ju0
c
500   continue
c
      ierr = 0
c
cC     Support for Petsc timing
c      call PLogEventEnd(ILUT_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c***********************************************************************
      subroutine qsplit  (a, ind, n, ncut)
      implicit none
      integer n, ind(n), ncut
      real*8 a(n)
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a.
c-----------------------------------------------------------------------
      real*8 tmp, abskey
      integer j, itmp, first, last, mid
c
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1    mid = first
      abskey = abs(a(mid))
      do j=first+1, last
         if (abs(a(j)) .gt. abskey) then
            mid = mid+1
c           interchange
            tmp = a(mid)
            itmp = ind(mid)
            a(mid) = a(j)
            ind(mid) = ind(j)
            a(j)  = tmp
            ind(j) = itmp
         endif
      enddo
c
c     interchange
c
      tmp = a(mid)
      a(mid) = a(first)
      a(first)  = tmp
c
      itmp = ind(mid)
      ind(mid) = ind(first)
      ind(first) = itmp
c
c     test for while loop
c
      if (mid .eq. ncut) return
      if (mid .gt. ncut) then
         last = mid-1
      else
         first = mid+1
      endif
      goto 1
      end
c
c
c***********************************************************************
      subroutine amux (y, n,x,a,ja,ia) 
c-----------------------------------------------------------------------
c     Y = A * X
c     input:
c       n     = row dimension of A
c       x     = array of length equal to the column dimension of matrix A
c       a, ja, ia = input matrix in compressed sparse row format.
c     output:
c       y     = real array of length n, containing the product y=Ax
c-----------------------------------------------------------------------
      integer n, ja(*), ia(*)
      real*8  x(*), y(*), a(*), tmp
      integer i, k
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(MATMULT_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      do i= 1,n
        tmp = 0.0
        do k=ia(i),ia(i+1)-1
          tmp = tmp + a(k)*x(ja(k))
        enddo
        y(i) = tmp
      enddo
c
cC     Support for Petsc timing
c      call PLogFlops(ia(n+1)*2)
c      call PLogEventEnd(MATMULT_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
      return
      end
c
c
c***********************************************************************
      subroutine filter ( n,drptol, a,ja,ia, len, ierr)
c
      real*8 a(*),drptol
      integer ja(*),ia(n+1),n,len,ierr
c-----------------------------------------------------------------------
c     This module removes any elements whose absolute value is small
c     from an input matrix A and puts the resulting matrix back in A.
c-----------------------------------------------------------------------
c on entry:
c---------
c  n	  = integer. row dimension of matrix
c  drptol = real. drop tolerance used for dropping strategy.
c  a, ja, ia = input matrix in compressed sparse format
c  len	 = integer. the amount of space available in arrays b and jb.
c
c on return:
c---------- 
c  a, ja, ia = resulting matrix in compressed sparse format.
c 
c  ierr	= integer. containing error message.
c         ierr .eq. 0 indicates normal return
c         ierr .gt. 0 indicates that there is'nt enough
c         space is a and ja to store the resulting matrix.
c         ierr then contains the row number where filter stopped.
c  len   = actual size of a, ja
c
c       NOT in a vector or parallel mode! <- Howard Hu
c----------------------------------------------------------------------c
c           contributed by David Day,  Sep 19, 1989.                   c
c----------------------------------------------------------------------c
c local variables
      integer index,row,k,k1,k2 
c
      ierr =0
      index = 1
      do row=1,n
c
        k1 = ia(row)
        k2 = ia(row+1) - 1
        ia(row) = index
        do k=k1,k2
          if ( dabs(a(k)).gt.drptol ) then 
            a(index) =  a(k)
            ja(index) = ja(k)
            index = index + 1
          endif
        enddo
c
      enddo
      ia(n+1) = index
      len = index
c
      return
      end
