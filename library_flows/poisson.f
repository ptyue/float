c***********************************************************************
      subroutine Ssolver (sol, n,rhs,ioptp,maxits,epsm,smax,
     &                    sx,isx,jsx,slu,jslu,jsu,nA,a,ja,ia,nB,b,jb,ib,
     &                    rwork,lenrw,iout )
c
c     this routine solves the approximate Schur complement
c
c     INPUT:
c       rhs(n) = right-hand-side vector of the equation
c       ioptp = option for the solver
c             = x1  solving S = precon of B(F^-1)B^t
c             = x2  solving for S=B(F^-1)B^t - CG
c             = x3  solving for S=BB^t - CG
c             = x5  solving for S=B(F^-1)B^t - CR
c             = x6  solving for S=B(F^-1)B^t - CR, S formed
c           x=1 -> diagonal preconditioned
c           x=2 -> ilu0 preconditioned with simplerified S
c           x=3 -> ilu0 preconditioned with full S
c       maxits = maximum number of iterations allowed in the solver
c       epsm   = tolerance for stopping criterion
c       smax   = magnitudes of the variables
c       sx,isx,jsx   = S maxtrix in CRS format
c       slu,jslu,jsu = preconditioner for S in MSR format
c       nA,a,ja,ia   = matrix A in CRS
c       nB,b,jb,ib   = matrix B in CSR
c       iout = write unit number
c
c     OUTPUT:
c       sol(n)  = solution vector
c       ierr    = error index
c
c     WORKING VARIABLES
c       rwork(lenrw) 
c         if ( ioptp = x1 )       lenrw = 0
c         if ( ioptp = x2/x3 )    lenrw = 3*n+nA
c         if ( ioptp = x5/x6 )    lenrw = 5*n+nA
c
c     Howard Hu    August 10, 1999
c***********************************************************************
      implicit none
      integer lenrw,iout,ioptp,maxits,n
      integer nA,ia(nA+1),ja(*),nB,ib(nB+1),jb(*)
      double precision a(*),b(*)
      double precision rhs(n),sol(n), rwork(lenrw),epsm,smax
      integer isx(*),jsx(*), jslu(*),jsu(*)
      real*8 sx(*),slu(*)
c
c     local variables
      integer ierr, iopt,iopt1
      double precision eps
      integer locrw,locr,locz,locp,locs,locq, locx,i
c
      ierr = 0
      iopt = mod(ioptp,10)
c
      if ( iopt.eq.1 ) then
c
c        direct solver - with preconditioner
c        -----------------------------------
         iopt1 = ioptp/10
         call Fsolver(sol, n,rhs,slu,jslu,jsu, iopt1 )

      elseif ( iopt.eq.2 .or. iopt.eq.3 ) then
c
c        solver conjudate gradient
c        -------------------------
         locr   = 1
         locz   = locr   + n
         locp   = locz   + n
         locx   = locp   + n
         locrw  = locx   + nA
         if ( locrw-1.ne.lenrw ) then
            write(iout,*)  'Ssolver: rwork too short, current size=',
     &                     lenrw,locrw-1,' is needed'
            stop
         endif
         eps = epsm*smax
c
         do i=1,n
            sol(i) = rhs(i)*slu(i)
         enddo
         call dcgS(sol, n,rhs,nA,a,ja,ia,nB,b,jb,ib,
     &             eps, ioptp,maxits, iout,ierr, 
     &             rwork(locr),rwork(locz),rwork(locp),rwork(locx),
     &             sx,isx,jsx, slu,jslu,jsu )

      elseif ( iopt.eq.5 .or. iopt.eq.6 ) then
c
c        solver conjudate residual
c        -------------------------
         locr   = 1
         locz   = locr   + n
         locp   = locz   + n
         locs   = locp   + n
         locq   = locs   + n
         locx   = locq   + n
         locrw  = locx   + nA
         if ( locrw-1.ne.lenrw ) then
            write(iout,*) 'Ssolver: rwork too short, current size=',
     &                    lenrw,locrw-1,' is needed'
            stop
         endif
         eps = epsm*smax
c
         do i=1,n
            sol(i) = 0.d0
         enddo
         call dcrS(sol, n,rhs,nA,a,ja,ia,nB,b,jb,ib,
     &             eps, ioptp,maxits, iout,ierr,
     &             rwork(locr),rwork(locz),rwork(locp),
     &             rwork(locs),rwork(locq),rwork(locx),
     &             sx,isx,jsx, slu,jslu,jsu )

      endif
      if (ierr.ne.0 ) then
         write(*,*) 'error in Ssolver, ierr=', ierr
         stop
      endif
c
      return
      end

c*********************************************************************
      subroutine dcgS (sol, n,rhs,nA,a,ja,ia,nB,b,jb,ib,
     &                 eps, ioptp,maxits, iout,ierr, r,z,p,x,
     &                 sx,isx,jsx,slu,jlu,ju )
c
c     preconditioned conjugate gradient method for symmetric systems
c
c INPUT 
c    n      = order of the matrix.
c    rhs    = right-hand side vector.
c    ioptp  = index for the solver
c    eps    = convergence criterion
c    maxits = max. number of iterations
c    iout   = unit number
c    r(n),z(n),p(n),x(nA) = work arrays
c
c OUTPUT  
c    sol    = final solution (initial guess on input)
c    ierr   = error index
c --------------------------------------------------------------------
      implicit none
      integer maxits, iout, ierr
      integer nA,ia(*),ja(*),nB,ib(*),jb(*)
      real*8 rhs(*), sol(*), slu(*), eps,res
      real*8 r(*), z(*), p(*), x(nA)
      integer n,iter, i,ioptp,iopt1
      real*8 dnrm2,ddot,err0, err, rnew,rold, beta, rq,alpha
      real*8 a(*),b(*)
      integer isx(*),jsx(*), jlu(*),ju(*)
      real*8 sx(*)
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(CG_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      iopt1 = ioptp/10
      iter = 0
      if ( mod(ioptp,10).eq.6 ) then
        call amux (r, n,sol,sx,jsx,isx )
      else
        call amuxS (r, sol,nA,a,ja,ia, nB,b,jb,ib, ioptp, x)
      endif

      do i=1,n
         r(i) = rhs(i) - r(i)
      enddo
      err0 = dnrm2(n,r,1)
      if ( err0.eq.0.d0 ) goto 200
c
      do iter=1,maxits-1

         call Fsolver(z, n,r,slu,jlu,ju,iopt1 )
c
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
         if ( mod(ioptp,10).eq.6 ) then
            call amux (z, n,p,sx,jsx,isx )
         else
            call amuxS (z, p,nA,a,ja,ia, nB,b,jb,ib, ioptp, x)
         endif

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
200   if ( iout.gt.0 ) write(iout, 199) iter, err0, err
199   format(' dcgS: iters=',i3,' norms=',1pe12.4,' ->',1pe12.4)
      res = err
c
cC     Support for Petsc timing
c      call PLogEventEnd(CG_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c
c*********************************************************************
      subroutine dcrS (sol, n,rhs,nA,a,ja,ia,nB,b,jb,ib,
     &                 eps, ioptp,maxits, iout,ierr, r,z,p,s,q,x,
     &                 sx,isx,jsx,slu,jlu,ju )
c
c     preconditioned conjugate residual method for symmetric systems
c
c INPUT 
c    n      = order of the matrix.
c    rhs    = right-hand side vector.
c    ioptp  = index for the solver
c    eps    = convergence criterion
c    maxits = max. number of iterations
c    iout   = unit number
c    r(n),z(n),p(n),s(n),q(n),x(nA) = work arrays
c
c OUTPUT
c    sol    = final solution (initial guess on input)
c    ierr   = error index
c --------------------------------------------------------------------
      implicit none
      integer maxits, iout, ierr
      integer nA,ia(*),ja(*),nB,ib(*),jb(*)
      real*8 rhs(*), sol(*), slu(*), eps,res
      real*8 r(*), z(*), p(*), s(*), q(*), x(nA)
      integer n,iter, i,ioptp,iopt1
      real*8 dnrm2,ddot,err0, err, rnew,rold,beta,rq,alpha

      real*8 a(*),b(*)
      integer isx(*),jsx(*),jlu(*),ju(*)
      real*8 sx(*)
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(CG_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      iopt1 = ioptp/10
      iter = 0
      if ( mod(ioptp,10).eq.6 ) then
        call amux (r, n,sol,sx,jsx,isx )
      else
        call amuxS (r, sol,nA,a,ja,ia, nB,b,jb,ib, ioptp, x)
      endif

      do i=1,n
         r(i) = rhs(i) - r(i)
      enddo
      err0 = dnrm2(n,r,1)
      if ( err0.eq.0.d0 ) goto 200
c
      do iter=1,maxits-1

         call Fsolver(z, n,r,slu,jlu,ju,iopt1 )

         if ( mod(ioptp,10).eq.6 ) then
            call amux (s, n,z,sx,jsx,isx )
         else
            call amuxS (s, z,nA,a,ja,ia, nB,b,jb,ib, ioptp, x)
         endif

         rnew = ddot(n,r,1,s,1)
c
         if ( iter.eq.1) then
            do i=1,n
               p(i) = z(i)
               q(i) = s(i)
            enddo
         else
            beta = rnew/rold
            do i=1,n
               p(i) = z(i) + beta*p(i)
               q(i) = s(i) + beta*q(i)
            enddo
         endif
         rold = rnew
c
         rq    = ddot(n,q,1,q,1)
         alpha = rnew/rq
         do i=1,n
            sol(i) = sol(i) + alpha*p(i)
            r(i)   = r(i)   - alpha*q(i)
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
200   if ( iout.gt.0 ) write(iout, 199) iter, err0, err
199   format(' dcrS: iters=',i3,' norms=',1pe12.4,' ->',1pe12.4)
      res = err
c
cC     Support for Petsc timing
c      call PLogEventEnd(CG_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      return
      end
c
c***********************************************************************
      subroutine amuxS (y, x,nA,a,ja,ia,nB,b,jb,ib, ioptp, r )
c
c     Y = (B*BT) * X
c
c     input:
c       n     = row dimension of A
c       x     = array of length equal to the column dimension of matrix A
c       nA,a,ja,ia = input matrix in compressed sparse row format.
c       nB,b,jb,ib
c
c     output:
c       y     = real array of length n, containing the product y=Ax
c-----------------------------------------------------------------------
      implicit none
      integer nA,ja(*),ia(*),nB,ib(*),jb(*), ioptp
      real*8  x(*), y(*), a(*), b(*), tmp, r(*)
      integer i,j, k
c
cC     Support for Petsc timing
c      include '../include/petsc.h'
cc
c      call PLogEventBegin(MATMULT_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
c       r=(B^t)*x
c       --------
        do i=1,nA
           r(i) = 0.d0
        enddo
        do i=1,nB
           tmp = x(i)
           do k=ib(i),ib(i+1)-1
              j = jb(k)
              r(j) = r(j) + b(k)*tmp
           enddo
        enddo
c
c       r=r/diag(A)
c       -----------
        if ( mod(ioptp,10).ne.3 ) then
           do i=1,nA
              r(i) = r(i)/a(ia(i))
           enddo
        endif
c
c       y=B*r
c       -----
        do i=1,nB
           tmp = 0.d0
           do k=ib(i),ib(i+1)-1
              j = jb(k)
              tmp = tmp + b(k)*r(j)
           enddo
           y(i) = tmp
        enddo

c
cC     Support for Petsc timing
c      call PLogEventEnd(MATMULT_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
      return
      end
c
c
c***********************************************************************
      subroutine FormS (nA,a,ja,ia, nB,b,jb,ib, sx,jsx,isx, ioptp, r)
c
c     this routine forms the approximate Schur complement
c
c-----------------------------------------------------------------------
      implicit none
      integer nA,ja(*),ia(*),nB,jb(*),ib(*),jsx(*),isx(*),ioptp
      double precision a(*),b(*),sx(*),r(nA)
c
      integer i,j,k,ik,jk
      double precision tmp
c
      do i=1,nA
         r(i) = 0.d0
      enddo
c
      do i=1,nB

         do k=ib(i),ib(i+1)-1
            r(jb(k)) = b(k)
         enddo

         do jk=isx(i),isx(i+1)-1
            j = jsx(jk)
c
            tmp = 0.d0
            do k=ib(j),ib(j+1)-1
               ik = jb(k)
               if ( mod(ioptp,10).eq.3 ) then
                  tmp = tmp + b(k)*r(ik)
               else
                  tmp = tmp + b(k)*r(ik)/a(ia(ik))
               endif
            enddo
            sx(jk) = tmp
c
         enddo

         do k=ib(i),ib(i+1)-1
            r(jb(k)) = 0.d0
         enddo
      enddo
c
c      write(20,'(10i7)') (isx(i),i=1,nB+1)
c      do i=1,nB
c        write(20,'(4(i7,e12.5))') (jsx(k),sx(k), k=isx(i),isx(i+1)-1)
c      enddo
      return
      end
