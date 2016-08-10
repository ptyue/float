c***********************************************************************
      subroutine collision ( nbrigid,bdseg,eps,dt,pdiameter,pmass,
     &                       pposition,pvelocity,pacceleration, pwall,
     &                       pcollision,ncontmax, pnew,g,icont)
c
c     this routine implements colllision models
c
c     INPUT: 
c        nbrigid  = number of rigid particles
c        bdseg    = number of total boundary segments
c        eps    = relative size of the collision zone
c                 (relative to the diameter of the particles)
c        pdiameter(nbrigid)      = particle pdiametereter
c        pmass(3,nbrigid)        = particle mass
c        pposition(3,nbrigid)    = particle position
c        pvelocity(3,nbrigid)    = particle velocity
c        pacceleration(3,nbrigid)= particle acceleration
c        pwall(5,bdseg)        = information about boundary segments
c
c     OUTPUT:
c        pcollision(3,nbrigid)  = particle acceleration due to the
c                               collision forces
c
c     WORKING ARRAYS
c        pnew(2,nbrigid),icont(4,ncontmax),g(5,ncontmax)
c           g(1,ncount) = Gij
c           g(2,ncount) = (Ri+Rj+delta)**2
c           g(3,ncount) = 1/(1/mi+1/mj))
c           g(4,ncount) = nijx
c           g(5,ncount) = nijy
c-----------------------------------------------------------------------
c
      implicit none
      integer nbrigid,bdseg,ncontmax
      double precision eps,dt,pdiameter(nbrigid),pmass(3,nbrigid),
     &                 pposition(3,nbrigid),pvelocity(3,nbrigid),
     &                 pacceleration(3,nbrigid),pcollision(3,nbrigid),
     &                 pwall(5,bdseg)
      integer icont(3,ncontmax)
      double precision pnew(2,nbrigid),g(5,ncontmax)
c
c     local variables
      logical check,lwall
      integer n,i,j,iter,ncount,ncound,nj,ncod(20)
      double precision tmp,dis,res1,res2,res,eps1,pij1,pij2,delta
      double precision xcoor,xploc
c
c     initialize the index for colliding particles
c     --------------------------------------------
      if(nbrigid.ge.1) delta = 0.001*eps*(2.d0*pmass(1,1))
      eps1 = 0.5625d0
c
      do i=1,nbrigid
        pcollision(1,i) = 0.d0
        pcollision(2,i) = 0.d0
        pcollision(3,i) = 0.d0
      enddo
      do n=1,ncontmax
        icont(1,n) = 0
        icont(2,n) = 0
        icont(3,n) = 1
        g(1,n)     = 0.d0
      enddo
      ncount = 0
      ncound = 0
c
c     projected particle position
c     ---------------------------
      tmp = 0.5d0*dt*dt
      do i=1,nbrigid
        do j=1,2
           pnew(j,i) = pposition(j,i) + dt*pvelocity(j,i) +
     &                 tmp*pacceleration(j,i)
        enddo
        pnew(1,i) = xploc(pnew(1,i))
      enddo
c
c     detecting collisions between particles and particle-wall
c     --------------------------------------------------------
10    do i=1,nbrigid
c
        nj = 0
        do n=1,ncound
          if ( icont(1,n).eq.i ) then
            nj = nj+1
            ncod(nj) = icont(2,n)
          endif
        enddo
c
c       betweend particle-wall
        do j=1,bdseg
c
          check = .true.
          do n=1,nj
            if ( ncod(n).eq.j+nbrigid ) check = .false.
          enddo
          if ( check ) then
            call wcollision (i,j,nbrigid,bdseg,eps,eps1,pdiameter,pmass,
     &                    pwall,pcollision,pnew,g,icont,ncount,ncontmax)
          endif
        enddo
c
c       between particle-particle
        do j=i+1,nbrigid
          check = .true.
          do n=1,nj
            if ( ncod(n).eq.j ) check = .false.
          enddo
          if ( check) then
            tmp = pnew(1,j)
            pij1 = xcoor(tmp,pnew(1,i)) - tmp
            pij2 = pnew(2,i) - pnew(2,j)
            res1 = pij1 + pcollision(1,i) - pcollision(1,j)
            res2 = pij2 + pcollision(2,i) - pcollision(2,j)
            tmp = (0.5*pdiameter(i)+0.5*pdiameter(j)+
     &             eps*max(pdiameter(i),pdiameter(j)) )**2
            res =  res1*res1 + res2*res2 - tmp
            if ( res.lt.tmp*eps1 ) then
               ncount = ncount + 1
               g(2,ncount) = tmp
               g(3,ncount) = 1.d0/(1.d0/pmass(1,i)+1.d0/pmass(1,j))
               g(4,ncount) = pij1
               g(5,ncount) = pij2
               icont(1,ncount) = i
               icont(2,ncount) = j
            endif
          endif
        enddo
c
      enddo
c
      if ( ncount.gt.ncontmax ) stop 'colllision: ncontmax is too small'
c
c      write(*,*) 'ncount=',ncount
c      do n=1,ncound
c        write(*,*) 'n=',n,icont(1,n),icont(2,n)
c      enddo
c      do n=1+ncound,ncount
c        write(*,*) 'after recheck n=',n,icont(1,n),icont(2,n)
c      enddo
c
      if ( ncount.eq.ncound ) goto 80
c
c     calculate the contact forces
c     ----------------------------
      iter = 1
50    dis  = 0.d0
      do n=1,ncount
         i = icont(1,n)
         j = icont(2,n)
         if ( j.le.nbrigid ) then
            res1 = g(4,n) + pcollision(1,i) - pcollision(1,j)
            res2 = g(5,n) + pcollision(2,i) - pcollision(2,j)
         else
            res1 = g(4,n) + pcollision(1,i)
            res2 = g(5,n) + pcollision(2,i)
         endif
c
         res = res1*res1 + res2*res2 - g(2,n)
         res2 = res1*g(4,n) + res2*g(5,n)
         if ( abs(res2).gt.1.d-15 ) then
            res1 = -res*g(3,n)/res2
         else
            res1 = 0.d0
         endif
         if ( abs(res1).ge.delta ) then
c
           lwall = .false.
           if ( icont(3,n).ge.0 ) then
             if ( res.lt.0.d0 .or. g(1,n).ge.delta ) lwall = .true.
           else
             if ( res.gt.0.d0 .or. g(1,n).ge.delta ) lwall = .true.
           endif
           if ( lwall ) then
c
c            update contact forces Gij
c            -------------------------
             res2 = g(1,n) + res1
             if ( icont(3,n).ge.0 ) then
                if ( res2.gt.0.d0 ) then
                   g(1,n) = res2
                   res    = abs(res1)
                else
c                  particle i and j are not touching
                   res1   = -g(1,n)
                   g(1,n) = 0.d0
                   res    = 0.d0
                endif
             else
                if ( res2.lt.0.d0 ) then
                   g(1,n) = res2
                   res    = abs(res1)
                else
c                  particle i and j are not touching
                   res1   = -g(1,n)
                   g(1,n) = 0.d0
                   res    = 0.d0
                endif
             endif
c
c            update the contact forces, x and y components
c            ---------------------------------------------
             tmp = 0.5d0*res1*g(4,n)
             pcollision(1,i) = pcollision(1,i) + tmp/pmass(1,i)
             if ( j.le.nbrigid ) 
     &       pcollision(1,j) = pcollision(1,j) - tmp/pmass(1,j)
             tmp = 0.5d0*res1*g(5,n)
             pcollision(2,i) = pcollision(2,i) + tmp/pmass(2,i)
             if ( j.le.nbrigid ) 
     &       pcollision(2,j) = pcollision(2,j) - tmp/pmass(2,j)
             dis = max(dis,res)
          endif
c
        endif
c
      enddo
c
c     convergence test
c     ----------------
      if ( dis.gt.delta ) then
         iter = iter + 1
c        write(*,*) 'iter=',iter,dis
         if ( iter.gt.5000 ) stop 'collision does not converge!'
         go to 50
      endif
c
c     check new contacts
c     ------------------
      eps1 = 0.d0
      ncound = ncount
      goto 10
c
c     return with particle acceleration due to collision
c     --------------------------------------------------
 80   continue
      tmp = 2.d0/(dt*dt)
      do i=1,nbrigid
         pcollision(1,i) = pcollision(1,i)*tmp
         pcollision(2,i) = pcollision(2,i)*tmp
      enddo
c
      return
      end 
c
c
c***********************************************************************
      subroutine wcollision (i,j,nbrigid,bdseg,eps,eps1,pdiameter,pmass,
     &                    pwall,pcollision,pnew,g,icont,ncount,ncontmax)
c
c     this routine detects the colllisions between particle-wall
c***********************************************************************
      implicit none
      integer i,j,nbrigid,bdseg,ncount,ncontmax,icont(3,ncontmax)
      double precision eps,eps1,pdiameter(nbrigid),pmass(3,nbrigid),
     &                 pcollision(3,nbrigid),pwall(5,bdseg)
      double precision pnew(2,nbrigid),g(5,ncontmax)
c
c     local variables
      logical lwall
      integer ix,jm1,jp1
      double precision tmp,res1,res2,res, pij1,pij2,
     &                 x21,y21,x1,y1,x2,y2,xn1,xn2,yn1,yn2
c     --------------------------------------------------------
c
c     no checking for periodic boundary segement
c     ------------------------------------------
      if ( abs(pwall(1,j)).le.0.001) return
c
      if ( int(pwall(1,j))/100.eq.1 ) then
c     collision between particle i and circular segment j
c     ---------------------------------------------------
          pij1 = pnew(1,i) - pwall(2,j)
          pij2 = pnew(2,i) - pwall(3,j)
          res1 = pij1 + pcollision(1,i)
          res2 = pij2 + pcollision(2,i)
          tmp = (0.5d0*pdiameter(i)+pwall(4,j)+eps*pdiameter(i))**2
          res = res1*res1 + res2*res2 - tmp
          lwall = .false.
          if ( pwall(4,j).gt.0 ) then
             if ( res.lt.tmp*eps1 ) then
                lwall = .true.
                ix = 1
             endif
          else
             if ( res.gt.-0.3*eps1*pdiameter(i)**2 ) then
                lwall = .true.
                ix = -1
             endif
          endif
          if ( lwall ) then
             ncount = ncount + 1
             g(2,ncount) = tmp
             g(3,ncount) = pmass(1,i)
             g(4,ncount) = pij1
             g(5,ncount) = pij2
             icont(1,ncount) = i
             icont(2,ncount) = j+nbrigid
             icont(3,ncount) = ix
          endif
      else
c
c     collision between particle i and linear segment j
c     -------------------------------------------------
         x21  = pwall(4,j) - pwall(2,j)
         y21  = pwall(5,j) - pwall(3,j)
c
         x1  = pnew(1,i)  - pwall(2,j)
         y1  = pnew(2,i)  - pwall(3,j)
         x2  = pnew(1,i)  - pwall(4,j)
         y2  = pnew(2,i)  - pwall(5,j)
c
         xn1 = x1 + pcollision(1,i)
         xn2 = x2 + pcollision(1,i)
         yn1 = y1 + pcollision(2,i)
         yn2 = y2 + pcollision(2,i)
c
c        distance/gap bewteen the particle to the wall
         tmp  = (x1*x21 + y1*y21)/(x21**2+y21**2)
         pij1  = x1 - tmp*x21
         pij2  = y1 - tmp*y21
         res1 = pij1 + pcollision(1,i)
         res2 = pij2 + pcollision(2,i)
         tmp  = ( 0.5d0*pdiameter(i)+eps*pdiameter(i) )**2
         res  = res1*res1 + res2*res2 - tmp
c
c        particle is located at the left-hand-side of the boundary
         if( (xn1*yn2-xn2*yn1) .gt. 0.d0 ) then
c
            if( res .lt. tmp*eps1 ) then
c
               if( ((xn1*x21+yn1*y21) .ge. 0) .and. 
     &             ((xn2*x21+yn2*y21) .le. 0) ) then
                   ncount = ncount + 1
                   g(2,ncount) = tmp
                   g(3,ncount) = pmass(1,i)
                   g(4,ncount) = pij1
                   g(5,ncount) = pij2
                   icont(1,ncount) = i
                   icont(2,ncount) = j+nbrigid
               else
                  if ( (xn1*xn1+yn1*yn1-tmp) .lt. tmp*eps1 ) then
                     jm1 = j-1
                     if ( j.eq.1 ) jm1 = bdseg
                     if( abs(pwall(1,jm1)).gt.0.001 .and.
     &                   y21*(pwall(4,jm1)-pwall(2,jm1))-
     &                   x21*(pwall(5,jm1)-pwall(3,jm1)).le.0.d0 .and.
     &                   xn1*(pwall(4,jm1)-pwall(2,jm1))+
     &                   yn1*(pwall(5,jm1)-pwall(3,jm1)).gt.0.d0) then
                         ncount = ncount + 1
                         g(2,ncount) = tmp
                         g(3,ncount) = pmass(1,i)
                         g(4,ncount) = x1
                         g(5,ncount) = y1
                         icont(1,ncount) = i
                         icont(2,ncount) = j+nbrigid
                     endif
                  endif
c
                  if ( (xn2*xn2+yn2*yn2-tmp) .lt. tmp*eps1 ) then
                     jp1 = j+1
                     if ( j.eq.bdseg ) jp1 = 1 
                     if( abs(pwall(1,jp1)).gt.0.001 .and.
     &                   y21*(pwall(4,jp1)-pwall(2,jp1)) -
     &                   x21*(pwall(5,jp1)-pwall(3,jp1)).gt.0.d0 .and.
     &                   xn2*x21+yn2*y21.lt.0.d0 ) then
                         ncount = ncount + 1
                         g(2,ncount) = tmp
                         g(3,ncount) = pmass(1,i)
                         g(4,ncount) = x2
                         g(5,ncount) = y2
                         icont(1,ncount) = i
                         icont(2,ncount) = j+nbrigid
                     endif
                  endif
c
               endif
            endif
         else
c            if( (xn1*x21+yn1*y21.ge.0).and.(xn2*x21+yn2*y21.le.0) ) then
c              write(*,'(A,i3,A,A,i3)') 'STOP: particle ',i,'is located '
c     &               ,'at the right-hand-side of the boundary segment',j
c              stop
c            endif
         endif
      endif
      return
      end  
