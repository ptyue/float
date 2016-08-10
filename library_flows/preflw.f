c***********************************************************************
      subroutine PREFLW (itime, itordr,dt,dt1,dtin,endtim,tsc,tsc1,
     &                   ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi,
     &                   ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                   indgu,indgp,indgtm,indgr,indgnb,indgphi,
     &                    indgpsi,indgstrm,
     &                   z,zs, u,uold,p,dudt,ncpg,ncps,ncpb,
     &                   selas,desdt,selasd,phi,phiold,psi,dphidt,
     &                   nbrigid,uvom,uvod1,duvot )
c
c     This routine prepares for the flow solver
c***********************************************************************
      implicit none
      integer ncpvelo,ncppres,ncptemp,ncpelas,
     &        ndgvelo,ndgpres,ndgtemp,ndgelas,
     &        indgu,indgp,indgtm,indgr,indgnb,indgstrm
      integer itime,nbrigid,itordr,ncpg,ncps,ncpb
      real*8 u(ncpg,ndgvelo),dudt(ncpg,ndgvelo),uold(ncpg,ndgvelo),
     &       p(ndgpres),selas(ncps,ndgelas),selasd(ncps,ndgelas),
     &       desdt(ncps,ndgelas),
     &       uvom(ncpb,nbrigid),uvod1(ncpb,nbrigid),duvot(ncpb,nbrigid)
      real*8 z(*),zs(*),dtin,endtim,tsc,tsc1,dt,dt1
c*CH start
      integer ndgphi, ncpphi, indgphi, indgpsi
      real*8 phi(ndgphi),phiold(ndgphi),psi(ndgphi),dphidt(ndgphi)
c*CH end
c
      integer i,j,n,np
c
c     constants for transient terms
c     -----------------------------
      if ( itordr.eq.1 ) then
         tsc1 = 0.d0
         tsc = 1.d0/dt
      elseif ( itordr.eq.2 ) then
         tsc1 = -dt/(dt+dt1)
         tsc = (1.d0-tsc1)/dt
      elseif ( itordr.eq.3 ) then
         tsc1 = -1.d0
         tsc = 2.d0/dt
      endif
c
c     time step
c     ---------
      dtin = dt
      endtim = dt
c
c     indx for global variables
c     -------------------------
c     Note: indgp should be the largest among indgng...indgpsi, 
c           in order to use matglb=.f.
      indgnb = 0
      indgr  = indgnb + ncpb*nbrigid
      indgu  = indgr  + ndgelas*ncpelas
      indgtm = indgu  + ndgvelo*ncpvelo
c*CH start
      indgphi= indgtm + ndgtemp*ncptemp
      indgpsi= indgphi+ ndgphi*ncpphi
c*CH end
      indgp  = indgpsi+ ndgphi*ncpphi

      indgstrm=indgp  + ndgpres*ncppres
c
c     stress components
c     -----------------
      n = indgr
      do j=1,ncpelas
         if ( itordr.ne.2 ) then
            do i=1,ndgelas
               zs(n+i) = -tsc*selas(j,i) + tsc1*desdt(j,i)
               z(n+i) = selas(j,i)
            enddo
         else
            do i=1,ndgelas
               zs(n+i) = -tsc*selas(j,i) + tsc1*desdt(j,i)
               z(n+i) = selas(j,i) + dt*desdt(j,i)
            enddo
         endif
         n = n + ndgelas
      enddo
c
c     velocities
c     ----------
      n = indgu
      do j=1,ncpvelo
         if ( itordr.ne.2 ) then
            do i=1,ndgvelo
               zs(n+i) = -tsc*u(j,i) + tsc1*dudt(j,i)
               z(n+i) = u(j,i)
            enddo
         else
            do i=1,ndgvelo
               zs(n+i) = -tsc*u(j,i) + tsc1*dudt(j,i)
               z(n+i) = u(j,i) + dt*dudt(j,i)
            enddo
         endif
         n = n + ndgvelo
      enddo
c
c     pressure
c     --------
      n = indgp
      do i=1,ndgpres
         zs(n+i) = 0.d0
         z(n+i) = p(i)
      enddo
c
c     particle velocities
c     -------------------
      n = indgnb
      do np=1,nbrigid
        do i=1,ncpb
           zs(n+1) = -tsc*uvod1(i,np) + tsc1*duvot(i,np)
           z(n+1)  = uvom(i,np)
           n = n + 1
        enddo
      enddo
c*CH start
c
c     Cahn-Hillard equation
c     ----------
      if(ncpphi.eq.1)then
         n = indgphi
         if ( itordr.ne.2 ) then
           zs(n+1:n+ndgphi) = -tsc*phi(1:ndgphi) + tsc1*dphidt(1:ndgphi)
            z(n+1:n+ndgphi) = phi(1:ndgphi)
         else
           zs(n+1:n+ndgphi) = -tsc*phi(1:ndgphi) + tsc1*dphidt(1:ndgphi)
            z(n+1:n+ndgphi) = phi(1:ndgphi) + dt*dphidt(1:ndgphi)
         endif
         n = indgpsi
            zs(n+1:n+ndgphi) = 0.d0
             z(n+1:n+ndgphi) = psi(1:ndgphi)
      endif
c*CH end
c
c     update
c     ------
      do j=1,ncpvelo
        do i=1,ndgvelo
           uold(j,i) = u(j,i)
        enddo
      enddo
      do j=1,ncpelas
        do i=1,ndgelas
           selasd(j,i) = selas(j,i)
        enddo
      enddo
c*CH start
      if(ncpphi.eq.1)phiold(1:ndgphi)=phi(1:ndgphi)
c*CH end
c
      return
      end
