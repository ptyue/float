c***********************************************************************
      subroutine PSTFLW ( itordr, dt,tsc, z,zs,
     &                    ndgvelo,ndgpres,ndgelas,ndgstrm,
     &                    indgu,indgp,indgr,indgnb,indgstrm,
     &                    ncpvelo,ncpelas,ncppres,u,uold,p,dudt,
     &                    ncpg,ncps,ncpb,selas,desdt,selasd,strm,
     &                    nbrigid,uvom,duvot,uvod1,
     &                    ncpphi,ndgphi, indgphi,indgpsi,phi,phiold,psi,
     &                     dphidt)
c
c     This routine recovers forces, velocities, pressure and stesses
c***********************************************************************
      implicit none
      integer ndgvelo,ndgpres,ndgelas,ndgstrm,
     &        indgu,indgp,indgr,indgnb,indgstrm,ncppres,
     &        ncpvelo,ncpelas,nbrigid, itordr, ncpg,ncps,ncpb
      real*8 u(ncpg,ndgvelo),dudt(ncpg,ndgvelo),
     &       uold(ncpg,ndgvelo),p(ndgpres),strm(ndgstrm),
     &       selas(ncps,ndgelas),desdt(ncps,ndgelas),
     &       selasd(ncps,ndgelas),
     &       uvom(ncpb,nbrigid),duvot(ncpb,nbrigid),uvod1(ncpb,nbrigid)
      real*8 z(*),zs(*),tsc,dt
      integer ndgphi,indgphi,indgpsi,ncpphi
      real*8 phi(ndgphi),psi(ndgphi),phiold(ndgphi),dphidt(ndgphi)
c
      integer i,j,n,np
c
c     stresses
c     --------
      n = indgr
      do j=1,ncpelas
         if ( itordr.eq.2 ) then
            do i=1,ndgelas
               selas(j,i) = z(n+i)
               desdt(j,i) = (selas(j,i)-selasd(j,i))/dt
            enddo
         else
            do i=1,ndgelas
               selas(j,i) = z(n+i)
               desdt(j,i) = z(n+i)*tsc+zs(n+i)
            enddo
         endif
         n = n + ndgelas
      enddo
c
c     velocities
c     ----------
      n = indgu
      do j=1,ncpvelo
         if ( itordr.eq.2 ) then
            do i=1,ndgvelo
               u(j,i) = z(n+i)
               dudt(j,i) = (u(j,i)-uold(j,i))/dt
            enddo
         else
            do i=1,ndgvelo
               u(j,i) = z(n+i)
               dudt(j,i) = z(n+i)*tsc+zs(n+i)
            enddo
         endif
         n = n + ndgvelo
      enddo
c
c     pressure
c     --------
      do j=1,ncppres
      n = indgp
      do i=1,ndgpres
         p(i) = z(n+i)
      enddo
      enddo
c
c     particle velocity
c     -----------------
      n = indgnb
      do np=1,nbrigid
         if ( itordr.eq.2 ) then
            do i=1,ncpb
               uvom(i,np) = z(n+i)
               duvot(i,np) = (uvom(i,np)-uvod1(i,np))/dt
            enddo
         else
            do i=1,ncpb
               uvom(i,np) = z(n+i)
               duvot(i,np) = z(n+i)*tsc+zs(n+i)
            enddo
         endif
         n = n + ncpb
      enddo
c
c     stream function
c     ---------------
      n   = indgstrm
      do i=1,ndgstrm
         strm(i) = z(n+i)
      enddo
c
c     phase function and chemical potential
c     ------------------------------------
      if(ncpphi.eq.1)then
      phi(1:ndgphi)=z(indgphi+1:indgphi+ndgphi)
      psi(1:ndgphi)=z(indgpsi+1:indgpsi+ndgphi)
      n = indgphi
      if(itordr.eq.2)then
         dphidt(1:ndgphi)=(phi(1:ndgphi)-phiold(1:ndgphi))/dt
      else
         dphidt(1:ndgphi)=z(n+1:n+ndgphi)*tsc+zs(n+1:n+ndgphi)
      endif
      endif
      return
      end
