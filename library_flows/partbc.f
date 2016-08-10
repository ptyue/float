c***********************************************************************
      subroutine partbc ( indgnb,neq,nelt,nbrigid,ncpb, z,zs,tsc,
     &                    fxym,fpart,amas, a,ia,arhs,jacoct,ifixpt )
c
c     this routine modifies the global matrix and vector due to moving
c     particles
c
c     ifixpt(3) = .true. the particle velocity (U,V, or Omega) is fixed
c
c     called by : front
c***********************************************************************
      implicit none
      integer nbrigid,neq,nelt,indgnb,ia(neq+1), ncpb
      real*8 fxym(ncpb,nbrigid),fpart(ncpb,nbrigid),amas(ncpb,nbrigid),
     &       a(nelt),arhs(neq),z(neq),zs(neq),tsc
      logical jacoct,ifixpt(3)
c
      integer np,nv,i,ir
c
      do np=1,nbrigid
        do i=1,ncpb
          nv = indgnb + ncpb*(np-1) + i
          fxym(i,np) = fpart(i,np) + arhs(nv)
          arhs(nv) = arhs(nv)+fpart(i,np)-amas(i,np)*(z(nv)*tsc+zs(nv))
          if ( ifixpt(i) ) arhs(nv) = 0.0
          if ( .not.jacoct ) then
             ir = ia(nv)
             a(ir) = a(ir) + amas(i,np)*tsc
             if ( ifixpt(i) ) a(ir) = 1.d50
          endif
        enddo
      enddo
c
      return
      end

