c***********************************************************************
      subroutine boundpre(kaxi,nic,nbd,ic,ibdnod,nbound,nside,
     &                    nrdgeom,ndggeom,x,y,nrdvelo,ndgvelo,indgu,
     &                    nvar,iperm,imp,neq,arhs,
     &                    ibdn,densn,shg1,shg1x,we1)
     
c     This subroutine implements the normal stress boundary condition
c     
c     periodic bc not considered yet! Highest order of interpolation=2
c
c     input:
c        ibdn(nbd): to be used in conjunction with ibpd=0
c              if ibdn=1/0, Y/N stress boundary condition
c        densn(nbd): normal stress on the boundary. densn<0 if p>0.   
c
c     Dec 14, 2007, Pengtao Yue
c        adopted from rhsbblfs in the foaming code     
c     Note:  Only  for the N-S solver, the i in imp(i) is the order of 
c        elimination! In other solvers, i is the global order.
c***********************************************************************
      implicit none
      integer kaxi,nic,nbd,nbound,nrdgeom,nrdvelo,ndggeom,ndgvelo,indgu,
     &        nvar,neq
      integer ic(nic+1),nside(nbound),ibdnod(nbd),iperm(nvar),imp(nvar),
     &        ibdn(nbd)
      real(8) arhs(neq),densn(nbd),x(ndggeom),y(ndggeom)
      real(8) shg1(3,3,2),shg1x(3,3,2),we1(3,2)
c     local variables
      integer i,i1,i2,k1,k2,n,l,ndlg,ndlv,nb,ne,nel
      real(8) xn(3),yn(3),pn(3),pu(3),rhu(3),rhv(3),
     &        pl,xxil,yxil,tmp,tmp1,tmp2,xl
      integer iu,iv,nu,nv,in(3),sec
      logical stressb
c
      ndlg=nrdgeom+1
      ndlv=nrdvelo+1
      if(ndlg.gt.3.or.ndlv.gt.3)then
         write(*,*)'ndl>3 in subroutine rhsbbl'
         stop
      endif
c
      sec=1         
      do nb=1,nbound
         i1=ic(sec)
         sec=sec+nside(nb)
         i2=ic(sec)
         nel=(i2-i1)/nrdgeom
         do ne=1,nel
c     calculate the local nodal values
            k1=i1+nrdgeom*(ne-1)
            if(ibdn(k1)/=1)cycle
            n=ibdnod(k1)
            xn(1)=x(n)
            yn(1)=y(n)
            pn(1)=densn(k1)
            in(1)=n
            stressb=.false.
            do i=2,ndlg
               k2=k1+i-1
               if(k2.eq.i2)then
                  k2=i1
               endif
               if(ibdn(k2)/=1)then
                  stressb=.true.
                  exit
               endif
               n=ibdnod(k2)
               xn(i)=x(n)
               yn(i)=y(n)
               pn(i)=densn(k2)
               in(i)=n
            enddo
            if(stressb)cycle
c     loop over integration points
            rhu(1:ndlv)=0.
            rhv(1:ndlv)=0.
            do l=1,3
               pu(1:ndlv)=shg1(1:ndlv,l,nrdvelo)
               pl=sum(pn(1:ndlg)*shg1(1:ndlg,l,nrdgeom))
               xxil=sum(xn(1:ndlg)*shg1x(1:ndlg,l,nrdgeom))
               yxil=sum(yn(1:ndlg)*shg1x(1:ndlg,l,nrdgeom))
               xl  =sum(xn(1:ndlg)*shg1 (1:ndlg,l,nrdgeom))

c
               tmp =pl*we1(l,nrdvelo)*(kaxi*(xl-1.)+1.)
               tmp1=tmp*xxil
	         tmp2=tmp*yxil
               do i=1,ndlv
                  rhu(i)=rhu(i)+tmp2*pu(i)
                  rhv(i)=rhv(i)-tmp1*pu(i)
               enddo
            enddo
c     assemble the global RHS
            if(ndlv.lt.ndlg) in(2)=in(3)
c      write(*,*)'ndlv',ndlv,ndlg
            do i=1,ndlv
               iu = indgu + in(i)
               iv = indgu + ndgvelo + in(i)
               nu = iperm(iu)
               nv = iperm(iv)
c               write(*,*)iu,nu,imp(iu)
               if(imp(nu).eq.0) arhs(nu)=arhs(nu)+rhu(i)
               if(imp(nv).eq.0) arhs(nv)=arhs(nv)+rhv(i)
            enddo
         enddo
      enddo            
      end
