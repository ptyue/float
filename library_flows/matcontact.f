c***********************************************************************
      subroutine contact(kaxi,nic,nbd,ic,nbound,ibdnod,nside,
     &                   nrdgeom,nrdvelo,nrdphi,ndggeom,ndgvelo,ndgphi,
     &                   indgu,indgphi,x,y,
     &                   neq,nvar,nelt,z,zs,a,ia,ja,iperm,arhs,
     &                   wallrela,wallener,tsc,phlam,pheps,
     &                   shg1,shg1x,we1,jwk,jacoct)
c     This subroutine modifies the b in Ax=b to account for stress 
c     at the bubble surface and free surface
c     
c     periodic bc not considered yet! Highest order of interpolation=2:
c        nrdgveo,nrdphi <=nrdgeom<=2
c
c     in/output:  a(nelt), arhs(neq)    
c        neq:  number of unknowns
c        nvar: number of variables (including boundaries)
c
c     Working array: jwk(neq)
c
c     Aug,06 2006, Pengtao Yue
c        adopted from rhsbbl, the old rhsbbl is obselete
c     Note:  Only  for the N-S solver, the i in imp(i) is the order of 
c        elimination! In other solvers, i is the global order.
c***********************************************************************
      implicit none
      integer kaxi,nic,nbd,nbound,neq,nvar,nelt,
     &        nrdphi,nrdgeom,nrdvelo,ndgphi,ndggeom,ndgvelo,
     &        indgphi,indgu
      integer ic(nic+1),nside(nbound),ibdnod(nbd),iperm(nvar),
     &        ia(neq+1),ja(nelt),jwk(neq)
      real(8) arhs(neq),x(ndggeom),y(ndggeom),z(nvar),zs(nvar),a(nelt)
      real(8) shg1(3,3,2),shg1x(3,3,2),we1(3,2)
      real(8) phlam,pheps,wallrela(nic),wallener(nic),tsc
      logical jacoct
c     local variables
      integer i,j,i1,i2,ik,jk,n,l,k,k1,k2,sec,sec1,sec2,
     &        nb,ne,nel,ndlg,ndlv,ndlphi,ndlmax,in(3),irowst,ilast
      real(8) xn(3),yn(3),un(3),vn(3),fn(3),fsn(3),pf(3),
     &        dpfds(3),rh(3),st(3,3),
     &        dxxi,dyxi,dsxi,dxis,xl,ul,vl,vel,fl,fsl,wel,dfdsl,epslam,
     &        wlre,wlen
      real(8) tmp
c
      ndlg=nrdgeom+1
      ndlv=nrdvelo+1
      ndlphi=nrdphi+1
      ndlmax=max(ndlg,ndlv,ndlphi)
c
      epslam=pheps**2/phlam
c
      if(ndlmax>3)then
         write(*,*)'ndl>3 in subroutine rhsbbl'
         stop
      endif
c
      sec1=1         
      do nb=1,nbound
        sec2=sec1+nside(nb)
        i1=ic(sec1)
        i2=ic(sec2)
        do sec=sec1,sec2-1
         nel=(ic(sec+1)-ic(sec))/nrdgeom
         wlre=wallrela(sec)
         wlen=wallener(sec)
         do ne=1,nel
c     calculate the local nodal values
            k1=ic(sec)+nrdgeom*(ne-1)
            n=ibdnod(k1)
            xn(1)=x(n)
            yn(1)=y(n)
            un(1)=z(indgu+n)
            vn(1)=z(indgu+n+ndgvelo)
            fn(1)=z(indgphi+n)
            fsn(1)=zs(indgphi+n)
            in(1)=n
            do i=2,ndlmax
               k2=k1+i-1
               if(k2.eq.i2)then
                  k2=i1
               endif
               n=ibdnod(k2)
               if(n<=ndggeom)then
                  xn(i)=x(n)
                  yn(i)=y(n)
               endif
               if(n<=ndgvelo)then
                  un(i)=z(indgu+n)
                  vn(i)=z(indgu+n+ndgvelo)
               endif
               if(n<=ndgphi)then
                  fn(i)=z(indgphi+n)
                  fsn(i)=zs(indgphi+n)
               endif
               in(i)=n
            enddo
c
            if(nrdgeom<ndlmax)then
               xn(ndlg)=xn(ndlmax)
               yn(ndlg)=yn(ndlmax)
            endif
            if(nrdvelo<ndlmax)then
               un(ndlv)=un(ndlmax)
               vn(ndlv)=vn(ndlmax)
            endif
            if(nrdphi<ndlmax)then
                fn(ndlphi)= fn(ndlmax)
               fsn(ndlphi)=fsn(ndlmax)
                in(ndlphi)= in(ndlmax)
            endif
c
c     loop over integration points
            rh(1:ndlphi)=0.
            st(1:ndlphi,1:ndlphi)=0.
            do l=1,3
c
               dxxi=sum(xn(1:ndlg) *shg1x(1:ndlg,l,nrdgeom))
               dyxi=sum(yn(1:ndlg) *shg1x(1:ndlg,l,nrdgeom))
               dsxi=sqrt(dxxi**2+dyxi**2)
               dxis=1./dsxi
c
               pf(1:ndlphi) = shg1(1:ndlphi,l,nrdphi)
               dpfds(1:ndlphi) = shg1x(1:ndlphi,l,nrdphi)*dxis
c
               ul    = sum(un(1:ndlv)* shg1(1:ndlv,l,nrdvelo))
               vl    = sum(vn(1:ndlv)* shg1(1:ndlv,l,nrdvelo))
               xl    = sum(xn(1:ndlg)* shg1(1:ndlg,l,nrdgeom))
               fl    = sum( fn(1:ndlphi)* pf  (1:ndlphi))
               fsl   = sum(fsn(1:ndlphi)* pf  (1:ndlphi))
               dfdsl = sum( fn(1:ndlphi)*dpfds(1:ndlphi))
               vel   = (ul*dxxi+vl*dyxi)/dsxi
c
               wel=we1(l,nrdphi)*epslam*dsxi
               if(kaxi==1)wel=wel*xl
c               
               tmp=wel*(wlre*(tsc*fl+fsl+vel*dfdsl)
     &                 +0.75d0*wlen*(1.d0-fl**2))
               do i=1,ndlphi
                  rh(i)=rh(i)+tmp*pf(i)
               enddo
               if(jacoct)cycle
               do j=1,ndlphi
                  tmp=wel*(wlre*(tsc*pf(j)+vel*dpfds(j))
     &               -1.5d0*wlen*fl*pf(j))
                  do i=1,ndlphi
                     st(i,j)=st(i,j)-tmp*pf(i)
                  enddo
               enddo
            enddo
c     assemble the global RHS and Matrix
            do i=1,ndlphi
               ik=iperm(indgphi+ndgphi+in(i))
               arhs(ik)=arhs(ik)+rh(i)
c
               if(jacoct)cycle
               irowst = ia(ik)
               ilast  = ia(ik+1)-1
               do k=irowst,ilast
                  jwk(ja(k)) = k
               enddo
c
               do j=1,ndlphi
                  jk=iperm(indgphi+in(j))
                  k=jwk(jk)
                  a(k)=a(k)+st(i,j)
               enddo
            enddo
c
         enddo
        enddo
        sec1=sec2
      enddo            
c     pause
      end

