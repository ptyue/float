c***********************************************************************
      subroutine rstextra (ndgelas,ndgphi,ncpelas,ncpphi,ncps,
     &                     selas,desdt,phi)
c
c     This routine reset the extra stress to zero in the Newtonian 
c        region
c     Note:
c        the order of basis function for phi must be equal to or larger
c        than that of extra stres (i.e. nrdelas<=nrdphi or ndgelas<=
c        ndgphi, which is always the case)
c     Pengtao Yue, May 13, 2005
c***********************************************************************
      implicit none
      intent(in)   :: ndgelas,ndgphi,ncpelas,ncpphi,ncps,phi
      intent(inout):: selas,desdt
      integer ncpelas,ncpphi,ndgelas,ndgphi,ncps
      real(8) selas(ncps,ndgelas),desdt(ncps,ndgelas),phi(ndgphi)
c
      integer i      
      real(8),parameter:: gate1=-0.99,gate2=-0.90
      real(8),parameter:: delta=gate2-gate1
      real(8) f
      if(ncpphi.ne.1.or.ncpelas.eq.0)return
      if(ndgphi.lt.ndgelas)then
         write(*,*)'Error in rstextra: ndgphi<ndgelas!'
         stop
      endif
      do i=1,ndgelas
         f=phi(i)
         if(f.lt.gate1)then
            selas(1:ncpelas,i)=0.d0
            desdt(1:ncpelas,i)=0.d0
         elseif(f.lt.gate2)then
            f=(f-gate1)/delta
            selas(1:ncpelas,i)=selas(1:ncpelas,i)*f
            desdt(1:ncpelas,i)=desdt(1:ncpelas,i)*f
         endif
      enddo
      write(*,'(A,F8.4)')'extra stress reset at points where phi<',gate2
      end         
                        