cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mshreorder(nelem,ngmx,ncmx,inod,nec,iwork,maxiwk)
c-----------------------------------------------------------------------
c   This subroutine reorders the elements using Cuthill & McKee's 
c       algorithm
c
c   input:      
c       nelem:      number of elements
c       ngmx:       number of nodes in each element
c       mcmx:       number of edges (neighboring elm.) in each element
c       inod:       element description table (old)
c       nec:        neighboring element table (old)
c       maxiwk:     maximum size of iwork (>=(2+ncmx+ngmx)*nelem)
c   output:
c       inod,nec    new inod and nec
c   working arrays:
c       iwork    
c   
c   1/6/2015, P. Yue
c-----------------------------------------------------------------------
      implicit none
      intent(in)::  nelem,ngmx,ncmx
      intent(inout):: inod,nec
      integer nelem,ngmx,ncmx,maxiwk
      integer inod(ngmx,nelem),nec(ncmx,nelem),iwork(maxiwk)
c   local variables
      integer lelmold, lelmnew, linod, lnec, leniwk
      lelmold=1
      lelmnew=lelmold+nelem
      linod=lelmnew+nelem
      lnec=linod+ngmx*nelem
      leniwk=lnec+ncmx*nelem-1
      if(leniwk>maxiwk) then
        write(*,*)'gmshreorder, maxiwk should >=',leniwk
        stop
      endif
      call mshnumnel(nelem,ncmx,nec,iwork(lelmold),iwork(lelmnew))
      call mshrrdr(nelem,ngmx,ncmx,inod,nec,iwork(lelmnew),
     &             iwork(linod),iwork(lnec))
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mshrrdr(nelem,ngmx,ncmx,inod,nec,elmnew,inodd,necd)
c-----------------------------------------------------------------------
c   input:
c       inod:   old element description table
c       nec:    old neighboring element table
c       elmnew: new element index
c   ouput:
c       inod,nec:   new arrays
c   working arrays:
c       inodd,necd: to store the old inod and nec
c   
c   1/6/2015,   P. Yue
c-----------------------------------------------------------------------
      implicit none
      intent(in)::  nelem,ngmx,ncmx,elmnew
      intent(inout):: inod, nec
      integer nelem,ngmx,ncmx
      integer inod(ngmx,nelem),nec(ncmx,nelem),elmnew(nelem),
     &        inodd(ngmx,nelem),necd(ncmx,nelem)
c   local variables
      integer n,ne,i
c      
      inodd=inod
      necd=nec
      do n=1,nelem
        ne=elmnew(n)
        inod(1:ngmx,ne)=inodd(1:ngmx,n)
        do i=1,ncmx
          if(necd(i,n)<=0)then
            nec(i,ne)=necd(i,n)
          else
            nec(i,ne)=elmnew(necd(i,n))
          endif
        enddo
      enddo
      end subroutine
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mshnumnel (nelem,ncmx,nec,rmhnel,rmhlen)
c-----------------------------------------------------------------------      
c     The routine renumbers the mesh elements using Cuthill & McKee's 
c       algorithm
c
c     OUTPUT 
c       rmhnel(new elem numb) = old elem numb
c       rmhlen(old elem numb) = new elem numb
c     Modified by Pengtao Yue on Mar 8,2005
c        from Howard's subroutine RNBNEL
c     reused for gmsh, P. Yue, 1/6/2015
c     NOTE: nec<=0 for bounary edges
c***********************************************************************
      implicit none
      intent  (in):: nelem,ncmx,nec
      intent (out):: rmhnel, rmhlen
      integer nelem,ncmx,nec(ncmx,nelem),rmhnel(nelem),rmhlen(nelem)
c
      integer i,ie,nact1,nact2,mact,ne,nne
c
      rmhnel(2:nelem) = 0
      rmhlen(2:nelem) = 0
      rmhnel(1) = 1
      rmhlen(1) = 1
c
      nact1 = 1
      nact2 = 1
      do while (nact2.lt.nelem )
c
         mact = 0
         do ie=nact1,nact2
            ne = rmhnel(ie)
            do i=1,3
               nne = nec(i,ne)
               if ( nne.gt.0 ) then
                  if ( rmhlen(nne).eq.0 ) then
                     mact = mact + 1
                     rmhnel(nact2+mact) = nne
                     rmhlen(nne) = nact2+mact
                  endif
               endif
            enddo
         enddo
c
         nact1 = nact2 + 1
         nact2 = nact2 + mact
c
      enddo
c
      return
      end
