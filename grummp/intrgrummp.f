c***********************************************************************
      subroutine GRunimsh(mxnbd,mxvert,mxelem,ngmx,ncmx,nbd,nvert,nelem,
     &                    refseg,bdseg,inod,nec,reft,x,y,h)
c
c-----------------------------------------------------------------------
c     This subroutine generate a uniform triangluar mesh using GRUMMP
c
c     input:
c        maximum numbers: (for array bound checking)
c           mxnbd,mxvert,mxelem                 
c        amphi.bdry file specifying the domain
c        h:       mesh size 
c    output:
c        nbd:     number of boundary nodes
c        nvert:   number of vertices
c        nelem:   number of elements(cells)
c        bdseg:   end nodes of boundary segments
c        refseg:  reference number of boundary segments
c        inod:    element description table
c        nec:     neighbouring elements info
c        x,y:     coordinates of vertices.
c        reft:    reference number of elements
c     Note: here we mean the same thing by vertex and node.
c-----------------------------------------------------------------------
c                              n3
c     element info:           /  \
c                          e3/    \e2
c                           /      \
c                          n1------n2
c                              e1
c     nec(i,ne)=0 when the edge is a boundary edge
c     boundary nodes info:
c           --bdnode i------------bdnode j--->
c                      bdseg n         
c           bdseg(1,n)=i,bdseg(2,n)=j
c                
c-----------------------------------------------------------------------
      implicit none
      
      
      intent (in)  :: h,mxnbd,mxvert,mxelem,ngmx,ncmx
      intent (out) :: nbd,nvert,nelem,bdseg,inod,nec,x,y,reft
      integer nbd,nvert,nelem,mxnbd,mxvert,mxelem,ngmx,ncmx
      integer bdseg(2,mxnbd),refseg(mxnbd),
     &        inod(ngmx,mxelem),nec(ncmx,mxelem),reft(mxelem)
      real*8  x(mxvert),y(mxvert),h
c-----------------------------------------------------------------------
      real*8 dG
      dG=100.0D0
      call GRUMMP_INIT(mxnbd,mxvert,mxelem,ngmx,     ncmx,
     &                 nbd,  nvert, nelem, refseg(1),bdseg(1,1),
     &                 inod(1,1),nec(1,1),reft(1),x(1),y(1),h,dG)
c-----------------------------------------------------------------------  
      return    
      end
c***********************************************************************
      subroutine GRrfnmsh(mxnbd,mxvert,mxelem,ngmx,ncmx,nvertd,phid,
     &                    nbd,nvert,nelem,bdseg,refseg,
     &                    inod,nec,reft,x,y,h1,h2,h3,dG,lrmsh)
c
c-----------------------------------------------------------------------
c     This subroutine returns the new refined mesh given an old mesh 
c        with a phi field defined on it.
c
c     input:
c        amphi.bdry: file specifying the domain
c        old mesh info: presaved in data file or memory by GRUMMP
c        maximum numbers: (for array bound checking)
c           mxnbd,mxvert,mxelem                 
c        old mesh related variables: 
c           nvertd:  number of vertices
c           phid:    phi field
c
c        mesh size:
c           h1:   mesh size along the interface(phi=0 contour)
c           h2:   mesh size at phi=1
c           h3:   mesh size at phi=-1
c           dG:   grading control in GRUMMP. Length scale in the mesh
c                 can change by no more than d/G as you move a distance
c                 of d.
c    output:
c        new mesh variables:
c           nbd:     number of boundary nodes
c           nvert:   number of vertices
c           nelem:   number of elements(cells)
c           bdseg:   end nodes of boundary segments
c           refseg:  reference number of boundary segments
c           inod:    element description table
c           nec:     neighbouring elements info
c           x,y:    coordinates of vertices.
c        lrmsh:   whether remesh has been performed
c                 T, remeshed; F, the old mesh is good enough, remesh
c                 not performed.
c-----------------------------------------------------------------------
      implicit none
      
      
      intent (in) :: mxnbd,mxvert,mxelem,ngmx,ncmx,nvertd,phid,h1,h2,h3,
     &               dG
      intent (out)   :: nbd,nvert,nelem,inod,nec,bdseg,refseg,x,y,lrmsh
      integer nvertd,nbd,nvert,nelem,mxnbd,mxelem,mxvert,ngmx,ncmx
      integer bdseg(2,mxnbd),refseg(mxnbd),inod(ngmx,mxelem),
     &        nec(ncmx,mxelem),reft(mxelem)
      real(8)  x(mxvert),y(mxvert),phid(nvertd),h1,h2,h3
      logical lrmsh
      real(8) dG
c
      integer lmsh,lwrmesh
      lwrmesh = 0
c
      call GRUMMP_ADP(mxnbd,     mxvert,    mxelem,    ngmx,    ncmx,
     &                nvertd,    phid(1),   nbd,       nvert,   nelem,
     &                bdseg(1,1),refseg(1), inod(1,1), nec(1,1),reft(1),
     &                x(1),      y(1),      h1,    h2, h3,     lmsh,
     &                dG,lwrmesh) 
c
      if(lmsh.eq.1)then
         lrmsh=.true.
      else
         lrmsh=.false.
      endif
c      
      end
      
C      
      subroutine WRAMPHIMESH()
      implicit none
      integer  :: lwrmesh
     	lwrmesh=1
     	
      call GRUMMP_ADP(0, 0, 0, 0, 0, 
     & 		    0, 0, 0, 0, 0,
     & 		    0, 0, 0, 0, 0,     	
     & 		    0, 0, 0, 0, 0, 0,     	
     & 		    0, lwrmesh)
          	
      end subroutine 
