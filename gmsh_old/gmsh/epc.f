cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epc(nvert, nnode, nelem, ngmx, maxtime,
     &               inod, phi, x, y, dist,
     &               iwork,rwork,maxiwk,maxrwk)
c-----------------------------------------------------------------------
c   This subroutine generates a signed distance function (dist) based 
c   on the zero level set of phi. The explicit positive coefficient (EPC)
c   is used.
c   Ref: Barth & Sethian, JCP, 145, 1-40, 1998
c   
c   input:
c       nvert:  number of vertices
c       nnode:  number of nodes
c       nelem:  number of elements
c       ngmx:   number of nodes in each element
c       maxtime:computatinoal time (usually the domain size)
c       maxiwk: length of the integer working array, >=nelem+nvert
c       maxrwk: length of the real(8) working array, >=3*nvert
c       inod:   element description table
c       phi:    field function 
c       x,y:    coordinates of nodes
c       dist:   signed distance function, measured from the zero level 
c               set of phi
c   working arrays:
c       iwork, rwork:   integer and real(8) working arrays
c
c   12/2014, P. Yue     
c-----------------------------------------------------------------------

      implicit none
      intent(in):: nvert, nnode, nelem, ngmx,  phi, x, y, maxtime
      integer nvert, nnode,nelem, ngmx, maxiwk, maxrwk
      real(8) maxtime, phi(nnode),x(nnode),y(nnode),dist(nvert)
      integer inod(ngmx,nelem)
      integer iwork(maxiwk)
      real(8) rwork(maxrwk)
c   local variables
      integer leflag, lvflag, leniwk, lenrwk, lfstar, lwstar, ldist2
      leflag=1
      lvflag=leflag+nelem
      leniwk=lvflag+nvert-1 !nelem+nvert
      lfstar=1
      lwstar=lfstar+nvert
      ldist2=lwstar+nvert
      lenrwk=ldist2+nvert-1 !3*nvert
      if(leniwk>maxiwk)then
        print*, 'subroutine epc, leniwk>maxiwk, increase size of iwork!'
        stop
      endif
      if(lenrwk>maxrwk)then
        print*, 'subroutine epc, leniwk>maxiwk, increase size of rwork!'
        stop
      endif
c   label the elements and vertices, and compute the distance function
c   of interfacial vertices.
      call epc_label(nvert,nnode,nelem, ngmx,
     &               inod,phi,x,y,iwork(leflag),iwork(lvflag),dist)
      call epc_smooth(nvert,nnode, nelem,ngmx,
     &                inod, x, y, iwork(leflag), iwork(lvflag),
     &                dist,rwork(lfstar),rwork(lwstar))     
      call epc_update(nvert,nnode, nelem,ngmx,maxtime,
     &                inod, x, y, iwork(leflag), iwork(lvflag),
     &                dist,rwork(lfstar),rwork(lwstar),rwork(ldist2))
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epc_update(nvert,nnode, nelem,ngmx,maxtime,
     &                      inod, x, y, eflag, vflag, dist,fstar,wstar,
     &                      dist2)
c-----------------------------------------------------------------------
c   This subroutine deals with the time integration in EPC
c   
c   input:
c       nvert, nnode, nelem, ngmx:  mesh parameters
c       maxtime:    computational time
c       inod:   element description table
c       x,y:    coordinates of mesh nodes
c       eflag:  element flag, see epc_label
c       vflag:  vertex flag, see epc_label
c   in/out:   
c       dist:   only the zero level set of the input is meaningful. At 
c               output, dist is a signed distance function
c   working arrays:
c       fstar,wstar,dist2   
c
c   12/2014, P. Yue
c-----------------------------------------------------------------------
      implicit none
      intent(in):: nvert, nelem, ngmx, x, y, eflag, vflag,maxtime
      integer nvert, nnode,nelem, ngmx 
      real(8) x(nnode),y(nnode),dist(nvert),fstar(nvert),wstar(nvert),
     &        dist2(nvert),maxtime
      integer inod(ngmx,nelem),eflag(nelem),vflag(nvert)
c   local variables
      real(8) dt,time
      integer num,it,i

      call epc_timestep(nvert,nnode, nelem,ngmx, inod, x, y,dt)
      num=maxtime/dt
c      print*, 'dt,num=',dt,num
      do it=1,num
      dist2=dist
        call epc_fstar(nvert,nnode, nelem,ngmx,
     &                 inod, x, y, eflag, vflag, dist,fstar,wstar)      
        do i=1,nvert
cc          print*,i,vflag(i),wstar(i),fstar(i)
          if(vflag(i)==1.and.wstar(i)>1.d-10)then
            dist(i)=dist(i)-dt*fstar(i)/wstar(i)
          endif
        enddo
        
        call epc_fstar(nvert,nnode, nelem,ngmx,
     &                 inod, x, y, eflag, vflag, dist,fstar,wstar)      
        do i=1,nvert
ccc          print*,i,vflag(i),wstar(i),fstar(i)
          if(vflag(i)==1.and.wstar(i)>1.d-10)then
            dist(i)=0.5*(dist2(i)+dist(i))-0.5*dt*fstar(i)/wstar(i)
          endif
        enddo
      enddo
      end
                    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epc_fstar(nvert,nnode, nelem,ngmx,
     &                     inod, x, y, eflag, vflag, dist,fstar,wstar)      
c-----------------------------------------------------------------------
c   This subroutine computes fstar abd wstar based on the EPC scheme
c   for level set function.
c   Ref. Barth and Sethian, JCP 145, 1-40, 1998
c   
c   input:
c       nvert   -   total number of vertices
c       nnode   -   total number of nodes
c       nelem   -   total number of elements
c       ngmx    -   maximum number of nodes in each element
c       inod    -   element-node description table
c       x,y     -   coordinates of nodes
c       eflag   -   element flag. 1, phi>0; 0, interfacial element; 
c                                  -1, phi<0.
c       vflag   -   vertex flag. 0, vertices of interfacial elements; 
c                                1, vertices need to be computed by EPC
c       dist    -   distance function to the interface. Only for the 
c                   vertices of interfacial elements
c   output:
c       fstar   -   phi* in B&S's paper
c       wstar   -   w   in B&S's paper
c
c   Revision record
c   12/2014     P. Yue, original code
c-----------------------------------------------------------------------
      implicit none
      intent(in):: nvert, nelem, ngmx, x, y, eflag, vflag, 
     &             dist
      intent(out):: fstar, wstar
      real(8), parameter:: eps=1.d-15
      integer nvert, nnode,nelem, ngmx
      real(8) x(nnode),y(nnode),dist(nvert),fstar(nvert),wstar(nvert)
      integer inod(ngmx,nelem),eflag(nelem),vflag(nvert)
      real(8) hmin
c   external functions
      real(8), external:: epc_norm2
c   local variables
      real(8) xn(3),yn(3),fn(3),area,gradn(2,3),normal(2,3),gradf(2)
      real(8) K(3),Kplus(3),Kminus(3),alpha(3),df(3),dff,ss,d
      integer n,i,nd
      real(8) tmp1,tmp2,tmp3
c
c   d=n for n-dimensional problem      
      d=2.d0
      hmin=1.d10
      fstar(1:nvert)=0.d0
      wstar(1:nvert)=0.d0
      do n=1,nelem
        if(eflag(n)==0)cycle
        xn(1:3)=x(inod(1:3,n))
        yn(1:3)=y(inod(1:3,n))
        fn(1:3)=dist(inod(1:3,n))
        call epc_shape(xn,yn,area,gradn)
        normal(1:2,1:3)=d*area*gradn(1:2,1:3)
c
        gradf(1)=sum(fn(1:3)*gradn(1,1:3))
        gradf(2)=sum(fn(1:3)*gradn(2,1:3))     
        tmp1=epc_norm2(gradf)
c
        if(eflag(n)>0)then 
            ss=1.d0
        else
            ss=-1.d0
        endif
        if(tmp1>eps)then
c
            tmp1=ss/(tmp1*d)
            do i=1,3
                K(i)=sum(gradf(1:2)*normal(1:2,i))*tmp1
                Kplus(i)=max(0.d0,K(i))
                Kminus(i)=min(0.d0,K(i))          
            enddo
c        
            dff=sum(K(1:3)*fn(1:3))
c
            if(abs(dff)>eps)then
                tmp1=min(sum(Kminus(1:3)),-eps)
                do i=1,3
                    tmp2=sum(Kminus(1:3)*(fn(i)-fn(1:3)))  
                    tmp3=Kplus(i)/tmp1*tmp2
                    df(i)=max(0.d0,tmp3/dff)
                enddo
                if(sum(df(1:3))>eps)then
                    alpha(1:3)=df(1:3)/sum(df(1:3))
                else
                    alpha(1:3)=1./3.
                endif
            else
                alpha(1:3)=1./3.
            endif
        else
c       this happens if dist is a constant within the cell. 
            alpha(1:3)=1./3.
        endif
c            
        do i=1,3
          nd=inod(i,n)
          fstar(nd)=fstar(nd)+alpha(i)*(dff-ss*area)
          wstar(nd)=wstar(nd)+alpha(i)*area
        enddo
      enddo
      end  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epc_shape(x,y,area,gradn)
c-----------------------------------------------------------------------
c   This subroutine computes the area and gradient of shape functions of
c   P1 triangular elements
c
c   input:
c       x,y     -   x and y coordinates of the three vertices   
c   output:
c       area    -   element area
c       gradn   -   gradn(i,j) is the i-th component of the gradient of 
c                   the shape function defined at the j-th local vertex
c
c   revision record:
c   11/25/2014  -   P. Yue, original code
c-----------------------------------------------------------------------
      implicit none
      intent(in):: x,y
      intent(out):: area, gradn
      real(8) x(3),y(3),area, gradn(2,3)
c   local variables     
      real(8) x_xi,x_eta,y_xi,y_eta,xi_x,xi_y,eta_x,eta_y,jacobi
c      
      x_xi=x(2)-x(1)
      x_eta=x(3)-x(1)
      y_xi=y(2)-y(1)
      y_eta=y(3)-y(1)
c
      jacobi=x_xi*y_eta-x_eta*y_xi
      area=jacobi/2.d0
c
      xi_x=y_eta/jacobi
      xi_y=-x_eta/jacobi
      eta_x=-y_xi/jacobi
      eta_y=x_xi/jacobi
c     grad(N(j)), where N(j) is the shape function associated 
c     with the j-th node
      gradn(1,1)=-xi_x-eta_x
      gradn(2,1)=-xi_y-eta_y
      gradn(1,2)=xi_x
      gradn(2,2)=xi_y
      gradn(1,3)=eta_x
      gradn(2,3)=eta_y
c
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epc_label(nvert, nnode, nelem, ngmx,
     &                    inod,phi, x, y, eflag, vflag, dist)
c-----------------------------------------------------------------------
c   This subroutine labels the interfacial cells and vertices, where the 
c   interface is defined as the level set phi=0
c   
c   input:
c       nvert   -   total number of vertices
c       nnode   -   total number of nodes
c       nelem   -   total number of elements
c       ngmx    -   maximum number of nodes in each element
c       inod    -   element-node description table
c       phi     -   phase-field varible or level-set function
c       x,y     -   coordinates of nodes
c   output:
c       eflag   -   element flag. 1, phi>0; 0, interfacial element; 
c                                  -1, phi<0.
c       vflag   -   vertex flag. 0, vertices of interfacial elements; 
c                                1, vertices need to be computed by EPC
c       dist    -   distance function to the interface. Only for the 
c                   vertices of interfacial elements
c
c   Revision record:
c   10/28/2014, P. Yue, original code
c-----------------------------------------------------------------------        
      implicit none
      intent(in):: nvert, nnode, nelem, ngmx, phi, x, y
      integer nvert, nnode,nelem, ngmx
      real(8) phi(nnode),x(nnode),y(nnode),dist(nvert)
      integer inod(ngmx,nelem),eflag(nelem),vflag(nvert)
c   working variables
      integer n,i
      real(8) f(3),xn(2,3),un(3)
 
      vflag(1:nvert)=1
      dist(1:nvert)=5.
      
      do n=1,nelem
        f(1:3)=phi(inod(1:3,n))
        if(f(1)>0.and.f(2)>0.and.f(3)>0) then
          eflag(n)=1
        elseif(f(1)<0.and.f(2)<0.and.f(3)<0) then
          eflag(n)=-1
        else
          eflag(n)=0
          vflag(inod(1:3,n))=0
          xn(1,1:3)=x(inod(1:3,n))
          xn(2,1:3)=y(inod(1:3,n))
          call epc_elmdist(f,xn,un)
          do i=1,3
            dist(inod(i,n))=min(dist(inod(i,n)),un(i))
cc            print*, '***', un(i),dist(inod(i,n))
          enddo
        endif
      enddo
c   set the signs of u
      do i=1,nvert
        if(vflag(i)/=0)then
          dist(i)=phi(i)
        elseif (phi(i)<0) then
          dist(i)=-dist(i)
        endif
      enddo
      return
      end     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epc_elmdist(f,xn,un)
c-----------------------------------------------------------------------
c   This sbroutine work on a triangular element. It computes the 
c   distance of each vertex to the interface, which is represented by 
c   the zero level set of f.
c
c   input
c       f   -   f value at vertices (interface is represented by f=0)
c       xn  -   coordinates at vertices.
c   output
c       un  -   distance of each vertex to the interface
c   working (local) variables
c       fn(3)       -   absolution values of f(3)
c       line(2,2)   -   line(1:2,1) and line(1:2,2) are the starting and
c                       end points of the line segment f=0
c       dl(2)       -   vector connecting the two ends of line(2,2)
c
c   Revision record
c   10/31/2014,     P. Yue  orginal code
c-----------------------------------------------------------------------         
      implicit none
      real(8), parameter:: eps=1.d-10
      intent(in):: f,xn
      intent(out):: un
      real(8) f(3),xn(2,3),un(3)
c   working variables      
      real(8) fn(3),line(2,2),dl(2),dx1(2),dx2(2),dldl,t
      integer i,i1,i2,j
      real(8), external:: epc_norm2
      
      fn(1:3)=abs(f(1:3))
cc      print*, f(1:3)
cc      print*, xn(1,1:3)
cc      print*, xn(2,1:3)
cc      pause
c     find the node with fn=0. 
      i1=0
      if(f(1)*f(2)*f(3)==0)then
        do  i=1,3
          if(f(i)==0) then
            i1=i-1
            exit
          endif
        enddo
      endif
c      mi=1      
c      minf=fn(mi)
c      do i=2:3
c        if(fn(i)<minf)then
c          mi=i
c          minf=fn(i)
c        endif
c      enddo
c   Find the ending point from an edge with large fn first.
c   this ensures that if fn=0 at one vertex, that same point will not be
c   used as the starting and ending point simultaneously.
      j=0
      do i=1,3
        i1=mod(i1,3)+1
        i2=mod(i1,3)+1
        if(f(i1)*f(i2)<=0)then
          j=j+1
          if(fn(i1)+fn(i2)<=eps)then
            line(1:2,j)=0.5*(xn(1:2,i1)+xn(1:2,i2))
          else
            line(1:2,j)=(fn(i1)*xn(1:2,i2)+fn(i2)*xn(1:2,i1))
     &                 /(fn(i1)+fn(i2))
          endif
          if(j==2)exit
        endif
      enddo
      if(j/=2)then
        stop 'error in epc_elmdist, j /=2'
      endif 
c
      dl(1:2)=line(1:2,2)-line(1:2,1)
      dldl=sum(dl*dl)
cc      print*,'dldl',dldl
      do i=1,3
        dx1(1:2)=xn(1:2,i)-line(1:2,1)
        dx2(1:2)=xn(1:2,i)-line(1:2,2)
        t=sum(dx1*dl)/dldl
cc        print*, 't',t
        if(t<1.and.t>0)then
          un(i)=epc_norm2(dx1-t*dl)
        else
          un(i)=min(epc_norm2(dx1),epc_norm2(dx2))
        endif
      enddo
cc      print*,'u',un(1:3)
      return
      end  
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epc_smooth(nvert,nnode, nelem,ngmx,
     &                     inod, x, y, eflag, vflag,dist,fstar,wstar)
c-----------------------------------------------------------------------      
c   This subroutine smooths the intial distance function to avoids
c   islands of low |dist| generated by numerical oscillation. 
c
c   For each vertex, the smoothing is obtained by the area weighted 
c   average of elements sharing the same vertex.
c   
c   input:
c       nvert:  number of vertices
c       nnode:  number of nodes
c       nelem:  number of elements
c       ngmx:   number of nodes in each element
c       inod:   element description table
c       x,y:    nodal coordinates
c       eflag:  element flag, see subroutine epc_label
c       vflag:  vertex flag, see subroutine epc_label
c   in/out:
c       dist:   distance function
c   working:
c       fstar,wstar
c
c   12/2014, P. Yue
c-----------------------------------------------------------------------      

      intent(in):: nvert, nelem, ngmx, x, y, eflag, vflag
      intent(inout):: dist
      integer nvert, nnode,nelem, ngmx
      real(8) x(nnode),y(nnode),dist(nvert),fstar(nvert),wstar(nvert)
      integer inod(ngmx,nelem),eflag(nelem),vflag(nvert)
c   local variables
      real(8) xn(3),yn(3),fn,area,gradn(2,3)
      integer n,i,nd,it
c
      do it=1,10
        fstar=0.d0
        wstar=0.d0
        do n=1,nelem
          if(eflag(n)==0)cycle
          xn(1:3)=x(inod(1:3,n))
          yn(1:3)=y(inod(1:3,n))
          call epc_shape(xn,yn,area,gradn)
          fn=sum(dist(inod(1:3,n)))/3.d0
          wstar(inod(1:3,n))=wstar(inod(1:3,n))+area
          fstar(inod(1:3,n))=fstar(inod(1:3,n))+area*fn
        enddo
        do i=1,nvert
          if(vflag(i)/=0) then
            dist(i)=fstar(i)/wstar(i)
          endif
        enddo          
      enddo
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
      subroutine epc_timestep(nvert,nnode, nelem,ngmx, inod, x, y,dt)
c-----------------------------------------------------------------------      
c   This subroutine computes the time step for reinitialization of 
c   distance function.
c       dt=0.25*hmin
c
c   input:
c       nvert, nnode, nelem, ngmx:  mesh parameters
c       inod:   element description table
c       x,y:    nodal coordinates
c   output:
c       dt:     time step
c-----------------------------------------------------------------------      
      intent(in):: nvert, nnode, nelem, ngmx, inod, x, y
      intent(inout):: dt
      integer nvert, nnode,nelem, ngmx
      real(8) x(nnode),y(nnode),dist(nvert),dt
      integer inod(ngmx,nelem)
c   local variables
      real(8) xn(3),yn(3),area,gradn(2,3),hmin
      integer n
c
      hmin=1.d10
      do n=1,nelem
          xn(1:3)=x(inod(1:3,n))
          yn(1:3)=y(inod(1:3,n))
          call epc_shape(xn,yn,area,gradn)
          hmin=min(hmin,sqrt(area))
      enddo
      dt=0.25*hmin
      end       

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      function epc_norm2(x)
c-----------------------------------------------------------------------           
c   This function returns the two norm of a two-component vector      
c-----------------------------------------------------------------------           
      implicit none
      real(8) epc_norm2,x(2)
      epc_norm2=sqrt(x(1)*x(1)+x(2)*x(2))
      return
      end function
      
          
      
