c***********************************************************************
      subroutine mesh2d (coor,coorcp,h,nu,nbs,nbsmx,bdseg,nba,sd,nbsd,
     &                   reft,nbt,regul,h0, iopt,err,tri )
c-----------------------------------------------------------------------
c     purpose : to mesh a domain given information of its boundaries
c     Note : This program uses a double precision working array for the 
c            coordinates and works for periodic mesh too!
c-----------------------------------------------------------------------
c     input :
c     ------
c           coor(2,nbsmx)  coordinates of the boundary vertices
c                          mapped coordinates for periodic geometry
c           coorcp(2,nbs)  copy of the boundary nodes 
c                          original coordinates for periodic geometry
c           h(nbsmx)     mesh spacing near each boundary vertex
c           nbs          number of boundary vertices
c           nbsmx        maximum number of vertices 
c
c           bdseg(2,nba) specification of the boundary segments 
c                        = vertex # of the i-th endpoint (i=1,2)
c           nba          number of boundary segments
c
c           sd(2,nbsd)   specification of the subdomains to be meshed
c                          sd(1,i) = +e, or sd(1,i) = -e
c                          where i is the subdomain number
c                                e is the segment # of ANY segment which
c                                  which borders the subdomain
c                                +e if the orientation of the seg. is the
c                                   same as that of the boundary
c                                -e otherwise
c                          sd(2,i) is the REFERENCE number of the subdomain
c           nbsd         number of subdomains to be meshed
c
c           regul        if regul = true, the mesh will be smoothed
c           h0           background mesh size
c           iopt         option index, 
c                        = 0 normal call
c                        = 1 form a mesh connecting boundary nodes only
c                        = 2 refine an existing mesh
c
c        working arrays :
c        ----------------
c           tri(ltri)    integer array
c        
c        output :
c        --------
c         nbs         total number of vertices (boundary + interior)
c         nbt         number of triangles
c         coor(1:2,nbs) coordinates of the vertices (boundary + interior)
c         nu(1:3,nbt) node number of three vertices for each triangle
c                     (in counterclockwise order)
c         nu(4:6,nbt) neighboring element number on each side of a triangle
c                     = 8*t+a (t=element number, a=side number)
c                     = -(8*t+a) interface element
c                     = 0 neighboring a boundry segment
c         reft(1:nbt) reference number for each triangle
c
c         err    err = 0    => normal return
c                else       => nbt = 0 and bug in data 
c         
c     definitions of the parameters
c     -----------------------------
c     nbtmx = 2*(nbs-1)
c     ltri  = max(3*nbtmx+nbs+1,4*nbs+2*nbsd,nba)
c
c     dimensions of the arrays
c     ------------------------
c     integer : nu(6*nbtmx)  , reft(nbtmx) , tri(ltri) 
c     integer : bdseg(2,nba) , sd(2,nbsd)
c     double precision  : coor(2*nbsmx)  , h(nbsmx)
c
c----------------------------------------------------------------------
c   programmer :  H.H. Hu
c                 based on mshptg.f by F. Hecht, INRIA
c-----------------------------------------------------------------------
c    
      integer nbs,nbt,nba,nbsmx,nbsd
      integer tri(7*nbsmx-5),bdseg(2,nba),sd(2,nbsd)
      integer nu(6,2*(nbsmx-1)),reft(2*(nbsmx-1))
      double precision coor(2,nbsmx),coorcp(2,nbs),h(nbsmx),h0,puis
      integer err,iopt
c
c     for parameters of the regularization
      logical regul
      integer itermx
      double precision omega,eps
      parameter (itermx=20, omega = 1.0, eps = .005)
c
c     puis         controls the rate of change of the mesh size during
c                  the generation of interior vertices
c                  .25 => a good value
c                  .1  => gives more vertices
      parameter (puis=0.25 )
c
c     for local variables
      integer i,nbtgrn,nbsgrn
c----------------------------------------------------------------------
      err = 0
      if ( nbs.lt.3 .or. nbsmx.lt.nbs ) then
        err = 1
        print *,'fatal error mshpts : number of points ',nbs
     &         ,' is < 3 or > ',nbsmx,' maximum '
        return
      endif
c
      if ( iopt.eq.2 ) go to 40
c
c     1. input data preparation 
c     -------------------------
      call bdinit (coor,nbs,tri,h(nbs+1),err)
      if ( err.ne.0 ) return
c
c     2. meshing the boundary nodes
c     -----------------------------
      call mshbdn (coor,nu,tri,nbs,err)
      if ( err.ne.0 ) return
c   
c     3. reset the coordinates
c     ------------------------
      do i=1,nbs
         coor(1,i) = coorcp(1,i)
         coor(2,i) = coorcp(2,i)
      enddo
c
c     4. treatment of the front
c     -------------------------
      call frtref(coor,nu,nbs,nbt,bdseg,nba,sd,nbsd,reft,tri,err)
      if ( err.ne.0 ) return
c
c     edge swappig for the transformed mesh
c     -------------------------------------
      do i=1,nbt
         call diagswp (coor,nu,i,4,nbs,err)
         call diagswp (coor,nu,i,5,nbs,err)
         call diagswp (coor,nu,i,6,nbs,err)
      enddo
      if ( err.ne.0 ) return
c
c     5. generate interior nodes
c     --------------------------
40    nbsgrn = nbs
      nbtgrn = nbt
      if ( iopt.eq.1 ) goto 50
      call mshrfn (coor,nu,h,h0,reft,nbsgrn,nbsmx,nbtgrn,puis,err)
      if ( err.ne.0 ) goto 60
c
c     6. mesh smoothing
c     -----------------
50    if ( iopt.ne.1 .and. regul ) then
        call mshvoi (nu,tri,tri(nbsgrn+2),nbtgrn,nbsgrn)
        call mshsmth (coor,nbs,nbsgrn,nu,tri,tri(nbsgrn+2),nbtgrn,
     &                omega,itermx,eps)
      endif
c
60    nbt = nbtgrn
      nbs = nbsgrn
c
      return
      end
c
c**********************************************************************
      subroutine bdinit (c,nbs,tri,pat,err)
c
c     rearrangment of the supplied bd. nodes
c----------------------------------------------------------------------
      implicit none
      integer nbs,tri(nbs),err
      double precision c(2,nbs),trfri(4)
      integer iii,ic,ip,i,j,jp,k,trik,tri3,ierr
      double precision aa1,aa2,xmin,xmax,ymin,ymax
      double precision xx,det,pat(nbs)
      double precision precis
      parameter (precis=1.d0)
c   
c     1. initialization (needed, e.g. symmetric domain)
c     ------------------------------------------------- 
      ierr = 0
      iii = 1
      xmin = c(1,1)
      ymin = c(2,1)
      xmax = c(1,1)
      ymax = c(2,1)
      do ic=1,nbs
         xmin = min(c(1,ic),xmin)
         ymin = min(c(2,ic),ymin)
         xmax = max(c(1,ic),xmax)
         ymax = max(c(2,ic),ymax)
         if ( c(1,ic).lt.c(1,iii) ) iii=ic
      enddo
      aa1 = precis/(xmax-xmin)
      aa2 = precis/(ymax-ymin)
      aa1 = min(aa1,aa2)
      aa2 = aa1*(c(2,iii)-ymin)
      trfri(1) = 1.d0/aa1
      trfri(2) = c(1,iii)
      trfri(3) = c(2,iii)
      trfri(4) = aa2
      do ic=1,nbs
         tri(ic) = ic
         c(1,ic) = aa1*(c(1,ic)-trfri(2))
         c(2,ic) = aa1*(c(2,ic)-trfri(3))
         pat(ic)= c(1,ic)**2 + c(2,ic)**2
      enddo
c
c     2. rearrangment of the points
c     -----------------------------
      call heapor(pat,tri,nbs)
      ip = 1
      xx = pat(ip)
      do jp=1,nbs
        if ( pat(jp).gt.xx ) then
          call heapor (pat(ip),tri(ip),jp-ip)
          do i=ip,jp-2
            if ( pat(i).eq.pat(i+1) ) then
              ierr = ierr + 1
              print *,' error, points ',tri(i),tri(i+1),' are same'
            endif
          enddo
          xx = pat(jp)
          ip = jp
        endif 
        ic = tri(jp)
        pat(jp) = c(2,ic)
      enddo
      call heapor (pat(ip),tri(ip),nbs-ip)
      do i=ip,jp-2
        if ( pat(i).eq.pat(i+1) ) then
          ierr = ierr + 1
          print *,' error, points ',tri(i),tri(i+1),' are same'
        endif
      enddo
      if ( ierr.ne.0 ) stop ' fatal error bdinit: points collide'
c
c     3. find the third node for the first element
c     --------------------------------------------
      do k=3,nbs
        det = c(1,tri(2))*c(2,tri(k)) - c(2,tri(2))*c(1,tri(k))
        if ( abs(det).gt.1.d-15 ) then
c         k is the first point not on the line
          trik = tri(k)
          do j=k-1,3,-1
            tri(j+1) = tri(j)
          enddo
          tri(3) = trik
c
          if ( det.lt.0 ) then
c            reverse node 2 and 3 (1-2-3: counterclockwise)
             tri3   = tri(3)
             tri(3) = tri(2)
             tri(2) = tri3
          endif
          return
        endif
      enddo
c
      print *,'fatal error bdinit, all points are on along a line'
      print *,'tri =',(tri(k),k=1,nbs)
c
      stop
      end
c
c**********************************************************************
      subroutine heapor (ra,ia,n)
c
c     Sorts arrray ra(n) into ascending order using Heapsort algorithm
c     reference:  "Numerical Recipes" by W.H.Press et.al.
c
      implicit none
      integer i,l,ir,j,n
      integer ia(n)
      double precision ra(n)
      double precision rra
      integer iia
c     
      if ( n.le.1 ) return
      l = n/2+1
      ir = n
2     if ( l.gt.1 ) then
         l = l-1
         iia = ia(l)
         rra = ra(l)
      else
         iia = ia(ir)
         rra = ra(ir)
         ia(ir) = ia(1)
         ra(ir) = ra(1)
         ir = ir - 1
         if ( ir.eq.1 ) then
            ia(1) = iia
            ra(1) = rra
            return
         endif
      endif
      i = l
      j = l+l
 10   if ( j.le.ir ) then
         if ( j.lt.ir ) then
            if ( ra(j).lt.ra(j+1) ) j=j+1
         endif
         if ( rra.lt.ra(j) ) then
            ia(i) = ia(j)
            ra(i) = ra(j)
            i = j
            j = j+j
         else
            j = ir+1
         endif
         go to 10
      endif
      ia(i) = iia
      ra(i) = rra
      go to 2
      end
c
c**********************************************************************
      subroutine mshbdn (c,nu,tri,nbs,err)
c
c     generation of the initial mesh with supplied bd nodes
c       using front method
c 
      implicit none
      integer nbs,nu(6,2*nbs-2),tri(nbs),tete
      double precision c(2,nbs)
      integer frntnd,err
      integer i,s,t,pf,ppf,psf,npf,pp,ps,taf,iaf,free,ttaf
      parameter (pp=3,ps=4)
c
c------------------------------------------------------------------------
c     definition of nu(1:6,2*nbs-2)
c     if ( nu(5:6,ie) = (0,0) ) then 
c        ie = front node
c        nu(1,ie) = node number
c        nu(2,ie) = 8*t + a
c                   t = triangle element number
c                   a = segment starting with node nu(1,ie)
c        nu(pp,ie) = pointer to the previous front node
c        nu(ps,ie) = pointer to the subsequent front node
c     else 
c        ie = triangle element
c        nu(1:3,ie) = three nodes of the trangle element
c        nu(ai,ie)  = di
c                   ai = (4:6) three segments (sides) of the element
c                        connnecting  nu(i-3,ie) -> nu(mod(i,3)+1,ie)
c                   di = 8*t + a
c                      t = triangle adjacent to the segment
c                      a = segment number in element t
c                   if di<0: segment ai is on meshing front, 
c                            and -di = pointer to the front node i 
c------------------------------------------------------------------------
c
c     1. initialization
c     -----------------
      do i=1,2*nbs-2
         nu(1,i) = i+1
         nu(2,i) = 0
         nu(3,i) = 0
         nu(4,i) = 0
         nu(5,i) = 0
         nu(6,i) = 0
      enddo
      nu(1,nbs+nbs-2) = 0
c
c     2. generate the first triangle and the front list
c     -------------------------------------------------
      free = 1
      t    = free
      free = nu(1,free)
      tete = free
      pf   = free
      do i=1,3
         nu(i  ,t) = tri(i)
         nu(3+i,t) = -pf
         ppf       = pf
         free      = nu(1,pf)
         pf        = free
         if ( i.eq.3 ) pf=tete
         nu(1,ppf) = tri(i)
         nu(2,ppf) = 8*t + i+3
         nu(ps,ppf) = pf
         nu(pp,pf ) = ppf
      enddo
c
c     3. creation of new triangles modification of the front
c     ------------------------------------------------------
      do i=4,nbs
         s = tri(i)
         pf = frntnd(c,nu,tete,s,nbs)
c
         t = free
         free = nu(1,free)
         npf  = free
         free = nu(1,free)
         ppf  = nu(pp,pf)
         psf  = nu(ps,pf)
         ttaf  = nu(2,pf)
         taf  = ttaf / 8 
         iaf  = ttaf - 8 * taf
c
c                  npf 
c               1  x s               ---
c                 / \                ---
c              4 /   \ 6        ---  new  ---
c               /  t  \              ---
c            2 /   5   \ 3           ---
c------ --<---x---------x---------- front--<---
c          psf \  iaf  /  pf         ---
c               \ taf /         --- omega ---
c                \   /               ---
c                 \ /                ---
c                  x ppf             ---
c                                    ---
c       3.1 generation of element t
c       ---------------------------
        nu(1,t) = s
        nu(2,t) = nu(1,psf)
        nu(3,t) = nu(1,pf )
        nu(4,t) = -npf
        nu(5,t) = 8*taf + iaf
        nu(6,t) = -pf
        nu(iaf,taf) = 8*t + 5
c
c       3.2 update the front list
c       -------------------------
        nu(ps,npf) = psf
        nu(pp,npf) = pf
        nu(ps,pf ) = npf
        nu(pp,psf) = npf
        nu(1,npf)  = s
        nu(2,npf)  = 8*t + 4
        nu(2,pf )  = 8*t + 6
c
c       3.3 check the generation
c       ------------------------
        call diagswp1 (c,nu,t,5,nbs,err)
        if ( err.ne.0 ) return
        call elmord (.true. ,c,nu,npf,nbs,err)
        if ( err.ne.0 ) return
        call elmord (.false.,c,nu,npf,nbs,err)
        if ( err.ne.0 ) return
      enddo
      end
c
c**********************************************************************
      integer function frntnd (c,nu,tete,s,nbs)
c     
c     given node s find the visible side on the front
c
      implicit none
      integer nbs,nu(6,nbs+nbs-2),tete,s
      double precision c(2,nbs),det,x,y
      integer pt,ppt
      logical init
c
      x = c(1,s)
      y = c(2,s)
      init = .true.
      pt = tete
10    continue
      ppt = pt
      pt = nu(4,pt)
      if ( pt.ne.tete ) then
         det = x*c(2,nu(1,pt)) - y*c(1,nu(1,pt))
         if ( det.lt.-1.d-15 ) then
            init = .false.
            goto 10
         elseif ( init.and.det.lt.1.d-15 ) then
            goto 10
         endif
      endif
      frntnd = ppt
      end
c
c**********************************************************************
      subroutine diagswp1 (c,nu,t,a,nbs,err)
c
c     diagonal swapping to optimize triangles, with meshing front info.
c----------------------------------------------------------------------
      implicit none
      integer nbs,nu(6,nbs+nbs-2),t,a,err
      double precision c(2,nbs)
      integer mxpile
      parameter (mxpile=256)
      integer pile(2,mxpile)
      integer t1,t2,i,s1,s2,s3,s4, idex
      double precision sin1,cos1,sin2,cos2
      integer tt1,tt,i11,i12,i13,i21,i22,i23,a1,a2,aa,p3(3)
      data p3/2,3,1/
      double precision x1,x2,x3,x4,y1,y2,y3,y4,det
c
      double precision epsilon
      data epsilon/-1.d-20/
c
c     initization
c     -----------
      err = 0
      idex = 0
5     i=1
      pile(1,i) = t 
      pile(2,i) = a
10    continue
      if ( i.gt.0 ) then
c
c       load two elements
c       -----------------
        t1 = pile(1,i)
        a1 = pile(2,i)
        i = i-1
        if ( t1.le.0 ) goto 10
        tt1 = nu(a1,t1)
        if ( tt1.le.0 ) goto 10
        t2 = tt1/8
        a2 = tt1-t2*8
        i11 = a1 - 3
        i12 = p3(i11) 
        i13 = p3(i12)
        i21 = a2 - 3
        i22 = p3(i21)
        i23 = p3(i22)
        s1 = nu(i13,t1)
        s2 = nu(i11,t1)
        s3 = nu(i12,t1)
        s4 = nu(i23,t2)
c
c       optimization of the quadrilateral s1,s2,s3,s4
c       ---------------------------------------------
        x1 = c(1,s1)
        x2 = c(1,s2) - x1
        x3 = c(1,s3) - x1
        x4 = c(1,s4) - x1
        y1 = c(2,s1)
        y2 = c(2,s2) - y1
        y3 = c(2,s3) - y1
        y4 = c(2,s4) - y1
        sin1 = y3*x2 - x3*y2
        cos1 = x3*(x3-x2) + y3*(y3-y2)
        sin2 = x4*y2 - y4*x2
        cos2 = x4*(x4-x2) + y4*(y4-y2)
        det = cos2*sin1 + cos1*sin2
c        if ( dabs(sin1).le.1.d-15 ) then
c           print *,'fatal error mshpot: 3 points collide ',s1,s2,s3
c           err = 20
c           return
c        endif
c        print *, i,det
        if ( det.ge.epsilon ) goto 10
c
c       switch diagonals of the quadrilateral, update nodes
c       ---------------------------------------------------
        nu(i12,t1) = s4
        nu(i22,t2) = s1
c
c       update segments a1,a2
c       ---------------------
        tt1 = nu(i22+3,t2)
        nu(a1 ,t1) = tt1
        if ( tt1.gt.0 ) then
           tt = tt1/8
           aa = tt1-8*tt
           nu(aa,tt) = 8*t1 + a1
        elseif ( tt1.lt.0 ) then
           nu(2,-tt1) = 8*t1 + a1
        endif
c
        tt1 = nu(i12+3,t1)
        nu(a2 ,t2) = tt1
        if ( tt1.gt.0 ) then
           tt = tt1/8
           aa = tt1-8*tt
           nu(aa,tt) = 8*t2 + a2
        elseif ( tt1.lt.0 ) then
           nu(2,-tt1) = 8*t2 + a2
        endif
c
        nu(i12+3,t1) = 8*t2 + i22+3
        nu(i22+3,t2) = 8*t1 + i12+3
        if ( i+4.gt.mxpile ) then
           if ( idex.eq.0 ) then
             epsilon = -1.d-15
             idex = 1
             go to 5
           endif
           print *,'mshpot fatal error: mxpile too small',mxpile
           err = 21
           return
        endif
c
        i = i+1
        pile(1,i) = t1
        pile(2,i) = a1
        i = i+1
        pile(1,i) = t2
        pile(2,i) = a2
        i = i+1
        pile(1,i) = t1
        pile(2,i) = i13+3
        i = i+1
        pile(1,i) = t2
        pile(2,i) = i23+3
        goto 10
      endif
      end
c
c**********************************************************************
      subroutine elmord (direct,c,nu,pfold,nbs,err)
c
c     re-order the nodes of a element, counter-clock-wisely
c
      implicit none
      integer nbs,nu(6,nbs+nbs-2),pfold,err
      double precision c(2,nbs),det
      logical direct
      integer pp,ps,i1,i2,i3,i4,i5,i6
      integer pf,psf,ppf,s1,s2,s3,t,t4,t5,a4,a5,tt4,tt5
c
      if ( direct ) then
        pp = 3
        ps = 4
        i1 = 1
        i2 = 3
        i3 = 2
        i4 = 6
        i5 = 5
        i6 = 4
      else
        pp = 4
        ps = 3
        i1 = 1
        i2 = 2
        i3 = 3
        i4 = 4
        i5 = 5
        i6 = 6
      endif
10    continue
      ppf = pfold
      pf  = nu(ps,pfold)
      psf = nu(ps,pf)
      s1  = nu(1,ppf)
      s2  = nu(1,pf)
      s3  = nu(1,psf)
      det =   ( c(1,s2) - c(1,s1) ) * ( c(2,s3) - c(2,s1) )
     &      - ( c(2,s2) - c(2,s1) ) * ( c(1,s3) - c(1,s1) )
      if ( (.not.direct .and. det.gt.0) .or. 
     &     (direct .and. det.lt.0) ) then
c        print *,'elmord convexification:',s1,s2,s3,det,direct
c
c       add triangle t and destruct one segment
c       ---------------------------------------
        if ( direct ) then
          tt4 = nu(2,ppf)
          tt5 = nu(2,pf)
        else
          tt4 = nu(2,pf)
          tt5 = nu(2,psf)
        endif
        t4 = tt4/8
        t5 = tt5/8
        a4 = tt4 - 8*t4
        a5 = tt5 - 8*t5
c
c       destruction of the front segment in pf
c       --------------------------------------
        nu(ps,ppf) = psf
        nu(pp,psf) = ppf
c
c       replace the front segment by element generation
c       -----------------------------------------------
        t = pf
        if ( direct ) then
          nu(2,ppf) = 8*t + i6
        else
          nu(2,psf) = 8*t + i6
        endif
c
c       on the created element
c       ----------------------
        nu(i1,t) = s1
        nu(i2,t) = s2
        nu(i3,t) = s3
        nu(i4,t) = 8*t4 + a4
        nu(i5,t) = 8*t5 + a5
        if ( direct ) then
          nu(i6,t) = -ppf
        else
          nu(i6,t) = -psf
        endif
        nu(a4,t4) = 8*t + i4
        nu(a5,t5) = 8*t + i5
        call diagswp1 (c,nu,t5,a5,nbs,err)
        if ( err.ne.0 ) return
        goto 10
      endif
      end
c
c**********************************************************************
      subroutine frtref (c,nu,nbs,nbt,bdseg,nba,sd,nbsd,reft,w,err)
c
c     treat the meshing front and define the element reference number
c
c      on return: nu(4-6,ie) >0 normal elements
c                 nu(4-6,ie) =0 neighboring a boundary segement
c                 nu(4-6,ie) <0 neighboring an interfacial segment
c-----------------------------------------------------------------
      implicit none
      integer nbs,nbt,nu(6,nbs+nbs-2),reft(2*(nbs-1))
      double precision c(2,nbs)
      integer nba,nbsd,bdseg(2,nba),sd(2,nbsd),err,w(*)
      integer i,j,is,ie,nbac,nbacpp,err1,itera,ne,nelem
      integer ss1,s1,s2,t,ta,is1,s2t,s3t,a,ap,isd,jsd
      double precision det2,det3, xcoor,x0,x21,xt2,xt3
      integer p3(1:3)
      integer vide
      parameter (vide=-2**30)
      data p3/2,3,1/
c
      if ( nba.eq.0 ) return
c
c     initialization of the boundary segments
c      sequencing and ordering the segments
c        reft(bd node #) = segment #
c     ---------------------------------------
      nbt = nbs+nbs-2
      do i=1,nbs
         reft(i) = 0
      enddo
      do i=1,nba
         reft(bdseg(1,i)) = vide
         reft(bdseg(2,i)) = vide
      enddo
      nbac = 0
      do a=1,nba
         s1 = min(bdseg(1,a),bdseg(2,a))
         s2 = max(bdseg(1,a),bdseg(2,a))
         if ( s1.ne.s2 ) then
            i = reft(s1)
 14         continue
	    if ( i.eq.vide ) then
               w(a) = reft(s1)
               reft(s1) = a
            else
               if ( s2.eq.max(bdseg(1,i),bdseg(2,i)) ) then
                  print *,'frtref warning: segment',i,' =(',s1,s2,
     &                  ') is the same as segment',a,' =(',
     &                  bdseg(1,a),bdseg(2,a),')'
                  nbac = nbac + 1
               else
                  i = w(i)
                  goto 14
               endif
            endif
         else
            print *,'frtref warning: nodes',s1,s2,' of segment',a,
     &              ' are the same'
            nbac = nbac + 1
         endif
      enddo
c
c     set the neighboring element information
c     ---------------------------------------
      nbacpp = 1       
      itera = 0
      err1 = 0
50    continue
      itera = itera + 1
      if ( err1.ne.0 ) then
         err = err1
         return
      endif
      if ( nbac.lt.nba ) then
c
         if ( nbacpp.eq.0 ) then
           print *,'frtref fatal error: algorithm fails,',
     &             nba,nbac,' iteration =',itera
           err = 7
           return
         endif
c
         nbacpp = 0
         do 120 ie=1,nbt
           if ( nu(5,ie).ne.0 ) then
             do 110 is=1,3
                s1  = nu(   is ,ie)
                s2t = nu(p3(is),ie)
                ss1 = min(s1,s2t)
                ap = 0
                a = reft(ss1)
80              continue
                if ( a.gt.0 ) then
                   s2 = max(bdseg(1,a),bdseg(2,a))
                   t  = ie
                   ta = 0
                   if ( s2.eq.max(s1,s2t) ) then
                      if ( nu(is+3,ie).gt.0 ) then
c                        reset the elements sharing the segment
                         ta       = nu(is+3,ie)/8
                         i        = nu(is+3,ie)-8*ta
                         nu(i,ta)    = -nu(i,ta)
                         nu(is+3,ie) = -nu(is+3,ie)
                      else
                         nu(is+3,ie) = vide
                      endif
                      goto 100
                   endif
                   ap = a
                   a  = w(a)
                   goto 80
                endif
c
c               if the boundary segment is not in the mesh
c
                if ( itera.eq.1 ) goto 110
                ss1 = s1
                ap = 0
                a = reft(ss1)
90              continue
                if ( a.gt.0 ) then
                   s2 = max(bdseg(1,a),bdseg(2,a))
                   t  = ie
                   ta = 0
c                  search element sharing segment a
                   is1  = is
                   s3t  = nu(p3(p3(is)),t)
                   x0 = c(1,s1)
                   xt2 = xcoor(x0,c(1,s2t)) - x0
                   xt3 = xcoor(x0,c(1,s3t)) - x0
                   x21 = xcoor(x0,c(1,s2))  - x0
                   det2 = xt2*(c(2,s2)-c(2,s1))-(c(2,s2t)-c(2,s1))*x21
                   det3 = xt3*(c(2,s2)-c(2,s1))-(c(2,s3t)-c(2,s1))*x21
                   if ( det2.gt.0 .and. det3.lt.0 ) then
c                     node s2 is between s2t and s3t
                      call mshfr1 (c,nu,nbs,t,ta,is1,s2,err)
c                       if ( err.ne.0 ) return
                       if ( err.ne.0 ) goto 101
                      goto 100
                   elseif ( det2.eq.0 .and. reft(s2t).eq.0 ) then
c                     node s2 is along s2t, however s2t is blocked
                      print *,'frtref fatal error: point ',s2t,
     &                     ' is blocked from the mesh front'
                      err1 = 10
                   elseif ( det3.eq.0 .and. reft(s3t).eq.0 ) then
c                     node s2 is along s3t, however s3t is blocked
                      print *,'frtref fatal error: point ',s3t,
     &                     ' is blocked from the mesh front'
                      err1 = 10
                   endif
                   ap = a
                   a  = w(a)
                   goto 90
                endif
                goto 110
100             continue
                nbacpp = nbacpp + 1
                if ( ap.eq.0 ) then
                   reft(ss1) = w(a)
                else
                   w(ap) = w(a)
                endif
                if ( nbac+nbacpp.eq.nba ) goto 130
110          continue
          endif
120     continue
        nbac = nbac + nbacpp
        goto 50
      endif
c
130   continue
c
c     compute the element reference number
c     ------------------------------------
      do i=1,nbs+nbsd+nbsd
         w(i) = 0
      enddo
      do i=1,nbsd
         a = abs(sd(1,i))
         s1 = min(bdseg(1,a),bdseg(2,a))
         w(i+i) = w(s1+nbsd+nbsd)
         w(s1+nbsd+nbsd) = i
      enddo
c
c     identify the element sharing side with a specified BD segment
c
      do 180 t=1,nbt  
         reft(t) = vide
         if ( nu(6,t).ne.0 ) then
            do 170 i=1,3
               ss1 = min(nu(i,t),nu(p3(i),t))
               jsd = nbsd+nbsd+ss1
160            continue           
               isd = w(jsd)
               if ( isd.gt.0 ) then
                 a = sd(1,isd)
                 if ( a.gt.0 ) then
                    if ( nu(i,t).eq.bdseg(1,a) .and.
     &                   nu(p3(i),t).eq.bdseg(2,a) ) then
                       reft(t) = sd(2,isd)
                       w(isd+isd-1) = t
                       w(jsd)=w(isd+isd)
                       goto 170
                    endif
                 elseif ( a.lt.0 ) then
                    if ( nu(i,t).eq.bdseg(2,-a) .and. 
     &                   nu(p3(i),t).eq.bdseg(1,-a) ) then
                       reft(t) = sd(2,isd)
                       w(isd+isd-1) = t
                       w(jsd) = w(isd+isd)
                       goto 170
                    endif
                 else
                    print *,'fatal error in domain',isd,'zero segment'
                    err = 11
                 endif
                 jsd = isd+isd
                 goto 160
               endif
170         continue         
         endif
180   continue
c
      do isd=1,nbsd
         if ( w(isd+isd-1).eq.0 ) then
            print *,' fatal error in frtref: domain ',isd
     &             ,' does not border any element '
            err = 11
         else
            w(isd+isd) = 3
         endif
      enddo
      if ( err.ne.0 ) then
        print *,'fatal error frtref: definition of subdomain'
        return
      endif
c
c     special case: for isolated elements
c
      do t=1,nbt
         if (nu(4,t).lt.0 .and. nu(5,t).lt.0 .and. nu(6,t).lt.0 )
     &      nu(1,t) = -nu(1,t)
      enddo
c
c     set reference number for other elements
c
      i = nbsd+nbsd
      do while ( i.gt.0 )
         w(i) = w(i)+1
         if ( w(i).le.6 ) then
            ta = nu(w(i),w(i-1))
            if ( ta.gt.0 ) then
               ta = ta/8
               if ( nu(1,ta).gt.0 ) then
                  nu(1,ta) = -nu(1,ta)
                  if ( reft(ta).ne.reft(w(i-1)) ) then
                     if ( reft(ta).ne.vide ) then
                        print *,'frtref error: domain element',ta,
     &                  ' ref old=',reft(ta),' ref new=',reft(w(i-1))
                         err = 10
                     else
                        reft(ta) = reft(w(i-1))
                     endif
                     w(i+1) = ta
                     w(i+2) = 3
                     i = i+2
                  endif
               endif
            endif
         else
            i = i-2
         endif
      enddo
c
c     clean-up nu and reft
c     --------------------
101   ne = 0
      do ie=1,nbt
         if ( nu(1,ie).ge.0 .or. nu(5,ie).eq.0 .or.
     &        reft(ie).eq.vide ) then
            w(ie) = 0
         else
            ne =  ne+1
            w(ie) = ne
         endif
      enddo
      nelem = ne
c
      do ie=1,nbt
         ne = w(ie)
         if ( ne.ne.0 ) then
            nu(1,ne) = -nu(1,ie)
            do j=2,3
              nu(j,ne) = nu(j,ie)
            enddo
            do j=4,6
               if ( nu(j,ie).gt.0 ) then
                  t = nu(j,ie)/8
                  i = nu(j,ie) - t*8 
                  nu(j,ne) = 8*w(t)+i
               elseif ( nu(j,ie).eq.vide ) then
                  nu(j,ne) = 0
               else
                  t = abs(nu(j,ie))/8
                  i = abs(nu(j,ie)) - t*8
                  nu(j,ne) = -(8*w(t)+i)
                  if ( reft(t).eq.vide ) nu(j,ne) = 0
               endif
            enddo
c
            reft(ne) = reft(ie)
         endif
      enddo
c
      do ie=nelem+1,nbt
         do i=1,6
            nu(i,ie) = 0
         enddo
         reft(ie) = 0
      enddo
c
      nbt = nelem
      return
      end
c
c**********************************************************************
      subroutine mshfr1 (c,nu,nbs,it1,ita,is1,s2,err)
c
c     reconnect the mesh such that a given boundary segment is one of
c       the edge in the mesh
c 
      implicit none
      integer nbs,nu(6,nbs+nbs-2),is1,s2,err,it1,ita
      double precision c(2,nbs),x,y,det, x1,xcoor
      integer mxpile
      parameter (mxpile=256)
      integer lst(3,mxpile)
      integer s1,s3,nbac,s2t,s3t,t,ta
      integer l1,l2,l3,la,p3(1:5)
      data p3 /2,3,1,2,3/
c
      t = it1
      s1 = nu(is1,t)
      x1 = c(1,s1)
      x = xcoor(x1,c(1,s2)) - x1
      y = c(2,s2) - c(2,s1)
      nbac = 0
      l1 = is1
      l2 = p3(l1)
      l3 = p3(l2)
      s2t = nu(l2,t)
      s3t = nu(l3,t)
      la = l2 + 3
c      print *,'  mshfr1 :',it1,is1,s1,s2
20    continue
      nbac = nbac + 1
      if ( nbac.gt.mxpile ) then
         print *,'mshfr1 error: mxpile too small',nbac,mxpile
         err = 8
         return
      endif
      lst(2,nbac) = t
      lst(3,nbac) = la
      ta = nu(la,t)
      if ( ta.le.0 ) then
         print *,'mshfr1 error: front is blocked at',t
         err = 9
         return
      endif
      t  = ta/8
      la = ta-8*t
      s3 = nu(p3(la-2),t)
      if ( s3.ne.s2 ) then
         det = x*(c(2,s3)-c(2,s1)) - y*(xcoor(x1,c(1,s3))-x1)
         if ( det.gt.0 ) then
c           s2 is between s2t and s3
            la = 3 + p3(la-3)
         elseif ( det.lt.0 ) then
c           s2 is between s3 and s3t
            la = 3 + p3(la-2)
         else
            print *,'mshfr1 error: point',s3
     &             ,' is not in the right direction'
            err = 10
            return
         endif
         goto 20
      endif
c
c     reconstruct the mesh with s1-s2 as an edge 
      call mshfr2 (c,nu,nbs,lst,nbac,it1,ita,s1,s2,err)
      return
      end
c
c**********************************************************************
      subroutine mshfr2 (c,nu,nbs,lst,nbac,t,ta,ss1,ss2,err)
      implicit none
      integer nbs,nbac,nu(6,nbs+nbs-2),lst(3,nbac)
      double precision c(2,nbs)
      integer t,ta,ss1,ss2,err
      integer ptlst,ttlst,pslst,pplst,s1,s2,s3,s4
      integer i,t1,a1,tt1,t2,a2,tt,i11,i12,i13,i21,i22,i23,aas,aa
      double precision det1,det4,det2,det3,x,y,x41,y41
     &                ,xcoor,x0,x1
      integer p3(3)
      integer vide
      parameter (vide=-2**30)
      data p3/2,3,1/
c
      x0 = c(1,ss1)
      x  = x0 - xcoor(x0,c(1,ss2))
      y  = c(2,ss1) - c(2,ss2)
      do i=1,nbac-1
         lst(1,i) = i+1
      enddo
      lst(1,nbac) = 0
      ttlst = 1
20    continue
      ptlst = ttlst
      pplst = 0
30    continue
      if ( ptlst.gt.0 ) then
         t1 = lst(2,ptlst)
         a1 = lst(3,ptlst)
         tt1 = nu(a1,t1)
         t2 = tt1/8
         a2 = tt1-t2*8
         i11 = a1 -3
         i12 = p3(i11) 
         i13 = p3(i12)
         i21 = a2 -3
         i22 = p3(i21)
         i23 = p3(i22)
         s1 = nu(i13,t1)
         s2 = nu(i11,t1)
         s3 = nu(i12,t1)
         s4 = nu(i23,t2)
         x1  = xcoor(x0,c(1,s1))
         x41 = xcoor(x0,c(1,s4)) - x1
         y41 = c(2,s4) - c(2,s1)
         det2 = (xcoor(x0,c(1,s2))-x1)*y41 - (c(2,s2)-c(2,s1))*x41
         det3 = (xcoor(x0,c(1,s3))-x1)*y41 - (c(2,s3)-c(2,s1))*x41
         if ( det2.gt.0 .and. det3.lt.0 ) then
c
c           le quadrilataire est convexe on le retourne
c           update the nodes
c           -------------------------
            nu(i12,t1) = s4
            nu(i22,t2) = s1
c
c           update the previous pointer
c           ---------------------------
            pslst = lst(1,ptlst)
            if ( pslst.gt.0 ) then
               aas = lst(3,pslst)
               if ( aas.eq.i22+3 ) then
                  lst(2,pslst) = t1
                  lst(3,pslst) = i11 + 3
               endif
            endif
c
c           update the mesh
c           ---------------
            tt1 = nu(i22+3,t2)
            nu(a1 ,t1) = tt1
            tt=abs(tt1)/8
            aa = abs(tt1)-8*tt
            if ( tt1.gt.0 ) then
               nu(aa,tt) = a1 + 8*t1
            elseif ( tt1.ne.vide ) then
               nu(aa,tt) = -(a1 + 8*t1)
            endif
c
            tt1 = nu(i12+3,t1)
            nu(a2,t2) = tt1
            tt=abs(tt1)/8
            aa = abs(tt1)-8*tt
            if ( tt1.gt.0 ) then
               nu(aa,tt) = a2 + 8*t2
            elseif ( tt1.ne.vide .and. tt.gt.0) then
               nu(aa,tt) = -(a2 + 8*t2)
            endif
c
            nu(i12+3,t1) = i22+3 + 8*t2
            nu(i22+3,t2) = i12+3 + 8*t1
            det1 = (xcoor(x0,c(1,s1))-x0)*y - (c(2,s1)-c(2,ss1))*x
            det4 = (xcoor(x0,c(1,s4))-x0)*y - (c(2,s4)-c(2,ss1))*x
            if ( det1.lt.0 .and. det4.gt.0 ) then
c
c              node s4 is in omega
               lst(2,ptlst) = t2
               lst(3,ptlst) = i22+3
            elseif ( det1.gt.0 .and. det4.lt.0 ) then
c              node s1 is in omega
               lst(2,ptlst) = t1
               lst(3,ptlst) = i12+3
            else
c              print *,' on supprime l''segment dans  lst ',t1,a1,t2,a2
               if ( pplst.eq.0 ) then
                  ttlst = lst(1,ptlst)
                  ptlst = ttlst
               else
                  ptlst        = lst(1,ptlst)
                  lst(1,pplst) = ptlst
               endif
               goto 30
            endif
         endif
         pplst = ptlst
         ptlst = lst(1,ptlst)
         goto 30
      endif
      if ( ttlst.ne.0 ) goto 20
      nu(i12+3,t1) = -nu(i12+3,t1)
      nu(i22+3,t2) = -nu(i22+3,t2)
      t  = t2      
      ta = t1
      do i=1,nbac
         call diagswp (c,nu,lst(2,i),4,nbs,err)
         call diagswp (c,nu,lst(2,i),5,nbs,err)
         call diagswp (c,nu,lst(2,i),6,nbs,err)
      enddo
      end
c
c**********************************************************************
      subroutine mshrfn (c,nu,h,h0,reft,nbs,nbsmx,nbt,puis,err)
c
c     modify the mesh, if the mesh size is too large, refine the mesh!
c----------------------------------------------------------------------
      implicit none
      integer nbs,nbsmx,nbt,err
      integer nu(6,2*(nbsmx-1)),reft(2*(nbsmx-1))
      double precision c(2,nbsmx),h(nbsmx),puis,h0
      integer t,s1,s2,s3,itera,nbsold
      double precision x,y,aire,hs,det,x1,x2,x3,y1,y2,y3,xcoor,xploc
c
c     coef         controls the number of interior vertices which will
c                  be generated
c                  .75 => a good value
      double precision coef
      parameter (coef = 0.75)
c
      err = 0
      itera = 0
20    continue
      itera = itera + 1
      nbsold = nbs 
      do t=1,nbt
         s1 = nu(1,t)
         s2 = nu(2,t)
         s3 = nu(3,t)
         x1 = c(1,s1)
         x2 = xcoor(x1,c(1,s2)) - x1
         x3 = xcoor(x1,c(1,s3)) - x1
         y1 = c(2,s1)
         y2 = c(2,s2) - y1
         y3 = c(2,s3) - y1
c
c        area of the trangle
c        -------------------
         det = x2*y3 - y2*x3
         if ( det.le.0 ) then
            print *,'element ',t,' has negative area=',det
            print *,'s1,s2,s3,(x1,y1),(x2,y2),(x3,y3)=',s1,s2,s3,
     &               x1,y1,c(1,s2),c(2,s2),c(1,s3),c(2,s3)
            err = 5
            return
         endif
         aire = det*coef
         hs = ((h(s1)**puis+h(s2)**puis+h(s3)**puis)/3.)**(1./puis)
         if ( aire.gt.hs*hs ) then
           x = x1 + (x2+x3)*0.333333333333333d0
           y = y1 + (y2+y3)*0.333333333333333d0
           if ( nbs.ge.nbsmx ) then
              print *,'mshrfn: mesh points=',nbs,' > maximum=',nbsmx
              err = 10
              return
           endif
           nbs = nbs + 1
           c(1,nbs) = xploc(x)
           c(2,nbs) = y
           h(nbs) = max(hs, min(dsqrt(0.33333333333333d0*aire),h0))
           call elmdvd (t,nbs,c,nu,reft,nbt,nbs,err)
           if(err.ne.0) return
         endif
      enddo
      if ( nbsold.ne.nbs ) goto 20
      end
c
c**********************************************************************
      subroutine elmdvd (t,s,c,nu,reft,nbt,nbs,err)
c
c     divide the current triangle into three smaller ones
c----------------------------------------------------------------------
      implicit none
      integer t,s,nbt,nbs,err
      integer nu(6,2*(nbs-1)),reft(2*(nbs-1))
      double precision c(2,nbs)
      integer t1,t2,t3,ta2,ta3,ia2,ia3,tta
c
      t1 = t
      nbt = nbt+1
      t2 = nbt
      nbt = nbt+1
      t3 = nbt
c
      nu(1,t2) = s
      nu(2,t2) = nu(2,t)
      nu(3,t2) = nu(3,t)
      nu(4,t2) = 8*t1 + 5
      nu(5,t2) = nu(5,t)
      nu(6,t2) = 8*t3 + 5
c
      nu(1,t3) = nu(1,t)
      nu(2,t3) = s
      nu(3,t3) = nu(3,t)
      nu(4,t3) = 8*t1 + 6
      nu(5,t3) = 8*t2 + 6
      nu(6,t3) = nu(6,t)
c
      tta = nu(5,t)
      ta2 = abs(tta)/8
      ia2 = abs(tta) - 8*ta2
      if ( tta.gt.0 ) then
         nu(ia2,ta2) = 8*t2 + 5
      elseif ( tta.lt.0 .and. ta2.gt.0 ) then
         nu(ia2,ta2) = -(8*t2 + 5)
      endif
c
      tta = nu(6,t)
      ta3 = abs(tta)/8
      ia3 = abs(tta) - 8*ta3
      if ( tta.gt.0 ) then
         nu(ia3,ta3) = 8*t3 + 6
      elseif ( tta.lt.0 .and. ta3.gt.0 ) then
         nu(ia3,ta3) = -(8*t3 + 6)
      endif
c
      nu(3,t1) = s
      nu(5,t1) = 8*t2 + 4
      nu(6,t1) = 8*t3 + 4
c
      reft(t2) = reft(t)
      reft(t3) = reft(t)
c
      call diagswp (c,nu,t1,4,nbs,err)
      if ( err.ne.0 ) return
      call diagswp (c,nu,t2,5,nbs,err)
      if ( err.ne.0 ) return
      call diagswp (c,nu,t3,6,nbs,err)
      if ( err.ne.0 ) return
      end
c
c**********************************************************************
      subroutine diagswp (c,nu,t,a,nbs,err)
c
c     diagonal swapping to optimize triangles, element level
c----------------------------------------------------------------------
      implicit none
      integer nbs,nu(6,nbs+nbs-2),t,a,err
      double precision c(2,nbs)
      integer mxpile
      parameter (mxpile=256)
      integer pile(2,mxpile)
      integer t1,t2,i,s1,s2,s3,s4, idex
      double precision sin1,cos1,sin2,cos2
      integer tt1,tt,i11,i12,i13,i21,i22,i23,a1,a2,aa,p3(3)
      data p3/2,3,1/
      double precision x1,x2,x3,x4,y1,y2,y3,y4,det,xcoor
c
      double precision epsilon
      data epsilon/-1.d-20/
c
c     initization
c     -----------
      err = 0
      idex = 0
5     i=1
      pile(1,i) = t 
      pile(2,i) = a
10    continue
      if ( i.gt.0 ) then
c
c       load two elements
c       -----------------
        t1 = pile(1,i)
        a1 = pile(2,i)
        i = i-1
        if ( t1.le.0 ) goto 10
        tt1 = nu(a1,t1)
        if ( tt1.le.0 ) goto 10
        t2 = tt1/8
        a2 = tt1-t2*8
        i11 = a1 - 3
        i12 = p3(i11) 
        i13 = p3(i12)
        i21 = a2 - 3
        i22 = p3(i21)
        i23 = p3(i22)
        s1 = nu(i13,t1)
        s2 = nu(i11,t1)
        s3 = nu(i12,t1)
        s4 = nu(i23,t2)
c
c       optimization of the quadrilateral s1,s2,s3,s4
c       ---------------------------------------------
        x1 = c(1,s1)
        x2 = xcoor(x1,c(1,s2)) - x1
        x3 = xcoor(x1,c(1,s3)) - x1
        x4 = xcoor(x1,c(1,s4)) - x1
        y1 = c(2,s1)
        y2 = c(2,s2) - y1
        y3 = c(2,s3) - y1
        y4 = c(2,s4) - y1
        sin1 = y3*x2 - x3*y2
        cos1 = x3*(x3-x2) + y3*(y3-y2)
        sin2 = x4*y2 - y4*x2
        cos2 = x4*(x4-x2) + y4*(y4-y2)
        det = cos2*sin1 + cos1*sin2
c        print *, i,det
        if ( det.ge.epsilon ) goto 10
c
c       switch diagonals of the quadrilateral, update nodes
c       ---------------------------------------------------
        nu(i12,t1) = s4
        nu(i22,t2) = s1
c
c       update segments a1,a2
c       ---------------------
        tt1 = nu(i22+3,t2)
        nu(a1 ,t1) = tt1
        tt = abs(tt1)/8
        aa = abs(tt1)-8*tt
        if ( tt1.gt.0 ) then
           nu(aa,tt) = 8*t1 + a1
        elseif ( tt1.lt.0 .and. tt.gt.0 ) then
           nu(aa,tt) = -(8*t1 + a1)
        endif
c
        tt1 = nu(i12+3,t1)
        nu(a2 ,t2) = tt1
        tt = abs(tt1)/8
        aa = abs(tt1)-8*tt
        if ( tt1.gt.0 ) then
           nu(aa,tt) = 8*t2 + a2
        elseif ( tt1.lt.0 .and. tt.gt.0 ) then
           nu(aa,tt) = -(8*t2 + a2)
        endif
c
        nu(i12+3,t1) = 8*t2 + i22+3
        nu(i22+3,t2) = 8*t1 + i12+3
        if ( i+4.gt.mxpile ) then
           if ( idex.eq.0 ) then
             epsilon = -1.d-15
             idex = 1
             go to 5
           endif
           print *,'diagswp fatal error: mxpile too small',mxpile
           err = 21
           return
        endif
c
        i = i+1
        pile(1,i) = t1
        pile(2,i) = a1
        i = i+1
        pile(1,i) = t2
        pile(2,i) = a2
        i = i+1
        pile(1,i) = t1
        pile(2,i) = i13+3
        i = i+1
        pile(1,i) = t2
        pile(2,i) = i23+3
        goto 10
      endif
      end
c
c**********************************************************************
      subroutine mshvoi (nu,w1,w,nbt,nbs)
c
c     Find points connected to each current node
c     w1: total # of points connected, w: record of these points
c----------------------------------------------------------------------
      implicit none
      integer nbt,nbs,nu(6,nbt),w1(0:nbs),w(3*nbt),i,j,is
      do i=0,nbs
        w1(i) = 0
      enddo
      do i=1,nbt
        do j=1,3
          w1(nu(j,i)) = w1(nu(j,i)) + 1
        enddo
      enddo
      do i=1,nbs
        w1(i)= w1(i-1) + w1(i)
      enddo
      do i=1,nbt
        do j=1,3
          is = nu(j,i) -1
          w1(is)    = w1(is) + 1
          w(w1(is)) = 8*i+j
        enddo
      enddo
      do i=nbs,1,-1
        w1(i) = w1(i-1)
      enddo
      w1(0) = 0
      end
c
c**********************************************************************
      subroutine mshsmth (c,nbd,nbs,nu,w1,w,nbt,omega,itermx,eps)
c
c     Smoothing of generated nodes -> the centroid of the polygon
c       smooth the interior nodes only (nbd+1 -> nbs)
c----------------------------------------------------------------------
      implicit none
      integer itermx
      integer nbd,nbs,nbt,nu(6,nbt),w1(0:nbs),w(3*nbt)
      double precision bx,by,err,omega,depx,depy,eps,dx,ymin,ymax
      integer iter,i1,i2,i,is,k,ic,nd,ta,ia
      double precision c(2,nbs)
      integer p3(1:3)
      data p3/2,3,1/
c
      double precision xcoor,xploc
c
c     calculate the mesh size dx
c     --------------------------
      ymin = c(2,1)
      ymax = c(2,1)
      do ic=2,nbs
        ymin = min(c(2,ic),ymin)
        ymax = max(c(2,ic),ymax)
      enddo
      dx = ymax-ymin
c
c     update nodes to the centroid of the polygon
c     -------------------------------------------
      do iter=1,itermx
         err = 0
         do is=nbd+1,nbs
            i1 = w1(is-1) + 1
            i2 = w1(is)
            if ( i2.ge.i1 ) then
               bx = 0.d0
               by = 0.d0
               do i=i1,i2
                 k = w(i)
                 ta = k/8
                 ia = k - ta*8
                 nd = nu(p3(ia),ta)
                 bx = bx + xcoor(c(1,is),c(1,nd))
                 by = by + c(2,nd)
               enddo
               bx = bx/(i2-i1+1.d0)
               by = by/(i2-i1+1.d0)
               depx = omega*(c(1,is)-bx)
               depy = omega*(c(2,is)-by)
               c(1,is) = xploc(c(1,is) - depx)
               c(2,is) =       c(2,is) - depy
               err = max(dabs(depx),dabs(depy),err)
            endif
         enddo
         if ( err.le.eps*dx ) return
c
      enddo
      print *,'mshsmth: not converged aft.',itermx,' iterations'
      end
