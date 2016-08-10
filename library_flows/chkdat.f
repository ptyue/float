c***********************************************************************
      subroutine chkdat (iou, icase,iax,kaxi, neq,nvar,
     &                   resold,rstart,upwind,newton,inert,
     &                   trans,dtin,endtim,bdpr,iforce,itmax,zconv,
     &                   nvert,nnode,nelem,inod,ngmx,x,y,ishape,
     &                   nbd,ibdnod,nic,ic, ipnode,isnode,
     &                   ibdu,ibdt,ibds, densu,denst,denss,ncpu,ncps,
     &                   ifluid,ivisc,itemp,itrel,ivheat,istres,icoe,
     &                   ro,grav,pvis,pterm,pelas, iorder, ndgstrm,
     &                  maxnod,maxelt,maxnbd,maxbdp,maxcoe,maxvar,maxnz)
c
c     this routine reads and checks data 
c***********************************************************************
      implicit none
      integer iou,icase,iax,kaxi,nnode,nvert,nbd,nelem,nvar,nz,nic,
     &        ivisc,itrel,ipnode,isnode,maxnod,maxelt,
     &        maxnbd,maxvar,maxnz,maxbdp,maxcoe,iforce,itmax,
     &        ngmx,ncpu,ncps,ndgstrm
      integer inod(ngmx,nelem),ibdnod(nbd+1),ic(nic+1),ishape(nelem),
     &        ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd)
      real*8 x(nnode),y(nnode),dtin,endtim,zconv,
     &       densu(ncpu,nbd),denst(nbd),denss(ncps,nbd)
      logical bdpr,upwind,trans,newton,inert,resold,rstart
      integer neq
      integer ifluid,itemp,iorder,icoe,istres
      logical ivheat
      double precision ro(icoe),grav(2,icoe),pvis(10,icoe),
     &                 pterm(10,icoe),pelas(10,icoe)
c
c     local varaibles
      integer i,idp
      real*8 valdp
c
c     test problem data
c     -----------------
      kaxi = 1
      if ( iax.eq.0 ) kaxi = 0
      if ( icase.eq.3 )  istres = 0
      if ( icase.eq.3 .or. icase.eq.5 ) iforce = 0
c
      nz = neq + ndgstrm
c
      call tesdat ( iou,icase,ifluid,iax, nnode,maxnod,nelem,maxelt,nbd,
     &              maxnbd,nvar,maxvar,nz,maxnz, nic,maxbdp,icoe,maxcoe,
     &              ivisc,itemp,itrel, istres,iforce,iorder, pvis,pelas,
     &              pterm,ro,grav,zconv,dtin,endtim,ivheat,inert,trans,
     &              resold, rstart, upwind )
c
c     computation ishape
c     ------------------
c      ibdnod(nbd+1) = ibdnod(1)
      do i=1,nelem
        ishape(i) = 1
        if ( ngmx.eq.3 ) then
          ishape(i) = 0
        elseif ( ngmx.eq.4 ) then
          if ( inod(4,i).eq.0 ) ishape(i) = 0
        elseif ( ngmx.eq.6 ) then
           ishape(i) = 0
        elseif ( ngmx.gt.6 ) then
          if( inod(7,i).eq.0 ) ishape(i) = 0
        endif
      enddo
c
      call tesbnd ( iou,nbd,nic,ic,icase,istres,ipnode,isnode)
c
      write(iou,'(/A)') 
     &            'chkdat: all checked data from the data file are OK.'
c
c     print problem information
c     -------------------------
      call pridat (iou, nnode,nvert,nbd,nelem,nvar, inod,ngmx,
     &             istres, icase, ifluid,iax, itmax,ivisc,itemp,
     &             itrel, iorder,iforce,icoe, x,y, ro,grav, pvis,pelas,
     &             pterm,dtin,endtim, zconv, trans, newton,
     &             inert,upwind,resold,rstart,ivheat )
c
      call pribnd (iou,bdpr,nbd,icase,istres,nic,ipnode,isnode,
     &        ic, ibdu,ibdt,ibds,idp,densu,denst,denss,ncpu,ncps,valdp )
c
      return
      end
c
c
c***********************************************************************
      subroutine tesdat (iou,icase,ifluid,iax,nnode,maxnod,nelem,maxelt,
     &                   nbd,maxnbd, nvar,maxvar, nz,maxnz, nic, maxbdp,
     &                   icoe,maxcoe, ivisc,itemp,itrel, istres, iforce,
     &                   iorder,pvis,pelas,pterm, ro,grav, zconv, dtin, 
     &                   endtim,ivheat,inert,trans,resold,rstart,upwind)
c
c     this routine checks the input data 
c***********************************************************************
      implicit none
      integer nnode,maxnod,nelem,maxelt,nbd,maxnbd,nvar,maxvar,nz,maxnz,
     &        nic,maxbdp,icoe,maxcoe
      integer icase,ifluid,iax,iou
      integer ivisc,itemp,itrel,istres,iforce, iorder
      double precision ro(icoe),grav(2,icoe),pvis(10,icoe),
     &                 pterm(10,icoe),pelas(10,icoe)
      real*8 zconv,dtin,endtim
      logical ivheat,inert,trans,upwind,resold,rstart
c
c     local variables
      logical istop,range,odd,const2,const3,const4
      integer i,i1,i2,i3
      real*8 a,b,c,d
c
c     0. functions
c     ------------
      range(i1,i2,i3) = ( i2.ge.i1 .and. i2.le.i3 )
      odd(i) = ( ((i/2)*2).ne.i )
      const2(a,b) = ( a.eq.0.d0 .and. b.eq.0.d0 )
      const3(a,b,c) = ( a.eq.0.d0 .and. b.eq.0.d0 .and. c.eq.0.d0 )
      const4(a,b,c,d) = ( a.eq.0.d0 .and. b.eq.0.d0 .and. c.eq.0.d0
     &                    .and. d.eq.0.d0 )
c
c     1. initialization
c     -----------------
      istop = .false.
c
      write(iou,'(A)') 'tesdat:'
c
c     2. test of icase, ifluid and iax
c     --------------------------------
      if ( .not.range(1,icase,5) ) then
        write(iou,*) 'icase=',icase,' is out of range 1 - 5'
        istop = .true.
      endif
      if ( .not.range(0,ifluid,3) ) then
        write(iou,*) 'ifluid=',ifluid,' is out of range 0 - 3'
        istop = .true.
      endif
      if ( .not.range(0,iax,2) ) then
        write(iou,*) 'iax=',iax,' is out of range 0 - 2'
        istop = .true.
      endif
c
c     3. test of nnode, nelem and nbd
c     -------------------------------
      if ( .not.range(1,nnode,maxnod) ) then
        write(iou,*) 'number of nodes nnode=',nnode,' must < ',maxnod
        istop = .true.
      endif
      if ( .not.range(1,nelem,maxelt) ) then
        write(iou,*) 'number of elements nelem=',nelem,' must < ',maxelt
        istop = .true.
      endif
      if ( .not.range(2,nbd,maxnbd) ) then
        write(iou,*) 'number of bd. nodes nbd=',nbd,' must < ',maxnbd
        istop = .true.
      endif
      if ( .not.range(2,nic,maxbdp) ) then
        write(iou,*) 'number of bd. segments nic=',nic,' must < ',maxbdp
        istop = .true.
      endif
      if ( odd(nbd) ) then
        write(iou,*) 'boundary nodes must be even, nbd=',nbd
        istop = .true.
      endif
c
c     4. test of resold,rstart,zconv,iforce,istres
c     ------------------------------------------------------
      if ( .not.(resold.or.rstart) .and. icase.eq.5 ) then
        write(iou,*) 'a result file is required when icase = 5'
        istop = .true.
      endif
      if ( zconv.lt.0.d0 ) then
        write(iou,*) 'relative convergence value zconv =',zconv,
     &               ' must > 0'
        istop = .true.
      endif
      if ( .not.range(0,iforce,1) ) then
        write(iou,*) 'iforce =',iforce,' must be 0 or 1. it is set to 0'
        iforce = 0
      endif
      if ( .not.range(0,istres,2) ) then
        write(iou,*) 'istres=',istres,' must be 0, 1 or 2.'
        istop = .true.
      endif
      if ( istres.eq.2 .and. icase.eq.2 ) then
        write(iou,'(A,i4/A)')
     &         'viscoelastic stress postprocessor  istres =',istres,
     &         'is not allowed for non-isothermal flows icase = 2'
        istop = .true.
      endif
      if ( istres.eq.0 .and. icase.eq.5 ) then
        write(iou,'(A/A,i4,A)')
     &         'stress postprocessing is required for icase = 5',
     &       'istres =',istres,' must be 1 newtonian or 2 viscoelastic'
        istop = .true.
      endif
      if ( istres.eq.2 .and. icase.eq.1 ) then
        write(iou,'(A,i4/A/A/A/A)')
     &  'viscoelastic stress postprocessor  istres =',istres,
     &  'is not available immediately for isothermal flows icase=1',
     &  'a run of polyflow with istres=0 will solve velocity field',
     &  'a second run with icase=5&istres=2 solves viscoelastic stress',
     &  'execution continues with istres=0'
        istres = 0
      endif
c
c     5. test of physical data: icoe,pvis(1),ro,ivisc,itemp
c     ------------------------------------------------------
      if ( icoe.gt.maxcoe ) then
        write(iou,*) 'icoe=',icoe,' must < maxcoe=',maxcoe
        istop = .true.
      endif
      if ( icase.ne.3 ) then
        do i=1,icoe
          if ( pvis(1,i).lt.0.d0 ) then
            write(iou,*) '0-shear rate viscosity =',pvis(1,i),' < 0.'
            istop = .true.
            endif
        enddo
        if ( .not.range(1,ivisc,6) ) then
          write(iou,*) 'shear-viscosity law ivisc=',ivisc,
     &                 ' is out of range 1 - 6.'
          istop = .true.
        endif
        if ( icase.eq.2 ) then
          if ( .not.range(1,itemp,3) ) then
            write(iou,*) 'thermal viscosity law itemp=',itemp,
     &                   ' is out of range 1 - 3.'
            istop = .true.
          endif
        else
          ivheat = .false.
        endif
      endif
c
      do i=1,icoe
        if ( icase.ne.3 .and. ro(i).lt.0.d0 ) then
          write(iou,*) 'density ro < 0, ro=',ro(i)
          istop = .true.
        endif
      enddo
c
c     6. test of gravity field
c     ------------------------
      do i=1,icoe
      if ( iax.ne.0 .and. grav(1,i).ne.0.d0 ) then
        if ( icase.eq.3 ) then
          write(iou,*) 'axisymmetric problem with convection velocity',
     &                 ' which destroys the symmetry'
        else
          write(iou,*) 'axisymmetric/swirling flow with gravity (gx)',
     &                 ' which destroys the symmetry'
        endif
        istop = .true.
      endif
      enddo
c
c     8. test of data compatibility for streamline upwinding
c     ------------------------------------------------------
c      if ( upwind ) then
c        do i=1,icoe
c          if ( ro(i).eq.0.d0 .or. const2(grav(1,i),grav(2,i)) ) then
c            write(iou,'(A/A/A)')
c     &           'streamline upwinding only has sense in the cases:',
c     &           ' - ro not equal to zero with icase=1 or 2',
c     &           ' - vx,vy not equal to zero with icase=3'
c            stop
c          endif
c        enddo
c      endif
c
c     9. test of physical data: pelas(1),pelas(6),itrel
c     -------------------------------------------------
      if ( icase.eq.4 .or. istres.eq.2 ) then
        do i=1,icoe
          if ( pelas(1,i).lt.0.d0 ) then
            write(iou,*) '0-shear rate relaxation time =',pelas(1,i),
     &                   ' < 0'
            istop = .true.
          endif
          if ( pelas(6,i).ne.0.d0 ) istres = 1
        enddo
        if ( .not.range(1,itrel,3) ) then
          write(iou,*) 'relaxation-time law  itrel=',itrel,
     &                 ' is out of range 1 - 3'
          istop = .true.
        endif
      endif
c
c     10. tests for transient problem
c     -------------------------------
      if ( trans )then
        if ( .not.range(1,iorder,3) ) then
          write(iou,*) 'iorder=',iorder,' is not allowed in the program'
          istop = .true.
        endif
        if ( endtim.le.0.d0 ) then
          write(iou,*) 'endtim must be strictly positive'
          istop = .true.
        endif
      else
        if ( rstart) then
          write(iou,*) 'steady problem rstart must be false'
          istop = .true.
        endif
      endif
c
c     11. check nvar and nz
c     ---------------------
      if ( .not.range(1,nvar,maxvar) ) then
        write(iou,*) 'number of variables nvar=',nvar,' must < ',maxvar
        istop =.true.
      endif
      if ( .not.range(1,nz,maxnz) ) then
        write(iou,*) 'size of array z  nz=',nz,' must < ',maxnz
        istop = .true.
      endif
c
c     12. stop of the program if istop is true
c     ----------------------------------------
      if ( istop ) then
        stop 'tesdat: error in checking input data, stop!'
      else
        write(iou,*) 'data from the data- and mesh-files are OK'
      endif
c
      return
      end
c
c
c***********************************************************************
      subroutine tesbnd ( iou,nbd,nic,ic,icase,istres,ipnode,isnode)
c
c     this function checks the boundary conditions
c***********************************************************************
      implicit none
      integer nbd,nic,ic(nic),icase,istres,ipnode,isnode,iou
c
      logical odd,range,impres,istop
      integer i,i1,i2,i3,icorn,istart,ins
c
c     0. functions
c     ------------
      odd(i) = ( ((i/2)*2).ne.i )
      range(i1,i2,i3) = ( i2.ge.i1 .and. i2.le.i3 )
c
c     1. initialization
c     -----------------
      istop = .false.
      impres = .true.
c
      write(iou,*) 'msg from tesbnd:'
c
c     2. check of isnode
c     ------------------
      if ( .not.(range(1,isnode,nbd).or.icase.eq.3.or.icase.eq.5) ) then
        write(iou,*) 'vanishing stream-function boundary node =',isnode,
     &               ' must < nbd=',nbd
        istop = .true.
      endif
c
c     3. check of the corners
c     -----------------------
      do i=1,nic
        icorn = ic(i)
        if ( .not.range(1,icorn,nbd) ) then
          write(iou,*) 'corner',i,' is ',icorn,', must < ',nbd
          istop = .true.
        endif
        if ( .not.odd(icorn) ) then
          write(iou,*) 'corner',i,' is ',icorn,', has to be odd'
          istop = .true.
        endif
      enddo
c
      istart = ic(1)
      if ( istart.ne.1 ) then
        write(iou,*) 'first corner is',istart,'. it must be 1'
        istop = .true.
      endif
      do i=2,nic
        if ( ic(i).le.istart ) then
          write(iou,1300) i-1,istart,i,ic(i)
1300      format('corner ',i2,' is ',i4,';  corner ',i2,' is ',i4/
     &           'the last one has to be greater than the first one')
          istop = .true.
        endif
        istart = ic(i)
      enddo
c
c     10. check of ipnode
c     -------------------
      if ( icase.ne.3 ) then
        if ( .not.impres ) then
          if ( ipnode.ne.0 ) write(iou,'(A,i4,A)')
     &             ' warning: vanishing pressure node ipnode=',ipnode,
     &             ' may not be imposed'
          if ( ipnode.eq.-1)
     &      write(*,*)'Dirichlet condition is imposed for some boundary
     &sections'
               
          if ( .not.range(1,ipnode,nbd).and.ipnode.ne.-1) then
            write(iou,*) 'vanishing pressure vertex ipnode=',ipnode,
     &                   ' must < nbd=',nbd
            istop = .true.
          endif
          if ( .not.odd(ipnode) ) then
            write(iou,*) 'vanishing pressure vertex ipnode =',ipnode,
     &                   ' must be an odd number'
            istop = .true.
          endif
        endif
      endif
c
      if ( istop ) then
        stop 'tesbnd: error in checking boundary data, stop!'
      else
        write(iou,*) 'boundary data are OK'
      endif
c
      return
      end
c
c
c***********************************************************************
      subroutine pridat ( iou, nnode,nvert,nbd,nelem,nvar,inod,ngmx,
     &                    istres, icase, ifluid,iax,
     &                    itmax,ivisc,itemp,itrel, iorder,iforce,icoe,
     &                    x,y,ro,grav, pvis,pelas, pterm, dtin, endtim,
     &                    zconv,trans,newton,inert,upwind,resold,
     &                    rstart,ivheat )
c
c     routine prints data
c***********************************************************************
      implicit none
      integer iou,nnode,nvert,nbd,nelem,nvar,
     &        istres,icase,ifluid,iax,itmax,
     &        ivisc,itemp,itrel,iorder,iforce,icoe,ngmx
      integer inod(ngmx,nelem)
      real*8 x(nnode),y(nnode),dtin,endtim,zconv
      double precision ro(icoe),grav(2,icoe),pvis(10,icoe),
     &                 pterm(10,icoe),pelas(10,icoe)
      logical trans,newton,inert,upwind,resold,rstart,ivheat
c
      integer ivi,imod,i
      logical grav0
c
c     1. problem information and data
c     -------------------------------
      write(iou,*) '                                '
      write(iou,*) '********************************'
      write(iou,*) '*                              *'
      write(iou,*) '*   p r o b l e m    d a t a   *'
      write(iou,*) '*                              *'
      write(iou,*) '********************************'
c
      if ( icase.eq.1 )  
     &  write(iou,*) 'isothermal flow problem            icase =',icase
      if ( icase.eq.2 )  
     &  write(iou,*) 'non isothermal flow problem        icase =',icase
      if ( icase.eq.3 ) 
     &  write(iou,*) 'thermal problem                    icase =',icase
      if ( icase.eq.4 ) 
     &  write(iou,*) 'isothermal viscoelastic flow       icase =',icase
      if ( icase.eq.5 ) 
     &  write(iou,*) 'stress postprocessor problem       icase =',icase
c
      if ( iax.eq.0 ) 
     &  write(iou,*) 'plane geometry                       iax =',iax
      if ( iax.eq.1 ) 
     &  write(iou,*) 'axisymmetric geometry                iax =',iax
      if ( iax.eq.2 ) 
     &  write(iou,*) 'swirling flow                        iax =',iax
c
      if ( trans ) then
        write(iou,*) '    ------------------ '
        write(iou,*) '    transient  problem '
        write(iou,*) '    ------------------ '
        if ( .not.(resold.or.rstart) ) write(iou,*) 'number of steady',
     &                 ' iterations before the transient problem ',itmax
        write(iou,'(A,1pe12.5)')
     &         ' final time is                             ',endtim
        write(iou,'(A,1pe12.5)')
     &         ' the value of dtin  is                     ',dtin
        if ( iorder.eq.1 ) write(iou,*) 
     &         'program uses implicit Euler method,   iorder =',iorder
        if ( iorder.eq.2 ) write(iou,*) 
     &         'program uses Galerkin method,         iorder =',iorder
        if ( iorder.eq.3 ) write(iou,*)
     &         'program uses Crank-Nicholson method,  iorder =',iorder
        if ( .not.(resold.or.rstart) ) write(iou,*) 
     &                              'initialization of problem at zero'
        if ( resold ) write(iou,*) 'old file contains only value of z'
        if ( rstart ) write(iou,*) 'old file contains value and ',
     &                             'derivates of z'
      else
        write(iou,*) '    ------------------ '
        write(iou,*) '    steady     problem '
        write(iou,*) '    ------------------ '
        write(iou,*) 'maximum number of iterations    itmax =',itmax
        write(iou,*) 'convergence test        zconv = ',zconv
        if ( upwind )      write(iou,*) 'streamline upwinding'
        if ( .not.upwind ) write(iou,*) 'no streamline upwinding'
        if ( resold )      write(iou,*) 'initialize with old results'
        if ( .not.resold ) write(iou,*) 'initialize at zero'
      endif
c
      if ( icase.ne.3 .and. icase.ne.5 ) then
        if ( iforce.eq.1) then
          write(iou,*) 'calculation of bd. forces       iforce =',iforce
        else
          write(iou,*) 'no calculation of bd. forces    iforce =',iforce
        endif
      endif
c
      if ( inert )      write(iou,*) 'inertia terms are accounted'
      if ( .not.inert ) write(iou,*) 'inertia terms are neglected'
c
      if ( .not.newton .and. (itemp.ne.1 .or. ivisc.ne.1) )
     &    write(iou,*) 'Picard iteration is requested for the viscosity'
c
c     2. physical data
c     ----------------
      write(iou,*) '                                '
      write(iou,*) '********************************'
      write(iou,*) '*                              *'
      write(iou,*) '*   p h y s i c a l   d a t a  *'
      write(iou,*) '*                              *'
      write(iou,*) '********************************'
c
      do i=1,icoe

      if (i.eq.1) then
         write(iou,*) 'For the continuous phase:'
      else
         write(iou,'(//A,i2,A)') 'For the dispersed phase',i-1,':'
      endif
c
      ivi = 0
      if ( pelas(6,i).ne.0.d0 ) ivi = 1
      if ( pelas(10,i).gt.0.d0 ) ivi = 2
      imod = 1
      if ( icase.eq.4 ) imod = 2*ifluid + ivi
      if ( ivi.lt.2 .and. imod.eq.1 ) then
        write(iou,*) 'generalized newtonian fluid'
        if ( istres.eq.1 ) write(iou,*) 
     &      ' with calculation of viscous stress        istres =',istres
        if ( istres.eq.2 ) write(iou,*) 
     &      ' with calculation of viscoelastic stress   istres =',istres
      elseif ( ivi.eq.2 ) then
        write(iou,*) 'linear elastic solid'
      endif
      if ( imod.eq.2 ) 
     &  write(iou,*) 'Maxwell/White-Metzner fluid       ifluid =',ifluid
      if ( imod.eq.3 ) 
     &  write(iou,*) 'Oldroyd fluid                     ifluid =',ifluid
      if ( imod.eq.4 )
     &  write(iou,*) 'Phan Thien-Tanner fluid           ifluid =',ifluid
      if ( imod.eq.5 ) 
     &  write(iou,*) 'Phan Thien-Tanner fluid with viscous comp.',
     &                                                ' ifluid =',ifluid
      if ( imod.eq.6 )
     &  write(iou,*) 'Giesekus fluid                    ifluid =',ifluid
      if ( imod.eq.7 )
     &  write(iou,*) 'Giesekus fluid with viscous comp. ifluid =',ifluid
c
      if ( icase.eq.2 ) then
        write(iou,*) 'dependence of viscosity on temperature:',
     &               '  itemp = ',itemp
        if (itemp.eq.1) write(iou,1200) pvis(1,i)
1200      format(' temperature independent viscosity  fac = fac0 =',
     &           e12.5)
        if (itemp.eq.2) write(iou,1201) pvis(1,i),pterm(1,i),pterm(2,i)
1201      format('  fac = fac0 * exp( -alfa * (t-talfa) )'/9x,
     &             'fac0 = ',e12.5,', alfa = ',e12.5,', talfa = ',e12.5)
        if (itemp.eq.3) write(iou,1202) pvis(1,i),pterm(1,i),pterm(2,i)
     &                                                      ,pterm(3,i)
1202       format('  fac = fac0 * exp( alfa/(t0+t) - alfa/(to+talfa) )'/
     &          7x,'fac0 = ',e12.5,', alfa = ',e12.5,', talfa = ',e12.5,
     &            ', t0 = ',e12.5)
      endif
c
      if ( icase.ne.3 ) then
        if ( ivi.le.1 ) then
        if( ivisc.gt.1 )
     &  write(iou,*) '-shear rate dependent viscosity:   ivisc = ',ivisc
        if( ivisc.eq.1 ) write(iou,1208) pvis(1,i)
1208      format(' -shear-rate independent viscosity visc = ',e12.5)
        elseif ( ivi.eq.2 ) then
          write(iou,'(A,e12.5)') ' -Youngs modulus=',pelas(1,i)
          write(iou,'(A,e12.5)') ' -Poissons ratio=',pelas(2,i)
        endif
c
        if( ivisc.eq.2 ) write(iou,1203) pvis(1,i),pvis(2,i),pvis(3,i),
     &                                   pvis(4,i)
1203        format('   viscosity law : bird carreau :'/
     &           '   visc = facinf + (fac-facinf) * '/17x,
     &           '(1 +tnat*tnat*gamma*gamma)**((expo-1)/2)'/
     &           '   fac = ',e12.5,', tnat = ',e12.5,', expo = ',e12.5,
     &           ', facinf = ',e12.5)
        if( ivisc.eq.3 ) write(iou,1204) pvis(1,i),pvis(3,i)
1204        format('   viscosity law : power law :'/
     &             '   visc = fac * gamma**(expo-1)'/
     &             '        fac = ',1pe12.5,' ,expo = ',1pe12.5)
        if( ivisc.eq.4 ) write(iou,1205) pvis(1,i),pvis(2,i),pvis(3,i)
1205        format('   viscosity law : bingham'/3x,
     &             'visc = fac + ystr/gamma',17x,'if gamma > gcrit'/8x,
     &           '= fac + ystr/gcrit*(2-gamma/gcrit)  if gamma < gcrit'/
     &           6x,'fac = ',e12.5,', ystr = ',e12.5,', gcrit = ',e12.5)
        if( ivisc.eq.5 ) write(iou,1206) pvis(1,i),pvis(2,i),pvis(3,i),
     &                                   pvis(4,i)
1206        format('   viscosity law : herschel-bulkley'/3x,
     &             'gamma = gcrit',27x,'if gamma < gcrit'/
     &             '   visc  = fac / gamma + tnat*gamma**(expo-1)'/3x,
     &             'fac = ',e12.5,', tnat = ',e12.5,', expo = ',e12.5,
     &             ', gcrit = ',e12.5)
        if( ivisc.eq.6 ) write(iou,1207) pvis(1,i),pvis(2,i),pvis(3,i)
1207        format('   viscosity law : cross'/3x,'visc = fac /',
     &             '(1 + (tnat*gamma)**expom)'/3x,'fac = ',e12.5,
     &             ', tnat = ',e12.5,', expom = ',e12.5)
c
        if ( imod.eq.3 .and. pelas(6,i).ne.0.d0 ) 
     &                                        write(iou,1210) pelas(6,i)
1210         format(/' -viscous component for the extra-stress :'/
     &               '    visc = visc*(1-ratio) + visc*ratio '/
     &               '    visc = visc1          + visc2      '/
     &               '    visc2 = viscous component for extra-stress'/
     &               '    ratio = ',e14.7)
        if ( imod.eq.4 .or.imod.eq.5 ) 
     &                             write(iou,1211) pelas(4,i),pelas(5,i)
1211         format(/' -Phan Thien-Tanner specific coefficients : '/
     &            '    eps = ',e14.7,'    xi = ',e14.7)
        if ( imod.eq.6 .or.imod.eq.7 ) write(iou,1212) pelas(4,i)
1212         format(/' -Giesekus specific coefficients: alpha = ',e14.7)
c
        write(iou,1214) ro(i)
1214    format(' -density ro =',e14.5)
        grav0 = .false.
        if ( grav(1,i).eq.0.d0 .and. grav(2,i).eq.0.d0 ) grav0 = .true.
        if ( grav0 )      write(iou,*) '-gravity field is neglected'
        if ( .not.grav0 ) write(iou,1213) grav(1,i),grav(2,i)
1213      format(' -gravity field is accounted:  gx = ',e14.5,
     &           ' gy =',e14.5)
      endif
c
      if ( icase.eq.2 .or. icase.eq.3 ) then
        if ( icase.eq.3 .and. .not.grav0 ) write(iou,*)
     &  '-convection velocity is imposed: u =',grav(1,i),' v=',grav(2,i)
        write(iou,*) '-average initial temperature within the domain ',
     &               'is:  tinit = ',pterm(10,i)
c
        if ( icase.eq.3 ) goto 90
        if ( ivheat )      write(iou,*) '-viscous heating is accounted'
        if ( .not.ivheat ) write(iou,*) '-viscous heating is neglected'
      endif
c
      if ( icase.eq.4 .or. istres.eq.2 ) then
      write(iou,*) '-dependence of relaxation time on shear rate:',
     &             ' itrel =',itrel
      if ( itrel.eq.1 ) write(iou,1230) pelas(1,i)
1230    format('  shear-rate independent relaxation time: ',
     &          'trelax = facr =',e12.5)
      if ( itrel.eq.2 ) write(iou,1231) pelas(1,i),pelas(2,i),pelas(3,i)
1231       format('  relaxation law : bird carreau :'/
     &        '   trelax = facr*(1 + tnatr**2*gamma**2)**((expor-1)/2)'/
     &       '   trelax = ',e12.5,', tnatr = ',e12.5,', expor = ',e12.5)
      if ( itrel.eq.3 ) write(iou,1232) pelas(1,i),pelas(3,i)
1232       format('  relaxation law : power law :'/
     &            '   trelax = facr * gamma**(expor-1)'/
     &            '   trelax = ',e12.5,', expor = ',e12.5)
      endif
c
      enddo
c
c     3. mesh data
c     ------------
90    write(iou,*) '                                '
      write(iou,*) '********************************'
      write(iou,*) '*                              *'
      write(iou,*) '*      m e s h    d a t a      *'
      write(iou,*) '*                              *'
      write(iou,*) '********************************'
      write(iou,*) 'number of elements,        nelem = ',nelem
      write(iou,*) 'number of nodes,           nnode = ',nnode
      write(iou,*) 'number of vertices,        nvert = ',nvert
      write(iou,*) 'number of boundary nodes,  nbd   = ',nbd
      write(iou,*) 'number of variables,       nvar  = ',nvar
c
      return
      end
c
c
c***********************************************************************
      subroutine pribnd ( iou,bdpr,nbd,icase,istres,nic,ipnode,
     &                    isnode,ic,ibdu,ibdt,ibds,idp,
     &                    densu,denst,denss,ncpu,ncps,valdp )
c
c     this routine prints boundary conditions
c***********************************************************************
      implicit none
      integer iou,nbd,icase,istres,ipnode,isnode,nic,idp, ncpu,ncps
      integer ic(nic+1),ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd)
      real*8 densu(ncpu,nbd),denst(nbd),denss(ncps,nbd),valdp
      logical bdpr
c
c     local variables
      integer iside,ipoint,npoint,ica,icb,is,i,k,kk,itest,j
c
c     1. initialization
c     -----------------
      if ( .not.bdpr .or. icase.eq.5 .and. istres.eq.1 ) return
c
      write(iou,*) ' '
      write(iou,*) '***********************'
      write(iou,*) '*                     *'
      write(iou,*) '* boundary conditions *'
      write(iou,*) '*                     *'
      write(iou,*) '***********************'
c
      if ( ipnode.ne.0 ) write(iou,*) 'pressure vanishes at ',
     &                                'boundary node number=',ipnode
      if ( isnode.ne.0 ) write(iou,*) 'stream-function vanishes at ',
     &                                'boundary node number=',isnode
c
c     2. boundary conditions
c     ----------------------
      do 20 iside=1,nic
c
      ica = ic(iside)
      icb = ic(iside+1)-1
      write(iou,*) ' '
      write(iou,*) '**********************************'
      write(iou,1201) iside,ica,icb
1201  format(' * side :',i2,'   corners :',i4,' -',i4,' *')
      write(iou,*) '**********************************'
c
c     5. print boundary conditions
c     ----------------------------
      write(iou,1500) 
1500  format(/78('-')/'|bd node | ibd    densu   | ibd    densv   |',
     &              ' ibd    densw   | ibdt    denst  |'/78('-'))
      npoint = icb - ica + 1
      do ipoint=1,npoint
        k = ica + ipoint - 1
        if ( k.gt.nbd ) k=1
        write(iou,1504) k,(ibdu(j,k),densu(j,k),j=1,ncpu),
     &                  ibdt(k),denst(k)
1504    format('|  ',i4,'  |',4(i2,' ',1pe12.5,' |'))
      enddo
      write(iou,1505)
1505  format(78('-'))
c
c     6. stresses
c     -----------
      if ( icase.ne.4.and.istres.ne.2 ) goto 20
c
      write(iou,1600) 
1600  format(/78('-')/'|bd node|ibds|     denss1    |     denss2',
     &       '    |     denss3    |     denss4    |'/78('-'))
c
      npoint = icb - ica + 1
      do ipoint=1,npoint
        k = ica + ipoint - 1
        if ( k.gt.nbd ) k=1
        write(iou,1603) k,ibds(k),(denss(j,k),j=1,ncps)
1603    format('| ',i4,'  |',i3,' |',4(1pe14.7,' |'))
      enddo
      write(iou,1605)
1605  format(78('-'))
c
20    continue
c
      return
      end
