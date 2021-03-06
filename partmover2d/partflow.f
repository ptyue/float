      program PARTFLOW
c***********************************************************************
c     This program calculates motion of particles in Newtonian
c       and viscoelastic fluids
c
c     by Howard H. Hu,   April, 1991
c
c     modified:
c              Jan.26, 1995: periodic boundary in x-direction
c              Mar.24, 1995: New "ALE" scheme, semi-implicit
c              May 2, 1995: LOCATC updated & UPDATE modified to
c                           handle particle collision
c              Nov.21, 1995: mshchk is added (checking mesh quality)
c              Feb.19, 1996: second order (in time) scheme
c                            least-square projection
c              June 24, 1996: periodic and non-periodic codes are merged
c              Oct.20, 1996: new method of implementing boundary cond.
c              Oct, 1998: new flowtype = 11/12 is added
c              Dec, 1998: constant pressure element is added
c     Modified By Pengtao Yue,     Jan, 2005
c         Cahn-Hillard equation added
c              Feb 23, 2005: calculating drop deformation is added
c***********************************************************************
      implicit none
      include '../include/parameter2d.h'
c
c     Definitions
      integer nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdstrm,nrdmsh,
     &        ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgstrm,ndgmsh
c
      integer ioptu,imu,lfilu,maxitu, iopts,lfils,maxits,
     &        ioptm,lfilm,maxitm, iopti,lfili,maxiti
     &        ,ioptp,maxitp
      real*8 epsu,epss,epsm,epsi,epsp
c
      logical lrmsh,lfilter
      integer nnode,nvert,nelem,inod(ngmx,maxelt),nec(ncmx,maxelt),
     &        nnodd,nvrtd,nelmd,inodd(ngmx,maxelt),necold(ncmx,maxelt),
     &        nbd,ibdnod(maxnbd),nic,ic(maxbdp),nbound,nside(maxbdp),
     &        nbdd,ibddd(maxnbd),icold(maxbdp)
      real*8 x(maxnod),y(maxnod),area(maxelt),aspr(maxelt),
     &       xold(maxnod),yold(maxnod)
      integer nvar,neq, norder
      integer ncpvelo,ncppres,ncptemp,ncpelas,
     &        indgu,indgp,indgtm,indgr,indgnb,indgstrm
      logical SLVvelo,SLVpres,SLVtemp,SLVelas
      real*8 u(ncpg,maxnod),p(maxvrt),umesh(ncpg,maxnod),
     &       dudt(ncpg,maxnod),dudx(ncpg,maxnod),
     &       dudy(ncpg,maxnod),vort(maxvrt),strm(maxvrt),
     &       selas(ncps,maxvrt),desdt(ncps,maxvrt),
     &       uold(ncpg,maxnod),selasd(ncps,maxvrt)
      real*8 xpos(ncpb,maxptc),fpart(ncpb,maxptc),amas(ncpb,maxptc),
     &       uvom(ncpb,maxptc),fxym(ncpb,maxptc),xpod(ncpb,maxptc),
     &       uvod1(ncpb,maxptc),duvot(ncpb,maxptc)
     &      ,fbody(ncpb,maxptc),fcols(ncpb,maxptc)
      logical bdrefn,ifixpt(ncpb)
      integer nbrigid,nbfluid,nbelast,nbtotal,
     &        nbshp(maxptc), ncr,ny, flowtype
      real*8 diaa(maxptc),diab(maxptc)
      real*8 dt1,dt,dtmin,dtmax,dtrate,disinc,velinc,rinc, time
      real*8 ueta(4,2*maxnth+1)
      integer itime,itordr,ntime,iprint,nsavstt,nptim,nnew1,nnew2,npfmt
      integer ibdu(ncpu,maxnbd),ibdt(maxnbd),ibds(maxnbd),ibdp(maxnbd)
      real*8 densu(ncpu,maxnbd),denst(maxnbd),denss(ncps,maxnbd),
     &       densp(maxnbd)
      real*8 z(maxvar),zs(maxvar),dtin,endtim,tsc,tsc1
      real*8 fxnbd(maxnbd),fynbd(maxnbd)
      integer nbd0,nic0,nbound0
c*CH     start of Definitions related with Cahn-Hillard Equation
      real*8 phi(maxnod),phiold(maxnod),psi(maxnod),dphidt(maxnod)
      integer ncpphi,nrdphi,ndgphi
      integer indgphi, indgpsi
      logical SLVphi
	real*8 phga,pheps,phlamb,phsft
c*CH     end of Definitions related with Cahn-Hillard Equation
c       Gmsh
      real(8) dist(maxvrt),domsize
c     for initialization
      logical noflow
      integer ncpveln,ncpelan
c     drop deformation
      integer ishm
      logical dropone
c
      integer icollision, bdseg
      real*8 pvelbd(10),pwall(5,maxbdp-maxptc),epscollision
c
      logical noexec,inert,ivheat,newton,upwind
      logical meshpr,bdpr,resold,rstart,trans
      integer icase,iax,ifluid,ivisc(maxcoe),itemp,itrel,icoe
      integer iforce,istres,ipnode,isnode,itmax,iorder
      double precision ro(maxcoe),grav(2,maxcoe),pvis(10,maxcoe),
     &                 ptemp(10,maxcoe),pelas(10,maxcoe),zconv
c     moving contact line
      real(8) wallrela(maxbdp),wallener(maxbdp)
c     normal stress boundary conditions
      integer ibdn(maxnbd)
      real(8) densn(maxnbd),bdpre(maxbdp)

c
      integer, allocatable :: iwork(:)
      real*8,  allocatable :: rwork(:)
      integer  ierr
c      common /work/ rwork
c      common /iwork/ iwork
c
      integer reft(maxelt)
      character*10 outfile
      real*8 xcv(maxnbd),ycv(maxnbd),strsx(maxnbd),strsy(maxnbd)
      integer ndglb,lenstr,nk,k, i
      integer io5,io6
c     io5 for temperary file(open and closed with the subprograms)
c     io6 is open over the runtime of program
      data io5/25/,io6/24/ !(note that 26 is used by other files)
c      data io5/10/,io6/11/
      real*8 beta
      data beta/0.d0/
c
      logical xperid, mshslpp,mshslpw
      real*8 xlength
      common /xperid/ xlength,xperid
      common /mesh/ nnode,nvert,nelem
      logical lfrmt,local
c     mesh size control
      real(8) h1,h2,h3,dG
c     time step control
      real(8) cflflow,cflphi
c     rescaling factors
      real(8) refr,refu,refp
c     reset extra stress
      integer irstr
c
      character filenum*5, filename*11
      logical :: lwrmesh=.false.
c   meshgen =1, grummp; 2, gmsh      
      integer :: meshgen=2
      
      

c
C     Support for Petsc timing
c      include '../include/petsc.h'
c
c      integer ierr
c
c      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
c
C     Create a new user-defined event.
C      - Note that PLogEventRegister() returns to the user a unique
C        integer event number, which should then be used for profiling
C        the event via PLogEventBegin() and PLogEventEnd().
C      - The user can also optionally log floating point operations
C        with the routine PLogFlops().
C
c      call PLogEventRegister(PARTMOVER_2D_EVENT,'-PartMover2d-  ',
c     * 'Red:',ierr)
c      call PLogEventRegister(SOLVER_EVENT,    'Linear Solver  ',
c     * 'Red:',ierr)
c      call PLogEventRegister(ILU0_EVENT,      ' ILU(0) Factor ',
c     * 'Red:',ierr)
c      call PLogEventRegister(ILUT_EVENT,      ' ILUT Factor   ',
c     * 'Red:',ierr)
c      call PLogEventRegister(ILUTP_EVENT,     ' ILUTP Factor  ',
c     * 'Red:',ierr)
c      call PLogEventRegister(GMRES_EVENT,     ' GMRES         ',
c     * 'Red:',ierr)
c      call PLogEventRegister(BICG_EVENT,      ' BICG          ',
c     * 'Red:',ierr)
c      call PLogEventRegister(CG_EVENT,        ' CG            ',
c     * 'Red:',ierr)
c      call PLogEventRegister(MATMULT_EVENT,   '   MatMult     ',
c     * 'Red:',ierr)
c      call PLogEventRegister(LUSOL_EVENT,     '   LU Solve    ',
c     * 'Red:',ierr)
c      call PLogEventRegister(ELEM_MAT_EVENT,  'Elem Mat Gen   ',
c     * 'Red:',ierr)
c      call PLogEventRegister(MESH_MOVE_EVENT, 'Mesh Movement  ',
c     * 'Red:',ierr)
c      call PLogEventRegister(PRE_MAT_EVENT,   'Prepare Matrix ',
c     * 'Red:',ierr)
c      call PLogEventRegister(REMESH_EVENT,    'Mesh Generation',
c     * 'Red:',ierr)
c      call PLogEventRegister(PROJECTION_EVENT,'Projection     ',
c     * 'Red:',ierr)
c      call PLogEventRegister(MAT_PROJ_EVENT,  ' Mat&Interpolat',
c     * 'Red:',ierr)
c      call PLogEventRegister(SOLV_PROJ_EVENT, ' Solver        ',
c     * 'Red:',ierr)
c
c      call PLogEventBegin(PARTMOVER_2D_EVENT,
c     *                    PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c
      bdrefn = .true.
      local  = .false.
c
c     order of interpolation functions
c     --------------------------------
      nrdgeom = 2
      nrdvelo = 2
      nrdpres = 1
      nrdelas = 1
      nrdtemp = 1
      nrdstrm = 1
      nrdmsh  = 1
      nrdphi  = 2
c
      allocate (iwork(maxiwk), STAT=ierr)
      if(ierr.ne.0)then
         write(*,*)'allocating iwork error, STAT=',ierr
         stop
      endif
      allocate (rwork(maxrwk), STAT=ierr)
      if(ierr.ne.0)then
         write(*,*)'allocating rwork error, STAT=',ierr
         stop
      endif


c
c-----------------------------------------------------------------------
c     Initialization
      call INITIA(ntime,nptim,time,iprint,nsavstt,norder,
     &            maxvrt,maxnod,maxelt,maxnbd,maxbdp,maxptc,maxnth,
     &            maxcoe,ngmx,ncmx,ncpg,ncps,ncpb,
     &            nvert,nnode,nelem,inod,x,y,reft,nec,area,aspr,
     &            nbd,ibdnod,nic,ic,nbound,nside,
     &            u,p,dudt,umesh,selas,desdt,phi,psi,dphidt,
     &            nbrigid,xpos,uvom,duvot,fxym,amas,fpart,ifixpt,
     &            nbfluid,nbelast,nbtotal,icollision,epscollision,
     &            nbshp,diaa,diab,ncr,ny,bdrefn,
     &            flowtype,pwall,bdseg,
     &            noexec,inert,ivheat,newton,upwind,lfrmt,
     &            meshpr,bdpr,resold,rstart,trans,
     &            icase,iax,ifluid,ivisc,itemp,itrel,icoe,
     &            iforce,istres,ipnode,isnode,itmax,iorder,irstr,
     &            ro,grav,pvis,ptemp,pelas,zconv,
     &            pvelbd,ueta, nnew1,nnew2,
     &            dt,dtmin,dtmax,dtrate,cflflow,cflphi,
     &               disinc,velinc,rinc,
     &            phga,phlamb,pheps,phsft,
     &            ioptu,imu,lfilu,maxitu,epsu,lfilter,
     &            iopts,lfils,maxits,epss,
     &            ioptm,lfilm,maxitm,epsm,
     &            iopti,lfili,maxiti,epsi,
     &            ioptp,maxitp,epsp,
     &           nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdstrm,nrdmsh,
     &           nrdphi,
     &           ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgstrm,ndgmsh,
     &           ndgphi,
     &            SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &            ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &            iwork,maxiwk,rwork,maxrwk, io5,
     &            noflow,ncpveln,ncpelan,dropone,
     &            h1,h2,h3,dG,
     &            wallrela,wallener,bdpre,
     &            refr,refu,refp,
     &            dist,domsize,meshgen)

c
      write(*,'(A,i4,A,i4)')
     &           'Ntime=',ntime,'   starting from last itime=',nptim
c
      mshslpw = .true.
      if ( flowtype.eq.3 .or. flowtype.eq.4 ) mshslpw = .false.
      if ( flowtype.eq.11 .or. flowtype.eq.12 ) mshslpw = .false.
      if ( flowtype.eq.99 ) mshslpw = .false.
c
      mshslpp = .true.
      do k=1,nbrigid
        if ( nbshp(k).ne.1 .or. dabs(diaa(k)/diab(k)-1.d0).gt.1.d-5 )
     &     mshslpp = .false.
      enddo
c     output the initial data
      k=ncpelas
      if (nptim==0)
     &  call tecpost ( 0,time, nbrigid,xpos, xcv,ycv,
     &                 strsx,strsy,nnode,nvert,nelem,inod,x,y,
     &                 nbd,ibdnod,nic,ic, nbound,nside,
     &                 bdseg,pwall,
     &                 ngmx,ncpg,ncps,ncpb,ncpvelo,k,ncpphi,
     &                 u,p,strm,vort,selas,phi,psi,io5 )
c      pause

c     set mesh velocity and acceleration to zero for flowtype/100=1
      if(flowtype/100.eq.1) umesh(1:ncpg,1:maxnod)=0.d0
c
c-----------------------------------------------------------------------
c     Main loop, for 'ntime' time steps
c-----------------------------------------------------------------------
c
      ishm=0 !for triangular mesh
      do itime=nptim+1,nptim+ntime
c     debug
c	u=1.
c     end debug
c
        close ( io6 )
        open ( io6, file = 'hh.lst' )
c
c       order of time step schemes
        itordr = 1
        if ( itime.ge.norder ) itordr = iorder
        write(*,'(A)')
     +  '*************************************************************'
c
c       Update mesh for the current time
c	pause 'in update'
      if(flowtype/100.eq.0)then
        call UPDATE ( itime, itordr, dt,dt1,dtmin,dtmax,dtrate,
     &                nnode,nelem,inod,x,y,area,aspr,ngmx,reft,
     &                ncpg,xold,yold,umesh,lrmsh,
     &                nbd,ibdnod,nic,ic,mshslpp,mshslpw,ifixpt,
     &                nbrigid,nbfluid,nbelast,ncpb,xpos,xpod,
     &                 uvom,uvod1,fxym,amas,fpart,fbody,fcols,
     &                ioptm,lfilm,maxitm,epsm,
     &                diaa,diab, disinc,velinc,rinc,ncr,
     &                icollision,epscollision,flowtype,pwall,bdseg,
     &                nrdmsh,nrdgeom,nrdvelo,
     &                iwork,maxiwk,rwork,maxrwk,io6)
c       print*, 'fxym =',(fxym(k,1)*amas(k,1),k=1,3)
c       print*, 'fpart=',fpart(1,1),fpart(2,1),fpart(3,1)
c       print*, 'fbody=',fbody(1,1),fbody(2,1),fbody(3,1)
c
c       pause 'after update'
c
        if ( flowtype.eq.4 ) then
           call height(itime,time,nbrigid,xpos,pwall(4,2))
        endif
c
        if ( meshpr ) call PRTMSH ( 'mesh.old',6,6, nvert,nnode,nelem,
     &          inod,x,y,reft,nbd,ibdnod,nic,ic, nbound,nside,nec, 20)
        if ( meshpr .and. itime.gt.1 ) lrmsh = .true.
c
        if ( lrmsh ) then
c
c          Save old mesh
           call MSHOLD ( nnode,nvert,nelem,inod,x,y,nec, ngmx,ncmx,
     &                   nbd,ibdnod,nic,ic,
     &                   nnodd,nvrtd,nelmd,inodd,xold,yold,necold,
     &                   nbdd,ibddd,icold )
c
c          Rremesh if the mesh is too distorted
c	pause 'in rmesh'
           call RMESH (itime,flowtype,pwall,bdseg,bdrefn,ncr,ny,
     &                 nbrigid,nbfluid,nbelast,nbtotal,
     &                 nbshp,diaa,diab,xpos,
     &                 nvert,nnode,nelem,inod,x,y,reft,nec,area,aspr,
     &                 nbd,ibdnod,nic,ic,nbound,nside,
     &                 maxvrt,maxnod,maxelt,maxnbd,maxbdp,
     &                 iwork,maxiwk,rwork,maxrwk )
c	pause 'after rmesh'
           if ( meshpr ) call PRTMSH ('mesh.new',6,6, nvert,nnode,nelem,
     &          inod,x,y,reft,nbd,ibdnod,nic,ic, nbound,nside,nec, 20)
c
        endif
      elseif(flowtype/100.eq.1)then
        call timestep(itime,ishm,dt,dt1,dtmin,dtmax,dtrate,
     &                   cflflow,cflphi,
     &                   nelem,ngmx,ncpg,inod,
     &                   nrdgeom,nrdvelo,nrdphi,ndggeom,ndgvelo,ndgphi,
     &                   x,y,u,phi)

c           call timestep(itime-nptim,dt,dt1,dtmin,dtmax,dtrate)
c          Save old mesh
        if(meshgen==1)then
          call MSHOLD ( nnode,nvert,nelem,inod,x,y,nec, ngmx,ncmx,
     &                   nbd,ibdnod,nic,ic,
     &                   nnodd,nvrtd,nelmd,inodd,xold,yold,necold,
     &                   nbdd,ibddd,icold )
c           lrmsh=.false.
          call GRrmesh( maxvrt,maxnod,maxelt,maxnbd,maxbdp,
     &                   ngmx,ncmx,nvrtd,phi,
     &                   nvert,nnode,x,y,
     &                   nelem,inod,nec,area,reft,aspr,
     &                   nbd,ibdnod,nic,ic,nbound,nside,
     &                   maxiwk,iwork,h1,h2,h3,dG,lrmsh)
        elseif(meshgen==2)then
          call gmshlrmsh(nvert,phi,dist,h1,lrmsh)
          if(lrmsh)then
            call MSHOLD ( nnode,nvert,nelem,inod,x,y,nec, ngmx,ncmx,
     &                   nbd,ibdnod,nic,ic,
     &                   nnodd,nvrtd,nelmd,inodd,xold,yold,necold,
     &                   nbdd,ibddd,icold )          
            call gmshrmsh(maxvrt,maxnod,maxelt,maxnbd,maxbdp,ngmx,ncmx,
     &                    nvrtd,nnodd,nelmd,
     &                    inodd,necold,xold,yold,phi,dist,
     &                    nvert,nnode,nelem,inod,nec,x,y,
     &                    nbd,ibdnod,nic,ic,nbound,nside,
     &                    area,reft,aspr,
     &                    h1,h2,h3,dG,domsize,
     &                    maxiwk,maxrwk,iwork,rwork)            
          endif
        endif
c	debug
c	lrmsh=.true.
c
      endif
c
        if(lrmsh)then
c      write(*,*)'old mesh',nnodd,nelmd,nbdd
c      write(*,*)'new mesh',nnode,nelem,nbd
           lwrmesh = .true.
           ndgvelo = ndglb(nrdvelo)
           ndgpres = ndglb(nrdpres)
           ndgelas = ndglb(nrdelas)
           ndggeom = ndglb(nrdgeom)
           ndgtemp = ndglb(nrdtemp)
           ndgstrm = ndglb(nrdstrm)
           ndgmsh  = ndglb(nrdmsh)
           ndgphi  = ndglb(nrdphi)

c	print*,ndgvelo,ndgpres,ndgelas,ndggeom,ndgtemp,ndgstrm,ndgmsh
c
c          Project the flow field
           call INTPLT ( nrdgeom,nrdvelo,nrdpres,nrdelas,nrdphi,
     &                   ndggeom,ndgvelo,ndgpres,ndgelas,ndgphi,
     &                   maxnod,maxvrt,
     &                   nnode,nelem,inod,ngmx,x,y,
     &                   nbd,ibdnod,nbound,nic,ic,nside, nbelast,
     &                   nnodd,nvrtd,nelmd,inodd,xold,yold,necold,ncmx,
     &                   nbdd,ibddd,icold,
     &                   u,p,dudt,ncpelas,selas,desdt,ncpg,ncps,
     &                   ncpphi,phi,psi,dphidt,
     &                   iopti,lfili,maxiti,epsi, local,beta,
     &                   iwork,maxiwk,rwork,maxrwk,io6 )
           call epc(nvert, nnode, nelem, ngmx, 10*h1,
     &              inod, phi, x, y, dist,
     &              iwork,rwork,maxiwk,maxrwk)
        endif
c
        time = time + dt
        write(*,'(A,i4,A,1pe11.4,A,1pe11.4,A,i2,A,l2)')
     &          'itime=',itime,' time=',time,' dt=',dt,
     &          ' itordr=',itordr,' rmesh=',lrmsh
c
c       Apply boundary conditions
        call BOUNDC ( flowtype,iax,pwall,bdseg,nnode,x,y,nbd,ibdnod,nic,
     &                  ic,
     &                nbrigid,nbfluid,nbelast,xpos,uvom,ncpb,ncpu,ncps,
     &                icase,pvis,pelas, pvelbd,ueta,ivisc(1),
     &                maxnth,ibdu,ibdt,ibds,densu,denst,denss,
     &                ibdn,densn,bdpre,
     &                ndgphi,phi,time,ipnode,ibdp,densp)

c
c       Prepare input for flow solver
        call PREFLW ( itime, itordr,dt,dt1,dtin,endtim,tsc,tsc1,
     &                ndgvelo,ndgpres,ndgtemp,ndgelas,ndgphi,
     &                ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                 indgstrm,
     &                z,zs, u,uold,p,dudt,ncpg,ncps,ncpb,
     &                selas,desdt,selasd,phi,phiold,psi,dphidt,
     &                nbrigid,uvom,uvod1,duvot )
c
c       Solve flow field
        nbound0 = nbound - nbfluid - nbelast
        nic0    = nic    - nbfluid - nbelast
        nbd0    = ic(nic0+1) - 1
        call pflow ( io6, tsc,tsc1,dtin,endtim,flowtype,
     &               neq,nvar, mshslpw,mshslpp, nnew1,nnew2,
     &               nvert,nnode,nelem,inod,x,y,reft,
     &               nbd,ibdnod,nic,ic,nbound,nside,
     &               nbd0,nic0,nbound0,
     &               maxnod,maxelt,maxnbd,maxbdp,maxcoe,maxvar,maxnz,
     &               nbrigid,nbfluid,nbelast,
     &               xpos,fxym,amas,fbody,ifixpt,
     &               ioptu,imu,lfilu,maxitu,epsu,lfilter,
c     &               ioptphi,imphi,lfilphi,maxitphi,epsphi,lfilterphi,
     &               iopts,lfils,maxits,epss,
     &               ioptm,lfilm,maxitm,epsm, ioptp,maxitp,epsp,
     &               nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdstrm,
     &                nrdmsh,nrdphi,
     &               ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgstrm,
     &                ndgmsh,ndgphi,
     &               SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &               ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &               indgu,indgp,indgtm,indgr,indgnb,indgphi,indgpsi,
     &                indgstrm,
     &               ngmx,ncpu,ncpg,ncps,ncpb,
     &               fxnbd,fynbd,umesh,z,zs, selas,dt,
     &               ibdu,ibdt,ibds,ibdp,ibdn,
     &               densu,denst,denss,densp,densn,
     &               noexec,inert,ivheat,newton,upwind,
     &               bdpr,resold,rstart,trans,
     &               icase,iax,ifluid,ivisc,itemp,itrel,icoe,
     &               iforce,istres,ipnode,isnode,itmax,itordr,
     &               ro,grav, pvis,ptemp,pelas,zconv,
     &               phga,phlamb,pheps,phsft,
     &               wallrela,wallener,
     &               iwork,maxiwk,rwork,maxrwk,
     &               refr,refu,refp,dropone)
c
c       Convert results from flow solver
        call PSTFLW ( itordr, dt,tsc, z,zs,
     &                ndgvelo,ndgpres,ndgelas,ndgstrm,
     &                indgu,indgp,indgr,indgnb,indgstrm,
     &                ncpvelo,ncpelas,ncppres,u,uold,p,dudt,
     &                ncpg,ncps,ncpb,selas,desdt,selasd,strm,
     &                nbrigid,uvom,duvot,uvod1,
     &                ncpphi,ndgphi,indgphi,indgpsi,phi,phiold,psi,
     &                  dphidt)
c
c     reset the extra stress in the Newtonian region
      if(icase.eq.4.and.mod(itime,irstr).eq.0)then
        call rstextra (ndgelas,ndgphi,ncpelas,ncpphi,ncps,
     &                 selas,desdt,phi)
      endif

c       Calculate velocity gradient, vorticity and stress field
        call VORTCT ( vort, u,dudx,dudy,ndgvelo,ngmx,ncpg,
     &                nbrigid,nbfluid,nbelast,
     &                nvert,nnode,nelem,inod,reft,x,y,
     &                selas,ncps,ndgelas,pelas,icoe,dt,
     &                nbd,ibdnod,xcv,ycv,strsx,strsy, iwork,rwork )
c     output the drop deformation info
        if(flowtype.eq.102)then
        call DropDf(time,io5,io6,ishm,iax,nvert,nelem,ngmx,ncmx,
     &              inod,nec,nrdgeom,nrdphi,ndggeom,ndgphi,x,y,phi,
     &              iwork,rwork,dropone)
c     extension flow
        elseif(flowtype.eq.103.or.flowtype.eq.104.or.flowtype.eq.113.or.
     &         flowtype.eq.114)then
        call DropDfExt(time,io5,nrdgeom,ndggeom,x,y,nrdphi,ndgphi,phi,
     &                 nic,ic,nbd,ibdnod)
        elseif(flowtype.eq.111)then
        call DropLenWid(time,io5,iax,ngmx,nelem,ndggeom,
     &                  inod,x,y,nrdphi,ndgphi,phi)
c
        elseif(flowtype.eq.124)then
        call dropspread(time,itime,iprint,io5,io6,ishm,
     &                  nvert,nelem,ngmx,ncmx,
     &                  inod,nec,nrdgeom,nrdphi,ndggeom,ndgphi,x,y,phi,
     &                  nic,ic,nbd,ibdnod,iwork,rwork)
        endif
c
c       Write results in data file
        npfmt = max0(iprint,1)
        if ( mod(itime,npfmt).eq.0 ) then
c
           k = ncpelas
           if ( nbelast.gt.0 ) k = max(ncpelas,3)
c
           if ( .not.lfrmt ) then
              stop 'lfrmt must be .true. !'
              call cvnust(itime,5,outfile,lenstr)
              outfile = outfile(1:lenstr)//'\0'

c              call WRESLT ( itime,time, nbrigid,xpos, xcv,ycv,
c     &                      strsx,strsy,nnode,nvert,nelem,inod,x,y,
c     &                      nbd,ibdnod,nic,ic,nbound,nside,
c     &                      bdseg,pwall,2,k,
c     &                      u,p,strm,vort,selas, outfile )
           else
c              call WRESLT2 ( itime,time, nbrigid,xpos, xcv,ycv,
c     &                       strsx,strsy,nnode,nvert,nelem,inod,x,y,
c     &                       nbd,ibdnod,nic,ic, nbound,nside,
c     &                       bdseg,pwall,
c     &                       ngmx,ncpg,ncps,ncpb,ncpvelo,k,
c     &                       u,p,strm,vort,selas, io5 )
              call tecpost ( itime,time, nbrigid,xpos, xcv,ycv,
     &                       strsx,strsy,nnode,nvert,nelem,inod,x,y,
     &                       nbd,ibdnod,nic,ic, nbound,nside,
     &                       bdseg,pwall,
     &                       ngmx,ncpg,ncps,ncpb,ncpvelo,k,ncpphi,
     &                       u,p,strm,vort,selas,phi,psi,io5 )

           endif
        endif
c
c       Save information for restarting
        if ( mod(itime,nsavstt).eq.0 .or. itime.eq.(nptim+ntime))then
            call WRSTRT( io5, itime,time,dt, nbrigid,xpos,uvom,duvot,
     &                    fxym,ndgvelo,ndgpres,ndgelas,ndgphi,
     &                    ncpvelo,ncpelas, nbelast, ncpphi,
     &                    nvert,nnode,nelem,inod,x,y,nec,area,aspr,reft,
     &                    nbd,ibdnod,nic,ic,nbound,nside,pwall,bdseg,
     &                    u,p,dudt,umesh,selas,desdt,phi,psi,dphidt,
     &                    ngmx,ncmx,ncpg,ncps,ncpb,
     &                    noflow,ncpveln,ncpelan,
     &                    dist,meshgen)
          if(meshgen==1)then
            if(lwrmesh.or. itime.eq.(nptim+ntime))then
               call WRAMPHIMESH()
               lwrmesh=.false.
            endif
          endif
          write(filenum,'(i5.5)')itime
          filename='hh'//filenum//'.stt'
          call system('cp hh.stt '//filename)          
        endif
c
c       Print results on screen and into files
c
        open(unit=29,file='hh.tdt',access='append')
        if(itime.eq.1) rewind 29
        write(29,'(A,i4,A,1pe11.4,A,1pe11.4,A,l2)') 'itime=',itime,
     &                           ' time=',time,' dt=',dt,' lrmsh=',lrmsh
        close(unit=29)
c
        if ( nbrigid.gt.0 ) then
        nk = nbrigid
        if ( nbrigid.gt.2 ) nk=2
        write(*,'(A,2(3(1pe12.4),1x))')
     &       ' posit:',((xpos(i,k),i=1,3),k=1,nk)
        write(*,'(A,2(3(1pe12.4),1x))')
     &       ' veloc:',((uvom(i,k),i=1,3),k=1,nk)
        write(*,'(A,2(3(1pe12.4),1x))')
     &       ' force:',((fxym(i,k),i=1,3),k=1,nk)
c
        open(unit=26,file='hh.pos',access='append')
        open(unit=27,file='hh.vel',access='append')
        open(unit=28,file='hh.fcs',access='append')
        if(itime.eq.1) then
          rewind 26
          rewind 27
          rewind 28
        endif
        write(26,1015) time,((xpos(i,k),i=1,3),k=1,nbrigid)
        write(27,1015) time,((uvom(i,k),i=1,3),k=1,nbrigid)
        write(28,1015) time,((fxym(i,k),i=1,3),k=1,nbrigid)
1015    format(5001(1pe10.3,1x))
        close(unit=26)
        close(unit=27)
        close(unit=28)
        endif
c
      enddo
c
      if(SLVphi.and.noflow)then
      time=0.d0
      itime=0
             call WRSTRT( io5, itime,time,dt, nbrigid,xpos,uvom,duvot,
     &                    fxym,ndgvelo,ndgpres,ndgelas,ndgphi,
     &                    ncpvelo,ncpelas, nbelast, ncpphi,
     &                    nvert,nnode,nelem,inod,x,y,nec,area,aspr,reft,
     &                    nbd,ibdnod,nic,ic,nbound,nside,pwall,bdseg,
     &                    u,p,dudt,umesh,selas,desdt,phi,psi,dphidt,
     &                    ngmx,ncmx,ncpg,ncps,ncpb,
     &                    noflow,ncpveln,ncpelan)
      endif
c
C     Support for Petsc timing
c      call PLogEventEnd(PARTMOVER_2D_EVENT,
c     *                  PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL)
c      call PetscFinalize(ierr)
c
      deallocate(iwork,rwork)
      stop
      end
c
c
c***********************************************************************
      subroutine INITIA(ntime,nptim,time, iprint,nsavstt,norder,
     &                  maxvrt,maxnod,maxelt,maxnbd,maxbdp,maxptc,
     &                    maxnth,maxcoe,ngmx,ncmx,ncpg,ncps,ncpb,
     &                  nvert,nnode,nelem,inod,x,y,reft,nec,area,aspr,
     &                  nbd,ibdnod,nic,ic,nbound,nside,
     &                  u,p,dudt,umesh,selas,desdt,phi,psi,dphidt,
     &                  nbrigid,xpos,uvom,duvot,fxym,amas,fpart,ifixpt,
     &                  nbfluid,nbelast,nbtotal,icollision,epscollision,
     &                  nbshp,diaa,diab,ncr,ny,bdrefn,
     &                  flowtype,pwall,bdseg,
     &                  noexec,inert,ivheat,newton,upwind,lfrmt,
     &                  meshpr,bdpr,resold,rstart,trans,
     &                  icase,iax,ifluid,ivisc,itemp,itrel,icoe,
     &                  iforce,istres,ipnode,isnode,itmax,iorder,irstr,
     &                  ro,grav,pvis,ptemp,pelas,zconv,
     &                  pvelbd,ueta, nnew1,nnew2,
     &                  dt,dtmin,dtmax,dtrate,cflflow,cflphi,
     &                  disinc,velinc,rinc,
     &                  gamma,plambda,epsilon,shift,
     &                  ioptu,imu,lfilu,maxitu,epsu,lfilter,
     &                  iopts,lfils,maxits,epss,
     &                  ioptm,lfilm,maxitm,epsm,
     &                  iopti,lfili,maxiti,epsi,
     &                  ioptp,maxitp,epsp,
     &                  nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,
     &                    nrdstrm,nrdmsh,nrdphi,
     &                  ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,
     &                    ndgstrm,ndgmsh,ndgphi,
     &                  SLVvelo,SLVpres,SLVtemp,SLVelas,SLVphi,
     &                  ncpvelo,ncppres,ncptemp,ncpelas,ncpphi,
     &                  iwork,leniw,rwork,lenrw, iou,
     &                  noflow,ncpveln,ncpelan,dropone,
     &                  h1,h2,h3,dG,
     &                  wallrela,wallener,bdpre,
     &                  refr,refu,refp,
     &                  dist,domsize,meshgen)
c
c     This routine initilizes the velocity, pressure field,
c     particle positions, forces on particles and
c     generates a mesh / read information from a restart file
c
c      flowtype =  1: periodic flow in x-direction
c               =  2: flow in an infinite channel in x-direction
c               =  3: flow in a channel closed in x-direction
c               =  4: fluidization in x-direction
c               = 11: flow inside a circular cylinder
c               = 12: flow in an annulus between two cylinders
c               = 99: flow in a geometry specified by linear segments
c     flowtype  = 1xx: use GRUMMP for mesh generation
c               = 101: rectangular with solid walls
c               = 102: Simple shear flow (zero initial polymer stress)
c               = 103: planar extensional flow
c               = 104: drop retraction case(the same as 101 except for
c                      only 1/4 domain used. xd, yd are drop radii)
c               = 105: drop head-on collision(2D planar)
c               = 106: drop/interface impact(2D planar)
c               = 107: drop/interface impact(2D planar,
c                       upper bound open)
c               = 111: axisymmetric domain with solid walls
c                       (drop collision)
c               = 113: uniaxial extensional flow
c                       (zero initial polymer stress)
c               = 114: drop retraction case(xd, yd are drop radii,
c                       rd not used)
c               = 115: drop head-on collision(2D axisymmetric)
c               = 116: drop/interface impact(2D axisymmetric)
c               = 117: drop/interface impact(2D axisymmetric,
c                       upper bound open)
c     flowtype  = 12x:  moving contact line
c                 121:  couette device (iax=0)
c                 122:  plunging tape/capillary tube (iax=0/1)
c                 123:  drop in a channel (iax=0/1)
c                 124:  drop spreading (iax=0/1)
c                 129:  drop imbibition in a two staged tube. (iax=0/1)
c
c     OUTPUT:
c       ntime,nptime,time: time step information
c       xpos(3,nbrigid) : particle positions
c       u,p,dudt : velocity, pressure and acceleration field
c       nvert,nnode,nbd,nelem,inod,ibdnod,nbound,x,y,
c         nic,ic:  information about the mesh at t=0
c
c       pwall = information about each boundary segements
c            pwall(1,i) = 0xy -> line segment from (x1,y1)->(x2,y2)
c                   pwall(2,i) = x1 ; pwall(3,i) = y1
c                   pwall(4,i) = x2 ; pwall(5,i) = y2
c            pwall(1,i) = 1xy -> a circle at (x0,y0) with radius R
c                   pwall(2,i) = x0 ; pwall(3,i) = y0
c                   pwall(4,i) = R  ; pwall(5,i) = free
c              where xy is the boundary condition number
c                 xy = 00 -> periodic boundary segment
c                 xy = 10 -> x and y components of velocity are specified
c                            (inflow, generally non-zero)
c                 xy = 11 -> x and y components of velocity are specified
c                            (both are zero)
c                 xy = 12 -> x-component of velocity and
c                            y-component of stress(zero) are specified
c                 xy = 21 -> x-component of stress(zero) and
c                            y-component of velocity are specified
c                 xy = 22 -> x and y components of stress are specified
c                            (both are zero)
c
c     collision models
c       icollision = 0 no model
c       icollision = 1 modify the particle position only
c       icollision = 2 modify the particle position and velocity,
c                      and the body force on the particle!
c     Mar 2005, modified by Pengtao Yue
c        add the 1xx flowtypes
c
c***********************************************************************
      implicit none
      integer maxvrt,maxnod,maxelt,maxnbd,maxbdp,maxptc,maxnth,ntime,
     &        nptim,nbd,nnode,nvert,nelem,nbound,
     &        nic,ncr,ny,nbrigid,
     &        norder,iou,nnew1,nnew2, ngmx,ncmx,ncpg,ncps,ncpb,bdseg,
     &        ioptu,imu,lfilu,maxitu, iopts,lfils,maxits,
     &        ioptm,lfilm,maxitm, iopti,lfili,maxiti, iprint,nsavstt,
     &        ioptp,maxitp, maxcoe, nbfluid,nbelast,nbtotal,
     &        nrdgeom,nrdvelo,nrdpres,nrdtemp,nrdelas,nrdstrm,nrdmsh,
     &        ndggeom,ndgvelo,ndgpres,ndgtemp,ndgelas,ndgstrm,ndgmsh
      integer inod(ngmx,maxelt),ibdnod(maxnbd),nec(ncmx,maxelt),
     &        reft(maxelt),ic(maxbdp),nside(maxbdp),nbshp(maxptc)
      integer ncpvelo,ncppres,ncptemp,ncpelas
      logical SLVvelo,SLVpres,SLVtemp,SLVelas
      real*8 x(maxnod),y(maxnod),area(maxelt),aspr(maxelt),
     &       u(ncpg,maxnod),p(maxvrt),dudt(ncpg,maxnod),
     &       umesh(ncpg,maxnod),selas(ncps,maxvrt),desdt(ncps,maxvrt),
     &       xpos(ncpb,maxptc),uvom(ncpb,maxptc),duvot(ncpb,maxptc),
     &       fxym(ncpb,maxptc),amas(ncpb,maxptc),fpart(ncpb,maxptc),
     &       time,diaa(maxptc),diab(maxptc),
     &       disinc,velinc,rinc,dt,dtmax,dtmin,dtrate,cflflow,cflphi,
     &       epsu,epss,epsm,epsi,epsp,gapp
      integer leniw,lenrw,iwork(leniw)
      real*8 rwork(lenrw)
      logical ifixpt(ncpb),bdrefn,lfilter
      real*8 uwall1,uwall2,uin,upmax
      logical noexec,inert,ivheat,newton,upwind
      logical meshpr,bdpr,resold,rstart,trans
      integer icase,iax,ifluid,ivisc(maxcoe),itemp,itrel,icoe
      integer iforce,istres,ipnode,isnode,itmax,iorder
      double precision ro(maxcoe),grav(2,maxcoe),pvis(10,maxcoe),
     &                 ptemp(10,maxcoe),pelas(10,maxcoe),zconv
      real*8 pvelbd(10),pdmain(10)
      integer icollision,flowtype
      real*8 pwall(5,maxbdp-maxptc),epscollision
c     mesh size control for GRUMMP/Gmsh
      real(8) h1,h2,h3,dG,hh
c     Gmsh
      real(8) dist(maxvrt),domsize
      integer meshgen 
c     rescaling factors
      real(8) refr,refu,refp
c*CH start     start of Definitions related with Cahn-Hillard Equation
      real*8 phi(maxnod),psi(maxnod),dphidt(maxnod)
      integer ncpphi,nrdphi,ndgphi
      logical SLVphi
	real*8 gamma,epsilon,plambda,shift
c*CH end     end of Definitions related with Cahn-Hillard Equation
c     for initialization
      logical noflow
      integer ncpveln,ncpelan
c
      logical dropone
c
      real*8 xlength
      logical xperid
      common /xperid/ xlength,xperid
c
      integer irstr
c     moving contact line, (1/Gamma) and wall energy density
      real(8) wallrela(maxbdp),wallener(maxbdp)
      real(8) bdpre(maxbdp)
c     local variables
      character*1 title
      logical random
      integer i,k,j,k1, ndglb, idx, itime
      real*8 dens,gx,gy,dpdx,pi,pi2,vlm
      real*8 ybar,ybar1,dudy,ueta(4,2*maxnth+1),uthin,duth,etath
      real*8 lambd1,tmp(50),xmax,xmin,width
      logical lfrmt
c     local variables related with C-H equation
      real*8 rd,xd,yd  ! drop size and position
c
      data pi/3.14159265358979d0/, pi2/0.0174532925199433d0/
c
      open(unit=iou,file='partdata')
        read(iou,1001) title
        read(iou,1001) title
        read(iou,*) rstart,random,(ifixpt(i),i=1,3),noexec,lfrmt,noflow
        read(iou,1001) title
        read(iou,*) nbtotal,ncr,ny,h1,h2,h3,dG
        read(iou,1001) title
        read(iou,*) flowtype,(pdmain(i),i=1,3)
c     modifications related with C-H equation
        read(iou,1001) title
        read(iou,*) SLVphi, gamma, plambda, epsilon,shift
        read(iou,1001) title
        read(iou,*) rd,xd,yd,dropone
c
        read(iou,1001) title
        read(iou,*) dtmin,dtmax,dtrate,cflflow,cflphi,disinc,velinc,rinc
        read(iou,1001) title
        read(iou,*) ntime,iprint,nsavstt,nnew1,nnew2,irstr,norder,iorder
        read(iou,1001) title
        read(iou,*) (pvelbd(i),i=1,4)
        read(iou,1001) title
        read(iou,*) dens,dpdx,gx,gy
        read(iou,1001) title
        read(iou,*) refr,refu,refp
        read(iou,1001) title
        read(iou,*) icase,iax,icoe,ifluid,itemp,itrel
        read(iou,1001) title
        do k=1,icoe
           read(iou,*) ro(k),ivisc(k),(pvis(i,k),i=1,9)
        enddo
        read(iou,1001) title
        if(SLVphi.and.(icoe.eq.2))then
           j=1
        else
           j=icoe
        endif
        do k=1,j
           read(iou,*) (pelas(i,k),i=1,10)
        enddo
        read(iou,1001) title
        do k=1,j
           read(iou,*) (ptemp(i,k),i=1,10)
        enddo
        read(iou,1001) title
        read(iou,*) meshpr,bdpr,inert,ivheat,newton,upwind,resold,trans
        read(iou,1001) title
        read(iou,*) iforce,istres,ipnode,isnode,itmax,zconv
        read(iou,1001) title
        read(iou,*) ioptu,imu,lfilu,maxitu,epsu,lfilter
c     modifications related with C-H equation
c        read(iou,1001) title
c        read(iou,*) ioptphi,imphi,lfilphi,maxitphi,epsphi,lfilterphi
c
        read(iou,1001) title
        read(iou,*) ioptm,lfilm,maxitm,epsm
        read(iou,1001) title
        read(iou,*) iopts,lfils,maxits,epss
        read(iou,1001) title
        read(iou,*) iopti,lfili,maxiti,epsi
        read(iou,1001) title
        read(iou,*) ioptp,maxitp,epsp
        read(iou,1001) title
        read(iou,*) icollision,epscollision,gapp
        if ( flowtype.eq.99 ) then
           read(iou,1001) title
           read(iou,*) bdseg
           do k=1,bdseg
             read(iou,*) (pwall(i,k),i=1,5)
           enddo
        endif
        if ( .not.random ) then
          read(iou,1001) title
          do k=1,nbtotal
             read(iou,*) j,nbshp(k),diaa(k),diab(k),(xpos(i,k),i=1,3)
          enddo
        endif
c
      close(iou)
1001  format(a1)
c
c     check for fluid type (icase)
c     ----------------------------
      if ( icase.eq.1 ) then
        write(*,'(A,i4,A)')
     &    ' Motion of ',nbtotal,' particle(s) in a Newtonian fluid'
        do k=1,icoe
           pelas(6,k)= 1.d0
        enddo
      else if ( icase.eq.4 ) then
        write(*,'(A,i4,A)')
     &    ' Motion of ',nbtotal,' particle(s) in a viscoelastic fluid'
      else if ( icase.eq.2 .or. icase.eq.3 .or. icase.eq.5 ) then
        write(*,'(A,i1,A)')
     &    ' Warning: Fluid type (icase=',icase,') may not work fine!!!'
      else
        write(*,*) ' Error: invalid fluid type: use different icase!'
        stop
      endif
c
c     check for the solver options
c     ----------------------------
      k = mod(ioptu,10)
      if ( k.eq.1 .or. k.eq.2 ) lfilu=0
      k = mod(ioptu/10,10)
      if ( k.eq.0 .or. k.eq.3 ) imu=0
c
      k = mod(ioptm,10)
      if ( k.eq.1 .or. k.eq.2 ) lfilm=0
c
      k = mod(iopts,10)
      if ( k.eq.1 .or. k.eq.2 ) lfils=0
c
      k = mod(iopti,10)
      if ( k.eq.1 .or. k.eq.2 ) lfili=0
c
c     check the flow type
c     -------------------
      if ( flowtype.ne.1 .and. random )  then
         write(*,*) 'You can not use random particle distribution',
     &              ' for the flow type selected'
         stop
      endif
      if ( flowtype/100.eq.0)then
      if ( flowtype.eq.1 ) then
        write(*,*) 'periodic boundary in x-direction'
      elseif ( flowtype.eq.2 ) then
        write(*,*) 'infinite channel in x-direction'
      elseif ( flowtype.eq.3 ) then
        write(*,*) 'channel is closed in x-direction'
      elseif ( flowtype.eq.4 ) then
        write(*,*) 'fluidization in x-direction'
      elseif ( flowtype.eq.11 ) then
        write(*,*) 'flow in a circular cylinder'
      elseif ( flowtype.eq.12 ) then
        write(*,*) 'flow in an annulus between two circular cylinders'
      elseif ( flowtype.eq.99 ) then
        write(*,*) 'flow in a geometry specified by linear segments'
      else
         stop 'error in flowtype2 value'
      endif
      elseif( flowtype/100.eq.1)then
        if(meshgen==1)then
          write(*,*)'***GRUMMP is used for mesh generation***'
        elseif(meshgen==2)then
          write(*,*)'***Gmsh is used for mesh generation***'
        else
          stop 'error in specifying mesh generator'
        endif
      if(mod(flowtype,100).eq.1)then
         write(*,*)'rectangular with solid walls'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for planar flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.2)then
         write(*,*)'Couette flow'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for planar flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.3)then
         write(*,*)'Planar extensional flow'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for planar flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.4)then
         write(*,*)'Drop retraction(planar)'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for planar flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.5)then
         write(*,*)'Drop head-on collision(2D planar)'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for planar flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.6)then
         write(*,*)'Drop/interface impact(2D planar)'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for planar flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.7)then
         write(*,*)'Drop/interface impact(2D planar, upper boundary open
     &)'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for planar flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.11)then
         write(*,*)'axisymmetric flow'
         if(iax.ne.1)then
            write(*,*)'iax is enforced to be 1 for axisymmetric flow'
            iax=1
         endif
      elseif(mod(flowtype,100).eq.13)then
         write(*,*)'uniaxial extensional flow'
         if(iax.ne.1)then
            write(*,*)'iax is enforced to be 1 for axisymmetric flow'
            iax=1
         endif
      elseif(mod(flowtype,100).eq.14)then
         write(*,*)'Drop retraction(axisymmetric)'
         if(iax.ne.1)then
            write(*,*)'iax is enforced to be 1 for axisymmetric flow'
            iax=1
         endif
      elseif(mod(flowtype,100).eq.15)then
         write(*,*)'Drop head-on collision(2D axisymmetric)'
         if(iax.ne.1)then
            write(*,*)'iax is enforced to be 1 for axisymmetric flow'
            iax=1
         endif
      elseif(mod(flowtype,100).eq.16)then
         write(*,*)'Drop/interface impact(2D axisymmetric)'
         if(iax.ne.1)then
            write(*,*)'iax is enforced to be 1 for axisymmetric flow'
            iax=1
         endif
      elseif(mod(flowtype,100).eq.17)then
         write(*,*)'Drop/interface impact(2D axisymmetric, upper boundar
     &y open)'
         if(iax.ne.1)then
            write(*,*)'iax is enforced to be 1 for axisymmetric flow'
            iax=1
         endif
      elseif(mod(flowtype,100).eq.21)then
         write(*,*)'Moving contact line: couette device'
         if(iax.ne.0)then
            write(*,*)'iax is enforced to be 0 for axisymmetric flow'
            iax=0
         endif
      elseif(mod(flowtype,100).eq.22)then
         write(*,*)'Moving contact line: plunging tape/capillary tube'
      elseif(mod(flowtype,100).eq.23)then
         write(*,*)'Moving contact line: drop in a channel'
      elseif(mod(flowtype,100)==24)then
         write(*,*)'***Drop speading problem***'
         if(iax==0)then
            write(*,*)'**2D planar**'
         else
            write(*,*)'**2D axisymmetric**'
         endif
      elseif(mod(flowtype,100)==29)then
         write(*,*)'***Drop migration in a two-staged tube***'
         if(iax==0)then
            write(*,*)'**2D planar**'
         else
            write(*,*)'**2D axisymmetric**'
         endif
      else
         stop 'error in flowtype2 value'
      endif
      else
         stop 'error in flowtype1 value'
      endif
c
c     build pwall for boundary information
c     ------------------------------------
      xperid = .false.
      if(flowtype/100.eq.0)then
      if ( flowtype.eq.1 ) then
         bdseg = 4
c        note: first list the periodic boundary
         pwall(1,1) = 0
         pwall(2,1) = pdmain(2)
         pwall(3,1) = pdmain(1)
         pwall(4,1) = pdmain(2)
         pwall(5,1) = 0.d0
         pwall(1,2) = 11
         pwall(2,2) = pdmain(2)
         pwall(3,2) = 0.d0
         pwall(4,2) = pdmain(3)
         pwall(5,2) = 0.d0
         pwall(1,3) = 0
         pwall(2,3) = pdmain(3)
         pwall(3,3) = 0.d0
         pwall(4,3) = pdmain(3)
         pwall(5,3) = pdmain(1)
         pwall(1,4) = 11
         pwall(2,4) = pdmain(3)
         pwall(3,4) = pdmain(1)
         pwall(4,4) = pdmain(2)
         pwall(5,4) = pdmain(1)
         xperid = .true.
         xlength = pdmain(3)-pdmain(2)
      elseif ( flowtype.eq.2 .or. flowtype.eq.3 .or.
     &         flowtype.eq.4 ) then
         bdseg = 4
         pwall(1,1) = 11
         pwall(2,1) = pdmain(2)
         pwall(3,1) = pdmain(1)
         pwall(4,1) = pdmain(2)
         pwall(5,1) = 0.d0
         pwall(1,2) = 11
         pwall(2,2) = pdmain(2)
         pwall(3,2) = 0.d0
         pwall(4,2) = pdmain(3)
         pwall(5,2) = 0.d0
         pwall(1,3) = 11
         pwall(2,3) = pdmain(3)
         pwall(3,3) = 0.d0
         pwall(4,3) = pdmain(3)
         pwall(5,3) = pdmain(1)
         pwall(1,4) = 11
         pwall(2,4) = pdmain(3)
         pwall(3,4) = pdmain(1)
         pwall(4,4) = pdmain(2)
         pwall(5,4) = pdmain(1)
      elseif ( flowtype.eq.11 ) then
         bdseg = 1
         pwall(1,1) = 111
         pwall(2,1) = 0.d0
         pwall(3,1) = 0.d0
         pwall(4,1) = -0.5*pdmain(1)
      elseif ( flowtype.eq.12 ) then
         bdseg = 2
         pwall(1,1) = 111
         pwall(2,1) = 0.d0
         pwall(3,1) = 0.d0
         pwall(4,1) = -0.5*pdmain(1)
         pwall(1,2) = 111
         pwall(2,2) = pdmain(3)
         pwall(3,2) = 0.d0
         pwall(4,2) = 0.5*pdmain(2)
      elseif ( flowtype.eq.99 ) then
         do i=1,bdseg
            if ( abs(pwall(1,i)).lt.0.01) then
               xperid = .true.
               if ( i.ne.1 ) then
c              reordering pwall
                 do j=1,5
                    do k=1,bdseg
                       tmp(k) = pwall(j,k)
                    enddo
                    do k=1,bdseg
                       k1 = i+k-1
                       if ( k1.gt.bdseg ) k1=k1-bdseg
                       pwall(j,k) = tmp(k1)
                    enddo
                 enddo
               endif
               go to 10
            endif
         enddo
 10      if ( xperid ) then
            xmax = pwall(2,1)
            xmin = pwall(2,1)
            do i=2,bdseg
               xmax = max(xmax,pwall(2,i))
               xmin = min(xmin,pwall(2,i))
            enddo
            xlength = xmax-xmin
         endif
      endif
      elseif(flowtype/100.eq.1)then
c       if(flowtype.eq.101.or.flowtype.eq.102.or.flowtype.eq.111)then
         if(flowtype==129)then
            call flow129bdry(iou,pwall,xd,yd,rd)
            domsize=max(pwall(2,1),pwall(3,1)+pwall(4,1))
         else
            xmin=pdmain(2)
            xmax=pdmain(3)
            width=pdmain(1)
            pwall(1,1) = xmin
            pwall(2,1) = xmax
            pwall(3,1) = width
            domsize=max(xmax-xmin,width)
            if(dropone)then
               hh=h3
            else
               hh=h2
            endif
            open(iou,file='amphi.bdry')
            write(iou,*)8,8
            write(iou,1010)xmin,0.d0
            write(iou,1010)xmax,0.d0
            write(iou,1010)xmax,width
            write(iou,1010)xmin,width
            write(iou,1010)xmin+hh,0.d0
            write(iou,1010)xmax,hh
            write(iou,1010)xmax-hh,width
            write(iou,1010)xmin,width-hh
            write(iou,1020)'r',1,'b',1,3,0,4,1
            write(iou,1020)'r',1,'b',2,3,1,5,2
            write(iou,1020)'r',1,'b',3,3,2,6,3
            write(iou,1020)'r',1,'b',4,3,3,7,0
            close(iou)
         endif
1010  format(2(1x,g12.5))
1020  format(1x,'polyline',2(1x,a4,1x,i4),5(1x,i4))

c       elseif(flowtype.eq.113)then
c         xmin=pdmain(2)   !not used
c         xmax=pdmain(3)   !not used
c         width=pdmain(1)  !R
c         pwall(1,1) = width
c         pwall(2,1) = xmin
c         pwall(3,1) = xmax
c         open(iou,file='amphi.bdry')
c         write(iou,1010)4,4
c         write(iou,1010)0.0,0.0
c         write(iou,1010)width,0.d0
c         write(iou,1010)sqrt(2./3.)*width,sqrt(1./3.)*width
c         write(iou,1010)0.0,width
c         write(iou,1020)'r',1,'b',1,2,0,1
c         write(iou,1030)'r',1,'b',2,width,1,2
c         write(iou,1030)'r',1,'b',3,width,2,3
c         write(iou,1020)'r',1,'b',4,2,3,0
c         close(iou)
c1030  format(1x,'arc',2(1x,a4,1x,i4),1x,g12.5, 2(1x,i4))
c         endif
      endif
c
c     check the boundary condition
c     ----------------------------
      uwall1 = pvelbd(1)
      uwall2 = pvelbd(2)
      uin    = pvelbd(3)
      upmax  = pvelbd(4)
      if ( upmax.eq.0 .and. uin.eq.0 .and. uwall1.eq.0 .and.
     &     uwall2.eq.0 )  then
        if ( flowtype.ne.1 .or. dpdx.eq.0 ) then
          write(*,*) 'Sedimentation problem'
        else
          write(*,'(A,e12.5)') ' Poiseuille flow with dp/dx=',dpdx
        endif
      endif
      if ( upmax.ne.0 )
     &   write(*,'(A,e12.5)') ' Poiseuille inflow with Umax=',upmax
      if ( uin.ne.0 )
     &   write(*,'(A,e12.5)') ' uniform inflow with Uin=',uin
      if ( flowtype.lt.10 .and. (uwall1.ne.0.or.uwall2.ne.0) )
     &   write(*,'(A,e12.5,A,e12.5)')
     &   ' planar Couette flow with U1=',uwall1,' U2=',uwall2
      if ( flowtype.eq.11 .and. uwall1.ne.0 ) then
         write(*,'(A,e12.5)') ' angular velocity (rpm)=',pvelbd(1)
         pvelbd(1) = pvelbd(1)*pi/30.d0
      endif
      if ( flowtype.eq.12 .and. (uwall1.ne.0.or.uwall2.ne.0) ) then
         write(*,'(A,e12.5)') ' angular velocity (outer rpm)=',pvelbd(1)
         write(*,'(A,e12.5)') ' angular velocity (inner rpm)=',pvelbd(2)
         pvelbd(1) = pvelbd(1)*pi/30.d0
         pvelbd(2) = pvelbd(2)*pi/30.d0
      endif
c
      if ( icollision.gt.0 )
     &   write(*,'(A,i2)') ' Using collision model No.',icollision
c
c     generate initial particle positions
c     -----------------------------------
      if(flowtype/100.eq.1.and.nbtotal.ne.0)then
         write(*,*)"GRUMMP can't deal with particle. nbtotal enforced to
     & be zero!"
         nbtotal=0
      endif
c     bdseg is of no use for grummp, just to be compatible with old part
      bdseg=4
c
      if ( random ) then
         do i=1,nbtotal
           nbshp(i) = 1
           diaa(i) = 1.d0
           diab(i) = 1.d0
         enddo
         if (.not.rstart ) then
            call RANDIST(nbtotal,xlength,pdmain(1),xpos,diaa,gapp)
c            do k=1,nbtotal
c              write(42,'(2i4,5f10.5)') k,1,1.d0,1.d0,(xpos(j,k),j=1,3)
c            enddo
        endif
      endif
c
      do k=1,nbtotal
        xpos(3,k) = xpos(3,k)*pi2
      enddo
      disinc = disinc*diaa(1)
      rinc   = rinc*pi2
c
c     reordering the particles
c     ------------------------
      do k=1,nbtotal
         iwork(k)   = nbshp(k)
         fpart(1,k) = diaa(k)
         fpart(2,k) = diab(k)
         amas(1,k) = xpos(1,k)
         amas(2,k) = xpos(2,k)
         amas(3,k) = xpos(3,k)
      enddo
      nbrigid = 0
      do k=1,nbtotal
         if ( iwork(k)/100.eq.0 ) then
             nbrigid = nbrigid+1
             nbshp(nbrigid)  = iwork(k)
             diaa(nbrigid)   = fpart(1,k)
             diab(nbrigid)   = fpart(2,k)
             xpos(1,nbrigid) = amas(1,k)
             xpos(2,nbrigid) = amas(2,k)
             xpos(3,nbrigid) = amas(3,k)
         endif
      enddo
      nbfluid = 0
      do k=1,nbtotal
         if ( iwork(k)/100.eq.1 ) then
             nbfluid = nbfluid+1
             nbshp(nbrigid+nbfluid)  = iwork(k)
             diaa(nbrigid+nbfluid)   = fpart(1,k)
             diab(nbrigid+nbfluid)   = fpart(2,k)
             xpos(1,nbrigid+nbfluid) = amas(1,k)
             xpos(2,nbrigid+nbfluid) = amas(2,k)
             xpos(3,nbrigid+nbfluid) = amas(3,k)
         endif
      enddo
      nbelast = 0
      do k=1,nbtotal
         if ( iwork(k)/100.eq.2 ) then
             nbelast = nbelast+1
             nbshp(nbrigid+nbfluid+nbelast)  = iwork(k)
             diaa(nbrigid+nbfluid+nbelast)   = fpart(1,k)
             diab(nbrigid+nbfluid+nbelast)   = fpart(2,k)
             xpos(1,nbrigid+nbfluid+nbelast) = amas(1,k)
             xpos(2,nbrigid+nbfluid+nbelast) = amas(2,k)
             xpos(3,nbrigid+nbfluid+nbelast) = amas(3,k)
         endif
      enddo
c
c     volume, mass, weight, moment of rigid particles
c     (static pressure in the continuous phase and the resultant
c      buoyancy force are taken out of the equations)
c     ----------------------------------------------------------
      do k=1,nbrigid
         if ( nbshp(k).eq.1 ) then
            vlm = 0.25*pi*diaa(k)*diab(k)
            amas(1,k) = dens*vlm
            amas(2,k) = amas(1,k)
            amas(3,k) = amas(1,k)*(diaa(k)**2+diab(k)**2)/16.
         elseif( nbshp(k).eq.2 ) then
            vlm = diaa(k)*diab(k)
            amas(1,k) = dens*vlm
            amas(2,k) = amas(1,k)
            amas(3,k) = amas(1,k)*(diaa(k)**2+diab(k)**2)/12.
         endif
         fpart(1,k) = (dens-ro(1))*vlm*gx - vlm*ro(1)*dpdx
         fpart(2,k) = (dens-ro(1))*vlm*gy
         fpart(3,k) = 0.d0
      enddo
c
c     body force terms in the fluid momentum equations
c     ------------------------------------------------
	
	if(SLVphi)then
c     dpdx maynot work properly in two phase systems
c     also this dpdx is not the real one, it's dpdx/ro(1)
      if(dpdx.ne.0.d0)then
         write(*,*)'dpdx may not work properly in two phase system'
      endif
      do i=1,icoe
         grav(1,i) = gx -dpdx
         grav(2,i) = gy
      enddo
      else
      grav(1,1) = -dpdx
      grav(2,1) = 0.d0
      do i=2,icoe
         grav(1,i) = (ro(i)-ro(1))/ro(i)*gx - dpdx
         grav(2,i) = (ro(i)-ro(1))/ro(i)*gy
      enddo
      endif
c
c     calculation of SLVxxxx, ncpxxx
c     ------------------------------
      SLVvelo = .false.
      SLVpres = .false.
      SLVtemp = .false.
      SLVelas = .false.
      ncpvelo = 0
      ncppres = 0
      ncptemp = 0
      ncpelas = 0
      if ( icase.eq.1 ) then
         SLVvelo = .true.
         SLVpres = .true.
         ncpvelo = 2
         ncppres = 1
         if ( iax.eq.2 ) ncpvelo = ncpvelo + 1
      elseif ( icase.eq.2 ) then
         SLVvelo = .true.
         SLVpres = .true.
         SLVtemp = .true.
         ncpvelo = 2
         ncppres = 1
         ncptemp = 1
         if ( iax.eq.2 ) ncpvelo = ncpvelo + 1
      elseif ( icase.eq.3 ) then
         SLVtemp = .true.
         ncptemp = 1
      elseif ( icase.eq.4 ) then
         SLVvelo = .true.
         SLVpres = .true.
         SLVelas = .true.
         ncpvelo = 2
         ncppres = 1
         ncpelas = 3 + iax
      elseif ( icase.eq.5 ) then
         SLVelas = .true.
         ncpvelo = 2
         ncpelas = 3 + iax
      endif
c*CH start
      if(SLVphi)then
         ncpphi=1
      else
         ncpphi=0
      endif
c*CH end
c
c     velocity profile for Poiseuille flow of a shear-dependent fluid
c     ---------------------------------------------------------------
      if ( (flowtype.eq.2 .or. flowtype.eq.4) .and.
     &      ivisc(1).ne.1 .and. upmax.ne.0.d0 )
     &       call pthin(maxnth,ivisc(1),pvis,upmax,pdmain(1),ueta)
c
      if ( rstart ) then
c -----------------------------------------------------
c       Reading in restart information from 'hh.stt'
c -----------------------------------------------------
c        dtmin = 0.d0
        open(unit=iou,file='hh.stt',form='unformatted')
          read(iou) nptim,time,dt
          read(iou) ((xpos(i,k),i=1,ncpb),k=1,nbrigid)
          read(iou) ((uvom(i,k),i=1,ncpb),k=1,nbrigid)
          read(iou) ((fxym(i,k),i=1,ncpb),k=1,nbrigid)
          read(iou) ((duvot(i,k),i=1,ncpb),k=1,nbrigid)
          read(iou) nvert,nnode,nelem
          read(iou) ((inod(k,i),k=1,ngmx),i=1,nelem)
          read(iou) ((nec(k,i),k=1,ncmx),i=1,nelem)
          read(iou) (x(i),y(i),i=1,nnode)
          read(iou) (area(i),aspr(i),i=1,nelem)
          read(iou) (reft(i),i=1,nelem)
          read(iou) nbd,nic,nbound
          read(iou) (ibdnod(i),i=1,nbd)
          read(iou) (ic(i),i=1,nic+1)
          read(iou) (nside(i),i=1,nbound)
          read(iou) ndgvelo
          read(iou) ((u(j,i),i=1,ndgvelo),j=1,ncpvelo)
          read(iou) ((dudt(j,i),i=1,ndgvelo),j=1,ncpvelo)
          read(iou) ((umesh(j,i),i=1,ndgvelo),j=1,ncpvelo)
          read(iou) ndgpres
          read(iou) (p(i),i=1,ndgpres)
          read(iou) ndgelas
          do j=1,ncpelas
            read(iou) (selas(j,i),i=1,ndgelas)
            read(iou) (desdt(j,i),i=1,ndgelas)
          enddo
          if ( nbelast.gt.0 ) then
            do j=1,3
              read(iou) (selas(j,i),i=1,ndgelas)
            enddo
          endif
c*CH start     modifications related with C-H equation
          if(ncpphi.eq.1)then
          read(iou) ndgphi
          read(iou) phi(1:ndgphi)
          read(iou) psi(1:ndgphi)
          read(iou) dphidt(1:ndgphi)
          endif
c*CH end
          read(iou) bdseg
          do j=1,bdseg
            read(iou) (pwall(i,j),i=1,5)
          enddo
          if(meshgen==2) read(iou)dist(1:nvert)
        close(unit=iou)
      else
c -----------------------------------------------------
c       generate information from fresh
c -----------------------------------------------------
        dt = dtmin
        nptim = 0
        time = 0.d0
        itime = 0
c
c       write(*,*) '  Generating new mesh'
        if(flowtype/100.eq.0)then
        call RMESH (itime,flowtype,pwall,bdseg,bdrefn,ncr,ny,
     &              nbrigid,nbfluid,nbelast,nbtotal,
     &              nbshp,diaa,diab,xpos,
     &              nvert,nnode,nelem,inod,x,y,reft,nec,area,aspr,
     &              nbd,ibdnod,nic,ic,nbound,nside,
     &              maxvrt,maxnod,maxelt,maxnbd,maxbdp,
     &              iwork,leniw,rwork,lenrw )
        else
          if(meshgen==1)then
            call GRinitmsh(flowtype,maxvrt,maxnod,maxelt,maxnbd,maxbdp,
     &                     ngmx,ncmx,nvert,nnode,x,y,
     &                     nelem,inod,nec,area,reft,aspr,
     &                     nbd,ibdnod,nic,ic,nbound,nside,
     &                     leniw,lenrw,iwork,rwork,
     &                     rd,xd,yd,dropone,epsilon,h1,h2,h3,dG)
          elseif(meshgen==2)then
            call gmshinitmsh(maxvrt,maxnod,maxelt,maxnbd,maxbdp,
     &                       ngmx,ncmx,nvert,nnode,nelem,
     &                       inod,nec,x,y,dist,
     &                       nbd,ibdnod,nic,ic,nbound,nside,
     &                       area,reft,aspr,
     &                       flowtype,rd,xd,yd,dropone,epsilon,
     &                       h1,h2,h3,dG,domsize,
     &                       leniw,lenrw,iwork,rwork)
          endif
        endif
c
c       initializing particle variables
c       -------------------------------
        do j=1,ncpb
          do i=1,nbrigid
             uvom(j,i) = 0.d0
             fxym(j,i) = 0.d0
             duvot(j,i) = 0.d0
          enddo
        enddo
c
        ndgvelo = ndglb(nrdvelo)
        ndgpres = ndglb(nrdpres)
        ndgelas = ndglb(nrdelas)
c
c       initializing fluid velocities
c       -----------------------------
        do j=1,ncpg
          do i=1,ndgvelo
             u(j,i) = 0.d0
             dudt(j,i) = 0.d0
             umesh(j,i) = 0.d0
          enddo
        enddo
        if ( flowtype.eq.2 .or. flowtype.eq.4 ) then
          do i=1,ndgvelo
             if ( ivisc(1).ne.1 .and. upmax.ne.0.d0 ) then
                call thin(y(i),ueta,maxnth,u(1,i),duth,etath)
             else
                ybar = y(i)/pdmain(1)
                ybar1 = 1.d0-ybar
                idx = 1
                if ( ybar.lt.1.d-10 .or. ybar1.lt.1.d-10 ) idx=0
                u(1,i) = uin*idx + uwall1*ybar1 + uwall2*ybar
     &                 + 4.0*upmax*ybar*ybar1
             endif
          enddo
        endif
        if(flowtype.eq.102)then
c
c     special treatment for couette flow
c     ----------------------------------
         do i=1,ndgvelo
            u(1,i)=(y(i)/width)*(uwall2-uwall1)+uwall1
         enddo
        endif
        if(flowtype.eq.103)then
c
c     special treatment for planar extensional flow
c     ---------------------------------------------
            do i=1,ndgvelo
               u(1,i)=-uwall1*x(i)
               u(2,i)=uwall1*y(i)
            enddo
        endif
        if(flowtype.eq.113)then
c
c     special treatment for uniaxial extensional flow
c     -----------------------------------------------
            do i=1,ndgvelo
               u(1,i)=-0.5d0*uwall1*x(i)
               u(2,i)=uwall1*y(i)
            enddo
         endif
c
c       initializing fluid pressure
c       ---------------------------
        do i=1,ndgpres
           p(i) = 0.d0
        enddo
c
c     Initial Pressure for drop/interface impact (flowtype 106 & 116)
c     ---------------------------------------------------------------
         if(flowtype.eq.106.or.flowtype.eq.116.or.flowtype.eq.107.or.
     &      flowtype.eq.117)then
            if(dropone)then
            do i=1,ndgpres
               if(y(i).gt.xd)then
                  p(i)=0.d0
               else
                  p(i)=(ro(1)-ro(2))*gy*(y(i)-xd)
               endif
            enddo
            else
            do i=1,ndgpres
               if(y(i).gt.xd)then
                  p(i)=0.d0
               else
                  p(i)=(ro(2)-ro(1))*gy*(y(i)-xd)
               endif
            enddo
            endif
         endif

c
c       initializing fluid elastic stress
c       ---------------------------------
        do j=1,ncpelas
          do i=1,ndgelas
             selas(j,i) = 0.d0
             desdt(j,i) = 0.d0
          enddo
        enddo
        if ( (flowtype.eq.2 .or. flowtype.eq.4).and.ncpelas.gt.0 ) then
          lambd1= pelas(1,1)*(1.0-pelas(6,1))
          do i=1,ndgelas
            if ( ivisc(1).ne.1 .and. upmax.ne.0.d0 ) then
               call thin(y(i),ueta,maxnth,uthin,dudy,etath)
            else
               etath = pvis(1,1)
               ybar = y(i)/pdmain(1)
               dudy = (uwall2-uwall1+4.*upmax*(1-2*ybar))/pdmain(1)
            endif
            selas(1,i) = 2.0*lambd1*etath*dudy**2
            selas(3,i) = etath*dudy
          enddo
        endif
c
        if ( nbelast.gt.0 ) then
           do j=1,3
              do i=1,ndgelas
                 selas(j,i) = 0.d0
              enddo
           enddo
        endif
c
c*CH start
c
c     Initializing phase function
c     ---------------------------
        ndgphi  = ndglb(nrdphi)
        if(ncpphi.eq.1)then
          call initphi(ndgphi,x,y,phi,dropone,epsilon,rd,xd,yd,flowtype)
          psi(1:ndgphi)=0.
          if(meshgen==2)
     &       call epc(nvert, nnode, nelem, ngmx, 10.d0*h1,
     &                inod, phi, x, y, dist,
     &                iwork,rwork,leniw,lenrw)
        endif
c* CH end
      endif

c*CH start
      if(ncpphi.eq.0)then
        phi(1:ndgphi)=1.
        psi(1:ndgphi)=0.
      endif
c*CH end
      ndggeom = ndglb(nrdgeom)
      ndgtemp = ndglb(nrdtemp)
      ndgstrm = ndglb(nrdstrm)
      ndgmsh  = ndglb(nrdmsh)
c
c     velocity initialization and parameters input
c        for moving contact line problem
c
      wallrela(1:nic)=0.d0
      wallener(1:nic)=0.d0
      bdpre(1:nic)=0.d0
      if(flowtype/10==12)then
         if(flowtype==129)then
            open(iou,file='flw129data')
            read(iou,1001) title
            read(iou,1001) title
            read(iou,*)bdpre(1),bdpre(5)
            read(iou,1001) title
            read(iou,1001) title
            read(iou,1001) title
            read(iou,1001) title
         else
            open(iou,file='contdata')
            read(iou,1001) title
            read(iou,1001) title
            read(iou,*)pvelbd(1)
         endif
         read(iou,1001) title
         do i=1,nic
         read(iou,*)wallrela(i),wallener(i)
         enddo
         close(iou)

c           initializing velocity field
            if(.not.rstart)then
               if(flowtype==121)then
                  uwall1=pvelbd(1) !vel on left wall
                  uwall2=-uwall1
                  do i=1,ndgvelo
                     u(2,i)=uwall1+(uwall2-uwall1)*(x(i)/(xmax-xmin))
                  enddo
               elseif(flowtype==122)then
                  uwall1=-pvelbd(1)  !vel on right wall
                  if(iax==0)then
                     uwall2=-2.d0/3.d0*uwall1   !vel at center line
                  else
                     uwall2=-0.5d0*uwall1
                  endif
                  do i=1,ndgvelo
                     u(2,i)=uwall2+uwall1*(1.d0-(x(i)/(xmax-xmin))**2)
                  enddo
               else
                  pvelbd(1)=0.d0
               endif
            endif
      endif
c
      if(noflow)then
c     set velocity to zero and calculate the initial phase function.
         write(*,*)'calculating the initial phase function with no flow'
         ncpveln=ncpvelo
         ncpelan=ncpelas
         SLVvelo=.false.
         SLVpres=.false.
         SLVelas=.false.
         SLVtemp=.false.
         ncpvelo=0
         ncppres=0
         ncpelas=0
         ncptemp=0
      endif
      return
      end
c
c
c***********************************************************************
      subroutine BOUNDC (flowtype,iax,pwall,bdseg,nnode,x,y, nbd,ibdnod,
     &                   nic,ic,nbrigid,nbfluid,nbelast,xpos,uvom,
     &                   ncpb,ncpu,ncps, icase,pvis,pelas, pvelbd,ueta,
     &                   ivisc,maxnth, ibdu,ibdt,ibds,densu,denst,denss,
     &                   ibdn,densn,bdpre,
     &                   ndgphi,phi,time,ipnode,ibdp,densp)
c
c     This routine applies the boundary conditions and generates the
c     boundary data for flow solver
c
c     Input :
c        flowtype: index for flow types
c        pwall(5,*): boundary information
c        x(nnode), y(nnode) = coord. of the nodes
c        nbrigid,xpos(3,nbrigid) = positions of the particles
c        nbfluid  = number of deformable fluid particles
c        nbelast  = number of deformable elastic particles
c        nbd    = number of node on the boundaries
c        ibdnod(nbd) = boundary node description
c        nic = number of boundary sections
c        ic(nic+1) = index for boundary sections
c        uvom = velocities of the particles
c        pvelbd(10) = inflow conditions
c        time: time, used for the time dependent BC
c     Output :
c        ibdu,densu: velocity BC
c        ibds,denss: extra stress BC
c        ibdt,denst: temperature
c        ipnode: the boundary node index at which pressure is fixed
c                for flowtype 103&113, 104&114
c***********************************************************************
c     For flowtype 102,103,113, only valid for Oldroyd-B/Newontian,
c        or N/N systems
c     Modified by Pengtao Yue, Jul 28, 2005
c        ibdp and densp added, those two will be used whenever ipnode=-1
c***********************************************************************
      implicit none
      integer iax,nnode,nbd,ibdnod(nbd),nic,ic(nic+1), nbfluid,nbelast,
     &        nbrigid,maxnth,icase, ncpb,ncpu,ncps,flowtype,bdseg,ndgphi
      real*8 x(nnode),y(nnode),xpos(ncpb,nbrigid),uvom(ncpb,nbrigid)
      real*8 pvis(10),pelas(10),pvelbd(10),pwall(5,bdseg)
      integer ibdu(ncpu,nbd),ibdt(nbd),ibds(nbd),ipnode,ibdp(nbd),
     &        ibdn(nbd)
      real*8 densu(ncpu,nbd),denst(nbd),denss(ncps,nbd),densp(nbd),
     &        densn(nbd),bdpre(nic)
      real*8 phi(ndgphi),time
c
      integer nwall,i,k,n,ncur,ivisc
      real*8 uwall1,uwall2,uin,upmax,speedx,speedy
      real*8 ybar,ybar1,dudy, ueta(4,2*maxnth+1),etath,wdth
      real*8 visc,lambd1, x0,xcoor
      real(8) mup,eps,lam,taupxx,taupyy,taupxy,taupzz,f1,f2
      integer ftype1,ftype2
c
c     Initialization
c     --------------
      visc = pvis(1)
      lambd1= pelas(1)*(1.0-pelas(6))
c
      uwall1 = pvelbd(1)
      uwall2 = pvelbd(2)
      uin    = pvelbd(3)
      upmax  = pvelbd(4)
c
      call INIBD (nbd, ibdu,ibdt,ibds, ibdp, densu,denst,denss,densp,
     &            ncpu,ncps )
      ibdn(1:nbd)=0
      densn(1:nbd)=0.d0

c
c     boundary conditions on the outer-boundary
c     -----------------------------------------
      ftype1=flowtype/100
      ftype2=mod(flowtype,100)
      if (ftype1.eq.0)then
      if (ftype2.eq.1) then
c     flow in a periodic channel
c
c        boundary section 1: y=0
         do i=ic(1),ic(2)-1
            ibdu(1,i) = 1
            ibdu(2,i) = 1
            densu(1,i) = uwall1
            densu(2,i) = 0.0
         enddo
c
c        boundary section 2: y=wdth
         do i=ic(2),ic(3)-1
            ibdu(1,i) = 1
            ibdu(2,i) = 1
            densu(1,i) = uwall2
            densu(2,i) = 0.0
         enddo
c
      elseif ( ftype2.eq.2 .or. ftype2.eq.4 ) then
c     flow in an infinite / a fluidized channel
c
         wdth = pwall(3,1)
c
c        boundary section 2: y=0
         do i=ic(2),ic(3)
            ibdu(1,i) = 1
            ibdu(2,i) = 1
            densu(1,i) = uwall1
            densu(2,i) = 0.0
         enddo
c
c        boundary section 4: y=wdth
         do i=ic(4),ic(5)-1
            ibdu(1,i) = 1
            ibdu(2,i) = 1
            densu(1,i) = uwall2
            densu(2,i) = 0.0
         enddo
         i=ic(1)
         ibdu(1,i) = 1
         ibdu(2,i) = 1
         densu(1,i) = uwall2
         densu(2,i) = 0.0
c
c        boundary section 1: x=xmin ( outflow)
         if ( icase.eq.1 .and.
     &        dabs(uwall1)+dabs(uwall2)+dabs(upmax).gt.1.d-10 ) then
            do i=ic(1),ic(2)
               ibdu(2,i) = 1
            enddo
         endif
         if ( icase.eq.4 ) then
            do i=ic(1),ic(2)
               ncur   = ibdnod(i)
               ybar = y(ncur)/wdth
               ybar1 = 1.d0 - ybar
               speedx = uin + uwall1*ybar1 + uwall2*ybar +
     &                  4.0*upmax*ybar*ybar1
               dudy = (uwall2-uwall1)/wdth + 4.0*upmax/wdth*(1.-2.*ybar)
               etath= visc
               if ( ivisc.ne.1 .and. upmax.ne.0.d0 )
     &            call thin(y(ncur),ueta,maxnth,speedx,dudy,etath)
               ibdu(1,i) = 1
               ibdu(2,i) = 1
               densu(1,i) = speedx
               densu(2,i) = 0.0
            enddo
         endif
c
c        boundary section 3: x=xmax (inflow)
         speedy = 0.0
         do i=ic(3),ic(4)
            ncur = ibdnod(i)
            ybar = y(ncur)/wdth
            ybar1 = 1.d0 - ybar
            speedx = uin+uwall1*ybar1+uwall2*ybar+4*upmax*ybar*ybar1
            dudy = (uwall2-uwall1)/wdth + 4.0*upmax/wdth*(1.-2.*ybar)
            etath= visc
            if ( ivisc.ne.1 .and. upmax.ne.0.d0 )
     &         call thin(y(ncur),ueta,maxnth,speedx,dudy,etath)
c            ibdu(1,i) = 1
c            ibdu(2,i) = 1
c            densu(1,i) = speedx
c            densu(2,i) = speedy
            if ( ibds(i).eq.0 ) then
               ibds(i) = 1
               denss(1,i) = 2.0*lambd1*etath*dudy*dudy
               denss(2,i) = 0.0
               denss(3,i) = 0.0
            endif
         enddo
c
      elseif ( ftype2.eq.3 ) then
c     flow in a closed channel
c
         do i=ic(1),ic(5)-1
            ibdu(1,i) = 1
            ibdu(2,i) = 1
         enddo
c*CH cavity flow bc
c         do i=ic(4),ic(5)-1
c           densu(1,i)=uwall2
c         enddo
c*CH cavity flow bc
c
      elseif ( ftype2.eq.11 ) then
c     flow in a rotating cylinder
c
         do i=ic(1),ic(2)-1
            ibdu(1,i) = 1
            ibdu(2,i) = 1
            ncur = ibdnod(i)
            densu(1,i) = -uwall1*(y(ncur)-pwall(3,1))
            densu(2,i) =  uwall1*(x(ncur)-pwall(2,1))
         enddo
c
      elseif ( ftype2.eq.12 ) then
c     flow between two rotating cylinders
c
         do i=ic(1),ic(2)-1
            ibdu(1,i) = 1
            ibdu(2,i) = 1
            ncur = ibdnod(i)
            densu(1,i) = -uwall1*(y(ncur)-pwall(3,1))
            densu(2,i) =  uwall1*(x(ncur)-pwall(2,1))
         enddo
c
         do i=ic(2),ic(3)-1
            ibdu(1,i) = 1
            ibdu(2,i) = 1
            ncur = ibdnod(i)
            densu(1,i) = -uwall2*(y(ncur)-pwall(3,2))
            densu(2,i) =  uwall2*(x(ncur)-pwall(2,2))
         enddo
c
      elseif ( ftype2.eq.99 ) then
c     flow in an arbitrary geometry
c
         n = 0
         do k=1,bdseg
         if ( mod(int(pwall(1,k)),100).ne.0 ) then
             n = n+1
             if ( mod(int(pwall(1,k)),100).eq.10 ) then
c               boundary type 10 ( inflow)
                wdth = dsqrt( (pwall(4,k)-pwall(2,k))**2 +
     &                        (pwall(5,k)-pwall(3,k))**2 )
                do i=ic(n),ic(n+1)-1
                   ncur = ibdnod(i)
                   ybar = dsqrt( (x(ncur)-pwall(2,k))**2 +
     &                           (y(ncur)-pwall(3,k))**2 ) /wdth
                   ybar1 = 1.d0 - ybar
                   ibdu(1,i) = 1
                   ibdu(2,i) = 1
                   densu(1,i) = pvelbd(1)+4.0*pvelbd(2)*ybar*ybar1
                   densu(2,i) = pvelbd(3)+4.0*pvelbd(4)*ybar*ybar1
                enddo
             elseif ( mod(int(pwall(1,k)),100).eq.11 ) then
c               no-slip velocity  boundary condition
                do i=ic(n),ic(n+1)-1
                  ibdu(1,i) = 1
                  ibdu(2,i) = 1
                enddo
             elseif ( mod(int(pwall(1,k)),100).eq.22 ) then
c               boundary type 22 ( outflow - initialized)
             elseif ( mod(int(pwall(1,k)),100).eq.12 ) then
                do i=ic(n),ic(n+1)-1
                  ibdu(1,i) = 1
                enddo
             elseif ( mod(int(pwall(1,k)),100).eq.21 ) then
                do i=ic(n),ic(n+1)-1
                  ibdu(2,i) = 1
                enddo
             endif
         endif
         enddo
c
      endif
      elseif(ftype1.eq.1)then
c     GRUMMP mesh
      if(ftype2.eq.1)then
c     rectangular domain with solid walls
         do i=ic(1),ic(5)-1
            ibdu(1,i)=1
            ibdu(2,i)=1
         enddo
      elseif(ftype2.eq.2)then
         ipnode=(ic(2)+ic(3))/4*2+1
c     rectangular domain with shear flow
         wdth=pwall(3,1)
c     lower wall
         do i=ic(1),ic(2)-1
            ibdu(1:2,i)=1
            densu(1,i)=uwall1
            densu(2,i)=0.d0
         enddo
c     upper wall
         do i=ic(3),ic(4)-1
            ibdu(1:2,i)=1
            densu(1,i)=uwall2
            densu(2,i)=0.d0
         enddo
c     right boundary
         do i=ic(2),ic(3)-1
            ibdu(1:2,i)=1
            densu(1,i)=uwall1+(uwall2-uwall1)/wdth*y(ibdnod(i))
            densu(2,i)=0.d0
         enddo
c     left boundary
         do i=ic(4),ic(5)-1
            ibdu(1:2,i)=1
            densu(1,i)=uwall1+(uwall2-uwall1)/wdth*y(ibdnod(i))
            densu(2,i)=0.d0
         enddo
         if(icase.eq.4)then
c     impose bc for extra stress(Only works for Oldroyd-B fluid)
c     right boundary
            mup=pvis(1)*(1.-pelas(6))
            eps=abs(uwall2-uwall1)/wdth
            lam=pelas(1)
            taupxx=(lam*eps-eps*(time+lam)*exp(-time/lam))*2*mup*eps
            taupyy=0.d0
            taupxy=(1.-exp(-time/lam))*mup*eps
c     right boundary
            i=(ic(2)+ic(3))/2
            if(ibdnod(i).le.ndgphi)then
               f1=phi(ibdnod(i))
            else
               f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
            endif
c            f1=nint((f1+1.)/2.)
            if(f1.gt.0)then
               f1=1.0
            else
               f1=0.0
            endif
            do i=ic(2),ic(3)
               if(densu(1,i).le.0.0)then
c                  if(i.le.ndgphi)then
c                     f1=phi(ibdnod(i))
c                  else
c                     f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
c                  endif
c                  f1=(1.+f1)/2.
                  ibds(i)=1
                  denss(1,i)=taupxx*f1
                  denss(2,i)=taupyy*f1
                  denss(3,i)=taupxy*f1
               endif
            enddo
c     left boundary
            i=(ic(4)+ic(5))/2
            if(ibdnod(i).le.ndgphi)then
               f1=phi(ibdnod(i))
            else
               f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
            endif
c            f1=nint((f1+1.)/2.)
            if(f1.gt.0)then
               f1=1.0
            else
               f1=0.0
            endif
            do i=ic(4),ic(5)-1
               if(densu(1,i).ge.0.0)then
c                  if(i.le.ndgphi)then
c                     f1=phi(ibdnod(i))
c                  else
c                     if(i.eq.ic(5)-1)then
c                        f1=(phi(ibdnod(i-1))+phi(ic(1)))/2.
c                     else
c                     f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
c                  endif
c                  endif
c                  f1=(1.+f1)/2.
                  ibds(i)=1
                  denss(1,i)=taupxx*f1
                  denss(2,i)=taupyy*f1
                  denss(3,i)=taupxy*f1
               endif
            enddo
            i=ic(1)
            if(densu(1,i).ge.0.0)then
               ibds(i)=1
               denss(1,i)=taupxx*f1
               denss(2,i)=taupyy*f1
               denss(3,i)=taupxy*f1
            endif
	   endif
      elseif(ftype2.eq.3)then
c     planar extensional flow
         ipnode=(ic(3)+ic(2))/4*2+1
c     extensional rate
         eps=pvelbd(1)
c     lower symmetry line
         do i=ic(1),ic(2)-1
            ibdu(2,i)=1
         enddo
c     lower boundary
c         do i=ic(1),ic(2)
c            ibdu(2,i)=1
c            ibdu(2,i)=1
c            ncur = ibdnod(i)
c            densu(2,i)=eps*(y(ncur)-wdth/2.)
c         enddo
c     left symmetry line
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
         enddo
         ibdu(1,ic(1))=1
c     right boundary
         do i=ic(2),ic(3)
            ibdu(1:2,i)=1
            ncur = ibdnod(i)
            densu(1,i)=-eps*x(ncur)
            densu(2,i)=eps*y(ncur)
         enddo
c     top boundary
         do i=ic(3),ic(4)
            ibdu(1:2,i)=1
            ncur = ibdnod(i)
            densu(1,i)=-eps*x(ncur)
            densu(2,i)=eps*y(ncur)
         enddo
         if(icase.eq.4)then
c     extra stress BC(only works for Oldroyd-B fluid)
            mup=pvis(1)*(1.-pelas(6))
            lam=pelas(1)
            taupxx=-2.*(mup*eps)/(1.+2.*lam*eps)
     &             *(1.-exp(-time/lam*(1.+2.*lam*eps)))
            taupyy= 2.*(mup*eps)/(1.-2.*lam*eps)
     &             *(1.-exp(-time/lam*(1.-2.*lam*eps)))
            taupxy=0.d0
            i=(ic(2)+ic(3))/2
            if(ibdnod(i).le.ndgphi)then
               f1=phi(ibdnod(i))
            else
               f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
            endif
c            f1=nint((f1+1.)/2.)
            if(f1.gt.0)then
            do i=ic(2),ic(3)
c               if(i.le.ndgphi)then
c                  f1=phi(ibdnod(i))
c               else
c                  f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
c               endif
c               f1=(1.+f1)/2.
               ibds(i)   =1
               denss(1,i)=taupxx*f1
               denss(2,i)=taupyy*f1
               denss(3,i)=taupxy*f1
            enddo
            endif
         endif
      elseif(ftype2.eq.4)then
c     drop retraction (2D Planar)
         ipnode=ic(3)
c     non-slip walls and symmetric boundaries
c         do i=ic(1),ic(2)-1
c            ibdu(2,i)=1
c         enddo
c         do i=ic(2),ic(4)
c            ibdu(1:2,i)=1
c         enddo
c         do i=ic(4)+1,ic(5)-1
c            ibdu(1,i)=1
c         enddo
c         ibdu(1,ic(1))=1

c     zero normal stress at the open boundaries
         do i=ic(1),ic(2)
            ibdu(2,i)=1
         enddo
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
         enddo
         ibdu(1,ic(1))=1
         if(icase.eq.4)then
            ibds(ic(2):ic(4))=1
         endif
      elseif(ftype2.eq.11)then
c     axisymmetric domain(rectangular) with solid walls
         do i=ic(1),ic(4)-1
            ibdu(1:2,i)=1
         enddo
c     axis
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
c            ibdu(2,i)=0
         enddo
      elseif(ftype2.eq.13)then
c     uniaxial extensional flow(ipnode is odd, because nrdpres=1,nrdgeom=2)
         ipnode=(ic(3)+ic(2))/4*2+1
c     extensional rate
         eps=pvelbd(1)
c     lower symmetry line
         do i=ic(1),ic(2)-1
            ibdu(2,i)=1
         enddo
         ibdu(1,ic(1))=1
c     axis of axisymmetry
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
         enddo
c     right and top bc
c         do i=ic(2),ic(4)
c            ibdu(1:2,i)=1
c            ncur = ibdnod(i)
c            densu(1,i)=-eps/2.*x(ncur)
c            densu(2,i)=eps*y(ncur)
c         enddo
         do i=ic(2),ic(3)
            ibdu(1:2,i)=1
            ncur = ibdnod(i)
            densu(1,i)=-eps/2.*x(ncur)
            densu(2,i)=eps*y(ncur)
         enddo
         do i=ic(3),ic(4)
            ibdu(1:2,i)=1
            ncur = ibdnod(i)
            densu(1,i)=-eps/2.*x(ncur)
            densu(2,i)=eps*y(ncur)
         enddo
         if(icase.eq.4)then
c     extra stress BC(only works for Oldroyd-B fluid)
            mup=pvis(1)*(1.-pelas(6))
            lam=pelas(1)
            taupxx=-(mup*eps)/(1.+lam*eps)
     &             *(1.-exp(-time/lam*(1.+lam*eps)))
            taupyy=2.*(mup*eps)/(1.-2.*lam*eps)
     &             *(1.-exp(-time/lam*(1.-2.*lam*eps)))
            taupxy=0.d0
            taupzz=-(mup*eps)/(1.+lam*eps)
     &             *(1.-exp(-time/lam*(1.+lam*eps)))
            i=(ic(2)+ic(3))/2
            if(ibdnod(i).le.ndgphi)then
               f1=phi(ibdnod(i))
            else
               f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
            endif
c            f1=nint((f1+1.)/2.)
c            if(f1.gt.0)then
c               f1=1.0
c            else
c               f1=0.0
c            endif
            if(f1.gt.0.0)then
            do i=ic(2),ic(3)
c               if(i.le.ndgphi)then
c                  f1=phi(ibdnod(i))
c               else
c                  f1=(phi(ibdnod(i-1))+phi(ibdnod(i+1)))/2.
c               endif
c               f1=(1.+f1)/2.
               ibds(i)=1
               denss(1,i)=taupxx*f1
               denss(2,i)=taupyy*f1
               denss(3,i)=taupxy*f1
               denss(4,i)=taupzz*f1
            enddo
            endif
         endif
      elseif(ftype2.eq.14)then
c     drop retraction (2D axisymmetric)
         ipnode=ic(3)
c     non-slip walls and symmetric boundaries
c         do i=ic(1),ic(2)-1
c            ibdu(2,i)=1
c         enddo
c         do i=ic(2),ic(4)
c            ibdu(1:2,i)=1
c         enddo
c         do i=ic(4)+1,ic(5)-1
c            ibdu(1,i)=1
c         enddo
c         ibdu(1,ic(1))=1

c     zero normal stress at the open boundaries
         do i=ic(1),ic(2)
            ibdu(2,i)=1
         enddo
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
         enddo
         ibdu(1,ic(1))=1
         if(icase.eq.4)then
            ibds(ic(2):ic(4))=1
         endif
      elseif(ftype2.eq.5.or.ftype2.eq.15)then
c     drop head-on collision  (only one quadrant of the domain)
         ipnode=ic(3)
         do i=ic(1),ic(2)
            ibdu(2,i)=1
         enddo
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
         enddo
         ibdu(1,ic(1))=1
         do i=ic(2),ic(4)
            ibdu(1:2,i)=1
         enddo
c         if(icase.eq.4)then
c            ibds(ic(2):ic(4))=1
c         endif
      elseif(ftype2.eq.6.or.ftype2.eq.16)then
c     drop/interface impact
         ipnode=ic(4)
         do i=ic(1),ic(4)-1
            ibdu(1:2,i)=1
         enddo
c     line of symmetry(6) or axisymmetry(16)
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
c            ibdu(2,i)=0
         enddo
      elseif(ftype2.eq.7.or.ftype2.eq.17)then
c     drop/interface impact, upper bound is open
         ipnode=-1
         ibdp(ic(3):ic(4))=1
         densp(ic(3):ic(4))=0.d0
         do i=ic(1),ic(3)
            ibdu(1:2,i)=1
         enddo
c     line of symmetry(7) or axisymmetry(17)
         do i=ic(4),ic(5)-1
            ibdu(1,i)=1
c            ibdu(2,i)=0
         enddo
c     set zero polymer stress at the open boundary
         if(icase.eq.4)then
            ibds(ic(3):ic(4))=1
         endif
      elseif(ftype2==21)then
c     couette device
         ipnode=1
         uwall1=pvelbd(1) !vel on left wall
         uwall2=-uwall1
         wdth=pwall(2,1)-pwall(1,1)
         do i=ic(1),ic(2)-1
               ibdu(1,i)=1
               ibdu(2,i)=1
               densu(2,i)=uwall1+(uwall2-uwall1)*(x(ibdnod(i))/wdth)
         enddo
         do i=ic(2),ic(3)
               ibdu(1,i)=1
               ibdu(2,i)=1
               densu(2,i)=uwall2
         enddo
         do i=ic(3)+1,ic(4)-1
               ibdu(1,i)=1
               ibdu(2,i)=1
               densu(2,i)=uwall1+(uwall2-uwall1)*(x(ibdnod(i))/wdth)
         enddo
         do i=ic(4),ic(5)-1
               ibdu(1,i)=1
               ibdu(2,i)=1
               densu(2,i)=uwall1
         enddo
         if(icase==4)then
            mup=pvis(1)*(1.-pelas(6))
            lam=pelas(1)
            eps=(uwall2-uwall1)/wdth
            f1=(lam-(time+lam)*exp(-time/lam))*2.*mup*eps**2
            f2=(1-exp(-time/lam))*mup*eps
            if(phi(ibdnod(ic(1)))>0)then
c           boundary section one is viscoelastic
               do i=ic(1),ic(2)
                  ibds(i)=1
                  denss(2,i)=f1
                  denss(3,i)=f2
               enddo
            endif
            if(phi(ibdnod(ic(4)))>0)then
c           boundary section three is viscoelastic
               do i=ic(3),ic(4)
                  ibds(i)=1
                  denss(2,i)=f1
                  denss(3,i)=f2
               enddo
            endif
         endif
      elseif(ftype2==22)then
c     plunging tape or capillary tube
         ipnode=1
         uwall1=-pvelbd(1)  !vel on right wall
         wdth=pwall(2,1)-pwall(1,1)
         if(iax==0)then
            uwall2=-0.5d0*uwall1   !vel at center line
         else
            uwall2=-uwall1
         endif
c
         do i=ic(1),ic(2)-1
               ibdu(1,i)=1
               ibdu(2,i)=1
               densu(2,i)=uwall1
     &                   +(uwall2-uwall1)*(1.d0-(x(ibdnod(i))/wdth)**2)
         enddo
c
         do i=ic(2),ic(3)
               ibdu(1,i)=1
               ibdu(2,i)=1
               densu(2,i)=uwall1
         enddo
c
         do i=ic(3)+1,ic(4)
               ibdu(1,i)=1
               ibdu(2,i)=1
               densu(2,i)=uwall1
     &                   +(uwall2-uwall1)*(1.d0-(x(ibdnod(i))/wdth)**2)
         enddo
c
         do i=ic(4),ic(5)-1
               ibdu(1,i)=1
         enddo
         if(icase==4)then
            mup=pvis(1)*(1.-pelas(6))
            lam=pelas(1)
            f1=(lam-(time+lam)*exp(-time/lam))*2.*mup
            f2=(1-exp(-time/lam))*mup
            if(phi(ibdnod(ic(1)))>0)then
c           boundary section one is viscoelastic
               do i=ic(1),ic(2)
                  eps=-2.*(uwall2-uwall1)*x(ibdnod(i))/wdth**2
                  ibds(i)=1
                  denss(2,i)=f1*eps**2
                  denss(3,i)=f2*eps
               enddo
            endif
            if(phi(ibdnod(ic(4)))>0)then
c           boundary section three is viscoelastic
               do i=ic(3),ic(4)
                  eps=-2.*(uwall2-uwall1)*x(ibdnod(i))/wdth**2
                  ibds(i)=1
                  denss(2,i)=f1*eps**2
                  denss(3,i)=f2*eps
               enddo
            endif
         endif
      elseif(ftype2==23)then
         ipnode=0
         do i=ic(2),ic(3)
               ibdu(1,i)=1
               ibdu(2,i)=1
         enddo
c
         do i=ic(4),ic(5)-1
               ibdu(1,i)=1
         enddo
         ibdu(1,ic(1))=1
      elseif(ftype2==24)then
c     capillary spreading
         ipnode=0
c
         do i=ic(1),ic(2)
               ibdu(1,i)=1
               ibdu(2,i)=1
         enddo
c
         do i=ic(4),ic(5)-1
               ibdu(1,i)=1
         enddo
c
      elseif(ftype2==29)then
c
c     drop imbibition in a two-staged tube
c
         ipnode=0
c
         do i=ic(1),ic(2)
            ibdn(i)=1
            densn(i)=-bdpre(1)
         enddo
c
         do i=ic(2),ic(5)
            ibdu(1:2,i)=1
         enddo
c
         do i=ic(5),ic(6)
            ibdn(i)=1
            densn(i)=-bdpre(5)
         enddo
c
         do i=ic(6),ic(7)-1
            ibdu(1,i)=1
         enddo
         ibdu(1,1)=1

      endif

      endif
c
c     boundary conditions on rigid particle surfaces
c     ----------------------------------------------
      nwall = nic-nbrigid-nbfluid-nbelast
      do k=1,nbrigid
         x0 = xpos(1,k)
         do i=ic(k+nwall),ic(k+nwall+1)-1
            ncur = ibdnod(i)
            ibdu(1,i) = -k
            ibdu(2,i) = -k
            densu(1,i) = uvom(1,k)-uvom(3,k)*(y(ncur)-xpos(2,k))
            densu(2,i) = uvom(2,k)+uvom(3,k)*(xcoor(x0,x(ncur))-x0)
         enddo
      enddo
c
c       do i=1,nbd
c          write(*,'(i5,2e15.5)') i,densu(1,i),densu(2,i)
c       enddo
c       pause 'end of boundc'
c
      return
      end
c
c**********************************************************************
      subroutine bdmvel ( ubdmsh,vbdmsh,mshslpp,mshslpw,nbd,ncpb,
     &                    nbrigid,nbfluid,nbelast,
     &                    nic,ic,nnode,x,y,ibdnod,u,v,uvom,xpos)
c
c   this routine imposes boundary conditions for the mesh velocity
c**********************************************************************
      implicit none
      integer nbd,ncpb,nbrigid,nbfluid,nbelast,
     &        nic,nnode,ic(nic+1),ibdnod(nbd)
      double precision ubdmsh(nbd),vbdmsh(nbd),x(nnode),y(nnode)
      double precision u(nnode),v(nnode)
      double precision uvom(ncpb,nbrigid),xpos(ncpb,nbrigid)
      logical mshslpw,mshslpp,xperid

      double precision uxmsh,xlength,ufront,ubehind,posmax,posmin
      integer i,k,n,nwall,nd
      common /xperid/ xlength,xperid
c
      nwall = nic-nbrigid-nbfluid-nbelast
      if ( mshslpw ) then
c        slips on the wall
         if ( xperid ) then
            uxmsh = 0.d0
            do i=ic(nwall+1),ic(nic+1)-1
               uxmsh = uxmsh + u(ibdnod(i))
            enddo
            uxmsh = uxmsh/(ic(nic+1)-ic(nwall+1))
            do i=1,ic(nwall+1)-1
               ubdmsh(i) = uxmsh
               vbdmsh(i) = 0.d0
            enddo
         else
c           for infinite channel
            ufront  = uvom(1,1)
            ubehind = uvom(1,1)
            posmax  = xpos(1,1)
            posmin  = xpos(1,1)
            do i=2,nbrigid
               if ( xpos(1,i) .gt. posmax ) then
                  posmax = xpos(1,i)
                  ufront = uvom(1,i)
               endif
               if ( xpos(1,i) .lt. posmin ) then
                  posmin = xpos(1,i)
                  ubehind = uvom(1,i)
               endif
            enddo
c
            uxmsh = (ufront-ubehind)/(posmax-posmin+1.d-30)
            do i=1,ic(nwall+1)-1
               nd = ibdnod(i)
               if ( x(nd).ge.posmax ) then
                  ubdmsh(i) = ufront
               elseif ( x(nd).le.posmin ) then
                  ubdmsh(i) = ubehind
               else
                  ubdmsh(i) = uxmsh*(x(nd)-posmin) + ubehind
               endif
               vbdmsh(i) = 0.d0
            enddo
         endif
      else
c        no-slip on the wall
         do i=1,ic(nwall+1)-1
            ubdmsh(i) = 0.d0
            vbdmsh(i) = 0.d0
         enddo
      endif
c
      if ( mshslpp ) then
         do k=1,nbrigid
            do i=ic(k+nwall),ic(k+nwall+1)-1
               ubdmsh(i) = uvom(1,k)
               vbdmsh(i) = uvom(2,k)
            enddo
         enddo
      else
         do k=1,nbrigid
            do i=ic(k+nwall),ic(k+nwall+1)-1
               n = ibdnod(i)
               ubdmsh(i) = u(n)
               vbdmsh(i) = v(n)
            enddo
         enddo
      endif
c
      do i=ic(nwall+nbrigid+1),ic(nic+1)-1
         n = ibdnod(i)
         ubdmsh(i) = u(n)
         vbdmsh(i) = v(n)
      enddo
c
      return
      end
c
c**********************************************************************
      subroutine bdmacl ( ubdmsh,vbdmsh,mshslpp,mshslpw,fxym,nbd,ncpb,
     &                    nbrigid,nbfluid,nbelast,
     &                    nic,ic,ibdnod,uvom,xpos,x,y,nnode,ifixpt)
c
c
c     This routine define boundary conditions for the mesh acceleration
c**********************************************************************
      implicit none
      integer nbd,ncpb,nbrigid,nbfluid,nbelast,
     &        nic,nnode,ic(nic+1),ibdnod(nbd)
      double precision fxym(ncpb,nbrigid),ubdmsh(nbd),vbdmsh(nbd)
      double precision uvom(ncpb,nbrigid),xpos(ncpb,nbrigid),
     &                 x(nnode),y(nnode)
      logical ifixpt(ncpb),mshslpw,mshslpp,xperid
c
      double precision uxmsh,xlength,ufront,ubehind,posmax,posmin
      integer i,k,n,nwall,nd
      common /xperid/ xlength,xperid
      double precision dx,dy,tmp,xcoor
c
      nwall = nic-nbrigid-nbfluid-nbelast
      if ( mshslpw ) then
c        slips on the wall
         if ( xperid ) then
            uxmsh = 0.d0
            do i=1,nbrigid
               uxmsh = uxmsh + fxym(1,i)
            enddo
            if ( nbrigid.gt.0 ) then
               uxmsh = uxmsh/nbrigid
            else
               uxmsh = 0.d0
            endif
            do i=1,ic(nwall+1)-1
               ubdmsh(i) = uxmsh
               vbdmsh(i) = 0.d0
            enddo
         else
c           for infinite channel
            ufront  = fxym(1,1)
            ubehind = fxym(1,1)
            posmax  = xpos(1,1)
            posmin  = xpos(1,1)
            do i=2,nbrigid
               if ( xpos(1,i) .gt. posmax ) then
                  posmax = xpos(1,i)
                  ufront = fxym(1,i)
               endif
               if ( xpos(1,i) .lt. posmin ) then
                  posmin  = xpos(1,i)
                  ubehind = fxym(1,i)
               endif
            enddo
c
            uxmsh = (ufront-ubehind)/(posmax-posmin+1.d-30)
            do i=1,ic(nwall+1)-1
               nd = ibdnod(i)
               if ( x(nd).ge.posmax ) then
                  ubdmsh(i) = ufront
               elseif ( x(nd).le.posmin ) then
                  ubdmsh(i) = ubehind
               else
                  ubdmsh(i) = uxmsh*(x(nd)-posmin) + ubehind
               endif
               vbdmsh(i) = 0.d0
            enddo
         endif
      else
c        no-slip on the wall
         do i=1,ic(nwall+1)-1
            ubdmsh(i) = 0.d0
            vbdmsh(i) = 0.d0
         enddo
      endif
c
      if ( mshslpp ) then
         if ( .not.ifixpt(1) ) then
            do k=1,nbrigid
               do i=ic(k+nwall),ic(k+nwall+1)-1
                  ubdmsh(i) = fxym(1,k)
               enddo
            enddo
         else
            do k=1,nbrigid
               do i=ic(k+nwall),ic(k+nwall+1)-1
                  ubdmsh(i) = 0.d0
               enddo
            enddo
         endif
         if ( .not.ifixpt(2) ) then
            do k=1,nbrigid
               do i=ic(k+nwall),ic(k+nwall+1)-1
                  vbdmsh(i) = fxym(2,k)
               enddo
            enddo
         else
            do k=1,nbrigid
               do i=ic(k+nwall),ic(k+nwall+1)-1
                  vbdmsh(i) = 0.d0
               enddo
            enddo
         endif
      else
         do k=1,nbrigid
            do i=ic(k+nwall),ic(k+nwall+1)-1
               n = ibdnod(i)
               dx = xcoor(xpos(1,k),x(n)) - xpos(1,k)
               dy = y(n) - xpos(2,k)
               tmp = uvom(3,k)*uvom(3,k)
               ubdmsh(i) = fxym(1,k) - fxym(3,k)*dy - tmp*dx
               vbdmsh(i) = fxym(2,k) + fxym(3,k)*dx - tmp*dy
            enddo
         enddo
      endif
c
      do i=ic(nwall+nbrigid+1),ic(nic+1)-1
         ubdmsh(i) = 0.d0
         vbdmsh(i) = 0.d0
      enddo
c
      return
      end
c
c***********************************************************************
       subroutine height (itime,time,nbrigid,xpos,xmax)
c
c     this routine give the height of the fluidized bed
c***********************************************************************
       implicit none
       include '../include/parameter2d.h'
       integer itime,nbrigid,j,k,lim,loclg
       real*8 xpos(ncpb,nbrigid),time,xmax,hfludz
       real*8 tmpos(maxptc),dlarge
c
       lim = min0(nbrigid,10)
       do k=1,nbrigid
          tmpos(k) = xmax-xpos(1,k)
       enddo
c
       do k=1,lim
          dlarge = tmpos(k)
          loclg = k
          do j=k+1,nbrigid
             if (tmpos(j).gt.dlarge) then
                dlarge = tmpos(j)
                loclg = j
             endif
          enddo
          tmpos(loclg) = tmpos(k)
          tmpos(k) = dlarge
       enddo
c
       hfludz = 0.d0
       do k=1,lim
          hfludz = hfludz + tmpos(k)
       enddo
       hfludz = hfludz/lim
c
       open(unit=30,file='hh.hfd',access='append')
          if(itime.eq.1) rewind 30
          write(30,'(i4,2(1pe12.5))') itime-1,time,hfludz
       close(unit=30)
       return
       end

c***********************************************************************
      subroutine initphi(nnode,x,y,phi,dropone,eps,rd,xd,yd,flowtype)
c-----------------------------------------------------------------------
c     This subroutine intialize the phi field.
c        (This subroutine doesn't distinguish nodes and vetices)
c
c     input:
c        nnode:   number of nodes(vertices for p1 elements)
c        x,y:     coordiantes of nodes
c        dropone: whether phi=1 in the drop phase
c        eps:     capillary width
c        rd,xd,yd:drop radius and poistion
c     output:
c        phi:     phi field
c     Note:
c        For drop retraction cases(flowtype=104 or 114)
c           1. drop is centered at (0,0);
c           2. the radius in x, y direction are xd and yd respectively;
c           3. rd is not used
c        For drop/interface impact case(flowtype=106 or 116)
c           1. rd: drop radius
c           2. xd: liquid depth of the liquid reservior
c           3. yd: y coordinate of the drop center
c        For moving contact line problem (flowtype=121,122)
c           1. rd: the angle between interface and x axis (DEGREE)
c           2. xd,yd:   the point around which the interface may rotate
c           phi=1 above interface and -1 below
c     Mar 11, 2005, Pengtao Yue
c-----------------------------------------------------------------------
      implicit none
      real(8), parameter:: pi=3.141592653589793d0
      intent (in) :: nnode,x,y,dropone,eps,rd,xd,yd,flowtype
      intent (out):: phi
      integer nnode,flowtype
      real(8) x(nnode),y(nnode),phi(nnode),eps,rd,xd,yd
      logical dropone
      logical new129
c
      integer i
      real(8) rtmp,rr,rtmp1,rtmp2,nx,ny,h
c
      if(flowtype.eq.104.or.flowtype.eq.114)then
c
c     special settings for drop retraction cases
c
c
        do i=1,nnode
c          rr=sqrt(x(i)**2+y(i)**2)
          rr=sqrt(xd*yd)
          rtmp=(sqrt((x(i)/xd)**2+(y(i)/yd)**2)-1)*rr
          rtmp=rtmp/(sqrt(2.)*eps)
          if(rtmp.gt.100.)then
            phi(i)=1.
          elseif(rtmp.lt.-100.)then
            phi(i)=-1.
          else
            rtmp=exp(rtmp)
            phi(i)=(rtmp-1./rtmp)/(rtmp+1./rtmp)
          endif
        enddo
      elseif(flowtype.eq.106.or.flowtype.eq.116.or.flowtype.eq.107.or.
     &       flowtype.eq.117)then
c
c     special setting for drop/interface impact
c
        do i=1,nnode
          rtmp1=sqrt(x(i)**2+((y(i)-yd))**2)
          rtmp1=(rtmp1-rd)/(sqrt(2.)*eps)
          rtmp2=(y(i)-xd)/(sqrt(2.)*eps)
          if(rtmp1.lt.0.and.rtmp2.lt.0)then
            rtmp=-10.
	    elseif(rtmp1.lt.0)then
            rtmp=rtmp1
          elseif(rtmp2.lt.0)then
            rtmp=rtmp2
          else
            rtmp=min(rtmp1,rtmp2)
          endif
          if(rtmp.ge.10.)then
            phi(i)=1.
          elseif(rtmp.le.-10.)then
            phi(i)=-1.
          else
            rtmp=exp(rtmp)
            phi(i)=(rtmp-1./rtmp)/(rtmp+1./rtmp)
          endif
        enddo
      elseif(flowtype==121.or.flowtype==122)then
         nx=-sin(rd/180*pi)
         ny= cos(rd/180*pi)
         do i=1,nnode
            rtmp=((x(i)-xd)*nx+(y(i)-yd)*ny)/(sqrt(2.d0)*eps)
c          rtmp=max(-100.d0,min(100.d0,rtmp))
          if(rtmp.gt.100.)then
            phi(i)=1.d0
          elseif(rtmp.lt.-100.)then
            phi(i)=-1.d0
          else
            rtmp=exp(rtmp)
            phi(i)=(rtmp-1.d0/rtmp)/(rtmp+1.d0/rtmp)
          endif
         enddo
      elseif(flowtype==123)then
         nx=-sin(rd/180*pi)
         ny= cos(rd/180*pi)
         do i=1,nnode
            rtmp1=(x(i)*nx+(y(i)-xd)*ny)/(sqrt(2.d0)*eps)
            rtmp2=-(-x(i)*nx+(y(i)-yd)*ny)/(sqrt(2.d0)*eps)

c          rtmp=max(-100.d0,min(100.d0,rtmp))
          if(rtmp1>=3. .and. rtmp2>=3.)then
            phi(i)=1.d0
          elseif(rtmp1<=-3. .or. rtmp2<=-3.)then
            phi(i)=-1.d0
          elseif((rtmp1<3..and.rtmp1>-3.).and.
     &           (rtmp2<3..and.rtmp2>-3.))then
            if(abs(rtmp1)<abs(rtmp2))then
               rtmp=rtmp1
            else
               rtmp=rtmp2
            endif
            rtmp=exp(rtmp)
            phi(i)=(rtmp-1.d0/rtmp)/(rtmp+1.d0/rtmp)
          elseif(rtmp1<3. .and. rtmp1>-3.)then
            rtmp=exp(rtmp1)
            phi(i)=(rtmp-1.d0/rtmp)/(rtmp+1.d0/rtmp)
          elseif(rtmp2<3. .and. rtmp2>-3.)then
            rtmp=exp(rtmp2)
            phi(i)=(rtmp-1.d0/rtmp)/(rtmp+1.d0/rtmp)
          endif
        enddo
      elseif(flowtype==129)then
        new129=.true. !new129=.false. may not work properly,
                      !recheck before using
        if(new129)then
        do i=1,nnode
          rtmp=sqrt(x(i)**2+((y(i)-yd))**2)
          rtmp=(rtmp-rd)/(sqrt(2.)*eps)
          if(rtmp.gt.100.)then
            phi(i)=1.
          elseif(rtmp.lt.-100)then
            phi(i)=-1.
          else
            rtmp=exp(rtmp)
            phi(i)=(rtmp-1./rtmp)/(rtmp+1./rtmp)
          endif
        enddo
        else
        h=xd
        do i=1,nnode
          rtmp=sqrt(x(i)**2+((y(i)-yd))**2)
          rtmp=(rtmp-rd)/(sqrt(2.)*eps)
          if(y(i)<h-1.d-10)then
            phi(i)=-1.
          elseif(rtmp.gt.100.)then
            phi(i)=1.
          elseif(rtmp.lt.-100)then
            phi(i)=-1.
          else
            rtmp=exp(rtmp)
            phi(i)=(rtmp-1./rtmp)/(rtmp+1./rtmp)
          endif
c        another interface at y=-1.
          rtmp=-(y(i)-(h-1.))/(sqrt(2.)*eps)
          if(y(i)<h-0.5)then
            if(y(i)<h-1.5)then
               phi(i)=1.
            elseif(rtmp.gt.100.)then
               phi(i)=1.
            elseif(rtmp.lt.-100)then
               phi(i)=-1.
            else
               rtmp=exp(rtmp)
               phi(i)=(rtmp-1./rtmp)/(rtmp+1./rtmp)
            endif
          endif
        enddo
        endif
      else
        do i=1,nnode
c          rtmp=sqrt(((x(i)-xd)*4./5.)**2+((y(i)-yd)*5./4.)**2)
          rtmp=sqrt(((x(i)-xd))**2+((y(i)-yd))**2)
          rtmp=(rtmp-rd)/(sqrt(2.)*eps)
c          rtmp=max(-100.d0,min(100.d0,rtmp))
          if(rtmp.gt.100.)then
            phi(i)=1.
          elseif(rtmp.lt.-100.)then
            phi(i)=-1.
          else
            rtmp=exp(rtmp)
            phi(i)=(rtmp-1./rtmp)/(rtmp+1./rtmp)
          endif
        enddo
      endif
      if(dropone)phi(1:nnode)=-phi(1:nnode)
      return
      end
