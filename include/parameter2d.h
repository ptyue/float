

c     parameters in for partmover2d
c
c     ngmx - maximum number of nodes for geometry in each element
c     ncmx - maximum number of faces in each element
c     nsmx - maximum number of nodes for extra-stress
c     numx - maximum number of nodes for velocity
c     npmx - maximum number of nodes for pressure
c     ntmx - maximum number of nodes for temperature
c     nint - maximum number of numerical integeration nodes
c
      integer ngmx,ncmx,nsmx,numx,npmx,ntmx,nint
      parameter ( ngmx = 6 )
      parameter ( ncmx = 3 )
      parameter ( nsmx = 3 )
      parameter ( numx = 6 )
      parameter ( npmx = 3 )
      parameter ( ntmx = 6 )
      parameter ( nint = 9  )
c
c     ncpu - max. number of velocity componenets
c     ncps - max. number of stress components
c     ncpb - max. number of unknows for each particle
c     ncpg - space dimension
c
      integer ncpu,ncps,ncpb,ncpg
      parameter ( ncpu = 3 )
      parameter ( ncps = 4 )
      parameter ( ncpb = 3 )
      parameter ( ncpg = 2 )
c
      integer maxvrt,maxnbd,maxptc,maxbdp,maxelt,maxnod,maxcoe,
     &        maxvar,maxnz,maxiwk,maxrwk,maxnth
c
      parameter ( maxvrt = 80000 )
      parameter ( maxnbd = 10000 )
      parameter ( maxptc = 101 )
      parameter ( maxcoe = 3 )
      parameter ( maxbdp = maxptc+20)
      parameter ( maxelt = 2*(maxvrt-1) )
      parameter ( maxnod = 4*maxvrt )
c
c     for icase=1, Newtonian fluid
c      parameter ( maxvar = 2*maxnod + maxvrt )
c      parameter ( maxiwk = 40*maxvar )
c      parameter ( maxrwk = 40*maxvar )
c
c     for icase=1, Newtonian fluid with interface
c      parameter ( maxvar = 2*maxnod + maxvrt + maxptc*ncpb + 2*maxnod)
c      parameter ( maxiwk = 140*maxvar )
c      parameter ( maxrwk = 140*maxvar )
c
c     for icase=4, viscoelastic fluid with interface
      parameter ( maxvar = 2*maxnod + maxvrt + 3*maxvrt +maxptc*ncpb
     &                   + 2*maxnod)
      parameter ( maxiwk = 180*maxvar )
      parameter ( maxrwk = 180*maxvar )
c
      parameter ( maxnz = maxvar+maxvrt) !2*maxvar + maxvrt )
      parameter ( maxnth = 500 )
c
c     AFTER changes in this file, please recompile partflow.f!!!


