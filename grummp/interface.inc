      INTERFACE TO subroutine GRUMMP_INIT 
     & [C,ALIAS:'GRUMMP_INIT__FPiN12_0PdN313_']
     &        (mxnbd,mxvert,mxelem,ngmx,ncmx,
     &         nbd,nvert,nelem,refseg,bdseg,
     &         inod,nec,reft,coorx,coory,
     &         h,dG)
      integer mxnbd    [reference]        
      integer mxvert   [reference]        
      integer mxelem   [reference]        
      integer ngmx     [reference]        
      integer ncmx     [reference]        
      integer nbd      [reference]        
      integer nvert    [reference]        
      integer nelem    [reference]        
      integer refseg   [reference]
      integer bdseg    [reference]        
      integer inod     [reference]        
      integer nec      [reference]        
      integer reft     [reference]
      real*8  coorx    [reference]        
      real*8  coory    [reference]        
      real*8  h        [reference]        
      real*8  dG       [reference]
      END
      
      INTERFACE TO subroutine GRUMMP_ADP
     & [C,ALIAS:'GRUMMP_ADP__FPiN50PdT0N70T6N46T0T6']
     &           (mxnbd,mxvert,mxelem,ngmx,ncmx,
     &            nvertd,phid, nbd,nvert,nelem,
     &            bdseg,refseg, inod,nec,reft,
     &            coorx,coory,h1,h2,h3,lrmsh,dG)
      integer mxnbd    [reference]        
      integer mxvert   [reference]        
      integer mxelem   [reference]        
      integer ngmx     [reference]        
      integer ncmx     [reference]        
      integer nvertd   [reference]        
      real*8  phid     [reference]        
      integer nbd      [reference]        
      integer nvert    [reference]        
      integer nelem    [reference]        
      integer refseg   [reference]
      integer bdseg    [reference]        
      integer inod     [reference]        
      integer nec      [reference]        
      integer reft     [reference]        
      real*8  coorx    [reference]        
      real*8  coory    [reference]        
      real*8  h1       [reference]        
      real*8  h2       [reference]      
      real*8  h3       [reference]      
      integer lrmsh    [reference]      
      real*8  dG       [reference]      

      END