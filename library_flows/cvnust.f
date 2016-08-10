c***********************************************************************
      subroutine cvnust ( intval,ndigit,strval,lenstr )
c
c     This routine builds a string representation of an integer value.
c
c     Input arguments :
c       intval : integer value to represent
c       ndigit : requested number of digits in the string representation
c                (automatically adjusted if zero is specified)
c
c     Output arguments :
c       strval : string representation of the value ( 'XXXX' is re-
c                returned if the string is too short)
c       lenstr : length of the built string
c               (equal to ndigit if successful
c                greater if the given value was higher than expected
c                but small enough with respect to the length of strval)
c***********************************************************************
      integer intval,ndigit,lenstr
      character*(*) strval
c
      integer i,resval,auxval,digcnt,maxdgt,iaux
      character*1 auxchr,auxch2
      character*10 ordstr
      data ordstr/'0123456789'/
c
      do i=1,len(strval)
        strval(i:i) = '0'
      enddo
c
      resval = abs(intval)
      maxdgt = len(strval)
      digcnt = 0
c
      do while ( resval.ne.0 ) 
c
          if (digcnt.ge.maxdgt) then
c             The representation of the spec. value has too many digits
              do i=1,len(strval)
                strval(i:i) = 'X'
              enddo
              lenstr = 0
              return
          endif
c
c         Extraction of the rightmost digit of the residual value
          auxval = resval/10
          digcnt = digcnt + 1
          iaux = resval - auxval*10 + 1
          strval(digcnt:digcnt) = ordstr(iaux:iaux)
          resval = auxval
c
      enddo
c
c     Adaptation of the string length
      if (digcnt.le.ndigit) then
c         The resulting length is the one requested by the user
          lenstr = max(1,ndigit)
      else
c         The requested length is insufficient but may be increased
          lenstr = digcnt
      endif
c
c     Inversion of the built string.
      do i=1,lenstr/2
        auxchr = strval(i:i)
        iaux = lenstr + 1 - i
        auxch2 = strval(lenstr+1-i:lenstr+1-i)
        strval(i:i) = auxch2
        strval(iaux:iaux) = auxchr
      enddo
c
      return
      end
