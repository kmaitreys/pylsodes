      subroutine dewset (n, itol, rtol, atol, ycur, ewt)
c***begin prologue  dewset
c***subsidiary
c***purpose  set error weight vector.
c***type      double precision (sewset-s, dewset-d)
c***author  hindmarsh, alan c., (llnl)
c***description
c
c  this subroutine sets the error weight vector ewt according to
c      ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
c  with the subscript on rtol and/or atol possibly replaced by 1 above,
c  depending on the value of itol.
c
c***see also  dlsode
c***routines called  (none)
c***revision history  (yymmdd)
c   791129  date written
c   890501  modified prologue to slatec/ldoc format.  (fnf)
c   890503  minor cosmetic changes.  (fnf)
c   930809  renamed to allow single/double precision versions. (ach)
c***end prologue  dewset
c**end
      integer n, itol
      integer i
      double precision rtol, atol, ycur, ewt
      dimension rtol(*), atol(*), ycur(n), ewt(n)
c
c***first executable statement  dewset
      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
      return
c----------------------- end of subroutine dewset ----------------------
      end