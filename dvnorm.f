      double precision function dvnorm (n, v, w)
c***begin prologue  dvnorm
c***subsidiary
c***purpose  weighted root-mean-square vector norm.
c***type      double precision (svnorm-s, dvnorm-d)
c***author  hindmarsh, alan c., (llnl)
c***description
c
c  this function routine computes the weighted root-mean-square norm
c  of the vector of length n contained in the array v, with weights
c  contained in the array w of length n:
c    dvnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
c
c***see also  dlsode
c***routines called  (none)
c***revision history  (yymmdd)
c   791129  date written
c   890501  modified prologue to slatec/ldoc format.  (fnf)
c   890503  minor cosmetic changes.  (fnf)
c   930809  renamed to allow single/double precision versions. (ach)
c***end prologue  dvnorm
c**end
      integer n,   i
      double precision v, w,   sum
      dimension v(n), w(n)
c
c***first executable statement  dvnorm
      sum = 0.0d0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      dvnorm = sqrt(sum/n)
      return
c----------------------- end of function dvnorm ------------------------
      end