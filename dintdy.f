      subroutine dintdy (t, k, yh, nyh, dky, iflag)
c***begin prologue  dintdy
c***subsidiary
c***purpose  interpolate solution derivatives.
c***type      double precision (sintdy-s, dintdy-d)
c***author  hindmarsh, alan c., (llnl)
c***description
c
c  dintdy computes interpolated values of the k-th derivative of the
c  dependent variable vector y, and stores it in dky.  this routine
c  is called within the package with k = 0 and t = tout, but may
c  also be called by the user for any k up to the current order.
c  (see detailed instructions in the usage documentation.)
c
c  the computed values in dky are gotten by interpolation using the
c  nordsieck history array yh.  this array corresponds uniquely to a
c  vector-valued polynomial of degree nqcur or less, and dky is set
c  to the k-th derivative of this polynomial at t.
c  the formula for dky is:
c               q
c   dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c              j=k
c  where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c  the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c  communicated by common.  the above sum is done in reverse order.
c  iflag is returned negative if either k or t is out of bounds.
c
c***see also  dlsode
c***routines called  xerrwd
c***common blocks    dls001
c***revision history  (yymmdd)
c   791129  date written
c   890501  modified prologue to slatec/ldoc format.  (fnf)
c   890503  minor cosmetic changes.  (fnf)
c   930809  renamed to allow single/double precision versions. (ach)
c   010418  reduced size of common block /dls001/. (ach)
c   031105  restored 'own' variables to common block /dls001/, to
c           enable interrupt/restart feature. (ach)
c   050427  corrected roundoff decrement in tp. (ach)
c***end prologue  dintdy
c**end
      integer k, nyh, iflag
      double precision t, yh, dky
      dimension yh(nyh,*), dky(*)
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l,
     2   lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter,
     3   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      common /dls001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   iownd(6), iowns(6),
     3   icf, ierpj, iersl, jcur, jstart, kflag, l,
     4   lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      double precision c, r, s, tp
      character*80 msg
c
c***first executable statement  dintdy
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0d0*uround*sign(abs(tn) + abs(hu), hu)
      if ((t-tp)*(t-tn) .gt. 0.0d0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = ic
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = ic
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   msg = 'dintdy-  k (=i1) illegal      '
      call xerrwd (msg, 30, 51, 0, 1, k, 0, 0, 0.0d0, 0.0d0)
      iflag = -1
      return
 90   msg = 'dintdy-  t (=r1) illegal      '
      call xerrwd (msg, 30, 52, 0, 0, 0, 0, 1, t, 0.0d0)
      msg='      t not in interval tcur - hu (= r1) to tcur (=r2)      '
      call xerrwd (msg, 60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
c----------------------- end of subroutine dintdy ----------------------
      end