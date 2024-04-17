      subroutine dsrcms (rsav, isav, job)
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks dls001, dlss01, which are used
c internally by one or more odepack solvers.
c
c rsav = real array of length 224 or more.
c isav = integer array of length 71 or more.
c job  = flag indicating to save or restore the common blocks:
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ils, ilss
      integer i, lenils, leniss, lenrls, lenrss
      double precision rsav,   rls, rlss
      dimension rsav(*), isav(*)
      save lenrls, lenils, lenrss, leniss
      common /dls001/ rls(218), ils(37)
      common /dlss01/ rlss(6), ilss(34)
      data lenrls/218/, lenils/37/, lenrss/6/, leniss/34/
c
      if (job .eq. 2) go to 100
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 15 i = 1,lenrss
 15     rsav(lenrls+i) = rlss(i)
c
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      do 25 i = 1,leniss
 25     isav(lenils+i) = ilss(i)
c
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 115 i = 1,lenrss
 115     rlss(i) = rsav(lenrls+i)
c
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      do 125 i = 1,leniss
 125     ilss(i) = isav(lenils+i)
c
      return
c----------------------- end of subroutine dsrcms ----------------------
      end