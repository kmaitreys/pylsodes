      subroutine lsodes (f, neq, y, t, tout, itol, rtol, atol, itask,
     1                 istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, liw, mf
      double precision y, t, tout, rtol, atol, rwork
      dimension neq(*), y(*), rtol(*), atol(*), rwork(lrw), iwork(liw)
! This is the 12 November 2003 version of
! DLSODES: Livermore Solver for Ordinary Differential Equations
!          with general Sparse Jacobian matrix.
!
! This version is in double precision.
!
! DLSODES solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! DLSODES is a variant of the DLSODE package, and is intended for
! problems in which the Jacobian matrix df/dy has an arbitrary
! sparse structure (when the problem is stiff).
!
! Authors:       Alan C. Hindmarsh
!                Center for Applied Scientific Computing, L-561
!                Lawrence Livermore National Laboratory
!                Livermore, CA 94551
! and
!                Andrew H. Sherman
!                J. S. Nolen and Associates
!                Houston, TX 77084
!-----------------------------------------------------------------------
! References:
! 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
!     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
!     North-Holland, Amsterdam, 1983, pp. 55-64.
!
! 2.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
!     Yale Sparse Matrix Package: I. The Symmetric Codes,
!     Int. J. Num. Meth. Eng., 18 (1982), pp. 1145-1151.
!
! 3.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
!     Yale Sparse Matrix Package: II. The Nonsymmetric Codes,
!     Research Report No. 114, Dept. of Computer Sciences, Yale
!     University, 1977.
!-----------------------------------------------------------------------
! Summary of Usage.
!
! Communication between the user and the DLSODES package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  See also the example
! problem (with program and output) following this summary.
!
! A. First provide a subroutine of the form:
!               SUBROUTINE F (NEQ, T, Y, YDOT)
!               DOUBLE PRECISION T, Y(*), YDOT(*)
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  If the problem is nonstiff,
! use a method flag MF = 10.  If it is stiff, there are two standard
! choices for the method flag, MF = 121 and MF = 222.  In both cases,
! DLSODES requires the Jacobian matrix in some form, and it treats this
! matrix in general sparse form, with sparsity structure determined
! internally.  (For options where the user supplies the sparsity
! structure, see the full description of MF below.)
!
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 121), but if this is not feasible, DLSODES will
! compute it internally by difference quotients (MF = 222).
! If you are supplying the Jacobian, provide a subroutine of the form:
!               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
!               DOUBLE PRECISION T, Y(*), IAN(*), JAN(*), PDJ(*)
! Here NEQ, T, Y, and J are input arguments, and the JAC routine is to
! load the array PDJ (of length NEQ) with the J-th column of df/dy.
! I.e., load PDJ(i) with df(i)/dy(J) for all relevant values of i.
! The arguments IAN and JAN should be ignored for normal situations.
! DLSODES will call the JAC routine with J = 1,2,...,NEQ.
! Only nonzero elements need be loaded.  Usually, a crude approximation
! to df/dy, possibly with fewer nonzero elements, will suffice.
!
! D. Write a main program which calls Subroutine DLSODES once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages by
! DLSODES.  On the first call to DLSODES, supply arguments as follows:
! F      = name of subroutine for right-hand side vector f.
!          This name must be declared External in calling program.
! NEQ    = number of first order ODEs.
! Y      = array of initial values, of length NEQ.
! T      = the initial value of the independent variable t.
! TOUT   = first point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = relative tolerance parameter (scalar).
! ATOL   = absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution: actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional inputs used.
! RWORK  = real work array of length at least:
!             20 + 16*NEQ            for MF = 10,
!             20 + (2 + 1./LENRAT)*NNZ + (11 + 9./LENRAT)*NEQ
!                                    for MF = 121 or 222,
!          where:
!          NNZ    = the number of nonzero elements in the sparse
!                   Jacobian (if this is unknown, use an estimate), and
!          LENRAT = the real to integer wordlength ratio (usually 1 in
!                   single precision and 2 in double precision).
!          In any case, the required size of RWORK cannot generally
!          be predicted in advance if MF = 121 or 222, and the value
!          above is a rough estimate of a crude lower bound.  Some
!          experimentation with this size may be necessary.
!          (When known, the correct required length is an optional
!          output, available in IWORK(17).)
! LRW    = declared length of RWORK (in user dimension).
! IWORK  = integer work array of length at least 30.
! LIW    = declared length of IWORK (in user dimension).
! JAC    = name of subroutine for Jacobian matrix (MF = 121).
!          If used, this name must be declared External in calling
!          program.  If not used, pass a dummy name.
! MF     = method flag.  Standard values are:
!          10  for nonstiff (Adams) method, no Jacobian used
!          121 for stiff (BDF) method, user-supplied sparse Jacobian
!          222 for stiff method, internally generated sparse Jacobian
! Note that the main program must declare arrays Y, RWORK, IWORK,
! and possibly ATOL.
!
! E. The output from the first call (or any call) is:
!      Y = array of computed values of y(t) vector.
!      T = corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DLSODES was successful, negative otherwise.
!          -1 means excess work done on this call (perhaps wrong MF).
!          -2 means excess accuracy requested (tolerances too small).
!          -3 means illegal input detected (see printed message).
!          -4 means repeated error test failures (check all inputs).
!          -5 means repeated convergence failures (perhaps bad Jacobian
!             supplied or wrong choice of MF or tolerances).
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
!          -7 means a fatal error return flag came from sparse solver
!             CDRV by way of DPRJS or DSOLSS.  Should never happen.
!          A return with ISTATE = -1, -4, or -5 may result from using
!          an inappropriate sparsity structure, one that is quite
!          different from the initial structure.  Consider calling
!          DLSODES again with ISTATE = 3 to force the structure to be
!          reevaluated.  See the full description of ISTATE below.
!
! F. To continue the integration after a successful return, simply
! reset TOUT and call DLSODES again.  No other parameters need be reset.
!
!-----------------------------------------------------------------------
! Example Problem.
!
! The following is a simple example problem, with the coding
! needed for its solution by DLSODES.  The problem is from chemical
! kinetics, and consists of the following 12 rate equations:
!    dy1/dt  = -rk1*y1
!    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
!                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
!    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
!                + rk11*rk14*y4 + rk12*rk14*y6
!    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
!    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
!    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
!    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
!    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
!    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
!    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
!                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
!                - rk6*y10 - rk9*y10
!    dy11/dt = rk10*y8
!    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
!                - rk15*y2*y12 - rk17*y10*y12
!
! with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
!      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
!      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
!      rk15 = rk17 = 100.0.
!
! The t interval is from 0 to 1000, and the initial conditions
! are y1 = 1, y2 = y3 = ... = y12 = 0.  The problem is stiff.
!
! The following coding solves this problem with DLSODES, using MF = 121
! and printing results at t = .1, 1., 10., 100., 1000.  It uses
! ITOL = 1 and mixed relative/absolute tolerance controls.
! During the run and at the end, statistical quantities of interest
! are printed (see optional outputs in the full description below).
!
!     EXTERNAL FEX, JEX
!     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
!     DIMENSION Y(12), RWORK(500), IWORK(30)
!     DATA LRW/500/, LIW/30/
!     NEQ = 12
!     DO 10 I = 1,NEQ
! 10    Y(I) = 0.0D0
!     Y(1) = 1.0D0
!     T = 0.0D0
!     TOUT = 0.1D0
!     ITOL = 1
!     RTOL = 1.0D-4
!     ATOL = 1.0D-6
!     ITASK = 1
!     ISTATE = 1
!     IOPT = 0
!     MF = 121
!     DO 40 IOUT = 1,5
!       CALL DLSODES (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL,
!    1     ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
!       WRITE(6,30)T,IWORK(11),RWORK(11),(Y(I),I=1,NEQ)
! 30    FORMAT(//' At t =',D11.3,4X,
!    1    ' No. steps =',I5,4X,' Last step =',D11.3/
!    2    '  Y array =  ',4D14.5/13X,4D14.5/13X,4D14.5)
!       IF (ISTATE .LT. 0) GO TO 80
!       TOUT = TOUT*10.0D0
! 40    CONTINUE
!     LENRW = IWORK(17)
!     LENIW = IWORK(18)
!     NST = IWORK(11)
!     NFE = IWORK(12)
!     NJE = IWORK(13)
!     NLU = IWORK(21)
!     NNZ = IWORK(19)
!     NNZLU = IWORK(25) + IWORK(26) + NEQ
!     WRITE (6,70) LENRW,LENIW,NST,NFE,NJE,NLU,NNZ,NNZLU
! 70  FORMAT(//' Required RWORK size =',I4,'   IWORK size =',I4/
!    1   ' No. steps =',I4,'   No. f-s =',I4,'   No. J-s =',I4,
!    2   '   No. LU-s =',I4/' No. of nonzeros in J =',I5,
!    3   '   No. of nonzeros in LU =',I5)
!     STOP
! 80  WRITE(6,90)ISTATE
! 90  FORMAT(///' Error halt.. ISTATE =',I3)
!     STOP
!     END
!
!     SUBROUTINE FEX (NEQ, T, Y, YDOT)
!     DOUBLE PRECISION T, Y, YDOT
!     DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
!    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
!     DIMENSION Y(12), YDOT(12)
!     DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
!    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
!    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
!    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
!    4   RK19/50.0D0/, RK20/50.0D0/
!     YDOT(1)  = -RK1*Y(1)
!     YDOT(2)  = RK1*Y(1) + RK11*RK14*Y(4) + RK19*RK14*Y(5)
!    1           - RK3*Y(2)*Y(3) - RK15*Y(2)*Y(12) - RK2*Y(2)
!     YDOT(3)  = RK2*Y(2) - RK5*Y(3) - RK3*Y(2)*Y(3) - RK7*Y(10)*Y(3)
!    1           + RK11*RK14*Y(4) + RK12*RK14*Y(6)
!     YDOT(4)  = RK3*Y(2)*Y(3) - RK11*RK14*Y(4) - RK4*Y(4)
!     YDOT(5)  = RK15*Y(2)*Y(12) - RK19*RK14*Y(5) - RK16*Y(5)
!     YDOT(6)  = RK7*Y(10)*Y(3) - RK12*RK14*Y(6) - RK8*Y(6)
!     YDOT(7)  = RK17*Y(10)*Y(12) - RK20*RK14*Y(7) - RK18*Y(7)
!     YDOT(8)  = RK9*Y(10) - RK13*RK14*Y(8) - RK10*Y(8)
!     YDOT(9)  = RK4*Y(4) + RK16*Y(5) + RK8*Y(6) + RK18*Y(7)
!     YDOT(10) = RK5*Y(3) + RK12*RK14*Y(6) + RK20*RK14*Y(7)
!    1           + RK13*RK14*Y(8) - RK7*Y(10)*Y(3) - RK17*Y(10)*Y(12)
!    2           - RK6*Y(10) - RK9*Y(10)
!     YDOT(11) = RK10*Y(8)
!     YDOT(12) = RK6*Y(10) + RK19*RK14*Y(5) + RK20*RK14*Y(7)
!    1           - RK15*Y(2)*Y(12) - RK17*Y(10)*Y(12)
!     RETURN
!     END
!
!     SUBROUTINE JEX (NEQ, T, Y, J, IA, JA, PDJ)
!     DOUBLE PRECISION T, Y, PDJ
!     DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
!    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
!     DIMENSION Y(12), IA(*), JA(*), PDJ(12)
!     DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
!    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
!    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
!    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
!    4   RK19/50.0D0/, RK20/50.0D0/
!     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), J
! 1   PDJ(1) = -RK1
!     PDJ(2) = RK1
!     RETURN
! 2   PDJ(2) = -RK3*Y(3) - RK15*Y(12) - RK2
!     PDJ(3) = RK2 - RK3*Y(3)
!     PDJ(4) = RK3*Y(3)
!     PDJ(5) = RK15*Y(12)
!     PDJ(12) = -RK15*Y(12)
!     RETURN
! 3   PDJ(2) = -RK3*Y(2)
!     PDJ(3) = -RK5 - RK3*Y(2) - RK7*Y(10)
!     PDJ(4) = RK3*Y(2)
!     PDJ(6) = RK7*Y(10)
!     PDJ(10) = RK5 - RK7*Y(10)
!     RETURN
! 4   PDJ(2) = RK11*RK14
!     PDJ(3) = RK11*RK14
!     PDJ(4) = -RK11*RK14 - RK4
!     PDJ(9) = RK4
!     RETURN
! 5   PDJ(2) = RK19*RK14
!     PDJ(5) = -RK19*RK14 - RK16
!     PDJ(9) = RK16
!     PDJ(12) = RK19*RK14
!     RETURN
! 6   PDJ(3) = RK12*RK14
!     PDJ(6) = -RK12*RK14 - RK8
!     PDJ(9) = RK8
!     PDJ(10) = RK12*RK14
!     RETURN
! 7   PDJ(7) = -RK20*RK14 - RK18
!     PDJ(9) = RK18
!     PDJ(10) = RK20*RK14
!     PDJ(12) = RK20*RK14
!     RETURN
! 8   PDJ(8) = -RK13*RK14 - RK10
!     PDJ(10) = RK13*RK14
!     PDJ(11) = RK10
! 9   RETURN
! 10  PDJ(3) = -RK7*Y(3)
!     PDJ(6) = RK7*Y(3)
!     PDJ(7) = RK17*Y(12)
!     PDJ(8) = RK9
!     PDJ(10) = -RK7*Y(3) - RK17*Y(12) - RK6 - RK9
!     PDJ(12) = RK6 - RK17*Y(12)
! 11  RETURN
! 12  PDJ(2) = -RK15*Y(2)
!     PDJ(5) = RK15*Y(2)
!     PDJ(7) = RK17*Y(10)
!     PDJ(10) = -RK17*Y(10)
!     PDJ(12) = -RK15*Y(2) - RK17*Y(10)
!     RETURN
!     END
!
! The output of this program (on a Cray-1 in single precision)
! is as follows:
!
!
! At t =  1.000e-01     No. steps =   12     Last step =  1.515e-02
!  Y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
!                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
!                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
!
!
! At t =  1.000e+00     No. steps =   33     Last step =  7.880e-02
!  Y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
!                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
!                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
!
!
! At t =  1.000e+01     No. steps =   48     Last step =  1.239e+00
!  Y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
!                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
!                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
!
!
! At t =  1.000e+02     No. steps =   91     Last step =  3.764e+00
!  Y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
!                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
!                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
!
!
! At t =  1.000e+03     No. steps =  111     Last step =  4.156e+02
!  Y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
!               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
!                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
!
!
! Required RWORK size = 442   IWORK size =  30
! No. steps = 111   No. f-s = 142   No. J-s =   2   No. LU-s =  20
! No. of nonzeros in J =   44   No. of nonzeros in LU =   50
!
!-----------------------------------------------------------------------
! Full Description of User Interface to DLSODES.
!
! The user interface to DLSODES consists of the following parts.
!
! 1.   The call sequence to Subroutine DLSODES, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is a description of
!      optional inputs available through the call sequence, and then
!      a description of optional outputs (in the work arrays).
!
! 2.   Descriptions of other routines in the DLSODES package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      Common, and obtain specified derivatives of the solution y(t).
!
! 3.   Descriptions of Common blocks to be declared in overlay
!      or similar environments, or to be saved when doing an interrupt
!      of the problem and continued solution later.
!
! 4.   Description of two routines in the DLSODES package, either of
!      which the user may replace with his/her own version, if desired.
!      These relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part 1.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RWORK and IWORK are also used for conditional and
! optional inputs and optional outputs.  (The term output here refers
! to the return from Subroutine DLSODES to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 on input.
!
! The descriptions of the call arguments are as follows.
!
! F      = the name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT)
!               DOUBLE PRECISION T, Y(*), YDOT(*)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are arrays of length NEQ.
!          Subroutine F should not alter y(1),...,y(NEQ).
!          F must be declared External in the calling program.
!
!          Subroutine F may access user-defined quantities in
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!          (dimensioned in F) and/or Y has length exceeding NEQ(1).
!          See the descriptions of NEQ and Y below.
!
!          If quantities computed in the F routine are needed
!          externally to DLSODES, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use DINTDY instead.
!
! NEQ    = the size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may be decreased, but not increased, during the problem.
!          If NEQ is decreased (with ISTATE = 3 on input), the
!          remaining components of Y should be left undisturbed, if
!          these are to be accessed in F and/or JAC.
!
!          Normally, NEQ is a scalar, and it is generally referred to
!          as a scalar in this user interface description.  However,
!          NEQ may be an array, with NEQ(1) set to the system size.
!          (The DLSODES package accesses only NEQ(1).)  In either case,
!          this parameter is passed as the NEQ argument in all calls
!          to F and JAC.  Hence, if it is an array, locations
!          NEQ(2),... may be used to store other integer data and pass
!          it to F and/or JAC.  Subroutines F and/or JAC must include
!          NEQ in a Dimension statement in that case.
!
! Y      = a real array for the vector of dependent variables, of
!          length NEQ or more.  Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          on the first call, Y must contain the vector of initial
!          values.  On output, Y contains the computed solution vector,
!          evaluated at T.  If desired, the Y array may be used
!          for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F and JAC.  Hence its length may exceed NEQ, and locations
!          Y(NEQ+1),... may be used to store other real data and
!          pass it to F and/or JAC.  (The DLSODES package accesses only
!          Y(1),...,Y(NEQ).)
!
! T      = the independent variable.  On input, T is used only on the
!          first call, as the initial point of the integration.
!          on output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
!
! TOUT   = the next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial T, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal T interval, whose endpoints are
!          TCUR - HU and TCUR (see optional outputs, below, for
!          TCUR and HU).
!
! ITOL   = an indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = a relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = an absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!             The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector E = (E(i)) of estimated local errors
!          in y, according to an inequality of the form
!                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
!          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
!          and the RMS-norm (root-mean-square norm) here is
!          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation.  See Part 4 below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = an index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at t = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          On input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
!             the conditional inputs IA and JA,
!             and any of the optional inputs except H0.
!             In particular, if MITER = 1 or 2, a call with ISTATE = 3
!             will cause the sparsity structure of the problem to be
!             recomputed (or reread from IA and JA if MOSS = 0).
!          Note:  a preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful for the
!          purpose of outputting the initial conditions.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 on input.
!
!          On output, ISTATE has the following values and meanings.
!           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again
!              (the excess work step counter will be reset to 0).
!              In addition, the user may increase MXSTEP to avoid
!              this error return (see below on optional inputs).
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!          -7  means a fatal error return flag came from the sparse
!              solver CDRV by way of DPRJS or DSOLSS (numerical
!              factorization or backsolve).  This should never happen.
!              The integration was successful as far as T.
!
!          Note: an error return with ISTATE = -1, -4, or -5 and with
!          MITER = 1 or 2 may mean that the sparsity structure of the
!          problem has changed significantly since it was last
!          determined (or input).  In that case, one can attempt to
!          complete the integration by setting ISTATE = 3 on the next
!          call, so that a new structure determination is done.
!
!          Note:  since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other inputs, before
!          calling the solver again.
!
! IOPT   = an integer flag to specify whether or not any optional
!          inputs are being used on this call.  Input only.
!          The optional inputs are listed separately below.
!          IOPT = 0 means no optional inputs are being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means one or more optional inputs are being used.
!
! RWORK  = a work array used for a mixture of real (double precision)
!          and integer work space.
!          The length of RWORK (in real words) must be at least
!             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          LWM = 0                                    if MITER = 0,
!          LWM = 2*NNZ + 2*NEQ + (NNZ+9*NEQ)/LENRAT   if MITER = 1,
!          LWM = 2*NNZ + 2*NEQ + (NNZ+10*NEQ)/LENRAT  if MITER = 2,
!          LWM = NEQ + 2                              if MITER = 3.
!          In the above formulas,
!          NNZ    = number of nonzero elements in the Jacobian matrix.
!          LENRAT = the real to integer wordlength ratio (usually 1 in
!                   single precision and 2 in double precision).
!          (See the MF description for METH and MITER.)
!          Thus if MAXORD has its default value and NEQ is constant,
!          the minimum length of RWORK is:
!             20 + 16*NEQ        for MF = 10,
!             20 + 16*NEQ + LWM  for MF = 11, 111, 211, 12, 112, 212,
!             22 + 17*NEQ        for MF = 13,
!             20 +  9*NEQ        for MF = 20,
!             20 +  9*NEQ + LWM  for MF = 21, 121, 221, 22, 122, 222,
!             22 + 10*NEQ        for MF = 23.
!          If MITER = 1 or 2, the above formula for LWM is only a
!          crude lower bound.  The required length of RWORK cannot
!          be readily predicted in general, as it depends on the
!          sparsity structure of the problem.  Some experimentation
!          may be necessary.
!
!          The first 20 words of RWORK are reserved for conditional
!          and optional inputs and optional outputs.
!
!          The following word in RWORK is a conditional input:
!            RWORK(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot.  Required if ITASK is
!                       4 or 5, and ignored otherwise.  (See ITASK.)
!
! LRW    = the length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = an integer work array.  The length of IWORK must be at least
!             31 + NEQ + NNZ   if MOSS = 0 and MITER = 1 or 2, or
!             30               otherwise.
!          (NNZ is the number of nonzero elements in df/dy.)
!
!          In DLSODES, IWORK is used only for conditional and
!          optional inputs and optional outputs.
!
!          The following two blocks of words in IWORK are conditional
!          inputs, required if MOSS = 0 and MITER = 1 or 2, but not
!          otherwise (see the description of MF for MOSS).
!            IWORK(30+j) = IA(j)     (j=1,...,NEQ+1)
!            IWORK(31+NEQ+k) = JA(k) (k=1,...,NNZ)
!          The two arrays IA and JA describe the sparsity structure
!          to be assumed for the Jacobian matrix.  JA contains the row
!          indices where nonzero elements occur, reading in columnwise
!          order, and IA contains the starting locations in JA of the
!          descriptions of columns 1,...,NEQ, in that order, with
!          IA(1) = 1.  Thus, for each column index j = 1,...,NEQ, the
!          values of the row index i in column j where a nonzero
!          element may occur are given by
!            i = JA(k),  where   IA(j) .le. k .lt. IA(j+1).
!          If NNZ is the total number of nonzero locations assumed,
!          then the length of the JA array is NNZ, and IA(NEQ+1) must
!          be NNZ + 1.  Duplicate entries are not allowed.
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note:  The work arrays must not be altered between calls to DLSODES
! for the same problem, except possibly for the conditional and
! optional inputs, and except for the last 3*NEQ words of RWORK.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside DLSODES between calls, if
! desired (but not for use by F or JAC).
!
! JAC    = name of user-supplied routine (MITER = 1 or MOSS = 1) to
!          compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  It is to have the form
!               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
!               DOUBLE PRECISION T, Y(*), IAN(*), JAN(*), PDJ(*)
!          where NEQ, T, Y, J, IAN, and JAN are input, and the array
!          PDJ, of length NEQ, is to be loaded with column J
!          of the Jacobian on output.  Thus df(i)/dy(J) is to be
!          loaded into PDJ(i) for all relevant values of i.
!          Here T and Y have the same meaning as in Subroutine F,
!          and J is a column index (1 to NEQ).  IAN and JAN are
!          undefined in calls to JAC for structure determination
!          (MOSS = 1).  otherwise, IAN and JAN are structure
!          descriptors, as defined under optional outputs below, and
!          so can be used to determine the relevant row indices i, if
!          desired.
!               JAC need not provide df/dy exactly.  A crude
!          approximation (possibly with greater sparsity) will do.
!               In any case, PDJ is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Calls to JAC are made with J = 1,...,NEQ, in that order, and
!          each such set of calls is preceded by a call to F with the
!          same arguments NEQ, T, and Y.  Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user Common block by F and not recomputed by JAC,
!          if desired.  JAC must not alter its input arguments.
!          JAC must be declared External in the calling program.
!               Subroutine JAC may access user-defined quantities in
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
!          See the descriptions of NEQ and Y above.
!
! MF     = the method flag.  Used only for input.
!          MF has three decimal digits-- MOSS, METH, MITER--
!             MF = 100*MOSS + 10*METH + MITER.
!          MOSS indicates the method to be used to obtain the sparsity
!          structure of the Jacobian matrix if MITER = 1 or 2:
!            MOSS = 0 means the user has supplied IA and JA
!                     (see descriptions under IWORK above).
!            MOSS = 1 means the user has supplied JAC (see below)
!                     and the structure will be obtained from NEQ
!                     initial calls to JAC.
!            MOSS = 2 means the structure will be obtained from NEQ+1
!                     initial calls to F.
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on Backward
!                     Differentiation Formulas (BDFs).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      sparse Jacobian, given by Subroutine JAC.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) sparse Jacobian
!                      (using NGP extra calls to F per df/dy value,
!                      where NGP is an optional output described below.)
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!          If MITER = 1 or MOSS = 1, the user must supply a Subroutine
!          JAC (the name is arbitrary) as described above under JAC.
!          Otherwise, a dummy argument can be used.
!
!          The standard choices for MF are:
!            MF = 10  for a nonstiff problem,
!            MF = 21 or 22 for a stiff problem with IA/JA supplied
!                     (21 if JAC is supplied, 22 if not),
!            MF = 121 for a stiff problem with JAC supplied,
!                     but not IA/JA,
!            MF = 222 for a stiff problem with neither IA/JA nor
!                     JAC supplied.
!          The sparseness structure can be changed during the
!          problem by making a call to DLSODES with ISTATE = 3.
!-----------------------------------------------------------------------
! Optional Inputs.
!
! The following is a list of the optional inputs provided for in the
! call sequence.  (See also Part 2.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of these inputs requires IOPT = 1, and in that
! case all of these inputs are examined.  A value of zero for any
! of these optional inputs will cause the default value to be used.
! Thus to use a subset of the optional inputs, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! Name    Location      Meaning and Default Value
!
! H0      RWORK(5)  the step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  the maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  the minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! SETH    RWORK(8)  the element threshhold for sparsity determination
!                   when MOSS = 1 or 2.  If the absolute value of
!                   an estimated Jacobian element is .le. SETH, it
!                   will be assumed to be absent in the structure.
!                   The default value of SETH is 0.
!
! MAXORD  IWORK(5)  the maximum order to be allowed.  The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
!
! MXSTEP  IWORK(6)  maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!-----------------------------------------------------------------------
! Optional Outputs.
!
! As optional additional output from DLSODES, the variables listed
! below are quantities related to the performance of DLSODES
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of these outputs are defined
! on any successful return from DLSODES, and on any return with
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, outputs relevant to the error will be defined,
! as noted below.
!
! Name    Location      Meaning
!
! HU      RWORK(11) the step size in t last used (successfully).
!
! HCUR    RWORK(12) the step size to be attempted on the next step.
!
! TCUR    RWORK(13) the current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  On output, TCUR
!                   will always be at least as far as the argument
!                   T, but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! NST     IWORK(11) the number of steps taken for the problem so far.
!
! NFE     IWORK(12) the number of f evaluations for the problem so far,
!                   excluding those for structure determination
!                   (MOSS = 2).
!
! NJE     IWORK(13) the number of Jacobian evaluations for the problem
!                   so far, excluding those for structure determination
!                   (MOSS = 1).
!
! NQU     IWORK(14) the method order last used (successfully).
!
! NQCUR   IWORK(15) the order to be attempted on the next step.
!
! IMXER   IWORK(16) the index of the component of largest magnitude in
!                   the weighted local error vector ( E(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENRW   IWORK(17) the length of RWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(18) the length of IWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! NNZ     IWORK(19) the number of nonzero elements in the Jacobian
!                   matrix, including the diagonal (MITER = 1 or 2).
!                   (This may differ from that given by IA(NEQ+1)-1
!                   if MOSS = 0, because of added diagonal entries.)
!
! NGP     IWORK(20) the number of groups of column indices, used in
!                   difference quotient Jacobian aproximations if
!                   MITER = 2.  This is also the number of extra f
!                   evaluations needed for each Jacobian evaluation.
!
! NLU     IWORK(21) the number of sparse LU decompositions for the
!                   problem so far.
!
! LYH     IWORK(22) the base address in RWORK of the history array YH,
!                   described below in this list.
!
! IPIAN   IWORK(23) the base address of the structure descriptor array
!                   IAN, described below in this list.
!
! IPJAN   IWORK(24) the base address of the structure descriptor array
!                   JAN, described below in this list.
!
! NZL     IWORK(25) the number of nonzero elements in the strict lower
!                   triangle of the LU factorization used in the chord
!                   iteration (MITER = 1 or 2).
!
! NZU     IWORK(26) the number of nonzero elements in the strict upper
!                   triangle of the LU factorization used in the chord
!                   iteration (MITER = 1 or 2).
!                   The total number of nonzeros in the factorization
!                   is therefore NZL + NZU + NEQ.
!
! The following four arrays are segments of the RWORK array which
! may also be of interest to the user as optional outputs.
! For each array, the table below gives its internal name,
! its base address, and its description.
! For YH and ACOR, the base addresses are in RWORK (a real array).
! The integer arrays IAN and JAN are to be obtained by declaring an
! integer array IWK and identifying IWK(1) with RWORK(21), using either
! an equivalence statement or a subroutine call.  Then the base
! addresses IPIAN (of IAN) and IPJAN (of JAN) in IWK are to be obtained
! as optional outputs IWORK(23) and IWORK(24), respectively.
! Thus IAN(1) is IWK(IPIAN), etc.
!
! Name    Base Address      Description
!
! IAN    IPIAN (in IWK)  structure descriptor array of size NEQ + 1.
! JAN    IPJAN (in IWK)  structure descriptor array of size NNZ.
!         (see above)    IAN and JAN together describe the sparsity
!                        structure of the Jacobian matrix, as used by
!                        DLSODES when MITER = 1 or 2.
!                        JAN contains the row indices of the nonzero
!                        locations, reading in columnwise order, and
!                        IAN contains the starting locations in JAN of
!                        the descriptions of columns 1,...,NEQ, in
!                        that order, with IAN(1) = 1.  Thus for each
!                        j = 1,...,NEQ, the row indices i of the
!                        nonzero locations in column j are
!                        i = JAN(k),  IAN(j) .le. k .lt. IAN(j+1).
!                        Note that IAN(NEQ+1) = NNZ + 1.
!                        (If MOSS = 0, IAN/JAN may differ from the
!                        input IA/JA because of a different ordering
!                        in each column, and added diagonal entries.)
!
! YH      LYH            the Nordsieck history array, of size NYH by
!          (optional     (NQCUR + 1), where NYH is the initial value
!           output)      of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the solution,
!                        evaluated at t = TCUR.  The base address LYH
!                        is another optional output, listed above.
!
! ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
!                        corrections on each step, scaled on output
!                        to represent the estimated local error in y
!                        on the last step.  This is the vector E  in
!                        the description of the error control.  It is
!                        defined only on a successful return from
!                        DLSODES.
!
!-----------------------------------------------------------------------
! Part 2.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with DLSODES.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!
!     Form of Call                  Function
!   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
!                             output of messages from DLSODES, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!   CALL XSETF(MFLAG)         Set a flag to control the printing of
!                             messages by DLSODES.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!   CALL DSRCMS(RSAV,ISAV,JOB) saves and restores the contents of
!                             the internal Common blocks used by
!                             DLSODES (see Part 3 below).
!                             RSAV must be a real array of length 224
!                             or more, and ISAV must be an integer
!                             array of length 71 or more.
!                             JOB=1 means save Common into RSAV/ISAV.
!                             JOB=2 means restore Common from RSAV/ISAV.
!                                DSRCMS is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with DLSODES.
!
!   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
!        (see below)          orders, at a specified point t, if
!                             desired.  It may be called only after
!                             a successful return from DLSODES.
!
! The detailed instructions for using DINTDY are as follows.
! The form of the call is:
!
!   LYH = IWORK(22)
!   CALL DINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
!
! The input parameters are:
!
! T         = value of independent variable where answers are desired
!             (normally the same as the T last returned by DLSODES).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional outputs for TCUR and HU.)
! K         = integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (See optional outputs).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by DLSODES directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with DINTDY.
! LYH       = the base address of the history array YH, obtained
!             as an optional output as shown above.
! NYH       = column length of YH, equal to the initial value of NEQ.
!
! The output parameters are:
!
! DKY       = a real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part 3.  Common Blocks.
!
! If DLSODES is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in:
!   (1) the call sequence to DLSODES, and
!   (2) the two internal Common blocks
!         /DLS001/  of length  255  (218 double precision words
!                      followed by 37 integer words),
!         /DLSS01/  of length  40  (6 double precision words
!                      followed by 34 integer words),
!
! If DLSODES is used on a system in which the contents of internal
! Common blocks are not preserved between calls, the user should
! declare the above Common blocks in the calling program to insure
! that their contents are preserved.
!
! If the solution of a given problem by DLSODES is to be interrupted
! and then later continued, such as when restarting an interrupted run
! or alternating between two or more problems, the user should save,
! following the return from the last DLSODES call prior to the
! interruption, the contents of the call sequence variables and the
! internal Common blocks, and later restore these values before the
! next DLSODES call for that problem.  To save and restore the Common
! blocks, use Subroutine DSRCMS (see Part 2 above).
!
!-----------------------------------------------------------------------
! Part 4.  Optionally Replaceable Solver Routines.
!
! Below are descriptions of two routines in the DLSODES package which
! relate to the measurement of errors.  Either routine can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) DEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODES call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by DEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
! in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
! routine (see below), and also used by DLSODES in the computation
! of the optional output IMXER, the diagonal Jacobian approximation,
! and the increments for difference quotient Jacobians.
!
! In the user-supplied version of DEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! optional outputs.  In DEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of H**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in DEWSET the statements:
!     DOUBLE PRECISION RLS
!     COMMON /DLS001/ RLS(218),ILS(37)
!     NQ = ILS(33)
!     NST = ILS(34)
!     H = RLS(212)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!
! (b) DVNORM.
! The following is a real function routine which computes the weighted
! root-mean-square norm of a vector v:
!     D = DVNORM (N, V, W)
! where
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by Subroutine DEWSET.
!
! If the user supplies this function, it should return a non-negative
! value of DVNORM suitable for use in the error control in DLSODES.
! None of the arguments should be altered by DVNORM.
! For example, a user-supplied DVNORM routine might:
!   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of y.
!-----------------------------------------------------------------------
!
!***REVISION HISTORY  (YYYYMMDD)
! 19810120  DATE WRITTEN
! 19820315  Upgraded MDI in ODRV package: operates on M + M-transpose.
! 19820426  Numerous revisions in use of work arrays;
!           use wordlength ratio LENRAT; added IPISP & LRAT to Common;
!           added optional outputs IPIAN/IPJAN;
!           numerous corrections to comments.
! 19830503  Added routine CNTNZU; added NZL and NZU to /LSS001/;
!           changed ADJLR call logic; added optional outputs NZL & NZU;
!           revised counter initializations; revised PREP stmt. numbers;
!           corrections to comments throughout.
! 19870320  Corrected jump on test of umax in CDRV routine;
!           added ISTATE = -7 return.
! 19870330  Major update: corrected comments throughout;
!           removed TRET from Common; rewrote EWSET with 4 loops;
!           fixed t test in INTDY; added Cray directives in STODE;
!           in STODE, fixed DELP init. and logic around PJAC call;
!           combined routines to save/restore Common;
!           passed LEVEL = 0 in error message calls (except run abort).
! 20010425  Major update: convert source lines to upper case;
!           added *DECK lines; changed from 1 to * in dummy dimensions;
!           changed names R1MACH/D1MACH to RUMACH/DUMACH;
!           renamed routines for uniqueness across single/double prec.;
!           converted intrinsic names to generic form;
!           removed ILLIN and NTREP (data loaded) from Common;
!           removed all 'own' variables from Common;
!           changed error messages to quoted strings;
!           replaced XERRWV/XERRWD with 1993 revised version;
!           converted prologues, comments, error messages to mixed case;
!           converted arithmetic IF statements to logical IF statements;
!           numerous corrections to prologues and internal comments.
! 20010507  Converted single precision source to double precision.
! 20020502  Corrected declarations in descriptions of user routines.
! 20031105  Restored 'own' variables to Common blocks, to enable
!           interrupt/restart feature.
! 20031112  Added SAVE statements for data-loaded constants.
!
!-----------------------------------------------------------------------
! Other routines in the DLSODES package.
!
! In addition to Subroutine DLSODES, the DLSODES package includes the
! following subroutines and function routines:
!  DIPREP   acts as an iterface between DLSODES and DPREP, and also does
!           adjusting of work space pointers and work arrays.
!  DPREP    is called by DIPREP to compute sparsity and do sparse matrix
!           preprocessing if MITER = 1 or 2.
!  JGROUP   is called by DPREP to compute groups of Jacobian column
!           indices for use when MITER = 2.
!  ADJLR    adjusts the length of required sparse matrix work space.
!           It is called by DPREP.
!  CNTNZU   is called by DPREP and counts the nonzero elements in the
!           strict upper triangle of J + J-transpose, where J = df/dy.
!  DINTDY   computes an interpolated value of the y vector at t = TOUT.
!  DSTODE   is the core integrator, which does one step of the
!           integration and the associated error control.
!  DCFODE   sets all method coefficients and test constants.
!  DPRJS    computes and preprocesses the Jacobian matrix J = df/dy
!           and the Newton iteration matrix P = I - h*l0*J.
!  DSOLSS   manages solution of linear system in chord iteration.
!  DEWSET   sets the error weight vector EWT before each step.
!  DVNORM   computes the weighted RMS-norm of a vector.
!  DSRCMS   is a user-callable routine to save and restore
!           the contents of the internal Common blocks.
!  ODRV     constructs a reordering of the rows and columns of
!           a matrix by the minimum degree algorithm.  ODRV is a
!           driver routine which calls Subroutines MD, MDI, MDM,
!           MDP, MDU, and SRO.  See Ref. 2 for details.  (The ODRV
!           module has been modified since Ref. 2, however.)
!  CDRV     performs reordering, symbolic factorization, numerical
!           factorization, or linear system solution operations,
!           depending on a path argument ipath.  CDRV is a
!           driver routine which calls Subroutines NROC, NSFC,
!           NNFC, NNSC, and NNTC.  See Ref. 3 for details.
!           DLSODES uses CDRV to solve linear systems in which the
!           coefficient matrix is  P = I - con*J, where I is the
!           identity, con is a scalar, and J is an approximation to
!           the Jacobian df/dy.  Because CDRV deals with rowwise
!           sparsity descriptions, CDRV works with P-transpose, not P.
!  DUMACH   computes the unit roundoff in a machine-independent manner.
!  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
!           error messages and warnings.  XERRWD is machine-dependent.
! Note:  DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
! All the others are subroutines.
!
!-----------------------------------------------------------------------
      external dprjs, dsolss
      double precision dumach, dvnorm
!
! Declare all other variables.
      integer init, mxstep, mxhnil, nhnil, nslast, nyh, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l,
     2   lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter,
     3   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem,
     1   j, kgo, lenrat, lenhyt, leniw, lenrw, lf0, lia, lja,
     2   lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, mxstp0, ncolm
      double precision rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision con0, conmin, ccmxj, psmall, rbig, seth
      double precision atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0
      dimension mord(2)
      logical ihit
      character*60 msg
      save lenrat, mord, mxstp0, mxhnl0
!-----------------------------------------------------------------------
! The following two internal Common blocks contain
! (a) variables which are local to any subroutine but whose values must
!     be preserved between calls to the routine ("own" variables), and
! (b) variables which are communicated between subroutines.
! The block DLS001 is declared in subroutines DLSODES, DIPREP, DPREP,
! DINTDY, DSTODE, DPRJS, and DSOLSS.
! The block DLSS01 is declared in subroutines DLSODES, DIPREP, DPREP,
! DPRJS, and DSOLSS.
! Groups of variables are replaced by dummy arrays in the Common
! declarations in routines where those variables are not used.
!-----------------------------------------------------------------------
      common /ls001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   init, mxstep, mxhnil, nhnil, nslast, nyh, iowns(6),
     3   icf, ierpj, iersl, jcur, jstart, kflag, l,
     4   lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
!
      common /lss01/ con0, conmin, ccmxj, psmall, rbig, seth,
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
     
      data mord(1), mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
!-----------------------------------------------------------------------
! In the Data statement below, set LENRAT equal to the ratio of
! the wordlength for a real number to that for an integer.  Usually,
! LENRAT = 1 for single precision and 2 for double precision.  If the
! true ratio is not an integer, use the next smaller integer (.ge. 1).
!-----------------------------------------------------------------------
      data lenrat/2/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) return
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all inputs and various initializations.
! If ISTATE = 1, the final setting of work space pointers, the matrix
! preprocessing, and other initializations are done in Block C.
!
! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
! MF, ML, and MU.
!-----------------------------------------------------------------------
 20   if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      moss = mf/100
      mf1 = mf - 100*moss
      meth = mf1/10
      miter = mf1 - 10*meth
      if (moss .lt. 0 .or. moss .gt. 2) go to 608
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 3) go to 608
      if (miter .eq. 0 .and. miter .eq. 3) moss = 0
! Next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0D0
      hmxi = 0.0D0
      hmin = 0.0D0
      seth = 0.0D0
      go to 68
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout -t)*h0 .lt. 0.0D0) go to 614
 50   hmax = rwork(6)
      if  (hmax .lt. 0.0D0) go to 615
      hmxi = 0.0D0
      if (hmax .gt. 0.0D0) hmxi = 1.0D0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0D0) go to 616
      seth = rwork(8)
      if (seth .lt. 0.0D0) go to 609
! Check RTOL and ATOL for legality. ------------------------------------
 60   rtoli = rtol(1)
      atoli = atol(1)
      do 65 i = 1, n
          if (itol .ge. 3) rtoli = rtol(i)
          if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
          if (rtoli .lt. 0.0D0) go to 619
          if (atoli .lt. 0.0D0) go to 620
 65       continue
!-----------------------------------------------------------------------
! Compute required work array lengths, as far as possible, and test
! these against LRW and LIW.  Then set tentative pointers for work
! arrays.  Pointers to RWORK/IWORK segments are named by prefixing L to
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are denoted  WM, YH, SAVF, EWT, ACOR.
! If MITER = 1 or 2, the required length of the matrix work space WM
! is not yet known, and so a crude minimum value is used for the
! initial tests of LRW and LIW, and YH is temporarily stored as far
! to the right in RWORK as possible, to leave the maximum amount
! of space for WM for matrix preprocessing.  Thus if MITER = 1 or 2
! and MOSS .ne. 2, some of the segments of RWORK are temporarily
! omitted, as they are not needed in the preprocessing.  These
! omitted segments are: ACOR if ISTATE = 1, EWT and ACOR if ISTATE = 3
! and MOSS = 1, and SAVF, EWT, and ACOR if ISTATE = 3 and MOSS = 0.
!-----------------------------------------------------------------------
      lrat = lenrat
      if (istate .eq. 1) nyh = n
      lwmin = 0
      if (miter .eq. 1) lwmin = 4*n + 10*n/lrat
      if (miter .eq. 2) lwmin = 4*n + 11*n/lrat
      if (miter .eq. 3) lwmin = n + 2
      lenyh = (maxord+1)*nyh
      lrest = lenyh + 3*n
      lenrw = 20 + lwmin + lrest
      iwork(17) = lenrw
      leniw = 30
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)
     1   leniw = leniw + n + 1
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
      lia = 31
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)
     1   leniw = leniw + iwork(lia+n) - 1
      iwork(18) = leniw
      if (leniw .gt. liw) go to 618
      
      return
      end