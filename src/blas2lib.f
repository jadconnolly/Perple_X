c this file contains:

c 1) BLAS level 2 fortran subroutines. modified subroutines are named by 
c    appending a digit to the original name or in the case of 6 char
c    names by replacing the final character by a digit (1).
c    BLAS is a freely-available software package at www.netlib.org/blas/

c 2) lpsol a fortran subroutine, and any non-BLAS subroutines it invokes,
c    to solve linear programming problems; and nlpsol a fortran subroutine,
c    and any non-BLAS subroutines it invokes, to solve non-linear programming 
c    problems by succesive quadratic programming. Both after

c    Gill P E, Murray W, Saunders M A and Wright M H (1984) 
c    Procedures for optimization problems with a mixture of bounds and 
c    general linear constraints ACM Trans. Math. Software 10 282�298

c    Gill PE, Hammarling S, Murray W, Saunders MA and Wright MH (1986) 
c    User�s Guide for LSSOL (version 1.0) Report SOL 86�1 
c    Department of Operations Research, Stanford University.

      subroutine lpsol (n,nclin,a,lda,bl,bu,cvec,istate,x,iter,obj,ax,
     *                  clamda,iw,leniw,w,lenw,ifail,istart,tol,lpprob)
c----------------------------------------------------------------------
c     lpsol solves (may also be called as the fp problem)

c           minimize               c' x

c                                 ( x)
c           subject to    bl  .le.(  ).ge.  bu,
c                                 (ax)

c     where a is a constant nclin by n matrix.
c     the feasible region is defined by a mixture of linear equality or
c     inequality constraints on x.

c     n - the number of variables (dimension of x).
c     nclin - the number of general linear constraints (rows of a).
c----------------------------------------------------------------------
      implicit none

      character prbtyp*2, msg*6

      logical cset, done, found, halted, rowerr, rset, unitq, vertex

      integer istart, lpprob, ifail, iter, lda, leniw,
     *        lenw, n, nclin, istate(n+nclin), iw(leniw),
     *        ianrmj, inform, it, itmax, j, jinf, jmax,
     *        ldh, ldr, litotl, lwtotl, minact, minfxd,
     *        nact1, nactiv, nartif, ncnln, nctotl, 
     *        nfree, ngq, nmoved, nrejtd, nrz, numinf, nviol, nz

      double precision a(lda,*), ax(*), bl(n+nclin), bu(n+nclin), tol, 
     *                 clamda(n+nclin), cvec(*), w(lenw), x(n), obj, 
     *                 amin, condmx, errmax, feamax, feamin,
     *                 xnorm, dnrm2

      external dnrm2

      integer lkactv, lkx, lfeatu, lanorm, lad, ld, lgq, lcq, lrlam,
     *        lr, lt, lq, lwtinf, lwrk
      common/ cstlcl /lkactv,lkx,lfeatu,lanorm,lad,ld,lgq,lcq,lrlam,
     *                lr, lt, lq, lwtinf, lwrk

      double precision wmach
      common/ cstmch /wmach(10)

      integer ldt, ldq, ncolt
      common/ ngg004 /ldt, ncolt, ldq

      integer itnfix, kdegen, ndegen, nfix
      double precision tolinc, tolx0
      common/ ngg005 /tolx0, tolinc, kdegen, ndegen, itnfix, nfix(2)

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      logical prnt
      integer isdel, jadd, jdel
      double precision alfa, trulam
      common/ ngg007 /alfa, trulam, isdel, jdel, jadd, prnt

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin

      integer itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *        maxact, mxfree, maxnz
      common/ ngg010 /itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *                maxact, mxfree, maxnz

      double precision bigbnd, bigdx, bndlow, bndupp, tolact, tolfea, 
     *                 tolrnk
      common/ ngg011 /bigbnd, bigdx, bndlow, bndupp,
     *                tolact, tolfea, tolrnk
c----------------------------------------------------------------------
c                                 istart - 0 - cold start, 1 - warm start, 2 - hot
      lcrash = istart
c                                 cold start:  x may be provided, but 
c                                 in perplex it isn't and therefore 
c                                 x is initialized in cmcrsh.

c                                 warm start:  initial working set is in istate.

c                                 hot  start:  work arrays iw and w have 
c                                 been the first three elements of iw contain details
c                                 dimensions of the initial working set.

c                                 hot start is only used for static composition 
c                                 optimizations.

c                                 tol - feasibility tolerance
      tolfea = tol
c                                 problem type 1 - fp, 2 - lp
      lprob = lpprob
c                                 f(n) parameters
      nctotl = n + nclin
      itmax1 = max(50,5*(nctotl))
      maxact = min(n,nclin)
      maxnz = n
      mxfree = n

      if (nclin.lt.n) then
         mxfree = nclin + 1
         maxnz = mxfree
      end if

      inform = 0
      iter = 0
      condmx = 1d0/epspt5
c                                 if the matrix of general constraints
c                                 has fewer rows than columns (nclin < n),
c                                 a non-zero value is known for minfxd.
c                                 signal with vertex
      vertex = nclin.lt.n

      minfxd = n - mxfree
      minact = mxfree - maxnz

      ldt = max(maxnz,maxact)
      ncolt = mxfree

      ldr = ldt
      ncnln = 0
      ldh = 1
c                                 integer workspace pointers
      lkactv = 4
      lkx = 4 + n
c                                 real workspace pointers
      lfeatu = 1
c                                 for lpcore
      lanorm = 1 + nctotl
      lad = lanorm + nclin
      ld = lad + nclin
      lgq = ld + n
      lcq = lgq + n

      if (lprob.eq.1) then
         prbtyp = 'fp'
         cset = .false.
         lrlam = lcq
      else 
         prbtyp = 'lp'
         cset = .true.
         lrlam = lcq + n
      end if

      lr = lrlam + n
      lt = lr
      lq = lt + ldt*ncolt

      if (nclin.eq.0) then
         ldq = 1
         lwtinf = lq
      else
         ldq = max(1,mxfree)
         lwtinf = lq + ldq*ldq
      end if

      lwrk = lwtinf + nctotl

      litotl = 3 + 2*n
      lwtotl = lwrk + nctotl - 1

      nmoved = 0
c                                 nfix(j) counter of the number of times 
c                                 variables have been placed on the working set,
c                                 j=1 if infeasible, j=2 if feasible
      nfix = 0
c                                 ndegen: degenerate steps counter (incremented
c                                 by cmchzr)
      ndegen = 0
c                                 last iteration at which x was put on a constraint
      itnfix = 0

      if (istart.lt.2) then
c                                 cold or warm start. just about 
c                                 everything must be initialized.
c                                 the exception is istate on warm.
         clamda(1:nctotl) = tol/2d0
         w(1:nctotl) = tol

         ianrmj = lanorm

         do j = 1, nclin
            w(ianrmj) = dnrm2 (n,a(j,1),lda)
            ianrmj = ianrmj + 1
         end do

         if (nclin.gt.0) call scond (nclin,w(lanorm),1,asize,amin)

         call scond (nctotl,w(lfeatu),1,feamax,feamin)

         w(lwtinf:lwtinf+nctotl-1) = w(lfeatu:lfeatu+nctotl-1)/feamin

         call cmcrsh (istart,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
     *               lda,istate,iw(lkactv),iw(lkx),bigbnd,tolact,a,ax,
     *               bl,bu,clamda,x,w(lgq),w(lwrk))
c                                 get tq factorization of working set matrix
         unitq = .true.
         nz = nfree

         if (nactiv.gt.0) then
            it = nactiv + 1
            nact1 = nactiv
            nactiv = 0
            ngq = 0

            call rzadds (unitq,vertex,nact1,it,nactiv,nartif,nz,nfree,
     *                   nrejtd,ngq,n,ldq,lda,ldt,istate,iw(lkactv),
     *                   iw(lkx),condmx,a,w(lt),w(lgq),w(lq),w(lwrk),
     *                   w(ld),w(lrlam))
         end if

      else
c                                 hot can use previous clamda if
c                                 conditions are close, otherwise
c                                 better to steer cleer and reset
c                                 everything (but istate).
         clamda(1:nctotl) = tol/2d0
         x(1:n) = 0d0
         w(lgq:lgq+n-1) = 0d0
c                                  these counters reflect the values 
c                                  in istate, ergo better not to reset
c                                  istate.
         unitq = iw(1).eq.1
         nfree = iw(2)
         nactiv = iw(3)
         nz = nfree - nactiv

      end if

      if (cset) then
c                                 copy transformed linear term to cq.
         w(lcq:lcq+n-1) = cvec(1:n)

         call cmqmul (6,n,nz,nfree,ldq,unitq,iw(lkx),w(lcq),w(lq),
     *               w(lwrk))
      end if

      rset = .false.
      itmax = itmax2
      jinf = 0
c                                 minimizing the sum of infeasibilities:
c                                 nrz = nz steepest-descent in the 2-norm.
c                                 nrz = 0 steepest-descent in the infinity norm.
c                                 2-norm may be marginally faster?
      nrz = 0

      do
c                                 move x onto the constraints in the working set.
          call cmsetx (rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,ldt,
     *            istate,iw(lkactv),iw(lkx),jmax,errmax,xnorm,a,ax,bl,
     *            bu,w(lfeatu),w(lt),x,w(lq),w(ld),w(lwrk))

         if (rowerr) then
            msg = 'infeas'
            numinf = 1
            obj = errmax
            exit
         end if

         call lpcore (prbtyp,msg,cset,rset,unitq,iter,itmax,
     *             jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,nz,istate,
     *             iw(lkactv),iw(lkx),obj,numinf,xnorm,a,ax,bl,bu,
     *             cvec,clamda,w(lfeatu),x,w)

         found = msg.eq.'feasbl' .or. msg.eq.'optiml' .or. msg .eq.
     *        'weak  ' .or. msg.eq.'unbndd' .or. msg.eq.'infeas'
         halted = msg.eq.'itnlim'

         if (found) call cmdgen ('optimal',n,nclin,nmoved,iter,numinf,
     *                        istate,bl,bu,clamda,w(lfeatu),x)

         done = found .and. nviol.eq.0 .and. nmoved.eq.0

         if ( done .or. halted ) exit

      end do

      if (.not.rowerr) then 
c                                 set clamda for hot start
c                                 and or yclos routines 
         call cmprnt (nfree,n,nctotl,nactiv,iw(lkactv),iw(lkx),
     *             clamda,w(lrlam))
c                                 also set for hot start
         iw(1) = 0
         if (unitq) iw(1) = 1
         iw(2) = nfree
         iw(3) = nactiv

      end if 

      if (msg.eq.'optiml') then
         inform = 0
      else if (msg.eq.'feasbl') then
         inform = 0
      else if (msg.eq.'weak  ') then
         inform = 1
      else if (msg.eq.'unbndd') then
         inform = 2
      else if (msg.eq.'infeas') then
         inform = 3
      else if (msg.eq.'itnlim') then
         inform = 4
      else if (msg.eq.'errors') then
         inform = 6
      else if (msg.eq.'noprob') then
         inform = 7
      end if

      ifail = inform

      if (ifail.lt.3) ifail = 0

      if (ifail.lt.4) then
c                                 signal warm start possible, the calling
c                                 routine decides what to do.
         istart = 1
      else
         istart = 0
      end if
c                                 end of lpsol
      end

      subroutine nlpsol (n,nclin,lda,ldr,a,bl,bu,objfun,iter,istate,
     *                   clamda,objf,gradu,r,x,iw,leniw,w,lenw,idead)
c----------------------------------------------------------------------
c     nlpsol minimizes f(x)
c                                    (   x )
c            subject to    bl  .le.  ( a*x )  .le.  bu

c     by sequential quadratic programming algorithm, where f(x) is a smooth function,  
c     a is a constant matrix
c     n - the number of variables (dimension of  x)
c     nclin - the number of linear constraints (rows of the matrix  a)
c     ncnln - the number of nonlinear constraints (dimension of  c(x))
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical linobj, rowerr, bad

      integer iter, lda, ldr, leniw, lenw, n, nclin, idead, i, ianrmj,
     *        ikx, inform, maxnz, minact, itmxsv, itns, j, jinf, jmax, 
     *        lax, lclam, ldaqp, litotl, lwtotl, maxact, minfxd, mxfree,
     *        nact1, nctotl, ngq, nres, numinf, nlperr, 
     *        nmajor, nminor, nrank, nrejtd, nz1, istate(n+nclin), 
     *        iw(leniw), lent, lenzy, lfeatl, lgrad, lkx

      double precision a(lda,*), bl(n+nclin), bu(n+nclin), g0, x(n), 
     *                 clamda(n+nclin),  gradu(n), r(ldr,*), w(lenw), 
     *                 amin, condmx, ctx, errmax, objf, fdnorm, feamax,
     *                 feamin, obj, rootn, ssq1, suminf, xnorm, dnrm2

      external objfun, dnrm2

      integer lkactv, lanorm, lap, lqpdx, lres, lres0, lqphz,
     *        lqpgq, lgq, lrlam, lt, lq, lwtinf, lwrk1, lqptol
      common/ cstlnp /lkactv, lanorm, lap, lqpdx, lres, lres0, lqphz,
     *                lqpgq, lgq, lrlam, lt, lq, lwtinf, lwrk1, lqptol

      integer liperm , ladx, lbl, lbu, ldx, lgq1, lx1, lwrk2, lhctrl
      common/ cstln2 /liperm, ladx, lbl, lbu, ldx, lgq1, lx1,
     *                lwrk2, lhctrl

      double precision wmach
      common/ cstmch /wmach(10)

      integer ldt, ldq, ncolt
      common/ ngg004 /ldt, ncolt, ldq

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin

      integer itmax1, itmax2, lprob
      common/ ngg016 /itmax1, itmax2, lprob

      double precision rcndbd, rfrobn, drmax, drmin
      common/ ngg018 /rcndbd, rfrobn, drmax, drmin

      logical unitq
      integer nactiv, nfree, nz
      common/ ngg001 /nactiv, nfree, nz, unitq

      double precision cdint, ctol, dxlim, epsrf, eta, fdint, ftol,
     *                 hcndbd
      common/ ngg021 /cdint, ctol, dxlim, epsrf, eta,
     *                fdint, ftol, hcndbd

      equivalence (itmxnp,nmajor), (itmax2,nminor)
c----------------------------------------------------------------------
c                                 f(n) parameters 
      hcndbd = max(1d0/(1d2*wmach(3)*dble(n)),1d6)
      nctotl = n + nclin

      nmajor = max(50,3*nctotl)
      nminor = nmajor
      itmax1 = nmajor

      rootn = sqrt(dble(n))

      idead = -1

      inform = 0
c                                 constraint feasibility tolerance
      tolfea = ctol

      minact = 0
      minfxd = 0

      mxfree = n - minfxd
      maxact = max(1,min(n,nclin))
      maxnz = n - (minfxd+minact)

      if (nclin.eq.0) then
         ldq = 1
         ldt = 1
         ncolt = 1
      else
         ldq = max(1,mxfree)
         ldt = max(maxnz,maxact)
         ncolt = mxfree
      end if

      ldaqp = max(nclin,1)
      if (nclin.gt.0) ldaqp = lda

c                                 array lengths that depend on problem dimension
      if (nclin.eq.0) then
         lent = 0
         lenzy = 0
      else
         lent = ldt*ncolt
         lenzy = ldq*ldq
      end if
c                                 w(1:2*n) were reserved for fwd and ctl difference
c                                 increments, these are now free (set in cxt009).
      lkactv = 1
      lkx = 1 + n
      liperm = 1 + 2*n

      lhctrl = 1 + n
      lanorm = lhctrl + n
      lqpgq = lanorm + nclin 
      lgq = lqpgq + n
      lrlam = lgq + n
      lt = lrlam + n
      lq = lt + lent
c                                 addresses used by npiqp
      lap = lq + lenzy
      lqpdx = lap + nclin 
      lres = lqpdx + n
      lres0 = lres + n
      lqphz = lres0 + n
      lwtinf = lqphz + n
      lwrk1 = lwtinf + nctotl
      lqptol = lwrk1 + nctotl
c                                  addresses used by npcore
      ladx = lqptol + nctotl
      lbl = ladx + nclin 
      lbu = lbl + nctotl
      ldx = lbu + nctotl
      lgq1 = ldx + n
      lfeatl = lgq1 + n
      lx1 = lfeatl + nctotl
      lwrk2 = lx1 + n
      lgrad = lwrk2 + nctotl

      litotl = liperm + nctotl - 1
      lwtotl = lgrad + n - 1

      lax = lwtotl + 1
      lwtotl = lwtotl + nclin
      lax = min(lax,lwtotl)

      tolrnk = 0d0
      rcndbd = sqrt(hcndbd)
c                                 load feasibility tolerances.
      w(lfeatl:lfeatl+nctotl-1) = tolfea

      if (nclin.gt.0) then

         ianrmj = lanorm

         do j = 1, nclin
            w(ianrmj) = dnrm2 (n,a(j,1),lda)
            ianrmj = ianrmj + 1
         end do

         call scond (nclin,w(lanorm),1,asize,amin)

      end if

      call scond (nctotl,w(lfeatl),1,feamax,feamin)

      w(lwtinf:lwtinf+nctotl-1) = w(lfeatl:lfeatl+nctotl-1)/feamin
c                                 input values of x and (optionally) istate 
c                                 are used by lscrsh to define the initial working set.
      call lscrsh (nclin,nctotl,nactiv,nfree,n,lda,
     *            istate,iw(lkactv),tolact,a,w(lax),bl,bu,x,
     *            w(lwrk1))

      nres = 0
      ngq = 0
      condmx = 1d0/epspt5

      unitq = .true.
      iter = 0

      ikx = lkx
      do i = 1, n
         iw(ikx) = i
         ikx = ikx + 1
      end do

      call smload ('upper-triangular',n,n,0d0,1d0,r,ldr)
      rfrobn = rootn
      nrank = 0

      rhonrm = 0d0
      rhodmp = 1d0
c                                 re-order kx so that the free variables come first
      call lsbnds (unitq,inform,nz,nfree,nrank,nres,ngq,n,ldq,lda,ldr,
     *             ldt,istate,iw(lkx),condmx,a,r,w(lt),w(lres0),
     *             w(lgq),w(lq),w(lwrk1),w(lwrk2),w(lrlam))
c                                 factorize the linear constraints in the initial working set.
      if (nactiv.gt.0) then

         nact1 = nactiv
         nactiv = 0

         call lsadds (unitq,inform,nact1,nactiv,nz,nfree,
     *                nrank,nrejtd,nres,ngq,n,ldq,lda,ldr,ldt,istate,
     *                iw(lkactv),iw(lkx),condmx,a,r,w(lt),w(lres0),
     *                w(lgq),w(lq),w(lwrk1),w(lwrk2),w(lrlam))

      end if
c                                 move x on to the linear constraints and
c                                 find a feasible point.
      ssq1 = 0d0
      linobj = .false.

      call lssetx (linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,n,
     *            nctotl,ldq,lda,ldr,ldt,istate,iw(lkactv),iw(lkx),
     *            jmax,errmax,ctx,xnorm,a,w(lax),bl,bu,w(lgq),w(lres),
     *            w(lres0),w(lfeatl),r,w(lt),x,w(lq),w(lwrk1),w(lwrk2))
c                                 call lscore to find a feasible x.
c                                 use work2 as the multiplier vector.
      jinf = 0
      lclam = lwrk2

      itmxsv = itmax1
      itmax1 = nminor

      call lscore ('fp problem',linobj,unitq,nlperr,itns,
     *            jinf,nclin,nctotl,nactiv,nfree,nrank,nz,nz1,n,lda,
     *            ldr,istate,iw(lkactv),iw(lkx),ctx,obj,ssq1,suminf,
     *            numinf,xnorm,bl,bu,a,w(lclam),w(lax),w(lfeatl),r,x,w)

      itmax1 = itmxsv

      if (nlperr.gt.0) return
c                                get finite difference increments
c                                at initial point, no need to save
c                                the point as already been output
c                                by resub
      outrpc = .false.
c                                compute objective function
      call objfun (n,x,objf,gradu,bad)

      if (bad) then 
         return
      end if 

      g0 = objf

      cntrl = .false.

      if (numric.and.lopt(66)) then
c                                 get finite difference increments:
         call chfd (n,fdnorm,objf,objfun,bl,bu,gradu,x,bad)
c                                 chfd will in fact return 2nd order 
c                                 derivatives, but setting cntrl at this
c                                 point is too costly
      else if (numric) then 
c                                 recompute gradient at first order
         call numder (objf,objfun,gradu,x,fdnorm,bl,bu,n,bad)

         fdincs = .false.

      end if

      if (bad) then 
         return
      end if

      w(lgrad:lgrad+n-1) = gradu(1:n)
      w(lgq:lgq+n-1) = gradu(1:n)
c                                 transform grad (w(lgq))
      call cmqmul (6,n,nz,nfree,ldq,unitq,iw(lkx),w(lgq),w(lq),w(lwrk1))
c                                 solve the problem.
      call npcore (unitq,inform,iter,n,nclin,nctotl,nactiv,nfree,nz,
     *             ldaqp,ldr,istate,iw(lkactv),iw(lkx),objf,
     *             fdnorm,xnorm,objfun,a,w(lax),bl,bu,clamda,
     *             w(lfeatl),w(lgrad),gradu,r,x,iw,w,lenw)

      if (g0-objf.gt.ftol.or.inform.eq.-1) then 
         idead = inform
      else
         idead = -2
      end if

      end

      subroutine lsmove (hitcon,hitlow,linobj,unitgz,nclin,nrank,nrz,n,
     *                  ldr,jadd,numinf,alfa,ctp,ctx,xnorm,ap,ax,bl,bu,
     *                  gq,hz,p,res,r,x,work)
c----------------------------------------------------------------------
c     lsmove  changes x to x + alfa*p and updates ctx, ax, res and gq
c     accordingly.
c     if a bound was added to the working set,  move x exactly on to it,
c     except when a negative step was taken (cmalf may have had to move
c     to some other closer constraint.)
c----------------------------------------------------------------------
      implicit none

      logical hitcon, hitlow, linobj, unitgz

      integer jadd, ldr, n, nclin, nrank, nrz, numinf

      double precision ap(*), ax(*), bl(*), bu(*), gq(*), hz(*), p(n),
     *                 r(ldr,*), res(*), work(*), x(n), bnd, dnrm2,
     *                 alfa, ctp, ctx, xnorm

      external dnrm2
c----------------------------------------------------------------------
      call daxpy (n,alfa,p,1,x,1)
      if (linobj) ctx = ctx + alfa*ctp

      if (hitcon .and. jadd.le.n) then
         bnd = bu(jadd)
         if (hitlow) bnd = bl(jadd)
         if (alfa.ge.0d0) x(jadd) = bnd
      end if
      xnorm = dnrm2 (n,x,1)

      if (nclin.gt.0) call daxpy (nclin,alfa,ap,1,ax,1)

      if (nrz.le.nrank) then
         if (unitgz) then
            res(nrz) = res(nrz) - alfa*hz(nrz)
         else
            call daxpy (nrz,(-alfa),hz,1,res,1)
         end if

         if (numinf.eq.0) then

c           update the transformed gradient gq so that
c           gq = gq + alfa*r'(hz).
c                            (0 )
            if (unitgz) then
               call daxpy (n-nrz+1,alfa*hz(nrz),r(nrz,nrz),ldr,
     *                     gq(nrz),1)
            else
               call dcopy (nrz,hz,1,work,1)
               call dtrmv ('u','t','n',nrz,r,ldr,work)
               if (nrz.lt.n) call dgemv ('t',nrz,n-nrz,1d0,r(1,nrz+1),
     *                                  ldr,hz,0d0,work(nrz+1))
               call daxpy (n,alfa,work,1,gq,1)

            end if
         end if
      end if
c                                 end of lsmove
      end

      subroutine cmalf1 (firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)
c----------------------------------------------------------------------
c     cmalf1  finds steps palfa1, palfa2 such that
c        x + palfa1*p  reaches a linear constraint that is currently not
c                      in the working set but is satisfied.
c        x + palfa2*p  reaches a linear constraint that is currently not
c                      in the working set but is violated.
c     the constraints are perturbed by an amount featol, so that palfa1
c     is slightly larger than it should be,  and palfa2 is slightly
c     smaller than it should be.  this gives some leeway later when the
c     exact steps are computed by cmalf.
c
c     constraints in the working set are ignored  (istate(j) .ge. 1).
c
c     if negstp is true, the search direction will be taken to be  - p.

c     values of istate(j)....

c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c
c     the values  -2  and  -1  do not occur once a feasible point has
c     been found.
c----------------------------------------------------------------------
      implicit none
 
      logical firstv, negstp, lastv

      integer jadd1, jadd2, n, nctotl, i, j, js, istate(nctotl)

      double precision absatp, atp, atx, res, rownrm, bigalf, bigbnd, 
     *                 palfa1, palfa2, pnorm, anorm(*), ap(*), ax(*), 
     *                 bl(nctotl), bu(nctotl), featol(nctotl), p(n),x(n)

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9
c----------------------------------------------------------------------
      lastv = .not. firstv
      jadd1 = 0
      jadd2 = 0
      palfa1 = bigalf

      palfa2 = 0d0
      if (firstv) palfa2 = bigalf

      do j = 1, nctotl
         js = istate(j)
         if (js.le.0) then
            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               rownrm = 1d0
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               rownrm = 1d0 + anorm(i)
            end if
            if (negstp) atp = -atp

            if (abs(atp).le.epspt9*rownrm*pnorm) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.

               res = -1d0
            else if (atp.le.0d0 .and. js.ne.-2) then

c              a'x  is decreasing and the lower bound is not violated.

c              first test for smaller palfa1.

               absatp = -atp
               if (bl(j).gt.(-bigbnd)) then
                  res = atx - bl(j) + featol(j)
                  if (bigalf*absatp.gt.abs(res)) then
                     if (palfa1*absatp.gt.res) then
                        palfa1 = res/absatp
                        jadd1 = j
                     end if
                  end if
               end if

               if (js.eq.-1) then

c                 the upper bound is violated.  test for either larger
c                 or smaller palfa2, depending on the value of firstv.

                  res = atx - bu(j) - featol(j)
                  if (bigalf*absatp.gt.abs(res)) then
                     if (firstv .and. palfa2*absatp.gt.res .or.
     *                   lastv .and. palfa2*absatp.lt.res) then
                        palfa2 = res/absatp
                        jadd2 = j
                     end if
                  end if
               end if
            else if (atp.gt.0d0 .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.
c              test for smaller palfa1.

               if (bu(j).lt.bigbnd) then
                  res = bu(j) - atx + featol(j)
                  if (bigalf*atp.gt.abs(res)) then
                     if (palfa1*atp.gt.res) then
                        palfa1 = res/atp
                        jadd1 = j
                     end if
                  end if
               end if

               if (js.eq.-2) then

c                 the lower bound is violated.  test for a new palfa2.

                  res = bl(j) - atx - featol(j)
                  if (bigalf*atp.gt.abs(res)) then
                     if (firstv .and. palfa2*atp.gt.res .or. lastv .and.
     *                   palfa2*atp.lt.res) then
                        palfa2 = res/atp
                        jadd2 = j
                     end if
                  end if
               end if
            end if

         end if
      end do
c                                 end of cmalf1
      end

      subroutine cmr1md (n,nrank,ldr,lenv,lenw,r,v,w,c,s)
c----------------------------------------------------------------------
c     cmr1md  modifies the  nrank*n  upper-triangular matrix  r  so that
c     q*(r + v*w')  is upper triangular,  where  q  is orthogonal,
c     v  and  w  are vectors, and the modified  r  overwrites the old.
c     q  is the product of two sweeps of plane rotations (not stored).
c     if required,  the rotations are applied to the nu columns of
c     the matrix  u.
c
c     the matrix v*w' is an (lenv) by (lenw) matrix.
c     the vector v is overwritten.
c----------------------------------------------------------------------
      implicit none

      integer j, ldr, lenv, lenw, n, nrank

      double precision c(n), r(ldr,*), s(n), v(n), w(n)
c----------------------------------------------------------------------
      j = min(lenv,nrank)
      if (nrank.gt.0) then

c        reduce  v to beta*e(j)  using a backward sweep of rotations
c        in planes (j-1, j), (j-2, j), ..., (1, j).

         call ssrotg('fixed','backwards',j-1,v(j),v,1,c,s)

c        apply the sequence of rotations to r. this generates a spike in
c        the j-th row of r, which is stored in s.

         call sutsrs('left',n,1,j,c,s,r,ldr)

c        form  beta*e(j)*w' + r.  this a spiked matrix, with a row
c        spike in row j.

         call daxpy (min(j-1,lenw),v(j),w,1,s,1)
         call daxpy (lenw-j+1,v(j),w(j),1,r(j,j),ldr)

c        eliminate the spike using a forward sweep of rotations in
c        planes (1, j), (2, j), ..., (j-1, j).

         call susqr ('left',n,1,j,c,s,r,ldr)

      end if
c                                 end of cmr1md
      end

      subroutine lpcolr(nrz,ldr,r,rzz)
c----------------------------------------------------------------------
c     lpcolr  loads the last column of the  nrz x nrz  triangular factor
c     rz  with the multiple  rzz  of the  nrz-th unit vector.
c----------------------------------------------------------------------
      implicit none

      integer ldr, nrz
      double precision rzz, r(ldr,*)
c----------------------------------------------------------------------
      if (nrz.eq.0) return

      call sload (nrz-1,0d0,r(1,nrz),1)
      r(nrz,nrz) = rzz
c                                 end of lpcolr
      end

      subroutine npalf (inform,n,nclin,alfa,alfmin,alfmax,bigbnd,
     *                  dxnorm,anorm,adx,ax,bl,bu,dx,x)
c----------------------------------------------------------------------
c     npalf  finds a step alfa such that the point x + alfa*p reaches
c     one of the slacks or linear constraints.  the step alfa is the
c     maximum step that can be taken without violating one of the slacks
c     or linear constraints that is currently satisfied.
c----------------------------------------------------------------------
      implicit none

      integer i, j, inform, n, nclin

      double precision alfa, alfmax, alfmin, bigbnd, dxnorm,
     *                  adx(*), anorm(*), ax(*), bl(*), bu(*), 
     *                  dx(n), x(n), adxi, axi, res, rownrm

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9
c----------------------------------------------------------------------
      alfa = alfmax
      j = 1

c     +    while (j .le. n+nclin .and. alfa.gt.alfmin) do
   20 if (j.le.n+nclin .and. alfa.gt.alfmin) then

         if (j.le.n) then
            axi = x(j)
            adxi = dx(j)
            rownrm = 1d0
         else 

            i = j - n
            axi = ax(i)
            adxi = adx(i)
            rownrm = anorm(i) + 1d0

         end if

         res = -1d0
         if (adxi.le.-epspt9*rownrm*dxnorm) then

c           constraint decreasing.

            adxi = -adxi
            if (bl(j).gt.-bigbnd) res = axi - bl(j)
         else if (adxi.gt.epspt9*rownrm*dxnorm) then

c           constraint increasing.

            if (bu(j).lt.bigbnd) res = bu(j) - axi
         end if

         if (res.gt.0d0 .and. alfa*adxi.gt.res) alfa = res/adxi

         j = j + 1
         go to 20
c        +    end while
      end if

c     determine alfa, the bound on the step to be taken.

      alfa = max(alfa,alfmin)

      inform = 0
      if (alfa.ge.alfmax) inform = 1
c                                 end of npalf
      end

      subroutine cmtsol(mode,nrowt,n,t,y)
c----------------------------------------------------------------------
c     cmtsol  solves equations involving a reverse-triangular matrix  t
c     and a right-hand-side vector  y,  returning the solution in  y.
c----------------------------------------------------------------------
      implicit none

      integer j, jj, l, n1, mode, n, nrowt

      double precision t(nrowt,*), y(n), yj          
c----------------------------------------------------------------------
      n1 = n + 1
      if (mode.eq.1) then

c        mode = 1  ---  solve  t * y(new) = y(old).

         do j = 1, n
            jj = n1 - j
            yj = y(j)/t(j,jj)
            y(j) = yj
            l = jj - 1
            if (l.gt.0 .and. yj.ne.0d0) call daxpy (l,(-yj),t(j+1,jj),1,
     *          y(j+1),1)
         end do
      else

c        mode = 2  ---  solve  t' y(new) = y(old).

         do j = 1, n
            jj = n1 - j
            yj = y(j)/t(jj,j)
            y(j) = yj
            l = jj - 1
            if (l.gt.0 .and. yj.ne.0d0) call daxpy (l,(-yj),t(jj,j+1),
     *          nrowt,y(j+1),1)
         end do
      end if

c     reverse the solution vector.

      if (n.gt.1) then
         l = n/2
         do j = 1, l
            jj = n1 - j
            yj = y(j)
            y(j) = y(jj)
            y(jj) = yj
         end do
      end if
c                                 end of cmtsol
      end

      subroutine nggnbu(n,nu,nrank,ldr,i,j,r,u,c,s)
c----------------------------------------------------------------------
c     nggnbu  interchanges the  i-th  and  j-th  (i.lt.j)  columns of
c     an  nrank*n  upper-trapezoidal matrix  r   and restores the
c     resulting matrix to upper-trapezoidal form using two sweeps of
c     plane rotations applied on the left.  r is overwritten.

c     if nu.gt.0,  the rotations are applied to the  nu  columns of
c     the matrix  u.
c----------------------------------------------------------------------
      implicit none 

      integer i, j, ldr, n, nrank, nu, lenj

      double precision c(n), r(ldr,*), s(n), u(n,*)
c----------------------------------------------------------------------
c     swap the elements of the i-th and j-th columns of r on, or above,
c     the main diagonal.

      call dswap(min(i,nrank),r(1,i),r(1,j),1)
      lenj = min(j,nrank)

      if (lenj.gt.i) then

c        reduce elements  r(i+1,j), ..., r(lenj,j)  to  beta*e(lenj)
c        using a backward sweep in planes
c        (lenj-1,lenj), (lenj-2,lenj), ..., (i+1,lenj).
c        if required, apply the sequence of rotations to u.

         call ssrotg ('fixed','backwards',lenj-i-1,r(lenj,j),r(i+1,j),1,
     *               c(i+1),s(i+1))

         if (nu.gt.0) call sgesrc ('left','bottom','backwards',n,nu,i+1,
     *                            lenj,c,s,u,n)

c        put 0d0s into the j-th column of r in positions corresponding
c        to the sub-diagonals of the i-th column.

         s(i) = r(lenj,j)
         call sload (lenj-i,0d0,r(i+1,j),1)

c        apply the sequence of rotations to r.  this generates a spike
c        in the lenj-th row of r, which is stored in s.

         call sutsrs('left',n,i+1,lenj,c,s,r,ldr)

c        eliminate the spike using a forward sweep in planes
c        (i,lenj), (i+1,lenj), ..., (lenj-1,lenj).
c        if necessary, apply the sequence of rotations to u.

         call susqr ('left',n,i,lenj,c,s,r,ldr)
c
         if (nu.gt.0) call sgesrc('left','bottom','forwards',lenj,nu,i,
     *                            lenj,c,s,u,n)
      end if
c                                 end of nggnbu
      end

      subroutine lsdel (unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,res,r,t,gq,
     *                  zy,c,s)
c----------------------------------------------------------------------
c     lsdel  updates the least-squares factor r and the factorization
c     a(free) (z y) = (0 t) when a regular, temporary or artificial
c     constraint is deleted from the working set.
c----------------------------------------------------------------------
      implicit none

      logical unitq

      integer jdel, kdel, lda, ldr, ldt, ldzy, n, nactiv, nt,
     *        nfree, ngq, nrank, nres, nrz, nz, kactiv(n), kx(n),
     *        i, ir, itdel, jart, k, ka, ld, npiv, nrz1, nsup, idamax

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                 s(n), t(ldt,*), zy(ldzy,*), cs, sn

      external idamax

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin
c----------------------------------------------------------------------
      if (jdel.gt.0) then

c        regular constraint or temporary bound deleted.

         if (jdel.le.n) then

c           case 1.  a simple bound has been deleted.
c           columns nfree+1 and ir of r must be swapped.

            ir = nz + kdel

            itdel = 1
            nfree = nfree + 1

            if (nfree.lt.ir) then
               kx(ir) = kx(nfree)
               kx(nfree) = jdel
               if (nrank.gt.0) call nggnbu(n,nres,nrank,ldr,nfree,ir,r,
     *                                     res,c,s)
               call dswap(ngq,gq(nfree,1),gq(ir,1),n)
            end if

            if (.not. unitq) then

c              copy the incoming column of  a(free)  into the end of t.

               do ka = 1, nactiv
                  i = kactiv(ka)
                  t(ka,nfree) = a(i,jdel)
               end do

c              expand q by adding a unit row and column.

               if (nfree.gt.1) then
                  call sload (nfree-1,0d0,zy(nfree,1),ldzy)
                  call sload (nfree-1,0d0,zy(1,nfree),1)
               end if
               zy(nfree,nfree) = 1d0
            end if
         else

c           case 2.  a general constraint has been deleted.

            itdel = kdel
            nactiv = nactiv - 1

c           delete row  kdel  of t and move up the ones below it.
c           t becomes reverse lower hessenberg.

            do i = kdel, nactiv
               kactiv(i) = kactiv(i+1)
               ld = nfree - i
               call dcopy (i+1,t(i+1,ld),ldt,t(i,ld),ldt)
            end do
         end if

         nz = nz + 1

         if (nactiv.eq.0) then
            dtmax = 1d0
            dtmin = 1d0
         else

c           restore the nactiv x (nactiv+1) reverse-hessenberg matrix  t
c           to reverse-triangular form.  the last nactiv super-diagonal
c           elements are removed using a backward sweep of plane
c           rotations.  the rotation for the singleton in the first
c           column is generated separately.

            nsup = nactiv - itdel + 1

            if (nsup.gt.0) then
               npiv = nfree - itdel + 1
               if (nsup.gt.1) then
                  call dcopy (nsup-1,t(nactiv-1,nz+1),ldt-1,s(nz+1),1)
                  call nggqzz('remove',nactiv,1,nsup,c(nz+1),s(nz+1),
     *                        t(1,nz+1),ldt)
               end if

               call srotgc(t(nactiv,nz+1),t(nactiv,nz),cs,sn)
               t(nactiv,nz) = 0d0
               s(nz) = -sn
               c(nz) = cs

               call sgesrc('right','variable','backwards',nfree,nfree,
     *                     nz,npiv,c,s,zy,ldzy)
               call sgesrc('left ','variable','backwards',npiv,ngq,nz,
     *                     npiv,c,s,gq,n)

               nt = min(nrank,npiv)

               if (nt.lt.npiv .and. nt.gt.0) then

c                 r is upper trapezoidal, pretend r is (nt x n) and
c                 apply the rotations in columns  max(nt,nz)  thru npiv.

                  call sgesrc('right','variable','backwards',nt,n,
     *                        max(nt,nz),npiv,c,s,r,ldr)
               end if

c              apply the column transformations to the triangular part
c              of r.  a subdiagonal element is generated that must be
c              eliminated by a row rotation before the next column
c              transformation can be applied.

               if (nz.lt.nt) then
                  call sutsqr('right',nt,nz,nt,c,s,r,ldr)
               end if

c              apply the row rotations to the remaining rows of r.

               call sgesrc('left','variable','backwards',nt,n-nt,nz,nt,
     *                     c,s,r(1,min(nt+1,n)),ldr)

               if (nres.gt.0) call sgesrc('left','variable','backwards',
     *                                    nt,nres,nz,nt,c,s,res,n)

            end if
            call scond (nactiv,t(nactiv,nz+1),ldt-1,dtmax,dtmin)
         end if
      end if

      nrz1 = nrz + 1

      if (nz.gt.nrz) then
         if (jdel.gt.0) then
            jart = nrz1 - 1 + idamax(nz-nrz1+1,gq(nrz1,1))
         else
            jart = -jdel
         end if

         if (jart.gt.nrz1) then

c           swap columns nrz1 and jart of r.

            if (unitq) then
               k = kx(nrz1)
               kx(nrz1) = kx(jart)
               kx(jart) = k
            else
               call dswap(nfree,zy(1,nrz1),zy(1,jart),1)
            end if

            call dswap(ngq,gq(nrz1,1),gq(jart,1),n)
            if (nrank.gt.0) call nggnbu(n,nres,nrank,ldr,nrz1,jart,r,
     *                                  res,c,s)
         end if
      end if

      nrz = nrz1
c                                 end of lsdel
      end

      subroutine npsetx (unitq,nclin,nactiv,nfree,nz,n,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,bl,bu,rpq,rpq0,dx,gq,r,t,zy,work)
c----------------------------------------------------------------------
c     npsetx   defines a point which lies on the initial working set for
c     the qp subproblem.  this routine is similar to lssetx except
c     that advantage is taken of the fact that the initial estimate of
c     the solution of the least-squares subproblem is zero.
c----------------------------------------------------------------------
      implicit none

      logical unitq

      integer ldaqp, ldr, ldt, ldzy, n, nactiv, nctotl,
     *        nfree, nz, istate(nctotl), kactiv(n), kx(n),
     *        i, j, k, nfixed, nr, nclin

      double precision adx(*), aqp(ldaqp,*), bl(nctotl), bu(nctotl),
     *                 dx(n), gq(n), r(ldr,*), rpq(n), rpq0(n),
     *                 t(ldt,*), work(n), zy(ldzy,*), dxnorm, gdx, bnd,
     *                 ddot, dnrm2

      external ddot, dnrm2
c----------------------------------------------------------------------
      nfixed = n - nfree

      gdx = 0d0

      call sload (n,0d0,dx,1)
      call sload (n,0d0,rpq,1)
      call sload (n,0d0,rpq0,1)

      if (nactiv+nfixed.gt.0) then

c        set  work = residuals for constraints in the working set.
c        solve for  dx,  the smallest correction to  x  that gives a
c        point on the constraints in the working set.
c        set the fixed variables on their bounds,  solve the triangular
c        system  t*(dxy) = residuals,  and define  dx = y*(dxy).
c        use  (dxy)  to update  d(=pr)  as  d = d - r'( 0 ).
c                                                     (dxy)
c
         do i = 1, nfixed
            j = kx(nfree+i)
            if (istate(j).le.3) then
               bnd = bl(j)
               if (istate(j).eq.2) bnd = bu(j)
               dx(j) = bnd
               work(nfree+i) = bnd
            else
               work(nfree+i) = 0d0
            end if
         end do
c
         do i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(nz+i) = bnd - ddot (n,aqp(k,1),ldaqp,dx)
         end do

         if (nactiv.gt.0) call cmtsol (1,ldt,nactiv,t(1,nz+1),
     *                                work(nz+1))
         call dcopy (nactiv+nfixed,work(nz+1),1,dx(nz+1),1)
         if (nz.gt.0) call sload (nz,0d0,dx,1)

         gdx = ddot (nactiv+nfixed,gq(nz+1),1,dx(nz+1))

         if (nz.lt.n) then
            call dgemv ('n',nz,n-nz,-1d0,r(1,nz+1),ldr,dx(nz+1),1d0,
     *                 rpq)
            if (nz.lt.n) then
               nr = ldr
               if (nz+1.eq.n) nr = 1
               call dcopy (n-nz,dx(nz+1),1,rpq(nz+1),1)
               call dscal (n-nz,-1d0,rpq(nz+1),1)
               call dtrmv ('u','n','n',n-nz,r(nz+1,nz+1),nr,rpq(nz+1))

            end if
         end if

         call cmqmul (2,n,nz,nfree,ldzy,unitq,kx,dx,zy,work)

      end if

c     compute the 2-norm of  dx.
c     initialize  a*dx.

      dxnorm = dnrm2 (n,dx,1)
      if (nclin.gt.0) call dgemv ('n',nclin,n,1d0,aqp,ldaqp,dx,0d0,adx)
c                                 end of npsetx
      end

      subroutine lsmuls (n,nactiv,nfree,lda,ldt,numinf,nz,
     *                   nrz,istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,
     *                   jtiny,jbigst,kbigst,trulam,a,anorms,gq,rlamda,
     *                   t,wtinf)
c----------------------------------------------------------------------
c     lsmuls  first computes the lagrange multiplier estimates for the
c     given working set.  it then determines the values and indices of
c     certain significant multipliers.  in this process, the multipliers
c     for inequalities at their upper bounds are adjusted so that a
c     negative multiplier for an inequality constraint indicates non-
c     optimality.  all adjusted multipliers are scaled by the 2-norm
c     of the associated constraint row.  in the following, the term
c     minimum refers to the ordering of numbers on the real line,  and
c     not to their magnitude.

c     jsmlst  is the index of the minimum of the set of adjusted
c             multipliers with values less than  - dinky.  a negative
c             jsmlst defines the index in q'g of the artificial
c             constraint to be deleted.
c     ksmlst  marks the position of general constraint jsmlst in kactiv.

c     jbigst  is the index of the largest of the set of adjusted
c             multipliers with values greater than (1 + dinky).
c     kbigst  marks its position in kactiv.

c     on exit,  elements 1 thru nactiv of rlamda contain the unadjusted
c     multipliers for the general constraints.  elements nactiv onwards
c     of rlamda contain the unadjusted multipliers for the bounds.
c----------------------------------------------------------------------
      implicit none 

      double precision dinky, trulam
      integer jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *                  lda, ldt, n, nactiv, nfree, nrz, numinf,
     *                  nz

      double precision a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                  wtinf(*)
      integer istate(*), kactiv(n), kx(n)

      double precision anormj, biggst, blam, rlam, scdlam, smllst,
     *                  tinylm
      integer i, is, j, k, l, nfixed
c----------------------------------------------------------------------
      nfixed = n - nfree

      jsmlst = 0
      ksmlst = 0
      smllst = -dinky

      tinylm = dinky
      jtiny = 0

      jbigst = 0
      kbigst = 0
      biggst = 1d0 + dinky

      if (nrz.lt.nz) then

c        compute jsmlst for the artificial constraints.

         do j = nrz + 1, nz
            rlam = -abs(gq(j))
            if (rlam.lt.smllst) then
               smllst = rlam
               jsmlst = -j
            else if (rlam.lt.tinylm) then
               tinylm = rlam
               jtiny = j
            end if
         end do

      end if

c     compute jsmlst for regular constraints and temporary bounds.
c     first, compute the lagrange multipliers for the general
c     constraints in the working set, by solving  t'*lamda = y'g.

      if (n.gt.nz) call dcopy (n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call cmtsol (2,ldt,nactiv,t(1,nz+1),rlamda)
c     set elements nactiv, nactiv+1,... of  rlamda  equal to
c     the multipliers for the bound constraints.

      do l = 1, nfixed
         j = kx(nfree+l)
         blam = rlamda(nactiv+l)
         do k = 1, nactiv
            i = kactiv(k)
            blam = blam - a(i,j)*rlamda(k)
         end do
         rlamda(nactiv+l) = blam
      end do

c     find jsmlst and ksmlst.

      do k = 1, n - nz
         if (k.gt.nactiv) then
            j = kx(nz+k)
         else
            j = kactiv(k) + n
         end if

         is = istate(j)

         i = j - n
         if (j.le.n) anormj = 1d0
         if (j.gt.n) anormj = anorms(i)

         rlam = rlamda(k)

c        change the sign of the estimate if the constraint is in
c        the working set at its upper bound.

         if (is.eq.2) rlam = -rlam
         if (is.eq.3) rlam = abs(rlam)
         if (is.eq.4) rlam = -abs(rlam)

         if (is.ne.3) then
            scdlam = rlam*anormj
            if (scdlam.lt.smllst) then
               smllst = scdlam
               jsmlst = j
               ksmlst = k
            else if (scdlam.lt.tinylm) then
               tinylm = scdlam
               jtiny = j
            end if
         end if

         if (numinf.gt.0 .and. j.gt.jinf) then
            scdlam = rlam/wtinf(j)
            if (scdlam.gt.biggst) then
               biggst = scdlam
               trulam = rlamda(k)
               jbigst = j
               kbigst = k
            end if
         end if
      end do
c                                 end of lsmuls
      end

      subroutine cmfeas(n,nclin,istate,bigbnd,nviol,jmax,errmax,ax,bl,
     *                  bu,featol,x)
c----------------------------------------------------------------------
c     cmfeas  checks the residuals of the constraints that are believed
c     to be feasible.  the number of constraints violated by more than
c     featol is computed, along with the maximum constraint violation.
c----------------------------------------------------------------------
      implicit none

      double precision bigbnd, errmax
      integer jmax, n, nclin, nviol

      double precision ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), x(n)
      integer istate(n+nclin)

      double precision biglow, bigupp, con, feasj, res
      integer is, j
c----------------------------------------------------------------------
      biglow = -bigbnd
      bigupp = bigbnd

c     compute the number of constraints (nviol) violated by more than
c     featol and  the maximum constraint violation (errmax).
c     (the residual of a constraint in the working set is treated as if
c     it were an equality constraint fixed at that bound.)

      nviol = 0
      jmax = 0
      errmax = 0d0

      do j = 1, n + nclin
         is = istate(j)

         if (is.ge.0) then
            feasj = featol(j)

            if (j.le.n) then
               con = x(j)
            else
               con = ax(j-n)
            end if

c           check for constraint violations.

            if (bl(j).gt.biglow) then
               res = bl(j) - con
               if (res.gt.feasj) then
                  nviol = nviol + 1
                  go to 20
               end if
            end if

            if (bu(j).lt.bigupp) then
               res = bu(j) - con
               if (res.lt.(-feasj)) then
                  nviol = nviol + 1
                  res = -res
                  go to 20
               end if
            end if

c           this constraint is satisfied,  but count a large residual
c           as a violation if the constraint is in the working set.

            res = 0d0

            if (is.eq.1) then
               res = abs(bl(j)-con)

            else if (is.eq.2) then
               res = abs(bu(j)-con)

            else if (is.eq.3) then
               res = abs(bu(j)-con)
            end if

            if (res.gt.feasj) nviol = nviol + 1

   20       if (res.gt.errmax) then
               jmax = j
               errmax = res
            end if
         end if
      end do
c                                 end of cmfeas
      end

      subroutine srchc (first,done,imprvd,inform,maxf,numf,
     *                  alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *                  toltny,alfa,alfbst,fbest,gbest)
c----------------------------------------------------------------------
c     srchc finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     srchc   requires both  f(alpha)  and  f'(alpha) to be evaluated at
c     points in the interval.  estimates of the minimizer are computed
c     using safeguarded cubic interpolation.
c----------------------------------------------------------------------
      implicit none

      logical badfun, braktd, closef, crampd, extrap, fitok, done, 
     *        found, moved, quitf, quiti, setxw, wset, first, imprvd

      integer inform, maxf, numf, nsamea, nsameb

      double precision a, absr, artifa, artifb, b, daux, dtry, factor,
     *                 fw, gw, q, r, s, scale, tol, tolmax, truea,
     *                 trueb, xmidpt, xtry, xw, alfa, alfbst, alfmax, 
     *                 epsaf, fbest, ftry, g0, gbest, gtry, targtg, 
     *                 tolabs, tolrel, toltny

      save              braktd, crampd, extrap, moved, wset, nsamea,
     *                  nsameb, a, b, factor, xtry, xw, fw, gw, tolmax
c----------------------------------------------------------------------
      badfun = .false.
      quitf = .false.
      quiti = .false.
      imprvd = .false.

      if (first) then

c        first entry.  initialize various quantities, check input data
c        and prepare to evaluate the function at the initial alfa.

         first = .false.
         numf = 0
         alfbst = 0d0
         badfun = alfmax .le. toltny .or. g0 .ge. 0d0
         done = badfun
         moved = .false.

         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            wset = .false.
            nsamea = 0
            nsameb = 0

            tolmax = tolabs + tolrel*alfmax
            a = 0d0
            b = alfmax + tolmax
            factor = 5d0
            tol = tolabs
            xtry = alfa
         end if
      else

c        subsequent entries. the function has just been evaluated at
c        alfa = alfbst + xtry,  giving ftry and gtry.

         numf = numf + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b = alfmax - alfbst + tolmax
         end if

c        see if the new step is better.  if alfa is large enough that
c        ftry can be distinguished numerically from zero,  the function
c        is required to be sufficiently negative.

         closef = abs(ftry-fbest) .le. epsaf
         if (closef) then
            imprvd = abs(gtry) .le. abs(gbest)
         else
            imprvd = ftry.lt.fbest
         end if

         if (imprvd) then

c           an improvement.  the new point becomes the
c           origin and other points are shifted accordingly.

            fw = fbest
            fbest = ftry
            gw = gbest
            gbest = gtry
            alfbst = alfa
            moved = .true.
c
            a = a - xtry
            b = b - xtry
            xw = 0d0 - xtry
            wset = .true.
            extrap = xw.lt.0d0 .and. gbest.lt.0d0 .or. xw .gt.
     *               0d0 .and. gbest.gt.0d0

c           decrease the length of the interval of uncertainty.

            if (gtry.le.0d0) then
               a = 0d0
               nsamea = 0
            else
               b = 0d0
               nsameb = 0
               braktd = .true.
            end if
         else

c           the new function value is not better than the best point so
c           far.  the origin remains unchanged but the new point may
c           qualify as xw.  xtry must be a new bound on the best point.

            if (xtry.le.0d0) then
               a = xtry
               nsamea = 0
            else
               b = xtry
               nsameb = 0
               braktd = .true.
            end if

c           if xw has not been set or ftry is better than fw, update the
c           points accordingly.

            if (wset) then
               setxw = ftry.lt.fw .or. .not. extrap
            else
               setxw = .true.
            end if

            if (setxw) then
               xw = xtry
               fw = ftry
               gw = gtry
               wset = .true.
               extrap = .false.
            end if
         end if

c        check the termination criteria.  wset will always be .true..

         tol = tolabs + tolrel*alfbst
         truea = alfbst + a
         trueb = alfbst + b

         found = abs(gbest) .le. targtg
         quitf = numf .ge. maxf
         quiti = b - a .le. tol + tol

         if (quiti .and. .not. moved) then

c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.

            tol = tol/1d1
            tolabs = tol
            quiti = abs(fw) .le. epsaf .or. tol .le. toltny
         end if

         done = quitf .or. quiti .or. found

c        proceed with the computation of a trial steplength.
c        the choices are...
c        1. parabolic fit using derivatives only, if the f values are
c           close.
c        2. cubic fit for a minimizer, using both f and f'.
c        3. damped cubic or parabolic fit if the regular fit appears to
c           be consistently overestimating the distance to a minimizer.
c        4. bisection, geometric bisection, or a step of  tol  if
c           choices 2 or 3 are unsatisfactory.

         if (.not. done) then
            xmidpt = (a+b)/2d0
            s = 0d0
            q = 0d0

            if (closef) then

c              fit a parabola to the two best gradient values.

               s = gbest
               q = gbest - gw

            else

c              fit cubic through  fbest  and  fw.


               fitok = .true.
               r = 3d0 * (fbest-fw)/xw + gbest + gw
               absr = abs(r)
               s = sqrt(abs(gbest))*sqrt(abs(gw))

c              compute  q =  the square root of  r*r - gbest*gw.
c              the method avoids unnecessary underflow and overflow.

               if ((gw.lt.0d0 .and. gbest.gt.0d0)
     *             .or. (gw.gt.0d0 .and. gbest.lt.0d0)) then
                  scale = absr + s
                  if (scale.eq.0d0) then
                     q = 0d0
                  else
                     q = scale*sqrt((absr/scale)**2+(s/scale)**2)
                  end if
               else if (absr.ge.s) then
                  q = sqrt(absr+s)*sqrt(absr-s)
               else
                  fitok = .false.
               end if

               if (fitok) then

c                 compute a minimizer of the fitted cubic.

                  if (xw.lt.0d0) q = -q
                  s = gbest - r - q
                  q = gbest - gw - q - q
               end if
            end if

c           construct an artificial interval  (artifa, artifb)  in which
c           the new estimate of a minimizer must lie.  set a default
c           value of xtry that will be used if the polynomial fit fails.

            artifa = a
            artifb = b
            if (.not. braktd) then

c              a minimizer has not been bracketed.  set an artificial
c              upper bound by expanding the interval  xw  by a suitable
c              factor.

               xtry = -factor*xw
               artifb = xtry
               if (alfbst+xtry.lt.alfmax) factor = 5d0 * factor

            else if (extrap) then

c              the points are configured for an extrapolation.
c              set a default value of  xtry  in the interval  (a, b)
c              that will be used if the polynomial fit is rejected.  in
c              the following,  dtry  and  daux  denote the lengths of
c              the intervals  (a, b)  and  (0, xw)  (or  (xw, 0),  if
c              appropriate).  the value of  xtry is the point at which
c              the exponents of  dtry  and  daux  are approximately
c              bisected.

               daux = abs(xw)
               dtry = b - a
               if (daux.ge.dtry) then
                  xtry = 5d0 * dtry * (0.1d0 + dtry/daux)/11d0
               else
                  xtry = sqrt(daux)*sqrt(dtry)/2d0
               end if
               if (xw.gt.0d0) xtry = -xtry

c              reset the artificial bounds.  if the point computed by
c              extrapolation is rejected,  xtry will remain at the
c              relevant artificial bound.

               if (xtry.le.0d0) artifa = xtry
               if (xtry.gt.0d0) artifb = xtry
            else

c              the points are configured for an interpolation.  the
c              default value xtry bisects the interval of uncertainty.
c              the artificial interval is just (a, b).

               xtry = xmidpt

               if (nsamea.ge.3 .or. nsameb.ge.3) then

c                 if the interpolation appears to be overestimating the
c                 distance to a minimizer,  damp the interpolation.

                  factor = factor/5d0
                  s = factor*s
               else
                  factor = 1d0
               end if
            end if

c           the polynomial fits give  (s/q)*xw  as the new step.
c           reject this step if it lies outside  (artifa, artifb).

            if (q.ne.0d0) then
               if (q.lt.0d0) s = -s
               if (q.lt.0d0) q = -q
               if (s*xw.ge.q*artifa .and. s*xw.le.q*artifb) then

c                 accept the polynomial fit.

                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = 0d0
                  end if

               end if
            end if
         end if
      end if

      if (.not. done) then
         alfa = alfbst + xtry
         if (braktd .or. alfa.lt.alfmax-tolmax) then

c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)
c
            if (xtry.le.a+tol .or. xtry.ge.b-tol) then
               if (a+b .le. 0d0) then
                  xtry = -tol
               else
                  xtry = tol
               end if
               alfa = alfbst + xtry
            end if
         else

c           the step is close to, or larger than alfmax, replace it by
c           alfmax to force evaluation of  f  at the boundary.

            braktd = .true.
            xtry = alfmax - alfbst
            alfa = alfmax
         end if
      end if

c     exit.

      if (done) then
         if (badfun) then
            inform = 8
         else if (found) then
            if (alfbst.lt.alfmax) then
               inform = 1
            else
               inform = 2
            end if
         else if (moved) then
            inform = 3
         else if (quitf) then
            inform = 7
         else if (crampd) then
            inform = 4
         else
            inform = 6
         end if
      end if
c                                 end of srchc
      end

      subroutine lsadd (unitq,inform,ifix,iadd,jadd,nactiv,nz,nfree,
     *                  nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,a,r,
     *                  t,res,gq,zy,w,c,s)
c----------------------------------------------------------------------
c     lsadd updates the factorization,  a(free) * (z y) = (0 t), when
c     a constraint is added to the working set. if nrank.gt.0, the
c     factorization  (r) = pcq  is also updated, where c is the
c                    (0)
c     least squares matrix, r is upper-triangular, and p is an
c     orthogonal matrix. the matrices c and p are not stored.

c     there are three separate cases to consider (although each case
c     shares code with another)...

c     (1) a free variable becomes fixed on one of its bounds when there
c         are already some general constraints in the working set.

c     (2) a free variable becomes fixed on one of its bounds when there
c         are only bound constraints in the working set.

c     (3) a general constraint (corresponding to row iadd of a) is
c         added to the working set.

c     in cases (1) and (2), we assume that  kx(ifix) = jadd.
c     in all cases, jadd is the index of the constraint being added.

c     if there are no general constraints in the working set,  the
c     matrix  q = (z y) is the identity and will not be touched.

c     if  nres.gt.0,  the row transformations are applied to the rows
c     of the  (n by nres) matrix  res.
c     if ngq.gt.0, the column transformations are applied to the
c     columns of the (ngq by n) matrix  gq'.
c----------------------------------------------------------------------
      implicit none

      logical bound, overfl, unitq

      integer iadd, ifix, inform, jadd, lda, ldr, ldt, ldzy, n, kx(n), 
     *        i, nactiv, nfree, ngq, nrank, nres, nz, nanew, nfmin, 
     *        npiv, nt

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                 s(n), t(ldt,*), w(n), zy(ldzy,*), condmx, 
     *                 cond, condbd, dtnew, tdtmax, tdtmin, dnrm2, sdiv

      external dnrm2, sdiv 

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin
c----------------------------------------------------------------------
c     if the condition estimator of the updated factors is greater than
c     condbd,  a warning message is printed.

      condbd = 1d0/epspt9

      overfl = .false.
      bound = jadd .le. n

      if (bound) then

c        a simple bound has entered the working set.  iadd  is not used.

         nanew = nactiv

         if (unitq) then

c           q  is not stored, but kx defines an ordering of the columns
c           of the identity matrix that implicitly define  q.
c           define the sequence of pairwise interchanges p that moves
c           the newly-fixed variable to position nfree.
c           reorder kx accordingly.

            do i = 1, nfree - 1
               if (i.ge.ifix) then
                  w(i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if
            end do

         else

c           q  is stored explicitly.

c           set  w = the  (ifix)-th  row of  q.
c           move the  (nfree)-th  row of  q  to position  ifix.

            call dcopy (nfree,zy(ifix,1),ldzy,w,1)
            if (ifix.lt.nfree) then
               call dcopy (nfree,zy(nfree,1),ldzy,zy(ifix,1),ldzy)
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else

c        a general constraint has entered the working set.
c        ifix  is not used.

         nanew = nactiv + 1

c        transform the incoming row of  a  by  q'.  use c as workspace.

         call dcopy (n,a(iadd,1),lda,w,1)
         call cmqmul(8,n,nz,nfree,ldzy,unitq,kx,w,zy,c)

c        check that the incoming row is not dependent upon those
c        already in the working set.

         dtnew = dnrm2 (nz,w,1)
         if (nactiv.eq.0) then

c           this is the only general constraint in the working set.

            cond = sdiv (asize,dtnew,overfl)
            tdtmax = dtnew
            tdtmin = dtnew
         else

c           there are already some general constraints in the working
c           set. update the estimate of the condition number.

            tdtmax = max(dtnew,dtmax)
            tdtmin = min(dtnew,dtmin)
            cond = sdiv (tdtmax,tdtmin,overfl)
         end if
c
         if (cond.gt.condmx .or. overfl) go to 60
c
         if (unitq) then
c
c           first general constraint added.  set  q = i.
c
            call smload ('general',nfree,nfree,0d0,1d0,zy,ldzy)
            unitq = .false.
         end if
      end if
c
      if (bound) then
         npiv = nfree
      else
         npiv = nz
      end if
c
      nt = min(nrank,npiv)
c
      if (unitq) then

c        q (i.e., zy) is not stored explicitly.
c        apply the sequence of pairwise interchanges p that moves the
c        newly-fixed variable to position nfree.

         if (ngq.gt.0) call sgeapr('left','transpose',nfree-1,w,ngq,gq,
     *                             n)
c
         if (nrank.gt.0) then
c
c           apply the pairwise interchanges to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1), s(2), ..., s(nt-1).
c
            call sutsr1 ('right',n,ifix,nt,s,r,ldr)
c
            if (nt.lt.npiv) then
c
c              r is upper trapezoidal.  apply the interchanges in
c              columns  nt  thru  npiv.
c
               do i = ifix, nt - 1
                  w(i) = i
               end do
c
               call sgeapr('right','normal',nfree-1,w,nt,r,ldr)
            end if
c
c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.
c
            call suhqr('left ',n,ifix,nt,c,s,r,ldr)
            if (nres.gt.0) call sgesrc('left','variable','forwards',nt,
     *                                 nres,ifix,nt,c,s,res,n)
         end if
      else

c        full matrix q.  define a sweep of plane rotations p such that
c                           pw = beta*e(npiv).
c        the rotations are applied in the planes (1,2), (2,3), ...,
c        (npiv-1,npiv).  the rotations must be applied to zy, r, t
c        and gq'.

         call ssrotg('varble','forwrds',npiv-1,w(npiv),w,1,c,s)
c
         if (bound .and. nactiv.gt.0) then
c
            call dcopy (nactiv,s(nz),1,w(nz),1)
c
            s(nz) = s(nz)*t(nactiv,nz+1)
            t(nactiv,nz+1) = c(nz)*t(nactiv,nz+1)
c
            call nggqzz('create',nactiv,1,nactiv,c(nz+1),s(nz+1),
     *                  t(1,nz+1),ldt)
            call dcopy (nactiv,s(nz),1,t(nactiv,nz),ldt-1)
c
            call dcopy (nactiv,w(nz),1,s(nz),1)
         end if
c
         if (ngq.gt.0) call sgesrc('left ','variable','forwards',npiv,
     *                             ngq,1,npiv,c,s,gq,n)
         call sgesrc('right','variable','forwards',nfree,nfree,1,npiv,c,
     *               s,zy,ldzy)
c
         if (nrank.gt.0) then
c
c           apply the rotations to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nt-1).
c
            nt = min(nrank,npiv)
            call sutsrh('right',n,1,nt,c,s,r,ldr)
c
            if (nt.lt.npiv) then
c
c              r is upper trapezoidal.  pretend r is (nt x n) and
c              apply the rotations in columns  nt  thru  npiv.

               call sgesrc('right','variable','forwards',nt,n,nt,npiv,c,
     *                     s,r,ldr)
            end if

c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.

            call suhqr('left ',n,1,nt,c,s,r,ldr)
            if (nres.gt.0) call sgesrc('left','variable','forwards',nt,
     *                                 nres,1,nt,c,s,res,n)
         end if

         if (bound) then

c           the last row and column of zy has been transformed to plus
c           or minus the unit vector e(nfree).  we can reconstitute the
c           columns of gq and r corresponding to the new fixed variable.

            if (w(nfree).lt.0d0) then
               nfmin = min(nrank,nfree)
               if (nfmin.gt.0) call dscal (nfmin,-1d0,r(1,nfree),1)
               if (ngq.gt.0) call dscal (ngq,-1d0,gq(nfree,1),n)
            end if

c           the diagonals of t have been altered.  recompute the
c           largest and smallest values.

            if (nactiv.gt.0) then
               call scond (nactiv,t(nactiv,nz),ldt-1,tdtmax,tdtmin)
               cond = sdiv (tdtmax,tdtmin,overfl)
            end if
         else

c           general constraint.  install the new row of t.

            call dcopy (nanew,w(nz),1,t(nanew,nz),ldt)
         end if
      end if

c     prepare to exit.  check the magnitude of the condition estimator.

   60 if (nanew.gt.0) then
         if (cond.lt.condmx .and. .not. overfl) then

c           the factorization has been successfully updated.

            inform = 0
            dtmax = tdtmax
            dtmin = tdtmin

         else
c           the proposed working set appears to be linearly dependent.

            inform = 1

         end if
      end if
c                                 end of lsadd
      end



      subroutine lsgetp(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,ap,res,
     *                  hz,p,gq,cq,r,zy,work)
c----------------------------------------------------------------------
c     lsgetp  computes the following quantities for  lscore.
c     (1) the vector  (hz1) = (rz1)(pz1).
c         if x is not yet feasible,  the product is computed directly.
c         if  rz1 is singular,  hz1  is zero.  otherwise  hz1  satisfies
c         the equations
c                        rz1'hz1 = -gz1,
c         where  g  is the total gradient.  if there is no linear term
c         in the objective,  hz1  is set to  dz1  directly.
c     (2) the search direction p (and its 2-norm).  the vector p is
c         defined as  z*(pz1), where  (pz1)  depends upon whether or
c         not x is feasible and the nonsingularity of  (rz1).
c         if  numinf.gt.0,  (pz1)  is the steepest-descent direction.
c         otherwise,  x  is the solution of the  nrz*nrz  triangular
c         system   (rz1)*(pz1) = (hz1).
c     (3) the vector ap,  where a is the matrix of linear constraints.
c----------------------------------------------------------------------
      implicit none

      logical linobj, singlr, unitgz, unitq

      integer lda, ldr, ldzy, n, nclin, nfree, nrank, nrz, kx(n),
     *        numinf

      double precision a(lda,*), ap(*), cq(*), gq(n), hz(*), p(n),
     *                 r(ldr,*), res(*), work(n), zy(ldzy,*), gtp, 
     *                 ddot, dnrm2, ctp, pnorm

      external ddot, dnrm2
c----------------------------------------------------------------------
      if (singlr) then

c        the triangular factor for the current objective function is
c        singular,  i.e., the objective is linear along the last column
c        of z1.  this can only occur when unitgz is true.

         if (nrz.gt.1) then
            call dcopy (nrz-1,r(1,nrz),1,p,1)
            call dtrsv ('u','n','n',nrz-1,r,ldr,p)
         end if
         p(nrz) = -1d0

         gtp = ddot (nrz,gq,1,p)
         if (gtp.gt.0d0) call dscal (nrz,-1d0,p,1)

         if (nrz.le.nrank) then
            if (numinf.eq.0) then
               if (unitgz) then
                  hz(nrz) = r(nrz,nrz)*p(nrz)
               else
                  call sload (nrz,0d0,hz,1)
               end if
            else
               hz(1) = r(1,1)*p(1)
            end if
         end if
      else

c        the objective is quadratic in the space spanned by z1.

         if (linobj) then
            if (unitgz) then
               if (nrz.gt.1) call sload (nrz-1,0d0,hz,1)
               hz(nrz) = -gq(nrz)/r(nrz,nrz)
            else
               call dcopy (nrz,gq,1,hz,1)
               call dscal (nrz,-1d0,hz,1)
               call dtrsv ('u','t','n',nrz,r,ldr,hz)
            end if
         else
            call dcopy (nrz,res,1,hz,1)
         end if

c        solve  rz1*pz1 = hz1.

         call dcopy (nrz,hz,1,p,1)
         call dtrsv ('u','n','n',nrz,r,ldr,p)

      end if

c     compute  p = z1*pz1  and its norm.

      if (linobj) ctp = ddot (nrz,cq,1,p)
      pnorm = dnrm2 (nrz,p,1)

      call cmqmul (1,n,nrz,nfree,ldzy,unitq,kx,p,zy,work)

c     compute  ap.

      if (nclin.gt.0) call dgemv ('n',nclin,n,1d0,a,lda,p,0d0,ap)

c                                 end of lsgetp
      end

      subroutine lsgset(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,res,featol,gq,cq,r,
     *                  x,wtinf,zy,wrk)
c----------------------------------------------------------------------
c     lsgset  finds the number and weighted sum of infeasibilities for
c     the bounds and linear constraints.   an appropriate transformed
c     gradient vector is returned in  gq.
c
c     positive values of  istate(j)  will not be altered.  these mean
c     the following...
c
c               1             2           3
c           a'x = bl      a'x = bu     bl = bu
c
c     other values of  istate(j)  will be reset as follows...
c           a'x lt bl     a'x gt bu     a'x free
c              - 2           - 1           0
c
c     if  x  is feasible,  lsgset computes the vector q(free)'g(free),
c     where  g  is the gradient of the the sum of squares plus the
c     linear term.  the matrix q is of the form
c                    (q(free)  0      ),
c                    (  0      i(fixed))
c     where  q(free)  is the orthogonal factor of  a(free)  and  a  is
c     the matrix of constraints in the working set.  the transformed
c     gradients are stored in gq.
c----------------------------------------------------------------------
      implicit none

      character prbtyp*2

      logical linobj, singlr, unitgz, unitq

      integer lda, ldr, ldzy, n, nclin, nfree, nrank, nrz, numinf, nz,
     *        j, k, istate(*), kx(n), isrank

      double precision a(lda,*), bl(*), bu(*), cq(*), featol(*), gq(n),
     *                 r(ldr,*), res(*), wrk(n), wtinf(*), x(n), dnrm2,
     *                 zy(ldzy,*), bigbnd, suminf, tolrnk, biglow, 
     *                 bigupp, ctx, feasj, rownrm, s, weight, ddot

      external ddot, dnrm2, isrank
c----------------------------------------------------------------------
      bigupp = bigbnd
      biglow = -bigbnd

      numinf = 0
      suminf = 0d0
      call sload (n,0d0,gq,1)

      do j = 1, n + nclin

         if (istate(j).le.0) then
            feasj = featol(j)
            if (j.le.n) then
               ctx = x(j)
            else
               k = j - n
               ctx = ddot (n,a(k,1),lda,x)
            end if
            istate(j) = 0

c           see if the lower bound is violated.

            if (bl(j).gt.biglow) then
               s = bl(j) - ctx
               if (s.gt.feasj) then
                  istate(j) = -2
                  weight = -wtinf(j)
                  go to 20
               end if
            end if

c           see if the upper bound is violated.

            if (bu(j).ge.bigupp) cycle
            s = ctx - bu(j)
            if (s.le.feasj) cycle
            istate(j) = -1
            weight = wtinf(j)

c           add the infeasibility.

   20       numinf = numinf + 1
            suminf = suminf + abs(weight)*s
            if (j.le.n) then
               gq(j) = weight
            else
               call daxpy (n,weight,a(k,1),lda,gq,1)
            end if
         end if

      end do

c     install  gq,  the transformed gradient.

      singlr = .false.
      unitgz = .true.

      if (numinf.gt.0) then
         call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,gq,zy,wrk)
         unitgz = .true.
      else if (numinf.eq.0 .and. prbtyp.eq.'fp') then
         call sload (n,0d0,gq,1)
      else
c                                 ready for the optimality phase.
c                                 set nrz so that rz1 is nonsingular.
         if (nrank.eq.0) then
            if (linobj) then
               call dcopy (n,cq,1,gq,1)
            else
               call sload (n,0d0,gq,1)
            end if
            nrz = 0
         else
c                                 compute gq = - r' * (transformed residual)
            call sscmv (nrank,-1d0,res,gq)
            call dtrmv ('u','t','n',nrank,r,ldr,gq)
            if (nrank.lt.n) call dgemv ('t',nrank,n-nrank,-1d0,
     *                                 r(1,nrank+1),ldr,res,0d0,
     *                                 gq(nrank+1))
            if (linobj) call daxpy (n,1d0,cq,1,gq,1)
            unitgz = .false.
            rownrm = dnrm2 (n,r(1,1),ldr)
            if (rownrm.le.tolrnk .or. abs(r(1,1)).le.rownrm*tolrnk) then
               nrz = 0
            else
               nrz = isrank(min(nrank,nz),r,ldr+1,tolrnk)
            end if
         end if
         singlr = .false.
      end if
c                                 end of lsgset
      end

      subroutine lsfeas(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *                  ax,bl,bu,featol,x,work)
c----------------------------------------------------------------------
c     lsfeas  computes the number of constraints that are violated by more
c     than  featol  and the 2-norm of the constraint violations.
c----------------------------------------------------------------------
      implicit none

      integer jmax, n, nclin, nviol, istate(n+nclin), i, is, j, idamax

      double precision ax(*), bl(n+nclin), bu(n+nclin), featol(n+nclin),
     *                 work(n+nclin), x(n), bigbnd, cvnorm, errmax,
     *                 biglow, bigupp, con, feasj, res, tolj, dnrm2

      external dnrm2, idamax
c----------------------------------------------------------------------
      biglow = -bigbnd
      bigupp = bigbnd

c     compute nviol,  the number of constraints violated by more than
c     featol,  and cvnorm,  the 2-norm of the constraint violations and
c     residuals of the constraints in the working set.

      nviol = 0

      do j = 1, n + nclin

         feasj = featol(j)
         is = istate(j)
         res = 0d0

         if (is.ge.0 .and. is.lt.4) then
            if (j.le.n) then
               con = x(j)
            else
               i = j - n
               con = ax(i)
            end if

            tolj = feasj

c           check for constraint violations.

            if (bl(j).gt.biglow) then
               res = bl(j) - con
               if (res.gt.feasj) nviol = nviol + 1
               if (res.gt.tolj) go to 20
            end if

            if (bu(j).lt.bigupp) then
               res = bu(j) - con
               if (res.lt.(-feasj)) nviol = nviol + 1
               if (res.lt.(-tolj)) go to 20
            end if

c           this constraint is satisfied,  but count the residual as a
c           violation if the constraint is in the working set.

            if (is.le.0) res = 0d0
            if (is.eq.1) res = bl(j) - con
            if (is.ge.2) res = bu(j) - con
            if (abs(res).gt.feasj) nviol = nviol + 1
         end if

   20    work(j) = res

      end do

      jmax = idamax(n+nclin,work)
      errmax = abs(work(jmax))

      cvnorm = dnrm2 (n+nclin,work,1)
c                                 end of lsfeas
      end

      subroutine cmmul2 (n,nrz,nz,zerolm,notopt,numinf,trusml,
     *                   smllst,jsmlst,tinyst,jtiny,gq)
c----------------------------------------------------------------------
c     cmmul2  updates jsmlst and smllst when there are artificial
c     constraints.
c     on input,  jsmlst  is the index of the minimum of the set of
c     adjusted multipliers.
c     on output, a negative jsmlst defines the index in q'g of the
c     artificial constraint to be deleted.
c----------------------------------------------------------------------
      implicit none

      integer jsmlst, jtiny, n, notopt, nrz, numinf, nz, j

      double precision smllst, tinyst, trusml, zerolm, gq(n), rlam
c----------------------------------------------------------------------
      do j = nrz + 1, nz
         rlam = -abs(gq(j))

         if (rlam.lt.zerolm) then
            if (numinf.eq.0) notopt = notopt + 1

            if (rlam.lt.smllst) then
               trusml = gq(j)
               smllst = rlam
               jsmlst = -j
            end if

         else if (rlam.lt.tinyst) then
            tinyst = rlam
            jtiny = -j
         end if
      end do
c                                 end of cmmul2
      end

      subroutine nggnfm (n,k1,k2,s,a,lda)
c----------------------------------------------------------------------
c  nggnfm applies a  sequence  of  pairwise interchanges to either  the
c  left,  or the right,  of the  n by n  upper triangular matrix  u,  to
c  transform u to an  upper hessenberg matrix. the interchanges are
c  applied in planes k1 up to k2.

c  the upper hessenberg matrix, h, is formed as

c     h = u*p',

c  where p is a permutation matrix of the form

c     p = p(k2 - 1)*...*p(k1+1)*p(k1),

c  p(k) being a pairwise interchange for the  (k, k+1) plane.
c  the  two by two
c  interchange part of p(k), r(k), is assumed to have the form

c     r(k) = (0  1).
c              (1  0)

c  the matrix  u must be supplied in the n by n leading upper triangular
c  part of the array  a, and this is overwritten by the upper triangular
c  part of  h.

c  the  sub-diagonal elements of  h, h(k + 1, k),  are returned in the
c  elements s(k),  k = k1, k1 + 1, ..., k2 - 1.

c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, n, i, j

      double precision a(lda,*), s(*), temp
c----------------------------------------------------------------------
      if ((min(n,k1).lt.1) .or. (k2.le.k1) .or. (k2.gt.n)) return

c        apply  the  plane interchanges to  columns  k1  up to
c        (k2 - 1) and  form   the   additional  sub-diagonal
c        elements,   storing  h(j + 1, j) in s(j).

         do j = k1, k2 - 1

            do i = 1, j
               temp = a(i,j+1)
               a(i,j+1) = a(i,j)
               a(i,j) = temp
            end do

            s(j) = a(j+1,j+1)
            a(j+1,j+1) = 0d0

         end do

c                                 end of nggnfm.
      end

      subroutine srchq (first,done,imprvd,inform,maxf,numf,
     *                  alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,
     *                  tolrel,toltny,alfa,alfbst,fbest)
c----------------------------------------------------------------------
c     srchq  finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     srchq  requires  f(alpha) (but not f'(alpha)) to be evaluated
c     in the interval.  new estimates of a minimizer are computed using
c     safeguarded quadratic interpolation.

c     reverse communication is used to allow the calling program to
c     evaluate f.  some of the parameters must be set or tested by the
c     calling program.  the remainder would ordinarily be local
c     variables.

c     first         must be .true. on the first entry.
c                   it is subsequently altered by srchq.

c     maxf          is an upper limit on the number of times srchq is
c                   to be entered consecutively with done = .false.
c                   (following an initial entry with first = .true.).

c     alfa          is the first estimate of a minimizer.  alfa is
c                   subsequently altered by srchq (see below).

c     alfmax        is the upper limit of the interval to be searched.

c     alfsml        is intended to prevent inefficiency when a minimizer
c                   is very small, for cases where the calling program
c                   would prefer to redefine f'(alfa).  alfsml is
c                   allowed to be zero.  early termination will occur if
c                   srchq determines that a minimizer lies somewhere in
c                   the interval [0, alfsml) (but not if alfmax is
c                   smaller that alfsml).

c     epsaf         is an estimate of the absolute precision in the
c                   computed value of f(0).

c     ftry          the value of f at the new point
c                   alfa = alfbst + xtry.

c     g0            is the value of f'(0).  g0 must be negative.

c     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
c                   such that if f has already been evaluated at alfa,
c                   it will not be evaluated closer than tol(alfa).
c                   these values may be reduced by srchc .

c     targtg        is the target value of abs(f'(alfa)). the search
c                   is terminated when
c                    abs(f'(alfa)) le targtg and f(alfa) lt 0.

c     toltny        is the smallest value that tolabs is allowed to be
c                   reduced to.

c     output parameters (relevant to the calling program)
c     ---------------------------------------------------

c     imprvd        is .true. if the previous alfa was the best point so
c                   far.  any related quantities should be saved by the
c                   calling program (e.g., arrays) before paying
c                   attention to the variable done.

c     done = .false.  means the calling program should evaluate ftry
c                   for the new trial step alfa, and reenter srchq.

c     done = .true.   means that no new alfa was calculated.  the value
c                   of inform gives the result of the search as follows

c                   inform = 1 means the search has terminated
c                              successfully with alfbst < alfmax.

c                   inform = 2 means the search has terminated
c                              successfully with alfbst = alfmax.

c                   inform = 3 means that the search failed to find a
c                              point of sufficient decrease in maxf
c                              functions, but a lower point was found.

c                   inform = 4 means alfmax is so small that a search
c                              should not have been attempted.

c                   inform = 5 means that the search was terminated
c                              because of alfsml (see above).

c                   inform = 6 means the search has failed to find a
c                              useful step.  the interval of uncertainty
c                              is [0,b] with b < 2*tolabs. a minimizer
c                              lies very close to alfa = 0, or f'(0) is
c                              not sufficiently accurate.

c                   inform = 7 if no better point could be found after
c                              maxf  function calls.

c                   inform = 8 means the input parameters were bad.
c                              alfmax le toltny  or  g0 ge zero.
c                              no function evaluations were made.

c     numf          counts the number of times srchq has been entered
c                   consecutively with done = .false. (i.e., with a new
c                   function value ftry).

c     alfa          is the point at which the next function ftry must
c                   be computed.

c     alfbst        should be accepted by the calling program as the
c                   approximate minimizer, whenever srchq returns
c                   inform = 1, 2 or 3.

c     fbest         will be the corresponding value of f.

c     the following parameters retain information between entries
c     -----------------------------------------------------------
c
c     braktd        is .false. if f has not been evaluated at the far
c                   end of the interval of uncertainty.  in this case,
c                   the point b will be at alfmax + tol(alfmax).
c
c     crampd        is .true. if alfmax is very small (le tolabs).  if
c                   the search fails, this indicates that a zero step
c                   should be taken.
c
c     extrap        is .true. if alfbst has moved at least once and xv
c                   lies outside the interval of uncertainty.  in this
c                   case, extra safeguards are applied to allow for
c                   instability in the polynomial fit.
c
c     moved         is .true. if a better point has been found, i.e.,
c                   alfbst gt 0.

c     vset          records whether a third-best point has been defined.

c     wset          records whether a second-best point has been
c                   defined.  it will always be .true. by the time the
c                   convergence test is applied.

c     nsamea        is the number of consecutive times that the
c                   left-hand end point of the interval of uncertainty
c                   has remained the same.

c     nsameb        similarly for the right-hand end.

c     a, b, alfbst  define the current interval of uncertainty.
c                   a minimizer lies somewhere in the  interval
c                   [alfbst + a, alfbst + b].

c     alfbst        is the best point so far.  it lies strictly within
c                   [atrue,btrue]  (except when alfbst has not been
c                   moved, in which case it lies at the left-hand end
c                   point).  hence we have a .le. 0 and b.gt.0.

c     fbest         is the value of f at the point alfbst.

c     fa            is the value of f at the point alfbst + a.

c     factor        controls the rate at which extrapolated estimates of
c                   alfa  may expand into the interval of uncertainty.
c                   factor is not used if a minimizer has been bracketed
c                   (i.e., when the variable braktd is .true.).

c     fv, fw        are the values of f at the points alfbst + xv  and
c                   alfbst + xw.  they are not defined until  vset  or
c                   wset  are .true..

c     xtry          is the trial point within the shifted interval
c                   (a, b).  the new trial function value must be
c                   computed at the point alfa = alfbst + xtry.

c     xv            is such that alfbst + xv is the third-best point.
c                   it is not defined until vset is .true..

c     xw            is such that alfbst + xw is the second-best point.
c                   it is not defined until wset is .true..  in some
c                   cases,  xw will replace a previous xw that has a
c                   lower function but has just been excluded from
c                   (a,b).

c     closef     is .true. if the worst function fv is within epsaf of
c                fbest (up or down).
c
c     found      is .true. if the sufficient decrease conditions holds
c                at alfbst.
c
c     quitf      is .true. when  maxf  function calls have been made.
c
c     quitfz     is .true. when the three best function values are
c                within epsaf of each other, and the new point satisfies
c                fbest le ftry le fbest+epsaf.
c
c     quiti      is .true. when the interval of uncertainty is less than
c                2*tol.

c     quits      is .true. as soon as alfa is too small to be useful;
c                i.e., btrue le alfsml.

c     xinxw      is .true. if xtry is in (xw,0) or (0,xw).
c----------------------------------------------------------------------
      implicit none

      logical badfun, braktd, closef, crampd, extrap, found,
     *                  moved, quitf, quitfz, quiti, quits, setxv, vset,
     *                  wset, xinxw,done, first, imprvd

      integer inform, maxf, numf, nsamea, nsameb

      double precision alfa, alfbst, alfmax, alfsml, epsaf, fbest,
     *                 ftry, g0, targtg, tolabs, tolrel, toltny, a, 
     *                 artifa, artifb, b, daux, dtry, endpnt, fa,
     *                 factor, fv, fw, gv, gw, q, s, tol, tolmax,
     *                 truea, trueb, xmidpt, xtry, xv, xw

      save              braktd, crampd, extrap, moved, vset, wset,
     *                  nsamea, nsameb, a, b, fa, factor, xtry, xw, fw,
     *                  xv, fv, tolmax
c----------------------------------------------------------------------
      imprvd = .false.
      badfun = .false.
      quitf = .false.
      quitfz = .false.
      quits = .false.
      quiti = .false.

      if (first) then

c        first entry.  initialize various quantities, check input data
c        and prepare to evaluate the function at the initial step alfa.

         first = .false.
         numf = 0
         alfbst = 0d0
         badfun = alfmax .le. toltny .or. g0 .ge. 0d0
         done = badfun
         moved = .false.

         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            vset = .false.
            wset = .false.
            nsamea = 0
            nsameb = 0

            tolmax = tolrel*alfmax + tolabs
            a = 0d0
            b = alfmax + tolmax
            fa = 0d0
            factor = 5d0
            tol = tolabs
            xtry = alfa
         end if
      else

c        subsequent entries.  the function has just been evaluated at
c        alfa = alfbst + xtry,  giving ftry.

         numf = numf + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b = alfmax - alfbst + tolmax
         end if

c        check if xtry is in the interval (xw,0) or (0,xw).

         if (wset) then
            xinxw = 0d0.lt.xtry .and. xtry .le. xw .or. xw .le.
     *              xtry .and. xtry.lt.0d0
         else
            xinxw = .false.
         end if

         imprvd = ftry.lt.fbest
         if (vset) then
            closef = abs(fbest-fv) .le. epsaf
         else
            closef = .false.
         end if

         if (imprvd) then

c           we seem to have an improvement.  the new point becomes the
c           origin and other points are shifted accordingly.

            if (wset) then
               xv = xw - xtry
               fv = fw
               vset = .true.
            end if

            xw = 0d0 - xtry
            fw = fbest
            wset = .true.
            fbest = ftry
            alfbst = alfa
            moved = .true.

            a = a - xtry
            b = b - xtry
            extrap = .not. xinxw

c           decrease the length of (a,b).

            if (xtry.ge.0d0) then
               a = xw
               fa = fw
               nsamea = 0
            else
               b = xw
               nsameb = 0
               braktd = .true.
            end if
         else if (closef .and. ftry-fbest.lt.epsaf) then

c           quit if there has been no progress and ftry, fbest, fw
c           and fv are all within epsaf of each other.

            quitfz = .true.
         else

c           the new function value is no better than the current best
c           point.  xtry must an end point of the new (a,b).

            if (xtry.lt.0d0) then
               a = xtry
               fa = ftry
               nsamea = 0
            else
               b = xtry
               nsameb = 0
               braktd = .true.
            end if

c           the origin remains unchanged but xtry may qualify as xw.

            if (wset) then
               if (ftry.lt.fw) then
                  xv = xw
                  fv = fw
                  vset = .true.

                  xw = xtry
                  fw = ftry
                  if (moved) extrap = xinxw
               else if (moved) then
                  if (vset) then
                     setxv = ftry.lt.fv .or. .not. extrap
                  else
                     setxv = .true.
                  end if

                  if (setxv) then
                     if (vset .and. xinxw) then
                        xw = xv
                        fw = fv
                     end if
                     xv = xtry
                     fv = ftry
                     vset = .true.
                  end if
               else
                  xw = xtry
                  fw = ftry
               end if
            else
               xw = xtry
               fw = ftry
               wset = .true.
            end if
         end if

c        check the termination criteria.

         tol = tolabs + tolrel*alfbst
         truea = alfbst + a
         trueb = alfbst + b
c
         found = moved .and. abs(fa-fbest) .le. -a*targtg
         quitf = numf .ge. maxf
         quiti = b - a .le. tol + tol
         quits = trueb .le. alfsml

         if (quiti .and. .not. moved) then

c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.

            tol = tol/1d1
            tolabs = tol
            quiti = abs(fw) .le. epsaf .or. tol .le. toltny
         end if

         done = quitf .or. quitfz .or. quits .or. quiti .or. found

c        proceed with the computation of an estimate of a minimizer.
c        the choices are...
c        1. parabolic fit using function values only.
c        2. damped parabolic fit if the regular fit appears to be
c           consistently overestimating the distance to a minimizer.
c        3. bisection, geometric bisection, or a step of tol if the
c           parabolic fit is unsatisfactory.

         if (.not. done) then
            xmidpt = (a+b)/2d0
            s = 0d0
            q = 0d0

c           fit a parabola.
c           see if there are two or three points for the parabolic fit.

            gw = (fw-fbest)/xw
            if (vset .and. moved) then

c              three points available.  use fbest, fw and fv.

               gv = (fv-fbest)/xv
               s = gv - (xv/xw)*gw
               q = 2d0 * (gv-gw)

            else

c              only two points available.  use fbest, fw and g0.

               if (moved) then
                  s = g0 - 2d0 * gw
               else
                  s = g0
               end if
               q = 2d0 * (g0-gw)

            end if

c           construct an artificial interval (artifa, artifb) in which
c           the new estimate of the steplength must lie.  set a default
c           value of  xtry  that will be used if the polynomial fit is
c           rejected. in the following, the interval (a,b) is considered
c           the sum of two intervals of lengths  dtry  and  daux, with
c           common end point the best point (zero).  dtry is the length
c           of the interval into which the default xtry will be placed
c           and endpnt denotes its non-zero end point.  the magnitude of
c           xtry is computed so that the exponents of dtry and daux are
c           approximately bisected.

            artifa = a
            artifb = b
            if (.not. braktd) then

c              a minimizer has not yet been bracketed.
c              set an artificial upper bound by expanding the interval
c              xw  by a suitable factor.

               xtry = -factor*xw
               artifb = xtry
               if (alfbst+xtry.lt.alfmax) factor = 5d0 * factor
            else if (vset .and. moved) then

c              three points exist in the interval of uncertainty.
c              check if the points are configured for an extrapolation
c              or an interpolation.

               if (extrap) then

c                 the points are configured for an extrapolation.

                  if (xw.lt.0d0) endpnt = b
                  if (xw.gt.0d0) endpnt = a
               else

c                 if the interpolation appears to be overestimating the
c                 distance to a minimizer,  damp the interpolation step.

                  if (nsamea.ge.3 .or. nsameb.ge.3) then
                     factor = factor / 5d0
                     s = factor*s
                  else
                     factor = 1d0
                  end if

c                 the points are configured for an interpolation.  the
c                 artificial interval will be just (a,b).  set endpnt so
c                 that xtry lies in the larger of the intervals (a,b)
c                 and  (0,b).

                  if (xmidpt.gt.0d0) then
                     endpnt = b
                  else
                     endpnt = a
                  end if

c                 if a bound has remained the same for three iterations,
c                 set endpnt so that  xtry  is likely to replace the
c                 offending bound.

                  if (nsamea.ge.3) endpnt = a
                  if (nsameb.ge.3) endpnt = b
               end if

c              compute the default value of  xtry.

               dtry = abs(endpnt)
               daux = b - a - dtry
               if (daux.ge.dtry) then
                  xtry = 5d0 *dtry*(0.1d0 + dtry/daux) / 11d0
               else
                  xtry = sqrt(daux)*sqrt(dtry) / 2d0
               end if
               if (endpnt.lt.0d0) xtry = -xtry

c              if the points are configured for an extrapolation set the
c              artificial bounds so that the artificial interval lies
c              within (a,b).  if the polynomial fit is rejected,  xtry
c              will remain at the relevant artificial bound.

               if (extrap) then
                  if (xtry.le.0d0) then
                     artifa = xtry
                  else
                     artifb = xtry
                  end if
               end if
            else

c              the gradient at the origin is being used for the
c              polynomial fit.  set the default xtry to one tenth xw.

               if (extrap) then
                  xtry = -xw
               else
                  xtry = xw / 1d1
               end if

            end if

c           the polynomial fits give (s/q)*xw as the new step.  reject
c           this step if it lies outside (artifa, artifb).

            if (q.ne.0d0) then
               if (q.lt.0d0) s = -s
               if (q.lt.0d0) q = -q
               if (s*xw.ge.q*artifa .and. s*xw.le.q*artifb) then

c                 accept the polynomial fit.

                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = 0d0
                  end if

               end if
            end if
         end if
      end if

      if (.not. done) then
         alfa = alfbst + xtry
         if (braktd .or. alfa.lt.alfmax-tolmax) then

c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)

            xmidpt = (a+b) / 2d0
            if (xtry.le.a+tol .or. xtry.ge.b-tol) then
               if (xmidpt.le.0d0) then
                  xtry = -tol
               else
                  xtry = tol
               end if
            end if

            if (abs(xtry).lt.tol) then
               if (xmidpt.le.0d0) then
                  xtry = -tol
               else
                  xtry = tol
               end if
            end if
            alfa = alfbst + xtry
         else

c           the step is close to or larger than alfmax, replace it by
c           alfmax to force evaluation of the function at the boundary.

            braktd = .true.
            xtry = alfmax - alfbst
            alfa = alfmax
         end if
      end if

c     exit.

      if (done) then
         if (badfun) then
            inform = 8
         else if (quits) then
            inform = 5
         else if (found) then
            if (alfbst.lt.alfmax) then
               inform = 1
            else
               inform = 2
            end if
         else if (moved) then
            inform = 3
         else if (quitf) then
            inform = 7
         else if (crampd) then
            inform = 4
         else
            inform = 6
         end if
      end if
c                                 end of srchq
      end

      subroutine cmmul1 (n,lda,ldt,nactiv,nfree,nz,istate,
     *                   kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                   jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                   jbigst,kbigst,a,anorms,gq,rlamda,t,wtinf)
c----------------------------------------------------------------------
c     cmmul1 computes the lagrange multiplier estimates for the
c     given working set.  it then determines the values and indices of
c     certain significant multipliers.  in this process, the multipliers
c     for inequalities at their upper bounds are adjusted so that a
c     negative multiplier for an inequality constraint indicates non-
c     optimality.  all adjusted multipliers are scaled by the 2-norm
c     of the associated constraint row.  in the following, the term
c     minimum refers to the ordering of numbers on the real line,  and
c     not to their magnitude.
c
c     jsmlst          is the index of the constraint whose multiplier is
c                     the minimum of the set of adjusted multipliers
c                     with values less than  small.
c     rlamda(ksmlst)  is the associated multiplier.
c
c     jbigst          is the index of the constraint whose multiplier is
c                     the largest of the set of adjusted multipliers
c                     with values greater than (1 + small).
c     rlamda(kbigst)  is the associated multiplier.
c
c     on exit,  elements  1  thru  nactiv  of  rlamda  contain the
c     unadjusted multipliers for the general constraints.  elements
c     nactiv  onwards of  rlamda  contain the unadjusted multipliers
c     for the bounds.
c----------------------------------------------------------------------
      implicit none

      integer jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *        lda, ldt, n, nactiv, nfree, notopt, numinf, nz,
     *        istate(*), kactiv(n), kx(n), i, is, j, k, l, nfixed

      double precision biggst, smllst, tinyst, trubig, trusml, zerolm,
     *                 a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                 wtinf(*), anormj, blam, rlam, scdlam
c----------------------------------------------------------------------
      nfixed = n - nfree
      jtiny = 0
      jsmlst = 0
      ksmlst = 0
      jbigst = 0
      kbigst = 0

c     compute  jsmlst  for regular constraints and temporary bounds.

c     first, compute the lagrange multipliers for the general
c     constraints in the working set, by solving  t'*lamda = y'g.

      if (n.gt.nz) call dcopy (n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call dtrsv ('u','t','n',nactiv,t(1,nz+1),ldt,
     *                            rlamda)

c     set elements  nactiv, nactiv+1,... of  rlamda  equal to
c     the multipliers for the bound constraints.

      do l = 1, nfixed

         j = kx(nfree+l)
         blam = rlamda(nactiv+l)

         do k = 1, nactiv
            i = kactiv(k)
            blam = blam - a(i,j)*rlamda(nactiv-k+1)
         end do

         rlamda(nactiv+l) = blam

      end do

c     find  jsmlst  and  ksmlst.

      do k = 1, n - nz
         if (k.gt.nactiv) then
            j = kx(nz+k)
         else
            j = kactiv(nactiv-k+1) + n
         end if

         is = istate(j)
c
         i = j - n
         if (j.le.n) anormj = 1d0
         if (j.gt.n) anormj = anorms(i)
c
         rlam = rlamda(k)
c
c        change the sign of the estimate if the constraint is in
c        the working set at its upper bound.
c
         if (is.eq.2) rlam = -rlam
         if (is.eq.3) rlam = abs(rlam)
         if (is.eq.4) rlam = -abs(rlam)
c
         if (is.ne.3) then
            scdlam = rlam*anormj
c
            if (scdlam.lt.zerolm) then
               if (numinf.eq.0) notopt = notopt + 1
c
               if (scdlam.lt.smllst) then
                  smllst = scdlam
                  trusml = rlamda(k)
                  jsmlst = j
                  ksmlst = k
               end if
            else if (scdlam.lt.tinyst) then
               tinyst = scdlam
               jtiny = j
            end if
         end if
c
         scdlam = rlam/wtinf(j)
         if (scdlam.gt.biggst .and. j.gt.jinf) then
            biggst = scdlam
            trubig = rlamda(k)
            jbigst = j
            kbigst = k
         end if
      end do
c                                 end of cmmul1
      end

      subroutine npiqp (feasqp,unitq,nqperr,n,nclin,
     *              ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,
     *              numinf,istate,kactiv,kx,dxnorm,gdx,qpcurv,aqp,
     *              adx,ax,bl,bu,clamda, dx,qpbl,qpbu,qptol,r,x,wtinf,w)
c----------------------------------------------------------------------
c     npiqp    does the following:

c     (1)  generate the upper and lower bounds for the qp  subproblem.

c     (2)  compute the  tq  factors of the rows of  aqp  specified by
c          the array  istate.  the part of the factorization defined by
c          the first contiguous group of linear constraints does not
c          need to be recomputed.  the remaining rows (which could be
c          comprised of both linear and nonlinear constraints) are
c          included as new rows of the  tq  factorization stored in
c          t and zy.  note that if there are no nonlinear constraints,
c          no factorization is required.

c     (3)  solve the  qp  subproblem.
c                 minimize     1/2 (w p - d)'(wp - d) + g'p

c                 subject to   qpbl .le. ( p) .le. qpbu,
c                                        (ap)

c          where  w  is a matrix (not stored) such that  w'w = h  and
c          wq = r,  d  is the zero vector,  and  g  is the gradient.
c          if the subproblem is infeasible, compute the point which
c          minimizes the sum of infeasibilities.

c     (4)   find the value of each slack variable for which the merit
c          function is minimized.

c     (5)  compute  dx,  the search direction for
c          the slack variables, the multipliers and the variables.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical feasqp, unitq, linobj

      integer ldaqp, ldr, linact, minits, n, nactiv, nclin,
     *        nfree, nlnact, nqperr, numinf, nz, istate(*), kactiv(n), 
     *        kx(n), j, jinf, nctotl, nmajor, nminor, nrank,
     *        nrejtd, ntry, nviol, nz1

      double precision adx(*), aqp(ldaqp,*), ax(*), bl(*),
     *                 bu(*), clamda(*), dx(n), qpbl(*),
     *                 qpbu(*), qptol(*), r(ldr,*), 
     *                 w(*), wtinf(*), x(n), dxnorm, gdx, 
     *                 biglow, bigupp, blj, buj, con, condmx,
     *                 ssq, ssq1, suminf, weight, wscale,
     *                 wtmax, wtmin, qpcurv, ddot, dnrm2, sdiv

      external ddot, dnrm2, sdiv

      integer lkactv, lanorm, lcjdx, lqpdx, lrpq, lrpq0, lqphz,
     *        lhpq, lgq, lrlam, lt, lzy, lwtinf, lwrk1, lqptol
      common/ cstlnp /lkactv, lanorm, lcjdx, lqpdx, lrpq, lrpq0, lqphz,
     *                lhpq, lgq, lrlam, lt, lzy, lwtinf, lwrk1, lqptol

      integer ldt, ncolt, ldzy
      common/ ngg004 /ldt, ncolt, ldzy

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin

      integer itmax1, itmax2, lprob
      common/ ngg016 /itmax1, itmax2, lprob

      double precision rcndbd, rfrobn, drmax, drmin
      common/ ngg018 /rcndbd, rfrobn, drmax, drmin

      double precision cdint, ctol, dxlim, epsrf, eta, fdint, ftol,
     *                 hcndbd
      common/ ngg021 /cdint, ctol, dxlim, epsrf, eta,
     *                fdint, ftol, hcndbd

      equivalence (itmxnp,nmajor), (itmax2,nminor)
c----------------------------------------------------------------------
      feasqp = .true.
      linobj = .true.

      biglow = -bigbnd
      bigupp = bigbnd
      ssq1 = 0d0

      nctotl = n + nclin
      nrank = n
      nrejtd = 0

c     generate the upper and lower bounds upon the search direction, the
c     weights on the sum of infeasibilities and the nonlinear constraint
c     violations.

      wscale = -1d0

      do j = 1, nctotl

         if (j.le.n) then
            con = x(j)
         else
            con = ax(j-n)
         end if

         blj = bl(j)
         buj = bu(j)
         if (blj.gt.biglow) blj = blj - con
         if (buj.lt.bigupp) buj = buj - con

         weight = 1d0

         if (abs(blj).le.qptol(j)) blj = 0d0
         if (abs(buj).le.qptol(j)) buj = 0d0

         wtinf(j) = weight
         qpbl(j) = blj
         qpbu(j) = buj

      end do

      if (wscale.gt.0d0) then
         wscale = 1d0/wscale
         call dscal (nctotl,wscale,wtinf,1)
      end if

      call scond (nctotl,wtinf,1,wtmax,wtmin)
      wtmin = epspt9*wtmax

      do j = 1, nctotl
         wtinf(j) = max(wtinf(j),wtmin)
      end do 

c     set the maximum allowable condition estimator of the constraints
c     in the working set.  note that a relatively well-conditioned
c     working set is used to start the qp iterations.

      condmx = 1d0/epspt3

c     solve for dx, the vector of minimum two-norm that satisfies the
c     constraints in the working set.

      call npsetx(unitq,nclin,nactiv,nfree,nz,n,nctotl,ldzy,ldaqp,
     *            ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,adx,qpbl,qpbu,
     *            w(lrpq),w(lrpq0),dx,w(lgq),r,w(lt),w(lzy),w(lwrk1))

c     solve a quadratic program for the search direction  dx  and
c     multiplier estimates  clamda.

c     if there is no feasible point for the subproblem,  the sum of
c     infeasibilities is minimized subject to the linear constraints
c     (1  thru  jinf)  being satisfied.

      jinf = nctotl

      ntry = 1

      do

         call lscore('qp subproblem',linobj,unitq,nqperr,
     *            minits,jinf,nclin,nctotl,nactiv,nfree,nrank,nz,nz1,n,
     *            ldaqp,ldr,istate,kactiv,kx,gdx,ssq,ssq1,suminf,numinf,
     *            dxnorm,qpbl,qpbu,aqp,clamda,adx,qptol,r,dx,w)

         nviol = 0

         if (numinf.gt.0) then
c                                 count the violated linear constraints.
            do j = 1, nctotl
               if (istate(j).lt.0) nviol = nviol + 1
            end do

            if (nviol.gt.0) then
               ntry = ntry + 1
               unitq = .true.
               nactiv = 0
               nfree = n
               nz = n
               call iload (nctotl,0,istate,1)

               call npsetx (unitq,nclin,nactiv,nfree,nz,n,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,qpbl,qpbu,w(lrpq),w(lrpq0),dx,w(lgq),r,w(lt)
     *                  ,w(lzy),w(lwrk1))
            end if

         end if

         if (nviol.eq.0 .or. ntry.gt.2) exit

      end do

      nlnact = 0
      linact = nactiv - nlnact

c     compute  hpq = r'r(pq)  from the transformed gradient of the qp
c     objective function and  r(pq)  from the transformed residual.

      call dscal (n,-1d0,w(lrpq),1)
      call daxpy (n,-1d0,w(lgq),1,w(lhpq),1)
      qpcurv = 2d0 * ssq
c                                 end of npiqp
      end

      subroutine cmsinf (n,nclin,lda,istate,bigbnd,numinf,suminf,bl,bu,
     *                   a,featol,cvec,x,wtinf)
c----------------------------------------------------------------------
c     cmsinf  finds the number and weighted sum of infeasibilities for
c     the bounds and linear constraints.   an appropriate gradient
c     is returned in cvec.

c     positive values of  istate(j)  will not be altered.  these mean
c     the following...

c               1             2           3
c           a'x = bl      a'x = bu     bl = bu

c     other values of  istate(j)  will be reset as follows...
c           a'x lt bl     a'x gt bu     a'x free
c              - 2           - 1           0
c----------------------------------------------------------------------
      implicit none

      integer lda, n, nclin, numinf, istate(*), j, k

      double precision a(lda,*), bl(*), bu(*), cvec(n), featol(*),
     *                 wtinf(*), x(n), biglow, bigupp, ctx, feasj, s, 
     *                 weight, bigbnd, suminf, ddot
      external ddot
c----------------------------------------------------------------------
      bigupp = bigbnd
      biglow = -bigbnd

      numinf = 0
      suminf = 0d0
      call sload (n,0d0,cvec,1)

      do 40 j = 1, n + nclin
         if (istate(j).le.0) then
            feasj = featol(j)
            if (j.le.n) then
               ctx = x(j)
            else
               k = j - n
               ctx = ddot (n,a(k,1),lda,x)
            end if
            istate(j) = 0

c           see if the lower bound is violated.

            if (bl(j).gt.biglow) then
               s = bl(j) - ctx
               if (s.gt.feasj) then
                  istate(j) = -2
                  weight = -wtinf(j)
                  go to 20
               end if
            end if

c           see if the upper bound is violated.

            if (bu(j).ge.bigupp) go to 40
            s = ctx - bu(j)
            if (s.le.feasj) go to 40
            istate(j) = -1
            weight = wtinf(j)

c           add the infeasibility.

   20       numinf = numinf + 1
            suminf = suminf + abs(weight)*s
            if (j.le.n) then
               cvec(j) = weight
            else
               call daxpy (n,weight,a(k,1),lda,cvec,1)
            end if
         end if
   40 continue
c                                 end of cmsinf
      end

      subroutine rzdel (unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,t,gqm,q,c,s)
c----------------------------------------------------------------------
c     rzdel   updates the matrices  z, y, t, r  and  d  associated with
c     factorizations

c              a(free) * q(free)  = ( 0 t)
c                        q(free)  = ( z y)
c
c     when a regular, temporary or artificial constraint is deleted
c     from the working set.

c     the  nactiv x nactiv  upper-triangular matrix  t  is stored
c     with its (1,1) element in position  (it,jt)  of the array  t.
c----------------------------------------------------------------------
      implicit none

      logical unitq

      integer it, jdel, kdel, lda, ldq, ldt, n, nactiv, nfree,
     *        ngq, nrz, nz, kactiv(n), kx(n), i, ir, itdel, j, jart, jt,
     *        k, l, npiv, nrz1, nsup, idamax

      double precision a(lda,*), c(n), gqm(n,*), q(ldq,*), s(n),
     *                 t(ldt,*), cs, sn

      external idamax

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin
c----------------------------------------------------------------------
      jt = nz + 1

      if (jdel.gt.0) then

c        regular constraint or temporary bound deleted.

         if (jdel.le.n) then

c           case 1.  a simple bound has been deleted.
c           columns  nfree+1  and  ir  of gqm' must be swapped.

            ir = nz + kdel

            itdel = nactiv + 1
            nfree = nfree + 1
            if (nfree.lt.ir) then
               kx(ir) = kx(nfree)
               kx(nfree) = jdel
               call dswap(ngq,gqm(nfree,1),gqm(ir,1),n)
            end if

            if (.not. unitq) then

c              copy the incoming column of  a(free)  into the end of  t.

               do 20 k = 1, nactiv
                  i = kactiv(k)
                  t(nactiv-k+1,nfree) = a(i,jdel)
   20          continue

c              expand  q  by adding a unit row and column.

               if (nfree.gt.ldq) then
c debug debug 691
                  write (*,*) 'wtf nfree > ldq we are gonna crash'
               else
                  if (nfree.gt.1) then
                     call sload (nfree-1,0d0,q(nfree,1),ldq)
                     call sload (nfree-1,0d0,q(1,nfree),1)
                  end if
                  q(nfree,nfree) = 1d0
               end if
            end if
         else

c           case 2.  a general constraint has been deleted.
c           delete row  itdel  of  t  and move up the ones below it.
c           t  becomes lower hessenberg.

            itdel = kdel
            do k = itdel, nactiv
               j = jt + k - 1
               do l = itdel, k - 1
                  i = it + l - 1
                  t(i,j) = t(i+1,j)
               end do
            end do
c
            do 80 i = nactiv - itdel + 1, nactiv - 1
               kactiv(i) = kactiv(i+1)
   80       continue
            nactiv = nactiv - 1
         end if

         nz = nz + 1

         if (nactiv.eq.0) then
            dtmax = 1d0
            dtmin = 1d0
         else

c           restore the nactiv x (nactiv+1) upper-hessenberg matrix  t
c           to upper-triangular form.  the  nsup  super-diagonal
c           elements are removed by a backward sweep of rotations.
c           the rotation for the  (1,1)-th  element of  t  is generated
c           separately.

            nsup = itdel - 1

            if (nsup.gt.0) then
               npiv = jt + itdel - 1
               if (nsup.gt.1) then
                  call dcopy (nsup-1,t(it+1,jt+1),ldt+1,s(jt+1),1)
                  call suhqr('right',nactiv,1,nsup,c(jt+1),s(jt+1),
     *                        t(it,jt+1),ldt)
               end if

               call srotgc(t(it,jt+1),t(it,jt),cs,sn)
               t(it,jt) = 0d0
               s(jt) = -sn
               c(jt) = cs
               call sgesrc('right','variable','backwards',nfree,nfree,
     *                     nz,npiv,c,s,q,ldq)
               call sgesrc('left ','variable','backwards',npiv,ngq,nz,
     *                     npiv,c,s,gqm,n)
            end if

            jt = jt + 1
            call scond (nactiv,t(it,jt),ldt+1,dtmax,dtmin)
         end if
      end if

      nrz1 = nrz + 1

      if (nz.gt.nrz) then
         if (jdel.gt.0) then
            jart = nrz1 - 1 + idamax(nz-nrz1+1,gqm(nrz1,1))
         else
            jart = -jdel
         end if

         if (jart.gt.nrz1) then

c           swap columns  nrz1  and  jart  of  q  and  gqm.

            if (unitq) then
               k = kx(nrz1)
               kx(nrz1) = kx(jart)
               kx(jart) = k
            else
               call dswap(nfree,q(1,nrz1),q(1,jart),1)
            end if

            call dswap(ngq,gqm(nrz1,1),gqm(jart,1),n)
         end if
      end if

      nrz = nrz1
c                                 end of rzdel
      end

      subroutine cmchzr(firstv,n,nclin,istate,bigalf,bigbnd,pnorm,
     *                  hitlow,move,onbnd,unbndd,alfa,alfap,jhit,anorm,
     *                  ap,ax,bl,bu,featol,featlu,p,x)
c----------------------------------------------------------------------
c     cmchzr  finds a step alfa such that the point x + alfa*p reaches
c     one of the linear constraints (including bounds).

c     when x is infeasible, the number of
c     infeasibilities will never increase.  if the number stays the
c     same, the sum of infeasibilities will decrease.  if the number
c     decreases by one or more,  the sum of infeasibilities will usually
c     decrease also, but occasionally it will increase after the step
c     alfa  is taken.  (convergence is still assured because the number
c     has decreased.)

c     three possible steps are computed as:
c     alfaf = the maximum step that can be taken without violating
c              one of the constraints that are currently satisfied.
c     alfai = reaches a linear constraint that is currently violated.
c              usually this will be the furthest such constraint along
c              p, subject to the angle between the constraint normal and
c              p being reasonably close to the maximum value among
c              infeasible constraints,  but if firstv = .true. it will
c              be the first one along p.  the latter case applies only
c              when the problem has been determined to be infeasible,
c              and the sum of infeasibilities are being minimized.
c              (alfai is not defined when x is feasible.)
c     alfai is needed occasionally when infeasible, to prevent
c     going unnecessarily far when alfaf is quite large.  it will
c     always come into effect when x is about to become feasible.
c     (the sum of infeasibilities will decrease initially as alfa
c     increases from zero, but may start increasing for larger steps.
c     choosing a large alfai allows several elements of  x  to
c     become feasible at the same time.

c     in the end, alfa = alfaf  if x is feasible, or if
c     alfai > alfap (where  alfap  is the perturbed step from pass 1).
c     otherwise, alfa = alfai.

c     bigalf defines what should be treated as an unbounded step.
c     bigbnd provides insurance for detecting unboundedness.
c            if alfa reaches a bound as large as bigbnd, it is
c            classed as an unbounded step.
c     featol is the array of current feasibility tolerances used by
c            cmsinf.  typically in the range 0.5*tolx to 0.99*tolx,
c            where tolx is the featol specified by the user.
c     tolinc (in common) is used to determine stepmn (see below),
c            the minimum positive step.
c     istate is set as follows:
c            istate(j) = -2  if a'x.lt.bl - featol
c                      = -1  if a'x.gt.bu + featol
c                      =  0  if a'x is not in the working set
c                      =  1  if a'x is in the working set at bl
c                      =  2  if a'x is in the working set at bu
c                      =  3  if a'x is in the working set (an equality)
c                      =  4  if x(j) is temporarily fixed.
c            values -2 and -1 do not occur once feasible.
c     bl     the lower bounds on the variables.
c     bu     the upper bounds on ditto.
c     x      the values of       ditto.
c     p      the search direction.

c     hitlow  = true  if a lower bound restricted alfa.
c             = false otherwise.
c     move    = true  if  exact ge stepmn  (defined at end of code).
c     onbnd   = true  if  alfa = exact.  this means that the step  alfa
c                     moves x  exactly onto one of its constraints,
c                     namely  bound.
c             = false if the exact step would be too small
c                     (exact.lt.stepmn).
c               (with these definitions,  move = onbnd).
c     unbndd  = true  if alfa = bigalf.  jhit may possibly be zero.
c               the parameters hitlow, move, onbnd, bound and exact
c               should not be used.
c     jhit    = the index (if any) such that constraint jhit reaches
c               a bound.
c     bound   = the bound value bl(jhit) or bu(jhit) corresponding
c               to hitlow.
c     exact   = the step that would take constraint jhit exactly onto
c               bound.
c     alfa    = an allowable, positive step.
c               if unbndd is true,  alfa = stepmx.
c               otherwise,          alfa = max(stepmn, exact).

c     cmchzr is based on minos 5.2 routine m5chzr, which implements the
c     expand procedure to deal with degeneracy. the step alfaf is
c     chosen as in the two-pass approach of paula harris (1973), except
c     that this version insists on returning a positive step, alfa.
c     two features make this possible:

c        1. featol increases slightly each iteration.

c        2. the blocking constraint, when added to the working set,
c           retains the value ax(jhit) + alfa * ap(jhit),
c           even if this is not exactly on the blocking bound.

c     for infeasible variables moving towards their bound, we require
c     the rate of change of the chosen constraint to be at least gamma
c     times as large as the biggest available.  this still gives us
c     freedom in pass 2.
c     gamma = 0.1 and 0.01 seemed to inhibit phase 1 somewhat.
c     gamma = 0.001 seems to be safe.
c----------------------------------------------------------------------
      implicit none 

      double precision gamma
      parameter         (gamma=1.0d-3)

      logical firstv, hitlow, move, onbnd, unbndd, blockf, blocki

      integer jhit, n, nclin, istate(n+nclin), i, j, jhitf, jhiti, js,
     *        itnfix, kdegen, ndegen, nfix

      double precision alfa, alfap, bigalf, bigbnd, pnorm, anorm(*), 
     *                 ap(*), ax(*), bl(n+nclin), bu(n+nclin), 
     *                 featlu(n+nclin), featol(n+nclin), p(n), x(n), 
     *                 alfai, atp, atpabs, atpmxf, atpmxi, atpscd, atx,
     *                 biglow, bigupp, bound, delta, exact, res,
     *                 stepmn, tolpiv

      double precision tolinc, tolx0
      common/ ngg005 /tolx0, tolinc, kdegen, ndegen, itnfix, nfix(2)

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9
c----------------------------------------------------------------------

c     tolpiv is a tolerance to exclude negligible elements of a'p.

      biglow = -bigbnd
      bigupp = bigbnd
      tolpiv = epspt9*pnorm

c     first pass -- find steps to perturbed constraints, so that
c     alfap will be slightly larger than the true step.
c     in degenerate cases, this strategy gives us some freedom in the
c     second pass.  the general idea follows that described by p.m.j.
c     harris, p.21 of mathematical programming 5, 1 (1973), 1--28.

      atpmxi = 0d0
      alfap = bigalf

      do 20 j = 1, n + nclin
         js = istate(j)

         if (js.le.0) then
            delta = featol(j)

            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               atpabs = abs(atp)
               atpscd = atpabs
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               atpabs = abs(atp)
               atpscd = atpabs/(1d0+anorm(i))
            end if

            if (atpscd.le.tolpiv) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step. 

               res = -1d0

            else if (atp.le.0d0 .and. js.ne.-2) then

c              a'x  is decreasing and the lower bound is not violated.
c              test for smaller alfap.
c              if the upper bound is violated. test for bigger atp.

               if (bl(j).gt.biglow) then
                  res = atx - bl(j) + delta
                  if (res.lt.alfap*atpabs) alfap = res/atpabs
               end if

               if (js.eq.-1) atpmxi = max(atpmxi,atpscd)

            else if (atp.gt.0d0 .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.
c              test for smaller alfap.
c              if the lower bound is violated. test for bigger atp.

               if (bu(j).lt.bigupp) then
                  res = bu(j) - atx + delta

                  if (res.lt.alfap*atp) alfap = res/atp
               end if

               if (js.eq.-2) atpmxi = max(atpmxi,atpscd)
            end if
         end if
   20 continue

c     second pass.
c     for feasible variables, recompute steps without perturbation.
c     amongst constraints that are closer than alfap, choose the one
c     that makes the largest angle with the search direction.
c     for infeasible variables, find the largest step subject to a'p
c     being no smaller than gamma * max(a'p).

      if (firstv) then
         alfai = bigalf
      else
         alfai = 0d0
      end if

      atpmxf = 0d0
      atpmxi = gamma*atpmxi
      jhitf = 0
      jhiti = 0

      do 40 j = 1, n + nclin
         js = istate(j)

         if (js.le.0) then

            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               atpabs = abs(atp)
               atpscd = atpabs
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               atpabs = abs(atp)
               atpscd = atpabs/(1d0+anorm(i))
            end if

            if (atpscd.le.tolpiv) then

c              this constraint appears to be constant along p. it is
c              not used to compute the step.

               res = -1d0

            else if (atp.le.0d0 .and. js.ne.-2) then

c              a'x  is decreasing.
c              test for bigger a'p if the lower bound is satisfied.
c              test for smaller alfaf.

               if (atpscd.gt.atpmxf) then

                  if (bl(j).gt.biglow) then
                     res = atx - bl(j)

                     if (res.le.alfap*atpabs) then
                        atpmxf = atpscd
                        jhitf = j
                     end if
                  end if
               end if

               if (js.eq.-1) then

c                 the upper bound is violated.
c                 test for bigger or smaller alfai,  depending on the
c                 value of firstv.

                  if (firstv) then
                     res = atx - bu(j)

                     if (res.le.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if

                  else if (atpscd.ge.atpmxi) then
                     res = atx - bu(j)

                     if (res.gt.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if
                  end if
               end if

            else if (atp.gt.0d0 .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.
c              test for smaller alfap.

               if (atpscd.gt.atpmxf) then

                  if (bu(j).lt.bigupp) then
                     res = bu(j) - atx

                     if (res.le.alfap*atp) then
                        atpmxf = atpscd
                        jhitf = j
                     end if
                  end if
               end if

               if (js.eq.-2) then

c                 the lower bound is violated.
c                 test for bigger or smaller alfai,  depending on the
c                 value of firstv.

                  if (firstv) then
                     res = bl(j) - atx

                     if (res.le.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  else if (atpscd.ge.atpmxi) then
                     res = bl(j) - atx

                     if (res.gt.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  end if
               end if
            end if
         end if
   40 continue

c     see if a feasible and/or infeasible constraint blocks.

      blockf = jhitf.gt.0
      blocki = jhiti.gt.0
      unbndd = .not. (blockf .or. blocki)

      if (unbndd) go to 60

      if (blockf) then

c        a constraint is hit which is currently feasible.
c        the corresponding step alfaf is not used, so no need to get it,
c        but we know that alfaf .le. alfap, the step from pass 1.

         jhit = jhitf
         if (jhit.le.n) then
            atp = p(jhit)
         else
            atp = ap(jhit-n)
         end if
         hitlow = atp.lt.0d0
      end if

c     if there is a choice between alfaf and alfai, it is probably best
c     to take alfai.  however, we can't if alfai is bigger than alfap.

      if (blocki .and. alfai.le.alfap) then

c        an infeasible variable reaches its violated bound.

         jhit = jhiti
         if (jhit.le.n) then
            atp = p(jhit)
         else
            atp = ap(jhit-n)
         end if
         hitlow = atp.gt.0d0
      end if

c DEBUG DEBUG, seems to work

      if (jhit.eq.0) then 
         unbndd = .true.
         if (unbndd) go to 60
      end if 


      if (jhit.le.n) then
         atx = x(jhit)
      else
         atx = ax(jhit-n)
      end if


c     try to step exactly onto bound, but make sure the exact step
c     is sufficiently positive.  (exact will be alfaf or alfai.)
c     since featol increases by  tolinc  each iteration, we know that
c     a step as large as  stepmn  (below) will not cause any feasible
c     variables to become infeasible (where feasibility is measured
c     by the current featol).

      if (hitlow) then
         bound = bl(jhit)
      else
         bound = bu(jhit)
      end if

      unbndd = abs(bound) .ge. bigbnd
      if (unbndd) go to 60

      stepmn = tolinc*featlu(jhit)/abs(atp)
      exact = (bound-atx)/atp
      alfa = max(stepmn,exact)
      onbnd = alfa.eq.exact
      move = exact .ge. stepmn
      if (.not. move) ndegen = ndegen + 1

      return

c     unbounded.

   60 alfa = bigalf
      move = .true.
      onbnd = .false.
c                                 end of cmchzr
      end

      subroutine cmprt (nfree,n,nctotl,nactiv,kactiv,kx,clamda,rlamda)
c----------------------------------------------------------------------
c cmprt creates the expanded lagrange multiplier vector clamda.
c----------------------------------------------------------------------
      implicit none

      integer n, nactiv, nctotl, nfree, kactiv(n), kx(n), j, k, 
     *        nfixed, nz

      double precision clamda(nctotl), rlamda(n)
c----------------------------------------------------------------------
      nz = nfree - nactiv

      call sload (nctotl,0d0,clamda,1)
      nfixed = n - nfree

      do k = 1, nactiv + nfixed
         if (k.le.nactiv) j = kactiv(k) + n
         if (k.gt.nactiv) j = kx(nz+k)
         clamda(j) = rlamda(k)
      end do
c                                 end of cmprt
      end

      subroutine rzadd (unitq,rset,inform,ifix,iadd,jadd,it,nactiv,nz,
     *                  nfree,nrz,ngq,n,lda,ldq,ldr,ldt,kx,condmx,drzz,
     *                  a,r,t,gqm,q,w,c,s)
c----------------------------------------------------------------------
c     rzadd  updates the matrices  z, y, t, r  and  d  associated with
c     factorizations
c
c              a(free) * q(free)  = ( 0 t)
c                        q(free)  = ( z y)
c                      r' *d * r  =   hz
c
c     a) the matrices  r  and  t  are upper triangular.
c     b) the arrays  t  and  r  may be the same array.
c     c) the  nactiv x nactiv  upper-triangular matrix  t  is stored
c        with its (1,1) element in position  (it,jt) of the
c        array  t.   the integer  jt  is always  nz+1.  during regular
c        changes to the working set,  it = 1;  when several constraints
c        are added simultaneously,  it  points to the first row of the
c        existing  t.
c     d) the matrix  r  is stored in the first  nz x nz  rows
c        and columns of the  nfree x nfree  leading principal minor of
c        the array  r.
c     e) if  rset  is  false,   r  is not touched.
c
c     there are three separate cases to consider (although each case
c     shares code with another)...
c
c     (1) a free variable becomes fixed on one of its bounds when there
c         are already some general constraints in the working set.
c
c     (2) a free variable becomes fixed on one of its bounds when there
c         are only bound constraints in the working set.
c
c     (3) a general constraint (corresponding to row  iadd  of  a) is
c         added to the working set.
c
c     in cases (1) and (2), we assume that  kx(ifix) = jadd.
c     in all cases,  jadd  is the index of the constraint being added.
c
c     if there are no general constraints in the working set,  the
c     matrix  q = (z y)  is the identity and will not be touched.
c
c     if  ngq.gt.0,  the column transformations are applied to the
c     columns of the  (ngq x n)  matrix  gqm'.
c----------------------------------------------------------------------
      implicit none

      logical rset, unitq, bound, overfl


      integer iadd, ifix, inform, it, jadd, lda, ldq, ldr, nanew,
     *        ldt, n, nactiv, nfree, ngq, nrz, nz, kx(n), i, j, jt, k,
     *        npiv, nsup

      double precision a(lda,*), c(n), gqm(n,*), q(ldq,*), r(ldr,*),
     *                 s(n), t(ldt,*), w(n), condmx, drzz, cond, condbd,
     *                 dtnew, tdtmax, tdtmin, dnrm2, sdiv

      external dnrm2, sdiv 

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin
c----------------------------------------------------------------------
      condbd = 1d0/epspt9

      overfl = .false.
      bound = jadd .le. n
      jt = nz + 1

      if (bound) then
c                                 a simple bound has entered the working set. iadd is not used.
         nanew = nactiv

         if (unitq) then
c                                 q is not stored, but kx defines an ordering of the columns
c                                 of the identity matrix that implicitly define q.
c                                 define the sequence of pairwise interchanges p that moves
c                                 the newly-fixed variable to position  nfree.
c                                 reorder kx.
            do 20 i = 1, nfree - 1
               if (i.ge.ifix) then
                  w(i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if
   20       continue
         else
c                                 q is stored explicitly.
c                                 set w = (ifix)-th row of q.
c                                 move (nfree)-th row of q to position ifix.
            call dcopy (nfree,q(ifix,1),ldq,w,1)
            if (ifix.lt.nfree) then
               call dcopy (nfree,q(nfree,1),ldq,q(ifix,1),ldq)
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else
c                                 a general constraint has entered the working set.
c                                 ifix is not used.
         nanew = nactiv + 1
c                                 transform the incoming row of a by q'.
         call dcopy (n,a(iadd,1),lda,w,1)
         call cmqmul(8,n,nz,nfree,ldq,unitq,kx,w,q,c)
c                                 check that incoming row is not dependent upon those
c                                 already in the working set.
         dtnew = dnrm2 (nz,w,1)

         if (nactiv.eq.0) then
c                                 this is the only general constraint in the working set.
            cond = sdiv (asize,dtnew,overfl)
            tdtmax = dtnew
            tdtmin = dtnew
         else
c                                 there are already general constraints in the working
c                                 set. update the estimate of the condition number.
            tdtmax = max(dtnew,dtmax)
            tdtmin = min(dtnew,dtmin)
            cond = sdiv (tdtmax,tdtmin,overfl)
         end if

         if (cond.gt.condmx .or. overfl) go to 80

         if (unitq) then
c                                 first general constraint added.  set  q = i.
            call smload ('general',nfree,nfree,0d0,1d0,q,ldq)
            unitq = .false.
            it = 0
         end if
      end if

      if (bound) then
         npiv = nfree
      else
         npiv = nz
      end if

      if (unitq) then

c        the orthogonal matrix  q  (i.e.,  q) is not stored explicitly.
c        apply  p, the sequence of pairwise interchanges that moves the
c        newly-fixed variable to position  nfree.

         if (ngq.gt.0) call sgeapr ('left','transpose',nfree-1,w,ngq,
     *                              gqm,n)

         if (rset) then

c           apply the pairwise interchanges to  rz.
c           the subdiagonal elements generated by this process are
c           stored in  s(ifix), s(2), ..., s(nrz-1).

            nsup = nrz - ifix
            call nggnfm (nrz,ifix,nrz,s,r,ldr)
         end if
      else

c        the matrix  q  is stored explicitly.
c        define a sweep of plane rotations p such that
c                           pw = beta*e(npiv).
c        the rotations are applied in the planes (1, 2), (2, 3), ...,
c        (npiv-1, npiv).  the rotations must be applied to q, gqm', r
c        and t.

         call ssrotg('varble','forwrds',npiv-1,w(npiv),w,1,c,s)
c
         if (ngq.gt.0) call sgesrc('left ','variable','forwards',npiv,
     *                             ngq,1,npiv,c,s,gqm,n)
         call sgesrc('right','variable','forwards',nfree,nfree,1,npiv,c,
     *               s,q,ldq)

         if (rset) then

c           apply the rotations to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nrz-1).

            nsup = nrz - 1
            call sutsrh('right',nrz,1,nrz,c,s,r,ldr)
         end if
      end if

      if (rset) then

c        eliminate the  nsup  subdiagonal elements of  r  stored in
c        s(nrz-nsup), ..., s(nrz-1)  with a left-hand sweep of rotations
c        in planes (nrz-nsup, nrz-nsup+1), ..., (nrz-1, nrz).

         call suhqr('left ',nrz,nrz-nsup,nrz,c,s,r,ldr)

         if (nsup.gt.0 .and. drzz.ne.1d0) then
            drzz = c(nrz-1)**2 + drzz*s(nrz-1)**2
         end if
      end if

      if (.not. unitq) then
         if (bound) then

c           bound constraint added.   the rotations affect columns
c           nz+1  thru  nfree  of  gqm'  and  t.

c           the last row and column of  q  has been transformed to plus
c           or minus the unit vector  e(nfree).  we can reconstitute the
c           column of gqm' corresponding to the new fixed variable.
c
            if (w(nfree).lt.0d0) then
               if (ngq.gt.0) call dscal (ngq,-1d0,gqm(nfree,1),n)
            end if

            if (nactiv.gt.0) then
               t(it,jt-1) = s(jt-1)*t(it,jt)
               t(it,jt) = c(jt-1)*t(it,jt)

               if (nactiv.gt.1) then
                  call sutsrh('right',nactiv,1,nactiv,c(jt),s(jt),
     *                        t(it,jt),ldt)
                  call dcopy (nactiv-1,s(jt),1,t(it+1,jt),ldt+1)
               end if

               jt = jt - 1
               call scond (nactiv,t(it,jt),ldt+1,tdtmax,tdtmin)
               cond = sdiv (tdtmax,tdtmin,overfl)
            end if
         else

c           general constraint added.  install  w  at the front of  t.
c           if there is no room,  shift all the rows down one position.

            it = it - 1
            if (it.le.0) then
               it = 1
               do 60 k = 1, nactiv
                  j = jt + k - 1
                  do 40 i = k, 1, -1
                     t(i+1,j) = t(i,j)
   40             continue
   60          continue
            end if
            jt = jt - 1
            call dcopy (nanew,w(jt),1,t(it,jt),ldt)
         end if
      end if
c

c     prepare to exit.  check the magnitude of the condition estimator.

   80 if (nanew.gt.0) then
         if (cond.lt.condmx .and. .not. overfl) then

c           the factorization has been successfully updated.

            inform = 0
            dtmax = tdtmax
            dtmin = tdtmin

         else

c           the proposed working set appears to be linearly dependent.

            inform = 1

         end if
      end if
c                                 end of rzadd
      end


      subroutine npfeas(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,
     *                  nviol,ax,bl,bu,featol,x,work)
c----------------------------------------------------------------------
c     npfeas  computes the following...
c     (1)  the number of constraints that are violated by more
c          than  featol  and the 2-norm of the constraint violations.
c----------------------------------------------------------------------
      implicit none

      double precision bigbnd, cvnorm, errmax
      integer jmax, n, nclin, nviol

      double precision ax(*), bl(n+nclin), bu(n+nclin),
     *                 featol(n+nclin), work(n+nclin), x(n)

      integer istate(n+nclin)

      double precision biglow, bigupp, con, feasj, res, tolj
      integer is, j

      double precision dnrm2
      integer idamax
      external dnrm2, idamax
c----------------------------------------------------------------------
      biglow = -bigbnd
      bigupp = bigbnd

c     compute nviol, the number of constraints violated by more than
c     featol,  and cvnorm,  the 2-norm of the constraint
c     violations and residuals of the constraints in the qp working set.

      nviol = 0

      do j = 1, n + nclin

         feasj = featol(j)
         res = 0d0

            if (j.le.n) then
               con = x(j)
            else
               con = ax(j-n)
            end if
c
            tolj = feasj

c        check for constraint violations.

         if (bl(j).gt.biglow) then
            res = bl(j) - con
            if (res.gt.feasj) nviol = nviol + 1
            if (res.gt.tolj) go to 20
         end if

         if (bu(j).lt.bigupp) then
            res = bu(j) - con
            if (res.lt.(-feasj)) nviol = nviol + 1
            if (res.lt.(-tolj)) go to 20
         end if

c        this constraint is satisfied,  but count the residual as a
c        violation if the constraint is in the working set.

         is = istate(j)

         if (is.eq.0) then
            res = 0d0
         else if (is.eq.1 .or. is.le.-2) then
            res = bl(j) - con
         else if (is.ge.2 .or. is.eq.-1) then
            res = bu(j) - con
         end if

         if (abs(res).gt.feasj) nviol = nviol + 1

c        set the array of violations.

   20    work(j) = res

      end do

      jmax = idamax(n+nclin,work)
      errmax = abs(work(jmax))

      cvnorm = dnrm2 (n+nclin,work,1)
c                                 end of npfeas
      end

      subroutine npsrch (inform,n,objfun,alfa,alfmax,alfsml,
     *                  dxnorm,epsrf,eta,gdx,grdalf,glf,objf,
     *                  objalf,xnorm,dx,grad,gradu,x1,x)
c----------------------------------------------------------------------
c     npsrch finds the steplength alfa that gives sufficient decrease in
c     the augmented lagrangian merit function.
c
c     on exit, if inform = 1, 2 or 3,  alfa will be a nonzero steplength
c     with an associated merit function value  objalf  which is lower
c     than that at the base point. if  inform = 4, 5, 6, 7 or 8,  alfa
c     is zero and  objalf will be the merit value at the base point.

c     set the input parameters and tolerances for srchc  and srchq.

c     tolrx   is the tolerance on relative changes in dx resulting from
c             changes in alfa.
c
c     tolax   is the tolerance on absolute changes in dx resulting from
c             changes in alfa.
c
c     tolabs  is the tolerance on absolute changes in alfa.
c
c     tolrel  is the tolerance on relative changes in alfa.
c
c     toltny  is the magnitude of the smallest allowable value of alfa.
c             if m(tolabs) - m(0).gt.epsaf,  the linesearch tries
c             steps in the range toltny .le. alfa .le. tolabs.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical done, first, imprvd, bad

      integer inform, n, j, maxf, numf

      double precision alfa, alfmax, alfsml, dxnorm, epsrf,
     *                 eta, gdx, glf, grdalf, objalf, objf,
     *                 xnorm, dx(n), grad(n),
     *                 gradu(n), x(n), x1(n), alfbst, epsaf,
     *                 fbest, ftry, g0, gbest, gtry,
     *                 oldf, oldg, q, s, t, targtg, tgdx, tglf,
     *                 tobj, tolabs, tolax, tolrel, tolrx,
     *                 toltny, ddot

      external objfun, ddot

      double precision wmach
      common/ cstmch /wmach(10)

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9
c-----------------------------------------------------------------------
      if (numric) then
         maxf = 15
      else
         maxf = 10
      end if

      epsaf = epsrf*(1d0+abs(objalf))
      tolax = epspt8
      tolrx = epspt8

      if (tolrx*xnorm+tolax.lt.dxnorm*alfmax) then
         tolabs = (tolrx*xnorm+tolax)/dxnorm
      else
         tolabs = alfmax
      end if

      tolrel = max(tolrx,wmach(3))

      t = 0d0

      do j = 1, n
         s = abs(dx(j))
         q = abs(x(j))*tolrx + tolax
         if (s.gt.t*q) t = s/q
      end do
c
      if (t*tolabs.gt.1d0) then
         toltny = 1d0/t
      else
         toltny = tolabs
      end if
c
      oldf = objalf
      oldg = grdalf
      alfbst = 0d0
      fbest = 0d0
      gbest = (1d0 - 1d-4)*oldg
      targtg = (1d-4 - eta)*oldg
      g0 = gbest

      first = .true.

c     commence main loop, entering srchc  or srchq two or more times.
c     first = .true. for the first entry, .false. for subsequent entries
c     done  = .true. indicates termination, in which case the value of
c     inform gives the result of the search.
c     inform = 1 if the search is successful and alfa < alfmax.
c            = 2 if the search is successful and alfa = alfmax.
c            = 3 if a better point was found but too many functions
c                were needed (not sufficient decrease).
c            = 4 if alfmax < tolabs (too small to do a search).
c            = 5 if alfa < alfsml (srchq only -- maybe want to switch
c                to central differences to get a better direction).
c            = 6 if the search found that there is no useful step.
c                the interval of uncertainty is less than 2*tolabs.
c                the minimizer is very close to alfa = zero
c                or the gradients are not sufficiently accurate.
c            = 7 if there were too many function calls.
c            = 8 if the input parameters were bad
c                (alfmax le toltny  or  oldg ge 0).

      do

         if (numric) then
            call srchq (first,done,imprvd,inform,maxf,numf,
     *               alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest)
         else
            call srchc (first,done,imprvd,inform,maxf,numf,
     *               alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest,gbest)
         end if

         if (imprvd) then

            objf = tobj
            objalf = tobj

            if (.not.numric) then

               call dcopy (n,gradu,1,grad,1)
               gdx = tgdx
               glf = tglf

            end if

         end if

         if (.not. done) then
c                                 if ~done, compute objfun for srchc/q
            call dcopy (n,x1,1,x,1)
            call daxpy (n,alfa,dx,1,x,1)
c                                 hack for simplicial composition
            call badalf (alfa,n,x,x1,dx,'a')

            outrpc = .true.

            call objfun (n,x,tobj,gradu,bad)

            outrpc = .false.

            if (bad) then 
               inform = -1
               return
            end if

            ftry = tobj - oldf - 1d-4 * oldg * alfa
c                                 compute auxiliary gradient info 
            if (.not.numric) then

               gtry = ddot (n,gradu,1,dx)
               tgdx = gtry
               tglf = gtry
               gtry = gtry - 1d-4 * oldg

            end if

         else

            exit

         end if

      end do

      alfa = alfbst

      if (.not. imprvd) then

         call dcopy (n,x1,1,x,1)
         call daxpy (n,alfa,dx,1,x,1)
c                                 hack for simplicial composition
         call badalf (alfa,n,x,x1,dx,'b')

      end if
c                                 end of npsrch
      end

      subroutine npupdt (n,ldr,alfa,glf1,glf2,qpcurv,
     *                   gq1,gq2,hpq,rpq,r,wrk1,wrk2)
c----------------------------------------------------------------------
c     npupdt  computes the bfgs update for the approximate hessian of
c     the lagrangian.  if the approximate curvature of the lagrangian
c     function is negative,  a nonnegative penalty vector omega(i) of
c     minimum two norm is computed such that the approximate curvature
c     of the augmented lagrangian will be positive. if no finite penalty
c     vector exists,  the bfgs update is doed with the approximate
c     curvature modified to be a small positive value.
c
c     on entry,  gq1 and gq2 contain the transformed objective gradients
c     at x1 and x2,  hpq contains  r'r(pq), the transformed hessian
c     times the transformed search direction.  the vectors gq1 and hpq
c     are not saved.  if the regular bfgs quasi-newton update could not
c     be done, the first character of lsumry is loaded with 'm'.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character lsumry*5

      logical ssbfgs

      integer ldr, n, j, idamax

      double precision alfa, glf1, glf2, qpcurv,
     *                 hpq(n), r(ldr,*), rpq(n),
     *                 wrk1(n), wrk2(n), gq1(n),
     *                 curvl, eta, rtgtp, rtyts,
     *                 tinycl, trace1, trace2, gq2(n), dnrm2, sdiv

      external dnrm2, sdiv , idamax

      double precision rcndbd, rfrobn, drmax, drmin
      common/ ngg018 /rcndbd, rfrobn, drmax, drmin
c----------------------------------------------------------------------
c     set curvl = (g2 - g1)'dx,  the approximate curvature along dx of
c     the (augmented) lagrangian.  at first, the curvature is not scaled
c     by the steplength alfa.

      curvl = glf2 - glf1
      tinycl = qpcurv*1d-1
      ssbfgs = curvl .le. alfa*tinycl

c     test if curvl is sufficiently positive.  if there are no nonlinear
c     constraints,  no update can be doed.

      if (curvl.lt.tinycl) then 
         lsumry(1:1) = 'modified bfgs'
         curvl = tinycl
      end if

      do j = 1, n
         wrk2(j) = gq2(j) - gq1(j)
      end do
c
      rtgtp = sqrt(qpcurv)
      rtyts = sqrt(alfa*curvl)
      eta = 1d0
      if (ssbfgs) eta = rtyts/(rtgtp*alfa)
c
      trace1 = dnrm2 (n,hpq,1)/rtgtp
      trace2 = dnrm2 (n,wrk2,1)/(rtyts*eta)
      rfrobn = eta*sqrt(abs((rfrobn-trace1)*(rfrobn+trace1)+trace2**2))
c     update the cholesky factor of  q'hq.

c     normalize the vector  rpq (= r(pq)).
c
      call dscal (n,(1d0/rtgtp),rpq,1)
c
c     do the self-scaled or regular bfgs update.
c     form the vector wrk1 = gamma * (gq2 - gq1) - beta * r'r*pq,
c     where  gamma = 1/sqrt(curv) = 1/sqrt((gq2 - gq1)'sq)
c
      call dscal (n,(1d0/rtgtp),hpq,1)
c
      if (ssbfgs) then
         do 140 j = 1, n
            call dscal (j,eta,r(1,j),1)
            wrk1(j) = wrk2(j)/rtyts - eta*hpq(j)
  140    continue
      else
         do 160 j = 1, n
            wrk1(j) = wrk2(j)/rtyts - hpq(j)
  160    continue
      end if

c     do the update to  r = r + rpq*wrk1'.
c     rpq is overwritten. arrays gq1 and hpq are used to store the
c     sines and cosines defined by the plane rotations.

      call cmr1md (n,n,ldr,n,n,r,rpq,wrk1,gq1,hpq)
c                                 end of npupdt
      end

      subroutine cmalf (firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfa,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  anorm,ap,ax,bl,bu,featol,p,x)
c----------------------------------------------------------------------
c     cmalf finds a step alfa such that the point x + alfa*p reaches
c     one of the linear constraints (including bounds).  two possible
c     steps are defined as follows...
c
c     alfa1   is the maximum step that can be taken without violating
c             one of the linear constraints that is currently satisfied.
c     alfa2   reaches a linear constraint that is currently violated.
c             usually this will be the furthest such constraint along p,
c             but if firstv = .true. it will be the first one along p.
c             this is used only when the problem has been determined to
c             be infeasible, and the sum of infeasibilities are being
c             minimized.  (alfa2  is not defined if numinf = 0.)
c
c     alfa will usually be the minimum of alfa1 and alfa2.
c     alfa could be negative (since we allow inactive constraints
c     to be violated by as much as featol).  in such cases, a
c     third possible step is computed, to find the nearest satisfied
c     constraint (perturbed by featol) along the direction  - p.
c     alfa  will be reset to this step if it is shorter.  this is the
c     only case for which the final step  alfa  does not move x exactly
c     onto a constraint (the one denoted by jadd).
c
c     constraints in the working set are ignored  (istate(j) ge 1).
c
c     jadd    denotes which linear constraint is reached.
c
c     hitlow  indicates whether it is the lower or upper bound that
c             has restricted alfa.
c
c     values of istate(j)....
c
c     - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c
c     the values -2 and -1 do not occur once a feasible point has been
c     found.
c----------------------------------------------------------------------
      implicit none

      double precision alfa, atphit, bigalf, bigbnd, palfa, pnorm
      integer inform, jadd, n, nctotl, numinf
      logical firstv, hitlow

      double precision anorm(*), ap(*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), p(n), x(n)
      integer istate(nctotl)

      double precision absatp, alfa1, alfa2, apmax1, apmax2, atp, atp1,
     *                  atp2, atx, palfa1, palfa2, res, rownrm
      integer i, j, jadd1, jadd2, js, jsave1, jsave2
      logical hlow1, hlow2, lastv, negstp, step2

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9
c----------------------------------------------------------------------
      inform = 0

c     first pass -- find steps to perturbed constraints, so that
c     palfa1 will be slightly larger than the true step, and
c     palfa2 will be slightly smaller than it should be.
c     in degenerate cases, this strategy gives us some freedom in the
c     second pass.  the general idea follows that described by p.m.j.
c     harris, p.21 of mathematical programming 5, 1 (1973), 1--28.

      negstp = .false.
      call cmalf1 (firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,palfa1,
     *            palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,featol,p,x)
c
      jsave1 = jadd1
      jsave2 = jadd2
c

c     second pass -- recompute step-lengths without perturbation.
c     amongst constraints that are less than the perturbed steps,
c     choose the one (of each type) that makes the largest angle
c     with the search direction.

      alfa1 = bigalf
      alfa2 = 0d0
      if (firstv) alfa2 = bigalf

      apmax1 = 0d0
      apmax2 = 0d0
      atp1 = 0d0
      atp2 = 0d0
      hlow1 = .false.
      hlow2 = .false.
      lastv = .not. firstv

      do 20 j = 1, nctotl
         js = istate(j)
         if (js.le.0) then
            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               rownrm = 1d0
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               rownrm = anorm(i) + 1d0
            end if

            if (abs(atp).le.epspt9*rownrm*pnorm) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step. 

               res = -1d0

            else if (atp.le.0d0 .and. js.ne.-2) then

c              a'x  is decreasing.

c              the lower bound is satisfied.  test for smaller alfa1.

               absatp = -atp
               if (bl(j).gt.(-bigbnd)) then
                  res = atx - bl(j)
                  if (palfa1*absatp.ge.res .or. j.eq.jsave1) then
                     if (apmax1*rownrm*pnorm.lt.absatp) then
                        apmax1 = absatp/(rownrm*pnorm)
                        alfa1 = res/absatp
                        jadd1 = j
                        atp1 = atp
                        hlow1 = .true.
                     end if
                  end if
               end if

               if (js.eq.-1) then

c                 the upper bound is violated.  test for either a bigger
c                 or smaller alfa2,  depending on the value of firstv.

                  res = atx - bu(j)
                  if ((firstv .and. palfa2*absatp.ge.res .or.
     *                lastv .and. palfa2*absatp.le.res)
     *                .or. j.eq.jsave2) then
                     if (apmax2*rownrm*pnorm.lt.absatp) then
                        apmax2 = absatp/(rownrm*pnorm)
                        if (absatp.ge.1d0) then
                           alfa2 = res/absatp
                        else if (res.lt.bigalf*absatp) then
                           alfa2 = res/absatp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2 = j
                        atp2 = atp
                        hlow2 = .false.
                     end if
                  end if
               end if
            else if (atp.gt.0d0 .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.

c              test for smaller alfa1.

               if (bu(j).lt.bigbnd) then
                  res = bu(j) - atx
                  if (palfa1*atp.ge.res .or. j.eq.jsave1) then
                     if (apmax1*rownrm*pnorm.lt.atp) then
                        apmax1 = atp/(rownrm*pnorm)
                        alfa1 = res/atp
                        jadd1 = j
                        atp1 = atp
                        hlow1 = .false.
                     end if
                  end if
               end if

               if (js.eq.-2) then

c                 the lower bound is violated.  test for a new alfa2.

                  res = bl(j) - atx
                  if ((firstv .and. palfa2*atp.ge.res .or. lastv .and.
     *                palfa2*atp.le.res) .or. j.eq.jsave2) then
                     if (apmax2*rownrm*pnorm.lt.atp) then
                        apmax2 = atp/(rownrm*pnorm)
                        if (atp.ge.1d0) then
                           alfa2 = res/atp
                        else if (res.lt.bigalf*atp) then
                           alfa2 = res/atp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2 = j
                        atp2 = atp
                        hlow2 = .true.
                     end if
                  end if
               end if
            end if

         end if
   20 continue

c     determine alfa, the step to be taken.

c     in the infeasible case, check whether to take the step alfa2
c     rather than alfa1...
c
      step2 = numinf.gt.0 .and. jadd2.gt.0

c     we do so if alfa2 is less than alfa1 or (if firstv is false)
c     lies in the range  (alfa1, palfa1)  and has a smaller value of
c     atp.
c
      step2 = step2 .and. (alfa2.lt.alfa1 .or. lastv .and. alfa2.le.
     *        palfa1 .and. apmax2.ge.apmax1)
c
      if (step2) then
         alfa = alfa2
         palfa = palfa2
         jadd = jadd2
         atphit = atp2
         hitlow = hlow2
      else
         alfa = alfa1
         palfa = palfa1
         jadd = jadd1
         atphit = atp1
         hitlow = hlow1

c        if alfa1 is negative, the constraint to be added (jadd)
c        remains unchanged, but alfa may be shortened to the step
c        to the nearest perturbed satisfied constraint along  - p.

         negstp = alfa.lt.0d0
         if (negstp) then
            call cmalf1(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)

            alfa = -min(abs(alfa),palfa1)
         end if
      end if

c     test for undefined or infinite step.

      if (jadd.eq.0) then
         alfa = bigalf
         palfa = bigalf
      end if

      if (alfa.ge.bigalf) inform = 3
c                                 end of cmalf
      end

      subroutine lssetx(linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,
     *                  n,nctotl,ldzy,lda,ldr,ldt,istate,kactiv,kx,jmax,
     *                  errmax,ctx,xnorm,a,ax,bl,bu,cq,res,res0,featol,
     *                  r,t,x,zy,p,work)
c----------------------------------------------------------------------
c     lssetx  computes the point on a working set that is closest to the
c     input vector  x  (in the least-squares sense).  the norm of  x,
c     the transformed residual vector  pr - rq'x,  and the constraint
c     values ax  are also initialized.

c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable  rowerr
c     is set.
c----------------------------------------------------------------------
      implicit none

      logical linobj, rowerr, unitq

      integer jmax, lda, ldr, ldt, ldzy, n, nactiv, nclin, nctotl,
     *        nfree, istate(nctotl), kactiv(n), kx(n), nrank, nz, i, 
     *        is, j, k, ktry, idamax

      double precision a(lda,*), ax(*), bl(nctotl), bu(nctotl), cq(*),
     *                 featol(nctotl), p(n), r(ldr,*), res(*), res0(*),
     *                 t(ldt,*), work(nctotl), x(n), zy(ldzy,*), bnd,
     *                 ddot, dnrm2, ctx, errmax, xnorm

      external ddot, dnrm2, idamax
c----------------------------------------------------------------------
c     move  x  onto the simple bounds in the working set.

      do 20 k = nfree + 1, n
         j = kx(k)
         is = istate(j)
         bnd = bl(j)
         if (is.ge.2) bnd = bu(j)
         if (is.ne.4) x(j) = bnd
   20 continue

c     move  x  onto the general constraints in the working set.
c     we shall make  ntry  tries at getting acceptable row errors.

      ktry = 1
      jmax = 1
      errmax = 0d0

c     repeat
   40 if (nactiv.gt.0) then

c        set  work = residuals for constraints in the working set.
c        solve for p, the smallest correction to x that gives a point
c        on the constraints in the working set.  define  p = y*(py),
c        where  py  solves the triangular system  t*(py) = residuals.

         do 60 i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(i) = bnd - ddot (n,a(k,1),lda,x)
   60    continue

         call cmtsol(1,ldt,nactiv,t(1,nz+1),work)
         call sload (n,0d0,p,1)
         call dcopy (nactiv,work,1,p(nz+1),1)

         call cmqmul(2,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
         call daxpy (n,1d0,p,1,x,1)
      end if

c     compute the 2-norm of  x.
c     initialize  ax  for all the general constraints.

      xnorm = dnrm2 (n,x,1)
      if (nclin.gt.0) call dgemv ('n',nclin,n,1d0,a,lda,x,0d0,ax)

c     check the row residuals.

      if (nactiv.gt.0) then
         do 80 k = 1, nactiv
            i = kactiv(k)
            j = n + i
            is = istate(j)
            if (is.eq.1) work(k) = bl(j) - ax(i)
            if (is.ge.2) work(k) = bu(j) - ax(i)
   80    continue

         jmax = idamax(nactiv,work)
         errmax = abs(work(jmax))
      end if

      ktry = ktry + 1
c     until    (errmax .le. featol(jmax) .or. ktry.gt.ntry
      if (.not. (errmax.le.featol(jmax) .or. ktry.gt.5)) go to 40
c
      rowerr = errmax.gt.featol(jmax)

c     compute the linear objective value  c'x  and the transformed
c     residual  pr  -  rq'x = res0  -  rq'x.

      if (nrank.gt.0 .or. linobj) then
         call dcopy (n,x,1,p,1)
         call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
      end if

      ctx = 0d0
      if (linobj) ctx = ddot (n,cq,1,p)

      if (nrank.gt.0) then
         call dtrmv ('u','n','n',nrank,r,ldr,p)
         if (nrank.lt.n) call dgemv ('n',nrank,n-nrank,1d0,r(1,nrank+1),
     *                              ldr,p(nrank+1),1d0,p)
         call dcopy (nrank,res0,1,res,1)
         call daxpy (nrank,-1d0,p,1,res,1)
      end if
c                                 end of lssetx
      end

      subroutine lpcore (prbtyp,msg,cset,rset,unitq,iter,
     *                   itmax,jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,
     *                   nz,istate,kactiv,kx,obj,numinf,xnorm,a,
     *                   ax,bl,bu,cvec,featol,featlu,x,w)
c----------------------------------------------------------------------
c     lpcore  is a subroutine for linear programming.
c     on entry, it is assumed that an initial working set of
c     linear constraints and bounds is available.  the arrays  istate,
c     kactiv  and  kx  will have been set accordingly
c     and the arrays  t  and  q  will contain the tq factorization of
c     the matrix whose rows are the gradients of the active linear
c     constraints with the columns corresponding to the active bounds
c     removed.  the tq factorization of the resulting (nactiv by nfree)
c     matrix is  a(free)*q = (0 t),  where q is (nfree by nfree) and t
c     is upper-triangular.

c     over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.

c     after kdegen consecutive iterations, variables within featol of
c     their bounds are set exactly on their bounds and iterative
c     refinement is used to satisfy the constraints in the working set.
c     featol is then reduced to tolx0 for the next cycle of iterations.

c     values of istate(j) for the linear constraints.......

c     istate(j)
c     ---------
c          0    constraint j is not in the working set.
c          1    constraint j is in the working set at its lower bound.
c          2    constraint j is in the working set at its upper bound.
c          3    constraint j is in the working set as an equality.

c     constraint j may be violated by as much as featol(j).
c----------------------------------------------------------------------
      implicit none 

      character empty*6
      parameter (empty='      ')

      double precision obj, xnorm
      integer iter, itmax, jinf, lda, n, nactiv, nclin, nfree,
     *                  nrz, numinf, nviol, nz
      logical cset, rset, unitq
      character*2       prbtyp
      character*6       msg

      double precision a(lda,*), ax(*), bl(n+nclin), bu(n+nclin),
     *                  cvec(*), featlu(n+nclin), featol(n+nclin), w(*),
     *                  x(n)
      integer istate(n+nclin), kactiv(n), kx(n)

      double precision alfap, alfhit, bigalf, biggst, condmx, 
     *                  condt, dinky, dnorm, dzz, errmax, flmax,
     *                  gfnorm, grznrm, gznorm, objsiz, rtmax, smllst,
     *                  suminf, tinyst, trubig, trusml, wssize, zerolm

      integer iadd, ifix, inform, is, it, j, jbigst,
     *                  jmax, jsmlst, jtiny, kbigst, kdel, ksmlst,
     *                  ldr, ltmp,
     *                  nctotl, nfixed, ngq, nmoved, notopt, ntfixd

      integer         lkactv,lkx,lfeatu,lanorm,lad,ld,lgq,lcq,lrlam,
     *                lr, lt, lq, lwtinf, lwrk

      common/ cstlcl /lkactv,lkx,lfeatu,lanorm,lad,ld,lgq,lcq,lrlam,
     *                lr, lt, lq, lwtinf, lwrk

      logical firstv, fp, hitlow, lp, move, onbnd, overfl,
     *                  unbndd

      double precision ddot, dnrm2, sdiv 
      external ddot, dnrm2, sdiv 

      double precision wmach
      common/ cstmch /wmach(10)

      integer ldt, ldq, ncolt
      common/ ngg004 /ldt, ncolt, ldq

      integer itnfix, kdegen, ndegen, nfix
      double precision tolinc, tolx0
      common/ ngg005 /tolx0, tolinc, kdegen, ndegen, itnfix, nfix(2)

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      logical prnt
      integer isdel, jadd, jdel
      double precision alfa, trulam
      common/ ngg007 /alfa, trulam, isdel, jdel, jadd, prnt

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin

      integer itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *        maxact, mxfree, maxnz
      common/ ngg010 /itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *                maxact, mxfree, maxnz

      double precision bigbnd, bigdx, bndlow, bndupp, tolact, tolfea, 
     *                 tolrnk
      common/ ngg011 /bigbnd, bigdx, bndlow, bndupp,
     *                tolact, tolfea, tolrnk

      save firstv
c----------------------------------------------------------------------
c     specify the machine-dependent parameters.

      flmax = wmach(7)
      rtmax = wmach(8)

      if (cset) then
         ngq = 2
      else
         ngq = 1
      end if

      lp = prbtyp.eq.'lp'
      fp = .not. lp

      ldr = ldt
      it = 1

c     we need a temporary array when changing the active set.
c     use the multiplier array.

      ltmp = lrlam

      if (iter.eq.0) then

c        first entry.  initialize.

         jadd = 0
         jdel = 0
         isdel = 0
         firstv = .false.

         alfa = 0d0
         dzz = 1d0
      end if

      nctotl = n + nclin
      nviol = 0

      condmx = flmax

      call cmsinf (n,nclin,lda,istate,bigbnd,numinf,suminf,bl,bu,a,
     *             featol,w(lgq),x,w(lwtinf))

      if (numinf.gt.0) then
         call cmqmul(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),w(lwrk))
      else if (lp) then
         call dcopy (n,w(lcq),1,w(lgq),1)
      end if

      if (numinf.eq.0 .and. lp) then
         obj = ddot (n,cvec,1,x)
      else
         obj = suminf
      end if

      msg = empty

c*    ======================start of main loop==========================
c     +    do while (msg.eq.empty)
   20 if (msg.eq.empty) then

         gznorm = 0d0
         if (nz.gt.0) gznorm = dnrm2 (nz,w(lgq),1)

         if (nrz.eq.nz) then
            grznrm = gznorm
         else
            grznrm = 0d0
            if (nrz.gt.0) grznrm = dnrm2 (nrz,w(lgq),1)
         end if

         gfnorm = gznorm
         if (nfree.gt.0 .and. nactiv.gt.0) gfnorm = dnrm2 (nfree,w(lgq),
     *       1)

c        define small quantities that reflect the size of x, r and
c        the constraints in the working set.

         if (prnt) then
            condt = 1d0
            if (nactiv.gt.0) condt = sdiv (dtmax,dtmin,overfl)
            jdel = 0
            jadd = 0
            alfa = 0d0
         end if

         if (numinf.gt.0) then
            dinky = epspt8*abs(suminf)
         else
            objsiz = 1d0 + abs(obj)
            wssize = 0d0
            if (nactiv.gt.0) wssize = dtmax
            dinky = epspt8*max(wssize,objsiz,gfnorm)
         end if

c        if the reduced gradient z'g is small enough,
c        lagrange multipliers will be computed.

         if (numinf.eq.0 .and. fp) then
            msg = 'feasbl'
            nfixed = n - nfree
            call sload (nactiv+nfixed,0d0,w(lrlam),1)
            go to 20
         end if

         if (grznrm.le.dinky) then

c           the point  x  is a constrained stationary point.
c           compute lagrange multipliers.

c           define what we mean by 'tiny' and non-optimal multipliers.

            notopt = 0
            jdel = 0
            zerolm = -dinky
            smllst = -dinky
            biggst = dinky + 1d0
            tinyst = dinky

            call cmmul1 (n,lda,ldt,nactiv,nfree,nz,istate,
     *                   kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                   jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                   jbigst,kbigst,a,w(lanorm),w(lgq),w(lrlam),
     *                   w(lt),w(lwtinf))

            if (nrz.lt.nz) call cmmul2(n,nrz,nz,zerolm,notopt,
     *                                 numinf,trusml,smllst,jsmlst,
     *                                 tinyst,jtiny,w(lgq))

            if (abs(jsmlst).gt.0) then

c              delete a constraint.
c              cmmul1  or  cmmul2  found a non-optimal multiplier.

               trulam = trusml
               jdel = jsmlst

               if (jsmlst.gt.0) then

c                 regular constraint.

                  kdel = ksmlst
                  isdel = istate(jdel)
                  istate(jdel) = 0
               end if
            else
               if (numinf.gt.0 .and. jbigst.gt.0) then

c                 no feasible point exists for the constraints but
c                 the sum of the constraint violations can be reduced
c                 by moving off constraints with multipliers greater
c                 than 1.

                  jdel = jbigst
                  kdel = kbigst
                  isdel = istate(jdel)
                  if (trubig.le.0d0) is = -1
                  if (trubig.gt.0d0) is = -2
                  istate(jdel) = is
                  trulam = trubig
                  firstv = .true.
                  numinf = numinf + 1
               end if
            end if

            if (jdel.eq.0) then
               if (numinf.gt.0) then
                  msg = 'infeas'
               else
                  msg = 'optiml'
               end if
               go to 20
            else 
c debug debug
               if (jdel.gt.0.and.nfree.eq.ldq) then 
                  msg = 'infeas'
                  goto 20
               end if

            end if

c           constraint  jdel  has been deleted.
c           update the  tq  factorization.

            call rzdel (unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,w(lt),w(lgq),w(lq),
     *                  w(ld),w(lrlam))
            if (rset) call lpcolr(nrz,ldr,w(lr),1d0)

            prnt = .false.
         else

c           compute a search direction.

            if (iter.ge.itmax) then
               msg = 'itnlim'
               go to 20
            end if

            prnt = .true.
            iter = iter + 1

            call dcopy (nrz,w(lgq),1,w(ld),1)
            call dscal (nrz,-1d0,w(ld),1)

            dnorm = dnrm2 (nrz,w(ld),1)

            call cmqmul(1,n,nrz,nfree,ldq,unitq,kx,w(ld),w(lq),w(lwrk))
            call dgemv ('no transpose',nclin,n,1d0,a,lda,w(ld),0d0,
     *                 w(lad))

c           find the constraint we bump into along d.
c           update  x  and  ax  if the step alfa is nonzero.

c           alfhit is initialized to bigalf. if it remains that value
c           after the call to  cmchzr, it is regarded as infinite.

            bigalf = sdiv (bigdx,dnorm,overfl)

            call cmchzr(firstv,n,nclin,istate,bigalf,bigbnd,dnorm,
     *                  hitlow,move,onbnd,unbndd,alfhit,alfap,jadd,
     *                  w(lanorm),w(lad),ax,bl,bu,featol,featlu,w(ld),x)

            if (unbndd) then
               msg = 'unbndd'
               go to 20
            end if

            alfa = alfhit
            call daxpy (n,alfa,w(ld),1,x,1)

            if (nclin.gt.0) call daxpy (nclin,alfa,w(lad),1,ax,1)
            xnorm = dnrm2 (n,x,1)

c           add a constraint to the working set.
c           update the  tq  factors of the working set.
c           use  d  as temporary work space.

            if (bl(jadd).eq.bu(jadd)) then
               istate(jadd) = 3
            else if (hitlow) then
               istate(jadd) = 1
            else
               istate(jadd) = 2
            end if

            if (jadd.gt.n) then
               iadd = jadd - n
            else
               if (alfa.ge.0d0) then
                  if (hitlow) then
                     x(jadd) = bl(jadd)
                  else
                     x(jadd) = bu(jadd)
                  end if
               end if
               do 40 ifix = 1, nfree
                  if (kx(ifix).eq.jadd) go to 60
   40          continue
   60       end if

            call rzadd (unitq,rset,inform,ifix,iadd,jadd,it,nactiv,nz,
     *                  nfree,nrz,ngq,n,lda,ldq,ldr,ldt,kx,condmx,dzz,a,
     *                  w(lr),w(lt),w(lgq),w(lq),w(lwrk),w(lrlam),w(ld))

            nz = nz - 1
            nrz = nrz - 1

            if (jadd.le.n) then

c                                 a simple bound has been added.
               nfree = nfree - 1
            else
c                                a general constraint has been added.

               nactiv = nactiv + 1
               kactiv(nactiv) = iadd
            end if

c           increment featol.

            call daxpy (nctotl,tolinc,featlu,1,featol,1)

            if (mod(iter,kchk).eq.0) then

c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  set inform to force iterative
c              refinement and a switch to phase 1.

               call cmfeas(n,nclin,istate,bigbnd,nviol,jmax,errmax,ax,
     *                     bl,bu,featol,x)

            end if

            if (mod(iter,kdegen).eq.0) then

c              every  kdegen  iterations, reset  featol  and
c              move  x  on to the working set if it is close.

               call cmdgen('end of cycle',n,nclin,nmoved,iter,
     *                     numinf,istate,bl,bu,featol,featlu,x)

               nviol = nviol + nmoved
            end if

            if (nviol.gt.0) then
               msg = 'resetx'
               go to 20
            end if

            if (numinf.ne.0) then
               call cmsinf(n,nclin,lda,istate,bigbnd,numinf,suminf,bl,
     *                     bu,a,featol,w(lgq),x,w(lwtinf))

               if (numinf.gt.0) then
                  call cmqmul(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),
     *                        w(lwrk))
               else if (lp) then
                  call dcopy (n,w(lcq),1,w(lgq),1)
               end if
            end if

            if (numinf.eq.0 .and. lp) then
               obj = ddot (n,cvec,1,x)
            else
               obj = suminf
            end if
         end if
         go to 20
c        +    end while
      end if
c     ======================end of main loop============================

      if (msg.eq.'optiml') then
         if (lp) then
            if (nrz.lt.nz) then
               msg = 'weak  '
            else
               ntfixd = 0
               do 80 j = 1, n
                  if (istate(j).eq.4) ntfixd = ntfixd + 1
   80          continue
               if (ntfixd.gt.0) msg = 'weak  '
            end if
            if (abs(jtiny).gt.0) msg = 'weak  '
         end if
      else if (msg.eq.'unbndd' .and. numinf.gt.0) then
         msg = 'infeas'
      end if
c                                 end of lpcore
      end

      subroutine npcore (unitq,inform,majits,n,nclin,
     *                  nctotl,nactiv,nfree,nz,ldaqp,ldr,
     *                  istate,kactiv,kx,objf,fdnorm,xnorm,
     *                  objfun,aqp,ax,bl,bu,clamda,
     *                  featol,grad,gradu,r,x,iw,w,lenw)
c----------------------------------------------------------------------
c     npcore  is the core routine for nlpopt, a sequential quadratic
c     programming (sqp) method for nonlinearly constrained optimization.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical unitq, convpt, convrg, done, error, feasqp,
     *        goodgq, infeas, newgq, optiml, overfl, bad

      integer inform, ldaqp, ldr, lenw, majits, n, nactiv, nclin, 
     *        nctotl, nfree, nz, istate(*), iw(*), 
     *        kactiv(n), kx(n), isum, info, jmax, linact, majit0, 
     *        nqpinf, numinf, nviol,
     *        nlnact, nlserr, nmajor, nminor, nqperr


      double precision aqp(ldaqp,*), ax(*), bl(nctotl), bu(nctotl),
     *                 clamda(nctotl), featol(nctotl), grad(n), xnorm,
     *                 gradu(n), r(ldr,*), w(lenw), x(n), fdnorm, objf,
     *                 alfa, alfbnd, alfdx, alflim, alfmax, alfmin,
     *                 alfsml, cond, condh, condhz, condt, ddot, obj,
     *                 cvnorm, dinky, drzmax, drzmin, dxnorm, errmax,
     *                 flmax, gdx, gfnorm, glf1, glf2, glnorm, gltest,
     *                 grdalf, gtest, gznorm, objalf, objsiz,
     *                 qpcurv, rootn, rtftol, rtmax, xsize, dnrm2, sdiv 

      external ddot, dnrm2, sdiv, objfun

      integer lkactv, lanorm, lcjdxx, lqpdx, lrpq, lqrwrk, lqphz,
     *        lhpq, lgq, lrlam, lt, lzy, lwtinf, lwrk1, lqptol
      common/ cstlnp /lkactv, lanorm, lcjdxx, lqpdx, lrpq, lqrwrk,lqphz,
     *                lhpq, lgq, lrlam, lt, lzy, lwtinf, lwrk1, lqptol

      integer liperm , ladx, lbl, lbu, ldx, lgq1, lx1, lwrk2, lhctrl
      common/ cstln2 /liperm, ladx, lbl, lbu, ldx, lgq1, lx1,
     *        lwrk2, lhctrl

      double precision wmach
      common/ cstmch /wmach(10)

      integer ldt, ncolt, ldzy
      common/ ngg004 /ldt, ncolt, ldzy

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin

      integer itmax1, itmax2, lprob
      common/ ngg016 /itmax1, itmax2, lprob

      double precision rcndbd, rfrobn, drmax, drmin
      common/ ngg018 /rcndbd, rfrobn, drmax, drmin

      double precision cdint, ctol, dxlim, epsrf, eta, fdint, ftol,
     *                 hcndbd
      common/ ngg021 /cdint, ctol, dxlim, epsrf, eta,
     *                fdint, ftol, hcndbd

      equivalence (itmxnp,nmajor), (itmax2,nminor)
c----------------------------------------------------------------------
c                                 specify machine-dependent parameters.
      flmax = wmach(7)
      rtmax = wmach(8)
c                                 initialize
      nqpinf = 0
      majit0 = majits

      alfa = 0d0
      alfdx = 0d0
      rtftol = sqrt(ftol)
      rootn = sqrt(dble(n))
c                                 information from the feasibility phase will be used to generate a
c                                 hot start for the first qp subproblem.
      call dcopy (nctotl,featol,1,w(lqptol),1)

      objalf = objf
      newgq = .false.
      isum = 0

      do
c                                 loop to find good gradient
c        if (rids.eq.3) then
c           write (*,*) isum, objf
c           write (*,*) x
c        end if

         do
c                                 loop to follow the gradient
            if (newgq) then

               if (numric) then
c                                 in principle no call to objfun is necessary
c                                 because switching to 2nd order gradient
c              call objfun (n,x,obj,gradu)
c              if (obj.ne.objf) then
c                 write (*,*) 'wtf',objf-obj
c              end if 
c                                 compute derivatives
                  call numder (objf,objfun,grad,x,fdnorm,bl,bu,n,bad)

                  if (bad) then
                     inform = -1
                     return
                  end if

               end if

               w(lgq:lgq+n-1) = grad(1:n)

               call cmqmul (6,n,nz,nfree,ldzy,unitq,kx,w(lgq),w(lzy),
     *                      w(lwrk1))

               newgq = .false.

            end if
c                                 (1) solve iqp search direction and multiplier estimates.
c                                 (2) compute the search direction for the slack variables
c                                 and multipliers.
            call npiqp (feasqp,unitq,nqperr,n,nclin,
     *            ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,numinf,istate,
     *            kactiv,kx,dxnorm,gdx,qpcurv,aqp,w(ladx),ax,
     *            bl,bu,clamda,w(ldx),w(lbl),w(lbu),w(lqptol),r,
     *            x,w(lwtinf),w)

            if (feasqp) then
               nqpinf = 0
            else
               nqpinf = nqpinf + 1
            end if
c                                 convergence test:
c                                 compute norms of projected gradient and
c                                 gradient of the free variables.
            gznorm = 0d0
            if (nz.gt.0) gznorm = dnrm2 (nz,w(lgq),1)

            gfnorm = gznorm
            if (nfree.gt.0.and.nactiv.gt.0) 
     *         gfnorm = dnrm2 (nfree,w(lgq),1)
c                                 if the forward-difference estimates
c                                 are small, switch to central differences
c                                 and re-solve the qp.
            goodgq = .true.

            if (numric .and. .not. cntrl) then

               glnorm = dnrm2 (n,w(lhpq),1)
               gltest = (1d0+abs(objf))*epsrf/fdnorm

               if (glnorm.le.gltest) then
c                                 up the ante:
                  goodgq = .false.
                  cntrl = .true.
                  newgq = .true.

               end if

            end if

            if (goodgq) exit
c                                 end of inner loop
         end do
c                                 count violated constraints (> featol)
c                                 get 2-norm of the working set constraint residuals
         call npfeas(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *            ax,bl,bu,featol,x,w(lwrk2))
c                                 define small quantities that reflect the magnitude of objf and
c                                 the norm of grad(free).
         objsiz = 1d0 + abs(objf)

         xsize = 1d0 + xnorm
         gtest = max(objsiz,gfnorm)
         dinky = rtftol*gtest

         if (nactiv.eq.0) then
            condt = 0d0
         else if (nactiv.eq.1) then
            condt = dtmin
         else
            condt = sdiv (dtmax,dtmin,overfl)
         end if

         call scond (n,r,ldr+1,drmax,drmin)

         condh = sdiv (drmax,drmin,overfl)
         if (condh.lt.rtmax) then
            condh = condh*condh
         else
            condh = flmax
         end if

         if (nz.eq.0) then
            condhz = 1d0
         else if (nz.eq.n) then
            condhz = condh
         else
            call scond (nz,r,ldr+1,drzmax,drzmin)
            condhz = sdiv (drzmax,drzmin,overfl)
            if (condhz.lt.rtmax) then
               condhz = condhz*condhz
            else
               condhz = flmax
            end if
         end if
c                                 test for convergence.
c                                 convpt checks for a k-t point at the initial
c                                 point or after a large change in x.
         convpt = dxnorm .le. epspt8*gtest .and. nviol.eq.0 .and.
     *         nqperr .le. 1

         optiml = gznorm.lt.dinky .and. nviol.eq.0

         convrg = majits.gt.0 .and. alfdx .le. rtftol*xsize

         infeas = convrg .and. .not. feasqp .or. nqpinf.gt.7

         done = convpt .or. (convrg .and. optiml) .or. infeas

         objalf = objf
         grdalf = gdx
         glf1 = gdx

         alfa = 0d0
         error = majits .ge. nmajor

         if (.not. (done .or. error)) then

            majits = majits + 1
c                                 copy information needed for the bfgs update.
            call dcopy (n,x,1,w(lx1),1)
            call dcopy (n,w(lgq),1,w(lgq1),1)
c                                 compute parameters for \linesearch.
c                                 alfmin is the smallest allowable step predicted by the qp
c                                 subproblem.
            alfmin = 1d0
            if (.not. feasqp) alfmin = 0d0
c                                 alfmax is the largest feasible steplength subject to 
c                                 input specified alflim.
            alfmax = sdiv (bigdx,dxnorm,overfl)

            call npalf(info,n,nclin,alfa,alfmin,alfmax,bigbnd,
     *                 dxnorm,w(lanorm),w(ladx),ax,bl,bu,w(ldx),x)

            alfmax = alfa
            if (alfmax.lt.1d0+epspt3 .and. feasqp) alfmax = 1d0
c                                 alfbnd is a tentative upper bound on the steplength. if the
c                                 merit function is decreasing at alfbnd and certain
c                                 conditions hold, alfbnd will be increased in multiples of
c                                 two (subject to not being greater than alfmax).
            alfbnd = alfmax
c                                 alfsml trips computation of central differences. if a
c                                 trial steplength falls below alfsml, the linesearch is
c                                 terminated.
            alfsml = 0d0
            if (numric .and. cntrl) then
               alfsml = sdiv (fdnorm,dxnorm,overfl)
               alfsml = min(alfsml,alfmax)
            end if
c                                 compute the steplength
            alflim = sdiv ((1d0+xnorm)*dxlim,dxnorm,overfl)
            alfa = min(alflim,1d0)

            call npsrch (nlserr,n,objfun,alfa,alfmax,alfsml,
     *                   dxnorm,epsrf,eta,gdx,grdalf,glf2,objf,objalf,
     *                   xnorm,w(ldx),grad,gradu,w(lx1),x)

c                                 npsrch  sets nlserr to the following values...
c                                -1  objective function failure
c                                 1  if the search is successful and alfa < alfmax.
c                                 2  if the search is successful and alfa = alfmax.
c                                 3  if a better point was found but too many functions
c                                    were needed (not sufficient decrease).

c                                 values of nlserr occurring with a nonzero value of alfa.
c                                 4  if alfmax < tolabs (too small to do a search).
c                                 5  if alfa  < alfsml (srchq only -- maybe want to switch
c                                    to central differences to get a better direction).
c                                 6  if the search found that there is no useful step.
c                                    the interval of uncertainty is less than 2*tolabs.
c                                    the minimizer is very close to alfa = zero
c                                    or the gradients are not sufficiently accurate.
c                                 7  if there were too many function calls.
c                                 8  if the input parameters were bad
c                                    (alfmax le toltny or uphill).
            if (nlserr.lt.0) then
               inform = -1
               return
            end if

            error = nlserr .ge. 4

            if (error) then
c                                 linesearch failed to find a better point.
c                                 if exact or central difference gradients
c                                 are in use or the kt conditions are satisfied, 
c                                 stop.  otherwise, switch to central differences 
c                                 and solve the qp again.
               if (numric .and. cntrl) then

                  if (.not. optiml) then

                     error = .false.

                     newgq = .true.

                  end if

               end if

            else

               if (numric) then
c                                 compute the missing gradients.
c              call objfun (n,x,obj,gradu)
c              if (obj.ne.objf) then
c                 write (*,*) 'wtf 2',objf-obj
c              end if 

                  call numder (objf,objfun,grad,x,fdnorm,bl,bu,n,bad)

                  if (bad) then
                     inform = -1
                     return
                  end if

                  gdx = ddot (n,grad,1,w(ldx))
                  glf2 = gdx

               end if

               call dcopy (n,grad,1,w(lgq),1)

               call cmqmul (6,n,nz,nfree,ldzy,unitq,kx,w(lgq),w(lzy),
     *                      w(lwrk1))

               xnorm = dnrm2 (n,x,1)

               if (nclin.gt.0) call daxpy (nclin,alfa,w(ladx),1,ax,1)
               alfdx = alfa*dxnorm
c                                 update the factors of the approximate hessian of the
c                                 lagrangian function.
               call npupdt (n,ldr,alfa,glf1,glf2,qpcurv,w(lgq1),
     *                      w(lgq),w(lhpq),w(lrpq),r,w(lwrk2),w(lwrk1))

               call scond (n,r,ldr+1,drmax,drmin)
               cond = sdiv (drmax,drmin,overfl)

               if (cond.gt.rcndbd .or.rfrobn.gt.rootn*1d2*drmax) then
c                                 reset the condition estimator and range-space
c                                 partition of q'hq.
                  call nprset (unitq,n,nfree,nz,ldzy,ldr,iw(liperm),kx,
     *                     w(lgq),r,w(lzy),w(lwrk1),w(lqrwrk))

               end if
            end if
         end if

         if (done .or. error) exit

      end do

      if (done) then
         if (convrg .and. optiml) then
            inform = 0
         else if (convpt) then
            inform = 1
         else if (infeas) then
            inform = 3
         end if
      else
         if (majits.ge.nmajor) then
            inform = 4
         else if (optiml) then
            inform = 1
         else
            inform = 6
         end if
      end if
c                                 used to set clamda here, but no point
c                                 end of npcore
      end

      subroutine cmsetx (rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,
     *                  ldt,istate,kactiv,kx,jmax,errmax,xnorm,a,ax,bl,
     *                  bu,featol,t,x,q,p,work)
c----------------------------------------------------------------------
c     cmsetx  computes the point on a working set that is closest in the
c     least-squares sense to the input vector x.
c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable rowerr
c     is set.
c----------------------------------------------------------------------
      implicit none

      double precision errmax, xnorm
      integer jmax, lda, ldq, ldt, n, nactiv, nclin, nfree, nz
      logical rowerr, unitq

      double precision a(lda,*), ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), p(n), q(ldq,*), t(ldt,*),
     *                  work(n), x(n)
      integer istate(n+nclin), kactiv(n), kx(n)

      double precision bnd
      integer i, is, j, k, ktry

      double precision ddot, dnrm2
      integer idamax
      external ddot, dnrm2, idamax
c----------------------------------------------------------------------
c                                 move  x  onto the simple bounds in the working set.
      do 20 k = nfree + 1, n
         j = kx(k)
         is = istate(j)
         bnd = bl(j)
         if (is.ge.2) bnd = bu(j)
         if (is.ne.4) x(j) = bnd
   20 continue
c                                 move  x  onto the general constraints in the working set.
c                                 ntry  attempts are made to get acceptable row errors.
      ktry = 1
      jmax = 1
      errmax = 0d0

c     repeat
   40 if (nactiv.gt.0) then

c        set work = residuals for constraints in the working set.
c        solve for p, the smallest correction to x that gives a point
c        on the constraints in the working set.  define  p = y*(py),
c        where  py  solves the triangular system  t*(py) = residuals.

         do 60 i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(nactiv-i+1) = bnd - ddot (n,a(k,1),lda,x)
   60    continue

         call dtrsv ('u','n','n',nactiv,t(1,nz+1),ldt,work)
         call sload (n,0d0,p,1)
         call dcopy (nactiv,work,1,p(nz+1),1)

         call cmqmul(2,n,nz,nfree,ldq,unitq,kx,p,q,work)
         call daxpy (n,1d0,p,1,x,1)
      end if

c     compute the 2-norm of  x.
c     initialize  ax  for all the general constraints.

      xnorm = dnrm2 (n,x,1)
      if (nclin.gt.0) call dgemv ('n',nclin,n,1d0,a,lda,x,0d0,ax)

c     check the row residuals.

      if (nactiv.gt.0) then
         do 80 k = 1, nactiv
            i = kactiv(k)
            j = n + i
            is = istate(j)
            if (is.eq.1) work(k) = bl(j) - ax(i)
            if (is.ge.2) work(k) = bu(j) - ax(i)
   80    continue

         jmax = idamax(nactiv,work)
         errmax = abs(work(jmax))
      end if

      ktry = ktry + 1
c     until    (errmax .le. featol(jmax) .or. ktry.gt.ntry
      if (.not. (errmax.le.featol(jmax) .or. ktry.gt.5)) go to 40

      rowerr = errmax.gt.featol(jmax)
c                                 end of cmsetx
      end

      subroutine rzadds (unitq,vertex,k2,it,nactiv,nartif,nz,nfree,
     *                   nrejtd,ngq,n,ldq,lda,ldt,istate,kactiv,kx,
     *                   condmx,a,t,gqm,q,w,c,s)
c----------------------------------------------------------------------
c     rzadds includes general constraints 1:k2 as new rows of
c     the tq factorization:
c              a(free) * q(free)  = ( 0 t)
c                        q(free)  = ( z y)
c     a) the  nactiv x nactiv  upper-triangular matrix  t  is stored
c        with its (1,1) element in position  (it,jt)  of the array  t.
c----------------------------------------------------------------------
      implicit none

      logical unitq, vertex, overfl, rset

      integer it, k2, lda, ldq, ldt, n, nactiv, kx(n), k, l,
     *        nartif, nfree, ngq, nrejtd, nz, istate(*), kactiv(n),
     *        i, iadd, iartif, ifix, inform, iswap, j, jadd, jt, nzadd

      double precision a(lda,*), c(n), gqm(n,*), q(ldq,*), s(n), dnrm2, 
     *                 t(ldt,*), w(n), condmx, cndmax, cond, delta,sdiv,
     *                 drzz, dtnew, rnorm, rowmax, rtmax, tdtmax, tdtmin

      external dnrm2, sdiv 

      double precision wmach
      common/ cstmch /wmach(10)

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin
c----------------------------------------------------------------------
      rtmax = wmach(8)

      jt = nz + 1

c     estimate the condition number of the constraints already
c     factorized.

      if (nactiv.eq.0) then
         dtmax = 0d0
         dtmin = 1d0
         if (unitq) then

c           first general constraint added.  set  q = i.

            call smload ('general',nfree,nfree,0d0,1d0,q,ldq)
            unitq = .false.
         end if
      else
         call scond (nactiv,t(it,jt),ldt+1,dtmax,dtmin)
      end if

      do 20 k = 1, k2
         iadd = kactiv(k)
         jadd = n + iadd
         if (nactiv.lt.nfree) then

            overfl = .false.

c           transform the incoming row of  a  by  q'.

            call dcopy (n,a(iadd,1),lda,w,1)
            call cmqmul(8,n,nz,nfree,ldq,unitq,kx,w,q,s)

c           check that the incoming row is not dependent upon those
c           already in the working set.

            dtnew = dnrm2 (nz,w,1)
            if (nactiv.eq.0) then

c              this is the first general constraint in the working set.

               cond = sdiv (asize,dtnew,overfl)
               tdtmax = dtnew
               tdtmin = dtnew
            else

c              there are already some general constraints in the working
c              set. update the estimate of the condition number.

               tdtmax = max(dtnew,dtmax)
               tdtmin = min(dtnew,dtmin)
               cond = sdiv (tdtmax,tdtmin,overfl)
            end if

            if (cond.ge.condmx .or. overfl) then

c              this constraint appears to be dependent on those already
c              in the working set.  skip it.

               istate(jadd) = 0
               kactiv(k) = -kactiv(k)
            else
               if (nz.gt.1) then

                  delta = w(nz)
                  call sgrfg (nz-1,delta,w,1,0d0,w(nz))
                  if (w(nz).gt.0d0) then

                     call dgemv ('n',nfree,nz,1d0,q,ldq,w,0d0,s)
                     call dger (nfree,nz,-1d0,s,w,q,ldq)

                     if (ngq.gt.0) then
                        call dgemv ('t',nz,ngq,1d0,gqm,n,w,0d0,s)
                        call dger (nz,ngq,-1d0,w,s,gqm,n)
                     end if
                  end if

                  w(nz) = delta
               end if
               it = it - 1
               jt = jt - 1
               nactiv = nactiv + 1
               nz = nz - 1
               call dcopy (nactiv,w(jt),1,t(it,jt),ldt)
               dtmax = tdtmax
               dtmin = tdtmin
            end if
         end if
   20 continue

      if (nactiv.lt.k2) then

c        some of the constraints were classed as dependent and not
c        included in the factorization.  re-order the part of  kactiv
c        that holds the indices of the general constraints in the
c        working set.  move accepted indices to the front and shift
c        rejected indices (with negative values) to the end.

         l = 0
         do 40 k = 1, k2
            i = kactiv(k)
            if (i.ge.0) then
               l = l + 1
               if (l.ne.k) then
                  iswap = kactiv(l)
                  kactiv(l) = i
                  kactiv(k) = iswap
               end if
            end if
   40    continue

c        if a vertex is required,  add some temporary bounds.
c        we must accept the resulting condition number of the working
c        set.

         if (vertex) then
            rset = .false.
            cndmax = rtmax
            drzz = 1d0
            nzadd = nz
            do 80 iartif = 1, nzadd
               if (unitq) then
                  ifix = nfree
                  jadd = kx(ifix)
               else
                  rowmax = 0d0
                  do 60 i = 1, nfree
                     rnorm = dnrm2 (nz,q(i,1),ldq)
                     if (rowmax.lt.rnorm) then
                        rowmax = rnorm
                        ifix = i
                     end if
   60             continue
                  jadd = kx(ifix)

                  call rzadd (unitq,rset,inform,ifix,iadd,jadd,it,
     *                        nactiv,nz,nfree,nz,ngq,n,lda,ldq,ldt,ldt,
     *                        kx,cndmax,drzz,a,t,t,gqm,q,w,c,s)
               end if
               nfree = nfree - 1
               nz = nz - 1
               nartif = nartif + 1
               istate(jadd) = 4
   80       continue
         end if

         if (it.gt.1) then

c           if some dependent constraints were rejected,  move  t  to
c           the top of the array  t.

            do 120 k = 1, nactiv
               j = nz + k
               do 100 i = 1, k
                  t(i,j) = t(it+i-1,j)
  100          continue
  120       continue
         end if
      end if

      nrejtd = k2 - nactiv
c                                 end of rzadds
      end

      subroutine cmqmul(mode,n,nz,nfree,nq,unitq,kx,v,zy,wrk)
c----------------------------------------------------------------------
c     cmqmul  transforms the vector  v  in various ways using the
c     matrix  q = (z  y)  defined by the input parameters.

c        mode               result
c        ----               ------
c          1                v = z v
c          2                v = y v
c          3                v = q v

c     on input,  v  is assumed to be ordered as  (v(free)  v(fixed)).
c     on output, v  is a full n-vector.

c          4                v = z'v
c          5                v = y'v
c          6                v = q'v

c     on input,  v  is a full n-vector.
c     on output, v  is ordered as  (v(free)  v(fixed)).

c          7                v = y'v
c          8                v = q'v

c     on input,  v  is a full n-vector.
c     on output, v  is as in modes 5 and 6 except that v(fixed) is not
c     set.

c     modes  1, 4, 7 and 8  do not involve  v(fixed).
c----------------------------------------------------------------------
      implicit none

      integer mode, n, nfree, nq, nz
      logical unitq

      double precision v(n), wrk(n), zy(nq,*)
      integer kx(n)

      integer j, j1, j2, k, l, lenv, nfixed
c----------------------------------------------------------------------
      nfixed = n - nfree
      j1 = 1
      j2 = nfree
      if (mode.eq.1 .or. mode.eq.4) j2 = nz
      if (mode.eq.2 .or. mode.eq.5 .or. mode.eq.7) j1 = nz + 1
      lenv = j2 - j1 + 1
      if (mode.le.3) then

c        mode = 1, 2  or  3.

         if (nfree.gt.0) call sload (nfree,0d0,wrk,1)

c        copy  v(fixed)  into the end of  wrk.

         if (mode.ge.2 .and. nfixed.gt.0) call dcopy (nfixed,v(nfree+1),
     *       1,wrk(nfree+1),1)

c        set  wrk  =  relevant part of  zy * v.

         if (lenv.gt.0) then
            if (unitq) then
               call dcopy (lenv,v(j1),1,wrk(j1),1)
            else
               call dgemv ('n',nfree,j2-j1+1,1d0,zy(1,j1),nq,v(j1),
     *                    1d0,wrk)
            end if
         end if

c        expand  wrk  into  v  as a full n-vector.

         call sload (n,0d0,v,1)

         do k = 1, nfree
            j = kx(k)
            v(j) = wrk(k)
          end do

c        copy  wrk(fixed)  into the appropriate parts of  v.

         if (mode.gt.1) then
            do l = 1, nfixed
               j = kx(nfree+l)
               v(j) = wrk(nfree+l)
            end do
         end if

      else

c        mode = 4, 5, 6, 7  or  8.
c        put the fixed components of  v  into the end of  wrk.

         if (mode.eq.5 .or. mode.eq.6) then
            do l = 1, nfixed
               j = kx(nfree+l)
               wrk(nfree+l) = v(j)
            end do
         end if

c        put the free  components of  v  into the beginning of  wrk.

         if (nfree.gt.0) then
            do k = 1, nfree
               j = kx(k)
               wrk(k) = v(j)
            end do

c           set  v  =  relevant part of  zy' * wrk.

            if (lenv.gt.0) then
               if (unitq) then
                  call dcopy (lenv,wrk(j1),1,v(j1),1)
               else
                  call dgemv ('t',nfree,j2-j1+1,1d0,zy(1,j1),nq,wrk,
     *                       0d0,v(j1))
               end if
            end if
         end if

c        copy the fixed components of  wrk  into the end of  v.

         if (nfixed.gt.0 .and. (mode.eq.5 .or. mode.eq.6))
     *       call dcopy (nfixed,wrk(nfree+1),1,v(nfree+1),1)

      end if
c                                 end of cmqmul
      end

      subroutine lscore (prbtyp,linobj,unitq,inform,iter,
     *                   jinf,nclin,nctotl,nactiv,nfree,nrank,nz,nrz,n,
     *                   lda,ldr,istate,kactiv,kx,ctx,ssq,ssq1,suminf,
     *                   numinf,xnorm,bl,bu,a,clamda,ax,featol,r,x,w)
c----------------------------------------------------------------------
c     lscore is a subroutine for linearly constrained linear-least
c     squares. 

c     values of istate(j) for the linear constraints.......

c     istate(j)

c          0    constraint j is not in the working set.
c          1    constraint j is in the working set at its lower bound.
c          2    constraint j is in the working set at its upper bound.
c          3    constraint j is in the working set as an equality.

c     constraint j may be violated by as much as featol(j).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*2 prbtyp

      logical linobj, unitq , convrg, cyclin, error, firstv, hitcon, 
     *       hitlow, needfg, overfl, prnt, rowerr, singlr, stall,
     *       statpt, unbndd, uncon, unitgz, weak

      integer inform, iter, jinf, lda, ldr, n, nactiv, nclin,
     *        nctotl, nfree, nrank, nrz, numinf, nz, istate(nctotl), 
     *        kactiv(n), kx(n), iadd, ifix, irefn, is, isdel, itmax,
     *        jadd, jbigst, jdel, jmax1, jsmlst, jtiny, kbigst,
     *        kdel, ksmlst, ngq, nphase, nres, nstall, nviol

      double precision ctx, ssq, ssq1, suminf, xnorm, a(lda,*), ax(*), 
     *                 bl(nctotl), bu(nctotl), cnorm,
     *                 clamda(nctotl), featol(nctotl), r(ldr,*), w(*),
     *                 x(n), absrzz, alfa, alfhit, atphit, bigalf, 
     *                 condmx, condrz, condt, ctp, dinky, drzmax,
     *                 drzmin, err1, err2, flmax, gfnorm, grznrm,
     *                 gznorm, objsiz, palfa, pnorm, resnrm, rownrm,
     *                 trulam, wssize, dnrm2, sdiv

      external dnrm2, sdiv 

      integer lkactv, lanorm, lap, lpx, lres, lres0, lhz,
     *        lgq, lcq, lrlam, lt, lzy, lwtinf, lwrk, lqptol
      common/ cstlnp /lkactv, lanorm, lap, lpx, lres, lres0, lhz,
     *                lgq, lcq, lrlam, lt, lzy, lwtinf, lwrk, lqptol

      double precision wmach
      common/ cstmch /wmach(10)

      integer ldt, ncolt, ldzy
      common/ ngg004 /ldt, ncolt, ldzy

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin

      integer itmax1, itmax2, lprob
      common/ ngg016 /itmax1, itmax2, lprob
c----------------------------------------------------------------------
      flmax = wmach(7)

c     set up the adresses of the contiguous arrays  (res0, res)
c     and  (gq, cq).

      nres = 0
      if (nrank.gt.0) nres = 2
      ngq = 1
      if (linobj) ngq = 2

c     initialize.

      irefn = 0
      iter = 0
      jadd = 0
      jdel = 0
      nphase = 1
      nstall = 0
      numinf = -1
      nrz = 0

      if (prbtyp.eq.'fp') then
         itmax = itmax2
      else
         itmax = itmax1
      end if

      alfa = 0d0
      condmx = flmax
      drzmax = 1d0
      drzmin = 1d0
      ssq = 0d0

      cyclin = .false.
      error = .false.
      firstv = .false.
      prnt = .true.
      needfg = .true.
      stall = .true.
      uncon = .false.
      unitgz = .true.
      unbndd = .false.

c     =================== start of the main loop =======================
c     repeat
c        repeat
   20 if (needfg) then
         if (nrank.gt.0) then
            resnrm = dnrm2 (nrank,w(lres),1)
            ssq = (ssq1**2+resnrm**2)/2d0
         end if
c
         if (numinf.ne.0) then

c           compute the transformed gradient of either the sum of
c           infeasibilities or the objective.  initialize
c           singlr and unitgz.

            call lsgset(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,w(lres),featol,
     *                  w(lgq),w(lcq),r,x,w(lwtinf),w(lzy),w(lwrk))
            if (numinf.eq.0 .and. prbtyp.ne.'fp' .and. nphase.eq.1) then
               itmax = iter + itmax2
               nphase = 2
            end if
         end if
      end if

      gznorm = 0d0
      if (nz.gt.0) gznorm = dnrm2 (nz,w(lgq),1)

      if (nrz.eq.nz) then
         grznrm = gznorm
      else
         grznrm = 0d0
         if (nrz.gt.0) grznrm = dnrm2 (nrz,w(lgq),1)
      end if

      gfnorm = gznorm
      if (nfree.gt.0 .and. nactiv.gt.0) gfnorm = dnrm2 (nfree,w(lgq),1)

c        define small quantities that reflect the magnitude of  x,
c        r  and  the matrix of constraints in the working set.
c        use the largest and smallest diagonals of  r  to estimate
c        the condition number of  rz1.

      if (nrz.eq.0) then
         singlr = .false.
      else
         if (numinf.gt.0 .or. nrz.gt.nrank) then
            absrzz = 0d0
            singlr = .true.
         else
            call scond (nrz,r,ldr+1,drzmax,drzmin)
            absrzz = abs(r(nrz,nrz))
            rownrm = dnrm2 (n,r(1,1),ldr)
            singlr = absrzz .le. drzmax*tolrnk .or. rownrm .le.
     *               tolrnk .or. abs(r(1,1)) .le. rownrm*tolrnk
         end if

      end if

      condrz = sdiv (drzmax,drzmin,overfl)

      condt = 1d0
      if (nactiv.gt.0) condt = sdiv (dtmax,dtmin,overfl)

      if (prnt) then

         jdel = 0
         jadd = 0
         alfa = 0d0
      end if

      if (numinf.gt.0) then
         dinky = 0d0
      else
         objsiz = 1d0 + abs(ssq+ctx)


         wssize = 0d0
         if (nactiv.gt.0) wssize = dtmax
         dinky = epspt8*max(wssize,objsiz,gfnorm)
         if (uncon) then
            unitgz = grznrm .le. dinky
         end if
      end if

c     if the projected gradient  z'g  is small and rz is of full
c     rank, x is a minimum on the working set.  an additional
c     refinement step is allowed to take care of an inaccurate
c     value of dinky.

      statpt = .not. singlr .and. grznrm .le. dinky .or. irefn.gt.1
      if (.not. statpt) then

c        compute a search direction.

         prnt = .true.

         error = iter .ge. itmax
         if (.not. error) then

            irefn = irefn + 1
            iter = iter + 1

            call lsgetp(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,w(lap),
     *                  w(lres),w(lhz),w(lpx),w(lgq),w(lcq),r,w(lzy),
     *                  w(lwrk))

c           find the constraint we bump into along p.
c           update x and ax if the step alfa is nonzero.

c           alfhit is initialized to bigalf.  if it remains
c           that way after the call to cmalf, it will be
c           regarded as infinite.

            bigalf = sdiv (bigdx,pnorm,overfl)

            call cmalf(firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfhit,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  w(lanorm),w(lap),ax,bl,bu,featol,w(lpx),x)

c           if  rz1  is nonsingular,  alfa = 1.0  will be the
c           step to the least-squares minimizer on the
c           current subspace. if the unit step does not violate
c           the nearest constraint by more than featol,  the
c           constraint is not added to the working set.

            hitcon = singlr .or. palfa .le. 1d0
            uncon = .not. hitcon

            if (hitcon) then
               alfa = alfhit
            else
               jadd = 0
               alfa = 1d0
            end if

c           check for an unbounded solution or negligible step.

            unbndd = alfa .ge. bigalf
            stall = abs(alfa*pnorm) .le. epspt9*xnorm
            if (stall) then
               nstall = nstall + 1
               cyclin = nstall.gt.50
            else
               nstall = 0
            end if

            error = unbndd .or. cyclin
            if (.not. error) then

c              set x = x + alfa*p.  update ax, gq, res and ctx.

               if (alfa.ne.0d0) call lsmove(hitcon,hitlow,linobj,
     *                                unitgz,nclin,nrank,nrz,n,ldr,jadd,
     *                                numinf,alfa,ctp,ctx,xnorm,w(lap),
     *                                ax,bl,bu,w(lgq),w(lhz),w(lpx),
     *                                w(lres),r,x,w(lwrk))

               if (hitcon) then

c                 add a constraint to the working set.
c                 update the tq factors of the working set.
c                 use p as temporary work space.
c                 update  istate.

                  if (bl(jadd).eq.bu(jadd)) then
                     istate(jadd) = 3
                  else if (hitlow) then
                     istate(jadd) = 1
                  else
                     istate(jadd) = 2
                  end if
                  iadd = jadd - n
                  if (jadd.le.n) then

                     do 40 ifix = 1, nfree
                        if (kx(ifix).eq.jadd) go to 60
   40                continue
                  end if
   60             continue

                  call lsadd (unitq,inform,ifix,iadd,jadd,nactiv,nz,
     *                        nfree,nrank,nres,ngq,n,lda,ldzy,ldr,ldt,
     *                        kx,condmx,a,r,w(lt),w(lres),w(lgq),w(lzy),
     *                        w(lwrk),w(lrlam),w(lpx))

                  nrz = nrz - 1
                  nz = nz - 1

                  if (jadd.le.n) then

c                    a simple bound has been added.

                     nfree = nfree - 1
                  else

c                    a general constraint has been added.

                     nactiv = nactiv + 1
                     kactiv(nactiv) = iadd
                  end if

                  irefn = 0
               end if

c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  refine the current x and set inform so
c              that feasibility is checked in lsgset.

               call lsfeas (n,nclin,istate,bigbnd,cnorm,err1,jmax1,
     *                     nviol,ax,bl,bu,featol,x,w(lwrk))

               if (err1.gt.featol(jmax1)) then
                  call lssetx (linobj,rowerr,unitq,nclin,nactiv,nfree,
     *                        nrank,nz,n,nctotl,ldzy,lda,ldr,ldt,istate,
     *                        kactiv,kx,jmax1,err2,ctx,xnorm,a,ax,bl,bu,
     *                        w(lcq),w(lres),w(lres0),featol,r,w(lt),x,
     *                        w(lzy),w(lpx),w(lwrk))

                  if (rowerr) then

                     numinf = 1
                     error = .true.
                  else
                     numinf = -1
                     uncon = .false.
                     irefn = 0
                  end if
               end if
               needfg = alfa .ne. 0d0
            end if
         end if
      end if

c        until      statpt  .or.  error
      if (.not. (statpt .or. error)) go to 20

c     try and find the index jdel of a constraint to drop from
c     the working set.

      jdel = 0
c
      if (numinf.eq.0 .and. prbtyp.eq.'fp') then
         if (n.gt.nz) call sload (n-nz,0d0,w(lrlam),1)
         jtiny = 0
         jsmlst = 0
         jbigst = 0
      else

         call lsmuls (n,nactiv,nfree,lda,ldt,numinf,nz,nrz,
     *                istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,jtiny,
     *                jbigst,kbigst,trulam,a,w(lanorm),w(lgq),w(lrlam),
     *                w(lt),w(lwtinf))

      end if

      if (.not. error) then
         if (jsmlst.gt.0) then

c           lsmuls found a regular constraint with multiplier less
c           than (-dinky).

            jdel = jsmlst
            kdel = ksmlst
            isdel = istate(jdel)
            istate(jdel) = 0

         else if (jsmlst.lt.0) then

            jdel = jsmlst

         else if (numinf.gt.0 .and. jbigst.gt.0) then

c           no feasible point exists for the constraints but the
c           sum of the constraint violations may be reduced by
c           moving off constraints with multipliers greater than 1.

            jdel = jbigst
            kdel = kbigst
            isdel = istate(jdel)
            if (trulam.le.0d0) is = -1
            if (trulam.gt.0d0) is = -2
            istate(jdel) = is
            firstv = .true.
            numinf = numinf + 1
         end if
c
         if (jdel.ne.0 .and. singlr) then
c
c           cannot delete a constraint when rz is singular.
c           probably a weak minimum.
c
            jdel = 0
         else if (jdel.ne.0) then

c           constraint jdel has been deleted.
c           update the matrix factorizations.

            call lsdel(unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,w(lres),r,
     *                  w(lt),w(lgq),w(lzy),w(lwrk),w(lpx))
         end if
      end if

      irefn = 0
      convrg = jdel.eq.0

      prnt = .false.
      uncon = .false.
      needfg = .false.

c     until       convrg  .or.  error
      if (.not. (convrg .or. error)) go to 20

c     .....................end of main loop........................
      weak = jtiny.gt.0 .or. singlr
c
      if (error) then
         if (unbndd) then
            inform = 2
            if (numinf.gt.0) inform = 3
         else if (iter.ge.itmax) then
            inform = 4
         else if (cyclin) then
            inform = 5
         end if
      else if (convrg) then
         inform = 0
         if (numinf.gt.0) then
            inform = 3
         else if (prbtyp.ne.'fp' .and. weak) then
            inform = 1
         end if
      end if
c                                 set clamda

      call cmprt (nfree,n,nctotl,
     *            nactiv,kactiv,kx,clamda,w(lrlam))

c                                 end of lscore
      end

      subroutine cmprnt (nfree,n,nctotl,nactiv,kactiv,kx,clamda,rlamda)
c----------------------------------------------------------------------
c cmprnt creates the expanded lagrange multiplier vector clamda.
c----------------------------------------------------------------------
      implicit none

      integer n, nactiv, nctotl, nfree, kactiv(n), kx(n),
     *        j, k, nfixed, nz

      double precision clamda(nctotl), rlamda(n), rlam
c----------------------------------------------------------------------
      nz = nfree - nactiv

c     expand multipliers for bounds, linear and nonlinear constraints
c     into the  clamda  array.

      call sload (nctotl,0d0,clamda,1)
      nfixed = n - nfree

      do k = 1, nactiv + nfixed
         if (k.le.nactiv) then
            j = kactiv(k) + n
            rlam = rlamda(nactiv-k+1)
         else
            j = kx(nz+k)
            rlam = rlamda(k)
         end if
         clamda(j) = rlam
      end do
c                                 end of cmprnt
      end

      subroutine lsbnds (unitq,inform,nz,nfree,nrank,nres,ngq,n,ldzy,lda
     *                  ,ldr,ldt,istate,kx,condmx,a,r,t,res,gq,zy,w,c,s)
c----------------------------------------------------------------------
c     lsbnds updates the factor r as kx is reordered to reflect the
c     status of the bound constraints given by istate.  kx is reordered
c     so that the fixed variables come last.  one of two alternative
c     are used to reorder kx. one method needs fewer accesses to kx, the
c     other gives a matrix rz with more rows and columns.
c----------------------------------------------------------------------
      implicit none

      integer inform, lda, ldr, ldt, ldzy, n, nfree,
     *        ngq, nrank, nres, nz, istate(*), kx(n), iadd, 
     *        ifix, j, j2, jadd, k, l, lstart, nactv, nfixed

      logical unitq

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                 s(n), t(ldt,*), w(n), zy(ldzy,*), condmx
c----------------------------------------------------------------------
      nfixed = n - nfree

      if (nrank.lt.n .and. nrank.gt.0) then

c        r is specified but singular.  try and keep the dimension of rz
c        as large as possible.

         nactv = 0
         nfree = n
         nz = n
         j = n
c        +       while (j.gt.0  .and.  n-nfree.lt.nfixed) do
   20    if (j.gt.0 .and. n-nfree.lt.nfixed) then
            if (istate(j).gt.0) then
               jadd = j
               do 40 ifix = nfree, 1, -1
                  if (kx(ifix).eq.jadd) go to 60
   40          continue

c              add bound jadd.

   60          call lsadd (unitq,inform,ifix,iadd,jadd,nactv,nz,nfree,
     *                     nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,
     *                     a,r,t,res,gq,zy,w,c,s)

               nfree = nfree - 1
               nz = nz - 1
            end if
            j = j - 1
            go to 20
c           +       end while
         end if
      else

c        r is of full rank,  or is not specified.

         if (nfixed.gt.0) then

c           order kx so that the free variables come first.

            lstart = nfree + 1
            do 120 k = 1, nfree
               j = kx(k)
               if (istate(j).gt.0) then
                  do 80 l = lstart, n
                     j2 = kx(l)
                     if (istate(j2).eq.0) go to 100
   80             continue

  100             kx(k) = j2
                  kx(l) = j
                  lstart = l + 1

                  if (nrank.gt.0) call nggnbu(n,nres,nrank,ldr,k,l,r,
     *                                 res,c,s)
               end if
  120       continue

         end if
         nz = nfree
      end if
c                                 end of lsbnds
      end

      subroutine cmdgen (job,n,nclin,nmoved,iter,numinf,istate,
     *                   bl,bu,featol,featlu,x)
c----------------------------------------------------------------------
c     cmdgen does most of the manoeuvres associated with degeneracy.
c     the degeneracy-resolving strategy operates in the following way.

c     over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.

c     after kdegen consecutive iterations, variables within
c     featol of their bounds are set exactly on their bounds and x is
c     recomputed to satisfy the general constraints in the working set.
c     featol is then reduced to tolx0 for the next cycle of iterations.

c     featlu  is the array of user-supplied feasibility tolerances.
c     featol  is the array of current feasibility tolerances.

c     tolx0   is the minimum (scaled) feasibility tolerance.
c     tolx1   is the maximum (scaled) feasibility tolerance.
c     tolinc  is the scaled increment to the current featol.
c     kdegen  is the expand frequency (specified by the user).
c             it is the frequency of resetting featol to (scaled) tolx0.
c     ndegen  counts the number of degenerate steps (incremented
c             by cmchzr).
c     itnfix  is the last iteration at which a job = 'e' or 'o' entry
c             caused an x to be put on a constraint.
c     nfix(j) counts the number of times a job = 'o' entry has
c             caused the variables to be placed on the working set,
c             where j=1 if infeasible, j=2 if feasible.

c     tolx0*featlu and tolx1*featlu are both close to the feasibility
c     tolerance featlu specified by the user.  (they must both be less
c     than featlu.)

c     if job = 'e',  cmdgen has been called after a cycle of kdegen
c     iterations.  constraints in the working set are examined to see if
c     any are off their bounds by an amount approaching featol.  nmoved
c     returns how many.  if nmoved is positive,  x  is moved onto the
c     constraints in the working set.  it is assumed that the calling
c     routine will then continue iterations.

c     if job = 'o',  cmdgen is being called after a subproblem has been
c     judged optimal, infeasible or unbounded.  constraint violations
c     are examined as above.
c----------------------------------------------------------------------
      implicit none

      character job

      integer iter, n, nclin, nmoved, numinf, istate(n+nclin), is, j, 
     *        maxfix

      double precision bl(n+nclin), bu(n+nclin), d,
     *                 featlu(n+nclin), featol(n+nclin), x(n)

      double precision wmach
      common/ cstmch /wmach(10)

      integer itnfix, kdegen, ndegen, nfix
      double precision tolinc, tolx0
      common/ ngg005 /tolx0, tolinc, kdegen, ndegen, itnfix, nfix(2)
c----------------------------------------------------------------------
      nmoved = 0
c                                   job = 'end of cycle' or 'optimal'.
c                                   initialize local variables
      maxfix = 2

      if (job.eq.'o') then
c                                   job = 'optimal'.
c                                   return with nmoved = 0 if the last call was at the same
c                                   iteration,  or if there have already been maxfix calls with
c                                   the same state of feasibility.
         if (itnfix.eq.iter) return

         if (numinf.gt.0) then
            j = 1
         else
            j = 2
         end if

         if (nfix(j).ge.maxfix) return
         nfix(j) = nfix(j) + 1

      end if
c                                 reset featol to its minimum value.
      do j = 1, n + nclin
         featol(j) = tolx0*featlu(j)
      end do
c                                 count the number of times a variable is moved a nontrivial
c                                 distance onto its bound.
      itnfix = iter

      do j = 1, n
            is = istate(j)
            if (is.gt.0 .and. is.lt.4) then
               if (is.eq.1) then
                  d = abs(x(j)-bl(j))
               else
                  d = abs(x(j)-bu(j))
               end if

               if (d.gt.wmach(3)**0.6d0) nmoved = nmoved + 1
            end if
      end do
c                                 end of cmdgen
      end

      subroutine nprset (unitq,n,nfree,nz,nq,nrowr,iperm,kx,gq,r,zy,
     *                   work,qrwork)
c----------------------------------------------------------------------
c     nprset  bounds the condition estimator of the transformed hessian.
c     on exit, r is of the form
c                  (drz   0    )
c                  ( 0  sigma*i)
c     where d is a diagonal matrix such that drz has a bounded condition
c     number,  i is the identity matrix and sigma  is the geometric mean
c     of the largest and smallest elements of drz. the qr factorization
c     with interchanges is used to give diagonals of drz that are
c     decreasing in modulus.
c----------------------------------------------------------------------
      implicit none

      logical unitq

      integer n, nfree, nq, nrowr, nz, iperm(n), kx(n), j, jmax, jsave,
     *        nrank, isrank

      double precision gq(n), qrwork(2*n), r(nrowr,*), work(n),
     *                  zy(nq,*), drgm, drgs, gjmax, scle, sumsq, snorm

      double precision rcndbd, rfrobn, drmax, drmin
      common/ ngg018 /rcndbd, rfrobn, drmax, drmin
c----------------------------------------------------------------------
c                                 bound the condition estimator of q'hq.
      if (nz.gt.1) then
c                                 refactorize rz. interchanges give diagonals
c                                 of decreasing magnitude.
         do j = 1, nz - 1
            call sload (nz-j,0d0,r(j+1,j),1)
         end do

         call sgeqrp ('c',nz,nz,r,nrowr,work,iperm,qrwork)

         do j = 1, nz
            jmax = iperm(j)
            if (jmax.gt.j) then
               if (unitq) then
                  jsave = kx(jmax)
                  kx(jmax) = kx(j)
                  kx(j) = jsave
               else
                  call dswap(nfree,zy(1,jmax),zy(1,j),1)
               end if
c
               gjmax = gq(jmax)
               gq(jmax) = gq(j)
               gq(j) = gjmax
            end if
         end do

      end if

      drgm = 1d0

      if (nz.gt.0) then

         nrank = isrank(nz,r,nrowr+1,1d0/rcndbd)
         drgm = sqrt(abs(r(1,1)*r(nrank,nrank)))/2d0
         drgs = abs(r(1,1))/rcndbd

         if (nz.gt.nrank) then

            do j = nrank + 1, nz
               call sload (j-1,0d0,r(1,j),1)
            end do 

            call sload (nz-nrank,drgs,r(nrank+1,nrank+1),nrowr+1)

         end if

      end if
c                                 reset the range-space partition of the hessian.
      if (nz.lt.n) then

         do j = nz + 1, n
            call sload (j,0d0,r(1,j),1)
         end do

         call sload (n-nz,drgm,r(nz+1,nz+1),nrowr+1)

      end if
c                                 recompute the frobenius norm of r.
      scle = sqrt(dble(n-nz))*drgm
      sumsq = 1d0

      do j = 1, nz
         call sssq (j,r(1,j),1,scle,sumsq)
      end do

      rfrobn = snorm(scle,sumsq)
c                                 end of nprset
      end

      subroutine cmcrsh (istart,vertex,nclin,nctotl,nactiv,nartif,nfree,
     *                  n,lda,istate,kactiv,kx,bigbnd,tolact,a,ax,bl,bu,
     *                  featol,x,wx,work)
c-----------------------------------------------------------------------
c     cmcrsh  computes the quantities  istate (optionally),  kactiv,
c     nactiv,  nz  and  nfree  associated with the working set at x.
c
c     the computation depends upon the value of the input parameter
c     start,  as follows...
c
c     istart = 0      an initial working set will be selected. first,
c                     nearly-satisfied or violated bounds are added.
c                     next,  general linear constraints are added that
c                     have small residuals.
c
c     istart = 1      the quantities kactiv, nactiv and nfree are
c                     initialized from istate,  specified by the user.
c
c     if vertex is true, an artificial vertex is defined by fixing some
c     variables on their bounds.  infeasible variables selected for the
c     artificial vertex are fixed at their nearest bound.  otherwise,
c     the variables are unchanged.
c
c     values of istate(j)....
c
c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c-----------------------------------------------------------------------
      implicit none 

      double precision bigbnd, tolact
      integer lda, n, nactiv, nartif, nclin, nctotl, nfree, istart
      logical vertex

      double precision a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), work(n), wx(n), x(n)
      integer istate(nctotl), kactiv(n), kx(n)

      double precision b1, b2, biglow, bigupp, colmin, colsiz, flmax,
     *                  residl, resl, resmin, resu, tol, toobig
      integer i, imin, is, j, jfix, jfree, jmin, k

      double precision ddot
      external ddot

      double precision wmach
      common/ cstmch /wmach(10)
c-----------------------------------------------------------------------
      flmax = wmach(7)
      biglow = -bigbnd
      bigupp = bigbnd

      if (istart.eq.0) then
c                                 cold, initialize to lower bound
         x(1:n) = bl(1:n)

      else 
c                                 warm/hot move the variables inside their bounds.
         do j = 1, n

            b1 = bl(j)
            b2 = bu(j)
            tol = featol(j)

            if (b1.gt.biglow) then
               if (x(j).lt.b1-tol) x(j) = b1
            end if

            if (b2.lt.bigupp) then
               if (x(j).gt.b2+tol) x(j) = b2
            end if
         end do

      end if

      wx(1:n) = x(1:n)

      nfree = n
      nactiv = 0
      nartif = 0

      if (istart.eq.0) then
c                                 cold
         do j = 1, nctotl
            istate(j) = 0
            if (bl(j).eq.bu(j)) istate(j) = 3
         end do

      else
c                                warm (cmcrsh not called for hot)
         do j = 1, nctotl
            if (istate(j).gt.3 .or. istate(j).lt.0) istate(j) = 0
            if (bl(j).ne.bu(j) .and. istate(j).eq.3) istate(j) = 0
         end do

      end if

c     define nfree and kactiv.
c     ensure that the number of bounds and general constraints in the
c     working set does not exceed n.

      do j = 1, nctotl
         if (nactiv.eq.nfree) istate(j) = 0

         if (istate(j).gt.0) then
            if (j.le.n) then
               nfree = nfree - 1

               if (istate(j).eq.1) then
                  wx(j) = bl(j)
               else if (istate(j).ge.2) then
                  wx(j) = bu(j)
               end if
            else
               nactiv = nactiv + 1
               kactiv(nactiv) = j - n
            end if
         end if
      end do

c     if a cold start is required, attempt to add as many
c     constraints as possible to the working set.

      if (istart.eq.0) then

c        see if any bounds are violated or nearly satisfied.
c        if so, add these bounds to the working set and set the
c        variables exactly on their bounds.

         j = n
c        +       while (j .ge. 1  .and.  nactiv.lt.nfree) do
  100    if (j.ge.1 .and. nactiv.lt.nfree) then
            if (istate(j).eq.0) then
               b1 = bl(j)
               b2 = bu(j)
               is = 0
               if (b1.gt.biglow) then
                  if (wx(j)-b1.le.(1d0+abs(b1))*tolact) is = 1
               end if
               if (b2.lt.bigupp) then
                  if (b2-wx(j).le.(1d0+abs(b2))*tolact) is = 2
               end if
               if (is.gt.0) then
                  istate(j) = is
                  if (is.eq.1) wx(j) = b1
                  if (is.eq.2) wx(j) = b2
                  nfree = nfree - 1
               end if
            end if
            j = j - 1
            go to 100
c           +       end while
         end if

c        the following loop finds the linear constraint (if any) with
c        smallest residual less than or equal to tolact  and adds it
c        to the working set.  this is repeated until the working set
c        is complete or all the remaining residuals are too large.

c        first, compute the residuals for all the constraints not in the
c        working set.

         if (nclin.gt.0 .and. nactiv.lt.nfree) then

            do 120 i = 1, nclin
               if (istate(n+i).le.0) ax(i) = ddot (n,a(i,1),lda,wx)
  120       continue

            is = 1
            toobig = tolact + tolact

c           +          while (is.gt.0  .and.  nactiv.lt.nfree) do
  140       if (is.gt.0 .and. nactiv.lt.nfree) then
               is = 0
               resmin = tolact

               do 160 i = 1, nclin
                  j = n + i
                  if (istate(j).eq.0) then
                     b1 = bl(j)
                     b2 = bu(j)
                     resl = toobig
                     resu = toobig
                     if (b1.gt.biglow) resl = abs(ax(i)-b1)/(1d0+abs(b1)
     *                                       )
                     if (b2.lt.bigupp) resu = abs(ax(i)-b2)/(1d0+abs(b2)
     *                                       )
                     residl = min(resl,resu)
                     if (residl.lt.resmin) then
                        resmin = residl
                        imin = i
                        is = 1
                        if (resl.gt.resu) is = 2
                     end if
                  end if
  160          continue
c
               if (is.gt.0) then
                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j = n + imin
                  istate(j) = is
               end if
               go to 140
c              +          end while
            end if
         end if
      end if

      if (vertex .and. nactiv.lt.nfree) then

c        find an initial vertex by temporarily fixing some variables.

c        compute lengths of columns of selected linear constraints
c        (just the ones corresponding to variables eligible to be
c        temporarily fixed).

         do 200 j = 1, n
            if (istate(j).eq.0) then
               colsiz = 0d0
               do 180 k = 1, nclin
                  if (istate(n+k).gt.0) colsiz = colsiz + abs(a(k,j))
  180          continue
               work(j) = colsiz
            end if
  200    continue
c
c        find the  nartif  smallest such columns.
c        this is an expensive loop.  later we can replace it by a
c        4-pass process (say), accepting the first col that is within
c        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).

c        +       while (nactiv.lt.nfree) do
  220    if (nactiv.lt.nfree) then
            colmin = flmax
            do 240 j = 1, n
               if (istate(j).eq.0) then
                  if (nclin.eq.0) go to 260
                  colsiz = work(j)
                  if (colmin.gt.colsiz) then
                     colmin = colsiz
                     jmin = j
                  end if
               end if
  240       continue
            j = jmin

c           fix x(j) at its current value.

  260       istate(j) = 4
            nartif = nartif + 1
            nfree = nfree - 1
            go to 220
c           +       end while
         end if
      end if

      jfree = 1
      jfix = nfree + 1
      do 280 j = 1, n
         if (istate(j).le.0) then
            kx(jfree) = j
            jfree = jfree + 1
         else
            kx(jfix) = j
            jfix = jfix + 1
         end if
  280 continue
c                                 end of cmcrsh
      end

c START OF BLAS 2. 

      subroutine dtrmv (uplo, trans, diag, n, a, lda, x)
c-----------------------------------------------------------------------
c  dtrmv  does one of the matrix-vector operations
c     x := a*x,   or   x := a'*x,
c  where x is n element vector and a is an n by n unit, or non-unit,
c  upper or lower triangular matrix.
c-----------------------------------------------------------------------
      implicit none

      character diag*1, trans*1, uplo*1

      logical  nounit

      integer  lda, n, i, j

      double precision  a(lda,*), x(*), temp
c-----------------------------------------------------------------------
      if (n.eq.0) return

      nounit = (diag.eq.'n')

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.

      if (trans.eq.'n') then

c                                 form  x := a*x.

         if (uplo.eq.'u') then

               do j = 1, n
                  if (x(j).ne.0d0) then
                     temp = x(j)
                     do i = 1, j - 1
                        x(i) = x(i) + temp*a(i,j)
                     end do

                     if (nounit) x(j) = x(j)*a(j,j)

                  end if
               end do

         else

               do 60, j = n, 1, -1
                  if (x(j).ne.0d0) then
                     temp = x(j)
                     do 50, i = n, j + 1, -1
                        x(i) = x(i) + temp*a(i,j)
   50                continue
                     if (nounit) x(j) = x(j)*a(j,j)
                  end if
   60          continue

         end if

      else

c        form  x := a'*x.

         if (uplo.eq.'u') then

               do 100, j = n, 1, -1
                  temp = x(j)
                  if (nounit) temp = temp*a(j,j)
                  do 90, i = j - 1, 1, -1
                     temp = temp + a(i,j)*x(i)
   90             continue
                  x(j) = temp
  100          continue

         else

               do 140, j = 1, n
                  temp = x(j)
                  if (nounit) temp = temp*a(j,j)
                  do 130, i = j + 1, n
                     temp = temp + a(i,j)*x(i)
  130             continue
                  x(j) = temp
  140          continue

         end if

      end if
c                                 end of dtrmv.
      end

      integer function idamax (n, x)
c-----------------------------------------------------------------------
c  idamax returns the smallest value of i such that
c     abs(x(i)) = max(abs(x(j)))
c-----------------------------------------------------------------------
      implicit none

      integer n, i, imax

      double precision x(*), xmax
c----------------------------------------------------------------------
      if (n.gt.0) then
         imax = 1
         if (n.gt.1) then
            xmax = abs(x(1))

            do 10, i = 2, n

               if (xmax.lt.abs(x(i))) then
                  xmax = abs(x(i))
                  imax = i
               end if
   10       continue
         end if
      else
         imax = 0
      end if

      idamax = imax
c                                 end of idamax
      end

      subroutine daxpy (n, alpha, x, incx, y, incy)
c-----------------------------------------------------------------------
c  daxpy  does  y := alpha*x + y
c----------------------------------------------------------------------
      implicit none

      integer incx, incy, n, i, ix, iy

      double precision alpha, x(*), y(*)
c----------------------------------------------------------------------
      if (n.gt.0) then
         if (alpha.ne.0d0) then
            if ((incx.eq.incy).and.(incx.gt.0)) then
               do 10, ix = 1, 1 + (n - 1)*incx, incx
                  y(ix) = alpha*x(ix) + y(ix)
   10          continue
            else
               if (incy.ge.0) then
                  iy = 1
               else
                  iy = 1 - (n - 1)*incy
               end if
               if (incx.gt.0) then
                  do 20, ix = 1, 1 + (n - 1)*incx, incx
                     y(iy) = alpha*x(ix) + y(iy)
                     iy      = iy            + incy
   20             continue
               else
                  ix = 1 - (n - 1)*incx
                  do 30, i = 1, n
                     y(iy) = alpha*x(ix) + y(iy)
                     ix      = ix            + incx
                     iy      = iy            + incy
   30             continue
               end if
            end if
         end if
      end if
c                                 end of daxpy
      end

      subroutine dgemv (trans, m, n, alpha, a, lda, x, beta, y)
c----------------------------------------------------------------------
c  dgemv  does one of the matrix-vector operations
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.
c-----------------------------------------------------------------------
      implicit none
 
      integer lda, m, n, i, j, jx, lenx, leny, m4, n4

      character trans*1

      double precision a(lda,*), x(*), y(*), temp, temp1, temp2, temp3,
     *                 temp4, alpha, beta
c-----------------------------------------------------------------------
      if (m.eq.0.or.n.eq.0.or.(alpha.eq.0d0.and.beta.eq.1d0)) return

      if (trans.eq.'n') then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if

c DEBUG DEBUG avoids uninitialized ax/adx

      if (beta.eq.0d0) then 

         y(1:leny) = 0d0

      else if (beta.ne.1d0) then

         y(1:leny) = beta*y(1:leny)

      end if

      if (alpha.eq.0d0) return

      jx = 1

      if (trans.eq.'n') then
c                                 form  y := alpha*a*x + y.
c**** u n r o l l   t o   d e p t h   4 ********************************
            n4 = 4*(n/4)
            do 60, j = 1, n4, 4
               temp1 = alpha*x(jx)
               temp2 = alpha*x(jx + 1)
               temp3 = alpha*x(jx + 2)
               temp4 = alpha*x(jx + 3)
               if (temp1.ne.0d0.or.temp2.ne.0d0.or.temp3.ne.0d0.or.
     *             temp4.ne.0d0) then
                  do 50, i = 1, m
                     y(i) = ((((y(i) + temp1*a(i,j))
     *                        + temp2*a(i, j+1))
     *                        + temp3*a(i, j + 2))
     *                        + temp4*a(i, j + 3))
   50             continue
               end if
               jx = jx + 4
   60       continue
c**** clean-up loop ****************************************************
            do 80, j = n4 + 1, n
               temp = alpha*x(jx)
               if (temp.ne.0d0) then
                  do 70, i = 1, m
                     y(i) = y(i) + temp*a(i,j)
   70             continue
               end if
               jx = jx + 1
   80       continue

      else
c                                 form  y := alpha*a'*x + y.
c**** u n r o l l   t o   d e p t h   4 ********************************
            m4 = 4*(m/4)
            do 140, j = 1, m4, 4
               temp1 = alpha*x(jx)
               temp2 = alpha*x(jx + 1)
               temp3 = alpha*x(jx + 2)
               temp4 = alpha*x(jx + 3)
               if (temp1.ne.0d0.or.temp2.ne.0d0.or.temp3.ne.0d0.or.
     *             temp4.ne.0d0) then
                  do 130, i = 1, n
                     y(i) = ((((y(i) + temp1*a(j, i))
     *                        + temp2*a(j + 1, i))
     *                        + temp3*a(j + 2, i))
     *                        + temp4*a(j + 3, i))
  130             continue
               end if
               jx = jx + 4
  140       continue
c**** clean-up loop ****************************************************
            do 160, j = m4 + 1, m
               temp = alpha*x(jx)
               if (temp.ne.0d0) then
                  do 150, i = 1, n
                     y(i) = y(i) + temp*a(j, i)
  150             continue
               end if
               jx = jx + 1
  160       continue

      end if
c                                 end of dgemv.
      end

      subroutine dtrsv (uplo, trans, diag, n, a, lda, x)
c-----------------------------------------------------------------------
c  dtrsv  solves one of the systems of equations
c     a*x = b,   or   a'*x = b,
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c-----------------------------------------------------------------------
      implicit none

      integer lda, n, i, j

      character*1 diag, trans, uplo

      double precision a(lda,*), x(*), temp
    
      logical nounit
c-----------------------------------------------------------------------
      if (n.eq.0) return

      nounit = diag.eq.'n'

      if (trans.eq.'n') then
c                                 form  x := inv(a)*x.
         if (uplo.eq.'u') then

               do 20, j = n, 1, -1
                  if (x(j).ne.0d0) then
                     if (nounit) x(j) = x(j)/a(j,j)
                     temp = x(j)
                     do 10, i = j - 1, 1, -1
                        x(i) = x(i) - temp*a(i,j)
   10                continue
                  end if
   20          continue

         else

               do 60, j = 1, n
                  if (x(j).ne.0d0) then
                     if (nounit) x(j) = x(j)/a(j,j)
                     temp = x(j)
                     do 50, i = j + 1, n
                        x(i) = x(i) - temp*a(i,j)
   50                continue
                  end if
   60          continue

         end if

      else

c                                 form  x := inv(a')*x.

         if (uplo.eq.'u') then

               do 100, j = 1, n
                  temp = x(j)
                  do 90, i = 1, j - 1
                     temp = temp - a(i,j)*x(i)
   90             continue
                  if (nounit) temp = temp/a(j,j)
                  x(j) = temp
  100          continue

         else

               do 140, j = n, 1, -1
                  temp = x(j)
                  do 130, i = n, j + 1, -1
                     temp = temp - a(i,j)*x(i)
  130             continue
                  if (nounit) temp = temp/a(j,j)
                  x(j) = temp
  140          continue

         end if

      end if
c                                 end of dtrsv.
      end

      double precision function ddot (n, x, incx, y)
c-----------------------------------------------------------------------
c  ddot returns  ddot = x'y
c----------------------------------------------------------------------
      implicit none

      integer incx, n, i, ix, iy

      double precision x(*), y(*), sum
c----------------------------------------------------------------------
      sum = 0d0

      if (n.gt.0) then

         if (incx.eq.1) then

            do ix = 1, n
               sum = sum + x(ix)*y(ix)
            end do

         else

            iy = 1

            if (incx.gt.0) then
               do 20, ix = 1, 1 + (n - 1)*incx, incx
                  sum = sum + x(ix)*y(iy)
                  iy  = iy  + 1
   20          continue
            else
               ix = 1 - (n - 1)*incx
               do 30, i = 1, n
                  sum = sum + x(ix)*y(iy)
                  ix  = ix  + incx
                  iy  = iy  + 1
   30          continue
            end if
         end if
      end if

      ddot = sum
c                                 end of ddot
      end

      subroutine dger (m, n, alpha, x, y, a, lda)
c-----------------------------------------------------------------------
c  dger does a := alpha*x*y' + a,
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and a is an m by n matrix.
c-----------------------------------------------------------------------
      implicit none

      integer lda, m, n, i, j

      double precision alpha, a(lda,*), x(*), y(*), temp
c-----------------------------------------------------------------------
      if (m.eq.0.or.n.eq.0.or.alpha.eq.0d0) return

         do 20, j = 1, n
            if (y(j).ne.0d0) then
               temp = alpha*y(j)
               do 10, i = 1, m
                  a(i,j) = a(i,j) + x(i)*temp
   10          continue
            end if
   20    continue

c                                 end of dger.
      end

      double precision function dnrm2 (n, x, incx)
c-----------------------------------------------------------------------
c  dnrm2 returns the euclidean norm of a vector via the function
c  name, so that dnrm2 := sqrt(x'*x)
c-----------------------------------------------------------------------
      implicit none

      integer incx, n

      double precision x(*), norm, scale, ssq, snorm
c----------------------------------------------------------------------
      if (n.lt.1) then
         norm  = 0d0
      else if (n.eq.1) then
         norm  = abs(x(1))
      else
         scale = 0d0
         ssq   = 1d0
         call sssq (n, x, incx, scale, ssq)
         norm  = snorm(scale, ssq)
      end if

      dnrm2  = norm
c                                 end of dnrm2
      end

      subroutine sgeqrp (pivot,m,n,a,lda,zeta,perm,work)
c-----------------------------------------------------------------------
c  sgeqrp finds a qr factorization of the real m by n matrix a,
c  incorporating column interchanges, so that a is reduced to upper
c  triangular form by means of orthogonal transformations and column
c  permutations.
c-----------------------------------------------------------------------
      implicit none

      integer lda, m, n, perm(*), j, jmax, k, la

      character pivot*1

      double precision a(lda,*), work(*), zeta(*), eps, maxnrm, norm, 
     *                 temp, tol, dnrm2

      external dnrm2

      double precision wmach
      common/ cstmch /wmach(10)
c-----------------------------------------------------------------------
      if (min(m,n).eq.0) call errdbg ('sgeqrp')

      eps = wmach(3)
      do 20 j = 1, n
         work(j) = dnrm2 (m,a(1,j),1)
         work(j+n) = work(j)
   20 continue

      la = lda
      do 120 k = 1, min(m,n)

         maxnrm = 0d0
         jmax = k
         if (pivot.eq.'c') then
            do 40 j = k, n
               if (work(j+n).gt.maxnrm) then
                  maxnrm = work(j+n)
                  jmax = j
               end if
   40       continue
         else
            do 60 j = k, n
               if (work(j).gt.0d0) then
                  if (k.le.1) then
                     jmax = j
                     go to 80
                  else if ((work(j+n)/work(j)).gt.maxnrm) then
                     maxnrm = work(j+n)/work(j)
                     jmax = j
                  end if
               end if
   60       continue
   80       continue
         end if
         perm(k) = jmax
         if (jmax.gt.k) then
            call dswap(m,a(1,k),a(1,jmax),1)
            temp = work(k)
            work(k) = work(jmax)
            work(jmax) = temp
            work(jmax+n) = work(k+n)
         end if
         tol = eps*work(k)
         if (k.lt.m) then

            call sgrfg (m-k,a(k,k),a(k+1,k),1,tol,zeta(k))
            if (k.lt.n) then
               if (zeta(k).gt.0d0) then
                  if ((k+1).eq.n) la = m - k + 1

                  temp = a(k,k)
                  a(k,k) = zeta(k)

                  call dgemv ('transpose',m-k+1,n-k,1d0,a(k,k+1),la,
     *                       a(k,k),0d0,zeta(k+1))

                  call dger (m-k+1,n-k,-1d0,a(k,k),zeta(k+1),
     *                      a(k,k+1),la)

                  a(k,k) = temp
               end if

               do 100 j = k + 1, n
                  if (work(j+n).gt.0d0) then
                     temp = abs(a(k,j))/work(j+n)
                     temp = max((1d0+temp)*(1d0-temp),0d0)
                     norm = temp
                     temp = 1d0 + 1d-2*temp*(work(j+n)/work(j))**2
                     if (temp.gt.1d0) then
                        work(j+n) = work(j+n)*sqrt(norm)
                     else
                        work(j+n) = dnrm2 (m-k,a(k+1,j),1)
                     end if
                  end if
  100          continue
            end if
         end if
  120 continue

c     set the final  zeta  when  m.le.n.

      if (m.le.n) zeta(m) = 0d0
c                                 end of sgeqrp
      end

      subroutine sgeqr (m,n,a,lda,zeta)
c-----------------------------------------------------------------------
c  sgeqr   finds  the  qr factorization  of the real  m by n,  m .ge. n,
c  matrix a,  so that  a is reduced to upper triangular form by means of
c  orthogonal transformations.
c-----------------------------------------------------------------------
      implicit none

      integer lda, m, n, k, la

      double precision a(lda,*), zeta(*), temp
c----------------------------------------------------------------------
      if (n.eq.0) call errdbg ('sgeqr')

      la = lda
      do 20 k = 1, min(m-1,n)

         call sgrfg (m-k,a(k,k),a(k+1,k),1,0d0,zeta(k))
         if (zeta(k).gt.0d0 .and. k.lt.n) then
            if (k+1.eq.n) la = m - k + 1

            temp = a(k,k)
            a(k,k) = zeta(k)

            call dgemv ('transpose',m-k+1,n-k,1d0,a(k,k+1),la,a(k,k),
     *                 0d0,zeta(k+1))

            call dger (m-k+1,n-k,-1d0,a(k,k),zeta(k+1),a(k,k+1),la)

            a(k,k) = temp
         end if
   20 continue

c     set the final  zeta  when  m.eq.n.

      if (m.eq.n) zeta(n) = 0d0
c                                 end of sgeqr
      end

      double precision function dlantr (norm,uplo,diag,m,n,a,lda,work)
c-----------------------------------------------------------------------
c  dlantr  returns the value of the one norm,  or the frobenius norm, or
c  the  infinity norm,  or the  element of  largest absolute value  of a
c  trapezoidal or triangular matrix a.

c     dlantr = (max(abs(a(i,j))), norm = 'm' or 'm'
c              (
c              (norm1(a),         norm = '1', 'o' or 'o'
c              (
c              (normi(a),         norm = 'i' or 'i'
c              (
c              (normf(a),         norm = 'f', 'f', 'e' or 'e'

c  where  norm1  denotes the  one norm of a matrix (maximum column sum),
c  normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
c  normf  denotes the  frobenius norm of a matrix (square root of sum of
c  squares).  note that  max(abs(a(i,j)))  is not a  matrix norm.

c  -- lapack auxiliary routine
c-----------------------------------------------------------------------
      implicit none

      integer lda, m, n, i, j

      character diag*1, norm*1, uplo*1

      double precision a(lda,*), work(*), scale, sum, value

      logical udiag
c----------------------------------------------------------------------
      if (min(m,n).eq.0) then
         value = 0d0
      else if (norm.eq.'m') then

c        find max(abs(a(i,j))).

         if (diag.eq.'u') then
            value = 1d0
            if (uplo.eq.'u') then
               do 40 j = 1, n
                  do 20 i = 1, min(m,j-1)
                     value = max(value,abs(a(i,j)))
   20             continue
   40          continue
            else
               do 80 j = 1, n
                  do 60 i = j + 1, m
                     value = max(value,abs(a(i,j)))
   60             continue
   80          continue
            end if
         else
            value = 0d0
            if (uplo.eq.'u') then
               do 120 j = 1, n
                  do 100 i = 1, min(m,j)
                     value = max(value,abs(a(i,j)))
  100             continue
  120          continue
            else
               do 160 j = 1, n
                  do 140 i = j, m
                     value = max(value,abs(a(i,j)))
  140             continue
  160          continue
            end if
         end if
      else if (norm.eq.'o' .or. norm.eq.'1') then

c        find norm1(a).

         value = 0d0
         udiag = diag.eq.'u'
         if (uplo.eq.'u') then
            do 220 j = 1, n
               if ((udiag) .and. (j.le.m)) then
                  sum = 1d0
                  do 180 i = 1, j - 1
                     sum = sum + abs(a(i,j))
  180             continue
               else
                  sum = 0d0
                  do 200 i = 1, min(m,j)
                     sum = sum + abs(a(i,j))
  200             continue
               end if
               value = max(value,sum)
  220       continue
         else
            do 280 j = 1, n
               if (udiag) then
                  sum = 1d0
                  do 240 i = j + 1, m
                     sum = sum + abs(a(i,j))
  240             continue
               else
                  sum = 0d0
                  do 260 i = j, m
                     sum = sum + abs(a(i,j))
  260             continue
               end if
               value = max(value,sum)
  280       continue
         end if
      else if (norm.eq.'i') then

c        find normi(a).

         if (uplo.eq.'u') then
            if (diag.eq.'u') then
               do 300 i = 1, m
                  work(i) = 1d0
  300          continue
               do 340 j = 1, n
                  do 320 i = 1, min(m,j-1)
                     work(i) = work(i) + abs(a(i,j))
  320             continue
  340          continue
            else
               do 360 i = 1, m
                  work(i) = 0d0
  360          continue
               do 400 j = 1, n
                  do 380 i = 1, min(m,j)
                     work(i) = work(i) + abs(a(i,j))
  380             continue
  400          continue
            end if
         else
            if (diag.eq.'u') then
               do 420 i = 1, n
                  work(i) = 1d0
  420          continue
               do 440 i = n + 1, m
                  work(i) = 0d0
  440          continue
               do 480 j = 1, n
                  do 460 i = j + 1, m
                     work(i) = work(i) + abs(a(i,j))
  460             continue
  480          continue
            else
               do 500 i = 1, m
                  work(i) = 0d0
  500          continue
               do 540 j = 1, n
                  do 520 i = j, m
                     work(i) = work(i) + abs(a(i,j))
  520             continue
  540          continue
            end if
         end if
         value = 0d0
         do 560 i = 1, m
            value = max(value,work(i))
  560    continue
      else if (norm.eq.'f' .or. norm.eq.'e') then

c        find normf(a).

         if (uplo.eq.'u') then
            if (diag.eq.'u') then
               scale = 1d0
               sum = min(m,n)
               do 580 j = 2, n
                  call sssq (min(m,j-1),a(1,j),1,scale,sum)
  580          continue
            else
               scale = 0d0
               sum = 1d0
               do 600 j = 1, n
                  call sssq (min(m,j),a(1,j),1,scale,sum)
  600          continue
            end if
         else
            if (diag.eq.'u') then
               scale = 1d0
               sum = min(m,n)
               do 620 j = 1, n
                  call sssq (m-j,a(min(m,j+1),j),1,scale,sum)
  620          continue
            else
               scale = 0d0
               sum = 1d0
               do 640 j = 1, n
                  call sssq (m-j+1,a(j,j),1,scale,sum)
  640          continue
            end if
         end if
         value = scale*sqrt(sum)
      end if

      dlantr = value
c                                 end of dlantr
      end

      subroutine scsg (t, c, s)
c-----------------------------------------------------------------------
c  scsg  returns values c and s such that
c     c = cos(theta),   s = sin(theta)
c  for a given value of
c     t = tan(theta).
c  c is always non-negative and s has the same sign as t, so that
c     c = 1.0/sqrt(1.0 + t**2),   s = t/sqrt(1.0 + t**2).
c-----------------------------------------------------------------------
      implicit none

      double precision c, s, t, abst, eps, reps, rrteps, rteps

      logical first

      save first, eps, reps, rteps, rrteps
      data first/ .true. /

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if (first) then
         first  = .false.
         eps    =  wmach(3)
         reps   =  1/eps
         rteps  =  sqrt(eps)
         rrteps =  1/rteps
      end if
c
      abst = abs(t)
      if (abst.lt.rteps) then
         c = 1d0
         s = t
      else if (abst.gt.rrteps) then
         c = 1/abst
         s = sign(1d0, t)
      else
         c = 1/sqrt(1 + abst**2)
         s = c*t
      end if
c                                 end of scsg
      end

      subroutine sssq (n, x, incx, scale, sumsq)
c-----------------------------------------------------------------------
c  sssq  returns the values scl and smsq such that
c     (scl**2)*smsq = x(1)**2 +...+ x(n)**2 + (scale**2)*sumsq,
c  where x(i) = x(1 + (i - 1)*incx). the value of sumsq is assumed
c  to be at least unity and the value of smsq will then satisfy
c     1.0 .le. smsq .le. (sumsq + n) .
c  scale is assumed to be non-negative and scl returns the value
c     scl = max(scale, abs(x(i))) .
c-----------------------------------------------------------------------
      implicit none
 
      integer incx, ix, n

      double precision x(*), scale, sumsq, absxi
c----------------------------------------------------------------------
      if (n.gt.0) then
         do 10, ix = 1, 1 + (n - 1)*incx, incx
            if (x(ix).ne.0d0) then
               absxi = abs(x(ix))
               if (scale.lt.absxi) then
                  sumsq = 1     + sumsq*(scale/absxi)**2
                  scale = absxi
               else
                  sumsq = sumsq +       (absxi/scale)**2
               end if
            end if
   10    continue
      end if
c                                 end of sssq
      end

      subroutine scond (n, x, incx, xmax, xmin)
c-----------------------------------------------------------------------
c  scond  returns the values xmax and xmin given by
c     xmax = max(abs(x(i))),   xmin = min(abs(x(i))).
c             i                              i
c  if n is less than unity then xmax and xmin are returned as zero.
c----------------------------------------------------------------------
      implicit none

      double precision xmax, xmin, x(*)

      integer incx, n, ix
c----------------------------------------------------------------------
      if (n.lt.1) then
         xmax = 0d0
         xmin = 0d0
      else
         xmax = abs(x(1))
         xmin = xmax
         do 10 ix = 1 + incx, 1 + (n - 1)*incx, incx
            xmax = max(xmax, abs(x(ix)))
            xmin = min(xmin, abs(x(ix)))
   10    continue
      end if
c                                 end of scond
      end

      integer function isrank (n, x, incx, tol)
c----------------------------------------------------------------------
c  isrank finds the first element of the n element vector x for which
c     abs(x(k)).le.tol*max(abs(x(1)), ..., abs(x(k - 1)))
c  and returns the value (k - 1) in the function name isrank. if no
c  such k exists then isrank is returned as n.

c  if tol is supplied as less than zero then the value epsmch, where
c  epsmch is the relative machine precision, is used in place of tol.
c----------------------------------------------------------------------
      implicit none

      integer incx, n, ix, k

      double precision tol, x(*), tl, xmax          

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      k = 0
      if (n.ge.1) then

         ix = 1
         if (tol.lt.0d0) then
            tl = wmach(3)
         else
            tl = tol
         end if

         xmax = abs(x(ix))

c+       while(k.lt.n)loop
   10    if   (k.lt.n) then
            if (abs(x(ix)).le.tl*xmax) go to 20
            xmax = max(xmax, abs(x(ix)))
            k    = k  + 1
            ix   = ix + incx
            go to 10
         end if
c+       end while

      end if

   20 isrank = k
c                                 end of isrank
      end

      subroutine sutsqr (side, n, k1, k2, c, s, a, lda)
c-----------------------------------------------------------------------
c  sutsqr does the transformation
c     r := p*u*q'  when  side = 'l' or 'l'  ( left-hand side)
c     r := q*u*p'  when  side = 'r' or 'r'  (right-hand side),
c  where u and r are n by n  upper triangular  matrices, p is an
c  orthogonal matrix, consisting of a given sequence of plane rotations
c  to be  applied in planes  k1 to k2,  and q is a unitary matrix
c  consisting of a sequence of plane rotations, applied in planes  k1 to
c  k2, chosen to make r upper triangular.
c-----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, n, i, i1, j

      character side*1

      double precision a(lda,*), c(*), s(*), aij, ctemp, fill, stemp, 
     *                 temp     
c----------------------------------------------------------------------
      if ((min(n, k1).lt.1).or.(k2.le.k1).or.
     *   (k2.gt.n))return
      if (side.eq.'l') then

         do 20 j = k1 + 1, n

            aij = a(k1, j)
            do 10 i = k1, min(j - 1, k2 - 1)
               a(i,j) = s(i)*a(i + 1, j) + c(i)*aij
               aij = c(i)*a(i + 1, j) - s(i)*aij
   10       continue
            a(i,j) = aij
   20    continue

         do 40 j = k1, k2 - 1


            fill = -s(j)*a(j,j)
            a(j,j) = c(j)*a(j,j)

            call srotgc(a(j + 1, j+1), fill, ctemp, stemp)
            c(j) = ctemp
            s(j) = -stemp
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               stemp = -stemp
               do 30 i = 1, j
                  temp = a(i, j+1)
                  a(i, j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   30          continue
            end if
   40    continue
      else if (side.eq.'r') then

         do 60 j = k2 - 1, k1, -1

            if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
               ctemp = c(j)
               stemp = s(j)
               do 50 i = 1, j
                  temp = a(i, j+1)
                  a(i, j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   50          continue

               fill = s(j)*a(j + 1, j+1)
               a(j + 1, j+1) = c(j)*a(j + 1, j+1)

               call srotgc(a(j,j), fill, c(j), s(j))
            end if
   60    continue


         do 80 j = n, k1 + 1, -1
            i1 = min(k2, j)
            aij = a(i1, j)
            do 70 i = i1 - 1, k1, -1
               temp = a(i,j)
               a(i + 1, j) = c(i)*aij - s(i)*temp
               aij = s(i)*aij + c(i)*temp
   70       continue
            a(k1, j) = aij
   80    continue
      end if
c                                 end of sutsqr
      end

      subroutine sdscl (n, d, incd, x, incx)
c-----------------------------------------------------------------------
c  sdscl does x := diag(d)*x
c-----------------------------------------------------------------------
      implicit none

      integer incd, incx, n, i, id, ix

      double precision d(*), x(*)
c-----------------------------------------------------------------------
      if (n.gt.0) then
         if ((incd.eq.0).and.(incx.ne.0)) then
            call dscal (n, d(1), x, abs(incx))
         else if ((incd.eq.incx).and.(incd.gt.0)) then
            do 10, id = 1, 1 + (n - 1)*incd, incd
               x(id) = d(id)*x(id)
   10       continue
         else
            if (incx.ge.0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incd.gt.0) then
               do 20, id = 1, 1 + (n - 1)*incd, incd
                  x(ix) = d(id)*x(ix)
                  ix      = ix              + incx
   20          continue
            else
               id = 1 - (n - 1)*incd
               do 30, i = 1, n
                  x(ix) = d(id)*x(ix)
                  id      = id              + incd
                  ix      = ix              + incx
   30          continue
            end if
         end if
      end if
c                                 end of sdscl
      end

      subroutine ssrotg (pivot, direct, n, alpha, x, incx, c, s)
c-----------------------------------------------------------------------
c  ssrotg generates the parameters of an orthogonal matrix p such that

c     when   pivot = 'f' or 'f'   and   direct = 'f' or 'f'
c     or     pivot = 'v' or 'v'   and   direct = 'b' or 'b'

c        p*(alpha) = (beta),
c          (  x  )   (  0 )

c     when   pivot = 'f' or 'f'   and   direct = 'b' or 'b'
c     or     pivot = 'v' or 'v'   and   direct = 'f' or 'f'

c        p*(  x  ) = (  0 ),
c          (alpha) = (beta)

c  where alpha is a scalar and x is an n element vector.
c-----------------------------------------------------------------------
      implicit none 

      integer incx, n, i, ix

      character direct*1, pivot*1

      double precision c(*), s(*), x(*), alpha         
c----------------------------------------------------------------------
      if (n.gt.0) then
         if (direct.eq.'b') then
            ix = 1 + (n - 1)*incx
            if (pivot.eq.'v') then
               do 10, i = n, 2, -1
                  call srotgc(x(ix - incx), x(ix), c(i), s(i))
                  ix = ix - incx
   10          continue
               call srotgc(alpha, x(ix), c(1), s(1))
            else if (pivot.eq.'f') then

               do 20, i = n, 1, -1
                  call srotgc(alpha, x(ix), c(i), s(i))
                  s(i)  = -s(i)
                  x(ix) = -x(ix)
                  ix      =  ix      - incx
   20          continue
            end if
         else if (direct.eq.'f') then
            ix = 1
            if (pivot.eq.'v') then

               do 30, i = 1, n - 1
                  call srotgc(x(ix + incx), x(ix), c(i), s(i))
                  s(i)  = -s(i)
                  x(ix) = -x(ix)
                  ix      =  ix      + incx
   30          continue
               call srotgc(alpha, x(ix), c(n), s(n))
               s(n)  = -s(n)
               x(ix) = -x(ix)
            else if (pivot.eq.'f') then
               do 40, i = 1, n
                  call srotgc(alpha, x(ix), c(i), s(i))
                  ix = ix + incx
   40          continue
            end if
         end if
      end if
c                                 end of ssrotg
      end

      double precision function sdiv (a, b, fail)
c-----------------------------------------------------------------------
c  sdiv  returns the value div given by

c     div = (a/b                 if a/b does not overflow,
c           (
c           (0.0                 if a.eq.0.0,
c           (
c           (sign(a/b)*flmax   if a .ne. 0.0  and a/b would overflow,

c  where  flmax  is a large value, via the function name. in addition if
c  a/b would overflow then  fail is returned as true, otherwise  fail is
c  returned as false.

c  note that when  a and b  are both zero, fail is returned as true, but
c  div  is returned as  0.0. in all other cases of overflow  div is such
c  that  abs(div) = flmax.

c  when  b = 0  then  sign(a/b)  is taken as  sign(a).
c-----------------------------------------------------------------------
      implicit none
               
      logical fail, first

      double precision a, b, absb, div, flmax, flmin

      double precision wmach
      common/ cstmch /wmach(10)

      save                  first, flmin, flmax
      data                  first/ .true. /
c----------------------------------------------------------------------
      if (a.eq.0d0) then
         div = 0d0
         if (b.eq.0d0) then
            fail = .true.
         else
            fail = .false.
         end if
      else

         if (first) then
            first = .false.
            flmin =  wmach(10)
            flmax =  1/flmin
         end if

         if (b.eq.0d0) then
            div  =  sign(flmax, a)
            fail = .true.
         else
            absb = abs(b)
            if (absb.ge.1d0) then
               fail = .false.
               if (abs(a).ge.absb*flmin) then
                  div = a/b
               else
                  div = 0d0
               end if
            else
               if (abs(a).le.absb*flmax) then
                  fail = .false.
                  div  =  a/b
               else
                  fail = .true.
                  div  = flmax
                  if (((a.lt.0d0).and.(b.gt.0d0)).or.
     *                ((a.gt.0d0).and.(b.lt.0d0))    )
     *               div = -div
               end if
            end if
         end if
      end if

      sdiv  = div
c                                 end of sdiv
      end

      subroutine susqr (side, n, k1, k2, c, s, a, lda)
c-----------------------------------------------------------------------
c  susqr  restores an upper spiked matrix  h to upper triangular form by
c  applying a sequence of plane rotations, in planes  k1 up to k2,  from
c  either the left, or the right.
c-----------------------------------------------------------------------
      implicit none 

      integer k1, k2, lda, n, i, j

      character side*1

      double precision a(lda,*), c(*), s(*), aij, ctemp, spike, stemp, 
     *                 temp         
c----------------------------------------------------------------------
      if (min(n, k1).lt.1 .or. k2.le.k1 .or. k2.gt.n) return

      if (side.eq.'l') then

         do 20 j = k1, k2 - 1
            spike = s(j)
            do 10 i = k1, j - 1
               aij = a(i,j)
               a(i,j) = s(i)*spike + c(i)*aij
               spike = c(i)*spike - s(i)*aij
   10       continue

            call srotgc(a(j,j), spike, c(j), s(j))
   20    continue

         do 40 j = k2, n
            temp = a(k2, j)
            do 30 i = k1, k2 - 1
               aij = a(i,j)
               a(i,j) = s(i)*temp + c(i)*aij
               temp = c(i)*temp - s(i)*aij
   30       continue
            a(k2, j) = temp
   40    continue
      else if (side.eq.'r') then
         do 70 j = k2, k1 + 1, -1
            call srotgc(a(j,j), s(j - 1), ctemp, stemp)
            stemp = -stemp
            s(j - 1) = stemp
            c(j - 1) = ctemp
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               do 50 i = j - 1, k1 + 1, -1
                  spike = s(i - 1)
                  s(i - 1) = stemp*a(i,j) + ctemp*spike
                  a(i,j) = ctemp*a(i,j) - stemp*spike
   50          continue
               do 60 i = k1, 1, -1
                  temp = a(i, k1)
                  a(i, k1) = stemp*a(i,j) + ctemp*temp
                  a(i,j) = ctemp*a(i,j) - stemp*temp
   60          continue
            end if
   70    continue
      end if
c                                 end of susqr
      end

      subroutine sgesrc (side, pivot, direct, m, n, k1, k2, c, s, a, 
     *                   lda)
c-----------------------------------------------------------------------
c  sgesrc  does the transformation
c     a := p*a,   when   side = 'l' or 'l'  ( left-hand side)
c     a := a*p',  when   side = 'r' or 'r'  (right-hand side)
c  where a is an m by n matrix and p is an orthogonal matrix, consisting
c  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2.
c-----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, m, n, i, j

      character direct*1, pivot*1, side*1

      double precision a(lda,*), c(*), s(*), aij, ctemp, stemp, temp
          
      logical left, right
c----------------------------------------------------------------------
      left = (side.eq.'l')
      right = (side.eq.'r')
      if ((min(m, n, k1).lt.1).or.(k2.le.k1).or.
     *    ((left).and.(k2.gt.m)).or.
     *    ((right).and.(k2.gt.n)))return
      if (left) then
         if (pivot.eq.'v') then
            if (direct.eq.'f') then
               do 20 j = 1, n
                  aij = a(k1, j)
                  do 10 i = k1, k2 - 1
                     temp = a(i + 1, j)
                     a(i,j) = s(i)*temp + c(i)*aij
                     aij = c(i)*temp - s(i)*aij
   10             continue
                  a(k2, j) = aij
   20          continue
            else if (direct.eq.'b') then
               do 40 j = 1, n
                  aij = a(k2, j)
                  do 30 i = k2 - 1, k1, -1
                     temp = a(i,j)
                     a(i + 1, j) = c(i)*aij - s(i)*temp
                     aij = s(i)*aij + c(i)*temp
   30             continue
                  a(k1, j) = aij
   40          continue
            end if
         else if (pivot.eq.'t') then
            if (direct.eq.'f') then
               do 60 j = 1, n
                  temp = a(k1, j)
                  do 50 i = k1, k2 - 1
                     aij = a(i + 1, j)
                     a(i + 1, j) = c(i)*aij - s(i)*temp
                     temp = s(i)*aij + c(i)*temp
   50             continue
                  a(k1, j) = temp
   60          continue
            else if (direct.eq.'b') then
               do 80 j = 1, n
                  temp = a(k1, j)
                  do 70 i = k2 - 1, k1, -1
                     aij = a(i + 1, j)
                     a(i + 1, j) = c(i)*aij - s(i)*temp
                     temp = s(i)*aij + c(i)*temp
   70             continue
                  a(k1, j) = temp
   80          continue
            end if
         else if (pivot.eq.'b') then
            if (direct.eq.'f') then
               do 100 j = 1, n
                  temp = a(k2, j)
                  do 90 i = k1, k2 - 1
                     aij = a(i,j)
                     a(i,j) = s(i)*temp + c(i)*aij
                     temp = c(i)*temp - s(i)*aij
   90             continue
                  a(k2, j) = temp
  100          continue
            else if (direct.eq.'b') then
               do 120 j = 1, n
                  temp = a(k2, j)
                  do 110 i = k2 - 1, k1, -1
                     aij = a(i,j)
                     a(i,j) = s(i)*temp + c(i)*aij
                     temp = c(i)*temp - s(i)*aij
  110             continue
                  a(k2, j) = temp
  120          continue
            end if
         end if
      else if (right) then
         if (pivot.eq.'v') then
            if (direct.eq.'f') then
               do 140 j = k1, k2 - 1
                  if (c(j).ne.1d0 .or. s(j).ne.0d0) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 130 i = 1, m
                        temp = a(i, j+1)
                        a(i, j+1) = ctemp*temp - stemp*a(i,j)
                        a(i,j) = stemp*temp + ctemp*a(i,j)
  130                continue
                  end if
  140          continue
            else if (direct.eq.'b') then
               do 160 j = k2 - 1, k1, -1
                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 150 i = m, 1, -1
                        temp = a(i, j+1)
                        a(i, j+1) = ctemp*temp - stemp*a(i,j)
                        a(i,j) = stemp*temp + ctemp*a(i,j)
  150                continue
                  end if
  160          continue
            end if
         else if (pivot.eq.'t') then
            if (direct.eq.'f') then
               do 180 j = k1 + 1, k2
                  ctemp = c(j - 1)
                  stemp = s(j - 1)
                  if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
                     do 170 i = 1, m
                        temp = a(i,j)
                        a(i,j) = ctemp*temp - stemp*a(i, k1)
                        a(i, k1) = stemp*temp + ctemp*a(i, k1)
  170                continue
                  end if
  180          continue
            else if (direct.eq.'b') then
               do 200 j = k2, k1 + 1, -1
                  ctemp = c(j - 1)
                  stemp = s(j - 1)
                  if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
                     do 190 i = m, 1, -1
                        temp = a(i,j)
                        a(i,j) = ctemp*temp - stemp*a(i, k1)
                        a(i, k1) = stemp*temp + ctemp*a(i, k1)
  190                continue
                  end if
  200          continue
            end if
         else if (pivot.eq.'b') then
            if (direct.eq.'f') then
               do 220 j = k1, k2 - 1
                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 210 i = 1, m
                        temp = a(i,j)
                        a(i,j) = stemp*a(i, k2) + ctemp*temp
                        a(i, k2) = ctemp*a(i, k2) - stemp*temp
  210                continue
                  end if
  220          continue
            else if (direct.eq.'b') then
               do 240 j = k2 - 1, k1, -1
                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 230 i = m, 1, -1
                        temp = a(i,j)
                        a(i,j) = stemp*a(i, k2) + ctemp*temp
                        a(i, k2) = ctemp*a(i, k2) - stemp*temp
  230                continue
                  end if
  240          continue
            end if
         end if
      end if
c                                 end of sgesrc
      end

      double precision function snorm (scale, ssq)
c----------------------------------------------------------------------
c  snorm returns the value norm given by
c     norm = (scale*sqrt(ssq), scale*sqrt(ssq).lt.flmax
c            (
c            (flmax,             scale*sqrt(ssq) .ge. flmax
c----------------------------------------------------------------------
      implicit none

      double precision scale, ssq, flmax, norm, sqt

      logical first

      double precision wmach
      common/ cstmch /wmach(10)

      save first, flmax
      data first/ .true. /
c----------------------------------------------------------------------
      if (first) then
         first = .false.
         flmax =  1/wmach(10)
      end if

      sqt = sqrt(ssq)
      if (scale.lt.flmax/sqt) then
         norm = scale*sqt
      else
         norm = flmax
      end if

      snorm = norm
c                                 end of snorm.
      end

      subroutine sutsrh (side, n, k1, k2, c, s, a, lda)
c-----------------------------------------------------------------------
c  sutsrh applies a  given sequence  of  plane rotations  to either  the
c  left,  or the right,  of the  n by n  upper triangular matrix  u,  to
c  transform u to an  upper hessenberg matrix. the rotations are applied
c  in planes k1 up to k2.
c-----------------------------------------------------------------------
      implicit none

      integer  k1, k2, lda, n, i, j

      character side*1

      double precision a(lda,*), c(*), s(*), aij, ctemp, stemp, temp
c----------------------------------------------------------------------
      if ((min(n, k1).lt.1).or.(k2.le.k1).or.
     *   (k2.gt.n))return

      if (side.eq.'l') then

         do 20 j = n, k1, -1
            if (j.ge.k2) then
               aij = a(k2, j)
            else

               aij = c(j)*a(j,j)
               s(j) = -s(j)*a(j,j)
            end if
            do 10 i = min(k2, j) - 1, k1, -1
               temp = a(i,j)
               a(i + 1, j) = c(i)*aij - s(i)*temp
               aij = s(i)*aij + c(i)*temp
   10       continue
            a(k1, j) = aij
   20    continue
      else if (side.eq.'r') then

         do 40 j = k1, k2 - 1
            if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
               stemp = s(j)
               ctemp = c(j)
               do 30 i = 1, j
                  temp = a(i, j+1)
                  a(i, j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   30          continue
               s(j) = stemp*a(j + 1, j+1)
               a(j + 1, j+1) = ctemp*a(j + 1, j+1)
            end if
   40    continue
      end if
c                                 end of sutsrh.
      end

      subroutine sgeapr (side, trans, n, perm, k, b, ldb)
c-----------------------------------------------------------------------
c  sgeapr does one of the transformations
c     b := p'*b   or   b := p*b,   where b is an m by k matrix,
c  or
c     b := b*p'   or   b := b*p,   where b is a k by m matrix,
c  p being an m by m permutation matrix of the form
c     p = p(1, index(1))*p(2, index(2))*...*p(n, index(n))
c-----------------------------------------------------------------------
      implicit none 

      logical left, null, right, trnsp

      integer k, ldb, n, i, j, l

      character side*1, trans*1

      double precision perm(*), b(ldb,*), temp
c----------------------------------------------------------------------
      if (min(n, k).eq.0) return
      left = (side.eq.'l')
      right = (side.eq.'r')
      null = (trans.eq.'n')
      trnsp = (trans.eq.'t')
      if (left) then
         if (trnsp) then
            do 20 i = 1, n
               l = perm(i)
               if (l.ne.i) then
                  do 10 j = 1, k
                     temp = b(i,j)
                     b(i,j) = b(l, j)
                     b(l, j) = temp
   10             continue
               end if
   20       continue
         else if (null) then
            do 40 i = n, 1, -1
               l = perm(i)
               if (l.ne.i) then
                  do 30 j = 1, k
                     temp = b(l, j)
                     b(l, j) = b(i,j)
                     b(i,j) = temp
   30             continue
               end if
   40       continue
         end if
      else if (right) then
         if (trnsp) then
            do 60 j = n, 1, -1
               l = perm(j)
               if (l.ne.j) then
                  do 50 i = 1, k
                     temp = b(i,j)
                     b(i,j) = b(i, l)
                     b(i, l) = temp
   50             continue
               end if
   60       continue
         else if (null) then
            do 80 j = 1, n
               l = perm(j)
               if (l.ne.j) then
                  do 70 i = 1, k
                     temp = b(i, l)
                     b(i, l) = b(i,j)
                     b(i,j) = temp
   70             continue
               end if
   80       continue
         end if
      end if
c                                 end of sgeapr.
      end

      subroutine sutsr1 (side,n,k1,k2,s,a,lda)
c-----------------------------------------------------------------------
c  sutsr1 applies a  sequence  of  pairwise interchanges to either  the
c  left,  or the right,  of the  n by n  upper triangular matrix  u,  to
c  transform u to an  upper hessenberg matrix. the interchanges are
c  applied in planes k1 up to k2.
c----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, n, i, j

      character side*1

      double precision a(lda,*), s(*), aij, temp          
c----------------------------------------------------------------------
      if ((min(n,k1).lt.1) .or. (k2.le.k1) .or. (k2.gt.n)) return

      if (side.eq.'l') then

         do 40 j = n, k1, -1
            if (j.ge.k2) then
               aij = a(k2,j)
            else

               aij = 0d0
               s(j) = a(j,j)
            end if
            do 20 i = min(k2,j) - 1, k1, -1
               temp = a(i,j)
               a(i+1,j) = temp
               aij = aij
   20       continue
            a(k1,j) = aij
   40    continue

      else if (side.eq.'r') then

         do 80 j = k1, k2 - 1
            do 60 i = 1, j
               temp = a(i,j+1)
               a(i,j+1) = a(i,j)
               a(i,j) = temp
   60       continue
            s(j) = a(j+1,j+1)
            a(j+1,j+1) = 0d0
   80    continue
      end if
c                                 end of sutsr1
      end

      subroutine sgrfg (n, alpha, x, incx, tol, zeta)
c----------------------------------------------------------------------
c  sgrfg generates details of a generalized householder reflection such
c  that p*(alpha) = (beta),   p'*p = i.
c         (  x  )   (  0 )
c----------------------------------------------------------------------
      implicit none

      integer incx, n

      double precision alpha, tol, zeta, x(*), beta, scale, ssq

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if (n.lt.1) then
         zeta = 0d0
      else if (n.eq.1 .and. x(1).eq.0d0) then
         zeta = 0d0
      else

         if (n.eq.1) then

            if (alpha.eq.0d0) then
               zeta   =  1d0
               alpha  =  abs (x(1))
               x(1) = -sign(1d0, x(1))
            else if (abs(x(1)).le.max(wmach(3)*abs(alpha),tol))then
               zeta   =  0d0
            else
               if (abs(alpha).ge.abs(x(1))) then
                  beta = abs(alpha) *sqrt(1 + (x(1)/alpha)**2)
               else
                  beta = abs(x(1))*sqrt(1 + (alpha/x(1))**2)
               end if
               zeta = sqrt((abs(alpha) + beta)/beta)
               if (alpha.ge.0d0) beta = -beta
               x(1) = -x(1)/(zeta*beta)
               alpha  = beta
            end if
         else

            ssq   = 1d0
            scale = 0d0
            call sssq (n, x, incx, scale, ssq)

            if ((scale.eq.0d0).or.
     *          (scale.le.max(wmach(3)*abs(alpha), tol))) then
               zeta  = 0d0
            else if (alpha.eq.0d0) then
               zeta  = 1d0
               alpha = scale*sqrt(ssq)
               call dscal (n, -1/alpha, x, incx)
            else
               if (scale.lt.abs(alpha)) then
                  beta = abs(alpha)*sqrt(1 + ssq*(scale/alpha)**2)
               else
                  beta = scale       *sqrt(ssq +   (alpha/scale)**2)
               end if
               zeta = sqrt((beta + abs(alpha))/beta)
               if (alpha.gt.0d0) beta = -beta
               call dscal (n, -1/(zeta*beta), x, incx)
               alpha = beta
            end if
         end if
      end if
c                                 end of sgrfg
      end

      subroutine nggqzz (hess,n,k1,k2,c,s,a,lda)
c-----------------------------------------------------------------------
c  nggqzz either applies a given sequence of plane rotations to the
c  right of the n by n reverse lower triangular matrix t, to transform t
c  to a  reverse lower hessenberg matrix h, or restores a reverse lower
c  hessenberg matrix h to reverse lower triangular form t, by applying a
c  sequence of plane rotations from the right.
c----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, n, i, j

      character hess*1

      double precision a(lda,*), c(*), s(*), ctemp, stemp, suph, temp
c----------------------------------------------------------------------
      if (min(n,k1).lt.1 .or. k2.le.k1 .or. k2.gt.n) return

      if (hess.eq.'c') then

         do 40 j = k1, k2 - 1
            if ((c(j).ne.1d0) .or. (s(j).ne.0d0)) then
               stemp = s(j)
               ctemp = c(j)
               s(j) = stemp*a(n-j,j+1)
               a(n-j,j+1) = ctemp*a(n-j,j+1)
               do 20 i = n - j + 1, n
                  temp = a(i,j+1)
                  a(i,j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   20          continue
            end if
   40    continue

      else if (hess.eq.'r') then

         do 80 j = k2 - 1, k1, -1
            suph = s(j)
            call srotgc(a(n-j,j+1),suph,ctemp,stemp)
            stemp = -stemp
            s(j) = stemp
            c(j) = ctemp
            if ((ctemp.ne.1d0) .or. (stemp.ne.0d0)) then
               do 60 i = n - j + 1, n
                  temp = a(i,j+1)
                  a(i,j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   60          continue
            end if
   80    continue
      end if
c                                 end of nggqzz.
      end

      subroutine sutsrs (side, n, k1, k2, c, s, a, lda)
c-----------------------------------------------------------------------
c  sutsrs applies a  given sequence  of  plane rotations  to either  the
c  left,  or the right,  of the  n by n  upper triangular  matrix  u  to
c  transform  u  to an upper spiked matrix. the rotations are applied in
c  planes k1 up to k2.
c----------------------------------------------------------------------
      implicit none 

      integer k1, k2, lda, n, i, j

      character side*1

      double precision a(lda,*), c(*), s(*), aij, ctemp, spike, stemp, 
     *                 temp
c----------------------------------------------------------------------
      if ((min(n, k1).lt.1).or.(k2.le.k1).or.
     *   (k2.gt.n))return

      if (side.eq.'l') then

         do 20 j = n, k2, -1
            temp = a(k2, j)
            do 10 i = k2 - 1, k1, -1
               aij = a(i,j)
               a(i,j) = s(i)*temp + c(i)*aij
               temp = c(i)*temp - s(i)*aij
   10       continue
            a(k2, j) = temp
   20    continue

         do 40 j = k2 - 1, k1, -1
            spike = -s(j)*a(j,j)
            a(j,j) = c(j)*a(j,j)
            do 30 i = j - 1, k1, -1
               aij = a(i,j)
               a(i,j) = s(i)*spike + c(i)*aij
               spike = c(i)*spike - s(i)*aij
   30       continue
            s(j) = spike
   40    continue
      else if (side.eq.'r') then
         do 70 j = k1 + 1, k2
            ctemp = c(j - 1)
            stemp = s(j - 1)
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               do 50 i = 1, k1
                  temp = a(i, k1)
                  a(i, k1) = stemp*a(i,j) + ctemp*temp
                  a(i,j) = ctemp*a(i,j) - stemp*temp
   50          continue
               do 60 i = k1 + 1, j - 1
                  spike = s(i - 1)
                  s(i - 1) = stemp*a(i,j) + ctemp*spike
                  a(i,j) = ctemp*a(i,j) - stemp*spike
   60          continue
               s(j - 1) = stemp*a(j,j)
               a(j,j) = ctemp*a(j,j)
            end if
   70    continue
      end if
c                                 end of sutsrs
      end

      subroutine suhqr (side, n, k1, k2, c, s, a, lda)
c-----------------------------------------------------------------------
c  suhqr restores an upper hessenberg matrix h to upper triangular form
c  by  applying a sequence of  plane rotations  from either the left, or
c  the right. 
c----------------------------------------------------------------------
      implicit none 

      integer k1, k2, lda, n, i, j

      character side*1

      double precision a(lda,*), c(*), s(*), aij, ctemp, stemp, subh, 
     *                 temp            
c----------------------------------------------------------------------
      if (min(n, k1).lt.1 .or. k2.le.k1 .or. k2.gt.n)return

      if (side.eq.'l') then
         do 20 j = k1, n
            aij = a(k1, j)
            do 10 i = k1, min(j, k2) - 1
               temp = a(i + 1, j)
               a(i,j) = s(i)*temp + c(i)*aij
               aij = c(i)*temp - s(i)*aij
   10       continue
            if (j.lt.k2) then
               subh = s(j)
               call srotgc(aij, subh, c(j), s(j))
               a(j,j) = aij
            else
               a(k2, j) = aij
            end if
   20    continue
      else if (side.eq.'r') then
         do j = k2 - 1, k1, -1
            subh = s(j)
            call srotgc(a(j + 1, j+1), subh, ctemp, stemp)
            stemp = -stemp
            s(j) = stemp
            c(j) = ctemp
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               do i = j, 1, -1
                  temp = a(i, j+1)
                  a(i, j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
               end do
            end if
         end do 
      end if
c                                 end of suhqr
      end

      subroutine srotgc (a, b, c, s)
c----------------------------------------------------------------------
c  srotg blas except that c is always
c  returned as non-negative and b is overwritten by the tangent of the
c  angle that defines the plane rotation.
c----------------------------------------------------------------------
      implicit none

      double precision a, b, c, s, t, sdiv

      logical fail

      external sdiv 
c----------------------------------------------------------------------
      if (b.eq.0d0) then
         c  = 1d0
         s  = 0d0
      else
         t  = sdiv (b, a, fail)
         call scsg (t, c, s)
         a  = c*a + s*b
         b  = t
      end if
c                                 end of srotgc
      end

      subroutine dswap (n, x, y, incy)
c-----------------------------------------------------------------------
c dswap does the operations temp := x,   x := y,   y := temp.
c-----------------------------------------------------------------------
      implicit none

      integer incy, n, iy

      double precision x(*), y(*), temp
c----------------------------------------------------------------------
            do iy = 1, 1 + (n - 1)*incy, incy
               temp = x(iy)
               x(iy) = y(iy)
               y(iy) = temp
            end do
c                                 end of dswap
      end

      subroutine sload (n, const, x, incx)
c----------------------------------------------------------------------
c  sload does x = const*e,   e' = (1  1 ... 1).
c----------------------------------------------------------------------
      implicit none

      integer incx, n

      double precision const, x(*)
c----------------------------------------------------------------------
      x(1: 1 + (n- 1)*incx : incx) = const
c                                 end of sload
      end

      subroutine smcopy (matrix, m, n, a, lda, b, ldb)
c-----------------------------------------------------------------------
c  smcopy copies the m by n matrix a into the m by n matrix b.
c-----------------------------------------------------------------------
      implicit none

      character matrix*1

      integer m, n, lda, ldb, i, j

      double precision a(lda,*), b(ldb,*)
c----------------------------------------------------------------------
      if (matrix.eq.'g') then

         do j = 1, n
            b(1:m,j) = a(1:m,j)
         end do

      else if (matrix.eq.'u') then

         do j = 1, n
            i = min(m, j)
            b(1:i,j) = a(1:i,j)
         end do

      else if (matrix.eq.'l') then

         do j = 1, min(m, n)
            b(j:m,j) = a(j:m,j)
         end do

      end if
c                                 end of smcopy
      end

      subroutine smload (matrix, m, n, const, diag, a, lda)
c-----------------------------------------------------------------------
c  smload forms the m by n matrix a(i,j) = diag i.eq.j, const i.ne.j.
c-----------------------------------------------------------------------
      implicit none

      character matrix*1

      integer lda, m, n, i, j

      double precision const, diag, a(lda,*)
c-----------------------------------------------------------------------
      if (matrix.eq.'g') then

         do j = 1, n
            a(1:m,j) = const
         end do

      else if (matrix.eq.'u') then

         do j = 1, n
            i = min(m, j)
            a(1:i,j) = const
         end do

      else if (matrix.eq.'l') then

         do j = 1, min(m, n)
            a(j:m, j) = const
         end do

      end if

      if (const.ne.diag) then

         do i = 1, min(m, n)
            a(i, i) = diag
         end do

      end if
c                                 end of smload
      end

      subroutine iload (n, const, x, incx)
c-----------------------------------------------------------------------
c  iload  does x = const*e,   e' = (1  1 ... 1).
c-----------------------------------------------------------------------
      implicit none

      integer const, incx, n, x(*), ix
c----------------------------------------------------------------------
      ix = 1 + (n - 1)*incx

      x(1:ix:incx) = const
c                                 end of iload
      end

      subroutine dscal (n, alpha, x, incx)
c-----------------------------------------------------------------------
c  dscal does x := alpha*x
c----------------------------------------------------------------------
      implicit none

      integer incx, n, ix

      double precision  x(*), alpha
c----------------------------------------------------------------------
      ix = 1 + (n - 1)*incx

      x(1:ix:incx) = alpha * x(1:ix:incx)
c                                 end of dscal
      end

      subroutine dcopy (n, x, incx, y, incy)
c-----------------------------------------------------------------------
c  dcopy does y := x
c-----------------------------------------------------------------------
      implicit none

      integer incx, incy, n, ix, iy

      double precision x(*), y(*)
c-----------------------------------------------------------------------
      iy = 1 + (n - 1)*incy
      ix = 1 + (n - 1)*incx

      if (incx.eq.incy .and. incy.gt.0) then

         y(1:iy:incy) = x(1:iy:incy)

      else

         if (incx.ge.0) then
c                                 x ascending
            if (incy.gt.0) then 
c                                 y ascending
               y(1:iy:incy) = x(1:ix:incx)

            else 
c                                 y descending
               y(iy:1:incy) = x(1:ix:incx)

            end if

         else
c                                 x descending
            if (incy.gt.0) then 
c                                 y ascending
               y(1:iy:incy) = x(ix:1:incx)

            else 
c                                 y descending
               y(iy:1:incy) = x(ix:1:incx)

            end if

         end if

      end if
c                                 end of dcopy
      end

      subroutine sscmv (n, alpha, x, y)
c-----------------------------------------------------------------------
c  sscmv does y := alpha*x
c-----------------------------------------------------------------------
      implicit none

      double precision alpha, x(*), y(*)

      integer n
c----------------------------------------------------------------------
      if (alpha.eq.0d0) then

         call sload (n, 0d0, y, 1)

      else

         y(1:n) = alpha*x(1:n)

      end if
c                                 end of sscmv
      end

