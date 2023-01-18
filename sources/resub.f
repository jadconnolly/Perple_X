c routines only called by vertex/meemum, could be combined with nlib.f

      subroutine lpopt0 (idead)
c-----------------------------------------------------------------------
c lpopt0 - calls lp minimization after a call to initlp. lpopt0
c does the minimization, writes error messages if necessary.

c this is an utterly stupid formulation of the lp problem because i modified
c the lp code to impose the implicit constraint that phase amounts were between
c 0 and 1, but did not impose the constraint the sum of the amounts must be
c be 1 (which would be unwise in any case). this requires that compositions and
c g's must be normalized to the total number of moles (and therefore lots of 
c extra bookkeeping). 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,liw,lw,k,idead,inc,lphct,jter, lpprob

      parameter (liw=2*k1+3,lw=2*(k5+1)**2+7*k1+5*k5)  

      double precision ax(k5),x(k1),w(lw),oldt,oldp,gtot,
     *                 tol,oldx,clamda(k1+k5)

      integer iw(liw)

      logical quit, abort

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision g
      common/ cst2 /g(k1)

      integer jphct,istart
      common/ cst111 /jphct,istart

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision bl,bu
      common/ cstbup /bl(k1+k5),bu(k1+k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer tphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      logical abort1
      common/ cstabo /abort1

      double precision wmach
      common/ cstmch /wmach(10)

      save ax, x, clamda, w, iw
c-----------------------------------------------------------------------
      idegen = 0
      jdegen = 0
c                                 degeneracy test
      do k = 1, icp 
         if (b(k).eq.0d0) then 
            idegen = idegen + 1
            idg(idegen) = k
         else 
            jdegen = jdegen + 1
            jdg(jdegen) = k
         end if
      end do

      inc = istct - 1

      oldt = t
      oldp = p
      oldx = xco2
c                                logarithmic_p option
      if (lopt(14)) p = 1d1**p
c                                logarithmic_X option
      if (lopt(37)) xco2 = 1d1**xco2
c                                t_stop option
      if (t.lt.nopt(12)) t = nopt(12)

      if (lopt(61)) call begtim (1)

      call gall

      if (lopt(61)) call endtim (1,.false.,'Static GALL ')

      do k = 1, jphct
         c(k) = g(k+inc)/ctot(k+inc)
      end do
c                                 load the adaptive refinement cpd g's
      do k = 1, jpoint
         g2(k) = c(k)
      end do 
c                                 load the bulk into the constraint array
      bl(jphct+1:jphct+icp) = b(1:icp)
      bu(jphct+1:jphct+icp) = b(1:icp)

      lpprob = 2
      tol = wmach(4)

      if (lopt(61)) call begtim (13)

      call lpsol (jphct,hcp,a,k5,bl,bu,c,is,x,jter,gtot,ax,
     *            clamda,iw,liw,w,lw,idead,istart,tol,lpprob)
c                                 set istart according to static_LP_start
      if (istart.ne.0) istart = iopt(39)

      if (lopt(61)) call endtim (13,.false.,'Static optimization ')

      if (idead.gt.0) then
c                                 look for severe errors                                            
         call lpwarn (idead,'LPOPT ')
c                                 on severe error do a cold start.
c                                 necessary?
         istart = 0

      else if (isoct.eq.0) then 
c                                 no refinement, find the answer
         call yclos0 (x,is,jphct) 
c                                 final processing, .true. indicates static
         call rebulk (abort,.true.)

      else
c                                 save lphct to recover static solution if
c                                 no refinement 
         lphct = jphct 
c                                 find discretization points
c                                 for refinement
c        if (lopt(28)) call begtim (3)

         call yclos1 (x,clamda,jphct,quit)

c        if (lopt(28)) call endtim (3,.true.,'Static YCLOS1 ')

c                                 returns quit if nothing to refine
         if (quit) then 
c                                 final processing, .true. indicates static
            call rebulk (abort,.true.)

         else
c                                 initialize refinement point pointers
            do i = 1, ipoint
               hkp(i) = 0 
            end do 

c            if (lopt(28)) call begtim (4)
c                                 reoptimize with refinement
            call reopt (idead,gtot)

c            if (lopt(28)) call endtim (4,.true.,'Dynamic optimization ')

c                                 final processing, .false. indicates dynamic
             if (idead.eq.0) then 

               call rebulk (abort,.false.)

               if (abort) then
c                                 bad solution (lagged speciation) identified
c                                 in avrger
                  call lpwarn (102,'LPOPT0')
                  if (iopt(22).gt.2) idead = 102

               end if 

               if (lopt(32)) then 
c                                 if lagged speciation
                  if (abort1) idead = 104

               end if

            else if (idead.eq.-1) then
c                                 hail mary
               jphct = lphct
               idead = 0

               call yclos0 (x,is,jphct) 
               call rebulk (abort,.true.)

            end if 

         end if 

      end if 

      t = oldt
      p = oldp
      xco2 = oldx

      end 

      subroutine reopt (idead,ogtot)
c-----------------------------------------------------------------------
c reopt - given the results of an initial optimization for lpopt, reopt
c iteratively refines the solution by generating pseudocompounds in the
c neighborhood of the initial optimization.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw, lw, iter, idead, jstart, opt, i, j,
     *        idead1, jter, iprint, lpprob, xphct

      logical quit

      parameter (liw=2*k21+3,lw=2*(k5+1)**2+7*k21+5*k5)

      double precision ax(k5), clamda(k21+k5), w(lw), tot(k5), gtot,
     *                 ogtot, bl(k21+k5), bu(k21+k5), d2g(3), tol

      integer is(k21+k5), iw(liw)

      double precision wmach
      common/ cstmch /wmach(10)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision x
      common/ scrtch /x(k21)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer xis
      double precision xa,b,xc
      common/ cst313 /xa(k5,k1),b(k5),xc(k1),xis(k1+k5)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
c-----------------------------------------------------------------------
c                                 the pseudocompounds to be refined
c                                 are identified in jdv(1..npt)
      quit = .false.
      opt = npt
      idead1 = 0
      d2g(1) = ogtot

      jphct = jpoint
c                                 global composition coordinate counter
      zcoct = 0
c                                 --------------------------------------
c                                 generate pseudo compounds for the first 
c                                 iteration from static arrays
      call resub (1)

      if (jphct.eq.jpoint) then
c                                 if nothing to refine, set idead 
c                                 to recover previous solution,
c                                 DEBUG DEBUG set to error 102 because
c                                 likely failed aqlagd
         idead = 102
         return

      end if
c                                  initialization
      x(1:jphct) = 0d0
      bl(1:jphct) = 0d0
      bu(1:jphct) = 1d0
c                                  iopt(38) = 0, cold start, amounts 
c                                  and state are not set.
      if (iopt(38).eq.1) then 
c                                  iopt(38) = 1, set amounts, but not state.
c                                  stoichiometric phases
         do i = 1, lcpt
            x(ldv(i)) = lamt(i)
         end do
c                                  solution compositions
         do i = 1, lspt
            if (lsdv(i).eq.0) cycle
            x(lsdv(i)) = lsamt(i)
         end do

      else if (iopt(38).eq.2) then 
c                                  iopt(38) = 2, set amounts and state.
         is(jpoint+1:jphct) = 1
         is(1:jpoint) = xis(1:jpoint)
c                                  stoichiometric phases
         do i = 1, lcpt
            x(ldv(i)) = lamt(i)
         end do
c                                  solution compositions
         do i = 1, lspt
            if (lsdv(i).eq.0) cycle
            x(lsdv(i)) = lsamt(i)
            is(lsdv(i)) = lsst(i)
         end do
c                                  set constraint states
         is(jphct+1:jphct+icp) = 3

      end if

      iter = 2

      do
c                                 iter is incremented before the operations,
c                                 i.e., on the nth iteration, iter is n+1
         iter = iter + 1
c                                 set quit flag
         if (iter.gt.iopt(20)) quit = .true.
c                                 cold 0/warm 1 start
         if (iopt(38).eq.2) then
            jstart = 1
         else 
            jstart = 0
         end if

         iprint = 0
         lpprob = 2
         tol = wmach(4)
c                                  set constraint states
         is(jphct+1:jphct+icp) = 3
c                                  and bounds
         bl(jphct+1:jphct+icp) = b(1:icp)
         bu(jphct+1:jphct+icp) = b(1:icp)

         if (lopt(61)) call begtim (14)

         call lpsol (jphct,icp,cp2,k5,bl,bu,g2,is,x,jter,gtot,ax,
     *               clamda,iw,liw,w,lw,idead,jstart,tol,lpprob)

         if (lopt(61)) call endtim (14,.false.,'dynamic optimization ')

c                                 warn if severe error
         if (idead.gt.0) then

            if (idead.ne.3.or.quit) then 

               call lpwarn (idead,'REOPT')
               exit 

            end if 
c                                  the logic here, is that the composition
c                                  can't really be bad since we got through the 
c                                  initial minimization. do a mass balance check
c                                  just in case:
            tot(1:icp) = b(1:icp)

            do i = 1, jphct

               if (is(i).eq.1) cycle

               do j = 1, icp 
                  tot(j) = tot(j) - x(i)*cp2(j,i)
               end do

            end do

            do i = 1, icp

               if (dabs(tot(i)).gt.dsqrt(zero)) then

                  idead1 = 1

                  exit

               else if (dabs(tot(i)).gt.zero) then 

c                  write (*,'(/,a,/)') '**warning ver333** '//
c     *                   'You''ve got to ask yourself one '//
c     *                   'question: Do I feel lucky? Well, do ya, punk?'

                  idead1 = 3

               end if

            end do

            if (idead1.eq.1) then
c                                 let's blow this joint
c              write (*,'(/,a,/)') 'bad result on idead = 3, let''s '//
c    *                'blow this joint, the mass balance errors are:'
c              write (*,'(4(g14.6,2x))') (tot(i),i=1,icp)
               idead = 3

               call lpwarn (idead,'REOPT/MASS BALANCE')

               exit

c           else if (idead.eq.3) then

c              write (*,'(/,a,/)') '**warning ver333** '//
c    *                   'You''ve got to ask yourself one '//
c    *                   'question: Do I feel lucky? Well, do ya, punk?'

            end if

            idead = 0

         end if

         if (dabs(gtot-ogtot).lt.nopt(21)) then 
            quit = .true.
         else
            ogtot = gtot
         end if

c        if (iter.le.3) then 
c           d2g(iter) = gtot
c        else
c           d2g(1) = d2g(2)
c           d2g(2) = d2g(3)
c           d2g(3) = gtot
c        end if

c        if (iter.ge.3) then 
c           curve = d2g(3) + d2g(1) - 2d0*d2g(2)
c           write (*,'(g12.6,1x,i2,1x,g12.6)') curve, iter,gtot-ogtot
c           write (*,'(g12.6,1x,i2,1x,g12.6)') curve/dabs(gtot), iter,
c    *                                         (gtot-ogtot)/dabs(gtot)

c           if (dabs(curve/gtot).eq.0d0.or.
c    *          dabs((gtot-ogtot)/gtot).eq.0d0) then
c              quit = .true.
c           end if

c        end if

c        if (lopt(28)) call begtim (7)
c                                 analyze solution, get refinement points
         call yclos2 (clamda,x,is,iter,opt,idead,quit)

c        if (lopt(28)) call endtim (7,.true.,'YCLOS2 ')

         if (idead.gt.0) then 

            call lpwarn (idead,'REOPT')
            exit 

         end if 
c                                 save the id and compositions
c                                 of the refinement points, this
c                                 is necessary because resub rewrites
c                                 the zco array.
         call savpa

         if (quit) exit
c                                 save old counter 
         xphct = jphct
c                                 generate new pseudocompounds
         call resub (iter)
c                                 set the new values of is, x
         is(xphct+1:jphct) = 1
         x(xphct+1:jphct) = 0d0
         bl(xphct+1:jphct) = 0d0
         bu(xphct+1:jphct) = 1d0
c                                  save the old count
         xphct = jphct

      end do

c     if (count.gt.8000) then

c     write (*,*) ' '
c     write (*,*) 'function calls ',count
c     write (*,*) 'iterations ',rcount(1)
c     if (rcount(1).gt.0) 
c    *   write (*,*) 'function calls/iteration ',count/rcount(1)
c     write (*,*) 'good : bad ',rcount(2), rcount(3), 
c    *                          rcount(2) + rcount(3)
c     if (rcount(1).gt.0)  
c    *   write (*,*) 'iter/opt ', rcount(1)/(rcount(2) + rcount(3))
c     write (*,*) ' '

c     call prtptx

c     end if

c     rcount = 0

      
c     count = 0
c     write (*,*) count

c     write (*,*) 'end of reopt'

      end

      subroutine resub (iter)
c----------------------------------------------------------------------
c subroutine to generate new pseudocompounds around each refinenent 
c point during adaptive optimization. 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical swap, bad

      integer i, ids, lds, id, kd, iter, idif, ifail, help

      double precision gg, gsol1

      external gsol1

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ikp
      common/ cst61 /ikp(k1)

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer igood(h9), ibad(h9)
      data igood, ibad/h9*0,h9*0/
      save igood, ibad
c----------------------------------------------------------------------
c                                 reset refinement point flags
      do i = 1, jpoint
         hkp(i) = 0
      end do

c     if (.not.lopt(59)) then
c                                 if ~keep_all_rpcs, then reset counters
c        jphct = jpoint
c        zcoct = 0
c     end if
c                                 loop on previous stable phases
c                                 refine as necessay:
      lds = 0
      lsdv(1:npt) = 0

      do kd = 1, npt

         if (iter.eq.1) then 
c                                 static array pointer is
            id = jdv(kd) + istct - 1
c                                 solution model pointer is
            ids = ikp(id)
c                                 refine if a solution
            if (ids.eq.0) cycle 
c                                 get the refinement point composition
            if (id.gt.ipoint) then

               call setxyp (ids,id,bad)
c                                 save the composition for autorefine
               ststbl(id) = .true.
            else
               if (nrf(ids)) cycle
               call endpa (kd,id,ids)
            end if

            rkds = kd

         else
c                                 use pointer array lkp this uses 
c                                 negative values to index static
c                                 compounds and positive values to
c                                 point to solution models
            id = lkp(kd)
            rkds = mkp(kd)

            if (id.lt.0) then

               ids = ikp(-id)

               if (ids.eq.0) cycle 

               if (nrf(ids)) cycle

               rkds = id
c                                 endmember refinement point:
               call endpa (kd,-id,ids)

            else

               ids = id
c                                 solution refinement point:
               call getpa (ids,kd)

            end if

         end if
c                                 set local model and refinement point pointers
         rids = ids
c                                 set solution model parameters for
c                                 gsol1, don't call if the previous
c                                 refinement point was the same solution.
         if (ids.ne.lds) then
            call ingsol (ids)
            if (deriv(ids)) call ingend (ids)
         end if
c                                 save the previous solution, 
c                                 for lagged speciation with only
c                                 one solvent species this is
c                                 all that needs to be done.
         if (iter.eq.1) then
            gg = gsol1 (ids,.true.)
         else 
            gg = gsol1 (ids,.false.)
         end if
c                                 for electrolytic fluids set 
c                                 kwak0 to record the state of the
c                                 refinement point
         kwak0 = rkwak
         idif = 0

         if (nstot(ids).gt.1) then 

            call savrpc (gg,nopt(37),swap,idif)

            if (lopt(61)) call begtim (15)
c                                  normal solution
            call minfrc

            if (lopt(61)) call endtim (15,.false.,'minfrc')

         else 
c                                 don't save non-electrolytic pure fluids
            if (rkwak) cycle
c                                 save with 0-threshold
            call savkwk (gg,0d0,swap,idif)

         end if
c                                 save the location so that the 
c                                 amount can be initialized
         lsdv(kd) = idif

         lds = ids

      end do

c     write (*,*) 'end of resub'

      end

      subroutine sortin 
c-----------------------------------------------------------------------
c sort the first npt values of jdv
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, imin

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
c----------------------------------------------------------------------
      do j = 1, npt-1

         imin = jdv(j)

         do i = j+1, npt

            if (jdv(i).lt.imin) then 
               imin = jdv(i)
               jdv(i) = jdv(j)
               jdv(j) = imin
            end if
 
         end do 

      end do 

      end 

      subroutine savpa
c----------------------------------------------------------------------
c subroutine to save a copy of adaptive pseudocompound endmember fractions
c in the temporary array ycoor (also lcoor) used by resub to generate
c the new zcoor array for the subsequent iteration.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      integer i, kcoct, id, ids

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      kcoct = 0

      do i = 1, npt

         id = jdv(i)

         if (id.le.jpoint) then 
            lkp(i) = -(id + jiinc)
            cycle
         end if 

         ids = jkp(id)
         lkp(i) = ids
c                                 cycle on a compound
         if (ids.lt.0.or.id.le.jpoint) cycle

         lcoor(i) = kcoct
c                                 it's a solution:
         ycoor(kcoct+1:kcoct+nstot(ids)) = 
     *                            zco(icoz(id)+1:icoz(id)+nstot(ids))

         kcoct = kcoct + nstot(ids)

         if (lopt(58).and.(.not.refine.or.lopt(55))) then

            pa(1:nstot(ids)) = zco(icoz(id)+1:icoz(id)+nstot(ids))
c                                 only for pp comparison
            if (lorder(ids)) call makepp (ids)

            call savdyn (ids)

         end if

      end do 

      end

      subroutine getpa (ids,id)
c----------------------------------------------------------------------
c subroutine to recover independent endmember fractions from
c from the ycoor array loaded in savpa
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, ids, kcoor
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      kcoor = lcoor(id)

      pa(1:nstot(ids)) = ycoor(kcoor+1:kcoor+nstot(ids))

      call makepp (ids)

      end

      logical function solvs1 (id1,id2,ids)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds of solution
c ids, called only for final solution vales by avrger.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2, ids

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
c-----------------------------------------------------------------------
      solvs1 = .false.

      do i = 1, icp

         if (dcp(i,ids).lt.zero) cycle 

         if (dabs(cp3(i,id1)/cptot(id1) - cp3(i,id2)/cptot(id2))
     *                                    / dcp(i,ids).gt.soltol) then
            solvs1 = .true.
            exit 
         end if 
         
      end do 

      end

      logical function solvs4 (id1,id2)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds of a lagged
c aqueous solution based 

c ids, called only for final solution values by avrger.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c-----------------------------------------------------------------------
      solvs4 = .false.

      do i = 1, ns

         if (dabs(pa3(id1,i) - pa3(id2,i)).gt.soltol) then 
            solvs4 = .true.
            exit 
         end if 
         
      end do 

      end 

      subroutine avrger (abort)
c----------------------------------------------------------------------
c avrger combines discretization points into a single solution
c composition. on output

c     np  - is the number of solutions, 
c     ncpd - is the number of true compounds
c     ntot - np+ncpd

c this routine is unecessarily complicated, because it assumes
c pseudocompounds are not ordered by solution (but they are)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical check, bad, quit, notaq, abort

      integer idsol(k19),ksol(k19,k19),ids,xidsol,xksol,irep,jlist(k5),
     *        i,j,jdsol(k19,k19),jd,k,l,nkp(k19),xjdsol(k19),kk

      double precision bsol(k19,k19),cpnew(k19,k19),xx,xb(k19),msol,
     *                 bnew(k19),pnew(k19,m14),ncaq(k19,l10),ximp,sum

      logical solvs1, solvs4
      external solvs1, solvs4
c                                 -------------------------------------
c                                 global variables:
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                  x-coordinates for the final solution
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      logical abort1
      common/ cstabo /abort1

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c-----------------------------------------------------------------------
      abort = .false.
c                                first check if solution endmembers are
c                                among the stable compounds:
      do i = 1, ntot
c                                initialize ksol, this was not done 
c                                before nov 17, 2017; god only knows
c                                how it worked...
         ksol(1,i) = 0
c                                locate solution endmembers:
         if (kkp(i).lt.0) then 
            if (ikp(-kkp(i)).ne.0) then 
c                                we have an endmember
               nkp(i) = ikp(-kkp(i))
            else
               nkp(i) = kkp(i)
            end if 
         else 
            nkp(i) = kkp(i)
         end if 
      end do 
c                                check if any solutions
      do i = 1, ntot
         if (nkp(i).gt.0) goto 10
      end do 

      np = 0
      ncpd = ntot

      goto 99
c                                figure out how many solutions
c                                are present:
10    np = 0
      ncpd = 0
      soltol = nopt(8)

      do i = 1, ntot
          
         if (nkp(i).lt.0) then
c                                 the pseudocompound is a true compound
            ncpd = ncpd + 1 
            idsol(ntot) = ncpd
            bsol(ntot,ncpd) = amt(i)
            ksol(ntot,ncpd) = nkp(i)
            jdsol(ntot,ncpd) = i

         else 

            if (lopt(32).and.ksmod(nkp(i)).eq.39) then

               notaq = .false.
c                                 get lagged speciation
c                                 loaded into caq(i,1:ns+aqct)
               do k = 1, ns
                  pa(k) = pa3(i,k)
               end do 

               if (quack(jdv(i))) then 
c                                 pure solvent phase
                  msol = 0d0

                  do k = 1, ns

                     caq(i,k) = pa(k)
c                                 solvent molar weight
                     msol = msol + pa(k) * fwt(jnd(k))

                  end do 

                  do k = sn1, nat
                     caq(i,k) = 0d0
                  end do
c                                 total molality
                  caq(i,na2) = 1d0/msol
                  caq(i,na3) = msol

               else
c                                 impure solvent, get speciation
c                                 ximp, xb, sum, and msol are dummies
                  call gaqlgd (ximp,xb,sum,msol,i,bad,.true.)

               end if

            else

               notaq = .true.

            end if

            quit = .false.

            do j = 1, np
c                                 compare the compound to the np solutions 
c                                 identfied so far:        
               if (ksol(j,1).eq.nkp(i)) then

                  kk = jdsol(j,idsol(j))
c                                 if match check for a solvus
                  if (notaq) then

                     if (solvs1(i,kk,nkp(i))) cycle

                  else
c                                  special solvus test based on solvent 
c                                  speciation for lagged aq model.
                     if (solvs4(i,kk)) cycle

                     if (iopt(22).gt.2) then 
c                                  check pure and impure solvent coexist
                        if (caq(i,na1).eq.0d0.and.caq(kk,na1).ne.0d0.or.
     *                      caq(i,na1).ne.0d0.and.caq(kk,na1).eq.0d0) 
     *                                                              then 
c                                  pure solvent and impure solvent coexist
                            abort = .true.
                            return

                        end if

                     end if

                  end if 
c                                 the pseudocompound matches a solution
c                                 found earlier.
                  idsol(j) = idsol(j) + 1
                  bsol(j,idsol(j)) = amt(i)
                  jdsol(j,idsol(j)) = i
                  quit = .true.
                  exit

               end if

            end do

            if (quit) cycle 
c                                 the pseudocompound is a new solution 
c                                 phase.
            np = np + 1
            idsol(np) = 1
            ksol(np,1) = nkp(i)
            jdsol(np,1) = i
            bsol(np,1) = amt(i)

         end if

      end do
c                                 check if a solution occurs more than once
c                                 but the occurences are not sequential (this
c                                 can only occur if an endmember is immiscible 
c                                 with a general composition
      if (np.gt.2) then
 
         do i = 1, np

            check = .false.
            irep = 0

            do j = i+1, np
               if (ksol(j,1).ne.ksol(i,1)) then
                  check = .true.
               else 
                  irep = irep + 1
               end if 
            end do 

            if (check.and.irep.gt.0) then

               l = i + 1

               if (ksol(l,1).ne.ksol(i,1)) then 
c                                 not in sequence, find the next occurence
                  do j = i+2, np 
                     if (ksol(i,1).eq.ksol(j,1)) exit
                  end do 
c                                 swap phase at i+1 with the one at j
                  xidsol = idsol(l)
                  xksol = ksol(l,1)
                  do k = 1, xidsol
                     xb(k) = bsol(l,k)
                     xjdsol(k) = jdsol(l,k)
                  end do 

                  idsol(l) = idsol(j)
                  ksol(l,1) = ksol(j,1)
                  do k = 1, idsol(j)
                     bsol(l,k) = bsol(j,k)
                     jdsol(l,k) = jdsol(j,k)
                  end do 

                  idsol(j) = xidsol
                  ksol(j,1) = xksol
                  do k = 1, xidsol
                     bsol(j,k) = xb(k)
                     jdsol(j,k) = xjdsol(k)
                  end do 

               end if 
            end if 
         end do 
      end if 
c                                 if a solution is represented by
c                                 more than one pseudocompound get
c                                 the average composition
      do i = 1, np 
c                                 initialize
         bnew(i) = 0d0
         ids = ksol(i,1)

         cpnew(1:icomp,i) = 0d0

         pnew(i,1:nstot(ids)) = 0d0
c                               lagged speciation
         if (lopt(32).and.ksmod(ids).eq.39) ncaq(i,1:nat) = 0d0

         do j = 1, idsol(i)
            bnew(i) = bnew(i) + amt(jdsol(i,j))
         end do
c                                in case pure and impure solvent is going to be averaged
c                                count fraction of impure solvent
         ximp = 0d0
         sum = 0d0

         do j = 1, idsol(i)

            jd = jdsol(i,j)

            if (.not.refine.or.lopt(55)) then
c                                 savdyn with zero tolerance assures
c                                 final solution will be saved even if
c                                 lopt(58) is true.
c           if (.not.lopt(58).and.(.not.refine.or.lopt(55))) then 
c                                 load into pa and save for refinement
               pa(1:nstot(ids)) = pa3(jd,1:nstot(ids))
c                                 for pp comparison only
               if (lorder(ids)) call makepp (ids)

               call savdyn (ids)

            end if
c                                conditional for zero-mode stable phases
            if (bnew(i).gt.0d0) then 
               xx =  amt(jd)/bnew(i)
            else 
c                                this is gonna be a disaster if idsol(i) > 1
               xx = 1d0
            end if 
c                                save the new compositions
            do k = 1, icomp
               cpnew(k,i) = cpnew(k,i) + xx*cp3(k,jd)
            end do

            sum = sum + xx
            
            do k = 1, nstot(ids)
               pnew(i,k) = pnew(i,k) + xx*pa3(jd,k)
            end do

            if (lopt(32).and.ksmod(ids).eq.39) then

               if (caq(jd,na1).ne.0d0) ximp = ximp + xx
c                                lagged speciation (1:nsa), ionic strength (na1), total
c                                molality (na2), solvent mass (na3), err_log_kw (na4)
c                                pH, Delta_pH, solute molality, epsilon (nat)
               do k = 1, nat
                  ncaq(i,k) = ncaq(i,k) + xx*caq(jd,k)
               end do

            end if

         end do

         if (sum.lt.1d0-zero.or.sum.gt.1d0+zero) then 
            write (*,*) 'wugga'
         end if 

c         if (lopt(32).and.ksmod(ids).eq.39.and.ximp.gt.0d0) then
c                                 renomalize err_log_kw, pH, Delta_pH, epsilon
c            ncaq(i,na3+1) = ncaq(i,na3+1)/ximp
c            ncaq(i,na3+2) = ncaq(i,na3+2)/ximp
c            ncaq(i,na3+3) = ncaq(i,na3+3)/ximp
c            ncaq(i,nat) = ncaq(i,nat)/ximp

c         end if

      end do


      do i = 1, np
         sum = 0d0
        
            do k = 1, nstot(ksol(i,1))
               sum = sum + pnew(i,k)
            end do
c DEBUG691
            if (sum.lt.one.or.sum.gt.1d0+zero) then
               write (*,*) 'bad pa3 sum',ksol(i,1),sum
            end if

      end do
c                                 make a list of solutions as ordered 
c                                 in the input:
      if (np.gt.1) then

         call assort (jlist,ksol,np)

      else if (np.eq.1) then

         jlist(1) = 1

      end if
c                                 now reform the arrays kdv and b
      do i = 1, np

         kk = jlist(i)
         ids = ksol(kk,1)
         amt(i) = bnew(kk)
         kkp(i) = ksol(kk,1)

         cp3(1:icomp,i) = cpnew(1:icomp,kk)

         pa3(i,1:nstot(ids)) = pnew(kk,1:nstot(ids))
c                                 lagged speciation, ionic strength, tot molality
c                                 and solvent mass.
         if (lopt(32).and.ksmod(ids).eq.39) 
     *                            caq(i,1:nat) = ncaq(kk,1:nat)
c                                 if auto_refine is on:
c                                 check composition against solution model ranges
         if (lopt(11)) call sollim (ids,i)

      end do

9000  format (i2,3x,i8,3x,i8,3x,g14.6)

      do i = 1, ncpd

         k = np + i
         l = ksol(ntot,i)
         amt(k) = bsol(ntot,i)
         kkp(k) = l
c                                 for the sake of completeness load
c                                 compound composition into cp3 array
         cp3(1:icomp,k) = cp(1:icomp,-l)

      end do

      ntot = np + ncpd

99    end 

      subroutine sollim (ids,jd)
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c  ids - pointer to solution model
c  jd  - pointer to the phase in the final result
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical bad

      integer ids, ii, i, j, jd

      double precision stinc

      external stinc

      integer iam
      common/ cst4 /iam

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
c                                 load the indepedent endmeber fractions 
      pa(1:nstot(ids)) = pa3(jd,1:nstot(ids))
c                                 set stable flag
      stable(ids) = .true.
c                                 recover the prismatic composition
      call p2yx (ids,bad)

      if (bad) return
c                                 check x-ranges
      do ii = 1, pop1(ids)

         if (pwt(ii).le.zero.and.ii.lt.pop1(ids)) cycle

         do i = 1, istg(ids,ii)

            do j = 1, ispg(ids,ii,i)

               if (ksmod(ids).eq.20.and.j.eq.ns) cycle 
c                                 low limit:
               if (x(ii,i,j).lt.xlo(j,i,ii,ids)) then

                  xlo(j,i,ii,ids) = x(ii,i,j)
c                                 check if solution is at an unnatural limit
                  if (x(ii,i,j).ge.xmno(ids,ii,i,j).and.
     *                x(ii,i,j).lt.xmng(ids,ii,i,j)) then

                     if (.not.lopt(3)) then
c                                  error if composition has jumped too far
                        if (xmnh(ids,ii,i,j)-x(ii,i,j).gt.
     *                      xncg(ids,ii,i,j).and..not.refine) 
     *                                  call err993 (ids,ii,i,j,.false.)
c                                 relax limits according to subdivsion model
c                                 warn if MEEMUM
                        if (iam.eq.2) call meelim (x(ii,i,j),ids,ii,i,j)

                        if (imdg(j,i,ii,ids).eq.0) then 
c                                 cartesian
 
                           xmng(ids,ii,i,j) = xmng(ids,ii,i,j) 
     *                                             - xncg(ids,ii,i,j)

                        else 
c                                 assymmetric stretching towards xmin

                           xmng(ids,ii,i,j) =  stinc (x(ii,i,j),
     *                                    -xncg(ids,ii,i,j),ids,ii,i,j)

                        end if

                        if (xmng(ids,ii,i,j).lt.0d0) 
     *                                           xmng(ids,ii,i,j) = 0d0

                     end if

                     limit(ids) = .true.

                  end if

               end if 
c                                 high limit:
               if (x(ii,i,j).gt.xhi(j,i,ii,ids)) then

                  xhi(j,i,ii,ids) = x(ii,i,j)
c                                 check if solution is at an unnatural limit
                  if (x(ii,i,j).le.xmxo(ids,ii,i,j).and.
     *                x(ii,i,j).gt.xmxg(ids,ii,i,j)) then

                     if (.not.lopt(3)) then
c                                  error if composition has jumped too far
                        if (x(ii,i,j)-xmxh(ids,ii,i,j).gt.
     *                      xncg(ids,ii,i,j).and..not.refine) 
     *                                   call err993 (ids,ii,i,j,.true.)
c                                 relax limits according to subdivsion model
c                                 warn if MEEMUM
                        if (iam.eq.2) call meelim (x(ii,i,j),ids,ii,i,j)

                        if (imdg(j,i,ii,ids).eq.0) then 
c                                 cartesian
                           xmxg(ids,ii,i,j) = xmxg(ids,ii,i,j) 
     *                                            + xncg(ids,ii,i,j)

                        else 
c                                 assymmetric stretching 
                           xmxg(ids,ii,i,j) = stinc (x(ii,i,j),
     *                                     xncg(ids,ii,i,j),ids,ii,i,j)

                        end if 

                        if (xmxg(ids,ii,i,j).gt.1d0) 
     *                                           xmxg(ids,ii,i,j) = 1d0

                     end if

                     limit(ids) = .true.

                  end if

               end if

            end do

         end do

      end do

      end

      subroutine sorter (kdbulk,ic,jc)
c----------------------------------------------------------------------
c sorter compares assemblages to those already defined and reorders 
c the phases if the assemblage has been identified earlier
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,kdbulk,ic,jc,ids,ioct,inct

      logical reord, match, nomtch, ok 

      double precision cpt(k5,k5),pat(k5,m14),bt(k5),caqt(k5,l10)
c                                 x-coordinates for the final solution
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
c----------------------------------------------------------------------
c                                 look for a match with known assemblages
      match = .false.
 
      do i = 1, iasct

         if (np.ne.iavar(1,i).or.ncpd.ne.iavar(2,i)) cycle 

         nomtch = .false.

         do j = 1, ntot

            ok = .false.

            do k = 1, ntot 

               if (idasls(k,i).eq.kkp(j)) then 

                  ok = .true.
c                                 check that the phase occurs the same 
c                                 number of times in each assemblage:
                  inct = 0 
                  ioct = 0 

                  do l = 1, np
                     if (kkp(l).eq.kkp(j)) inct = inct + 1
                     if (idasls(l,i).eq.kkp(j)) ioct = ioct + 1
                  end do 

                  if (ioct.ne.inct) then
                     nomtch = .true.  
                     exit 
                  end if 

               end if 

            end do

            if (.not.ok) nomtch = .true. 
            if (nomtch) exit 

         end do 

         if (nomtch) cycle

         match = .true.
c                                 check if reordering is necessary
         reord = .false.

         do j = 1, ntot
            if (kkp(j).eq.idasls(j,i)) cycle
            reord = .true.
            exit
         end do

         if (reord) then 
c                                 reorder the result arrays of the
c                                 current occurence to match initial 
c                                 occurence:
            do j = 1, ntot

               do k = 1, ntot

                  ids = kkp(k)

                  if (ids.eq.idasls(j,i)) then
c                                 load temporary array
                     bt(j) = amt(k)

                     if (ids.gt.0) then 

                        do l = 1, icomp
                           cpt(l,j) = cp3(l,k)
                        end do

                        do l = 1, nstot(ids)
                           pat(j,l) = pa3(k,l)
                        end do

                        if (lopt(32).and.ksmod(ids).eq.39) then 

                           do l = 1, nat
                              caqt(j,l) = caq(k,l)
                           end do 

                        end if

                     end if 
c                                 this eliminates immiscible phases
                     kkp(k) = 0

                     exit 
 
                  end if 

               end do 

            end do
c                                 reload final arrays from temporary
            do j = 1, ntot

               amt(j) = bt(j)
               ids = idasls(j,i)
               kkp(j) = ids

               if (ids.gt.0) then 

                  do k = 1, icomp
                     cp3(k,j) = cpt(k,j)
                  end do

                  do k = 1, nstot(ids) 
                     pa3(j,k) = pat(j,k)
                  end do 

                  if (lopt(32).and.ksmod(ids).eq.39) then 

                     do l = 1, nat
                        caq(j,l) = caqt(j,l)
                     end do 

                  end if

               end if 
            end do 

         end if 

         if (ibulk.gt.k2) call error (183,0d0,k2,'SORTER')
         ibulk = ibulk + 1
         iap(ibulk) = i
         kdbulk = ibulk

         exit  

      end do 

      if (.not.match) then 
c                                 the assemblage is new:
         iasct = iasct + 1
         if (iasct.gt.k3) call error (184,0d0,k3,'SORTER')

         do i = 1, ntot
            idasls(i,iasct) = kkp(i)
         end do

         ibulk = ibulk + 1
         if (ibulk.gt.k2) call error (183,0d0,k2,'BLKMAT')
         kdbulk = ibulk 
         iap(ibulk) = iasct 

         iavar(1,iasct) = np
         iavar(2,iasct) = ncpd
         iavar(3,iasct) = np + ncpd

      end if 
                                
      if (outprt.or.iopt(34).ne.0) call outbl1 (ic,jc)
     
      end 

      subroutine outbl1 (ic,jc)
c----------------------------------------------------------------------
c output data for compositions and phases of assemblage ibulk
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ic,jc,i,j,ids
c                                 -------------------------------------
c                                 global variables
      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 x-coordinates for the final solution
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
c----------------------------------------------------------------------
c                                graphics output  
      write (n5,'(3(i8,1x))') ic,jc,iap(ibulk)
c                                phase molar amounts
      write (n5,1010) (amt(i),i=1,np+ncpd)
c                                solution phase compositions
      do i = 1, np

         ids = kkp(i)

         write (n5,1010) (pa3(i,j),j=1,nstot(ids))
c                                lagged speciation
         if (ksmod(ids).eq.39.and.lopt(32)) write (n5,1010) 
     *                                            (caq(i,j),j=1,nat)

      end do
c                                dependent potentials
      write (n5,1010) (mu(i),i=1,kbulk)

1010  format (10(g16.8,1x))

      end 

      subroutine yclos1 (x,clamda,jphct,quit)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement. this routine is only called as preparation
c for iterative refinement.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jphct, i, j, k, idsol(k5), kdv(h9), nsol, ids,
     *        mpt, iam, id, jdsol(k5,k5), ksol(k5), max

      external ffirst, degen 

      logical degen, solvus, quit, news, solvnt(k19)

      double precision x(*), slam(h9), clamda(*)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer ikp
      common/ cst61 /ikp(k1)

      integer tphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      npt = 0
      lcpt = 0
      lspt = 0
      nsol = 0
      quit = .true.
      soltol = nopt(25)

      do i = 1, jphct

         if (is(i).ne.1) then
c                                 make a list of found phases:
            id = i + jiinc

c           if (idegen.gt.0) then 
c              if (degen(id,1)) cycle
c           end if
c                                 currently endmember compositions are not 
c                                 refined (this is probably a mistake, but 
c                                 seems to work fine), so use id > ipoint
c                                 to identify a refineable point
            ids = ikp(id) 

            if (ids.ne.0) then 

               if (id.gt.ipoint.or.
     *             ksmod(ids).eq.39.and.lopt(32)) then 
c                                 a pseudocompound
                  quit = .false.
                  news = .true.

                  do j = 1, nsol

                     if (ids.ne.idsol(j)) cycle 
                     ksol(j) = ksol(j) + 1
                     jdsol(ksol(j),j) = id
                     news = .false.
                     exit 

                  end do

                  if (news) then 
c                                 new phase, add to list
                     nsol = nsol + 1
                     idsol(nsol) = ids
                     jdsol(1,nsol) = id
                     ksol(nsol) = 1

                  end if

               end if 

            end if 
c                                 new point, add to list
            npt = npt + 1
            jdv(npt) = i
            amt(npt) = x(i)
            if (id.gt.ipoint) then 
               jkp(i) = ikp(id)
            else if (x(i).gt.zero) then
c                                 save compound amts for starting guess
               lcpt = lcpt + 1
               ldv(lcpt) = i
               lamt(lcpt) = x(i)
            end if

            if (lopt(34)) then

               if (npt.eq.1) then
                  write (*,'(/,a)') 'Iteration dump at: '
                  write (*,'(a,'' = '',g12.6)') 
     *                  (vname(jv(j)), v(jv(j)), j = 1, ipot)
                  write (*,'(/,a,i2,a,i7)') 'iteration ',0,' jphct = ',
     *                  jphct
               end if 
               if (ikp(id).ne.0) then 
                  call dumper (1,i,0,ikp(id),x(i),clamda(i))
               else 
                  call dumper (1,i,0,-id,x(i),clamda(i))
               end if 
            end if

         end if 

      end do
c                                 amt < 0
      if (npt.gt.icp) call reject (is,1,solvnt)
c                                 amt < featol
      if (npt.gt.icp) call reject (is,2,solvnt)
c                                 is = 4
c     if (npt.gt.icp) call reject (is,3,solvnt)
c                                 get mus for lagged speciation
c     js(1:npt) = is(jdv(1:npt))

      call getmus (1,0,is,solvnt,.false.)

      if (.not.mus) then
         call muwarn (quit,0)
         return
      else
         xmu(1:icp) = mu(1:icp)
      end if

      do i = 1, isoct
         slam(i) = 1d99
         kdv(i) = 0
      end do 
c                                 perp 6.6.3, make a list of the least 
c                                 metastable composition of each solution.
      do 20 i = jpoint+1, jphct

c DEBUG why was this here? added ~6.7.6, removed april 21, 2017
c i think clamda(i).lt.0 allows degenerate compositions (and probably 
c therefore the 6.7.6 version may be better, on the bright side with
c it removed the solution composition gets refined (if endmember comps are
c being allowed, see ldsol code); restored again april 2017. removed sep 2018.

         if (is(i).ne.1) cycle

c        if (degen(i,1)) cycle

         id = i + jiinc 
         iam = ikp(id)

         if (lname(iam).eq.'liquid'.and.v(2).lt.nopt(20)) cycle 

         if (clamda(i).lt.slam(iam)) then
c                                the composition is more stable
c                                than the previous composition 
c                                of the solution. check if it's 
c                                one of the stable solutions
            do j = 1, nsol

                if (iam.eq.idsol(j)) then
c                                it's already stable, only accept
c                                it if its further than the solvus
c                                tolerance from any of the stable
c                                compositions.
                  do k = 1, ksol(j)
                     if (.not.solvus(jdsol(k,j)-jiinc,i,iam)) goto 20
                  end do

               end if
            end do
c                                the composition is acceptable
            slam(iam) = clamda(i)
            kdv(iam) = i

         end if

20    continue
c                                 load the metastable solutions into kdv
      mpt = 0 

      do i = 1, isoct

         if (kdv(i).eq.0) cycle
         mpt = mpt + 1
         kdv(mpt) = kdv(i)
         slam(mpt) = slam(i)
         quit = .false.

      end do 

      if (mpt.le.iopt(31)) then 
c                                 less metastable refinement points than
c                                 iopt(31)
         max = mpt

      else 
c                                 sort the metastable points to
c                                 find the most stable iopt(31) points
         max = iopt(31)

         call ffirst (slam,kdv,1,mpt,max,h9,ffirst)

      end if 
 
      do i = 1, max

         jdv(npt+i) = kdv(i)

         if (lopt(34)) then

            id = kdv(i)

            if (ikp(id+jiinc).ne.0) then 
               call dumper (1,id,0,ikp(id+jiinc),x(id),clamda(id))
            else 
               call dumper (1,id,0,-(id+jiinc),x(id),clamda(id))
            end if 

         end if

      end do

      if (quit) then 
c                                 zero mode filter and 
c                                 save amounts for final processing
         mpt = npt
         npt = 0

         do i = 1, mpt
            if (x(jdv(i)).lt.nopt(9)) cycle 
            npt = npt + 1
            jdv(npt) = jdv(i)
            amt(npt) = x(jdv(i))
         end do 

      else 

         npt = npt + max
c                                 sort the phases, why? don't know, but it's 
c                                 necessary
         call sortin
c                                 save the amounts for lp starting point
         do i = 1, npt
            lspt = lspt + 1
            lsamt(lspt) = x(jdv(i))
            lsst(lspt) = is(jdv(i))
         end do

      end if

      end

      logical function solvus (id1,id2,ids)
c-----------------------------------------------------------------------
c function to test if a solvus separates two static pseudocompounds of
c solution ids. called only by yclos1, modified in 688 to use the static
c composition matrix a, rather than cp.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id1, id2, ids

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c-----------------------------------------------------------------------
      solvus = .false.

      do i = 1, icp

         if (dcp(i,ids).eq.0d0) cycle

         if (dabs(a(i,id1)-a(i,id2))/dcp(i,ids).gt.soltol) then
            solvus = .true.
            exit
         end if

      end do

      end

      recursive subroutine ffirst (a, ind, left, right, k, n, dumsub)
c-----------------------------------------------------------------------
c find the k smallest values of array between indices left and right
c from http://en.wikipedia.org/wiki/Selection_algorithm
c-----------------------------------------------------------------------
      implicit none

      integer left, right, k, n, pivot, opivot, partit, ind(n)

      external dumsub

      double precision a(n)

      if (right.gt.left) then 

         opivot = left + (right-left)/2
         pivot = partit (a, ind, left, right, opivot, n)

         if (pivot.gt.k) then 
            call dumsub (a,ind,left,pivot-1,k,n,dumsub)
         else if (pivot.lt.k) then 
            call dumsub (a,ind,pivot+1,right,k-pivot,n,dumsub)
         end if 

      end if 

      end 

      integer function partit (a, ind, left, right, opivot, n)

      implicit none

      integer left, right, n, pivot, opivot, iold, ind(n), i

      double precision a(n), value, oldval

      value = a(opivot)
c                                 swap a(opivot) with a(right)
      iold = ind(opivot)
      a(opivot) = a(right)
      ind(opivot) = ind(right)
      a(right) = value
      ind(right) = iold

      pivot  = left

      do i = left, right-1

         if (a(i).le.value) then
 
            iold = ind(pivot)
            oldval = a(pivot)
            a(pivot) = a(i)
            ind(pivot) = ind(i)
            a(i) = oldval
            ind(i) = iold
            
            pivot = pivot + 1
         end if 
      end do 
c                                 swap a(right) with a(pivot)
      iold = ind(pivot)
      oldval = a(pivot)
      a(pivot) = a(right)
      ind(pivot) = ind(right)
      a(right) = oldval 
      ind(right) = iold
      
      partit = pivot

      end 

      subroutine subst1 (n)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x(k5), sum

      integer n, i, j, im1, ip1, nm1, ii, ip

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5) 

c                            solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)

      do i = 2, n
         sum = 0d0
         im1 = i-1
         do j = 1, im1
            sum = a(i,j)*x(j)+sum
         end do 
         ip = ipvt(i)
         x(i) = b(ip)-sum
      end do 
c                            solve ux = y for x:
      x(n) = x(n)/a(n,n)
      nm1 = n-1

      do ii = 1, nm1
         i = n-ii
         ip1 = i+1
         sum = 0d0
         do j = ip1, n
            sum = a(i,j)*x(j)+sum
         end do
         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)
      end do 

      b(n) = x(n)

      end

      subroutine yclos0 (x,is,jphct)
c----------------------------------------------------------------------
c subroutine to save optimization results for non-iterative refinement
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jphct, is(*)

      logical solvnt(1)

      double precision x(*)
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
      npt = 0

c     if (jphct.gt.jpoint) call errdbg ('jphct>jpoint, yclos0')

      do i = 1, jphct

         if (is(i).eq.1.or.x(i).lt.nopt(9)) cycle  
c                                 acceptable cases 0 active, between bounds
c                                                  2 active, upper bound 
         npt = npt + 1
         jdv(npt) = i
         amt(npt) = x(i)
c                                 this will be wrong if jphct > jpoint
         jkp(jdv(npt)) = -i - jiinc

      end do

      call getmus (1,0,is,solvnt,.false.)

      end

      double precision function ginc0 (dt,dp,id)
c-----------------------------------------------------------------------
c id indicates a static pseudocompound, use ginc for dynamic compound
c-----------------------------------------------------------------------
      implicit none

      double precision dt, dp, gphase

      external gphase

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      p = p + dp 
      t = t + dt 

      ginc0 = gphase (id)

      p = p - dp 
      t = t - dt

      end 

      subroutine yclos2 (clamda,x,is,iter,opt,idead,quit)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration > 1. quit is true for final
c iteration.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvus, degen

      external ffirst, solvus, degen

      integer i, id, jd, is(*), jmin(k19), opt, kpt, mpt, iter, tic, 
     *        idead, j, k

      double precision clamda(*), clam(k19), x(*)

      logical stabl(k19), solvnt(k19), quit, abort, test, good, bad

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save tic
      data tic/0/
c----------------------------------------------------------------------
c                                 npt is the number of points refined
c                                 from the previous cycle, opt is the
c                                 number of points in the original 
c                                 solution.
      do i = 1, k19
         jmin(i) = 0
         clam(i) = 1d99
         stabl(i) = .false.
      end do

      abort = .false.
      test = .false.

      npt = 0
      mpt = 0

      do i = 1, jphct

         if (is(i).eq.1.and.x(i).ne.0d0) then 
            write (*,*) ' is = 1, x = ',x(i),i,iter
         end if
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)
c                                 jd the solution model
         jd = jkp(i)
c                                 is = -1 violates lower bound (bad solution)
c                                 is = -2 violates upper bound (bad solution)
c                                 is = 0 not bounding
c                                 is = 1 lower bound
c                                 is = 2 upper bound
c                                 is = 3 input constraint
c                                 is = 4 temporary constraint (weak solution)
         if (is(i).ne.1) then
c                                 a stable point, add to list
            npt = npt + 1
            jdv(npt) = i
            amt(npt) = x(i)
            if (clamda(i).ne.0d0) is(i) = 4

            if (id.gt.0) stabl(id) = .true.

            if (lopt(32)) then
c                                 for lagged aq speciation
c                                 classify as solvent/solid
               if (jkp(i).lt.0) then
c                                 setting abort to true signals 
c                                 test in getmus:
                  if (quack(-jkp(i))) abort = .true.

                  if (ikp(-jkp(i)).eq.idaq) then
                     if (.not.quack(-jkp(i))) mpt = mpt + 1
                     solvnt(npt) = .true.
                     test = .true.
                  else 
                     solvnt(npt) = .false.
                  end if

               else if (jkp(i).eq.idaq) then

                  if (.not.quack(jkp(i))) mpt = mpt + 1
                  solvnt(npt) = .true.
                  test = .true.

               else

                  solvnt(npt) = .false.

               end if

            end if

            if (lopt(34)) then
c                                 dump iteration details
               if (npt.eq.1) write (*,'(/,a,i2,a,i7)') 'iteration ',
     *                      iter-1,' jphct = ',jphct
               call dumper (2,i,hkp(i),jkp(i),x(i),clamda(i))

            end if 

         else if (jkp(i).gt.0) then
c                                 a metastable solution cpd
            if (id.gt.0) then
c                                 and not an endmember
               if (clamda(i).lt.clam(id)) then
c                                 keep the least metastable point
                     jmin(id) = i
                     clam(id) = clamda(i)

               end if
            end if
         end if

      end do
c                                 amt < 0
      if (npt.gt.icp) call reject (is,1,solvnt)
c                                 amt < featol
      if (npt.gt.icp) call reject (is,2,solvnt)
c                                 is = 4
c     if (npt.gt.icp) call reject (is,3,solvnt)
c                                 at this point could signal bad result if
c                                 npt > icp, could also switch
c                                 icp for the non-degenerate 
c                                 component counter.

c                                 abort is set if pure solvent is stable, if 
c                                 it is the only solvent then reset abort
      if (abort.and.mpt.eq.0) abort = .false.
c                                 get mu's for lagged speciation
      if (abort.and.iopt(22).eq.0) then 

          quit = .true.
          idead = 103

      else 

         if (test) abort = .true.

         call getmus (iter,iter-1,is,solvnt,abort)

         if (abort) then 
c                                 undersaturated solute component
c                                 don't set idead if iopt =1, this
c                                 allows output of the result.
            if (lopt(32).and.iopt(22).ne.1.and.iopt(22).ne.99) then 
               idead = 101
            else 
               call lpwarn (101,'YCLOS2')
            end if

            if (.not.mus) then

               write (*,*) 'oink1',xmu(1:icp)
               call muwarn (quit,iter)
               mu(1:icp) = xmu(1:icp)

            end if

         else if (.not.mus) then 

            call muwarn (quit,iter)
            mu(1:icp) = xmu(1:icp)

         else

            xmu(1:icp) = mu(1:icp)

         end if

      end if 

      if (.not.quit) then
c                                 if not done iterating:
         kpt = 0

c                                 make a list of metastable refinement points
         do i = 1, opt
c                                 a refinement point can disappear if it 
c                                 converges on another, but new ones can't 
c                                 appear
            if (jmin(i).eq.0) then 
c                                 it's a compound
               cycle

            else if (stabl(hkp(jmin(i))).or.t.lt.nopt(20).and.
     *               lname(jkp(jmin(i))).eq.'liquid') then
c                                 contrary to what you might expect, the 1st condition
c                                 improves quality, because it stops the list 
c                                 from being clogged up with one phase
               cycle

            else
c                                 the refinement point is metastable, check that
c                                 it differs compositionally from the points 
c                                 already saved
               bad = .false.

               do j = 1, npt + kpt

                  if (jkp(jdv(j)).ne.jkp(jmin(i))) cycle

                  good = .false.
c                                 metastable point matches a refinement point, 
c                                 check composition
                  do k = 1, icp

                     if (dcp(k,jkp(jdv(j))).eq.0d0) cycle
c                                 accept if it differs by some tolerance
                     if (dabs( (cp2(k,jdv(j))-cp2(k,jmin(i)))
     *                     / dcp(k,jkp(jdv(j))) ).gt. 1d-4) then

                        good = .true.

                        exit

                     end if

                  end do

                  if (good) cycle

                  bad = .true.

                  exit

               end do

               if (bad) cycle

            end if 

            kpt = kpt + 1
            jmin(kpt) = jmin(i)
            clam(kpt) = clam(i)
c                                 temporarily use jdv(i>npt) for degeneracy test
            jdv(npt+kpt) = jmin(i)

         end do

         if (kpt.gt.iopt(31)) then 
c                                 sort if too many
            kpt = iopt(31)

            call ffirst (clam,jmin,1,kpt,iopt(31),k19,ffirst)

         end if

         do i = 1, kpt

           npt = npt + 1
           jdv(npt) = jmin(i)

           if (lopt(34)) call dumper (2,jdv(npt),hkp(jdv(npt)),
     *                       jkp(jdv(npt)),x(jdv(npt)),clamda(jdv(npt)))

         end do
c                                 sort the phases, this is only necessary if
c                                 metastable phases have been added
         call sortin
c                                 make a pointer to the original refinement 
c                                 point, this is used by resub
         do i = 1, npt
            mkp(i) = hkp(jdv(i))
         end do

      else
  
         mpt = npt
         npt = 0
c                                 check zero modes the amounts
         do i = 1, mpt

            if (x(jdv(i)).ge.nopt(9)) then

               npt = npt + 1
               amt(npt) = x(jdv(i))
               jdv(npt) = jdv(i)

            else if (x(jdv(i)).lt.-nopt(9).and.tic.lt.5) then 

               call warn (2,x(jdv(i)),i,'YCLOS2')
               tic = tic + 1

               if (tic.eq.5) call warn (49,x(1),2,'YCLOS2')

            end if 

         end do 

      end if 

      end

      subroutine nclos2 (clamda,x,is,iter,opt,idead,quit)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration > 1. quit is true for final
c iteration.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvus, degen

      external ffirst, solvus, degen

      integer i, id, jd, is(*), jmin(k19), opt, kpt, mpt, iter, tic, 
     *        idead, j, k, jmink(k19)

      double precision clamda(*), clam(k19), clamk(k19), x(*)

      logical stabl(k19), solvnt(k19), quit, abort, test, good, bad, 
     *        stablk(k19)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save tic
      data tic/0/
c----------------------------------------------------------------------
c                                 npt is the number of points refined
c                                 from the previous cycle, opt is the
c                                 number of points in the original 
c                                 solution.
      do i = 1, k19
         jmin(i) = 0
         jmink(i) = 0
         clam(i) = 1d99
         clamk(i) = 1d99
         stabl(i) = .false.
         stablk(i) = .false.
      end do

      abort = .false.
      test = .false.

      npt = 0
      mpt = 0

      do i = 1, jphct

         if (is(i).eq.1.and.x(i).ne.0d0) then 
            write (*,*) ' is = 1, x = ',x(i),i,iter
         end if
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)
c                                 jd the solution model
         jd = jkp(i)
c                                 is = -1 violates lower bound (bad solution)
c                                 is = -2 violates upper bound (bad solution)
c                                 is = 0 not bounding
c                                 is = 1 lower bound
c                                 is = 2 upper bound
c                                 is = 3 input constraint
c                                 is = 4 temporary constraint (weak solution)
         if (is(i).ne.1) then
c                                 a stable point, add to list
            npt = npt + 1
            jdv(npt) = i
            amt(npt) = x(i)
            if (clamda(i).ne.0d0) is(i) = 4

            if (id.gt.0) then 
               stabl(id) = .true.
               stablk(jd) = .true.
            end if

            if (lopt(32)) then
c                                 for lagged aq speciation
c                                 classify as solvent/solid
               if (jkp(i).lt.0) then
c                                 setting abort to true signals 
c                                 test in getmus:
                  if (quack(-jkp(i))) abort = .true.

                  if (ikp(-jkp(i)).eq.idaq) then
                     if (.not.quack(-jkp(i))) mpt = mpt + 1
                     solvnt(npt) = .true.
                     test = .true.
                  else 
                     solvnt(npt) = .false.
                  end if

               else if (jkp(i).eq.idaq) then

                  if (.not.quack(jkp(i))) mpt = mpt + 1
                  solvnt(npt) = .true.
                  test = .true.

               else

                  solvnt(npt) = .false.

               end if

            end if

            if (lopt(34)) then
c                                 dump iteration details
               if (npt.eq.1) write (*,'(/,a,i2,a,i7)') 'iteration ',
     *                      iter-1,' jphct = ',jphct
               call dumper (2,i,hkp(i),jkp(i),x(i),clamda(i))

            end if 

         else if (jkp(i).gt.0) then
c                                 a metastable solution cpd
            if (clamda(i).lt.clam(id)) then
c DEBUG DEBUG DEBUG               this is a cheap way of eliminating
c                                 compositionally redundant refinement 
c                                 points, the problem is that the critical 
c                                 value of clamda is problem/machine dependent. 
c              if (clamda(i).lt.1d-7) cycle 
c                                 keep the least metastable point
                  jmin(id) = i
                  clam(id) = clamda(i)

            end if

            if (clamda(i).lt.clamk(jd)) then
c                                 keep the least metastable point
               jmink(jd) = i
               clamk(jd) = clamda(i)

            end if

         end if

      end do

c                                 amt < 0
      if (npt.gt.icp) call reject (is,1,solvnt)
c                                 amt < featol
      if (npt.gt.icp) call reject (is,2,solvnt)
c                                 is = 4
c     if (npt.gt.icp) call reject (is,3,solvnt)
c                                 at this point could signal bad result if
c                                 npt > icp, could also switch
c                                 icp for the non-degenerate 
c                                 component counter.

c                                 abort is set if pure solvent is stable, if 
c                                 it is the only solvent then reset abort
      if (abort.and.mpt.eq.0) abort = .false.
c                                 get mu's for lagged speciation
      if (abort.and.iopt(22).eq.0) then 

          quit = .true.
          idead = 103

      else 

         if (test) abort = .true.

         call getmus (iter,iter-1,is,solvnt,abort)

         if (abort) then 
c                                 undersaturated solute component
c                                 don't set idead if iopt =1, this
c                                 allows output of the result.
            if (lopt(32).and.iopt(22).ne.1.and.iopt(22).ne.99) then 
               idead = 101
            else 
               call lpwarn (101,'YCLOS2')
            end if

            if (.not.mus) then

               write (*,*) 'oink2',xmu(1:icp)
               call muwarn (quit,iter)
               mu(1:icp) = xmu(1:icp)

            end if

         else if (.not.mus) then 

            call muwarn (quit,iter)
            mu(1:icp) = xmu(1:icp)

         else

            xmu(1:icp) = mu(1:icp)

         end if

      end if 

      if (.not.quit) then
c                                 if not done iterating:
         kpt = 0
c                                 look for metastable solutions
         do i = 1, isoct
c
            if (.not.stabl(i)) then
c                                 keep the least metastable 
c                                 composition
               kpt = kpt + 1
               jmin(kpt) = jmink(i)
               clam(kpt) = clam(i)
c                                 temporarily use jdv(i>npt) for degeneracy test
               jdv(npt+kpt) = jmink(i)

            end if

         end do
c                                 next check refinement points
         do i = 1, opt
c                                 a refinement point can disappear if it 
c                                 converges on another, but new ones can't 
c                                 appear
            if (jmin(i).eq.0) then 
c                                 it's a compound
               cycle

            else if (stabl(hkp(jmin(i))).or.t.lt.nopt(20).and.
     *               lname(jkp(jmin(i))).eq.'liquid') then
c                                 contrary to what you might expect, the 1st condition
c                                 improves quality, because it stops the list 
c                                 from being clogged up with one phase
               cycle

            else
c                                 the refinement point is metastable, check that
c                                 it differs compositionally from the points 
c                                 already saved
               bad = .false.

               do j = 1, npt + kpt

                  if (jkp(jdv(j)).ne.jkp(jmin(i))) cycle

                  good = .false.
c                                 metastable point matches a refinement point, 
c                                 check composition
                  do k = 1, icp

                     if (dcp(k,jkp(jdv(j))).eq.0d0) cycle
c                                 accept if it differs by some tolerance
                     if (dabs( (cp2(k,jdv(j))-cp2(k,jmin(i)))
     *                     / dcp(k,jkp(jdv(j))) ).gt. 1d-4) then

                        good = .true.

                        exit

                     end if

                  end do

                  if (good) cycle

                  bad = .true.

                  exit

               end do

               if (bad) cycle

            end if 

            kpt = kpt + 1
            jmin(kpt) = jmin(i)
            clam(kpt) = clam(i)
c                                 temporarily use jdv(i>npt) for degeneracy test
            jdv(npt+kpt) = jmin(i)

         end do

         do i = 1, kpt

           npt = npt + 1
           jdv(npt) = jmin(i)

           if (lopt(34)) call dumper (2,jdv(npt),hkp(jdv(npt)),
     *                       jkp(jdv(npt)),x(jdv(npt)),clamda(jdv(npt)))

         end do
c                                 sort the phases, this is only necessary if
c                                 metastable phases have been added
         call sortin
c                                 make a pointer to the original refinement 
c                                 point
         do i = 1, npt
            mkp(i) = hkp(jdv(i))
         end do

      else
  
         mpt = npt
         npt = 0
c                                 check zero modes the amounts
         do i = 1, mpt

            if (x(jdv(i)).ge.nopt(9)) then

               npt = npt + 1
               amt(npt) = x(jdv(i))
               jdv(npt) = jdv(i)

            else if (x(jdv(i)).lt.-nopt(9).and.tic.lt.5) then 

               call warn (2,x(jdv(i)),i,'YCLOS2')
               tic = tic + 1

               if (tic.eq.5) call warn (49,x(1),2,'YCLOS2')

            end if 

         end do 

      end if 

      end

      subroutine reject (is,icrit,solvnt)
c----------------------------------------------------------------------
c subroutine to reject bad results from lpopt by one of 3 criteria:
c  icrit = 1 -> amt(i) < 0
c  icrit = 2 -> amt(i) < featol
c  icrit = 3 -> is(i) = 4, weak solution
c cases 1 & 2 should be, but are not, flagged by the NAG version of 
c lpsol; this was not a problem with the pre-690 solver. This may be
c related to bugwandita, some forensic analysis is in order. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvnt(*)

      integer i, id, is(*), kpt, icrit, rej

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      kpt = 0
      rej = 0

      do i = 1, npt

         if (icrit.eq.1.and.amt(i).lt.0d0.or.
     *       icrit.eq.2.and.amt(i).lt.zero.or.
     *       icrit.eq.3.and.is(jdv(i)).eq.4) then 

c            write (*,'(a,i2,1x,g14.7,5(1x,i3))') 'rejecting ',
c    *                         icrit,amt(i),is(id),npt,npt-icp

             rej = rej + 1

             if (npt-rej.eq.icp) then

                do id = i + 1, npt
                   kpt = kpt + 1
                   jdv(kpt) = jdv(id)
                   amt(kpt) = amt(id)
                   solvnt(kpt) = solvnt(id)
                end do

                exit

            end if

         else

            kpt = kpt + 1
            jdv(kpt) = jdv(i)
            amt(kpt) = amt(i)
            solvnt(kpt) = solvnt(i)

         end if

      end do

      npt = kpt

      end

      subroutine rebulk (abort,stic)
c----------------------------------------------------------------------
c upon successful completion of an optimization with either static or
c dynamic pseudocompounds rebulk:
c     1) loads the generic arrays cp3, cptot, ctot3, kkp and x3
c     2) computes the amounts of saturated component phases.
c     4) checks for solvi and homogenizes miscible phases.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, id, jds, tictoc

      logical abort, stic, bad

      double precision c(k5), scp(k5)
 
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer ikp
      common/ cst61 /ikp(k1)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      character*5 cname
      common/ csta4 /cname(k5)
c                                 hcp is different from icp only if usv
      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      save tictoc
      data tictoc/0/
c----------------------------------------------------------------------
      do i = 1, npt
c                                 index for stoichiometric compounds, 
c                                 endmembers, and static compositions
         id = jdv(i) + jiinc

         if (stic) then
            jds = ikp(id)
         else
            jds = jkp(jdv(i))
         end if

         if (jdv(i).le.jpoint) then
c                                 load compositional data
            cptot(i) = ctot(id)

            do j = 1, icomp
               cp3(j,i) = cp(j,id)
            end do
c                                 set identifier flag
            if (ikp(id).ne.0) then

               jds = ikp(id)
c                                 an endmember
               kkp(i) = jds
c                                 load the solution composition in pa3
               pa3(i,1:nstot(jds)) = 0d0
c                                 figure out which endmember
               do j = 1, nstot(jds)
                  if (id.eq.jend(jds,2+j)) exit
               end do

               pa3(i,j) = 1d0

            else
c                                 a compound
               kkp(i) = -id

            end if

         else
c                                 save the solution model pointer
            kkp(i) = jds
c                                 get endmember fractions
            if (stic) then
               call setxyp (jds,id,bad)
            else
c                                 use getpa until it's decided 
c                                 whether to reset jphct after 
c                                 each iteration:
               call getpa (jds,i)

            end if
c                                 save endmember fractions
            pa3(i,1:nstot(jds)) = pa(1:nstot(jds))
c                                 get and save the composition
c                                 getscp uses the jdv pointer
c                                 only for lagged speciation
c DEBUG rkwak
            rkwak = .false.

            call getscp (scp,cptot(i),jds,jdv(i))

            cp3(1:icomp,i) = scp(1:icomp)

         end if
c                                 convert normalized amounts to molar 
c                                 amounts
         amt(i) = ctotal*amt(i)/cptot(i)

      end do 

      if (jbulk.gt.icp) then  
c                                 get the amounts of the saturated phases:
         do i = icp+1, jbulk
c                                 k is the saturated component pointer
            k = i - icp
c                                 initialize bulk
            c(k) = cblk(i)
c                                 subtract the amount of the component in the 
c                                 phases in the thermodynamic c-space
            do j = 1, npt 
               c(k) = c(k)- amt(j)*cp3(i,j)
            end do 

         end do 

         do i = jbulk, icp1, -1
c                                  cycle through the saturated phases
            npt = npt + 1
            id = idss(i-icp)
c                                  set case for solution in saturated component
c                                  space, the endmember composition is not set,
c                                  this is gonna cause problems, at least for 
c                                  meemum
            if (ikp(id).eq.0) then 
               kkp(npt) = -id
            else 
               write (*,1000) 
               stop 
            end if 
c                                  amount of the staurated phase
            amt(npt) = c(i-icp)/cp(i,id)
c                                  warn on undersaturation
            if (amt(npt).lt.nopt(9)) then 
               if (amt(npt).lt.-nopt(9).and.tictoc.lt.1) call warn (99,
     *             c(1),i,'the specified amount of saturated compon'//
     *                    'ent '//cname(i)//'is inadequate to saturat'//
     *                    'e the system at all conditions of interest')
               npt = npt - 1
               tictoc = tictoc + 1
               exit
            end if 
c                                  remove the saturated phase from 
c                                  the bulk composition.
            do j = icp+1, i - 1
               c(j) = c(j) - amt(npt)*cp(j,id)
            end do 
c                                  load the saturated phase composition 
            do j = 1, icomp
               cp3(j,npt) = cp(j,id)
            end do

         end do

      end if

      ntot = npt
c                                 test for solvi and average
c                                 homogeneous phases.
      call avrger (abort)

1000  format (/,'**error ver901** solutions not allowed in saturated ',
     *   'component composition space',/,'in adaptive optimization ',
     *   'calculations, this limitation can be removed upon request.')

      end

      subroutine getgc (lc,lg,ldc,iter) 
c----------------------------------------------------------------------
c iter is a flag indicating where the compositions are and is sort of 
c      related to the iteration count during optimization.
c jter is the true iteration count.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id, iter, ldc

      double precision lc(ldc,*), lg(*)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer hcp, idv
      common/ cst52  /hcp,idv(k7) 

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------
      do i = 1, npt

         id = jdv(i) 

         if (iter.gt.1) then

            lc(i,1:hcp) = cp2(1:hcp,id)

            lg(i) = g2(id)

         else

            lc(i,1:hcp) = a(1:hcp,id)

            lg(i) = c(id)

         end if

      end do

      end 

      subroutine getmus (iter,jter,is,solvnt,abort) 
c----------------------------------------------------------------------
c iter is a flag indicating where the compositions are and is sort of 
c      related to the iteration count during optimization.
c jter is the true iteration count.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvnt(*), bad, abort, cslut(k19), cslvt(k19), error, skip

      integer i, j, ier, ipvt(k8), iter, jter, imu(k8), kcp, lcp, mcp,
     *        inp(k8), is(*), kdv(k19)

      double precision comp(k8,k8), g, lc(k19,k19), lg(k19)

      character cname*5
      common/ csta4  /cname(k5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision x, y
      common/ cxt46 /x, y

      integer hcp, idv
      common/ cst52  /hcp,idv(k7) 

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      ier = 1
      error = .false.
c                                 load c and g into a local array to 
c                                 avoid a myriad of conditionals
      call getgc (lc,lg,k19,iter)
c                                 for lagged speciation 
      if (abort) then

         abort = .false.
c                                 isolc is the number of non-solvent components, 
c                                 the indices of which are in solc(isolc)
         do j = 1, isolc
            cslut(j) = .false.
            cslvt(j) = .false.
         end do 

         do i = 1, npt
            do j = 1, isolc

               if (lc(i,solc(j)).eq.0d0) cycle

               if (solvnt(i)) then 
                  cslvt(j) = .true.
               else 
                  cslut(j) = .true. 
               end if 

            end do
         end do

         do j = 1, isolc

            if (cslvt(j).and..not.cslut(j)) then

               bad = .false.
c                                 skip if the system is degenerate,
c                                 necessary only for oxides
               do i = 1, idegen
                  if (idg(i).eq.solc(j)) then 
                     bad = .true.
                     exit 
                  end if
               end do 

               if (bad) cycle 
c                                 a component is present only in the solvent
c                                 iteration will become unstable
               abort = .true.

               write (n13,'(i4,1x,4(g14.6,1x),a)') 1000+solc(j), 
     *                                             x, y, t, p,
     *                                'disolved_non-solvent_component'
               exit

            end if
         end do 

      end if 
c                                 look for a general solution if npt > icp
      if (npt.ge.hcp) then

         ier = 0 

         if (npt.eq.hcp) then

            do i = 1, npt

               comp(i,1:hcp) = lc(i,1:hcp)

               mu(i) = lg(i)

            end do

         else
c                                need to eliminate npt - hcp phases
c                                try filtering weak solution results
            kcp = npt - hcp
            lcp = 0
c                                first the obvious bad actors
            do i = 1, npt

               if (is(jdv(i)).eq.4.and.kcp.gt.0) then
                  kcp = kcp - 1
                  cycle
               end if

               lcp = lcp + 1

               kdv(lcp) = jdv(i)

            end do

            if (kcp.gt.0) then
c                                 next chop multiple instances of the 
c                                 same solution, this could be done a 
c                                 lot more cleverly
               mcp = 0

               do i = 1, lcp

                  skip = .false.

                  do j = 1, i - 1

                     if (jkp(kdv(j)).eq.jkp(kdv(i)).and.kcp.gt.0) then 
                        kcp = kcp - 1
                        skip = .true.
                        exit
                     end if

                  end do

                  if (skip) cycle

                  mcp = mcp + 1
                  kdv(mcp) = kdv(i)

               end do

               lcp = mcp

            end if


            if (lcp.ne.hcp) then

               ier = 2

            else

               lcp = 1

               do i = 1, npt

                  if (jdv(i).ne.kdv(lcp)) cycle

                  comp(lcp,1:hcp) = lc(i,1:hcp)

                  mu(lcp) = lg(i)

                  lcp = lcp + 1

               end do

            end if

         end if

         if (ier.eq.0) then

            call factor (comp,k8,hcp,ipvt,error)

            if (.not.error) call subst (comp,k8,ipvt,hcp,mu,error)

         end if

      end if 
c                                 ier ne 0 => look for degenerate solution
      if (idegen.gt.0.and.(ier.ne.0.or.error)) then

         ier = 0 
         kcp = 0

         do i = 1, npt
c                                 check if the phase contains a degenerate component
            bad = .false.

            do j = 1, idegen
               if (lc(i,idg(j)).lt.zero) cycle
c                                 reject
               bad = .true.
               exit
            end do

            if (bad) cycle
c                                 load the phase
            kcp = kcp + 1

            inp(kcp) = i

            comp(kcp,1:jdegen) = lc(i,jdg(1:jdegen))

            mu(kcp) = lg(i)

         end do

         if (kcp.gt.jdegen) then
c                                 over-determined, try eliminating
c                                 phases present in zero amount
            lcp = 0

            do i = 1, kcp

               if (amt(inp(i)).lt.zero) cycle

               lcp = lcp + 1

               mu(lcp) = mu(i)

               comp(lcp,1:jdegen) =  comp(i,1:jdegen)

            end do

            if (lcp.ne.jdegen) ier = 2

         else if (kcp.lt.jdegen) then

            ier = 3

         end if

         if (ier.eq.0) then 

            call factor (comp,k8,jdegen,ipvt,error)

            if (.not.error) call subst (comp,k8,ipvt,jdegen,mu,error)

            if (.not.error) then 
c                                 load the chemical potentials 
c                                 into their correct positions
               do i = jdegen, 1, -1
                  mu(jdg(i)) = mu(i)
               end do 
c                                 and NaN the missing values
c                                 in 691 for degenerate cases set
c                                 to a large positive number, this 
c                                 introduces a scale dependence
               do i = 1, idegen
                  mu(idg(i)) = nopt(7)
               end do

            end if

         end if 

      end if

      if (ier.ne.0) error = .true.

      if (.not.error) then 
c                                 output as requested:
         mus = .true. 

         if (lopt(33)) then 
c                                 output iteration bulk G and mu's
            g = 0d0

            do i = 1, hcp

               if (isnan(mu(i))) then

                  imu(i) = 0d0
                  cycle

               else 

                  g = g + cblk(i)*mu(i)
                  imu(i) = idint(mu(i))

               end if 

            end do

            if (jter.eq.0.or.lopt(34)) 
     *                    write (*,1000) (cname(i), i = 1, hcp)
            write (*,1010) jter, g, (imu(i), i = 1, hcp)
            if (lopt(34)) write (*,'(/)') 

         end if

      else

         mus = .false.

         if (lopt(33)) then 
c                                 output failure msg
            if (jter.eq.0.or.lopt(34)) 
     *                    write (*,1000) (cname(i), i = 1, hcp)
            write (*,1020) jter
            if (lopt(34)) write (*,'(/)') 

         end if
c                                 failed
         mu(1:hcp) = nopt(7)

      end if

1000  format (/,'Iteration',4x,'G(J/mol)',7x,20(5x,a))
1010  format (3x,i2,4x,f15.2,6x,20(i9,1x))
1020  format (3x,i2,4x,'chemical potential back-calculation failed.')

      end

      subroutine meemum (bad)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, idead

      logical nodata, bad

      integer itri(4),jtri(4),ijpt

      double precision wt(3), cum

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------- 
c                                 initialization
      rxn = .false.
c                                 normalize the composition vector, this 
c                                 is necessary for reasons of stupidity (lpopt0). 
      ctotal = 0d0

      do i = 1, icp
          ctotal = ctotal + cblk(i)
      end do 

      do i = 1, icp
         b(i) = cblk(i)/ctotal
      end do
c                                 set dependent variables
      call incdp0
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
      if (lopt(28)) call begtim(30)

      call lpopt0 (idead)

      if (lopt(28)) then 

         call endtim (30,.true.,'Total Opt ') 

         cum = 0d0 

         do i = 1, 29

            cum = cum + times(i)

         end do

         write (*,'(/,a,2x,g14.7,//,a)') 'sum of timed intervals ',cum,
     *                                 '----------------------------'
         write (666,'(/,a,2x,g14.7,//,a)') 'sum of timed intervals ',cum
     *                                ,'----------------------------'

      end if 

      if (idead.eq.0) then
c                                 compute derivative properties
         call getloc (itri,jtri,ijpt,wt,nodata)

         bad = .false.

      else 

         bad = .true.

      end if

      end 

      subroutine iniprp
c----------------------------------------------------------------------
c iniprp - read data files and initialization for meemum
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      logical first, err 
c-----------------------------------------------------------------------
c                                 version info
      call vrsion (6)
c                                 stage flag
      refine = .false.

      first = .true.
c                                 initialize outprt to .false. to force input1 to 
c                                 read input, subsequently outprt is set in setau2
      outprt = .false.
c                                 -------------------------------------------
c                                 open statements for units n1-n5 and n9
c                                 are in subroutine input1
      call input1 (first,err)
c                                 read thermodynamic data on unit n2:
      call input2 (first)
c                                 allow reading of auto-refine data 
      call setau1 
c                                 read data for solution phases on n9:
      call input9 (first)
c                                 load static compositions for manual autorefine
      if (refine) then
         call reload (refine)
      else
c                                 initialize the counter for future auto-refine
c                                 stage static compositions
         tpct = 0
         tcct = 0
         stpct = 1
c                                 call initlp to initialize arrays 
c                                 for optimization.
         call initlp

      end if
c                                 seismic data summary file
      if (lopt(50)) call outsei

c                                 unnecessary stuff for meemum
      call setau2

      end

      logical function degen (index,array)
c----------------------------------------------------------------------
c check compounds for degeneracy
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, index, array

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct
c----------------------------------------------------------------------

      degen = .false.
c                                 turn test off if aq_oxide_components, 
c                                 in principle this would allow disproportionation
c                                 for non-elemental components
      if (lopt(36)) return

      do i = 1, idegen

         if (array.eq.1) then 
            if (a(idg(i),index).gt.zero) then
               degen = .true.
               exit
            end if
         else if (array.eq.2) then 
            if (cp2(idg(i),index).gt.zero) then
               degen = .true.
               exit
            end if
         end if 

      end do

      end

