c routines only called by vertex/meemum

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

      integer k, idead, inc, lphct, jter, lpprob

c     parameter (liw=2*k1+3,lw=2*(k5+1)**2+7*k1+5*k5)

      double precision ax(k5),x(k1),oldt,oldp,gtot,
     *                 tol,oldx,clamda(k1+k5)

c     integer iwbig(liwbig)
c     double precision wbig(lwbig)

      logical quit, abort

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      double precision g
      common/ cst2 /g(k1)

      integer jphct,istart
      common/ cst111 /jphct,istart

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision bl,bu
      common/ cstbup /bl(k1+k5),bu(k1+k5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer tphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct

      double precision wmach
      common/ cstmch /wmach(10)

      save ax, x, clamda
c     save ax, x, clamda, w, iw
c-----------------------------------------------------------------------
      inc = istct - 1

      oldt = t
      oldp = p
      oldx = xco2
c                                logarithmic_X option
      if (lopt(37)) xco2 = 1d1**xco2
c                                t_stop/p_stop options
      if (lopt(46)) then
c                                pt_freeze
         t = nopt(12)
         p = nopt(3)
      else if (t.lt.nopt(12).and.nopt(3).eq.0d0) then
         t = nopt(12)
      else if (p.lt.nopt(3).and.nopt(12).eq.0d0) then
         p = nopt(3)
      else if (t.le.nopt(12).and.p.le.nopt(3)) then
         t = nopt(12)
         p = nopt(3)
      end if
c                                logarithmic_p option
      if (lopt(14)) p = 1d1**p

      if (lopt(61)) call begtim (1)

      call gall

      if (lopt(61)) call endtim (1,.false.,'Static GALL ')

      do k = 1, jphct
         c(k) = g(k+inc)/ctot(k+inc)
      end do
c                                 load the adaptive refinement cpd g's
      g2(1:jpoint) = c(1:jpoint)
c                                 load the bulk into the constraint array
      bl(jphct+1:jphct+icp) = b(1:icp)
      bu(jphct+1:jphct+icp) = b(1:icp)
c                                 set quack to .true. for all static cpds
      quack(1:jphct) = .true.

      lpprob = 2
      tol = wmach(4)

      if (lopt(61)) call begtim (13)

      call lpsol (jphct,hcp,a,k5,bl,bu,c,is,x,jter,gtot,ax,clamda,
     *            iwbig,liwbig,wbig,lwbig,idead,istart,tol,lpprob)
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
c                                 find discretization points for refinement
         call yclos1 (x,clamda,jphct,quit)
c                                 returns quit if nothing to refine
         if (quit) then 

c                                 final processing, .true. indicates static
            call rebulk (abort,.true.)

         else
c                                 initialize refinement point pointers
            hkp(1:ipoint) = 0 
c                                 reoptimize with refinement
            call reopt (idead,gtot)

            if (idead.eq.0) then 
c                                 final processing, .false. indicates dynamic
               call rebulk (abort,.false.)

            else if (idead.eq.-1) then
c                                 hail mary
               jphct = lphct
               idead = 0

               call yclos0 (x,is,jphct) 
               call rebulk (abort,.true.)

            end if

         end if 

      end if

      if (idead.eq.0.and.(abort.or.abort1)) then
c                                 optimization suceeded, but lagged speciation
c                                 failed.
         if (abort) then 
c                                 abort is set for bad lagged speciation 
c                                 solutions in avrger when pure and impure 
c                                 phases with same solvent composition coexist
            idead = 102

         else
c                                 gaqlagd couldn't speciate a previously 
c                                 speciated composition; this error probably 
c                                 doesn't occur, see comment in avrger on abort1
            idead = 104

         end if

         call lpwarn (idead,'LPOPT0')

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
     *                 ogtot, bl(k21+k5), bu(k21+k5), tol

      integer is(k21+k5), iw(liw)

      double precision wmach
      common/ cstmch /wmach(10)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision x
      common/ scrtch /x(k21)

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

      jphct = jpoint
c                                 global composition coordinate counter
      zcoct = 0
c                                 --------------------------------------
c                                 generate pseudo compounds for the first 
c                                 iteration from static arrays
      call resub (1,idead)
c                                 resub can set idead 103 for out-of-bounds
c                                 HKF-gfunc
      if (idead.gt.0) then
         call lpwarn (idead,'REOPT')
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

      iter = 1

      do
c                                 iter is incremented before the operations,
c                                 i.e., on the nth iteration, iter is n+1
         iter = iter + 1
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
c                                 set quit flag, the idead = 3 case 
c                                 resets quit to .false.
         if (iter.gt.iopt(20)) quit = .true.

         if (idead.gt.0) then

            if (idead.ne.3.or.idead.eq.3.and.idead1.ne.0) then 

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

                  write (*,'(/,a,/)') '**warning ver333** '//
     *                   'You''ve got to ask yourself one '//
     *                   'question: Do I feel lucky? Well, do ya, punk?'

                  idead1 = 3

               end if

            end do

            if (idead1.eq.1) then
c                                 let's blow this joint
               write (*,'(/,a,/)') 'bad result on idead = 3, let''s '//
     *                'blow this joint, the mass balance errors are:'
               write (*,'(4(g14.6,2x))') (tot(i),i=1,icp)

               call lpwarn (idead,'REOPT/MASS BALANCE')

               exit

            else if (idead1.eq.3) then
c                                 do another iteration
               if (quit)  quit = .false.

            end if

            idead = 0

         else

            idead1 = 0

         end if

         if (dabs(gtot-ogtot).lt.nopt(21).and.iter.gt.2) then 
            quit = .true.
         else
            ogtot = gtot
         end if
c                                 idead is zero coming into yclos2:
c                                 analyze solution, get refinement points
         call yclos2 (clamda,x,is,iter,opt,idead,quit)
c                                 yclos2 can set idead:
c                                 100 - pure and impure solvent coexist
c                                 101 - undersaturated solute
c                                 109 - invalid endmember
c                                 error trap 100 is definitely bad.
c                                 error trap 101 seems to be useless
         if (idead.gt.0) then 

            call lpwarn (idead,'REOPT')
            exit

         end if
c                                 save the id and compositions
c                                 of the refinement points, this
c                                 is necessary because resub rewrites
c                                 the zco array.
         call savpa (.false.)

         if (quit) exit
c                                 save old counter 
         xphct = jphct
c                                 generate new pseudocompounds
         call resub (iter,idead)
c                                 set the new values of is, x
         is(xphct+1:jphct) = 1
         x(xphct+1:jphct) = 0d0
         bl(xphct+1:jphct) = 0d0
         bu(xphct+1:jphct) = 1d0
c                                  save the old count
         xphct = jphct

      end do

      if (iter.gt.iopt(20).and.idead.eq.0) call lpwarn (108,'REOPT')

      end

      subroutine resub (iter,idead)
c----------------------------------------------------------------------
c subroutine to generate new pseudocompounds around each refinenent 
c point during adaptive optimization. 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical swap, bad, badsol

      integer i, ids, lds, id, kd, iter, idif, idead

      double precision gg, gsol1

      external gsol1, badsol

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
c                                 reject if bad endmember
            if (badsol(ids)) cycle
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
c                                 reject if bad endmember
               if (badsol(ids)) cycle

               if (nrf(ids)) cycle

               rkds = id
c                                 endmember refinement point:
               call endpa (kd,-id,ids)

            else

               ids = id
c                                 reject if bad endmember
               if (badsol(ids)) cycle
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
c                                 when initialized the dynamic list (jpoint+1:jphct)
c                                 does not include the static solution, technically
c                                 this could be recovered without calculation, but 
c                                 here it is simply recalculated and saved. additionally
c                                 for lagged speciation it is necessary to establish
c                                 whether the solvent is pure by calculation.
         if (iter.eq.1) then 

            gg = gsol1 (ids,.true.)

            if (lopt(32).and.ksmod(ids).eq.39.and.nstot(ids).eq.1.and.
     *          abort1) then
c                                 HKF g-func out-of-range error, only tested
c                                 for pure H2O solvent
               idead = 103
               return

            end if
c                                 electrolytic fluid, set kwak0 to record state
            kwak0 = rkwak

            if (nstot(ids).gt.1) then 
c                                 a normal solution or multicomponent solvent
               call savrpc (gg,nopt(37),swap,idif)

               if (lopt(61)) call begtim (15)

               call minfrc

               if (lopt(61)) call endtim (15,.false.,'minfrc')

            else if (.not.rkwak) then
c                                 a speciated electrolytic fluid, skip
c                                 pure 1-species solvent, this may cause
c                                 bad warm start behavior
               idif = jdv(kd)

               call savkwk (gg,0d0,swap,idif)

            end if

         else

            idif = jdv(kd)
c                                 iter > 1
            if (lopt(32).and.ksmod(ids).eq.39) then
c                                 for lagged speciation calculate and update
c                                 the previous solution
               gg = gsol1 (ids,.true.)
c                                 electrolytic fluid, set kwak0 to record state
               kwak0 = rkwak
c                                 save with 0-threshold, i.e., should always 
c                                 replace existing speciated points
               if (.not.rkwak) call savkwk (gg,0d0,swap,idif)

            end if

            if (nstot(ids).gt.1) then

               if (lopt(61)) call begtim (15)
c                                  normal solution
               call minfrc

               if (lopt(61)) call endtim (15,.false.,'minfrc')

            end if

         end if
c                                 save the location so that the 
c                                 amount can be initialized
         lsdv(kd) = idif

         lds = ids

      end do

      end

      logical function badsol (ids)
c-----------------------------------------------------------------------
c function to reject solutions with bad endmembers, currently only set
c by stix EoS
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, i

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c-----------------------------------------------------------------------
      do i = 1, lstot(ids)
         if (badend(jend(ids,2+i))) then 
            badsol = .true.
            return
         end if
      end do

      badsol = .false.

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

      subroutine savpa (statik)
c----------------------------------------------------------------------
c subroutine to save a copy of adaptive pseudocompound endmember fractions
c in the temporary array ycoor (also lcoor) used by resub to generate
c the new zcoor array for the subsequent iteration.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical statik, bad

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
         if (ids.lt.0) then 
            write (*,*) 'something molto rotten in denmark'
         end if 

         lcoor(i) = kcoct
c                                 it's a solution:
         if (.not.statik) then 
c                                 get the composition from zco and, optionally
c                                 save the composition in the auto-refine lst
c                                 (savdyn)
            ycoor(kcoct+1:kcoct+nstot(ids)) = 
     *                            zco(icoz(id)+1:icoz(id)+nstot(ids))

            if (lopt(58).and.(.not.refine.or.lopt(55))) then

               pa(1:nstot(ids)) = zco(icoz(id)+1:icoz(id)+nstot(ids))
c                                 only for pp comparison
               if (lorder(ids)) call makepp (ids)

               call savdyn (ids)

            end if

         else

            call setxyp (ids,id+jiinc,bad)

            ycoor(kcoct+1:kcoct+nstot(ids)) = pa(1:nstot(ids))

         end if

         kcoct = kcoct + nstot(ids)

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

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c-----------------------------------------------------------------------
      solvs4 = .false.

      do i = 1, ns

         if (dabs(pa3(id1,i) - pa3(id2,i)).gt.nopt(38)) then 
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

      logical check, quit, notaq, abort

      integer idsol(k19),ksol(k19,k19),ids,xidsol,xksol,irep,jlist(k5),
     *        i,j,jdsol(k19,k19),jd,k,l,nkp(k19),xjdsol(k19),kk

      double precision bsol(k19,k19),cpnew(k19,k19),xx,xb(k19),msol,
     *                 bnew(k19),pnew(k19,m14),ncaq(k19,l10),ximp

      logical solvs1, solvs4
      external solvs1, solvs4
c                                 -------------------------------------
c                                 global variables:
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

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c-----------------------------------------------------------------------
      abort = .false.
      abort1 = .false.
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
c                                figure out how many solutions
c                                are present:
      np = 0
      ncpd = 0
      soltol = nopt(8)

      do i = 1, ntot

         if (amt(i).lt.nopt(9)) cycle

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
c                                 pure solvent phase, if lopt(74) this is 
c                                 to be reported as a failed optimization.
                  if (lopt(74)) abort1 = .true.

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
c                                 solvent molar mass
                  caq(i,na3) = msol

               else
c                                 impure solvent, get speciation
c                                 ximp is a dummy
                  call gaqlgd (ximp,i,.true.)

                  if (rkwak.and.lopt(74)) then
c                                 it probably ONLY happened when no chemical
c                                 potentials were obtained after the initial 
c                                 optimization, so lagged speciation wasn't
c                                 possible and rkwak wasn't initialized. this
c                                 case was corrected May 21, 2025 (7.1.13+)

c                                 how/why this happens isn't clear to 
c                                 me, since the present aqlgd calculation
c                                 should be identical to one used to generate
c                                 the point? at least for pure water, for 
c                                 more complex solvents it's conceivable the
c                                 composition was generated with a different
c                                 set of chemical potentials.
                     abort1 = .true.
                     return

                  end if

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
c                                  in 691 this was iopt(22) < 2
                     if (lopt(72)) then 
c                                  check pure and impure solvent coexist
                        if (caq(i,na1).eq.0d0.and.caq(kk,na1).ne.0d0.or.
     *                      caq(i,na1).ne.0d0.and.caq(kk,na1).eq.0d0) 
     *                                                              then 
c                                  pure solvent and impure solvent coexist
c                                  signals error ver102
                            abort = .true.
                            return

                        end if

                     end if

                  end if 
c                                 the pseudocompound matches a solution
c                                 found earlier.
                  if (amt(i).gt.nopt(9)) then 
                     idsol(j) = idsol(j) + 1
                     bsol(j,idsol(j)) = amt(i)
                     jdsol(j,idsol(j)) = i
                  end if

                  quit = .true.

                  exit

               end if

            end do

            if (quit) cycle 
c                                 the pseudocompound is a new solution 
c                                 phase.
            if (amt(i).gt.nopt(9)) then 
               np = np + 1
               idsol(np) = 1
               ksol(np,1) = nkp(i)
               jdsol(np,1) = i
               bsol(np,1) = amt(i)
            end if

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

            do k = 1, nstot(ids)
               pnew(i,k) = pnew(i,k) + xx*pa3(jd,k)
            end do

            if (lopt(32).and.ksmod(ids).eq.39) then
c                                lagged speciation (1:nsa), ionic strength (na1), total
c                                molality (na2), solvent mass (na3), err_log_kw (na4)
c                                pH, Delta_pH, solute molality, epsilon (nat)
               do k = 1, nat
                  ncaq(i,k) = ncaq(i,k) + xx*caq(jd,k)
               end do

            end if

         end do

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

      end 

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
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

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
         if (iasct.gt.k3-1) call error (184,0d0,k3,'SORTER')

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

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
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
c                                for liquidus/solidus calcs output the
c                                the additional "dependent" potential
      if (icopt.eq.2) write (n5,1010) v(iv(1))

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
c                                 post processing assume dynamic
c                                 composition pointers, create those
c                                 here in case there's no further iteration.
            if (id.gt.ipoint) then

               jkp(i) = ikp(id)

            else

               jkp(i) = -id

            end if

            if (x(i).gt.zero) then
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

               call dumper (1,i,0,jkp(i),x(i),clamda(i))

            end if

         end if 

      end do
c                                 amt < 0
      if (npt.gt.icp) call reject (is,1,solvnt)
c                                 amt = 0
      if (npt.gt.icp) call reject (is,0,solvnt)
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

      logical solvus, degen, badsol

      external ffirst, solvus, degen, badsol

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
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)
c                                 jd the solution model/compound pointer
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
               if (jd.lt.0) then
c                                 setting abort to true signals 
c                                 test in getmus:
                  if (quack(-jd)) abort = .true.

                  if (ikp(-jd).eq.idaq) then

                     if (.not.quack(-jd)) mpt = mpt + 1
                     solvnt(npt) = .true.
                     test = .true.

                  else

                     solvnt(npt) = .false.

                  end if

               else if (jd.eq.idaq) then

                  if (.not.quack(jd)) mpt = mpt + 1
                  solvnt(npt) = .true.
                  test = .true.

               else

                  solvnt(npt) = .false.

               end if

            end if

            if (lopt(79)) then
c                                 check bad eos results
               if (jd.lt.0) then 
                  j = ikp(-jd)
               else
                  j = jd
               end if

               if (j.ne.0) then
c                                  badsol checks if a compound or
c                                  endmember of a solution has an invalid EoS
                  if (badsol(j)) then

                     idead = 109
                     return

                  end if

               end if

            end if

            if (lopt(34)) then
c                                 dump iteration details
               if (npt.eq.1) write (*,'(/,a,i2,a,i7)') 'iteration ',
     *                      iter-1,' jphct = ',jphct
               call dumper (2,i,id,jd,x(i),clamda(i))

            end if 

         else if (jd.gt.0) then
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
c                                 amt = 0
      if (npt.gt.icp) call reject (is,0,solvnt)
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
      if (abort.and.lopt(70)) then 

          idead = 100
          return

      else
c                                 test is only set T for aqueous fluid
         if (test.and.lopt(71)) then
            abort = .true.
         else 
            abort = .false.
         end if

         call getmus (iter,iter-1,is,solvnt,abort)
c                                 getmus sets abort T only if aqueous fluid
c                                 and lopt(71) and a solute component is 
c                                 undersaturated
         if (abort) then
c                                 report as error, no output
               idead = 101

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

      subroutine reject (is,icrit,solvnt)
c----------------------------------------------------------------------
c subroutine to reject bad results from lpopt by one of 3 criteria:
c  icrit = 1 -> amt(i) < 0
c  icrit = 0 -> amt(i) = 0
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

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      kpt = 0
      rej = 0

      do i = 1, npt

         if (icrit.eq.0.and.amt(i).eq.0d0.or.
     *       icrit.eq.1.and.amt(i).lt.0d0.or.
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

      integer i, j, id, jds, tictoc

      logical abort, stic, bad

      double precision c(k5), scp(k5)

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
c                                 hcp is different from icp only if usv
      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

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

         jds = jkp(jdv(i))

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

c           call chkpa (jds)
c                                 save endmember fractions
            pa3(i,1:nstot(jds)) = pa(1:nstot(jds))
c                                 get and save the composition
c                                 getscp uses the jdv pointer
c                                 only for lagged speciation

c DEBUG rkwak, prior to May 21, 2025 this DEBUG code set rkwak
c to .false. (meaning speciated), regardless of static/dynamic
c all satic compounds are not speciated, ergo it's not clear to 
c me why this did not cause problems. I doubt this is necessary.

            if (.not.stic) rkwak = .false.

            call getscp (scp,cptot(i),jds,jdv(i))

            cp3(1:icomp,i) = scp(1:icomp)

         end if
c                                 convert normalized amounts to molar 
c                                 amounts
         amt(i) = ctotal*amt(i)/cptot(i)

      end do 

      if (jbulk.gt.icp) then  
c                                 get the amounts of the saturated phases:
         do i = icp1, jbulk
c                                 corrected to use full component indexing
c                                 in C. JADC 9/13/2023
            c(i) = cblk(i)
c                                 subtract the amount of the component in the 
c                                 phases in the thermodynamic c-space
            do j = 1, npt 
               c(i) = c(i)- amt(j)*cp3(i,j)
            end do 

         end do 

         do i = jbulk, icp1, -1
c                                  cycle through the saturated phases
            npt = npt + 1
            id = idss(i-icp)

            if (mcfit) then
c                                  for mc_fit save the i'd of the saturated 
c                                  phases
               j = jbulk - i + 1
               skp(j) = -id
               xskp(j) = c(i)/cp(i,id)

            end if
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
            amt(npt) = c(i)/cp(i,id)
c                                  warn on undersaturation
            if (amt(npt).lt.nopt(9)) then 

               if (amt(npt).lt.-nopt(9).and.tictoc.le.isat) 
     *             call warn (99,c(1),i,
     *                    'the specified amount of saturated compon'//
     *                   'ent '//cname(i)//' is inadequate to saturat'//
     *                    'e the system at all conditions of interest')

               npt = npt - 1
               tictoc = tictoc + 1

               cycle

            end if 
c                                  remove the saturated phase from 
c                                  the temporary bulk composition
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
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
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

      logical solvnt(*), bad, abort, cslut(k19), cslvt(k19), error, 
     *        skip, readyn

      integer i, j, ier, ipvt(k8), iter, jter, imu(k8), kcp, lcp, mcp,
     *        inp(k8), is(*), kdv(k19)

      double precision comp(k8,k8), g, lc(k19,k19), lg(k19)

      external readyn

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
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

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


      subroutine chkblk (idead)
c-----------------------------------------------------------------------
c chkblk - checks that the bulk composition generated for (pseudo-)ternary 
c gridded minimization for degeneracy and bounds, if out of bounds the
c optimization will be counted as a bad result.
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer idead, k

      integer iam
      common/ cst4 /iam

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c------------------------------------------------------------------------
      idead = 0
c                                 bounds test
      do k = 1, hcp

         if (b(k).gt.0d0) then 

            cycle

         else if (dabs(b(k)).lt.zero) then

            b(k) = 0d0

         else 

            idead = 2
c                                 allow negative compositions in meemum
            if (iam.eq.1) return

         end if

      end do

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

      end

      subroutine meemum (bad)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, idead

      logical nodata, bad

      integer itri(4),jtri(4),ijpt

      double precision wt(3)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
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
c                                 check for degeneracy and out of bounds 
c                                 compositions (idead ~0, but nothing is done).
      call chkblk (idead)
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
c     if (lopt(28)) call begtim(30)

      call lpopt0 (idead)

c     if (lopt(28)) then 

c        call endtim (30,.true.,'Total Opt ') 

c        cum = 0d0 

c        do i = 1, 29

c           cum = cum + times(i)

c        end do

c        write (*,'(/,a,2x,g14.7,//,a)') 'sum of timed intervals ',cum,
c    *                                 '----------------------------'
c        write (666,'(/,a,2x,g14.7,//,a)') 'sum of timed intervals ',cum
c    *                                ,'----------------------------'

c     end if 

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
c iniprp - read data files and initialization for meemum and vertex
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, err

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c-----------------------------------------------------------------------
c                                 version info
      call vrsion (6)
c                                 iniprp resets refine if starting from 
c                                 autorefine stage. not clear if refine really
c                                 need to be initialized, probably not.
      refine = .false.
c                                 outprt must be initialized to .false. 
c                                 to force input1 to read project file name.
c                                 subsequently outprt is set in setau2
      outprt = .false.
c                                 first controls whether input1 reads the
c                                 option file, some initializations in 
c                                 input2, and whether input2 and input9
c                                 write summaries to console
      first = .true.
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

