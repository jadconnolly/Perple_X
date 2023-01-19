
      subroutine minfrc
c-----------------------------------------------------------------------
c minimize the omega function for the independent endmember fractions
c of solution ids subject to site fraction constraints

c     number of independent endmember fractions -> nstot-1 (<m19)
c     number of independent endmember fractions -> nz (<m20)
c     closure is forced in the objective function (gsol2)

c ingsol MUST be called prior to minfrc to initialize solution/p-t
c specific properties!

c endmember gibbs energies must be computed (presumably by gall)
c prior to the call to minfxc!
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical zbad, xref, swap

      integer i, nvar, iter, iwork(m22), idif, idead,
     *        istate(m21), nclin, ntot

      double precision ggrd(m19), lapz(m20,m19),gsol1, pinc,
     *                 bl(m21), bu(m21), gfinal, ppp(m19),
     *                 clamda(m21),r(m19,m19),work(m23),
     *                 yt(m4),zsite(m10,m11), sum

      external gsol2, gsol1

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision wmach
      common/ cstmch /wmach(10)
c-----------------------------------------------------------------------
      yt = pa

      nclin = nz(rids)
      ntot = nstot(rids)

      if (equimo(rids)) then
         nvar = ntot - 1
      else 
         nvar = ntot
      end if

      ppp(1:nvar) = pa(1:nvar)
c                                 initialize bounds
      if (boundd(rids)) then
c                                 the endmember fractions are bounded
         bu(1:nvar) = 1d0
         bl(1:nvar) = 0d0
         if (.not.lorder(rids)) nclin = 0

      else 
c                                 the model has site fractions
         bu(1:nvar) = 1d0
         bl(1:nvar) = -1d0

      end if 
c                                 load the local constraints 
c                                 from the global arrays
      lapz(1:nclin,1:nvar) = apz(rids,1:nclin,1:nvar)

      bl(nvar+1:nvar+nclin) = zl(rids,1:nclin)
      bu(nvar+1:nvar+nclin) = zu(rids,1:nclin)

      if (nvar.eq.ntot) then
c                                 closure for non-equimolar ordering
         nclin = nclin + 1
         bl(nvar+nclin) = 1d0
         bu(nvar+nclin) = 1d0
         lapz(nclin,1:nvar) = 1d0

      else if (nclin.eq.0) then 
c                                 closure for molecular models
         nclin = nclin + 1
         bl(nvar+nclin) = 0d0
         bu(nvar+nclin) = 1d0
         lapz(nclin,1:nvar) = 1d0

      end if

      if (deriv(rids)) then

         numric = .false.

      else

         numric = .true.

      end if
c                                 for lagged speciation calcultions
c                                 turn off lagged speciation if 
c                                 speciaition failed at the rpc 
c                                 composition
      if (lopt(32).and.kwak0) then
         swap= .true.
         lopt(32) = .false.
      else
         swap = .false.
      end if

      call nlpsol (nvar,nclin,m20,m19,lapz,bl,bu,gsol2,iter,istate,
     *            clamda,gfinal,ggrd,r,ppp,iwork,m22,work,m23,idead)

      if (swap) lopt(32) = swap

      rcount(1) = rcount(1) + iter
c                                 reject results outright if no improvement
c                                 (idead < 0) or infeasible (idead = 3)
      if (idead.lt.0.or.idead.eq.3) then 
         rcount(3) = rcount(3) + 1
         return
      end if 
c                                 reconstruct pa-array, this IS necessary.
      call ppp2pa (ppp,sum,nvar)
c                                 reject bad site populations, these may not
c                                 be useful
      if (boundd(rids)) then
         if (sum.gt.nopt(55)) then
c           write (*,*) 'oink 1',sum,rids
            rcount(3) = rcount(3) + 1
            return
         else if (sum.gt.1d0) then 
            pa(nstot(rids)) = 0d0
         end if
      end if
      
      if (zbad(pa,rids,zsite,fname(rids),.false.,fname(rids))) then
      
c        write (*,*) 'oink 3',rids
         rcount(3) = rcount(3) + 1
         return
      
      end if

      rcount(2) = rcount(2) + 1
c                                 save the final point, the point may have
c                                 already been saved by gsol2 but because
c                                 gsol2 uses a replicate threshold of nopt(37)
c                                 a near solution rpc would prevent gsol2 from 
c                                 saving the final composition. here the replicate
c                                 threshold is reduced to zero (sqrt(eps)).
      call makepp (rids)
c                                 if logical arg = T use implicit ordering
      gfinal = gsol1 (rids,.false.)
c                                 save the final QP result
      call savrpc (gfinal,0d0,swap,idif)

      if (lopt(54).and..not.swap) then

         yt(1:nstot(rids)) = pa(1:nstot(rids))
c                                 scatter in only for nstot-1 gradients
         pinc = 1d0 + nopt(48)
c                                 in case on 1st iteration set refine to 
c                                 to true and then reset to old value on 
c                                 exit
         xref = refine
         refine = .true.

         do i = 1, lstot(rids)

            pa(1:nstot(rids)) = yt(1:nstot(rids))/pinc

            pa(i) = pa(i) + (1d0 - 1d0/pinc)

            if (zbad(pa,rids,zsite,fname(rids),.false.,fname(rids))) 
     *                                                            cycle 
            call makepp (rids)
c                                 degeneracy test removed
c                                 if logical arg = T use implicit ordering
            gfinal = gsol1 (rids,.true.)
c                                 gsol1 computes rcp, rsum for ksmod(39)
c                                 but a direct call here to getscp will not.
            if (ksmod(rids).ne.39) call getscp (rcp,rsum,rids,rids)
c                                 save the scatter point
            call savrpc (gfinal,nopt(48)/2d0,swap,idif)

         end do
c                                 reset refine
         refine = xref

      end if

      end

      subroutine ppp2pa (ppp,sum,nvar)
c-----------------------------------------------------------------------
c reconstruct pa array from ppp array for minfxc/minfrc
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nvar

      double precision ppp(*), sum

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      sum = 0d0

      do i = 1, nvar

         sum = sum + ppp(i)
         pa(i) = ppp(i)

      end do

      if (nvar.lt.nstot(rids)) pa(nstot(rids)) = 1d0 - sum

      end

      subroutine gsol2 (nvar,ppp,gval,dgdp,bad)
c-----------------------------------------------------------------------
c function to evaluate gibbs energy of a solution for minfrc. can call 
c either gsol1 with order true or false, true seems to give better results
c presumably because it's using analytical gradients.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical zbad, saved, bad

      integer i, j, nvar, idif

      double precision ppp(*), gval, dgdp(*), psum, 
     *                 gsol1, g, bsum, zsite(m10,m11)

      external gsol1, zbad

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c-----------------------------------------------------------------------
      count = count + 1

      bad = .false.

      if (lopt(61)) call begtim (2)
c                                 reconstruct pa array
      call ppp2pa (ppp,psum,nvar)

      call makepp (rids)

      if (deriv(rids)) then
c                                 analytical derivatives:
         call getder (g,dgdp,rids)

         gval = g

         do j = 1, icp
c                                 degenerate sys, mu undefined:
            if (isnan(mu(j))) cycle
c                                 convert g to g'
            gval = gval - rcp(j)*mu(j)
c                                 convert dgdp to dg'dp
            do i = 1, nvar
               dgdp(i) = dgdp(i) - dcdp(j,i,rids)*mu(j)
            end do

         end do

      else
c                                 only numeric derivatives are
c                                 available, get g at the composition
         g = gsol1(rids,.false.)
c                                 level it
         call gsol5 (g,gval)

         if (lopt(32).and.rkwak) then 
            bad = .true.
         end if

      end if

      if (lopt(57).and.outrpc) then
c                                 try to eliminate bad results
         if (psum.lt.one.or.psum.gt.1d0+zero.or.bsum.lt.zero) return
         if (zbad(pa,rids,zsite,'a',.false.,'a')) return
c                                 save the composition
         call savrpc (g,nopt(37),saved,idif)

      end if

      if (lopt(61)) call endtim (2,.false.,'Dynamic G')

      end

      subroutine gsol5 (g,gval)
c-----------------------------------------------------------------------
c gsol5 is a shell called by minfrc/gsol2 that calls gsol1 to compute
c the true gibbs energy (g) and the leveled g (gval). both pa and pp
c must be set prior to calling gsol5.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j

      double precision gval, g

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c-----------------------------------------------------------------------
      gval = g

      do j = 1, icp
c                                 degenerate sys, mu undefined:
         if (isnan(mu(j))) cycle
c                                 convert g to g'
         gval = gval - rcp(j)*mu(j)

      end do

      end

      subroutine savrpc (g,tol,swap,idif)
c-----------------------------------------------------------------------
c save a dynamic composition/g for the lp solver
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical swap

      integer i, j, ntot, ltot, ttot, ipt, ist, idif

      double precision g, diff, tol, dtol

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c-----------------------------------------------------------------------
      ntot = nstot(rids)
      ltot = lstot(rids)
      ttot = tstot(rids)

      idif = 0

      if (.not.rkwak) then
c                                special case for electrolytic fluids
         call savkwk (g,tol,swap,idif)
         return

       end if
c                                o/d models use the pp array, which is 
c                                not normalized for non-equimolar o/d, do
c                                the normalization here
      if (.not.equimo(rids)) then

         diff = 0d0

         do j = 1, ltot
            diff = diff + pp(j)
         end do

         pp(1:ltot) = pp(1:ltot)/diff

      end if

      if (tol.eq.0d0) then
         dtol = zero
      else
         dtol = tol
      end if

      swap = .false.
c                                 degenerate bulk check is in earlier 
c                                 versions, probably was never done right
      do i = jpoint + 1, jphct
c                                 check if duplicate
         if (jkp(i).ne.rids.or..not.quack(i)) cycle

         if (lorder(rids)) then
            ist = icoz(i) + ntot
         else 
            ist = icoz(i)
         end if

         diff = 0d0

         do j = 1, ltot

            ipt = ist + j

            if (lorder(rids)) then
               diff = diff + dabs(pp(j) - zco(ipt))
            else
               diff = diff + dabs(pa(j) - zco(ipt))
            end if

         end do

         if (diff.gt.dtol) cycle

         if (diff.lt.zero) then
c                                 true zero difference, set swap to 
c                                 avoid scatter point replication.
            swap = .true.
c                                 if perfect replica swap lower g's
            if (diff.eq.0d0.and.g2(i).gt.g/rsum) g2(i) = g/rsum

            idif = i

            return

         else

            return

         end if

      end do
c                                 increment counters
      jphct = jphct + 1
      idif = jphct
      icoz(jphct) = zcoct
      zcoct = zcoct + ttot
c                                 lagged speciation quack flag
      quack(jphct) = rkwak
c                                 normalize and save the composition
      cp2(1:icomp,jphct) = rcp(1:icomp)/rsum
c                                 the solution model pointer
      jkp(jphct) = rids
c                                 the refinement point pointer
      hkp(jphct) = rkds
c                                 save the normalized g
      g2(jphct) = g/rsum
c                                 sum scp(1:icp)
      c2tot(jphct) = rsum
c                                 if it quacks like a duck then...
      quack(jphct) = rkwak
c                                 save the endmember fractions
      zco(icoz(jphct)+1:icoz(jphct)+ntot) = pa(1:ntot)
c                                 and normalized bulk fractions if o/d
      if (lorder(rids)) 
     *   zco(icoz(jphct)+ntot+1:icoz(jphct)+ttot) = pp(1:ltot)

      end 

      subroutine savkwk (g,tol,swap,idif)
c-----------------------------------------------------------------------
c save dynamic electrolytic fluid compositions/g for the lp solver.
c discriminate compositions on the basis of the solvent speciation
c under the assumption that for a given solvent speciation the latest
c result is always the best result.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical swap

      integer i, j, ltot, ist, idif

      double precision g, diff, tol, dtol

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c-----------------------------------------------------------------------
      ltot = lstot(rids)

      if (tol.eq.0d0) then
         dtol = zero
      else
         dtol = tol
      end if

      swap = .false.

      do i = jpoint + 1, jphct
c                                 check if duplicate
         if (quack(i)) cycle

         ist = icoz(i)
         diff = 0d0

c         do j = 1, ltot

c            diff = diff + dabs(pa(j) - zco(ist + j))

c         end do

         do j = 1, icomp

            diff = diff + dabs(cp2(j,i) - rcp(j)/rsum)

         end do

         if (diff.gt.dtol) then

            cycle

         else if (diff.le.zero) then

            swap = .true.
            idif = i
            exit

         end if

      end do
c                                 increment counters
      if (.not.swap) then
         jphct = jphct + 1
         idif = jphct
      end if

      icoz(idif) = zcoct
      zcoct = zcoct + ltot
c                                 lagged speciation quack flag
      quack(idif) = rkwak
c                                 normalize and save the composition
      cp2(1:icomp,idif) = rcp(1:icomp)/rsum
c                                 the solution model pointer
      jkp(idif) = rids
c                                 the refinement point pointer
      hkp(idif) = rkds
c                                 save the normalized g
      g2(idif) = g/rsum
c                                 renormalize the bulk to a mole of solvent
c                                 it's no longer clear to me why this is desireable.
      c2tot(idif) = rsum/rsmo
c                                 if it quacks like a duck...
      quack(idif) = rkwak
c                                 save the endmember fractions
      zco(icoz(idif)+1:icoz(idif)+ltot) = pa(1:ltot)

      end 

      subroutine numder (g,objfun,dgdp,ppp,fdnorm,bl,bu,nvar,bad)
c-----------------------------------------------------------------------
c subroutine to evaluate the gradient numerically for minfrc/minfxc
c on input sum is the total of the fractions, for bounded models this
c can be used to choose increment sign.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer i, nvar

      double precision ppp(*), dgdp(*), oldppp, 
     *            dpp, g, g1, g3, bl(*), bu(*), fdnorm

      double precision cdint, ctol, dxlim, epsrf, eta, fdint, ftol,
     *                 hcndbd
      common/ ngg021 /cdint, ctol, dxlim, epsrf, eta,
     *                fdint, ftol, hcndbd

      external objfun
c-----------------------------------------------------------------------
      fdnorm = 0d0

      do i = 1, nvar

         oldppp = ppp(i)

         if (cntrl) then
c                                 2nd order
            if (fdincs) then
               dpp = hctl(i)
            else
               dpp = cdint
            end if

         else
c                                1st order
            if (fdincs) then
               dpp = hfwd(i)
            else 
               dpp = fdint
            end if

         end if
c                                 rel/abs scaling
         dpp = dpp * (1d0 + dabs(ppp(i)))
c                                 first increment, doubled for central
         if (cntrl) dpp = 2d0*dpp
c                                 try to avoid invalid values (z<=0)
         if (ppp(i).gt.bu(i)-dpp) then 
c
            dpp = -dpp
c
         else if (ppp(i).gt.bl(i)+2d0*dpp) then
c                                 choose direction away from closest bound
           if (bl(i) + bu(i) - 2d0*ppp(i).lt.0d0) dpp = -dpp
c
         end if
c                                 apply the increment
         ppp(i) = oldppp + dpp

         if (dabs(dpp).gt.fdnorm) fdnorm = dabs(dpp)

         if (cntrl) then
c                                 g at the double increment
            call objfun (nvar,ppp,g3,dgdp,bad)
c                                 single increment
            ppp(i) = oldppp + dpp/2d0
c                                 g at the single increment
            call objfun (nvar,ppp,g1,dgdp,bad)

            dgdp(i) = (4d0*g1- 3d0*g-g3)/dpp

         else
c                                 g at the single increment
            call objfun (nvar,ppp,g1,dgdp,bad)

            dgdp(i) = (g1 - g)/dpp

         end if
c                                 reset ppp
         ppp(i) = oldppp

      end do

      end

      subroutine gsol4 (nvar,ppp,gval,dgdp)
c-----------------------------------------------------------------------
c gsol4 - a shell to call gsol1 from minfxc, ingsol must be called
c         prior to minfxc to initialize solution specific paramters. only
c         called for implicit o/d models. 

c         returns the p0 normalized g for non-equimolar o/d.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error 

      integer ids, i, nvar

      double precision ppp(*), gval, dgdp(*), d2s(j3,j3), gord, ddq(j3)

      double precision zz, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      external gord
c-----------------------------------------------------------------------
      ids = rids
c                                   ppp(1:nord) contains the 
c                                   proportions of the ordered species
c                                   pa(lstot+1:nstot).
c                                   -----------------------------------
c                                   set the remaining proportions
      call ppp2p0 (ppp,ids)

      if (.not.maxs) then

         if (equimo(ids)) then
c                                   analytical derivatives, equimolar o/d
            call gderiv (ids,gval,dgdp,.true.,error)

         else 
c                                   analytical derivatives, non-equimolar
           do i = 1, nvar
              ddq(i) = ppp(i)-p0a(lstot(ids)+i)
           end do

           call gpderi (ids,ddq,gval,dgdp,.true.,error)

        end if

      else
c                                   negentropy minimization:
c                                   will only be called for analytical
c                                   dnu = 0 case.
         call sderiv (ids,gval,dgdp,d2s)

         if (.not.equimo(ids)) call errdbg ('piggy wee, piggy waa')

      end if

      end

      subroutine ppp2p0 (ppp,ids)
c-----------------------------------------------------------------------
c set pa from p0a given current proportions of the ordered species in
c ppp, used by minfxc
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, jd, k

      double precision ppp(*), norm

      logical pin
      common/ cyt2 /pin(j3)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      pa(1:nstot(ids)) = p0a(1:nstot(ids))
c                                   update pa for the change in the 
c                                   proportions of the ordered species
      do k = 1, nord(ids)

         if (.not.pin(k)) cycle

         jd = lstot(ids) + k

         call dpinc (ppp(k)-p0a(jd),k,ids,jd)

      end do

      if (equimo(ids)) return
c                                 non-equimolar normalization
      norm = 1d0

      do k = 1, nord(ids)
         norm = norm +  dnu(k,ids) * (ppp(k)-p0a(lstot(ids)+k))
      end do

      pa(1:nstot(ids)) = pa(1:nstot(ids)) / norm

      end

      subroutine p2yx (id,bad)
c-----------------------------------------------------------------------
c converts the independent endmember fractions to 0-1 bounded barycentric 
c coordinates:

c     number of bounding vertices -> mstot (<m4)
c     number of independent fractions -> nstot (<m14)
c     number of linear constraints -> the number of independent
c        site fractions + closure (<m20)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad, site, comp, clos, inv, zbad

      integer liw, lw, mvar, mcon, nvar, i, jter, iwarn,
     *        iwarn1, iwarn2, lpprob

      double precision scp(k5), tol

      parameter (mvar=m4, mcon=m20, liw=2*mvar+3, 
     *           lw=2*(mcon+1)**2+7*mvar+5*mcon)

      integer ncon, id, is(mvar+mcon), iw(liw), idead, istart

      double precision ax(mcon), clamda(mvar+mcon), wrk(lw), c(mvar),
     *                 a(mcon,mvar), bl(mvar+mcon), bu(mvar+mcon), 
     *                 gopt, sum, b(mcon)

      double precision wmach
      common/ cstmch /wmach(10)

      double precision ayz
      common/ csty2z /ayz(h9,m20,m4)

      double precision ayc
      common/ csty2c /ayc(h9,k5,m4)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      external zbad

      save iwarn, iwarn1, iwarn2

      data iwarn, iwarn1, iwarn2/3*0/
c-----------------------------------------------------------------------
      bad = .false.
      inv = .false.

      tol = 1d2*zero
c                                 prismatic, need to invert to vertex fractions
      if (lstot(id).lt.mstot(id)) inv = .true.
c                                 choose constraints:
      if (lorder(id)) then
c                                 decompose to stoichiometric equivaluents
         call makepp (id)

         if (inv) then
c                                 prism
            site = .true.
            comp = .false.
c                                 explicit closure definitely helps
            clos = .true.

            if (.not.equimo(id)) 
     *         call errdbg ('unanticipated prism/non-eq molar/py2x')

c                                 get the disordered p's
            call minfxc (gopt,id,.true.)

         else
c                                get sum (needed for non-eq molar case):
            sum = 0d0

            do i = 1, lstot(id)
c DEBUG691
               if (pp(i).lt.-1d-2) then 
                  write (*,*) 'wtf, p2yx 2',fname(id),' pp ',
     *                        pp(1:lstot(id))
                  bad = .true.
                  return
               end if

               if (pp(i).lt.0d0) pp(i) = 0d0

               sum = sum + pp(i)

            end do

            x(1,1,1:lstot(id)) = pp(1:lstot(id))/sum

            if (pop1(id).gt.1) 
     *         call errdbg ('houston we have a problem, p2yx 1')

         end if

      else

         if (inv) then
c                                 reciprocal and/or relict
c                                 equipartition
            comp = .true.
            clos = .false.
            site = .false.

         else

            x(1,1,1:lstot(id)) = pa(1:lstot(id))

            if (pop1(id).gt.1) 
     *         call errdbg ('houston we have a problem, p2yx 1')

         end if

      end if

      if (.not.inv) return

      nvar = mstot(id)
      ncon = 0
c                                 dummy objective function coefficients
c                                 (only 1 feasible point?)
      c(1:nvar) = 1d0
      bl(1:nvar) = 0d0
      bu(1:nvar) = 1d0

      if (site) then 
c                                 get the site fraction constraints
         call p2zind (pa,b,ncon,id)
c                                 load the fractions
         bl(nvar+1:nvar+ncon) = b(1:ncon)
         bu(nvar+1:nvar+ncon) = b(1:ncon)
c                                 load the ayz constraint matrix
         a(1:ncon,1:nvar) = ayz(id,1:ncon,1:nvar)

      end if

      if (comp) then 
c                                 load the ayc constraint matrix
         a(ncon+1:ncon+icp,1:nvar) = ayc(id,1:icp,1:nvar)
c                                 get the bulk 
         call getscp (scp,sum,id,1)
c
         bl(nvar+ncon+1:nvar+ncon+icp) = scp(1:icp)
         bu(nvar+ncon+1:nvar+ncon+icp) = scp(1:icp)
         ncon = ncon + icp

      end if

      if (clos) then 
c                                 add the closure constraint
         ncon = ncon + 1
         a(ncon,1:nvar) = 1d0
         bl(nvar+ncon) = 1d0
         bu(nvar+ncon) = 1d0

      end if
c                                 cold start
      istart = 0
c                                 feasible point
      lpprob = 1

c     if (lopt(28)) call begtim (9)

      call lpsol (nvar,ncon,a,mcon,bl,bu,c,is,y,jter,gopt,ax,
     *            clamda,iw,liw,wrk,lw,idead,istart,tol,lpprob)

c     if (lopt(28)) call endtim (9,.true.,'p2y inversion')

      if (idead.gt.0) then
c                                 really bad inversion result
         if (iwarn.lt.11) then

            write (*,1010) fname(id),idead

            call prtptx

            if (iwarn.eq.10) call warn (49,0d0,202,'P2YX')

            iwarn = iwarn + 1

         end if

         badinv(id,1) = badinv(id,1) + 1

         bad = .true.

         return

      end if
c                                 the inversion is generally weak, take any answer
c                                 within 10% of closure or positivity
      sum = 0d0

      do i = 1, mstot(id)

         sum = sum + y(i)

      end do

      if (sum.gt.1.1.or.sum.lt.0.9) then
c                                 closure violation
         if (iwarn1.lt.11) then

            write (*,1000) fname(id),(sum-1d0)*1d2

            call prtptx

            if (iwarn1.eq.10) call warn (49,0d0,201,'P2YX')
            
            iwarn1 = iwarn1 + 1

         end if

         bad = .true.

         badinv(id,1) = badinv(id,1) + 1

         return

      end if

      sum = 0d0

      do i = 1, mstot(id)

         if (y(i).lt.0d0) then
c                                 could do another inversion without
c                                 positivity constraint to see if the
c                                 answer really is outside the prism.
            if (y(i).lt.-0.05) bad = .true.

            if (iwarn2.lt.11.and.y(i).lt.-tol) then

                write (*,1020) i,y(i),fname(id)

                if (bad) then
                   write (*,1040)
                else
                   write (*,1030) i
                end if

                call prtptx

                if (iwarn2.eq.10) call warn (49,0d0,203,'P2YX')

                iwarn2 = iwarn2 + 1

            end if

            if (bad) then

               badinv(id,1) = badinv(id,1) + 1 

               return

            end if

            y(i) = 0d0

         else 

            sum = sum + y(i)

         end if

      end do
c                                 renormalize
      y(1:mstot(id)) = y(1:mstot(id))/sum

      badinv(id,2) = badinv(id,2) + 1
c                                 convert the y's to x's
      call sety2x (id)

1000  format (/,'**warning ver201** p2y inversion for ',a,' violates ',
     *       'closure by ',f5.1,'%, the result',/,'will not be used t',
     *       'o compute compositional ranges, large violations may ind',
     *       'icate that',/,'the compositional polyhedron for the mode',
     *       'l does not span all possible model compositions.',/)
1010  format (/,'**warning ver202** p2y inversion for ',a,' failed, ',
     *       'idead = ',i2,', the result',/,'will not be used t',
     *       'o compute compositional ranges.',/)
1020  format (/,'**warning ver203** negative vertex fraction y(',i2,
     *       ') = ',g8.1,' for ',a,'.',/,'Large negative values may ',
     *       'indicate that the compositional polyhedron for the model',
     *     /,'does not span all possible model compositions.',/)
1030  format ('y(',i2,') will be zeroed for computing compositional ',
     *       'ranges.',/)
1040  format ('The composition will not be used to compute',
     *       ' compositional ranges.',/)

      end

      subroutine minfxc (gfinal,ids,mxs)
c-----------------------------------------------------------------------
c optimize solution gibbs energy or configurational entropy at constant 
c composition subject to site fraction constraints.

c returns the p0 normalized g for non-equimolar o/d.

c     number of independent endmember fractions -> <j3
c     number of constraints -> < j3*j5 * 2
c     requires that pp has been loaded in cxt7

c ingsol MUST be called prior to minfxc to initialize solution/p-t
c specific properties!

c this version uses only the ordered species proportions as variables.
c the original version used numeric derivatives with all endmember proportions
c as variables, it persisted as xmnfxc until 16/12/20.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical mxs

      integer ids, i, j, k, nvar, iter, iwork(m22), itic,
     *        istate(m21), nclin, lord, idead

      double precision ggrd(m19), gordp0, g0, 
     *                 bl(m21), bu(m21), gfinal, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),
     *                 lapz(m20,m19),xp(m14)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      logical pin
      common/ cyt2 /pin(j3)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision wmach
      common/ cstmch /wmach(10)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      external gsol4, gordp0
c-----------------------------------------------------------------------
c                                 compute the disordered g for bailouts
      g0 = gordp0(ids)
      nvar = nord(ids)
      maxs = mxs

      if (equimo(ids)) then 
c                                 initialize limit expressions from p0
         call p0limt (ids)
c                                 set initial p values and count the
c                                 the number of non-frustrated od
c                                 variables.
         call pinc0 (ids,lord)

         if (icase(ids).eq.0) then 
c                                 o/d reactions are independent and
c                                 pin settings from pinc0 are valid
c                                 regardless if whether p0 is fully 
c                                 disorderd
            if (lord.eq.0) then 
               gfinal = g0
               return
            end if

         else if (maxs) then
c                                 if maxs then p0 is likely partially 
c                                 ordered, the pin settings from pinc0
c                                 can't be relied upon, a routine could 
c                                 be added to check, but given that the 
c                                 maxs inversion is mostly likely to be
c                                 called for a general composition, the
c                                 lazy solution here is to keep everything
c                                 in:
            pin = .true.
            lord = nvar

         else if (icase(ids).eq.1) then 
c                                 p0 for ~maxs will always correspond to
c                                 the disordered limit, in this case 
c                                 pin for uncorrelated o/d reactions can 
c                                 be relied up. currently the only partly
c                                 correlated case (Omph(GHP)) can also be
c                                 relied upon, but this may change. a special
c                                 test for this and fully correlated cases 
c                                 could be made, but since ~maxs calls are
c                                 only for backup from specis the lazy solution
c                                 is adopted here for the fully correlated case.
            pin = .true.
            lord = nvar

         end if
c                                 variable bounds and local (ppp) variable
c                                 initialization
         do k = 1, nord(ids)

            if (pin(k)) then
               bu(k) = 1d0
               bl(k) = -1d0
            else
               bu(k) = pa(lstot(ids)+k)
               bl(k) = pa(lstot(ids)+k)
            end if

         end do
c                                constraints
         nclin = 0

         do k = 1, nord(ids)
c                                for each constraint
            do i = 1, ln(k,ids)

               nclin = nclin + 1
c                                bounds
               bu(nvar+nclin) = -tsum(i,k)
               bl(nvar+nclin) = -tsum(i,k) - l0c(2,i,k,ids)
c                                coefficients
               lapz(nclin,1:nvar) = 0d0

               do j = 1, jt(i,k,ids)

                  lapz(nclin,jid(j,i,k,ids)-lstot(ids)) = jc(j,i,k,ids)

               end do

               lapz(nclin,k) = -1d0

            end do

         end do
c                                 initialize ppp
         ppp(1:nvar) = pa(lstot(ids)+1:lstot(ids)+nvar)

      else
c                                 not equimolar, ok for HP melt models
c                                 as the endmember fractions are generally
c                                 related to a site fraction
         nclin = 0

         call qlim (bl,bu,lord,ids) 

         if (lord.eq.0) then 
            gfinal = g0
            return
         end if
c                                 initialize ppp
         do i = 1, nvar
            ppp(i) = (bl(i)+bu(i))/2d0
         end do
c                                 need to extract sderivs from gpderi
         if (maxs) call errdbg ('oink di oink oink!!')

      end if
c                                 solution model index
      rids = ids

      itic = 0

      xp(1:nvar) = ppp(1:nvar)

      numric = .false.

      call nlpsol (nvar,nclin,m20,m19,lapz,bl,bu,gsol4,iter,istate,
     *            clamda,gfinal,ggrd,r,ppp,iwork,m22,work,m23,idead)

      if (.not.maxs) then 
         if (idead.lt.0.or.idead.eq.3) then 
            gfinal = g0
            pa = p0a
            return
         end if
      end if
c                                 if ok:
c                                 set pa to correspond to the final 
c                                 values in ppp.
      call ppp2p0 (ppp,ids)

      end

      subroutine chfd (n,fdnorm,fx,objfun,bl,bu,grad,x,bad)
c----------------------------------------------------------------------
c chfd  computes difference intervals for gradients of f(x). intervals 
c are computed using a procedure that usually requires about two 
c function evaluations if the function is well scaled. fdest and cdest 
c are the 1st and 2nd order forward differences, sdest is the 2nd order
c centered difference.

c because npcore expects 1st order numerics, grad is set to fdest here.

c assumes the function at x (objf) has been computed on entry.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical done, first, lbad1, lbad2, bad

      integer n, info, iter, itmax, j

      double precision fdnorm, bl(n), bu(n),  grad(n), x(n), cdest,
     *                 epsa, errbnd, errmax, errmin, f1, f2, fdest, fx,
     *                 h, hcd, hfd, hmax, hmin, hopt, hphi, sdest, 
     *                 sumeps, sumsd, xj

      external objfun

      double precision cdint, ctol, dxlim, epsrf, eta, fdint, ftol,
     *                 hcndbd
      common/ ngg021 /cdint, ctol, dxlim, epsrf, eta,
     *                fdint, ftol, hcndbd

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9
c----------------------------------------------------------------------
      itmax = 3

      fdnorm = 0d0

      epsa = epsrf*(1d0+abs(fx))

      do j = 1, n

         xj = x(j)

         sumsd = 0d0
         sumeps = 0d0
         hfd = 0d0
         hcd = 0d0
         hmax = 0d0
         hmin = 1d0/epspt3
         errmax = 0d0
         errmin = 0d0
         hopt = 2d0*(1d0+abs(xj))*sqrt(epsrf)
c                                 set step direction away from the nearest
c                                 bound (without regard to other constraints,
c                                 numder may do this better),
         if (bu(j) + bl(j)- 2d0*xj.lt.0d0) then
            h = -1d1*hopt
         else
            h = 1d1*hopt
         end if

         iter = 0
         cdest = 0d0
         sdest = 0d0
         first = .true.

         do

            x(j) = xj + h
            call objfun (n,x,f1,grad,bad)
            lbad1 = bad

            x(j) = xj + h + h
            call objfun (n,x,f2,grad,bad)
            lbad2 = bad

            call chcore (done,first,epsa,epsrf,fx,info,iter,itmax,cdest,
     *                   fdest,sdest,errbnd,f1,f2,h,hopt,hphi)

            if (done) exit

         end do

         if (ksmod(rids).eq.39.and.lopt(32)) then
            if (lbad1.or.lbad2) then
               bad = .true.
               return
            end if
         end if

         grad(j) = cdest

         sumsd = sumsd + abs(sdest)
         sumeps = sumeps + epsa

         if (hopt.gt.hmax) then
            hmax = hopt
            errmax = errbnd
         end if

         if (hopt.lt.hmin) then
            hmin = hopt
            errmin = errbnd
         end if

         if (info.eq.0) hcd = max(hcd,hphi)

         if (hmin.gt.hmax) then
            hmin = hmax
            errmin = errmax
         end if

         if (4d0*sumeps.lt.hmin*hmin*sumsd) then
            hfd = hmin
            errmax = errmin
         else if (4d0*sumeps.gt.hmax*hmax*sumsd) then
            hfd = hmax
         else
            hfd = 2d0*sqrt(sumeps/sumsd)
            errmax = 2d0*sqrt(sumeps*sumsd)
         end if

         if (hcd.eq.0d0) hcd = 1d1*hfd

         fdnorm = max(fdnorm,hfd)
         hfwd(j) = hfd/(1d0+abs(xj))
         hctl(j) = hcd/(1d0+abs(xj))

         x(j) = xj

      end do
c                                 signal individual increments available:
      fdincs = .true.

      end

      subroutine chcore (done,first,epsa,epsr,fx,inform,iter,itmax,
     *                   cdest,fdest,sdest,errbnd,f1,f2,h,hopt,hphi)
c----------------------------------------------------------------------
c     chcore  implements algorithm  fd, the method described in
c     gill, p.e., murray, w., saunders, m.a., and wright, m. h.,
c     computing forward-difference intervals for numerical optimization,
c     siam journal on scientific and statistical computing, vol. 4,
c     pp. 310-321, june 1983.

c     the procedure is based on finding an interval (hphi) that
c     produces an acceptable estimate of the second derivative, and
c     then using that estimate to compute an interval that should
c     produce a reasonable forward-difference approximation.

c     one-sided difference estimates are used to ensure feasibility with
c     respect to an upper or lower bound on x. if x is close to an upper
c     bound, the trial intervals will be negative. the final interval is
c     always positive.

c     chcore has been designed to use a reverse communication
c     control structure, i.e., all evaluations of the function occur
c     outside this routine. the calling routine repeatedly calls  chcore
c     after computing the indicated function values.

c     bndlo, bndup, and rho control the logic of the routine.
c     bndlo and bndup are the lower and upper bounds that define an
c     acceptable value of the bound on the relative condition error in
c     the second derivative estimate.

c     the scalar rho is the factor by which the interval is multiplied
c     or divided, and also the multiple of the well-scaled interval
c     that is used as the initial trial interval.
c----------------------------------------------------------------------
      implicit none

      double precision bndlo, bndup

      parameter (bndlo=1.0d-3,bndup=1.0d-1)

      logical done, first, ce1big, ce2big, overfl, te2big

      integer inform, iter, itmax

      double precision cdest, epsa, epsr, errbnd, f1, f2, fdest, fx, h,
     *                 hopt, hphi, sdest, afdmin, cdsave, err1, err2, 
     *                 fdcerr, fdest2, fdsave, hsave, oldcd, oldh, oldsd
     *               , rho, sdcerr, sdsave, sdiv

      external sdiv

      save              cdsave, fdsave, hsave, oldh, rho, sdsave,
     *                  ce1big, ce2big, te2big
c----------------------------------------------------------------------
      iter = iter + 1
c                                 compute forward, backward, central and second-order
c                                 difference estimates.
      fdest = sdiv (f1-fx,h,overfl)
      fdest2 = sdiv (f2-fx, 2d0*h,overfl)

      oldcd = cdest
      cdest = sdiv (4d0*f1- 3d0*fx-f2, 2d0*h,overfl)

      oldsd = sdest
      sdest = sdiv (fx- 2d0*f1+f2, h*h, overfl)
c                                 compute  fdcerr  and  sdcerr,  bounds on relative condition
c                                 errors in first and second derivative estimates.
      afdmin = min(abs(fdest),abs(fdest2))
      fdcerr = sdiv (epsa, abs(h)/2d0 *afdmin, overfl)
      sdcerr = sdiv (epsa, abs(sdest)/4d0 *h*h, overfl)
c                                 select the correct case.
      if (first) then
c                                 first time through.
c                                 check that sdcerr is in the acceptable range.
         first = .false.
         done = sdcerr .ge. bndlo .and. sdcerr .le. bndup
         te2big = sdcerr.lt.bndlo
         ce2big = sdcerr.gt.bndup
         ce1big = fdcerr.gt.bndup

         if (.not. ce1big) then

            hsave = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest

         end if

         rho = epsr**(-0.16d0)/4d0

         if (te2big) then
c                                 truncation error may be too big 
c                                 (sdcerr is too small). decrease trial interval.
            rho = 1d1*rho
            oldh = h
            h = h/rho

         else if (ce2big) then
c                                 sdcerr is too large. increase trial interval.
            oldh = h
            h = h*rho

         end if

      else if (ce2big) then
c                                 during the last iteration, the trial interval was
c                                 increased in order to decrease sdcerr.
         if (ce1big .and. fdcerr.le.bndup) then

            ce1big = .false.
            hsave = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest

         end if
c                               if sdcerr is small, accept h. otherwise,
c                               increase h
         done = sdcerr .le. bndup

         if (.not. done) then
            oldh = h
            h = h*rho
         end if

      else if (te2big) then
c                                 in last iteration, interval was decreased in order
c                                 to reduce truncation error.
         done = sdcerr.gt.bndup

         if (done) then
c                                 sdcerr jumped from too small to too
c                           l     large. accept the previous value of h.
            h = oldh
            sdest = oldsd
            cdest = oldcd

         else
c                                 test whether fdcerr is sufficiently small.
            if (fdcerr.le.bndup) then

               ce1big = .false.
               hsave = h
               fdsave = fdest
               cdsave = cdest
               sdsave = sdest

            end if
c                                 check whether sdcerr is in range.
            done = sdcerr .ge. bndlo

            if (.not. done) then
c                                 sdcerr is still too small, decrease h again.
               oldh = h
               h = h/rho

            end if

         end if

      end if
c                                 either finished or have a new estimate of h.
      if (done) then
c                                 good second-derivative estimate found.
c                                 compute optimal interval.
         hphi = abs(h)
         hopt = 2d0*sqrt(epsa)/sqrt(abs(sdest))
c                                 err1 is the error bound on the forward estimate
c                                 with the final value of h. err2 is the difference of fdest
c                                 and the central-difference estimate with hphi.
         err1 = hopt*abs(sdest)
         err2 = abs(fdest-cdest)
         errbnd = max(err1,err2)
c                                 inform = 4 if the forward- and central-difference
c                                 estimates are not close.
         inform = 0
         if (errbnd.gt.abs(fdest)/2d0) inform = 4

      else

         done = iter .ge. itmax

         if (done) then

            if (ce1big) then
c                                 fdcerr was never small.  
c                                 probably a constant function.
               inform = 1
               hphi = hopt
               fdest = 0d0
               cdest = 0d0
               sdest = 0d0
               errbnd = 0d0

            else if (ce2big) then
c                                  fdcerr was small,  but sdcerr was 
c                                  never small probably a linear or odd function.
               inform = 2
               hphi = abs(hsave)
               hopt = hphi
               fdest = fdsave
               cdest = cdsave
               sdest = 0d0
               errbnd = 2d0*epsa/hopt

            else
c                                  the remaining case is the second
c                                  derivative changes too rapidly for an 
c                                  adequate interval to be found (sdcerr 
c                                  remained small as h was decreased itmax times).
               inform = 3
               hphi = abs(hsave)
               hopt = hphi
               fdest = fdsave
               cdest = cdsave
               sdest = sdsave
               errbnd = hopt*abs(sdest)/2d0 + 2d0*epsa/hopt

            end if

         end if

      end if
c                                 end of chcore
      end

      subroutine lscrsh (nclin,nctotl,nactiv,nfree,n,lda,istate,kactiv,
     *                   tolact,a,ax,bl,bu,x,wx)
c----------------------------------------------------------------------
c     lscrsh  computes the quantities istate, kactiv, nactiv, and nfree 
c     associated with the working set at x.

c     an initial working set is selected. nearly-satisfied or violated 
c     bounds are added, and lastly general linear constraints are added.

c     values of istate(j)
c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c----------------------------------------------------------------------
      implicit none

      integer lda, n, nactiv, nclin, nctotl, nfree, istate(nctotl), 
     *        kactiv(n), i, imin, is, j, nfixed

      double precision a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                 wx(n), x(n), tolact, residl, resl, resmin, 
     *                 resu, toobig, ddot

      external ddot

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      call dcopy (n,x,1,wx,1)
c                                 initialize variables
      nfixed = 0
      istate(1:nctotl) = 0
c                                 initialize constraints
      nactiv = 0

      do j = n+1, nctotl

         if (bl(j).eq.bu(j)) then
            istate(j) = 3
            nactiv = nactiv + 1
            kactiv(nactiv) = j - n
         end if

      end do
c                                 attempt to add as many constraints as possible to 
c                                 the working set.
c                                 -----------------------
c                                 check if any bounds are violated or nearly 
c                                 satisfied. if so, add these bounds to the working set and set the
c                                 variables exactly on their bounds.
      do j = n, 1, -1

         if (nfixed+nactiv.lt.n) then

               is = 0

               if (wx(j)-bl(j).le.(1d0+abs(bl(j)))*tolact) is = 1

               if (bu(j)-wx(j).le.(1d0+abs(bu(j)))*tolact) is = 2

               if (is.gt.0) then

                  istate(j) = is

                  if (is.eq.1) wx(j) = bl(j)
                  if (is.eq.2) wx(j) = bu(j)

                  nfixed = nfixed + 1

               end if

         else

            exit

         end if

      end do
c                                 find the linear constraint with
c                                 smallest residual <= tolact and add it
c                                 to the working set. repeat until the working set
c                                 is complete or all remaining residuals are too large.
      if (nclin.gt.0 .and. nactiv+nfixed.lt.n) then
c                                 compute residuals for all constraints not in
c                                 working set.
         do i = 1, nclin
            if (istate(n+i).le.0) ax(i) = ddot (n,a(i,1),lda,wx)
         end do

         is = 1

         toobig = tolact + tolact

         do

            if (is.gt.0 .and. nfixed+nactiv.lt.n) then

               is = 0
               resmin = tolact

               do i = 1, nclin

                  j = n + i

                  if (istate(j).eq.0) then

                     resl = abs(ax(i)-bl(j))/(1d0+abs(bl(j)))
                     resu = abs(ax(i)-bu(j))/(1d0+abs(bu(j)))

                     residl = min(resl,resu)

                     if (residl.lt.resmin) then

                        resmin = residl
                        imin = i
                        is = 1
                        if (resl.gt.resu) is = 2

                     end if

                  end if

               end do

               if (is.gt.0) then

                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j = n + imin
                  istate(j) = is

               end if

            else

               exit

            end if

         end do

      end if

      nfree = n - nfixed
c                                 end of lscrsh
      end

      subroutine lsadds (unitq,inform,k2,nactiv,nz,nfree,nrank,nrejtd,
     *                   nres,ngq,n,ldzy,lda,ldr,ldt,istate,kactiv,kx,
     *                   condmx,a,r,t,res,gq,zy,w,c,s)
c----------------------------------------------------------------------
c     lsadds  includes general constraints 1 thru k2 as new rows of
c     the tq factorization stored in t, zy. if nrank is nonzero, the
c     changes in q are reflected in nrank by n triangular factor r such
c     that
c                         c  =  p (r) q,
c                                 (0)
c     where  p  is orthogonal.
c----------------------------------------------------------------------
      implicit none

      logical unitq

      integer inform, k2, lda, ldr, ldt, ldzy, n, nz, istate(*), k,
     *        nactiv, nfree, ngq, nrank, nrejtd, nres, jadd, l,
     *        i, iadd, ifix, iswap, kactiv(n), kx(n)

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                 s(n), t(ldt,*), w(n), zy(ldzy,*), condmx, 
     *                 dnrm2

      external dnrm2

      double precision wmach
      common/ cstmch /wmach(10)

      double precision asize, dtmax, dtmin
      common/ ngg008 /asize, dtmax, dtmin
c----------------------------------------------------------------------
c                                 estimate the condition number of the constraints that are not
c                                 to be refactorized.
      if (nactiv.eq.0) then
         dtmax = 0d0
         dtmin = 1d0
      else
         call scond (nactiv,t(nactiv,nz+1),ldt-1,dtmax,dtmin)
      end if

      do k = 1, k2

         iadd = kactiv(k)
         jadd = n + iadd

         if (nactiv.lt.nfree) then

            call lsadd (unitq,inform,ifix,iadd,jadd,nactiv,nz,nfree,
     *                  nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,a,r,
     *                  t,res,gq,zy,w,c,s)

            if (inform.eq.0) then
               nactiv = nactiv + 1
               nz = nz - 1
            else
               istate(jadd) = 0
               kactiv(k) = -kactiv(k)
            end if
         end if
      end do

      if (nactiv.lt.k2) then

c        some of the constraints were classed as dependent and not
c        included in the factorization. re-order the part of  kactiv
c        that holds the indices of the general constraints in the
c        working set. move accepted indices to the front and shift
c        rejected indices (with negative values) to the end.

         l = 0

         do k = 1, k2
            i = kactiv(k)
            if (i.ge.0) then
               l = l + 1
               if (l.ne.k) then
                  iswap = kactiv(l)
                  kactiv(l) = i
                  kactiv(k) = iswap
               end if
            end if
         end do

      end if

      nrejtd = k2 - nactiv
c                                 end of lsadds
      end

      subroutine badalf (alfa,n,x,x1,dx,char)
c---------------------------------------------------------------------
c hack to prevent linesearch from generating bad alfa values for 
c simplicial solution models
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j, n

      character char*1

      double precision alfa, x(*), x1(*), dx(*), newa, sum
c---------------------------------------------------------------------
      return
      if (.not.boundd(rids)) return

      newa = alfa
      sum = 0d0

      do j = 1, n

         sum = sum + x(j)

         if (x(j).lt.nopt(50)) then

            if (x(j).gt.-nopt(50)) then 
               x(j) = 0d0
            else if (dx(j).ne.0d0) then
               if (x1(j)/dx(j).lt.newa) newa = -x1(j)/dx(j)
            end if

            end if

         end do

         if (newa.lt.alfa.and.newa.ge.0d0) then

c           write (*,*) char,' from ',alfa,'to',newa,' rids',rids

            alfa = newa

            call dcopy (n,x1,1,x,1)
            call daxpy (n,alfa,dx,1,x,1)

         end if

          if (sum.gt.1d0) then
c             write (*,*) 'uh oh?',sum
          end if
c                                 end of badalf
      end 