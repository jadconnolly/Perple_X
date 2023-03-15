      subroutine aqsol2 (g0,gso,mo,mu,is,gamm0,lnkw,bad)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, it, iwarn

      logical bad, kill

      double precision c(l9), mo(*), mu(*), dg, dn,
     *                 d(l9), is, gamm0, g0(*), lnkw, 
     *                 gso(*),qb(l9), 
     *                 cst, mi, chg, dchg, dgam, dis, dm

      double precision gcpd, solve, aqact

      external gcpd, solve, aqact

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      logical abort
      common/ cstabo /abort

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save iwarn
      data iwarn /0/
c----------------------------------------------------------------------
      if (epsln.lt.nopt(34).or.abort.or.yf(1).eq.0d0) then
c                                  vapor, same as checking lnkw
         bad = .true.
         return

      end if
c                                  set default charge balance ion (aq_ion_H+, lopt(44)
      if (lopt(44)) then
         ion = ihy
      else
         ion = ioh
      end if
c
      cst = 0.2d0
c                                  if default choice fails switch to back-up choice
      do k = 1, 2
c                                 set up coefficients for mo(ion) equation
         g0(ion) = gcpd (aqst+ion,.false.)
c                                 compute solute properties
         do i = 1, aqct
c                                 dg is the solvent oxide potentials - g
            g0(i) = gcpd (aqst+i,.false.)
            qr(i) = q(i)/q(ion)
            qb(i) = (q(ion)-q(i))*q(i)
            dg = -g0(i) + qr(i)*g0(ion)

            kill = .false.

            do j = 1, kbulk

               dn = aqcp(j,i) - qr(i)*aqcp(j,ion)

               if (dn.eq.0d0) cycle

               if (isnan(mu(j))) then
c                                 kill any species that depends on
c                                 an undetermined chemical potential
c                                 unless NOT lopt(36), this allows
c                                 oxide components without redox because
c                                 g(H+) is not a function of mu(O2) [but
c                                 is not fail-safe].
                  if (aqcp(j,i).ne.0d0.and..not.lopt(36)) then
                     kill = .true.
                     exit
                  else
                     cycle
                  end if

               else if (cblk(j).eq.0d0.and..not.lopt(36)) then
c                                 this check is necessary because lp may
c                                 give a zero-amount solution for the chemical
c                                 potential of an absent component. the test
c                                 cannot be made with oxide components.
                  if (aqcp(j,i).ne.0d0.and.j.le.jbulk) then
                     kill = .true.
                     exit
                  end if

               end if

               dg = dg + dn * mu(j)

            end do
c                                 normalize by RT
            if (kill) then
               dg = 0d0
            else if (dg/rt.gt.nopt(57)) then
               bad = .true.
               return
            else
               dg = dexp(dg/rt)
            end if

            if (q(i).ne.0d0) then
c                                  this is now c(i)*a(ion)^(q(i)) = mo(i)*gamma(i)*q(i)
c                                  the rhs q(i) is because the eq to be solved will be
c                                  sum (q(i)*m(i)) = 0.
               d(i) = q(i)*dg

            else
c                                  neutral species assumed to be ideal, molality is
               mo(i) = dg

            end if

         end do
c                                  initialize iteration loop
         lnkw = (gso(ns) - g0(ioh))/rt

         if (d(ion).eq.0d0) then
c                                 no hydrogen or no oxygen
            bad = .true.
            return

         end if
c                                 starting guess for the reference ion
         mi = dexp(lnkw/2d0)
c                                 ideal limit
         gamm0 = 1d0
c                                 guess derivative = 0
         dgam = 0d0
c
         it = 0
c                                 newton raphson iteration loop
         do
c                                 sum of charges
            chg = 0d0
            is = 0d0

            do i = 1, ichg

               j = jchg(i)

               mo(j) = d(j) * gamm0**qb(j) / q(j) * mi**qr(j)

               chg = chg + mo(j) * q(j)

               is = is + mo(j) * q2(j)

            end do

            is = is/2d0
c                                 update gamma0
            gamm0 = dexp(adh*dsqrt(is)/(1d0 + dsqrt(is)) + 0.2d0*is)
c                                 derivative of the ionic strength
            dis = 0d0

            do i = 1, ichg

               j = jchg(i)

               dis = dis + q2(j)*d(j)*(gamm0**(qb(j)-1d0)*mi**qr(j)
     *                                *qb(j)*dgam + mi**(qr(j)-1d0)
     *                                *gamm0**qb(j)*qr(j))/q(j)
            end do

            dis = dis/2d0
c                                 the derivative of gamm0
            dgam = dis*(2d0*cst*is**(1.5d0) + 4d0*cst*is 
     *                 + 2d0*cst*dsqrt(is) + adh)
     *                 / (2d0*dsqrt(is)*(1d0+dsqrt(is))**2)
     *                 * dexp((cst*is**(1.5) + adh*dsqrt(is) * cst*is)
     *                       / (1d0 + dsqrt(is)))
c                                 finally the derivative of the net chg
            dchg = 0d0

            do i = 1, ichg

               j = jchg(i)

               dchg = dchg + d(j)*(gamm0**(qb(j)-1d0)*mi**(qr(j))
     *                            *qb(j)*dgam
     *                           + mi**(qr(j)-1d0)*gamm0**(qb(j))*qr(j))

            end do

            dm = -chg/dchg

            gamm0 = gamm0 + dgam*dm


            if (mi+dm.lt.0d0) then 
               dm = -mi/2d0
            end if

            mi = mi + dm

            it = it + 1

            if (mi.le.1d-50.or.mi.gt.1d3.or.it.gt.iopt(21)) then
               bad = .true.
               exit
            else if (dabs(chg).lt.nopt(50)) then
               bad = .false.
               return
            end if

         end do

      end do
c                                 failure is the only path here
      if (bad.and.iwarn.lt.11) then

         call warn (64,is,it,' ')

         call prtptx

         if (iwarn.eq.10) call warn (49,0d0,64,'AQSOLV')

         iwarn = iwarn + 1

      end if

      end

