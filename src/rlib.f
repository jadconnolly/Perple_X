c-----------------------------------------------------------------------

c RLIB - a library of subprograms called by VERTEX, FRENDLY, ACTCOR
c BUILD, ADIABAT, and ISOCHOR.

c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.

c Please do not distribute this source.

c-----------------------------------------------------------------------

      character*10 function gname (id)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c-----------------------------------------------------------------------
      if (id.lt.0) then
         gname = names(-id)
      else if (id.gt.0) then
         gname = fname(id)
      end if

      end

      double precision function gexces (id)
c-----------------------------------------------------------------------
c gexces evaluates the contribution to the gibbs energy of a pseudocompound
c arising from configurational entropy, excess properties, and dqf corrections
c for special cases (internal EoS's or Van Laar gexces does not return
c the excess energy, see routines fexces, gvlaar, and toop).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision exces
      common/ cst304 /exces(m3,k1)
c-----------------------------------------------------------------------

      gexces = exces(1,id) + t * exces(2,id) + p * exces(3,id)

      end

      double precision function gfluid (y)
c-----------------------------------------------------------------------
c gfluid returns the fugacities computed from a binary phase described
c by an internal EoS referenced by cfluid, the composition of the phase
c is y (the amount of component 2).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision fo2, fs2, y

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)
c-----------------------------------------------------------------------
      xco2 = y

      call cfluid (fo2,fs2)

      gfluid = r*t*((1d0-y)*f(1) + y*f(2))

      end


      subroutine fexces (id,dg)
c-----------------------------------------------------------------------
c gexces evaluates the contribution to the gibbs energy of a pseudocompound
c arising from dqf corrections, and excess properties computed from an
c internal fluid EoS
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision dg,f,fo2,fs2

      common/ cst11 /f(3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision exces
      common/ cst304 /exces(m3,k1)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      dg = exces(1,id) + t * exces(2,id) + p * exces(3,id)
c                                 fexces is only called for solution model
c                                 type 0 - a binary solution in which the
c                                 first endmember is always the second 
c                                 special component. this model cannot be
c                                 called if there are no special components?
      xco2 = pa(1)

      call cfluid (fo2,fs2)

      dg = dg + r*t*(pa(2)*f(1) + pa(1)*f(2))

      end

      recursive double precision function gcpd (id,proj)
c-----------------------------------------------------------------------
c gcpd computes the gibbs free energy of a compound identified by
c the arguement 'id' from the thermochemical parameters stored
c in the array 'thermo' which is located in common block 'cst1'.
c the parameters are: g(pr,tr) [thermo(1,id)], s(pr,tr) [thermo
c (2,id)], v(pr,tr) [thermo(3,id)], the coefficients of the extended
c maier-kelley heat capacity equation:

c           cp(pr,t)=a+b*t+c/(t*t)+d*t*t+e/t**(1/2)+f/t+g/t**3

c [thermo(4-10,id)], and the coefficients of the volumetric
c equations given in the program documentation (Eqs 2.1-2.3):

c        b8 = 0 =>

c           v(p,t) = v(pr,tr) + b2*(t-tr) + b4*(p-pr)
c                           + b6*(p-pr)**2 + b7*(t-tr)**2

c        -3d0 < b8 < 0 =>

c           v(p,t) = v(pr,tr) * exp[ b3*(t-tr) + b8*(p-pr) ]

c        b8 < -3d0 =>

c           p = 3 * k * f * (2f + 1)^(5/2) * (1 - 3/2 * f * (4 + b8))
c           f = 1/2 * ((v0/v)^(2/3) - 1)
c           v0 = v(pr,tr) + b1 * (t - tr)
c           k = -v0 / (b2 + b3 * (t - tr))

c        b8 > 0 =>

c           v(p,t) = v(pr,t) * [1 - b8*p/(b8*p+Kt)]^(1/b8)
c               Kt = b6 + b7 * t
c            alpha = b1 + b2*t + b3/t - b4/t^2 + b5/t^(1/2)

c Ghiorso parameters:

c v(J/bar) dv/dt(J/K/bar) dv/dp(J/bar^2) dv/dp/dt(J/bar^2/K) -K'

c Perple_X parameters (1st index in array thermo):

c    3            11           12             13              18

c N.B. for some reason Ghiorso chooses a reference T = 1673 for tbe
c mechanical EoS parameters, this is hardwired as variable trv.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, iwarn, oldid, j

      logical proj

      double precision ialpha, vt, trv, pth, vdp, vdpbm3, gsixtr,
     *                 gstxgi, fs2, fo2, kt, gval, gmake, gkomab, kp,
     *                 a, b, c, gstxlq, glacaz, v1, v2, gmet, gmet2,
     *                 gterm2, km, kmk, lnfpur, gaq, ghkf, lamla2

      external vdpbm3, gsixtr, gstxgi, gmake, gkomab, gstxlq, glacaz,
     *         gaq,    lnfpur, gmet, gmet2, gterm2, ghkf, lamla2

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer make
      common / cst335 /make(k10)

      integer eos
      common/ cst303 /eos(k10)

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      double precision f
      common/ cst11 /f(3)

      integer iam
      common/ cst4 /iam

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      double precision mmu
      common/ cst39 /mmu(i6)

      save kt,trv,iwarn,oldid
      data kt,trv,iwarn,oldid/0d0,1673.15d0,0,0/
c---------------------------------------------------------------------

      if (make(id).ne.0) then
c                                 the phase is a made phase, compute
c                                 and sum the component g's.
         gval = gmake (id)
         goto 999

      else if (eos(id).eq.5) then
c                                 sixtrude 05 JGR EoS
         gval = gsixtr (id)
         goto 999

      else if (eos(id).eq.6) then
c                                 stixrude JGI '05 Eos
         gval = gstxgi (id)
c                                 landau O/D
         if (ltyp(id).eq.4) then 
c                                 in the 2011 data this is only qtz, 
c                                 but neglects the effect of the clapeyron 
c                                 slope on the transition T. This gives 
c                                 nonsensical results if extrapolated to high
c                                 pressure, therefore the transition effect
c                                 was commented out Feb 9, 2022. Apparently
c                                 the effect was not accounted for from the 
c                                 initial implementation in perple_X and was
c                                 added April 3, 2021.
c           call lamla4 (dg,lmda(id))
c           gval = gval + dg

         else if (ltyp(id).eq.7) then 
c                                 in the 2021 relative to the low T phase,
c                                 used pointlessly for magnetic entropy of
c                                 almost all Fe-bearing endmembers. 
            gval = gval + lamla2(lmda(id))

         end if

         goto 999

      else if (eos(id).eq.11) then
c                                 stixrude EPSL '09 Liquid Eos
         gval = gstxlq (id)
         goto 999

      else if (eos(id).eq.12) then
c                                read SGTE data and evaluate EOS after Brosh '07,'08:
c                                Nastia's implementation; see also eos(id) 17
         gval = gmet (id)
         goto 999

      else if (eos(id).eq.14) then
c                                read SGTE data and evaluate EOS after Brosh'07,'08
c                                (modified for liquid carbon)
         gval = gterm2 (id)
         goto 999

      else if (eos(id).eq.15) then
c                                Anderson density extrapolation aqueous species EoS
         gval = gaq (id)
         goto 999

      else if (eos(id).eq.16) then
c                                DEW/HKF aqueous species formulation
         gval = ghkf (id)
         goto 999

      else if (eos(id).eq.17) then
c                                read SGTE data and evaluate EOS after Brosh '07,'08:
c                                ecrg's implementation after Saxena & Eriksson 2015
c                                (quasi-harmonic terms are dodgy)
         gval = gmet2 (id)
         goto 999

      end if
c                                 all other EoS's with Cp function
      gval = thermo(1,id)
c                                 -sdt
     *      + t * (thermo(2,id) - thermo(4,id) * dlog(t)
     *      - t * (thermo(5,id) + (thermo(7,id) - thermo(24,id)*t) * t))
     *      - (thermo(6,id) + thermo(10,id) / t) / t
     *      + thermo(8,id) * dsqrt(t) + thermo(9,id)*dlog(t)
c                                 vdp-ndu term:
      if (eos(id).eq.8) then
c                                 HP Tait EoS, einstein thermal pressure
         pth = thermo(11,id)*(1d0/(dexp(thermo(15,id)/t)-1d0)
     *         -thermo(19,id))

         v1 = 1d0 + (p -pth)*thermo(17,id)
         v2 = 1d0 + (pr-pth)*thermo(17,id)

         if (v1.lt.0d0) then
c                                 destabilize the phase
            vdp = thermo(3,id)**2*p

            if (iwarn.le.5.and.oldid.ne.id) then
               call warn (60,t,1,names(id))
               iwarn = iwarn + 1
               oldid = id
               if (iwarn.eq.5) call warn (49,t,60,'GCPD_HP_Tait_I')
            end if

          else if (v2.lt.0d0) then
c                                 continue on the assumption that
c                                 v2 is small zero it and continue
             vdp = (thermo(16,id)*(v1**thermo(18,id)
     *             /thermo(20,id)-p+pr)+p-pr)*thermo(3,id)

            if (iwarn.le.5.and.oldid.ne.id) then
               call warn (60,t,2,names(id))
               iwarn = iwarn + 1
               oldid = id
               if (iwarn.eq.5) call warn (49,t,60,'GCPD_HP_Tait_II')
            end if

          else
c                                 int(vdp,p=Pr..Pf)
             vdp = (thermo(16,id)*(
     *            (v1**thermo(18,id) - v2**thermo(18,id))
     *            /thermo(20,id)-p+pr)+p-pr)*thermo(3,id)

          end if

      else if (eos(id).eq.9) then
c                                 True tait used for melt endmembers
c                                 kt = b6 + b5*(T-Tr)
c                                 vt = v0*exp(b1*(T-Tr))
          kt = thermo(16,id) + thermo(15,id) * (t-tr)
c                                 tait parameters, "c" is 1 - tait c
          a = thermo(19,id)/(thermo(19,id)+kt*thermo(17,id))
          b = thermo(18,id)/kt-thermo(21,id)
          c = 1d0 - (thermo(19,id)+kt*thermo(17,id))
     *             /(thermo(20,id)-kt*thermo(17,id))
c                                 int(vdp,p=Pr..Pf)
          vdp = ((((p*b+1d0)**c-(Pr*b+1d0)**c)/b/c+pr-p)*a-pr+p)*
c                                 vt
     *          thermo(3,id)*dexp(thermo(11,id)*(t-tr))

      else if (eos(id).eq.10) then
c                                 ideal gas EoS
         vdp = r*t*dlog(p/pr)

      else if (eos(id).eq.13) then
c                                 komabayashi/omori polynomials, murnaghan EoS
         vt = thermo(3,id) * dexp((thermo(11,id) + thermo(12,id)*t)*t +
     *        thermo(13,id)*dlog(t) + thermo(14,id)/t + thermo(23,id))
         kt = 1d0 / (thermo(15,id) + t*(thermo(16,id) + t*(thermo(17,id)
     *                             + thermo(18,id)*t)))
c                                 k'(T)
         kp = thermo(19,id) + thermo(20,id)*(t-tr)*dlog(t/tr)

         km = kp - 1d0
         kmk = km/kp

         vdp = vt * kt**(1d0/kp)/km
     *            * ((kt + kp*p )**kmk - (kt + kp*pr)**kmk)

      else if (thermo(18,id).eq.0d0) then
c                                 normal polynomial:
          vdp =  p * (thermo(3,id)
     *               + (thermo(17,id) * t + thermo(12,id)) * t
     *               + (thermo(16,id) * p + thermo(14,id)) * p)

      else if (thermo(18,id).gt.0d0) then
c                                 murnaghan EoS:
c                                 int(alpha,T=Tr..T)
         ialpha = (thermo(11,id) + thermo(12,id)*t)*t +
     *            thermo(13,id)*dlog(t) + thermo(14,id)/t +
     *            thermo(15,id)*dsqrt(t) + thermo(23,id)

         if (lopt(8)) then
c                                 use holland & powell's approximate form
            vt = thermo(3,id)*(1d0 + ialpha)

         else
c                                 v(t,pr), correct form
            vt = thermo(3,id)*dexp(ialpha)

         end if
c
         if (lopt(4)) then
c                                 compute kt using Anderson-Gruneisen parameter
c                                 and expansivity ala Helffrich & Connolly 2009.
            kt = thermo(16,id)*dexp(-thermo(21,id)*ialpha)

         else

            kt = thermo(16,id) + thermo(17,id)*t
c                                 a ****wit has entered a ridiculous
c                                 temperature
            if (kt.lt.0d0) then

               if (iwarn.lt.5.and.id.ne.oldid) then
                  call warn (46,t,id,names(id))
                  iwarn = iwarn + 1
                  oldid = id
                  if (iwarn.eq.5) call warn (49,t,46,'GCPD_Murnaghan')
               end if
c                                 destabalize the phase
               gcpd = thermo(3,id)**2*p

               return

            end if

         end if
c                                 Murnaghan EoS:
c                                 vdp = V(1,T)*KT**(1/K')/(K'-1)
c                                 * [(KT+K'*p)**(1-1/K') -
c                                 (KT+K'*pr)**(1-1/K')]
         vdp = vt * kt**(1d0/thermo(18,id))/thermo(22,id)
     *            *((kt+thermo(18,id)*p)**thermo(19,id)
     *            -(kt+thermo(20,id))**thermo(19,id))

      else if (thermo(18,id).lt.-3d0) then
c                                 3rd order Birch-Murnaghan
c                                 Ghirso Eos is a special case
c                                 indicated by thermo(16,id) = K0 = 0
         if (thermo(16,id).eq.0d0) then
c                                 assume ghiorso's expressions for v0 and k (KT)
c                                 vt = Volume at 1 bar, T
            vt = thermo(3,id) + thermo(11,id) * (t - trv)

            kt = -vt / (thermo(12,id) + thermo(13,id) * (t - trv) )

         else
c                                 int(alpha,T=Tr..Tf)
            ialpha = (thermo(11,id) + thermo(12,id)*t)*t +
     *               thermo(13,id)*dlog(t) + thermo(14,id)/t +
     *               thermo(15,id)*dsqrt(t) + thermo(23,id)
c                                 v(t,pr)
            vt = thermo(3,id)*dexp(ialpha)

            if (lopt(4)) then
c                                 compute kt using Anderson-Gruneisen parameter
c                                 and expansivity ala Helffrich & Connolly 2009.
               kt = thermo(16,id)*dexp(-thermo(21,id)*ialpha)

            else

               kt = thermo(16,id) + thermo(17,id)*t

            end if

         end if
c                                 a ****wit has entered a ridiculous
c                                 temperature
         if (kt.lt.0d0.or.vt.lt.0d0) then

            if (iwarn.lt.5.and.oldid.ne.id) then
               call warn (46,t,id,names(id))
               iwarn = iwarn + 1
               oldid = id
               if (iwarn.eq.5) call warn (49,t,46,'GCPD_BM3')
            end if
c                                 destabilize the phase
            vdp = thermo(3,id)**2*p

         else

            vdp = vdpbm3 (vt,kt,thermo(18,id))

         end if

      else
c                                 gottschalk.
         vdp = thermo(11,id)*dexp(thermo(13,id)*t)*
     *               (1d0 - dexp(thermo(18,id)*(p - pr)))

      end if

      gval = gval + vdp

c                                 check for transitions and Landau O/D:
      if (ltyp(id).ne.0) call mtrans (gval,vdp,id)
c                                 check for BERMAN temperature dependent
c                                 order/disorder:
      if (idis(id).ne.0) call disord (gval,idis(id))

c                                 fluids in the saturated
c                                 or thermodynamic composition space, these
c                                 are not used for mixture props.
      if (eos(id).gt.100) then

         if (eos(id).eq.201.or.eos(id).eq.202) then
c                                 species has been identified as a special composant
c                                 and eos is set by ifug
            if (eos(id).eq.201) then

               if (iam.ne.5) xco2 = 0d0
               call cfluid (fo2,fs2)
               gval = gval + r*t*f(1)

            else

               if (iam.ne.5) xco2 = 1d0
               call cfluid (fo2,fs2)
               gval = gval + r*t*f(2)

            end if

         else if (eos(id).lt.118) then
c                                 call appropriate pure fluid EoS
            gval = gval + r*t * lnfpur(eos(id))

         else if (eos(id).ge.600.and.eos(id).le.603) then
c                                 komabayashi & fei (2010) EoS for Fe
            gval = gkomab(eos(id),id,vdp)

         else if (eos(id).eq.605) then
c                                 Stoichiometic O rk fluid
            xco2 = 0d0
            call cfluid (fo2,fs2)
c                                 real O fluid (O and O2 species)
            gval = gval + r*t*f(1)
c                                 this is -RT(lnk2+lnk3)/2 (rksi5 k's)
c    *         -0.3213822427D7 / t + 0.6464888248D6 - 0.1403012026D3*t

         else if (eos(id).ge.610.and.eos(id).le.637) then
c                                 lacaze & Sundman (1990) EoS for Fe-Si-C alloys and compounds
c                                 Xiong et al., 2011 for Fe-Cr alloys
            gval = gval + glacaz(eos(id)) + vdp + thermo(1,id)

c        else if (eos(id).eq.800) then

c           pgpa = p/1d4
c           call interpolator (pgpa, t, gval, vt, a, b, err)

c           if (err) then 
c              write (*,*) 'interpolator set err true'
c              call prtptx
c           else 
c              write (*,'(a,g14.6)') 'from interpolator, g = ',gval
c              write (*,'(a,g14.6)') 'from interpolator, v = ',vt
c              write (*,'(a,g14.6)') 'from interpolator, g = ',a
c              write (*,'(a,g14.6)') 'from interpolator, g = ',b
c           end if

         end if

      end if
c                                 kill melt endmembers if T < T_melt
999   if (ifp(id).lt.0.and.t.lt.nopt(20)) gval = gval + 1d6
c                                 do legendre transform if proj
      if (proj) then
c                                 mobile components
         do j = 1, jmct
            gval = gval - vnumu(j,id) * mmu(j)
         end do

      end if

      gcpd = gval

      end

      subroutine zeroi (iarray,index,ivalue,n)

      implicit none

      integer n, iarray(n), index, ivalue, i

      do i = 1, index
         iarray(i) = ivalue
      end do

      end

      subroutine xchk (xmin, xmax, xinc, tname)

      implicit none

      double precision xmin, xmax, xinc

      character tname*10

      if (xmax.gt.1d0) then
         call warn (21,xmax,1,tname)
         xmax = 1d0
      end if

      if (xmin.lt.0d0) then
         call warn (22,xmin,1,tname)
         xmin = 0d0
      end if

      if (xmax.lt.xmin) then
         call warn (23,xmax,1,tname)
         xmax = 1d0
         xmin = 0d0
      end if

      if (xinc.le.0d0) then
         call warn (23,xinc,1,tname)
         xinc = 1d0
      end if

      end

      logical function badz (z)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision z
c----------------------------------------------------------------------
      if (z.gt.-nopt(50).and.z.le.nopt(55)) then
         badz = .false.
      else
         badz = .true.
      end if

      end

      subroutine loadit (id,make,nchk)
c---------------------------------------------------------------------
c loadit loads descriptive data for phases and species (name,comp,
c and therm) into the appropriate arrays (names,comps,thermo,vf,
c and vs).  the arguement 'id' indexes the phase in the arrays.
c note that loadit also computes reference state constants which
c are dependent on the state function being used and on its
c analytical expression.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,i,j,k

      logical make, nchk

      double precision gzero
      external gzero

      double precision z(14),smax,t0,qr2,vmax,dt,g1,g2

      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ipoint,kphct,imyn
      common/ cst60  /ipoint,kphct,imyn

      character*8 name
      common/ csta6 /name

      double precision emodu
      common/ cst318 /emodu(k15)

      integer ic
      common/ cst42 /ic(k0)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ikp
      common/ cst61 /ikp(k1)

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer eos
      common/ cst303 /eos(k10)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer junk
      double precision del, rand
      common/ cst321 /del(11,k10),rand(12,k10),junk

      double precision delta
      common/ cst325 /delta(11)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      character*8 eoscmp
      common/ cst98 /eoscmp(2)

      integer iam
      common/ cst4 /iam

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision fwt
      common/ cst338 /fwt(k10)
c---------------------------------------------------------------------

      if (id+1.gt.k10) call error (1,0d0,id+1,'k10')

      ipoint = iphct
c                               check for duplicates
      if (nchk) then
         do i = jmct + 1, iphct
            if (name.ne.names(i)) cycle
            call error (73,g1,i,name)
         end do
      end if
c                               load name and phase flag
      names(id) = name
c                               moduli flags, indicate availability of
c                               bulk and shear modulus

c                               if ikind = 0, no explicit moduli
c                                        = 1, just shear
c                                        = 2, shear and bulk
c                                        = 3, both
      iemod(id) = ikind

      ikp(id) = 0
      ifp(id) = 0
      idis(id) = 0
      lmda(id) = 0

      eos(id) = ieos

      if (lopt(7)) then

         do k = 1, ispec

            if (name.ne.cmpnt(idspe(k)).and.name.ne.eoscmp(k)) cycle
c                                 this is an awful mess, if there is a saturated
c                                 phase (ifct > 0) then ufluid will call the eos
c                                 identified by ifug irrespective of the eos value.
            if (ifct.eq.0.or.iam.eq.5) then

               if (ieos.gt.100.and.(iam.lt.5.or.iam.eq.15)) 
     *            call warn (56,r,k,name)
c                                 there is no saturated phase
c                                 assign it the default molecular fluid eos
               eos(id) = 200 + k

            else if (k.eq.1.and.idfl.ne.2.or.
     *               k.eq.2.and.idfl.ne.1) then
c                                 saturated phase, and it's not component(k), ergo
c                                 will be computed by ufluid. this only will work
c                                 for ispec < 3.
               eos(id) = ieos

            end if
c                                 set fluid flag, this flag is
c                                 used only to match fluid endmembers
c                                 with fluid pseudocompounds
            ifp(id) = 1
c                                 gflu used to indicate whether a fluid is
c                                 in the calculation. not clear why gflu
c                                 is set if saturated phase (formerly it
c                                 was set only if saturated phase was
c                                 ABSENT, 1/5/2017).
            gflu = .true.

            exit

         end do

      end if
c                                 use ieos flag to signal melt endmembers
c                                 in ifp array, this is only used by gcpd.
      if (eos(id).eq.3.or.eos(id).eq.9.or.eos(id).eq.11) then
c                                 liquid
         ifp(id) = -1

      else if (eos(id).eq.10.or.eos(id).gt.100.and.eos(id).le.202.or.
     *         eos(id).eq.605) then
c                                 fluid
         gflu = .true.
         ifp(id) = 1

      end if
c                                 load stoichiometry of components.
      fwt(id) = 0

      do i = 1, icomp
         cp(i,id) = comp(ic(i))
         fwt(id) = fwt(id) + cp(i,id)*atwt(i)
      end do
c                               convert to kg/mol
      fwt(id) = fwt(id)/1d3
c                               compositional array for frendly
      if (iam.eq.5.and.id.le.k5) then
         do i = 1, k0
            cp0(i,id) = comp(i)
         end do
      end if
c                               and just mobile components
      do i = 1, jmct
         vnumu(i,id) = comp(ic(i+jprct))
      end do

      if (make) return
c                               if aqueous solute species store name and
c                               compositional data in special arrays (in
c                               principle may need vnumu as well).
      if (ieos.eq.15.or.ieos.eq.16) then
c                                aqst is only properly initialized by programs that call input2,
c                                build sets it to -1 as a flag.
         if (aqst.eq.-1) aqst = iphct - 1

         aqct = iphct - aqst

         if (aqct.gt.l9) call error (1,r,aqct,'l9 (max aq species)')

         aqnam(aqct) = name
         aqtot(aqct) = 0d0

         do k = 1, icomp
            aqcp(k,aqct) = comp(ic(k))
            if (k.le.icp) aqtot(aqct) = aqtot(aqct) + comp(ic(k))
         end do
c                               locate H+/OH-, at this point HOH is in thermo(13), after
c                               conver it's in thermo(21)
         if (thermo(13,k10).eq.1d0) then
            ihy = aqct
         else if (thermo(13,k10).eq.2d0) then
            ioh = aqct
         end if

      end if
c                               load elastic props if present
      if (iemod(id).ne.0) then

         do i = 1, k15

            if (i.eq.3.and.emodu(i).gt.0d0) then
 
               call warn (99,0d0,i,'The T derivative of K for '//name//
     *                   ' is > 0, this is possible, but anomalous.')

            else if (i.eq.6.and.emodu(i).gt.0d0) then 

               call warn (99,0d0,i,'The T derivative of mu for '//name//
     *                   ' is > 0, this is possible, but anomalous.')

            else if (i.eq.2.and.emodu(i).lt.0d0) then 

               call warn (99,0d0,i,'The P derivative of mu for '//name//
     *                   ' is < 0, this is improbable. ')

            else if (i.eq.5.and.emodu(i).lt.0d0) then

               call warn (99,0d0,i,'The P derivative of K for '//name//
     *                   ' is < 0, this is improbable.')

            end if

            emod(i,id) = emodu(i)

         end do

      end if
c                               compute reference state constants,
c                               etc..
      do i = 1, k4
         thermo(i,id) = thermo(i,k10)
      end do
c                               load errors for MC calculations
      do i = 1, 11
         del(i,id) = delta(i)
      end do

      call conver (
c                               g0, s0, v0
     *             thermo(1,id),thermo(2,id),thermo(3,id),
c                               c1-c8
     *             thermo(4,id),thermo(5,id),thermo(6,id),
     *             thermo(7,id),thermo(8,id),thermo(9,id),
     *             thermo(10,id),thermo(24,id),
c                               b1-b12
     *             thermo(11,id),thermo(12,id),thermo(13,id),
     *             thermo(14,id),thermo(15,id),thermo(16,id),
     *             thermo(17,id),thermo(18,id),thermo(19,id),
     *             thermo(20,id),thermo(21,id),thermo(22,id),
c                               b13 on return
     *             thermo(23,id),
c                               ref stuff
     *             tr,pr,r,eos(id))

      if (tr.eq.0d0) then
         thermo(1,id) = thermo(1,k10)
         thermo(2,id) = thermo(2,k10)
      end if
c                              lmda transitions:
      if (ilam.ne.0) then

         lamin = lamin + 1

         if (lamin.gt.k9) call error (1,0d0,lamin,'k9')

         if (jlam.eq.5) then
c                                 holland and powell, bragg-williams model:
c                                 enthalpy change of disordering
            therlm(1,1,lamin) = tm(1,1) - pr*tm(2,1)
c                                 volume change of disordering
            therlm(2,1,lamin) = tm(2,1)
c                                 excess enthalpy
            therlm(3,1,lamin) = tm(3,1)
c                                 excess volume
            therlm(4,1,lamin) = tm(4,1)
c                                 n
            therlm(5,1,lamin) = tm(5,1)
c                                 fac - unused?
            therlm(6,1,lamin) = tm(6,1)
c                                 n+1
            therlm(7,1,lamin) = tm(5,1) + 1d0
c                                 f
            therlm(8,1,lamin) = tm(5,1)/(tm(5,1) + 1d0)

         else if (jlam.eq.4.or.jlam.eq.7) then
c                                 holland and powell, landau model:
c                                 4 - relative to the high T phase
c                                 7 - relative to the low T phase (stixrude 2021).
            do j = 1, ilam

               smax = tm(2,j)
               t0 = tm(1,j)
               vmax = tm(3,j)

               therlm(1,j,lamin) = t0
               therlm(2,j,lamin) = smax
c                                 this makes therlm(3) dt/dp
               therlm(3,j,lamin) = vmax/smax

               if (jlam.eq.4) then 
                  qr2 = dsqrt (1d0 - tr/t0)
c                                 PX ds5 landau
                  therlm(4,j,lamin) = (2d0*t0 + tr)*qr2/3d0
c                                 TC ds6 landau
                  therlm(7,j,lamin) = t0*(qr2 - qr2**3/3d0)
                  therlm(8,j,lamin) = qr2
c                                 Vdp coefficient
                  therlm(6,j,lamin) = vmax*qr2/thermo(3,k10)
               end if

            end do

         else if (jlam.eq.1) then
c                              ubc:
            do j = 1, ilam

               therlm(1,j,lamin)=tm(1,j)*tm(1,j)
               therlm(2,j,lamin)=tm(2,j)*tm(2,j)
               therlm(9,j,lamin)=tm(1,j)*tm(2,j)

               do k = 3, 8
                  therlm(k,j,lamin)=tm(k,j)
               end do
            end do


         else if (jlam.eq.2.or.jlam.eq.3) then
c                              helgeson:
            p = pr
            lmda(id) = lamin
            ltyp(id) = jlam
c                              set special eos flag to zero to prevent
c                              gzero from including special terms (currently
c                              only ieos 605).
            eos(id) = 0
c                                 now convert paramters:
            do k = 1, ilam
c                                 load into therlm:
               therlm(1,k,lamin) = tm(1,k)
               therlm(2,k,lamin) = tm(2,k)
c                                 c1-c7 cp coefficients
               do j = 4, 10
                  therlm(j+1,k,lamin) = tm(j,k)
               end do
c                                 the c8 heat capacity coefficient
               therlm(13,k,lamin) = tm(11,k)

               t = tm(1,k)
c                                 temporary counter value
               lct(id)  = k-1
c                                 g at trt:
               therlm(12,k,lamin) = gzero(id)
c                                 delta v trans:
               therlm(4,k,lamin) = 0d0
               if (tm(2,k).ne.0d0) therlm(4,k,lamin) = tm(3,k)/tm(2,k)
c                                 s + st at trt:
               if (t*1d-3.lt.1d0) then
                  dt = 1d0
               else
                  dt = 1d-3*t
               end if

               t = tm(1,k) + dt
               g1 = gzero(id)

               t = tm(1,k) - dt
               g2 = gzero(id)

               therlm(3,k,lamin) = tm(3,k) - (g1 - g2)/2d0/dt
c                              streamline the eos:
               do j = 1, 13
                  z(j) = 0d0
               end do

               call conver (
c                                 g0,s0, dummy
     *                     therlm(12,k,lamin),therlm(3,k,lamin), z(1),
c                                 c1-c8
     *                     therlm(5,k,lamin), therlm(6,k,lamin),
     *                     therlm(7,k,lamin), therlm(8,k,lamin),
     *                     therlm(9,k,lamin), therlm(10,k,lamin),
     *                     therlm(11,k,lamin),therlm(13,k,lamin),
c                                 dummies (b1-b12)
     *                     z(14),z(2),z(3),z(4),z(5),z(6),z(7),z(8),
     *                     z(9),z(10),z(11),z(12),z(13),
c                                 ref stuff
     *                     tm(1,k),pr,r,0)

            end do

            eos(id) = ieos

        else if (jlam.eq.6) then

           do k = 1, ilam
c                              SGTE format
c                              load into therlm:
              therlm(1,k,lamin) = tm(1,k)
              therlm(2,k,lamin) = tm(2,k)

              do j = 4, 14
c                              nastia originally had t12-t15 here, but there
c                              doesn't appear to be any data with more than
c                              the first 11 parameters
c                              ----------------------------------------------
c                              added t12 (t^3), t13 (sqrt(t)), t14 (ln(t)) for eleanor.
c                              note these DO NOT (and DID NOT) correspond to the coefficients
c                              in Eq 4 as currently stated the thermo data file html.
c                              JADC, 12/3/2017.
                 therlm(j+1,k,lamin) = tm(j,k)

              end do

           end do

         else

            write (*,*) 'no such transition model'
            call errpau

         end if

         lmda(id) = lamin
         lct(id)  = ilam
         ltyp(id) = jlam

      end if
c                              t dependent order: load berman and brown
c                              parameters, (this should be cleaned up)
      if (idiso.ne.0) then
         idsin = idsin + 1
         idis(id) = idsin
         do j = 1, m8
            therdi(j,idsin) = td(j)
         end do
      end if

      end

      subroutine univeq (i,ier)
c---------------------------------------------------------------------
c univeq computes the equilibrium condition of a univariant reaction.
c using a finite difference newton-raphson approximation. the method
c will fail when the derivative of the state function, with respect to
c the state variable v(i), goes to infinity.
c---------------------------------------------------------------------
      implicit none

      double precision vi, delv, del, u, gr, b, xg

      integer i,j,ier

      include 'perplex_parameters.h'

      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision ddv,vmax,vmin
      common/ cst9 /vmax(l2),vmin(l2),ddv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      ier = 0

      vi = v(i)
      del = delt(i)
      b = blim(i)
      u = ulim(i)
c                                 phase composition special case:
      if (i.eq.3) then
         if (vi.lt.1d1*del) then
            del = dabs(vi)/1d1
         else if (1d0-vi.lt.1d1*del) then
            del = dabs(1d0-vi)/1d1
         end if
      end if

      if (vi+dabs(del).gt.u.or.vi-dabs(del).lt.b) goto 30

      do j = 1, 100

         call grxn (gr)

         v(i) = vi + del
         call incdep (i)

         call grxn (dgr)
         xg = dgr
         dgr = dgr - gr

         if (dgr.eq.0d0) exit

         delv = gr*del/dgr

         if (dabs(delv/ddv(i)).gt.1d0) then

c            v(i) = vi - del
c            del = del*ddv(i)/delv/2d0
c            if (dabs(del/delt(i)).lt.1d-6) goto 30
c            cycle
c                                  changed 7/13/2014
             delv = dabs(delv)/delv * ddv(i)

         end if

         vi = vi - delv

         if (vi+dabs(del).gt.u.or.vi-dabs(del).lt.b) goto 30

         v(i) = vi

         call incdep (i)
c                                 used to be on the basis of utol,
c                                 but this allows uniform relative error.
         if (dabs(delv).lt.del) return

      end do
c                                 error code 1: iv and dv must be
c                                 switched
      ier = 1
      return
c                                 error code 2: value too far out of
c                                 range, refine iv increment.
30    ier = 2
      end

      subroutine unver (g,s,v,a,b,c,d,e,f,gg,c8,b1,b2,b3,b4,b5,b6,b7,b8,
     *                  b9,b10,b11,tr,pr,ieos)
c----------------------------------------------------------------------
c convert thermodynamic equation of state from a 0-0 reference state
c to a pr-tr reference state.

c corrections made corresponding to PJ's corrections in conver,
c June 16, 2004.
c----------------------------------------------------------------------
      implicit none

      integer ieos

      include 'perplex_parameters.h'

      double precision v,gg,a,b,c,d,e,f,g,b1,b2,b4,b5,b6,b7,b8,pr,tr,s,
     *                 c8,b3,b9,b10,b11
c----------------------------------------------------------------------
c                               Stixrude's EoS, Aq, CALPHAD exit without
c                               doing anything
      if (ieos.eq. 5.or.ieos.eq. 6.or.ieos.eq.11.or.ieos.eq.12.or.
     *    ieos.eq.14.or.ieos.eq.15.or.ieos.eq.17) then

          return

      else if (ieos.eq.16) then
c                                 HKF electrolyte model.
          b3 = b11
          return

      end if

      c8 = 12d0 * c8
      gg = 6d0 * gg
      e  = e / 4d0
      d  = 6d0 * d
      c  = 2d0 * c

      if (b8.eq.0d0) then
c                                normal vdp term:
         b6 = 3d0 * b6
         b4 = 2d0 * (b4 + b6 * pr)
         b2 = b2 + 2d0 * b7 * tr
         b  = 2d0 * (b - b7*pr)
         v  = v + b2 * tr + b4 * pr - b6 * pr * pr - b7 * tr * tr
         s  = -1d0 * ( s - ( a - b2 * pr  + a * dlog(tr) + b * tr
     *        - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 / 3d0
     *        + b7 * 2d0 * pr * tr ) )
         g  = g - ( s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *        - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *        - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *        - c8 * tr**4 / 4d0
     *        - v * pr + b2 * tr * pr + b4 * pr * pr / 2d0
     *        - b6 * pr**3/3d0 - b7 * tr * tr * pr )
      else

         b  = 2d0 * b

         s  = -1d0 * ( s - ( a + a * dlog(tr) + b * tr
     *        - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 / 3d0 ) )
         g  = g - ( s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *        - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *        - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *        - c8 * tr**4 / 4d0)

         if (ieos.eq.13) then 
c                                 komabayashi 2006
            b4 = -b4
            b2 = 2d0*b2

         else if (b8.gt.0d0.or.(b8.le.-3d0.and.b6.ne.0d0)) then
c                                 murnaghan or bm3:
            b2 = 2d0 * b2
            b4 = -b4
            b5 = b5 / 2d0
c                                 convert b6 back to K(Tr)
            b6 = b6 - b7*tr

         else if (b8.le.-3d0.and.b6.eq.0d0) then
c                                 ghirso, do nothing.

         else
c                                 v = f(exponential(beta))
c                                 zero b1, so users don't get confused.
            b1 = 0d0

         end if

      end if

      end

c----------------------------------------------------------------------

c PLIB -  subroutines common to FRENDLY and CONVEX.

c----------------------------------------------------------------------

      subroutine grxn (gval)
c----------------------------------------------------------------------
c compute free energy change of a stoichiometric rxn
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer j

      double precision gval, gphase, gproj

      external gphase, gproj

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iam
      common/ cst4 /iam
c---------------------------------------------------------------------
      gval = 0d0

      if (iam.eq.5) then
c                                 FRENDLY:
         do j = 1, iphct
            gval = gval + vnu(j) * (gphase(j) + r * t * dlog(act(j)))
         end do

      else
c                                 CONVEX:
c                                 no saturated phase components and no
c                                 saturated components:
         if (iffr.ne.1.or.isr.ne.1) call uproj
c                                 note that this call to uproj makes a
c                                 subsequent call in gall redundant if
c                                 sfol1 is used to trace a univariant
c                                 curve.
         do j = 1, ivct
            gval = gval + vnu(j) * gproj (idr(j))
         end do 

      end if

      end

      subroutine slope (iv1,iv2,s)
c---------------------------------------------------------------------
c slope computes the slope (div(1)/div(2)) of a univariant
c equilibria by approximating the derivatives dg/div(1) and
c dg/div(2) by finite differences.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer iv(2),iv1,iv2,i

      double precision dg(2),gr,s,gval

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      iv(1) = iv1
      iv(2) = iv2
      call grxn (gr)

      do i = 1, 2
         v(iv(i)) = v(iv(i)) + delt(iv(i))
         call incdep (iv(i))
         call grxn (gval)
         dg(i) = (gval-gr)/delt(iv(i))
c                                 the possibility exists here
c                                 for the value of an intensive
c                                 parameter to excede allowed
c                                 limits.  although this is
c                                 unlikely a limit test like the
c                                 one done in 'univeq' could be
c                                 used.
         v(iv(i)) = v(iv(i)) - delt(iv(i))
         call incdep (iv(i))
      end do

      s = -dg(2) / dg(1)

      end

      subroutine switch (div,ivi,ivd,jer)
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ivi,ivd,jer,iovd

      double precision div, s

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr
c---------------------------------------------------------------------
c                                 reset intensive variables
      call reptx
c                                 determine the sign for the
c                                 ivd increment
      call slope(ivd,ivi,s)
      jer = 0
c
      if (s.eq.0d0) then
c                                 a zero slope shouldn't occur
c                                 someplace is a bug
         jer = 1
         return
      end if

      div = s * div

      iovd = ivd
      ivd = ivi
      ivi = iovd

      end

      subroutine reptx
c-----------------------------------------------------------------------
c reset - resets indendent potentials to the last known stable condition
c along a univariant curve
c-----------------------------------------------------------------------
      implicit none

      double precision ptx

      integer ipt2

      include 'perplex_parameters.h'

      common/ cst32 /ptx(l5),ipt2

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      v(iv(1)) = ptx(ipt2-1)
      v(iv(2)) = ptx(ipt2)
      call incdp0

      end

      subroutine unlam (tm,id)
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ilam,id,jd,i,j,k

      double precision tm(m7,m6), z(12), g1, g0, s0, gcpd

      external gcpd

      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer eos
      common/ cst303 /eos(k10)
c-----------------------------------------------------------------------
      if (ltyp(id).eq.0) return

      jd = lmda(id)

      do i = 1, m7
         do j = 1, m6
            tm(i,j) = 0d0
         end do
      end do

      if (ltyp(id).eq.5) then
c                                 Bragg-Williams model:
          do j = 1, 6
             tm(j,1) = therlm(j,1,jd)
          end do

          tm(1,1) = tm(1,1) + pr*tm(2,1)

      else if (ltyp(id).eq.4) then
c                                 Landau model
c                                 HP Landau model
         do j = 1, lct(id)
            tm(1,j) = therlm(1,j,jd)
            tm(2,j) = therlm(2,j,jd)
            tm(3,j) = therlm(3,j,jd) * tm(2,j)
         end do

      else if (ltyp(id).eq.1) then
c                                 UBC
         do j = 1, lct(id)
            tm(1,j) = dsqrt (therlm(1,j,jd))
            tm(2,j) = dsqrt (therlm(2,j,jd))
         end do

      else if (ltyp(id).eq.2.or.ltyp(id).eq.3) then
c                                 Helgeson generic, maybe q/coe too.
         p = pr
c                                 temporary counter
         ilam = lct(id)

         do i = ilam, 1 , -1
c                                 get s transition:
c                                 load into therlm:
            tm(1,i) = therlm(1,i,jd)
            tm(2,i) = therlm(2,i,jd)

            do j = 5, 11
               tm(j-1,i) = therlm(j,i,jd)
            end do

            tm(11,i) = therlm(13,i,jd)

            t = tm(1,i)
c                              set transition type to null
c                              for call to gphase
            lct(id) = i - 1
c                             -s at trt, this should be
c                             changed to centered rel diff:
            g1 = gcpd (id,.false.)
            t = t + 1d-3

            tm(3,i) =  (gcpd (id,.false.) - g1)/1d-3

            g0 = therlm(12,i,jd)
            s0 = therlm(3,i,jd)

            do k = 1, 9
               z(k) = 0d0
            end do

            call unver (
c                                 g0, s0, dummy
     *                  g0,s0,z(1),
c                                 c1-c8
     *                  tm(4,i),tm(5,i),tm(6,i),tm(7,i),
     *                  tm(8,i),tm(9,i),tm(10,i),tm(13,i),
c                                 dummies
     *                  z(1),z(2),z(3),z(5),z(6),z(7),z(8),z(9),z(10),
     *                  z(11),z(12),
c                                 ref stuff
     *                  tm(1,i),pr,eos(id))

            tm(3,i) = s0 + tm(3,i)

         end do

         lct(id) = ilam

      end if

      end

      subroutine incdep (ind)
c-----------------------------------------------------------------------
c either indep or incdp0 are called whenever any primary potential
c variables are changed to reevaluate secondary variables. this
c cumbersome structure is necessitated by the fact that computational
c variables were mapped directly to the thermodynamic variables in
c cst5. a more rational strategy, which is used in werami/pssect, is
c to separate the computational and thermodynamic variables.

c incdep:

c  1) if ind = iind, computes the dependent variable if P(T) or T(P)

c  ind  - is the index of the current variable.
c  iind - is the index of the independent variable that v(idep) is a
c         function of.
c  idep - is the index of the dependent variable.

c  2) if jmct > 0, assigns the independent chemical potentials
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ind

      double precision var

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct
c----------------------------------------------------------------------
      if (ind.eq.iind.and.idep.ne.0) then
         var = v(iind)
         v(idep) = c0 + var*(c1 + var*(c2
     *                          + var*(c3 + var*c4)))
c    *                           var*(c3 + var*(c4 + c5*var))))
      end if

      if (jmct.gt.0) call subinc

      end

      subroutine subinc
c-----------------------------------------------------------------------
c assigns the independent chemical potentials, called by incdep and
c incdp0
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision gref, xp, gcpd

      external gcpd

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision mmu
      common/ cst39 /mmu(i6)
c----------------------------------------------------------------------
      do i = 1, jmct

            if (imaf(i).eq.1) then
c                                 the primary variable is a chemical
c                                 potential.
               mmu(i) = v(3+i)

            else
c                                 the primary variable is a fugacity or
c                                 an activity.
               if (imaf(i).eq.2) then
c                                 fugacity
                  xp = v(1)
                  v(1) = pr
                  gref = gcpd (idaf(i),.false.)
                  v(1) = xp

               else
c                                 activity
                  gref = gcpd (idaf(i),.false.)

               end if

               mmu(i) = gref + r*v(2)*v(3+i)*2.302585093d0

             end if

      end do

      end

      subroutine incdp0
c----------------------------------------------------------------------
c incdep 1) conditionally computes the dependent variable if one exists
c (idep ne 0), 2) assigns independent chemical potentials.

c  idep - is the index of the dependent variable
c  iind - is the index of the independent variable that v(idep) is a
c         function of.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision var

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct
c----------------------------------------------------------------------
      if (idep.ne.0) then
         var = v(iind)
         v(idep) = c0 + var*(c1 + var*(c2
     *                          + var*(c3 + var*c4)))
c    *                            var*(c3 + var*(c4 + c5*var))))
      end if

      if (jmct.gt.0) call subinc

      end


      integer function match (idim,ier,name)
c----------------------------------------------------------------------
c find the index of endmember identified by name
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer idim, ier

      character*8 name

      character mname*8
      common/ cst18a /mname(m4)

      ier = 0

      do match = 1, idim
         if (name.eq.mname(match)) exit
      end do

      if (match.gt.idim) ier = 1

      end

      subroutine readn (i,idim,tname)
c----------------------------------------------------------------------
c readn - read idim endmember names expected format into mname(i+1..i+idim)

c  name_1 name_2 ... | comment
c  name_n+1....      | comment

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, ier, idim, ict, i

      character name*8, tname*(*)

      character mname*8
      common/ cst18a /mname(m4)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ier = 0

      call readcd (n9,ier,.false.)
      if (ier.ne.0) goto 90

      ibeg = 1
      ict = i

      do while (ict-i.lt.idim)
c                                 find the name
         call readnm (ibeg,iend,com,ier,name)
         if (ier.ne.0) goto 90
         ict = ict + 1
         if (ict.gt.m4) call error (1,0d0,ict,
     *                 'm4 (maximum number of endmembers)')

         mname(ict) = name

         if (ibeg.ge.com.and.ict-i.lt.idim) then
            call readcd (n9,ier,.false.)
            ibeg = 1
            if (ier.ne.0) goto 90
         end if

      end do

      return

90    write (*,1000) tname,chars(1:com),name

      call errpau

1000  format ('**error ver200** READN bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a,/,
     *        'last name read was: ',a,/)

      end


      subroutine readda (rnums,idim,tname)
c----------------------------------------------------------------------
c readda - reads idim numbers, discarding all data in records
c          beyond the comment marker "|" or the end of the record.
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      integer idim, kdim, jdim, i, isnum, ier

      character tname*10, nums*(lchar)

      double precision rnums(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
c                                 read card scans for non blank data
c                                 card:
      kdim = 1
      jdim = 0
      isnum = 0
      com = 0
      ier = 1

      do while (jdim.lt.idim)

         call readcd (n9,ier,.true.)
         if (ier.ne.0) exit
c                                 got data, count how many "numbers"
         do i = 1, com
            if (chars(i).ne.' '.and.isnum.eq.0) then
               jdim = jdim + 1
               isnum = 1
            else if (chars(i).eq.' ') then
               isnum = 0
            end if
         end do
c                                 if there is no comment
c                                 marker there can be more
c                                 data than anticipated:
         if (jdim.gt.idim) jdim = idim
c                                 write numbers to string
         write (nums,*) chars(1:com),' '

         read (nums,*,iostat=ier) (rnums(i), i = kdim, jdim)
         if (ier.ne.0) exit

         kdim = jdim + 1

      end do

      if (ier.gt.0) then

         write (*,1000) tname, chars(1:com)
         write (*,1020)

         call errpau

      else if (ier.lt.0) then

         write (*,1010) tname
         write (*,1020)

         call errpau

      end if

1000  format ('**error ver209** READDA bad data, currently',
     *        ' reading solution model: ',/,a,/,'data was:',/,400a)
1010  format ('**error ver210** READDA read to end of file',
     *        ' reading solution model: ',/,a)
1020  format ('READDA was expecting numeric data.',/)

      end

      subroutine readx (idim,tname)
c----------------------------------------------------------------------
c readx - read excess function for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        W(name-name-...) number number number

c end_of_data is either a "|" or the end of the record.

c iord - reinstated as the max order of an excess term, JADC, 10/25/2021.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, jend, ier, iscan, lord, imax, match, idim,
     *        i, iscnlt

      double precision nums(m3)

      character name*8, begin*5, eod*3, tname*10, values*80, key*22

      logical ok

      external iscnlt, iscan

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)
c----------------------------------------------------------------------
      iterm = 0
      iord = 0
      xtyp = 0

      call readcd (n9,ier,.true.)

      write (begin,'(5a)') chars(1:5)

      if (begin.eq.'ideal') then
         return
      else if (begin.ne.'begin') then
         goto 90
      end if

      isub(1:m1,1:m2) = 0

      eod = ' '

      do while (eod.ne.'end')

         call readcd (n9,ier,.true.)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') chars(1:3)
c                                 find expansion type
         if (chars(2).eq.'k'.or.chars(2).eq.'K') then
            xtyp = 1
         else
            xtyp = 0
         end if
c                                 find brackets
         ibeg = iscan (1,com,'(') + 1
         imax = iscan (1,com,')') - 1

         if (ibeg.gt.com.or.imax.gt.com) cycle
c                                 data found
         iterm = iterm + 1
         if (iterm.gt.m1) call error (48,wg(1,1),m1,tname)

         lord = 0

         do while (ibeg.lt.imax)

            call readnm (ibeg,jend,imax,ier,name)
            if (ier.ne.0) goto 90

            lord = lord + 1
            if (lord.gt.m2) call error (49,wg(1,1),m2,tname)

            isub(iterm,lord) = match (idim,ier,name)

            if (ier.ne.0) goto 90

         end do

         if (xtyp.eq.0) then

            ibeg = imax + 2
c                                 read standard form margules pt functions
            rkord(iterm) = lord

            call redlpt (nums,ibeg,jend,ier)

            if (ier.ne.0) goto 90

            do i = 1, m3
               wg(iterm,i) = nums(i)
            end do

         else
c                                 initialize order
            rkord(iterm) = 0
c                                 rk form, read a new card for each term
            wk(1:m16,1:m17,iterm) = 0d0

            do

               ibeg = 1
c                                 check for end of data
               call readcd (n9,ier,.true.)

               write (begin,'(3a)') chars(1:3)

               if (begin.eq.'end') then
                  return
               else if (begin.eq.'Wk(') then
                  exit
               end if
c                                 we have a data card
               rkord(iterm) = rkord(iterm) + 1

               do
c                                 locate end of keyword
                  if (ibeg.ge.com) exit

                  jend = iscan (ibeg,com,'=') - 1
                  if (jend.ge.com) exit
c                                 write keyword
                  write (key,'(22a)',iostat=ier) chars(ibeg:jend)
                  if (ier.ne.0) call error (23,wg(1,1),ier,key)
c                                 locate data
                  ibeg = iscnlt (jend+2,com,' ')
                  jend = iscan (ibeg,com,' ')
c                                 write data
                  write (values,'(80a)',iostat=ier) chars(ibeg:jend)
                  if (ier.ne.0) call error (23,wg(1,1),ier,key)
c                                 shift pointer to next key
                  ibeg = iscnlt(jend,com,' ')
c                                 assign data
                  ok = .false.

                  do i = 1, m16
                     if (key.eq.wstrg(i)) then
                        read (values,*,iostat=ier)
     *                                          wk(i,rkord(iterm),iterm)
                        if (ier.ne.0) call error (23,wg(1,1),ier,key)
                        ok = .true.
                        exit
                     end if
                  end do

                  if (ok) cycle

                  call error (9,wg(1,1),i,key)

               end do

            end do

         end if

         if (rkord(iterm).gt.iord) iord = rkord(iterm)

      end do

      return

90    write (*,1000) tname,chars(1:com)
      write (*,1010) name

      call errpau

1000  format ('**error ver200** READX bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a)
1010  format ('last name read was: ',a,/,
     *        'usually this error is due to a mispelled ',
     *        'endmember name.',/)

      end

      subroutine readop (idim,kstot,tname)
c----------------------------------------------------------------------
c readop - tail of solution model to find optional dqf,
c          van laar size parameters, flagged endmembers, or 
c          reach_increment

c readop - reads data until it finds an     "end_of_model" record

c          van laar data is identified by a "begin_van_la" record
c          dqf data is identified by a      "begin_dqf_co" record
c          endmember flags by a             "begin_flagge" record
c          or the reach factor is           "reach_increm" record


c readop returns:

c          laar if van laar data found.
c          idqf  > 0 if dqf data found (else 0).
c          reach, set = 0 if no reach factor found.
c          and sets endmember flags of indicated endmembers
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, idim, kstot

      character tname*(*), key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c----------------------------------------------------------------------

      idqf = 0
      laar = .false.
      stck = .true.
      norf = .true.
      badx = .false.
      modres = .false.
      unbd = .false.

      do

         call redcd1 (n9,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'end_of_model') then

            exit

         else if (key.eq.'begin_model ') then
c                              found new model, current model
c                              lacks end_of_model keyword
            write (*,1000) tname,chars(1:length)

            call errpau

         else if (key.eq.'begin_van_laar_sizes') then
c                              read van laar data:
            laar = .true.
            call readvl (idim,kstot,tname)

         else if (key.eq.'begin_dqf_corrections') then
c                              read dqf data:
            call readdq (idim,tname)

         else if (key.eq.'reach_increment') then 

c           write (*,*) 'reach_increment obsolete, 6.9.1+ '//tname

         else if (key.eq.'low_reach') then

c           write (*,*) 'low_reach obsolete, 6.9.1+ '//tname

         else if (key.eq.'use_model_resolution') then

            modres = .true.
c           write (*,*) 'set low res for '//tname

         else if (key.eq.'reject_bad_composition') then

            badx = .true.

         else if (key.eq.'begin_flagged_endmembe') then

            call readef (idim,tname)

         else if (key.eq.'site_check_override') then

            stck = .false.

         else if (key.eq.'refine_endmembers') then

            norf = .false.

         else if (key.eq.'unbounded_composition') then 

            unbd = .true.

         else

            write (*,1010) tname,chars(1:length)
            write (*,1020)

            call errpau

         end if

      end do

1000  format (/,'**error ver200** READOP missing "end_of_model"',
     *          ' keyword at end',' of solution model:',a,/)

1010  format (/,'**error ver210** READOP bad data, currently',
     *          ' reading solution model: ',a,' data was:',/,400a)

1020  format (/,'This error is most probably due to an out-of-date ',
     *          'solution model file.',//,
     *          'Copy the current version from:',//,
     *          'www.perplex.ethz.ch/perplex/datafiles/',
     *          'solution_model.dat',//)

      end

      subroutine readvl (idim,kstot,tname)
c----------------------------------------------------------------------
c readvl - read van laar volumes for a solution models endmembers, assumes
c data on one line of less than 240 characters, the expected format

c        alpha(name) number number number

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, ier, iscan, imax, match, idim, index

      external iscan, match

      character name*8, eod*3, tname*10

      double precision nums(m3)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer kstot,jend,i,ict
c----------------------------------------------------------------------

      ict = 0

      eod = ' '

      do while (eod.ne.'end')

         call readcd (n9,ier,.true.)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') chars(1:3)
c                                 find brackets
         ibeg = iscan (1,com,'(') + 1
         imax = iscan (1,com,')') - 1

         if (ibeg.gt.com.or.imax.gt.com) cycle
c                                 data found
         ict = ict + 1
         if (ict.gt.m4) goto 91

         call readnm (ibeg,jend,imax,ier,name)
         if (ier.ne.0) goto 90

         index = match (idim,ier,name)
         if (ier.ne.0) goto 90

         ibeg = imax + 2

         call redlpt (nums,ibeg,jend,ier)

         if (ier.ne.0) goto 90

         do i = 1, m3
            vlaar(i,index) = nums(i)
         end do

      end do

      if (ict.lt.kstot) goto 91

      return

90    write (*,1000) tname,chars(1:com),vlaar(i,index)
      write (*,1001)
      call errpau

91    write (*,1010) tname
      call errpau

1000  format ('**error ver200** READVL bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

1001  format (/,'usually this error is caused by a mispelled ',
     *          'endmember name.',/)

1010  format (' **error ver201** READVL bad data, currently',
     *        ' reading solution model: ',a,/,
     *        ' this model requires 1 size parameter for',
     *        ' each independent endmember, READVL found ',i2,
     *        ' parameters.',/)

      end

      subroutine readdq (idim,tname)
c----------------------------------------------------------------------
c readvl - read dqf corrections for solution models endmembers, assumes
c data on one line of less than 240 characters, the expected format

c        dqf(name) number number number

c        output:

c          idqf           - the number of endmembers with dqf corrections
c          indq(idqf)     - pointer to corrected endmember in the phase name array
c          dqf(1..3,idqf) - the dqf parameters for the correction

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, ier, iscan, imax, match, idim

      external iscan, match

      character name*8, eod*3, tname*10

      double precision nums(m3)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer jend,i,idqf,indq
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c----------------------------------------------------------------------

      eod = ' '

      do while (eod.ne.'end')

         call readcd (n9,ier,.true.)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') chars(1:3)
c                                 find brackets
         ibeg = iscan (1,com,'(') + 1
         imax = iscan (1,com,')') - 1

         if (ibeg.gt.com.or.imax.gt.com) cycle
c                                 data found
         idqf = idqf + 1

         call readnm (ibeg,jend,imax,ier,name)
         if (ier.ne.0) goto 90

         indq(idqf) = match (idim,ier,name)
         if (ier.ne.0) goto 90

         ibeg = imax + 2

         call redlpt (nums,ibeg,jend,ier)
         if (ier.ne.0) goto 90

         do i = 1, m3
            dqf(i,idqf) = nums(i)
         end do

      end do

      return

90    write (*,1000) tname,chars(1:com),dqf(i,idqf)
      write (*,1001)

      call errpau

1000  format ('**error ver200** READDQ bad data, currently',
     *        'reading solution model: ',a,' data was:',/,400a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)
1001  format (/,'usually this error is caused by a mispelled ',
     *          'endmember name.',/)

      end

      subroutine readr (coeffs,enth,inds,idim,nreact,tname,eor)
c----------------------------------------------------------------------
c readr - read stoichiometric reaction data for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        name "=" (acoef(i), mame(i), i= 2..nreact) Enthalpy_value

c nreact >= 4 + nord for reciprocal reactions
c nreact = -1 on input for ordered/disorder reactions
c enthalpy_value is only read if on input nreact = -1

c end_of_data is either a "|" or the end of the record.

c eor - indicates end of data for 688 format solution models
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, ier, iscan, nreact, inds(k7),
     *        idim, match, iscnlt, i

      logical eor

      double precision coeffs(k7), enth(3), rnum

      character name*8, tname*10, tag*3

      external iscan, iscnlt, match

      character mname*8
      common/ cst18a /mname(m4)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ier = 0

      call readcd (n9,ier,.true.)
      if (ier.ne.0) goto 90

      ibeg = 1

      write (tag,'(3a)') chars(1:3)

      if (tag.eq.'end') then 
         eor = .true.
         return
      else 
         eor = .false.
      end if 
c                                 first name
      call readnm (ibeg,iend,com,ier,name)

      if (ier.ne.0) goto 90

      if (nreact.eq.-1) then
c                                 if nreact = -1, new name
         idim = idim + 1
         mname(idim) = name
         inds(1) = idim

      else
c                                 reciprocal sol, get index
         inds(1) = match(idim,ier,name)

         if (ier.ne.0) then
            rnum = 1d0
            goto 90
         end if

      end if
c                                 find marker '='
      ibeg = iscan (1,com,'=') + 1

      i = 2

      do
c                                 find a stoichiometric coeff
         call readfr (rnum,ibeg,iend,com,ier)
         if (ier.ne.0) exit

         coeffs(i) = rnum
c                                 find the name
         call readnm (ibeg,iend,com,ier,name)

         if (ier.ne.0) goto 90

         if (i.gt.k7) call error (1,0d0,i,'k7')

         inds(i) = match(idim,ier,name)
         if (ier.ne.0) goto 90

         if (nreact.gt.0.and.i.eq.nreact) exit
c                                 increment counter and find
c                                 the next coeff+name
         i = i + 1

      end do

      if (nreact.eq.-1) then
c                                 ordered compound, read
c                                 enthalpy, find marker '='
         ibeg = iscan (ibeg,com,'=') + 2
         call redlpt (enth,ibeg,iend,ier)

         nreact = i - 2

         if (ier.ne.0) goto 90

      else if (i.lt.3) then
c                                 a reaction with only 2 species
c                                 is unexpected, write error
         goto 90

      else

         nreact = i - 1

      end if

      return

90    write (*,1000) tname,chars(1:com),name,rnum

      call errpau

1000  format ('**error ver200** READR bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a,
     *        'last name read was: ',a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

      end

      subroutine readz (coeffs,inds,ict,idim,tname,tag)
c----------------------------------------------------------------------
c readz - read site fraction data for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        comments "=" a0 (acoef(i), mame(i), i= 1...end_of_data)

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, ier, iscan, ict, inds(k7),
     *        match, idim, i, iscnlt, jend, jbeg

      external iscan, iscnlt

      double precision rnum, coeffs(k7)

      character name*8, tname*10, tag*3

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ict = 0
      do i = 1, k7
         inds(i) = 0
         coeffs(i) = 0
      end do

      call readcd (n9,ier,.true.)
      if (ier.ne.0) goto 90
c                                 this first segment is only done for
c                                 the readlm routine:
      ibeg = 1
c                                 read the fist word
      call readnm (ibeg,iend,com,ier,name)
      read (name,'(a)') tag

      if (tag.eq.'end') then
c                                 if called from readlm tag may be the
c                                 end of data marker
         return

      else

         i = match(idim,ier,name)

         if (ier.eq.0) then
c                                 the "comment" is a name for
c                                 readlm
            ict = ict + 1
            inds(ict) = i

         end if

      end if
c                                 extract the species name for 688 models and return it as tag
      jbeg = iscan(1,iend,'(') 
      jend = iscan(1,iend,',')
c                                 if no comma, find right )
      if (jend.gt.iend) jend = iscan(1,iend,')')
      if (jend-jbeg.gt.4) jend = jbeg + 4
      write (tag,'(3a)') chars(jbeg+1:jend-1)
c                                 find start of data marker '='
      ibeg = iscan (iend,com,'=') + 1
      ict = ibeg
c                                 find a number
      call readfr (rnum,ibeg,iend,com,ier)
      if (ier.ne.0) goto 90
c                                 find the next non-blank character
c                                 if it's text, then the expression
c                                 has no constant, else save the
c                                 constant.
      if (chars(iscnlt(iend+1,lchar,'/')).lt.'A') then
c                                 assuming ascii collating sequence,
c                                 the next character is numeric, so
c                                 the number read previously is the a0 constant
         coeffs(1) = rnum

      else
c                                 no constant, reset character pointer.
         coeffs(1) = 0d0
         ibeg = ict

      end if
c                                 the rest of the data should
c                                 consist of coefficients followed
c                                 by names
      ict = 1

      do while (ibeg.lt.com)
c                                 find the number
         call readfr (rnum,ibeg,iend,com,ier)

         if (ier.ne.0) then
c                                 if called from readlm may be the
c                                 legitimate text string "delta"
            call readnm (ibeg,iend,com,ier,name)

            if (name.eq.'delta') then

               ibeg = iscan (iend,com,'=') + 1
               call readfr (rnum,ibeg,iend,com,ier)
               if (ier.ne.0) goto 90
               coeffs(ict+1) = rnum
               exit

            else
c                                 invalid data
               goto 90

            end if

         end if

         if (ier.ne.0) goto 90
c                                 find the name
         call readnm (ibeg,iend,com,ier,name)
c                                 for constant bounds, read delta
         if (name.eq.'delta') then
c                                 the lower bound should be the last number read
            coeffs(ict) = rnum
c                                 next find the delta
            ibeg = iscan (iend,com,'=') + 1
            call readfr (rnum,ibeg,iend,com,ier)
            if (ier.ne.0) goto 90
            coeffs(ict+1) = rnum

            exit

         end if


         if (ier.ne.0) goto 90

         ict = ict + 1
         coeffs(ict) = rnum
         inds(ict) = match(idim,ier,name)

         if (ier.ne.0) then
            write (*,1010) name,tname,chars(1:com)

            call errpau

         end if

      end do

      return

90    write (*,1000) tname,chars(1:com),name,rnum

      call errpau

1010  format (/,'**error ver201** invalid name: ',a,' in an expression',
     *        ' for solution model: ',a,/,' data was:',/,400a)
1000  format (/,'**error ver200** READZ bad data, currently',
     *        ' reading solution model: ',a,' data was:',/,400a,/,
     *        'last name read was: ',a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

      end

      double precision function gsixtr (id)
c-----------------------------------------------------------------------
c gsixtr computes G from the EoS formulated by Sixtrude & Bukowski '90,
c JGR 95:19311-19325.

c Sixtrude parameters:

c    F0    n  -v(J/bar) K0(bar) K0'   theta0 gamma0 q etaS0 Smag

c Perple_X parameters (1st index in array thermo):

c    1     2     3        4     5      6      7    8    9    10

c and for shear modulii

c    G0(bar)   G0'

c are in the Perple_X array emod(i,id) with i

c      1        2
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, itic, izap

      double precision a0, v0, v, df, f, dv, root, nr9t0, nr9t,
     *                 gamma0, k00, plg, c1, c2, c3, f1, gvq,
     *                 q, vq, v23, theta0, tol, k0p, a, ethv

      double precision nr9, qm1, d2f, tht, tht0, etht, etht0, ltht, df1,
     *                 ltht0, dtht, dtht0, d2tht, d2tht0, g, g0, dg,
     *                 dg0, d2g, d2g0, dfc, d2fc, dft, d2ft, dft0, d2ft0

      double precision smu
      common/ cst323 /smu

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save izap
      data izap /0/
c----------------------------------------------------------------------
c                                 assign local variables:
      a0     =  thermo(1,id)
      v0     = -thermo(3,id)

      theta0 =  thermo(6,id)
      gamma0 =  thermo(7,id)
      q      =  thermo(8,id)

      nr9    = thermo(11,id)
      c1     = thermo(12,id)
      c2     = thermo(13,id)
      c3     = thermo(14,id)
      qm1    = q - 1d0
      nr9t   = nr9*t
      nr9t0  = thermo(20,id)
c                                 initial guess for volume:
c                                 taylor(diff(FTH,v),v=v0,1)
c                                 JADC Feb 26, 2008.
c                                 the dfth0 could be loaded as a
c                                 constant.
      tht    = theta0/t
      tht0   = theta0/tr
      k00    = thermo(4,id)
      k0p    = thermo(5,id)

      dft   = nr9t*gamma0/v0*(3d0*plg(tht)/tht**3
     *        - dlog(1d0-exp(-tht)))
      dft0  = nr9t0*gamma0/v0*(3d0*plg(tht0)/tht0**3
     *        - dlog(1d0-exp(-tht0)))
c                                 taylor(diff(FC,v),v=v0,2)
c     v       = (k00-dft+dft0-p)/k00*v0
c                                 taylor(diff(FC,v),v=v0,3)
      root = k00*((2d0+2d0*k0p)*(p+dft-dft0) + k00)

      if (root.gt.0d0) then
         v = ((2d0+k0p)-dsqrt(root)/k00)*v0/(1d0+k0p)
         if (v.lt.v0/1d1.or.v.gt.v0*1d1) v = v0
      else
         v = v0
      end if

      f1 = 1d9
      itic = 0
c                                 change to use relative tolerance
c                                 JADC March 1, 2005. formerly 1d-1 bar.
      tol = 1d-6*p

      do while (dabs(f1).gt.tol)

         itic = itic + 1

         vq = (v/v0)**q
         gvq = gamma0*vq
         v23 = (v0/v)**r23
         f = 0.5d0*v23 - 0.5d0
         df = -v23/v/3d0
         d2f = 5d0/9d0*v23/v**2

         tht  =  theta0*dexp(-gamma0*((v/v0)**q-1d0)/q)/t
         if (tht.lt.1d-10) goto 90
         tht0 =  tht*t/tr
         etht  = dexp(-tht )
         etht0 = dexp(-tht0)
         ltht  = dlog(1d0 - etht )
         ltht0 = dlog(1d0 - etht0)
c                                 diff(theta/T,v)
         dtht  = -gvq/v*tht
         dtht0 = -gvq/v*tht0
c                                 diff(theta/T,v,v)
         d2tht  = gvq*tht /v**2*(gvq - qm1)
         d2tht0 = gvq*tht0/v**2*(gvq - qm1)

         g   = plg(tht )
         g0  = plg(tht0)
         dg  = tht **2*ltht *dtht
         dg0 = tht0**2*ltht0*dtht0
         d2g  = ((2d0*ltht  + tht *etht /(1d0-etht ))*dtht **2 +
     *          tht *ltht *d2tht )*tht
         d2g0 = ((2d0*ltht0 + tht0*etht0/(1d0-etht0))*dtht0**2 +
     *          tht0*ltht0*d2tht0)*tht0

         dfc = (c3*f+c1)*f*df
         d2fc = (2d0*c3*f+c1)*df**2+(c3*f+c1)*f*d2f

         dft  = nr9t /tht **3*(dg  -3d0/tht *g *dtht )
         dft0 = nr9t0/tht0**3*(dg0 -3d0/tht0*g0*dtht0)

         d2ft =  nr9t /tht **3*(3d0/tht *(dtht *(4d0/tht *g *dtht
     *                        - 2d0*dg ) - g *d2tht ) + d2g )
         d2ft0 = nr9t0/tht0**3*(3d0/tht0*(dtht0*(4d0/tht0*g0*dtht0
     *                        - 2d0*dg0) - g0*d2tht0) + d2g0)

         f1  = -dfc - dft + dft0 - p

         df1 = -d2fc - d2ft + d2ft0

         dv = f1/df1
         v = v - dv

         if (v.le.0d0.or.v/v0.gt.2d1.or.
     *       itic.gt.100.or.dabs(f1).gt.1d40) goto 90

      end do
c                                 if everything is ok, now get
c                                 helmoltz energy:
      goto 10
c                                 if we get here, failed to converge
90    if (izap.lt.10) then
         write (*,1000) t,p,names(id)
         izap = izap + 1
         if (izap.eq.10) call warn (49,r,369,'GETLOC')
      end if
c                                 destabilize the phase:
      gsixtr = 1d2*p

      return

10    vq = (v/v0)**q
      f = 0.5d0*(v0/v)**r23 - 0.5d0
      tht  =  theta0*dexp(-gamma0*(vq-1d0)/q)/t
      tht0 =  tht*t/tr

      a = a0 + c1*f**2*(0.5d0 + c2*f)
     *  + nr9*(t/tht**3*plg(tht ) -tr/tht0**3*plg(tht0))

      gsixtr = a + p*v - t*thermo(10,id)
c                                 thermal energy/v
      ethv = (dft0-dft)/gamma0/vq
c                                 etas0 = thermo(9,id)
c                                 g0 = emod(1,id)
c                                 g0p = emod(2,id)
c                                 etas = thermo(9,id)*v/v0
c                                 adiabatic shear modulus
      smu = (1d0 + 2d0*f)**(2.5d0)*
     *      (emod(1,id)*(1d0 - 5d0*f) + f*emod(2,id)*3d0*k00)
     *    -  thermo(9,id)*v/v0*ethv

1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude EoS.',
     *        ' Phase ',a,' will be destabilized.',/)

      end

      double precision function ghkf (id)
c-----------------------------------------------------------------------
c ghkf computes apparent G for aqueous species HKF formulation

c assumes perimittivity (epsln) and HKF g-function (gf) have been computed
c in common cxt37 (by slvnt2 or aqrxdo).

c HKF parameters are loaded into thermo as:

c thermo(1 ,id) = G0
c thermo(2 ,id) = S0
c thermo(5 ,id) = w (omega0)
c thermo(6 ,id) = q (charge)
c thermo(7 ,id) = a1
c thermo(8 ,id) = a2
c thermo(9 ,id) = a3
c thermo(10,id) = a4
c thermo(11,id) = c1
c thermo(12,id) = c2
c thermo(13,id) = -s + c1*dlog(tr) + c1 + w*yr + dlog(tr/(tr-theta))*c2/theta**2 => b8 in HKF_G.mws
c thermo(14,id) = (-w*yr-c1+s)*tr + w - a1*pr - a2*ln(psi+pr) + g + c2/theta => b9
c thermo(15,id) = -a3*pr-a4*ln(psi+pr) => b10
c thermo(16,id) = -c2/(tr-theta)/theta => b11
c thermo(17,id) = c2/theta^2 => b12
c thermo(18,id) = -c1-c2/theta^2
c thermo(19,id) = reference born radius 5d10*eta*q^2/(1622323167*eta*q+5d10*omega0)
c thermo(20,id) = q^2
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision ft, fp, omega, psi, theta, z, eta

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      integer iam
      common/ cst4 /iam

      save psi, theta, eta
      data psi, theta, eta/2600d0, 228d0, 694656.968d0/
c-----------------------------------------------------------------------
      if (id.eq.aqst+ihy) then
c                                 assumes proton is the only species
c                                 with zero G0, return G_H+(P,T) = 0.
         ghkf = 0d0
         return

      else if (iam.eq.5) then
c                                 needs solvent properties if frendly
        call slvnt0 (fp,ft)

      end if

      z = thermo(6,id)

      if (z.ne.0d0) then
c                                 ionic species
         omega = eta * z * (z/(thermo(19,id) + dabs(z)*gf)
     *                      - 1d0/(3.082d0 + gf))

      else
c                                 neutral species
         omega = thermo(5,id)

      end if

      ft = t - theta
      fp = dlog(psi+p)

      ghkf = thermo(14,id) + (thermo(13,id) + thermo(17,id)*dlog(ft)
     *                                      + thermo(18,id)*dlog(t))*t
     *     + thermo(16,id)*ft
     *     + thermo(7,id)*p + thermo(8,id)*fp
     *     + (thermo(9,id)*p + thermo(10,id)*fp + thermo(15,id))/ft
     *     + omega*(1d0/epsln - 1d0) - thermo(5,id)/epsln0

      end

      double precision function gfunc (rho)
c----------------------------------------------------------------------
c Shock et al 1992 dielectic g function [HKF_g_function.mws].
c g is in angstrom. rho is the CGS solvent density.

c hacked to use solvent density, this will, at a minimum violate the
c conditions based on p/t
c----------------------------------------------------------------------
      implicit none

      integer iwarn

      double precision g, tf, psat2, rho

      external psat2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      logical abort1
      common/ cstabo /abort1

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      save iwarn
      data iwarn/0/
c---------------------------------------------------------------------
      abort1 = .false.

      if (rho.gt.1d0) then
c                                 region III, rho = 1 g/cm3, g = 0
         g = 0d0

      else
c                                 region I function
         g = ((-6.557892d-6*t + 9.3295764d-3)*t
     *       -4.096745422)*((1d0 - rho)) **
     *       ((1.268348e-5*t - 1.767275512e-2)*t + 9.98834792)

         if (t.gt.428.15.and.p.lt.1d3) then
c
            tf = (t/300d0 - 1.427166667d0)
c                                 add region II perturbation term
            g  = g -
     *           (tf**4.8d0 + 0.366666D-15*tf**16)
     *         * ((((5.01799d-14*p - 5.0224d-11)*p - 1.504074d-7)*p
     *               + 2.507672d-4)*p - 0.1003157d0)

         end if
c                                 check on physical conditions
         if (rho.lt.0.35d0.or.
     *       (t.gt.623.15.and.p.lt.500d0).or.
     *       (t.le.623.15.and.p.lt.psat2(t))) then
c                                 warn

            if (iwarn.lt.10) then
               write (*,1000) t, p
               if (ns.eq.1) write (*,'(/,a,/)') 
     *                      'No result will be output.'

               iwarn = iwarn + 1
               if (iwarn.eq.10) call warn (49,r,277,'GFUNC')
            end if

            if (ns.eq.1) abort1 = .true.

            g = 0d0

         end if

      end if

      gfunc = g

1000  format (/,'**warning ver277** T= ',f8.2,' K P=',f9.1,' bar ',
     *       'is beyond the limits for',/,'the HKF g function. The ',
     *       'function will be zeroed.',/)

      end


      double precision function epsh2o (v)
c----------------------------------------------------------------------
c Sverjenski 2014 dielectic constant for pure water, v is the molar
c volume of water in j/bar [HKF_epsilon.mws].
c----------------------------------------------------------------------
      implicit none

      double precision v

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      epsh2o = dexp(-0.8016651D-4 * t + 0.4769870482D1 - 0.6871618D-1 *
     *         dsqrt(t - 0.27315D3)) * (0.1801526833D1 / v) **
     *         (-0.1576377D-2 * t + 0.1185462878D1 + 0.6810288D-1 *
     *         dsqrt(t - 0.27315D3))

      end

      double precision function gaq (id)
c-----------------------------------------------------------------------
c gaq computes apparent G for aqueous species with the Anderson et al.
c (GCA 1991) density model as extended by Holland & Powell (JMG 1998)
c and modified to account for P-Pr integration.

c parameters, 0 indicates property at Tr,Pr:

c t10 := -s0+Tr*b0+(-Tr*b0+c0)/Tr/dadt0*alpha0;
c t11 := (-Tr*b0+c0)/Tr/dadt0;
c t12 := -1/2*b0;
c t13 := Tr*s0+g0-pr*v0-1/2*b0*Tr^2+(-Tr*b0+c0)/Tr/dadt0*(-Tr*alpha0+beta0*pr);
c t14 := v0-(-Tr*b0+c0)/Tr/dadt0*beta0;

c vh2o = water volume
c vh2o0 = water volume
c alpha0 = water expansivity
c beta0 = water compressibility
c dadT0 = temperature derivative of alpha
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision vh2o, vh2o0, fh2o, tp

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save vh2o0
      data vh2o0/18.7231148995863/
c----------------------------------------------------------------------
      if (thermo(21,id).eq.1d0) then
c                                 thermo(21) flags hydronium. cause hp
c                                 use HSC convention, gaq(H+) is not zero.
         gaq = thermo(1,id)
         return

      end if
c                                 compare to last p-t, to save expensive
c                                 calls to the h2o eos.
      call pseos (vh2o,fh2o,1)
c                                 tp instead of t is a h-p innovation
      if (t.lt.500d0) then
         tp = t
      else
         tp = 500d0
      end if

      gaq = thermo(13,id)
     *      + t*(thermo(10,id) + thermo(11,id)*dlog(vh2o0/vh2o)/tp
     *         + thermo(12,id)*t) + thermo(14,id)*p

      end

      double precision function gstxlq (id)
c-----------------------------------------------------------------------
c gstxlq computes G from the liquid EoS formulated by Sixtrude et al 2009

c Stxrude parameters:

c   thermo(1)  = F0
c   thermo(2)  = S0 - Cv
c   thermo(3)  = V0
c   thermo(4)  = Cv
c   thermo(5)  = 4.5d0*K0*V0
c   thermo(6)  = 4.5d0*K0*V0*(K'-4)
c   thermo(7)  = y0 - y'
c   thermo(8)  = y'
c   thermo(9)  = T0
c   --- dependent ---
c   thermo(10) = (S0-Cv-Cv*y0)*T0
c   thermo(11) = Cv*(y0+ln(T0))-S0+Cv
c   thermo(12) = ln(v0)

c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, itic, izap

      logical bad

      double precision a5, a6, a7, a8, a9, a10, a11, v0, v, v2, df, f,
     *                 d2a, v23, tol, d2f, df2, da, a1

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save izap
      data izap /0/
c----------------------------------------------------------------------
c                                 assign local variables:
      v0 = thermo(3,id)

      a10 = thermo(4,id)*(thermo(9,id)-t)*thermo(7,id)
      a7  = (thermo(11,id)-thermo(4,id)*dlog(t))*t+thermo(10,id)
     *      - a10*thermo(12,id)
      a8  = thermo(5,id)
      a9  = thermo(6,id)
      a11 = thermo(4,id)*(thermo(9,id)-t)*thermo(8,id)/v0
      a6  = 3d0*a9
      a5  = 2d0*a8
c                                 initial guess for volume, taylor(diff(a,v),v=v0,3)
      a1 = (p+a11)*v0
      v = v0 + (9d0*(3d0*a8+a9)/(a5+a1*9d0)**2*(a1+a10) - 1d0)
     *         *9d0*v0*(a10+a1)/(a5+a1*9d0)

      if (v.lt.v0/1d1.or.v.gt.v0*1d1) v = v0

      itic = 0
      tol = 1d-6*p

      do

         itic = itic + 1
c                                 f, and derivatives
         v23 = (v0/v)**r23
         v2  = v**2
         f   = 0.5d0*v23 - 0.5d0
         df  = -v23/v/3d0
         df2 = df*df
         d2f = r59*v23/v2
c                                 a 1st and 2nd derivatives
c                                 da is actually diff(a,v) + p
         da  = (a5 + a6*f)*f*df + a10/v + a11 + p
         d2a = (df2 + f*d2f)*a5 + (2d0*df2 + f*d2f)*a6*f - a10/v2

         v = v - da/d2a

         if (v.le.0d0.or.itic.gt.100.or.dabs(da).gt.1d40) then
            bad = .true.
            exit
         else if (dabs(da).lt.tol) then
            bad = .false.
            exit
         end if

      end do

      if (bad) then
c                                 if we get here, failed to converge
         if (izap.lt.10) then
            write (*,1000) t,p,names(id)
            izap = izap + 1
            if (izap.eq.10) call warn (49,r,369,'GSTXLQ')
         end if
c                                 destabilize the phase.
         gstxlq  = 1d2*p

      else
c                                 everything ok, final f:
         f = 0.5d0*(v0/v)**r23 - 0.5d0
c                                 g = helmholtz enery + pv
         gstxlq  = (a8 + a9*f)*f**2 + a7 + a10*dlog(v) + a10 + a11*v
     *             + p*v + thermo(1,id)

      end if

1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude Liq EoS.',
     *        ' Phase ',a,' will be destabilized.',/)

      end

      double precision function gstxgi (id)
c-----------------------------------------------------------------------
c gstxgi computes G from the EoS formulated by Sixtrude & Lithgow-B

c Stxrude parameters:

c    F0    -n  -v(J/bar) K0(bar) K0'   theta0 gamma0 q etaS0 Smag

c Perple_X parameters (1st index in array thermo):

c    1     2     3        4     5      6      7    8    9    10

c and for shear modulii

c    G0(bar)   G0'

c are in the Perple_X array emod(i,id) with i

c      1        2
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, itic, izap

      logical bad

      double precision v0, v, df, f, dv, gamma0, k00, k0p,
     *           plg, c1, c2, c3, f1, aiikk, aiikk2, nr9t,
     *           root, aii, etas, a, ethv, gamma, da, nr9t0,
     *           fpoly, fpoly0, letht, letht0, z, aii2,
     *           v23, t1, t2, a2f

      double precision nr9, d2f, tht, tht0, etht, etht0, df1,
     *                 dtht, dtht0, d2tht, d2tht0,
     *                 dfc, d2fc, dfth, d2fth, dfth0, d2fth0

      double precision smu
      common/ cst323 /smu

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save izap
      data izap /0/
c----------------------------------------------------------------------
c     call begtim(2)
c                                 assign local variables:
      v0     = -thermo(3,id)
      nr9    = thermo(11,id)
      c1     = thermo(12,id)
      c2     = thermo(13,id)
      c3     = thermo(14,id)
      aii    = thermo(15,id)
      aiikk  = thermo(16,id)
      aiikk2 = thermo(18,id)
      aii2   = thermo(19,id)
      nr9t0  = thermo(20,id)

      t1     = thermo(6,id)/t
      t2     = t/tr
      nr9t   = nr9*t
c                                 initial guess for volume:
c                                 taylor(diff(FTH,v),v=v0,1)
c                                 JADC Feb 26, 2008.
c                                 the dfth0 could be loaded as a
c                                 constant.
      tht    = t1
      tht0   = tht*t2
      gamma0 = thermo(7,id)
      k00    = thermo(4,id)
      k0p    = thermo(5,id)

      dfth   = nr9t*gamma0/v0*(3d0*plg(tht)/tht**3
     *         - dlog(1d0-exp(-tht)))
      dfth0  = nr9t0*gamma0/v0*(3d0*plg(tht0)/tht0**3
     *         - dlog(1d0-exp(-tht0)))
c                                 taylor(diff(FC,v),v=v0,2)
c     v       = (k00-dfth+dfth0-p)/k00*v0
c                                 taylor(diff(FC,v),v=v0,3)
      root = k00*((2d0+2d0*k0p)*(p+dfth-dfth0) + k00)

      if (root.gt.0d0) then
         v = ((2d0+k0p)-dsqrt(root)/k00)*v0/(1d0+k0p)
         if (v.lt.v0/1d1.or.v.gt.v0*1d1) v = v0
      else
         v = v0
      end if

      itic = 0

      do

         itic = itic + 1
c                                 f, and derivatives
         v23 = (v0/v)**r23
         f = 0.5d0*v23 - 0.5d0
         df = -v23/v/3d0
         d2f = r59*v23/v**2
c                                 cold part derivatives
         dfc = (c3*f+c1)*f*df
         d2fc = (2d0*c3*f+c1)*df**2+(c3*f+c1)*f*d2f
c                                 debye T/T (tht)
         z  = 1d0+(aii+aiikk2*f)*f

         if (z.lt.0d0) then
            bad = .true.
            exit
         end if

         root = dsqrt(z)

         tht   = t1*root
         tht0  =  tht*t/tr
c                                 tht derivatives
         a2f   = aii2+aiikk2*f
         da    = a2f/root
         dtht  = t1*da*df
         d2tht = t1*((aiikk2/root-a2f**2/z**1.5d0)*df**2
     *               + da*d2f)

         dtht0 = dtht*t2
         d2tht0 = d2tht*t2
c                                 polylog functions:
         fpoly   = 3d0*plg(tht )/tht**3
         fpoly0  = 3d0*plg(tht0)/tht0**3
c                                 thermal part derivatives:
         etht  = dexp(-tht )

         if (1d0-etht.lt.0d0) then
            bad = .true.
            exit
         end if

         letht = dlog(1d0-etht)

         dfth = (letht-fpoly)*nr9t*dtht/tht
         d2fth = ((4d0*dtht**2/tht-d2tht)*(fpoly-letht)
     *         + dtht**2*etht/(1d0-etht))*nr9t/tht

         etht0 = dexp(-tht0)

         if (1d0-etht0.lt.0d0) then
            bad = .true.
            exit
         end if

         letht0 = dlog(1d0-etht0)

         dfth0 = (letht0-fpoly0)*nr9t0*dtht0/tht0
         d2fth0 = ((4d0*dtht0**2/tht0-d2tht0)*(fpoly0-letht0)
     *          + dtht0**2*etht0/(1d0-etht0))*nr9t0/tht0

         f1  = -dfc - dfth + dfth0 - p

         df1 = -d2fc - d2fth + d2fth0

         dv = f1/df1

         if (v - dv.lt.0d0) dv = v/2d0

         v = v - dv

         if (itic.gt.iopt(21).or.dabs(f1).gt.1d40) then
            bad = .true.
            exit
         else if (dabs(dv/(1d0+v)).lt.nopt(50)) then
            bad = .false.
            exit
         end if

      end do

      if (bad) then
c                                 if we get here, failed to converge
         if (izap.lt.10) then
            write (*,1000) t,p,names(id)
            izap = izap + 1
            if (izap.eq.10) call warn (49,r,369,'GSTX')
         end if
c                                 destabilize the phase.
         gstxgi  = 1d2*p

      else

c                                 everything is ok, now get
c                                 helmoltz energy:
         f = 0.5d0*(v0/v)**r23 - 0.5d0
         z = 1d0+(aii+aiikk2*f)*f
         root = dsqrt(z)
c                                 final estimate for tht
         tht   = t1*root
         tht0  = tht*t2
c                                 helmholtz enery
         a = thermo(1,id) + c1*f**2*(0.5d0 + c2*f)
     *     + nr9*(t/tht**3*plg(tht ) -tr/tht0**3*plg(tht0))

         gstxgi = a + p*v - t*thermo(10,id)
c                                 z = (theta/theta0)^2
         gamma = (2d0*f+1d0)*(aii+aiikk*f)/6d0/z
         etas = - gamma - thermo(17,id)/z*(2d0*f + 1d0)**2
c                                 thermal energy/V, based on
c                                 previous v estimate
         if (gamma.ne.0d0) then
            ethv = (dfth0-dfth)/gamma
         else
            ethv = 0d0
         end if
c                                 adiabatic shear modulus
         smu = (1d0 + 2d0*f)**(2.5d0)*(
     *         emod(1,id) + f*(thermo(21,id) + thermo(22,id)*f))
     *       - etas*ethv

      end if

1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude GI EoS.',
     *        ' Phase ',a,' will be destabilized.',/)

c     call endtim (1,.false.,'tot')
c     call endtim (2,.false.,'stx')
c     call begtim (1)

      end

      double precision function plg (t)
c-----------------------------------------------------------------------
c evaluates debye integral: int((ln(1-exp(-t))*t^2),t=0..t)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision t, p0, p1, p2, p3, p4, dinc
c-----------------------------------------------------------------------
c     call begtim(3)
      p0 = dexp(-t)
      p1 = 1d0
      p2 = t*t
      p3 = 2d0*t

      plg = -2.1646464674223d0

      do i = 1, 100000

         p4 = dfloat(i)

         p1 = p1 * p0

         dinc = (p2 + (p3 + 2d0/p4)/p4)*p1/p4/p4

         plg = plg + dinc

         if (dabs(dinc/(1d0+dabs(plg))).lt.nopt(50)) exit

      end do

c     call endtim (3,.false.,'plg')
      end

      double precision function vdpbm3 (vt,k,kprime)
c-----------------------------------------------------------------------
c vdpbm3 computes the vdp integral of a compound identified by id
c that is described by Birch-Murnaghan 3rd order EoS.
c    vt - is the volume at Pr & T
c    k  - is the bulk modulus at Pr & T
c    kprime - is -K' at Pr and T
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer itic, jerk

      double precision k, vt, rat, rat2, c0, c1, c2,
     *                 c3, c4, c5, a0, a1, v, df, f, dv, kprime

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save jerk
      data jerk /0/
c----------------------------------------------------------------------
c                                 constants:
      a0 = 0.375d0 * vt * k
      a1 = -0.125d0 * vt**2 * k
      c0 = (-28d0 -6d0 * kprime) * vt * a0
      c1 = (12d0 + 3d0 * kprime) * vt**2 * a0
      c2 = (16d0 + 3d0 * kprime) * a0
      c3 = a1 * vt * (-196d0 - 42d0 * kprime)
      c4 = a1 * (80d0 + 15d0 * kprime)
      c5 = a1 * vt * (108d0 + 27d0 * kprime)
c                                 use murnaghan guess for volume. GH, 6/23/16
c                                 initial guess for volume:
      dv = 1d0
      v = vt * (1d0 - kprime*p/k)**(dv/kprime)
      itic = 0

      do while (dabs(dv/(1d0+v)).gt.nopt(50))

         itic = itic + 1
         rat = (vt/v)**r13
         rat2 = rat**2
         f = p  + ((c0*v*rat+c1+c2*v**2*rat2)/v**3)
         df = (c3/rat2+c4*v/rat+c5)/v**4
         dv = f/df
         v = v - dv

         if (v.le.0d0.or.v.gt.1d6.or.itic.gt.20) then

            if (jerk.lt.10) then
               jerk = jerk + 1
               write (*,1000) t,p
               if (jerk.eq.10) call warn (49,r,369,'VDPBM3')
            end if

            vdpbm3 = 1d12*p

            return

         end if

      end do
c                                 and the vdp integral is:
      f = 0.5d0*((vt/v)**r23-1d0)
c                                 checked in BM3_integration.mws
      vdpbm3 = p*v - vt*(pr-4.5d0*k*f**2*(1d0-f*(4d0+kprime)))

1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Birch-Murnaghan ',
     *        'EoS, probably for Ghiorso et al. MELTS/PMELTS endmember',
     *        ' data.',/,
     *        'The affected phase will be destabilized.',/)

      end

      subroutine cartes (wt,ksite,lpoly,ids)
c---------------------------------------------------------------------
c subroutine to cartesian or transform subdivision on site ksite of
c solution ids (or sname). called by subdiv.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jsp, ksite, lpoly, ids

      double precision ycum, wt

      integer ntot,npairs
      common/ cst86 /ntot,npairs
c----------------------------------------------------------------------
      ycum = 0d0
      jsp = ndim(ksite,lpoly,ids)

      if (jsp.eq.0) then
c                                 a relict site with only one species
c                                 left over from reform, i have know idea
c                                 if this works. i doubt it does, until
c                                 jan 8, 2016 it was setting an array (y)
c                                 that was not returned to the calling
c                                 program.
         simp(1) = pxmn(lpoly,ksite,1)
         npairs = 1
         return

      end if

      call chopit (ycum,wt,0,jsp,ksite,lpoly,ids,0,.false.)

      end

      subroutine blanko (text,chars,nchar,ilen)
c-------------------------------------------------------------------
c unblnk - find the last non-blank character in text,
c     text - character string
c-------------------------------------------------------------------
      implicit none

      integer nchar, ilen

      character text*(*), chars(*)*1

      read (text,1000) chars(1:ilen)
c                                 scan for blanks:
      do nchar = ilen, 1, -1
         if (chars(nchar).gt.' ') exit
      end do

1000  format (400a)

      end

      subroutine sattst (ifer,make,good)
c----------------------------------------------------------------------
c sorts phases into the appropriate saturated phase list called by
c input2. returns good if data is valid
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j,ifer,idc

      logical good, make

      character name*8
      common/ csta6 /name

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ic
      common/ cst42 /ic(k0)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ifct,idfl
      common/ cst208 /ifct,idfl
c-----------------------------------------------------------------------

      good = .false.

      if (ifct.gt.0) then
c                               check for fluid species data
         do j = 1, ispec

            if (name.ne.cmpnt(idspe(j))) cycle
            ifer = ifer + 1
            good = .true.
            call loadit (j,.false.,.true.)
            return

         end do

      end if

      if (isat.gt.0) then
c                               check for saturated composants:
c                               reject the phase if it contains
c                               a thermodynamic component:
         do j = 1, icp
            if (comp(ic(j)).ne.0d0) return
         end do
c                               now load the phase if it has
c                               the saturated component idc:
         do j = isat, 1, -1
            idc = icp + j
            if (comp(ic(idc)).ne.0d0) then
               isct(j) = isct(j) + 1
               if (isct(j).gt.h6) call error (17,1d0,h6,'SATTST')
               iphct = iphct + 1
               if (iphct.gt.k1) call error (72,1d0,k1,
     *                            'SATTST increase parameter k1')
               ids(j,isct(j)) = iphct
               call loadit (iphct,make,.true.)
c                                set ltemp1 if a GFSM endmember
               if (ieos.gt.100.and.ieos.lt.200) ltemp1 = .true.
               good = .true.
               return
            end if
         end do

      end if

      end

      double precision function gmake (id)
c-----------------------------------------------------------------------
c gmake computes and sums the component g's for a make definition.
c the component g's may be calculated redundantly because gmake is
c called by gcpd, which in turn may be called by routines that call
c for a single g (e.g., gphase).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id, jd

      double precision g, gcpd

      external gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer mknum, mkind, meos
      double precision mkcoef, mdqf
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer make
      common / cst335 /make(k10)
c-----------------------------------------------------------------------

      jd = make(id)

      g = 0d0
c                                compute the sum of the component g's
      do i = 1, mknum(jd)

         g = g + mkcoef(jd,i) * gcpd (mkind(jd,i),.false.)

      end do
c                                add the dqf correction
      gmake = g + mdqf(jd,1) + t*mdqf(jd,2) + p*mdqf(jd,3)

      end

c       the x-files library for adaptive minimization

c local variables:

c insp(jstot/jstot+1) - pointer to original index of first kstot independent
c         endmembers, followed by jstot-kstot dependent endmembers, followed
c         by the ordered species, if present (jstot+1).
c istot - total number of endmembers excluding ordered species used
c         to formulate the solution model in the input file
c jmsol(jstot,msp) - chemical species present on m'th site of jstot'th
c          endmember, original indexing.
c jstot - total number of endmembers (excluding ordered species) used
c         in the computation (i.e., excluding missing endmembers), but
c         including dependent endmembers.
c kdsol(istot/istot+nord) - endmember index, 0 if missing, -2 if dependent,
c          -1 or 0 if ordered species, reordered for missing species.
c kstot - total number of endmembers (excluding ordered species)
c         used in the computation (i.e., excluding missing endmembers)

c global variables

c indx(i,j,k,l) - for solution i, pointer to the l'th original endmember
c                 index with species k on site j. eliminated in 6.9.1.
c mstot(i) - istot globally
c jgsol(i,j,k) - k species indices of endmember j in solution i (jmsol globally)
c knsp(i=1:mstot,ids) - points to the original (?) index of endmember i in ids

      subroutine reform (im,first)
c---------------------------------------------------------------------
c reform - counts the number of species that can be respresented for a
c solution given the present endmembers.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first

      integer kill,ikill,jkill,kill1,i,j,kosp(mst,msp),kill2,
     *        k,im,idsp,ksp(mst)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot
c----------------------------------------------------------------------

      if (jsmod.eq.20) then

         call reaqus

         kstot = jstot
         istot = jstot

         return

      else if (jsmod.eq.688) then 

         call reforn (im,first)

         return

      end if

      kill = 1

      if (first.and.isimp(1).gt.1) call warn (50,wg(1,1),isimp(1),tname)

      do

         do i = 1, isimp(1)
            ksp(i) = 1
            do j = 1, ivert(1,i)
               kosp(i,j) = 0
            end do
            do j = 1, isimp(1)
c                                 the number of endmembers that
c                                 should be present to represent
c                                 site i if all endmembers are present.
               if (i.ne.j) ksp(i) = ksp(i)*ivert(1,j)
            end do
         end do
c                                 count the number of species
c                                 missing on each site
         do i = 1, istot
            if (kdsol(i).eq.0) then
               do j = 1, isimp(1)
                  kosp(j,jmsol(i,j)) = kosp(j,jmsol(i,j)) + 1
               end do
            end if
         end do
c                                 find the species that is missing
c                                 the most from the model
         kill = 99
         kill1 = 0

         do i = 1, isimp(1)
            if (ivert(1,i).gt.1) then
               do j = 1, ivert(1,i)
c                                 idsp is the the number of species
c                                 possible - the number of missing
c                                 endmembers
                  idsp = ksp(i) -kosp(i,j)

                  if (idsp.lt.kill) then
                     ikill = i
                     jkill = j
                     kill = idsp
                  else if (idsp.eq.kill.and.kosp(i,j).gt.0) then
c                                 kill the species that will kill the
c                                 most dependent endmembers
                     kill2 = 0

                     do k = 1, istot
c                                 count the number that will be killed
                        if (jmsol(k,ikill).eq.jkill.and.kdsol(k).eq.-2)
     *                     kill2 = kill2 + 1
                     end do

                     if (kill2.gt.kill1) then
c                                 this is more than before (kill1)
                        kill1 = kill2
                        ikill = i
                        jkill = j
                        kill = idsp
                     end if
                  end if
               end do
            end if
         end do
c                                 kill the species jkill on site ikill
c                                 and reformulate the model (this is
c                                 inefficient, but we don't care here).
         call killsp (ikill,jkill)
c                                 check exit conditions
         if (istot.lt.2) then
c                                 failed, rejected too many endmembers

            im = im - 1
            if (first) call warn (25,wg(1,1),jstot,tname)
            jstot = 0

            exit

         else if (istot.le.jstot.or.kill.eq.99) then
c                                 changed from eq to le. JADC Nov 22, 2016
c                                 succeeded
            exit
c                                 reorder the insp array so that the
c                                 the first kstot endmembers are
c                                 independent, followed by the
c                                 dependent endmember, followed by
c                                 the ordered species.
         end if

      end do
c                                 eliminate sites with only one
c                                 species
      if (isimp(1).gt.1) call dedsit

      end

      subroutine dedsit
c---------------------------------------------------------------------
c dedsit - eliminates chemical mixing sites with only one species
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,itic,iwas(mst)
c                                 local input variables
      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *      nr(j3)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot
c----------------------------------------------------------------------

      itic = 0

      do i = 1, isimp(1)

         if (ivert(1,i).gt.1) then
            itic = itic + 1
            iwas(itic) = i
         end if

      end do

      if (itic.eq.isimp(1)) return

      isimp(1) = itic
c                                 a dead site, shift the counters and limits:
      do i = 1, isimp(1)

         ivert(1,i) = ivert(1,iwas(i))
c                                 shift subdivision ranges
         do j = 1, ivert(1,i) - 1
            pxmn(1,i,j) = pxmn(1,iwas(i),j)
            pxmx(1,i,j) = pxmx(1,iwas(i),j)
            pxnc(1,i,j) = pxnc(1,iwas(i),j)
            pimd(1,i,j) = pimd(1,iwas(i),j)
         end do

      end do

      do i = 1, istot + norder
c                                 reset the species pointers (jmsol)
         do j = 1, isimp(1)

            jmsol(i,j) = jmsol(i,iwas(j))

         end do

      end do

      if (isimp(1).eq.1) recip = .false.

      if (ordmod) then

         if (isimp(1).eq.1) jsmod = 6

      else if (depmod) then

         jsmod = 7

      else

         jsmod = 2

      end if

      end

      subroutine killsp (ikill,jkill)
c---------------------------------------------------------------------
c killsp - eliminates species jkill from site ikill in a solution model
c and reformulates the model accordingly
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical skip, bad

      integer jsp,jtic,morder,jend,
     *        i,j,ikill,jkill,kill,kdep,kdqf,ktic,jold,
     *        i2ni(m4),kwas(m4),k,l,itic,ijkill(m4),
     *        j2oj(msp),j2nj(msp),i2oi(m4),maxord
c                                 dqf variables
      integer indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c                                 local input variables
      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer nsub,nterm
      double precision acoef
      common/ cst107 /acoef(m10,m11,0:m0),
     *                nterm(m10,m11),nsub(m10,m11,m0)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)
c----------------------------------------------------------------------

      jend = isimp(1)

c                                 was 1,.site
      do i = 1, isimp(1)

         if (i.ne.ikill) then
c                                 nothing happens
            do j = 1, ivert(1,i)
              j2oj(j) = j
            end do

         else
c                                 on a site where we kill a species
            jsp = ivert(1,i) - 1
c                                 should also check and store subdivsion
c                                 model here (i.e., some ternary model
c                                 choices are not valid as binary):

c                                 make a pointer from the new species
c                                 index to the old index
            jtic = 0

            do j = 1, ivert(1,i)

               if (j.ne.jkill) then
                  jtic = jtic + 1
c                              pointer from new j to old j
                  j2oj(jtic) = j
c                              pointer from old j to new j
                  j2nj(j) = jtic
               end if

            end do
c                              now reload
            ivert(1,i) = jsp

            if (jsp.gt.1) then
c                              now shift subdivision ranges
               do j = 1, jsp - 1

                  pxmn(1,i,j) = pxmn(1,i,j2oj(j))
                  pxmx(1,i,j) = pxmx(1,i,j2oj(j))
                  pxnc(1,i,j) = pxnc(1,i,j2oj(j))
                  pimd(1,i,j) = pimd(1,i,j2oj(j))

               end do

            else

               pxmn(1,i,j) = 1d0
               pxmx(1,i,j) = 1d0
               pxnc(1,i,j) = 1d0

            end if
         end if
      end do

      kdep = 0

      if (depmod) then
c                                create an array which gives the
c                                original locations of the dependent
c                                endmembers, need this to be able to
c                                reorder the y2p array:
         do i = 1, istot
            if (kdsol(i).eq.-2) kdep = kdep + 1
         end do

      end if

      do i = 1, istot
c                                 kill endmembers with the species
c                                 to be deleted:
         if (jmsol(i,ikill).eq.jkill) kdsol(i) = -3
      end do
c                                 check for dependent endmembers
      call redep (-3)
c                                 at this point all ordered endmembers to
c                                 be killed are flagged by kdsol = -3.

c                                 now check the ordered species
      morder = 0

      if (ordmod) then
c                                 first check if the ordered endmember
c                                 may be stable
         do k = 1, norder
c                                 check if a missing constituent
            bad = .false.

            do j = 1, nr(k)
               if (kdsol(iddeps(j,k)).eq.-3) then
                  bad = .true.
                  exit
               end if
            end do

            if (bad) then
c                                 add species to the kill list
               kdsol(istot+k) = -3

            else

               morder = morder + 1
               kdsol(istot+k) = -1
               kwas(morder) = k

            end if

         end do

      end if

      itic = 0
      jtic = 0
      ktic = 0
      kill = 0

      do i = 1, istot + norder

         if (kdsol(i).ge.-2) then
c                                 replacement for istot (itic)
            itic = itic + 1
            if (i.le.istot) then
               ktic = ktic + 1
            end if
c                                 pointers from new to old endmember index (i2oi)
            i2oi(itic) = i
c                                 pointers from new to old endmember index (i2ni)
            i2ni(i) = itic
c                                 pointer to original species index
            iorig(itic) = iorig(i)
c                                 number of missing endmembers (jtic)
            if (kdsol(i).eq.0) jtic = jtic + 1
c                                reset the kdsol array
            kdsol(itic) = kdsol(i)

         else
c                                 kill records the killed endmembers
            kill = kill + 1
            ijkill(kill) = i

         end if

      end do

      do i = 1, itic

         if (i2oi(i).gt.istot) cycle
c                                 reset the species pointers (jmsol)
         do j = 1, jend
            if (j.eq.ikill) then
               jmsol(i,j) = j2nj(jmsol(i2oi(i),j))
            else
               jmsol(i,j) = jmsol(i2oi(i),j)
            end if
         end do
      end do
c                                reset total and present counters
      istot = ktic

      jstot = istot - jtic
c                                --------------------------------------
c                                excess terms:
      itic = 0
      maxord = 0

      do i = 1, iterm
c                                check for forbidden terms (i.e., terms
c                                with a missing endmember
         skip = .false.
c                                 macroscopic formulation
         do j = 1, kill
c                                 check if subscript points to a killed
c                                 endmember
            do k = 1, rkord(i)

               if (isub(i,k).eq.ijkill(j)) then
                  skip = .true.
                  exit
               end if

            end do

            if (skip) exit

         end do

         if (skip) cycle
c                               the term is acceptable
         itic = itic + 1

         rkord(itic) = rkord(i)

         isub(itic,1:rkord(i)) = i2ni(isub(i,1:rkord(i)))

         if (xtyp.eq.0) then
c                                save the coefficient
            do j = 1, m3
               wg(itic,j) = wg(i,j)
            end do

         else
c                                 redlich kistler
            do j = 1, rkord(itic)
               do k = 1, m16
                  wk(k,j,itic) = wk(k,j,i)
               end do
            end do

         end if

      end do
c                                reset counter
      iterm = itic
c                                --------------------------------------
c                                van laar volume functions
      if (laar) then
         do i = 1, istot + morder
            do j = 1, m3
               vlaar(j,i) = vlaar(j,i2oi(i))
            end do
         end do
      end if
c                                 --------------------------------------
c                                 dqf corrections, this is sloppy since
c                                 uses istot instead of kstot
      if (idqf.gt.0) then

         kdqf = 0
c                                 check if a retained species has a dqf
c                                 correction
         do j = 1, idqf
c                                 the itoi index must be in the inner loop
c                                 in case the values of indq are not sequential
            do i = 1, istot
               if (indq(j).eq.i2oi(i)) then
c                                 found a dqf'd endmember
                  kdqf = kdqf + 1
                  indq(kdqf) = i
                  do k = 1, m3
                     dqf(k,kdqf) = dqf(k,j)
                  end do
                  exit
               end if
            end do

            if (kdqf.eq.idqf) exit

         end do

         idqf = kdqf

      end if
c                                 --------------------------------------
c                                 configurational entropy model

c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      do i = 1, msite(h0)
c                                 for each species, read function to define
c                                 the site fraction of the species and eliminate
c                                 killed species

c                                 species counter is incremented in advance
c                                 and must be decremented before saving the
c                                 final value:
         jtic = 1

         do j = 1, zsp(h0,i)

            ktic = 0
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 macroscopic formulation:
               dead = .false.
               do l = 1, kill
                  if (nsub(i,j,k).eq.ijkill(l)) then
                     dead = .true.
                     exit
                  end if
               end do

               if (.not.dead) then
c                                 the term has survived (and therefore
c                                 also the species):
                  ktic = ktic + 1
c                                 but my dear peanut brained friend, do
c                                 not forget to move the pointer:
                  nsub(i,jtic,ktic) = i2ni(nsub(i,j,k))
                  acoef(i,jtic,ktic) = acoef(i,j,k)

               end if

            end do
c                                 ktic is the number of terms representing
c                                 the jth species, we won't count species
c                                 with no terms because the endmember configurational
c                                 entropy is assumed to be implicit.
            if (ktic.gt.0) then
c                                 increment the species counter
               znames(h0,i,jtic) = znames(h0,i,j)
               nterm(i,jtic) = ktic
               acoef(i,jtic,0) = acoef(i,j,0)
               jtic = jtic + 1

            end if

         end do

         zsp(h0,i) = jtic - 1
         zsp1(h0,i) = zsp(h0,i)

      end do
c                                 ---------------------------------------
c                                 ordered species:
      if (ordmod) then

         norder = morder

         if (morder.eq.0) then
c                                 there are no ordered species left
            ordmod = .false.

            if (depmod) then

               jsmod = 7

            else
c                                 why jsmod = 2?
               jsmod = 2

            end if

         else
c                                 shift the ordered species pointers
c                                 and data to eliminate kill ordered
c                                 species.
            do j = 1, morder

               jold = kwas(j)

               do i = 1, 3
                  denth(j,i) = denth(jold,i)
               end do

               nr(j) = nr(jold)

               do i = 1, nr(j)
                  iddeps(i,j) = i2ni(iddeps(i,jold))
                  depnu(i,j) = depnu(i,jold)
               end do

            end do

         end if

      end if
c                                 --------------------------------------
c                                 dependent endmember properties, the
      if (depmod) then
c                                 dependent endmembers have been reordered
c                                 in redep, but are still expressed in
c                                 terms of the old indices, so reset the
c                                 indices:
         do i = 1, mdep
            jdep(i) = i2ni(jdep(i))
            do j = 1, ndph(i)
               idep(i,j) = i2ni(idep(i,j))
            end do
         end do

      end if

      end

      subroutine kill01 (site)
c---------------------------------------------------------------------
c reformulates orphan vertex model (jsmod = 9) for missing orphans.
c site is the simplex with the orphan vertices.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ksp, site

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot
c----------------------------------------------------------------------
c                                 count the number of species
c                                 missing on site
      ksp = 0

      do

         do i = 1, ivert(1,site)
            if (kdsol(istot+i).eq.0) then
               ksp = ksp + 1
               call killsp (site,i)
               exit
            end if
         end do

         if (i.gt.ivert(1,site)) exit

      end do

      end

      subroutine reaqus
c---------------------------------------------------------------------
c reaqus - reformulates aqueous solution models to eliminate missing
c endmembers and shift ranges.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jq, jn, js, lm

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------

      jq = 0
      jn = 0
      js = 0
      lm = 0

      do i = 1, ns

         if (kdsol(i).ne.0) then
            js = js + 1
            iorig(js) = i
            kdsol(js) = kdsol(i)
            lm = lm + 1
            pxmn(1,1,lm) = pxmn(1,1,i)
            pxmx(1,1,lm) = pxmx(1,1,i)
            pxnc(1,1,lm) = pxnc(1,1,i)
            pimd(1,1,lm) = pimd(1,1,i)
         end if

      end do

      do i = ns+1, ns+nn

         if (kdsol(i).ne.0) then
            jn = jn + 1
            iorig(js+jn) = i
            kdsol(js+jn) = kdsol(i)
            lm = lm + 1
            pxmn(1,1,lm) = pxmn(1,1,i)
            pxmx(1,1,lm) = pxmx(1,1,i)
            pxnc(1,1,lm) = pxnc(1,1,i)
            pimd(1,1,lm) = pimd(1,1,i)
         end if

      end do

      do i = ns+nn+1, ns+nn+nq

         if (kdsol(i).ne.0) then
            jq = jq + 1
            iorig(js+jn+jq) = i
            kdsol(js+jn+jq) = kdsol(i)

            if (i.eq.ns+nn+nq) cycle
            lm = lm + 1
            pxmn(1,1,lm) = pxmn(1,1,i)
            pxmx(1,1,lm) = pxmx(1,1,i)
            pxnc(1,1,lm) = pxnc(1,1,i)
            pimd(1,1,lm) = pimd(1,1,i)

         end if

      end do

      ns = js
      nn = jn
      nq = jq

      if (ns.eq.0) then

         call warn (99,0d0,0,'rejecting '//tname//' because no solvent'
     *             //' species were identified')
         jstot = 0
         return

      else if (nq.eq.1) then

         call warn (99,0d0,0,'eliminating ions from '//tname//' because'
     *             //' only one charged species was identified')
         nq = 0

      end if

      jstot = ns + nn + nq

      end

      subroutine cmodel (im,idsol,first,found)
c---------------------------------------------------------------------
c cmodel - checks to see if solution models contain valid endmembers.
c modified to allow saturated phase/component endmembers, 10/25/05.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, ok, found

      character missin(m4)*8

      integer imiss, im, idsol, i, j, h, ineg, ipos

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer eos
      common/ cst303 /eos(k10)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character mname*8
      common/ cst18a /mname(m4)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer ikp
      common/ cst61 /ikp(k1)

      integer iam
      common/ cst4 /iam

      integer ixct,ifact
      common/ cst37 /ixct,ifact

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------
      jstot = 0
      ineg = 0
      ipos = 0
      ok = .false.
      found = .false.
c                              if called by build (iam = 4) skip the
c                              name check:
      if (iam.ne.4) then
c                              check that solution name matches a
c                              name requested in input from n1
         do i = 1, isoct

           if (tname.eq.fname(i)) then
c                              got a match, exit
               idsol = i
               im = im + 1
               ok = .true.
               found = ok
               exit

            end if
c
         end do
c                                 didn't find a match, read a new name:
         if (.not.ok) return

      end if

      do i = 1, istot

         kdsol(i) = 0
         ok = .false.
c                              solution with dependent endmembers, if endmember i
c                              is a dependent endmember and not excluded, flag it by setting kdsol(i) = -2
         do j = 1, mdep

            if (jdep(j).eq.i) then
c                              check against exclude list
               do h = 1, ixct
                  if (mname(i).eq.exname(h)) then
                     ok = .true.
                     exit
                  end if
               end do

               if (ok) exit

               kdsol(i) = -2
               ok = .true.
               exit

            end if

         end do

         if (ok) cycle

         if (jsmod.eq.20.and.i.gt.ns) then
c                                 aqueous solute, test against aqnam
            do h = 1, aqct

               if (aqnam(h).eq.mname(i)) then
c                                 got a valid endmember, count and
                  jstot = jstot + 1
                  kstot = jstot
c                                 create arrays of convenience, where j = 1, jstot
                  kdsol(i) = h + aqst
c                                 tests for aq solvent models
                  if (i.gt.nn+ns) then
c                                 must be a charged solute species
                     if (thermo(6,h+aqst).eq.0) then
                        write (*,*) aqnam(h),
     *                  ' is not a charged species '
     *                  ,'remove it from the list of charged species in'
     *                  ,'solution model:',tname
                        call errpau
                        stop
                     else if (thermo(6,h+aqst).gt.0) then
                        ipos = ipos + 1
                     else
                        ineg = ineg + 1
                     end if

                  else if (i.gt.nq) then
c                                  must be neutral species
                     if (thermo(6,h+aqst).ne.0) then
                        write (*,*) names(h),' is a charged species'
     *                            ,'remove it from the list of neutral'
     *                            ,' species in solution model:',tname
                        call errpau
                        stop
                     end if

                  end if

                  exit

               end if

            end do

         else
c                                 must be a real enitity:
            do h = 1, ipoint

               if (names(h).eq.mname(i)) then
c                                 got a valid endmember, count and
                  jstot = jstot + 1
                  kstot = jstot
c                                 create arrays of convenience, where j = 1, jstot
                  kdsol(i) = h
c                                 set kill flag:
c                                 don't reset ikp if it has been set for
c                                 another model.
                  if (iend(i).ne.0.and.ikp(h).eq.0) ikp(h) = -1

                  if (eos(h).eq.15.or.eos(h).eq.16) then

                     write (*,'(a,/)') names(h)//' is not described '//
     *                               'by a solvent EoS remove it from',
     *                               'the list of solvent species in '//
     *                               'solution model:'//tname

                     call errpau

                     stop

                  end if

                  exit

               end if

            end do

         end if
c                                 found all endmembers:
         if (jstot.eq.istot) exit

      end do

      if (jsmod.eq.20) then
c                                 for solvent models, check that charge balance
c                                 is possible
         if ((ipos.eq.0.and.ineg.gt.0).or.
     *       (ipos.gt.0.and.ineg.eq.0)) then

             write (*,'(/,a,/,a)')
     *                         'models with charged species must incl'//
     *                         'ude both postive and negative species',
     *                         'solution model: '//tname//
     *                         ' will be rejected'

             jstot = 1

         end if

      end if

      call redep (0)
c                                done if nothing is missing:
      if (jstot.eq.istot) then
         if (jstot.eq.1.and..not.lopt(32)) call warn (99,wg(1,1),1,
     *                tname//'will be rejected because '//
     *                'aq_lagged_speciation = F')
         return
      end if 
c                                missing endmember warnings:
      if (jstot.lt.2) then

         im = im - 1
         if (first) call warn (25,wg(1,1),jstot,tname)

      else

         imiss = 0

         do i = 1, istot
            if (kdsol(i).eq.0) then
               imiss = imiss + 1
               missin(imiss) = mname(i)
            end if
         end do

         if (first) then 
            write (*,1000) tname,(missin(i), i = 1, imiss)
            write (*,'(/)')
         end if

      end if

1000  format (/,'**warning ver114** the following endmembers',
     *          ' are missing for ',a,//,4(8(2x,a)))

      end

      subroutine redep (jkill)
c----------------------------------------------------------------------
c redep does reordering of dependent endmembers

c jkill is the value of kdsol that indicates invalid (missing) endmembers,
c this is 0 from cmodel and -3 from reform.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,l,ldep,k,jkill

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)
c----------------------------------------------------------------------
c                                check for dependent endmembers, necessary?
      if (depmod) then

         ldep = 0

         do 100 i = 1, mdep

            do j = 1, ndph(i)

               if (idep(i,j).le.istot) then
c                                look for excluded dependent endmembers
                  if (kdsol(jdep(i)).eq.jkill) goto 100
c                                looking for a dependent endmember component:
                  if (kdsol(idep(i,j)).eq.jkill.and.
     *                kdsol(jdep(i)).ne.-3) then
c                                a component is missing in cmodel, reset kdsol to 0
                     kdsol(jdep(i)) = 0
                     goto 100
                  end if

               else
c                                looking for an ordered species
                  do l = 1, norder
                     do k = 1, nr(l)
                        if (kdsol(iddeps(k,l)).eq.jkill) then
                           kdsol(jdep(i)) = 0
                           goto 100
                        end if
                     end do
                  end do

               end if
            end do
c                                dependent endmember is ok
            ldep = ldep + 1
            jdep(ldep) = jdep(i)
            ndph(ldep) = ndph(i)

            do j = 1, ndph(i)
               idep(ldep,j) = idep(i,j)
               nu(ldep,j) = nu(i,j)
            end do

            jstot = jstot + 1

100      continue

         mdep = ldep

         if (mdep.eq.0) depmod = .false.

      end if

      end

      subroutine rmodel (tn1,tn2)
c---------------------------------------------------------------------
c rmodel - reads solution models from LUN n9.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character tn1*6, tn2*22, tag*3, char*1, key*22, val*3,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40,
     *          estrg*80

      integer nreact,i,j,k,l,m,idim

      logical eor

      double precision coeffs(k7), rnums(m4), enth(3), sum

      integer ijk(mst),inds(k7),ict

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      character mname*8
      common/ cst18a /mname(m4)

      integer nsub,nterm
      double precision acoef
      common/ cst107 /acoef(m10,m11,0:m0),
     *                nterm(m10,m11),nsub(m10,m11,m0)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------
      mdep = 0
      norder = 0
      istot = 0

      do i = 1, m4
         kdsol(i) = 0
      end do
c                               set logical variables
      ordmod = .false.
      depmod = .false.
      recip = .false.

      do
          read (n9, '(3(a,1x))', iostat = i) tname
          if (i.ne.0) return
          read (tname,'(a)') char
          if (char.eq.' '.or.char.eq.'|') cycle
          if (tname.ne.'begin_mode') exit
      end do
c                             initialize strings
      tn1 = tname
      tn2 = 'unclassified'
c                             for 1-polytope models use the solution name
c                             for the polytope.
      poname(h0,1,1,1) = tname
c                             next look for optional abbreviation and full_name
      do
         call redcd1 (n9,i,key,val,nval1,nval2,nval3,strg,strg1)

         read (key,'(i3)', iostat = i) jsmod

         if (i.eq.0) then
c                             read jsmod, model type flag
            exit

         else if (key.eq.'abbreviation') then

            tn1 = strg1

         else if (key.eq.'full_name') then

            tn2 = strg1

         else

            write (estrg,'(4(a))') 'unrecognized text: ', key,
     *                             ' reading solution model ', tname
            call error (72, enth(1), i, estrg)

          end if

      end do

      if (jsmod.eq.20) then
c                                 aqueous model reads to different
c                                 arrays
         call raqmod
         istot = nq + nn + ns

         return

       else if (jsmod.eq.688) then 

         call rmoden

         return

       end if
c                                 check for disabled model types
      if (jsmod.eq.1.or.jsmod.eq.3.or.jsmod.eq.5)
     *                        call error (68,enth(1),jsmod,tname)
      if (jsmod.eq.6.or.jsmod.eq.8.or.jsmod.eq.9.or.
     *                                jsmod.eq.27) ordmod = .true.
      if (jsmod.ge.7.and.jsmod.le.10) depmod = .true.
      if (jsmod.ge.7.and.jsmod.le.10) recip = .true.
c                                 assign non-default props to
c                                 special models:
      if (jsmod.ge.30.and.jsmod.le.31) recip = .true.
c                                 read number of independent sites:
      if (recip) then
         call readda (rnums,1,tname)
         isimp(1) = idint(rnums(1))
      else
         isimp(1) = 1
      end if
c                                 this reads pre 6.8.8 format models
      poly(h0) = 1
c                                 total number of endmembers that define the 
c                                 polytope (does not include ordered endmembers):
      istot = 1
c                                 read number of species on each site:
      call readda (rnums,isimp(1),tname)
      do i = 1, isimp(1)
         ivert(1,i) = idint(rnums(i))
         istot = istot*ivert(1,i)
      end do
c                               counter for character read routines
c                               and starting index
      idim = istot

      if (istot.gt.m4) call error (1,rnums(1),idim,
     *                 'm4 (maximum number of endmembers)')

      i = 0
c                               read endmember names:
      call readn (i,idim,tname)
c                               compound formation models
      if (ordmod) then
c                               get the number of ordered species
         call readda (rnums,1,tname)
         norder = idint(rnums(1))

         if (istot+norder.gt.m4) call error (1,rnums(1),istot+norder,
     *                           'm4 (maximum number of endmembers)')

         if (norder.gt.j3) call error (5,rnums(1),norder,tname)
c                               get ordering reaction and name
c                               of ordered species:
         do i = 1, norder

            nreact = -1

            call readr (coeffs,enth,inds,idim,nreact,tname,eor)

            do j = 1, 3
               denth(i,j) = enth(j)
            end do

            nr(i) = nreact

            do j = 1, nreact
               depnu(j,i) = coeffs(j+1)
               iddeps(j,i) = inds(j+1)
            end do

         end do
c                               dummy routine for backward compatability:
c                               read the limit equations for the
c                               amount of the ordered endmembers
         call readlm (idim,tname)

      end if
c                               read dependent endmembers
      if (depmod) then
c                               number of dependent endmembers
         call readda (rnums,1,tname)
         mdep = idint(rnums(1))
         if (mdep.gt.m15) call error (1,enth(1),mdep,'m15')

         do i = 1, mdep
c                               nreact is returned by readr
            nreact = 0

            call readr (coeffs,enth,inds,idim,nreact,tname,eor)

            jdep(i) = inds(1)
            ndph(i) = nreact - 1
            if (ndph(i).gt.j4) call error (1,enth(1),ndph(i),'j4')

            sum = 0d0

            do j = 1, ndph(i)
               sum = sum + coeffs(j+1)
               nu(i,j) = coeffs(j+1)
               idep(i,j) = inds(j+1)
            end do

            if (dabs(sum-1d0).gt.nopt(50)) then

               write (*,'(/,a,g12.6,/)') 'coefficient sum = ', sum

               call error (72,sum,i,'dependent endmember '//
     *              mname(inds(1))//' definition coefficients do not'
     *              //'sum to 1, solution model:'//tname)

            end if

         end do

      end if
c                               read endmember flags:
      call readda (rnums,istot,tname)

      do i = 1, istot
         iend(i) = idint(rnums(i))
      end do
c                               read composition limits, subdivision type:
      do i = 1, isimp(1)

         do j = 1, ivert(1,i) - 1

            call readda (rnums,4,tname)

            pxmn(1,i,j) = rnums(1)
            pxmx(1,i,j) = rnums(2)
            pxnc(1,i,j) = rnums(3)
            pimd(1,i,j) = idint(rnums(4))

         end do

      end do
c                                create bragg-williams indexes
      do i = 2, isimp(1)
         ijk(i) = 1
      end do

      ijk(1) = 0

      do 20 l = 1, istot

         do m = 1, isimp(1)

            if (ijk(m).lt.ivert(1,m)) then

               ijk(m) = ijk(m) + 1
c                                increment only one index per endmember
               do i = 1, isimp(1)
                  jmsol(l,i) = ijk(i)
               end do

               if (ijk(m).eq.ivert(1,m).and.
     *             l.lt.istot-ivert(1,1)+1) then
c                                saturated counter, increment first
c                                unsaturated counter and reset all
c                                lower site counters to first index,
                  do j = 2, isimp(1)
                     if (ijk(j).lt.ivert(1,j)) then
                        ijk(j) = ijk(j) + 1
                        do k = 2, j-1
                           ijk(k) = 1
                        end do
                        ijk(1) = 0
                        goto 20
                     end if
                  end do
               end if

               goto 20

            end if
         end do
20    continue
c                                 read excess function
      call readx (idim,tname)
c                                 expansion for S(configurational)
      call readda (rnums,1,tname)

      msite(h0) = idint(rnums(1))
      if (msite(h0).gt.m10) call error (1,acoef(1,1,1),msite(h0),'m10')
c                                 for each site
      do i = 1, msite(h0)
c                                 read # of species, and site
c                                 multiplicty.
         call readda (rnums,2,tname)

         zmult(h0,i) = rnums(2)
c                                 if multiplicity is 0, the model
c                                 site has variable multiplicity
c                                 and molar site population expressions
c                                 are read rather than site fractions
c                                 in which case we need as many expressions
c                                 as species. nspm1 is the counter for the
c                                 number of expressions
         zsp1(h0,i) = idint(rnums(1))

         if (zmult(h0,i).gt.0) then
            zsp(h0,i) = zsp1(h0,i) - 1
            zsp1(h0,i) = zsp(h0,i)
         else
            zsp(h0,i) = zsp1(h0,i)
         end if

         if (zsp(h0,i).gt.m11) call error (1,0d0,zsp(h0,i),'m11')
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, zsp(h0,i)
c                                 read expression for site
c                                 fraction of species j on
c                                 site i.
            call readz (coeffs,inds,ict,idim,tname,tag)

            acoef(i,j,0) = coeffs(1)
            nterm(i,j) = ict - 1
            if (nterm(i,j).gt.m0) call error (33,0d0,m0,tname)
c                                 for each term:
            do k = 2, ict
c                                 all terms 1 order type, this
c                                 saved for compatability with
c                                 old versions:
               acoef(i,j,k-1)   = coeffs(k)
               nsub(i,j,k-1)    = inds(k)

            end do

         end do

      end do
c                                 look for van laar and/or dqf parameters
c                                 or the end of model marker
      call readop (idim,istot-mdep,tname)

      if (laar.and.iord.gt.2) call error (999,coeffs(1),800,'RMODEL')
c                                 save original indices, need this for
c                                 melt models etc that have species specific
c                                 equations of state.
      do i = 1, istot + norder
         iorig(i) = i
      end do

      end

      subroutine raqmod
c---------------------------------------------------------------------
c raqmod - reads aq electrolyte solution model, jsmod = 20.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j

      double precision rnums(10)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------
c                               read number of solvent species:
      call readda (rnums,1,tname)
      ns = idint(rnums(1))
c                               read solvent species names:
      i = 0
      if (ns.gt.0) call readn (i,ns,tname)
c                               read number of neutral species:
      call readda (rnums,1,tname)
      nn = idint(rnums(1))
c                               read neutral species names:
      i = ns
      if (nn.gt.0) call readn (i,nn,tname)
c                               read number of charged species:
      call readda (rnums,1,tname)
      nq = idint(rnums(1))
c                               read charged species names:
      i = nn + ns
      if (nq.gt.0) call readn (i,nq,tname)
c                               read composition limits, subdivision type
c                               for (nn + ns + nq) - 2 species
      if (i+nq.eq.2) i = i + 1
c
      do j = 1, i + nq - 1
c                               dummy for the ns'th species
         if (j.eq.ns) then

            pxmn(1,1,j) = 0d0
            pxmx(1,1,j) = 1d0
            cycle

         end if

         call readda (rnums,4,tname)

         pxmn(1,1,j) = rnums(1)
         pxmx(1,1,j) = rnums(2)
         pxnc(1,1,j) = rnums(3)
         pimd(1,1,j) = idint(rnums(4))

      end do
c                              look for van laar and/or dqf parameters
c                              or the end of model marker
      call readop (j,j,tname)

      do i = 1, nq+nn+ns
         iorig(i) = i
      end do

      end

      subroutine nmodel
c---------------------------------------------------------------------
c nmodel - creates some counters and conversion arrays and stores
c local solution model parameters in global arrays.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer itic,i,j,k

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)
c----------------------------------------------------------------------
c                                 make the insp arrays,

c                                 insp orders endmembers; first independent
c                                 then dependent, last ordered species.

c                                 jnsp points only to independent+ordered
c                                 species (used for reciprocal solutions
c                                 with ordering).

c                                 iy2p points from an endmember in the y
c                                 array (1..istot+norder) to the endmember
c                                 in the p array (1..kstot+norder)

      if (depmod) then

         itic = 0

         do i = 1, istot
            if (kdsol(i).gt.0) then
               itic = itic + 1
               insp(itic) = i
               jnsp(itic) = i
               iy2p(i) = itic
            end if
         end do

         do i = 1, mdep
            insp(itic+i) = jdep(i)
         end do

         do i = 1, norder
            insp(istot+i) = istot + i
            jnsp(itic+i) = istot + i
            iy2p(istot+i) = itic + i
         end do

      else

         do i = 1, istot + norder
            insp(i) = i
            jnsp(i) = i
            iy2p(i) = i
         end do
c                                added for reformulated reciprocal
c                                solutions with no dependent endmembers
c                                June 12, 2012, JADC.
c        if (recip) kstot = istot
c                                JADC Nov 22, 2016:
c                                changed to general case, i.e., if there
c                                are no dependent endmembers then
         kstot = istot

      end if

      if (depmod) then
c                                 make y2p array
         do i = 1, kstot + norder

            do j = 1, mdep

               y2p(i,j) = 0d0

               do k = 1, ndph(j)
                  if (idep(j,k).eq.jnsp(i)) y2p(i,j) = nu(j,k)
               end do

            end do

         end do

      end if

      end

      recursive double precision function gproj (id)
c-----------------------------------------------------------------------
c gproj - computes projected free energy of phase id. uproj must be
c called prior to any call to gproj.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,j

      double precision gphase, gcpd

      external gphase, gcpd

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision cp
      common/ cst12 /cp(k5,k10)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn
c---------------------------------------------------------------------
      if (id.le.ipoint) then

         gproj = gcpd (id,.true.)
c                                 if istct > 0 must be some saturated
c                                 components
         if (istct.gt.1) then
c                                 this is a screw up solution
c                                 necessary cause uf(1) and uf(2)
c                                 are allocated independent of ifct!
            if (ifct.gt.0) then
               do j = 1, 2
                  if (iff(j).ne.0) gproj = gproj - cp(iff(j),id)*uf(j)
               end do
            end if

            do j = icp1, icp + isat
               gproj = gproj - cp(j,id) * mu(j)
            end do

         end if

      else

         gproj = gphase(id)

      end if

      end

      double precision function gsol1 (id,order)
c-----------------------------------------------------------------------
c gsol1 computes the complete gibbs energy of a solution identified by 
c index ids and endmember composition pa input from cxt7. 

c if order is true, then the speciation of order disorder models is 
c computed; otherwise the speciation is taken as input.

c summary of the gsol functions/subroutines: 

c gsol  - assumes the endmember g's have not been calculated by gall and is
c         only called by WERAMI/MEEMUM/FRENDLY via function ginc. may be
c         called for any type of solution model.

c gsol1 - identical to gsol but can only been called after gall and is
c         only called by minfrc. ingsol must be called prior to
c         gsol1 to initialize p-t dependent model parameters. may be
c         called for any type of solution model. if gsol1 is called 
c         with order = .true. it returns the p0 normalized g, otherwise
c         p0 = pa. minfrc sets order = F except for scatter points.

c gsol2 - a shell to call gsol1 from minfrc, ingsol must be called
c         prior to minfrc to initialize solution specific paramters. may be
c         called for any type of solution model. gsol2 sets order = F.

c gsol4 - a shell to call gsol1 from minfxc, ingsol must be called
c         prior to minfxc to initialize solution specific paramters. only
c         called for implicit o/d models. p0 normalization

c gord  - a function to compute the excess + enthalpic + entropic effects of 
c         mixing for an implicit o/d solution model. called by gall (via specis), 
c         gsol4. pa normalization.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      logical order, bad

      double precision gg

      double precision omega, gfluid, gzero, gmech, gord, gdqf, gmech0,
     *                 gex, gfesi, gfesic, gfecr1, gerk, ghybrid, gfes

      external omega, gfluid, gzero, gex, gfesi, gmech, gord, gmech0,
     *         gfesic, gfecr1, gerk, ghybrid, gfes, gdqf

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision g
      common/ cst2 /g(k1)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c----------------------------------------------------------------------
      gg = 0d0
c                                 quack flag needed to distinguish phases
c                                 computed by lagged speciation
      rkwak = .true.

      if (specil(id)) then
c                                 special is reserved for special models 
c                                 that set standard flags (lorder and/or lrecip)
c                                 currently only Nastia's version of BCC/FCC Fe-Si-C Lacaze and Sundman
         gg =  gfesic (pa(1),pa(3),pa(4),
     *                 g(jend(id,3)),g(jend(id,4)),
     *                 g(jend(id,5)),g(jend(id,6)),ksmod(id))

      else if (simple(id)) then
c                                 -------------------------------------
c                                 macroscopic formulation for normal solutions.
         gg = gdqf (id) - t * omega (id,pa) + gex (id,pa) + gmech (id)

      else if (lorder(id).and.order) then
c                                 get the speciation, excess and entropy effects.
         if (.not.noder(id)) then

            call specis (gg,id)

         else

            call minfxc (gg,id,.false.)

         end if
c                                 for non-equimolar all the terms here
c                                 are computed for the pp mass, this is 
c                                 ok because the g will be normalized 
c                                 by the molar mass computed from the p0
c                                 array.
         gg = gg + gdqf(id) + gmech(id)

      else if (lorder(id)) then
c                                 add entropic + enthalpic + excess o/d effect
c                                 gdqf and gmech are for the p0 mass, gord is for
c                                 the pa mass, however this case is only called
c                                 when p0 = pa
         gg = gdqf(id) + gmech(id) + gord(id)

      else if (ksmod(id).eq.0) then
c                                 ------------------------------------
c                                 internal fluid eos
         gg = gfluid(pa(1)) + gmech0(id)

      else if (ksmod(id).eq.20) then
c                                 electrolytic solution, need to check
c                                 that this thing is getting the right
c                                 partial molar volumes.
         call slvnt1 (gg)

         call slvnt2 (gg)

      else if (ksmod(id).eq.26) then
c                                 ------------------------------------
c                                 andreas salt model
         call hcneos (gg,pa(1),pa(2),pa(3))

         gg = gg + gmech(id)

      else if (ksmod(id).eq.29) then
c                                 -------------------------------------
c                                 BCC Fe-Si Lacaze and Sundman
         gg =  gfesi(pa(1),g(jend(id,3)),g(jend(id,4)))

      else if (ksmod(id).eq.32) then
c                                 -------------------------------------
c                                 BCC Fe-Cr Andersson and Sundman
         gg =  gfecr1(pa(1),g(jend(id,3)),g(jend(id,4)))

      else if (ksmod(id).eq.39) then

         bad = .true.
c                                 -------------------------------------
c                                 generic hybrid EoS
         if (lopt(32)) then 
c                                 lagged speciation
c                                 the last argument cancels recalc, in
c                                 which case id is a dummy. smo the total
c                                 species molality it is necessary for 
c                                 renormalization.
            call gaqlgd (gg,rcp,rsum,rsmo,id,bad,.false.)

            if (.not.bad) rkwak = .false.

         end if
c                                 if bad OR molecular:
         if (bad) gg = ghybrid(pa) + gmech(id)

      else if (ksmod(id).eq.41) then
c                                 hybrid MRK ternary COH fluid
         call rkcoh6 (pa(2),pa(1),gg)

         gg = gg + gmech(id)

      else if (ksmod(id).eq.40) then
c                                 MRK silicate vapor
         gg = gmech0(id) + gerk(pa)

      else if (ksmod(id).eq.42) then
c                                 ------------------------------------
c                                 Fe-S fluid (Saxena & Eriksson 2015)
         gg =  gfes(pa(2),g(jend(id,3)),g(jend(id,4)))

      else

         write (*,*) 'what the **** am i doing here?'
         call errpau

      end if
c                                 set the bulk composition for non-lagged
c                                 speciation:
      if (rkwak) call getscp (rcp,rsum,rids,rids)

      gsol1 = gg

      end

      double precision function gord (id)
c-----------------------------------------------------------------------
c gord computes the excess + configurational + enthalpy_of_ordered_endmember 
c free energy of a solution identified by index ids for the speciation input
c via cxt7.

c ingsol MUST be called prior to gord to initialize solution
c specific parameters! 

c gord assumes the endmember g's have been calculated by gall.

c gord is called for explicit order-disorder models by minfxc, gsol1, and 
c specis.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      double precision omega, gex

      external omega, gex

      double precision enth
      common/ cxt35 /enth(j3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
c                                 entropic + excess o/d effect
      gord = - t * omega(id,pa) + gex(id,pa)
c                                 enthalpic effect o/d effct
      do k = 1, nord(id)
         gord = gord + pa(lstot(id)+k)*enth(k)
      end do

      end

      double precision function gmech (id)
c-----------------------------------------------------------------------
c gmech computes the mechanical mixture gibbs energy of a solution, if 
c the solution is an explicit o/d solution the sum does not include the 
c enthalpies of ordering for the ordered endmembers. 

c gmech assumes endmember g have been computed by a prior call to gall.
c gmchpt or gmchpr should be called when this is not the case.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision g
      common/ cst2 /g(k1)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      gmech = 0d0

      do k = 1, lstot(id)
         gmech = gmech + g(jend(id,2+k)) * pp(k)
      end do

      end

      double precision function gmchpt (id)
c-----------------------------------------------------------------------
c gmchpt computes the mechanical mixture gibbs energy of a solution, if 
c the solution is an explicit o/d solution the sum does not include the 
c enthalpies of ordering for the ordered endmembers. 

c gmchpt projects through mobile component potentials, but does not
c project through saturated components. use gmchpr for full projection.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      double precision gcpd

      external gcpd

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      gmchpt = 0d0

      do k = 1, lstot(id)
         gmchpt = gmchpt + gcpd (jend(id,2+k),.true.) * pp(k)
      end do

      end

      subroutine ingmfx (id) 
c-----------------------------------------------------------------------
c initialize endmember g's for minfxc if it is to be called as a bail-
c out routine for specis without prior call to gall
c-----------------------------------------------------------------------
      integer id
c-----------------------------------------------------------------------
c                                 projected endmember g's at p-t
      call geeend (id)
c                                 endmember dqfs at p-t
      call setdqf (id)
c                                 load into the endmember array
      call ingend (id)

      end 


      subroutine geeend (id)
c-----------------------------------------------------------------------
c geeend updates the projected free energies of solution id, it is 
c only called for output purposes (gsol/ginc/getphp).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      double precision gproj

      external gproj

      double precision g
      common/ cst2 /g(k1)

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c----------------------------------------------------------------------

      do k = 1, lstot(id)
         g(jend(id,2+k)) = gproj (jend(id,2+k))
      end do

      end

      double precision function gmchpr (id)
c-----------------------------------------------------------------------
c gmech computes the mechanical mixture gibbs energy of a solution, if 
c the solution is an explicit o/d solution the sum does not include the 
c enthalpies of ordering for the ordered endmembers. 

c gmchpr projects through mobile component potentials and saturated 
c components, called only for static cpds by CONVEX (gphase).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      double precision gproj

      external gproj

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      gmchpr = 0d0

      do k = 1, lstot(id)
         gmchpr = gmchpr + gproj (jend(id,2+k)) * pp(k)
      end do

      end

      double precision function gmech0 (id)
c-----------------------------------------------------------------------
c gmech computes the mechanical mixture gibbs energy of a solution at
c the refernce pressure.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      double precision gzero

      external gzero

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      gmech0 = 0d0

      do k = 1, lstot(id)
         gmech0 = gmech0 + gzero (jend(id,2+k)) * pp(k)
      end do

      end

      subroutine ingsol (id)
c-----------------------------------------------------------------------
c ingso1 initializes p-t dependent solution model id parameters for gsol1
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id
c----------------------------------------------------------------------
c                                 evaluate margules coefficients
      call setw (id)
c                                 evaluate dqf coefficients
      call setdqf (id)
c                                 load enthalpies of ordering
c                                 and load gend on the off-chance
c                                 that minfxc is used
      if (lorder(id)) call oenth (id)

      end

      logical function zbad (y,ids,z,text,endtst,text1)
c----------------------------------------------------------------------
c subroutine to compute site fractions computed from equations entered by
c user for configurational entropy (macroscopic form). with range checks.
c
c non-temkin models:
c     z(i,j) - molar site fraction of species j on site i.
c temkin models:
c     z(i,j) - molar amount of species j on site i.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical badz, bad, endtst

      external badz

      double precision y(*), zt, z(m10,m11)

      integer i,j,k,ids

      character text*(*), text1*(*)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
      bad = .false.
c                                 for each site
      do i = 1, msite(ids)

         zt = 0d0

         if (zmult(ids,i).ne.0d0.and.ksmod(ids).ne.688) then
c                                 get site fractions
            do j = 1, zsp(ids,i)

               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) +
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do

               if (badz(z(i,j))) bad = .true.

               zt = zt + z(i,j)

            end do

            z(i,j) = 1d0 - zt

            if (badz(z(i,j))) bad = .true.

         else if (zsp1(ids,i).gt.1) then
c                                 temkin or 688 model format, all species fractions are available
            do j = 1, zsp1(ids,i)
c                                 molar site population
               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) + 
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do
c                                 non-temkin (688)
               if (zmult(ids,i).gt.0d0.and.badz(z(i,j))) then 

                  if (endtst) then

                     write (*,1000) text1,z(i,j),text

                     call warn (72,zt,i,' z('//
     *                       znames(ids,i,j)//') on '//znames(ids,i,0)//
     *                       ' in '//text//' is invalid.')

                  end if 

                  bad = .true.

               end if

               zt = zt + z(i,j)

            end do

            if (ksmod(ids).eq.688.and.zmult(ids,i).gt.0d0) then 
c                                 non-temkin, fractions must sum to 1
               if (dabs(zt-1d0).gt.nopt(50)) then

c                 write (*,'(/,a,g14.6)') 'site fraction sum = ',zt

                  if (endtst) write (*,1000) text1,zt,text

                  bad = .true.

               end if

            else if (zt.gt.0d0) then
c                                 temkin, if site exists, check fractions
               do j = 1, zsp(ids,i)

                  if (badz(z(i,j)/zt)) bad = .true.

               end do

            else if (zt.lt.-nopt(50)) then
c                                 negative site?
               bad = .true.

            end if

         end if

      end do

      if (boundd(ids)) then 
c                                 molecular entropy, check fractions
          do i = 1, nstot(ids)
             if (y(i).lt.-nopt(50)) then
                bad = .true.
                exit
             else if (y(i).lt.0d0) then
                y(i) = 0d0
             end if
          end do

      end if

      zbad = bad

1000  format (/,'**error ver071** during testing of dependent endmember'
     *       ,' ',a,' the following invalid site fraction (z = ',g12.6,
     *        ')',/,'was found. The cause of this error may be either ',
     *       'the dependent endmember definition or invalid site',/,
     *       'fraction expressions for one or more of the independent ',
     *       'endmembers of ',a,/)

      end

      double precision function omega (id,y)
c----------------------------------------------------------------------
c subroutine to evaluate the delta configurational entropy of a solution
c with composition y, including the correction for endmember
c configurational negentropies. reciprocal end-member composition version.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision z,zt,dlnw,dlnz,y(*),n(m11)

      integer i,j,k,id
c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
      dlnw = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         dlnz = zt

         if (zmult(id,i).ne.0d0) then
c                                 standard model with fixed site multiplicity
c                                 get site fractions
            do j = 1, zsp(id,i)

               z = dcoef(0,j,i,id)
c                                 for each term:
               do k = 1, lterm(j,i,id)
                  z = z + dcoef(k,j,i,id) * y(ksub(k,j,i,id))
               end do

               call ckzlnz (z,dlnz)

               zt = zt + z

            end do

            z = 1d0 - zt
            call ckzlnz (z,dlnz)

            dlnw = dlnw - zmult(id,i)*dlnz

         else if (zsp(id,i).gt.1) then
c                                 variable site multiplicities
c                                 get site fractions
            do j = 1, zsp(id,i)

               n(j) = dcoef(0,j,i,id)
c                                 for each term:
               do k = 1, lterm(j,i,id)
c                                 n(j) is molar site population
                  n(j) = n(j) + dcoef(k,j,i,id) * y(ksub(k,j,i,id))
               end do
c                                 zt is the multiplicity
               zt = zt + n(j)

            end do

            if (zt.gt.0d0) then
c                                 if site has non-zero multiplicity
               do j = 1, zsp(id,i)
c                                 z is site fraction
                  z = n(j)/zt
                  call ckzlnz (z,dlnz)

               end do

            end if

            dlnw = dlnw - r*zt*dlnz

         end if

      end do
c                                 endmember corrections
      do i = 1, nstot(id)
         dlnw = dlnw - y(i)*scoef(i,id)
      end do

      omega = dlnw

      end

      double precision function omega0 (id,y)
c----------------------------------------------------------------------
c subroutine to evaluate the absolute configurational entropy of a solution
c with composition y. used only for p2yx inversion.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision z,zt,dlnw,dlnz,y(*),n(m11)

      integer i,j,k,id
c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
      dlnw = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         dlnz = zt

         if (zmult(id,i).ne.0d0) then
c                                 standard model with fixed site multiplicity
c                                 get site fractions
            do j = 1, zsp(id,i)

               z = dcoef(0,j,i,id)
c                                 for each term:
               do k = 1, lterm(j,i,id)
                  z = z + dcoef(k,j,i,id) * y(ksub(k,j,i,id))
               end do

               call ckzlnz (z,dlnz)
               zt = zt + z

            end do

            z = 1d0 - zt
            call ckzlnz (z,dlnz)

            dlnw = dlnw - zmult(id,i)*dlnz

         else if (zsp(id,i).gt.1) then
c                                 variable site multiplicities
c                                 get site fractions
            do j = 1, zsp(id,i)

               n(j) = dcoef(0,j,i,id)
c                                 for each term:
               do k = 1, lterm(j,i,id)
c                                 n(j) is molar site population
                  n(j) = n(j) + dcoef(k,j,i,id) * y(ksub(k,j,i,id))
               end do
c                                 zt is the multiplicity
               zt = zt + n(j)

            end do

            if (zt.gt.0d0) then
c                                 if site has non-zero multiplicity
               do j = 1, zsp(id,i)
c                                 z is site fraction
                  z = n(j)/zt

                  call ckzlnz (z,dlnz)

               end do

            end if

            dlnw = dlnw - r*zt*dlnz

         end if

      end do

      omega0 = dlnw

      end

      subroutine snorm0 (id,tname)
c------------------------------------------------------------------------
c compute endmember configurational entropies
c site fractions expressed as a function of end-member proportions.
c see deleted snorm0 for bragg-williams site fraction version.
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 tname

      logical zbad

      integer h,id

      double precision omega, zsite(m10,m11)

      external omega, zbad
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 get normalization constants
c                                 for each endmember
      do h = 1, nstot(id)
c                                 zero y-array
         pa(1:nstot(id)) = 0d0

         pa(h) = 1d0
c                                 check for valid site fractions
         if (stck.and.zbad(pa,id,zsite,tname,.false.,tname)) 
     *                                     call error (125,z(1),h,tname)
c                                 evaluate S
         scoef(h,id) = omega(id,pa)

      end do

      end

      subroutine setdqf (id)
c---------------------------------------------------------------------
c setdqf - evaluates dqf coefficients at p and t
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id

      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      do i = 1, jdqf(id)
c                                 index points to the endmember in the full
c                                 model:
         iq(i) = jndq(i,id)

         dq(i) = dqfg(1,i,id) + t*dqfg(2,i,id) + p*dqfg(3,i,id)

      end do

      end

      double precision function gdqf (id)
c----------------------------------------------------------------------
c subroutine to evaluate the excess G of solution id with endmember
c composition pp. assumes a previous call to setdqf.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      gdqf = 0d0 

      do i = 1, jdqf(id)

         gdqf = gdqf + pp(iq(i))*dq(i)

      end do

      end

      subroutine setw (id)
c---------------------------------------------------------------------
c setw - evaluates margules coeffiecients and, for laar models, size
c parameters for solution model id at p and t
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, k, l, i1, i2, id, j

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      double precision wgl, wkl, vlar
      common/ cxt2r /wgl(m3,m1,h9),wkl(m16,m17,m18,h9),vlar(m3,m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c                                 local alpha
      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 bookkeeping variables
      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)
c----------------------------------------------------------------------
      if (extyp(id).eq.1) then
c                                 redlich kistler is a special case
c                                     wk(1) = w0 cst
c                                     wk(2) = wT coefficient on T
c                                     wk(3) = wP some term in brosh's murnaghan-like excess term
c                                     wk(4) = wP1 some term in brosh's murnaghan-like excess term
c                                     wk(5) = wP2 some term in brosh's murnaghan-like excess term
c                                     wk(6) = wP0 coefficient on P
         do i = 1, jterm(id)
            do j = 1, rko(i,id)

               if (wkl(3,j,i,id).eq.0d0.or.wkl(4,j,i,id).eq.0d0.or.
     *             wkl(5,j,i,id).eq.0d0) then

                   wl(j,i) = wkl(1,j,i,id) + t*wkl(2,j,i,id) +
     *                                wkl(6,j,i,id)*p
               else
                   wl(j,i) = wkl(1,j,i,id) + t*wkl(2,j,i,id) +
     *                   4d0*((-1d0*wkl(5,j,i,id) -
     *                   dsqrt((2d0*wkl(5,j,i,id)*p +
     *                   wkl(4,j,i,id))/wkl(4,j,i,id)))*
     *                   wkl(3,j,i,id)*wkl(4,j,i,id)*dexp(-(-1d0 +
     *                   dsqrt((2d0*wkl(5,j,i,id)*p +
     *                   wkl(4,j,i,id))/wkl(4,j,i,id)))/
     *                   wkl(5,j,i,id)) +
     *                   wkl(3,j,i,id)*wkl(4,j,i,id)*
     *                   (wkl(5,j,i,id)+1d0))
               end if

            end do
         end do

         return

      end if

      do i = 1, jterm(id)
         w(i) = wgl(1,i,id) + t*wgl(2,i,id) + p*wgl(3,i,id)
      end do

      if (llaar(id)) then

         do i = 1, nstot(id)

            alpha(i) = vlar(1,i,id)
     *               + t * vlar(2,i,id) + p * vlar(3,i,id)

         end do

         do i = 1, jterm(id)
            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)
            w(i) = 2d0 * w(i)
     *                 * alpha(i1)*alpha(i2) / (alpha(i1) + alpha(i2))
         end do

      end if

      if (lorder(id)) then
c                                 set higher order derivatives for
c                                 speciation
         dt(1:nord(id)) = 0d0
         d2gx(1:nord(id),1:nord(id)) = 0d0
c                                 for both laar and regular need
c                                 the d(gex)/dy(k)/dy(l)
         do i = 1, jterm(id)
            do k = 1, nord(id)
               do l = 1, nord(id)
                  d2gx(l,k) = d2gx(l,k) + w(i) * dppp(l,k,i,id)
               end do
            end do
         end do


         if (llaar(id)) then
c                                 for laar also need:
            do i = 1, nstot(id)
               do k = 1, nord(id)
c                                 dt, derivative of sum(phi)
                  dt(k) = dt(k) + alpha(i)*dydy(i,k,id)
               end do
            end do

         end if

      end if

      end

      double precision function gex (ids,y)
c-----------------------------------------------------------------------
c evaluate the excess function for solution model ids. assumes prior
c call to setw

c input:
c      ids - solution pointer
c      y - composition array
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,ids

      double precision y(*), tphi, xpr, lex(m17,m18)

      double precision z, pa, p0a, x, w, yy, wl, pp
      common/ cxt7 /yy(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)
c                                 local alpha
      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c----------------------------------------------------------------------
      gex = 0d0

      if (extyp(ids).eq.1) then
c                                 redlich kistler; expand polynomial
         do i = 1, jterm(ids)

            lex(1:rko(i,ids),i) = 0d0

            do j = 1, rko(i,ids)
c                                 interchanged subscripts, G Helffrich, 4/8/16.
               lex(j,i) = lex(j,i) + wl(j,i)*
     *                        (y(jsub(1,i,ids))-y(jsub(2,i,ids)))**(j-1)
            end do

         end do

         do i = 1, jterm(ids)
            do j = 1, rko(i,ids)
               gex = gex + lex(j,i)*y(jsub(1,i,ids))*y(jsub(2,i,ids))
            end do
         end do

      else if (lexces(ids)) then

         if (llaar(ids)) then
c                                 holland & powell's version of the van laar
c                                 first compute "volumes"
            tphi = 0d0

            do i = 1, nstot(ids)
c                                 tphi is the sum of holland & powell's
c                                 phis
               tphi = tphi + alpha(i) * y(i)

             end do
c                                 dg is initialized as gph in the calling
c                                 program
            do i = 1, jterm(ids)
c                                 holland powell form, all terms regular
               gex = gex + w(i) * y(jsub(1,i,ids)) * y(jsub(2,i,ids))

            end do

            gex = gex/tphi

         else
c                                 macroscopic margules formulation by default
            do i = 1, jterm(ids)

               xpr = 1d0

               do j = 1, rko(i,ids)
                  xpr = xpr * y(jsub(j,i,ids))
               end do

               gex = gex + xpr * w(i)

            end do

         end if

      end if

      end

      subroutine gmodel (im,wham)
c---------------------------------------------------------------------
c gmodel - stores ALL solution model parameters in global arrays
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character text*80, abc*1, sname*10

      logical add, wham, zbad, bad

      integer im, nloc, i, j, id, jd, k, m, n, ii, killct, killid(20)

      double precision dx, gcpd, stinc, getstr, zsite(m10,m11), dinc

      external gcpd, zbad, stinc, getstr

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer eos
      common/ cst303 /eos(k10)

      integer nsub,nterm
      double precision acoef
      common/ cst107 /acoef(m10,m11,0:m0),
     *                nterm(m10,m11),nsub(m10,m11,m0)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)
c                                 GLOBAL SOLUTION PARAMETERS:

c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      double precision wgl, wkl, vlar
      common/ cxt2r /wgl(m3,m1,h9),wkl(m16,m17,m18,h9),vlar(m3,m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m14+2)
c                                 special model endmember indexing
      integer jspec
      common/ cxt8 /jspec(h9,m4)

      double precision cp
      common/ cst12 /cp(k5,k10)
c                                 dqf parameters
      integer idqf,indq
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      integer grid
      double precision rid
      common/ cst327 /grid(6,2),rid(5,2)

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)

      character specie*4
      integer jsp, ins
      common/ cxt33 /jsp,ins(nsp),specie(nsp)

      character mname*8
      common/ cst18a /mname(m4)

      integer iam
      common/ cst4 /iam

      integer tnq,tnn,tns
      common/ cxt337 /tnq,tnn,tns

      double precision stch
      common/ cst47 /stch(h9,h4,mst,msp,4)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
c----------------------------------------------------------------------
c                                 auto_refine changes
      if (refine.and.iam.eq.15) then
c                                 check for consistent auto-refine data
         read (n10,'(a)') sname
         if (tname.ne.sname) call error (63,r,i,'GMODEL')

      end if
c                                 read switch to make GALL use MINFXC
      read (tname,'(a)') abc

      if (abc.eq.'X') then
         noder(im) = .true.
         if (iam.lt.3) write (*,*) 'using MINFXC for ',tname
      else
         noder(im) = .false.
      end if 
c                                 initialize autorefine arrays
      stable(im) = .false.
      limit(im) = .false.
      lorch(im) = modres
c                                 initialize compositional distances
      do i = 1, icp
         dcp(i,im) = 0d0
      end do
c                                 check endmember counters:
      if (im.gt.h9) call error (52,dq(1),idqf,'GMODEL')
c                                 check for inconsistent model reformation
      if (kstot+mdep.gt.istot) call error (76,dq(1),idqf,tname)
c                                  out-of-date melt model special cases:
      if (jsmod.eq.24.or.jsmod.eq.25.or.jsmod.eq.27.or.jsmod.eq.28) then
         call error (72,r,i,'the solution model file contains an out-of'
     *        //'-date version of '//tname//' update the model or file')

      end if
c                                 set up simple counters for
c                                 charge balance models
      if (jsmod.eq.20) then

         isimp(1) = 1

         nq = tnq
         ns = tns
         nn = tnn

         nqs = nn + nq + ns
         nqs1 = nqs - 1
         nq1 = nq - 1
         sn = nn + ns
         ns1 = sn - 1
         sn1 = ns + 1
         qn = nn + nq

         ivert(1,1) = nqs

         do i = 1, nqs
            jmsol(i,1) = i
         end do

      end if
c                                 number of dependent + independent - ordered endmembers
c                                 prismatic space
      mstot(im) = istot
c                                 number of independent + ordered endmembers
      nstot(im) = kstot + norder
c                                 number of independent disordered endmembers
      lstot(im) = kstot
c                                 counter for o/d model coordinates
      tstot(im) = nstot(im)
      if (ordmod) tstot(im) = tstot(im) + kstot
c                                 reject bad compositions, only relevant
c                                 for relict equipartition models
      bdx(im) = badx
c                                 override norf if refine_endmembers option is set (default is false)
      if (lopt(39)) norf = .false.
c                                 refine endmembers if norf is false (default is true). since this
c                                 was running with the wrong default value of norf for more than a
c                                 year it seems doubtful that nrf is really useful (lopt(39) IS
c                                 important at least for WADDAH). 11/18.
      nrf(im) = norf
c                                 number of ordered species
      nord(im) = norder
c                                 number of species and multiplicity and
c                                 site ranges
      mcoor(im) = 0
      ncoor(im) = 0
      nsum(im) = 0
c                                 composition space:
c                                 number of polytopes:
      poly(im) = poly(h0)

      if (poly(h0).gt.1) then

         pop1(im) = poly(h0) + 1

      else

         pop1(im) = poly(h0)
         ipvert(1) = mstot(im)

      end if
c                                 polytope and composition 
c                                 variable names.
      do i = 1, pop1(im)
         do j = 1, isimp(i)
            do k = 1, ivert(i,j)
               poname(im,i,j,k) = poname(h0,i,j,k)
            end do
         end do
      end do

      k = 0
      bad = .false.

      do ii = 1, pop1(im)
c                                 number of chemical mixing sites, i.e., simplices, 
c                                 in polytope ii
         istg(im,ii) = isimp(ii)
c                                 total number of simplicies in composition space
         nsum(im) = nsum(im) + isimp(ii)
c                                 index of first and last vertices of the polytope
c                                 in array of vertices
         pvert(im,ii,1) = k + 1
         pvert(im,ii,2) = k + ipvert(ii)
         k = k + ipvert(ii)

         do i = 1, istg(im,ii)
c                                 number of species, i.e., vertices, in each simplex
            ispg(im,ii,i) = ivert(ii,i)
c                                 dimension of each simplex
            ndim(i,ii,im) = ivert(ii,i) - 1

            ncoor(im) = ncoor(im) + ispg(im,ii,i)
            mcoor(im) = mcoor(im) + ndim(i,ii,im)

            do j = 1, ndim(i,ii,im)
c                                 check for old subdivision schemes
               if (pimd(ii,i,j).gt.1) 
     *            call error (62,nopt(13),pimd(ii,i,j),tname)
c                                 allow inc to be either the number of points
c                                 or the resolution
               if (pxnc(ii,i,j).gt.1d0)
     *            pxnc(ii,i,j) = 1d0/pxnc(ii,i,j)
c                                 subdivision override (iopt(13))
               if (iopt(13).eq.1) then
c                                 make linear
                  pimd(ii,i,j) = 0

               else if (iopt(13).eq.2) then
c                                 make all stretch (if not already)
                  if (pimd(ii,i,j).eq.0) pimd(ii,i,j) = 1

               end if

               if (lorch(im)) then 

                  pxnc(ii,i,j) = pxnc(ii,i,j)

               else if (nopt(13).gt.0d0) then

                  if (pimd(ii,i,j).ne.0) then
c                                  set XINC for non-linear subdivision:
c                                  non_linear_switch toggles from default
c                                  i.e., adaptive optimization uses initial_resolution
c                                  and convex uses solution model value
                     if (lopt(38).and.iam.ne.15.or..not.lopt(38).and.
     *                               iam.eq.15) pxnc(ii,i,j) = nopt(13)

                  else

                     pxnc(ii,i,j) = nopt(13)
c                                 and for convexhull: perturb xmn by a per-mil scale increment to
c                                 reduce compositional degeneracies.
                     if (iam.eq.15) pxnc(ii,i,j) = pxnc(ii,i,j)
     *                                    * (1d0 + nopt(15)*float(im-5))

                  end if

               end if
c                                 set stretch parameters according to xmn specified
c                                 in the solution model:
               if (pimd(ii,i,j).ne.0) then

                  dx = pxnc(ii,i,j)
                  if (pxmn(ii,i,j).eq.0d0) pxmn(ii,i,j) = nopt(14)
                  stch(im,ii,i,j,1) = getstr (dx,pxmn(ii,i,j),bad)

                  if (bad) then 
                     call warn (99,dx,i,'GETSTR not converging on'
     *                       //'specified resolution for '//tname//
     *                       ' stretch_factor value will be used')
                     stch(im,ii,i,j,1) = nopt(14)
                  end if

                  stch(im,ii,i,j,2) = stch(im,ii,i,j,1) + 2d0
                  stch(im,ii,i,j,3) = 
     *                              stch(im,ii,i,j,2)/stch(im,ii,i,j,1)
                  stch(im,ii,i,j,4) = dlog(stch(im,ii,i,j,3))
                  pxmn(ii,i,j) = 0d0

               end if
c                                 save solution model values as hard limits
c                                 xmxo/xmxo is reset to 1/0 if hard limits not
c                                 set
               xmno(im,ii,i,j) = pxmn(ii,i,j)
               xmxo(im,ii,i,j) = pxmx(ii,i,j)
c                                 true hard limit record
               xmnh(im,ii,i,j) = pxmn(ii,i,j)
               xmxh(im,ii,i,j) = pxmx(ii,i,j)

               if (refine.and.iam.eq.15) then
c                                 new values from autorefine file
                  read (n10,*) pxmn(ii,i,j),pxmx(ii,i,j)

                  if (icopt.lt.4) then
c                                 set slop to the initial spacing
                     dinc = pxnc(ii,i,j)

                  else
c                                 adaptive minimization:
c                                 set slop to the final compositional
c                                 resolution of the exploratory stage
                     dinc = rid(4,1)

                  end if
c                                 widen the range by the exploratory resolution
                  if (pimd(ii,i,j).eq.0) then
                     pxmx(ii,i,j) = pxmx(ii,i,j) + dinc
                     pxmn(ii,i,j) = pxmn(ii,i,j) - dinc
                  else
                     pxmn(ii,i,j) = 
     *                             stinc (pxmn(ii,i,j),-dinc,im,ii,i,j)
                     pxmx(ii,i,j) = 
     *                             stinc (pxmx(ii,i,j),dinc,im,ii,i,j)
                  end if

                  if (pxmx(ii,i,j).gt.1d0) pxmx(ii,i,j) = 1d0
                  if (pxmn(ii,i,j).lt.0d0) pxmn(ii,i,j) = 0d0

                  if (lopt(3)) then
c                                 hard_limit test
                     if (pxmx(ii,i,j).gt.xmxo(im,ii,i,j)) 
     *                                  pxmx(ii,i,j) = xmxo(im,ii,i,j)
                     if (pxmn(ii,i,j).lt.xmno(im,1,i,j)) 
     *                                  pxmn(ii,i,j) = xmno(im,ii,i,j)
                  end if

                  pxnc(ii,i,j) = pxnc(ii,i,j)/nopt(17)

               end if

               imdg(j,i,ii,im) = pimd(ii,i,j)
               xmng(im,ii,i,j) = pxmn(ii,i,j)
               xmxg(im,ii,i,j) = pxmx(ii,i,j)
               xncg(im,ii,i,j) = pxnc(ii,i,j)

            end do

         end do
c                                 initialize high/low ranges
         do i = 1, istg(im,ii)
            do j = 1, ispg(im,ii,i)

               xlo(j,i,ii,im) = 1d0
               xhi(j,i,ii,im) = 0d0

            end do
         end do
      end do
c                                 -------------------------------------
c                                 classify the model
      ksmod(im) = jsmod
c                                 initialize/set type flags, presently no provision for
c                                 bw summation
      llaar(im) = .false.
      lexces(im) = .false.
      lorder(im) = .false.
      lrecip(im) = .false.
      specil(im) = .false.
      simple(im) = .false.
      extyp(im) = xtyp

      if (jsmod.eq.2.or.(jsmod.eq.688.or.jsmod.eq.7).and.
     *                 .not.ordmod) simple(im) = .true.

      if (jsmod.eq.30.or.jsmod.eq.31) specil(im) = .true.

c                                 this looks like bad news, for laar/recip
c                                 or laar/order, but appears to be overridden
c                                 by use of logical classification variables,
c                                 in which case, why is it here????
      if (laar.and.iterm.eq.0) laar = .false.

      if (iterm.gt.0) then
         lexces(im) = .true.
         if (laar) then
            llaar(im) = .true.
            extyp(im) = 2
         end if
      end if

      if (ordmod) lorder(im) = .true.
c                                 the ksmod(im) test is made because
c                                 reform may dump the dependent endmembers
c                                 setting depend = .false., while retaining
c                                 a dummy site with no mixing. reform should
c                                 be redone to truly reformulate multiple
c                                 models to single site models.
c                                 a non-reciprocal model (ksmod=5) with
c                                 dependent endmembers is also classified
c                                 as lrecip.
      if (recip.or.depmod) lrecip(im) = .true.
c                                 -------------------------------------
c                                 save the excess terms
      jterm(im) = iterm

      do i = 1, iterm

         rko(i,im) = rkord(i)

         if (xtyp.eq.0) then
c                                 van-laar implementation
c                                 probably assumes quadratic...
            if (laar.and.rkord(i).gt.2) 
     *         call error (999,dx,800,'RMODEL')
c                                 arbitrary expansion
            do j = 1, m3
               wgl(j,i,im) = wg(i,j)
            end do

         else
c                                 redlich-kistler
            do k = 1, rkord(i)
               do j = 1, m16
                  wkl(j,k,i,im) = wk(j,k,i)
               end do
            end do

         end if

         do j = 1, rkord(i)
c                                 isub points to the position in the list
c                                 of endmembers potentially including dependent
c                                 species. use iy2p to convert to independent
c                                 endmember pointers.
            jsub(j,i,im) = iy2p(isub(i,j))

            if (kdsol(isub(i,j)).eq.-2) then

               call error (72,r,i,'dependent endmember '
     *              //mname(iorig(isub(i,j)))//' in solution '
     *              //tname//'appears in an excess term.')

            end if

         end do

      end do
c                                 save global copy of kdsol
      ldsol(1:mstot(im),im) = kdsol(1:mstot(im))
c                                 insp points to the original position
c                                 of endmember i in the solution model input:
      knsp(1:mstot(im),im) = insp(1:mstot(im))
c                                 kmsol points to the species on the j'th site
c                                 of the i'th endmember, used for the xtoy
c                                 conversion
       do ii = 1, poly(h0)

          kmsol(im,pvert(im,ii,1):pvert(im,ii,2),1:istg(im,ii)) 
     *    = jmsol(pvert(im,ii,1):pvert(im,ii,2),1:istg(im,ii))

       end do
c                                 ----------------------------------------------
c                                 configurational entropy models

c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      nloc = 0

      if (zuffix(h0).eq.'none') zuffix(h0) = ' '
      zuffix(im) = zuffix(h0)

      do i = 1, msite(h0)
c                                 eliminate sites with 1 species
         if (zmult(h0,i).gt.nopt(50)) then
c                                 non-temkin
            if (zsp(h0,i).lt.nopt(50)) then
c                                 pad zuffix with the remaining species
               if (tzmult(h0,i).eq.1d0) then 
                  znames(h0,i,2) = ' '
               else 
                  write (znames(h0,i,2),'(i2)') idint(tzmult(h0,i))
               end if

               text = znames(h0,i,0)//'['//znames(h0,i,1)
     *                                   //znames(h0,i,2)//']'

               call unblnk (text)
               call mertxt (text,text,zuffix(h0),0)
               zuffix(h0) = text
               zuffix(im) = text
               cycle

            end if
         else
            if (zsp(h0,i).eq.1d0) then
c                                 pad zuffix with the remaining species
               text = znames(h0,i,0)//'['//znames(h0,i,1)//']'
               call unblnk (text)
               call mertxt (text,text,zuffix(h0),0)
               zuffix(h0) = text
               zuffix(im) = text
               cycle

            end if
         end if

         nloc = nloc + 1
c                                 # of species, and site r*multiplicty.
         zmult(im,nloc) = r*zmult(h0,i)
         tzmult(im,nloc) = tzmult(h0,i)
         znames(im,nloc,0) = znames(h0,i,0)
         zsp1(im,nloc) = zsp1(h0,i)
         zsp(im,nloc) = zsp(h0,i)
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, zsp1(h0,i)
c                                 # of terms in the
c                                 site fraction function and a0.
            lterm(j,nloc,im) = nterm(i,j)
            dcoef(0,j,nloc,im) = acoef(i,j,0)
            znames(im,nloc,j) = znames(h0,i,j)
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 term coefficient amd species index:
               dcoef(k,j,nloc,im) = acoef(i,j,k)
               ksub(k,j,nloc,im) = iy2p(nsub(i,j,k))

               if (kdsol(nsub(i,j,k)).eq.-2) then

                  call error (72,r,k,'dependent endmember '
     *              //mname(iorig(nsub(i,j,k)))//' in solution '
     *              //tname//' appears in a site fraction expression.')

               end if

            end do
         end do
      end do
c                                 number of distinct identisites for entropy
      msite(im) = nloc
c                                 --------------------------------------
c                                 van laar volumes, and pointers for "in" endmembers
      do i = 1, nstot(im)
c                                 if the solution is van laar save
c                                 the "volume" function.
         if (laar) vlar(1:m3,i,im) = vlaar(1:m3,jnsp(i))
c                                 initialize scoef's to zero for config
c                                 entropy calculation (done by snorm).
         scoef(i,im) = 0d0

      end do
c                                 -------------------------------------
      if (depmod) then
c                                 march, 2017: deleted y2p4z routine that converted
c                                 z(y) expressions to z(p), i.e., by simply
c                                 eliminating dependent endmembers.

c                                 save y -> p array
         ndep(im) = mdep

         do i = 1, nstot(im)
            y2pg(1:mdep,i,im) = y2p(i,1:mdep)
         end do

         do j = 1, mdep

            y(1:mstot(im)) = 0d0

            y(knsp(lstot(im)+j,im)) = 1d0

            call y2p0 (im)
c                                 check for invalid site fractions, this is only necessary
c                                 for H&P models that assume equipartition (which is not
c                                 implemented).
            if (zbad(pa,im,zsite,tname,.true.,
     *                     mname(iorig(knsp(lstot(im)+j,im))))) then

               if (iam.lt.3.or.iam.eq.4.or.iam.eq.15)
     *            call warn (59,y(1),i,
     *            mname(iorig(knsp(lstot(im)+j,im)))
     *            //' in solution model '//tname)

               if (stck) call error (78,y(1),i,tname)

            end if

         end do

      else

         ndep(im) = 0

      end if
c                                 -------------------------------------
c                                 relict equipartion warning:
      if (.not.stck.and..not.refine.and.iam.lt.3) then

         call warn (17,r,i,tname)

         if (lopt(56)) call wrnstp


      end if
c                                 -------------------------------------
c                                 dqf parameters
      jdqf(im) = idqf

      do i = 1, idqf

c                                 a dqf to an ordered species
c                                 or a dependent endmember
         if (kdsol(indq(i)).lt.0)
     *                       call error (227,r,indq(i),fname(im))
c                                 shift pointer from y array to p array
         jndq(i,im) = iy2p(indq(i))

         do j = 1, m3
            dqfg(j,i,im) = dqf(j,i)
         end do

      end do
c                                 -------------------------------------
c                                 if msite(h0) ne 0 get "normalization" constants (endmember
c                                 configurational entropies) for entropy model:
      if (msite(im).ne.0) call snorm0 (im,tname)
c                                 -------------------------------------
c                                 organize O/D model parameters
      call setord (im)

      if (.not.equimo(im)) then
c                                 non-equimolar restrictions:
         if (lrecip(im)) call error (72,r,i,'prismatic composition spa'/
     *        /'ce not anticipated for non-equimolar ordering: '//tname)
         if (ksmod(im).ne.688) call error (72,r,i,'non-equimolar order'/
     *               /'ing only allowed for 688 format solution models')

      end if
c                                 ----------------------------------------------
c                                 models with special endmember indexing:
      if (jsmod.eq.0) then
c                                 fluid eos, make pointer to co2
         do i = 1, 2
            id = kdsol(insp(i))
            if (cp(2,id).ne.0d0) then
               jspec(im,1) = i
               exit
            end if
         end do

      else if (jsmod.eq.39.or.jsmod.eq.20) then
c                                 generic hybrid fluid EoS, convert endmember EoS
c                                 flag to fluid species indices
         if (jsmod.eq.20) then
            j = ns
         else
            j = istot
         end if

         do i = 1, j
            k = eos(kdsol(insp(i)))
            if (k.gt.200) then
               write (*,1000) tname, names(kdsol(insp(i)))

1000  format (/,'**error ver888** a special component endmember cannot',
     *' be used in solution model ',a,/,
     *'This error can be corrected by any of the following actions:',/,
     *'1) Setting the GFSM option to True in perplex_option.dat',/,
     *'2) Deleting ',a,' from the special component section of the ',
     *'thermodynamic data file',/,
     *'3) Choosing a different solution model',/)

               call errpau

            else if (k.gt.100) then
               k = k - 100
            else
               call error (72,0d0,0,'invalid endmember EoS '//
     *               'specification in solution model '//tname)
            end if

            jspec(im,i) = k

         end do

      else
c                                 save original endmember indexes for hard-wired
c                                 solution models
         do i = 1, istot
            jspec(im,i) = iorig(i)
         end do

      end if
c                                  set fluid flags for non-special case melts
      if (lname(im).eq.'liquid'.or.lname(im).eq.'fluid') then
c                                  set ifp for t_melt and melt_is_fluid
         do i = 1, lstot(im)
c                                 of endmember i in the solution model input:
            if (lname(im).eq.'liquid') then
               ifp(kdsol(knsp(i,im))) = -1
            else
               ifp(kdsol(knsp(i,im))) = 1
            end if

         end do

      end if

      if (istot+norder.gt.m4) call error (39,0d0,m4,'INPUT9')
c                                 smod is used for all shear modulus calculations, it is
c                                 set to false if a) the EoS for an endmembers doesn't 
c                                 provide for its shear modulus or b) the endmember lacks
c                                 an explicit shear modulus function.
      smod(im) = .true.
c                                 pmod is used only for the explicit_bulk_modulus option
c                                 it is set to true if explicit functions for the bulk
c                                 modulus are present for all endmembers 
      pmod(im) = .true.

      killct = 0
c                                 set fluid flag from the full name of the
c                                 phase:
      if (lname(im).eq.'fluid') then
         fp(im) = .true.
      else
         fp(im) = .false.
      end if
c                                 classify liquid model as fluid/not fluid
c                                 according to the melt_is_fluid option, this
c                                 is only relevant for WERAMI and MEEMUM
      if (lname(im).eq.'liquid'.and.(iam.eq.3.or.iam.eq.2).and.lopt(6))
     *                                                 fp(im) = .true.

      do i = 1, lstot(im)
c                                 pointer to endmember
         id = kdsol(insp(i))
c                                 save the total mols of the endmember in a solution
c                                 specific array, this is done so the ordered 
c                                 endmembers do not need to be treated as a special case
         endt(im,i) = ctot(id)
         endc(im,i,1:icomp) = cp(1:icomp,id)
c                                 figure out the compositional distance between
c                                 the endmembers, this is used to scale the solvus
c                                 tolerance
         do j = i+1, lstot(im)

            jd = kdsol(insp(j))
            if (ctot(id)*ctot(jd).eq.0d0) cycle
c                                 switched from molar to mole fraction in 688
            do k = 1, icp

               dx = dabs(cp(k,id)/ctot(id) - cp(k,jd)/ctot(jd))

               if (dx.lt.nopt(50)) then
                  cycle
               else if (dcp(k,im).lt.dx) then
                  dcp(k,im) = dx
               end if

            end do

         end do
c                                 set ifp for models with an
c                                 endmember identified as "fluid"
         if (ifp(id).gt.0) fp(im) = .true.
c
         jend(im,2+i) = id
c                                 set shear/bulk moduli flags
         if (iemod(id).eq.0) then 
            smod(im) = .false.
            pmod(im) = .false.
         else if (iemod(id).eq.1) then 
            pmod(im) = .false.
         else if (iemod(id).eq.3) then 
            smod(im) = .false.
         end if
c                                 look for endmembers to be killed
         if (iend(insp(i)).ne.2) cycle

         add = .true.

         do j = 1, killct

            if (killid(j).eq.id) then
               add = .false.
               exit
            end if
         end do

         if (add) then

            killct = killct + 1
            if (killct.gt.10) call error (999,wg(1,1),killct,tname)
            killid(killct) = id

         end if

      end do
c                                 compute and save the total moles for the 
c                                 ordered endmembers
      if (lstot(im).lt.nstot(im)) then 
c                                  if o/d add the ordered endmembers
         do m = 1, norder
c                                  the ordered endmember compositions are 
c                                  a stoichiometric combination of the 
c                                  disordered endmembers
            j = lstot(im) + m

            endt(im,j) = 0d0
            endc(im,j,:) = 0d0

            do n = 1, nr(m)

               jd = jend(im,2+ideps(n,m,im))

               endt(im,j) = endt(im,j) + depnu(n,m) * ctot(jd)

               do i = 1, icomp
                  endc(im,j,i) = endc(im,j,i) + depnu(n,m) * cp(i,jd)
               end do

            end do

         end do

      end if
c                                 set pmod to false if explicit_bulk_modulus is not T
      if (.not.lopt(17)) pmod(im) = .false.

      if (.not.lopt(3)) then
c                                 hard limits are off, set limits to 0/1
         do ii = 1, pop1(im)
            do i = 1, isimp(ii)
               do j = 1, ivert(ii,i) - 1

                  xmxo(im,ii,i,j) = 1d0
                  xmno(im,ii,i,j) = 0d0

               end do 
            end do
         end do
      end if

      if (ksmod(im).eq.0. or.ksmod(im).eq.20.or.ksmod(im).eq.39.or.
     *    ksmod(im).eq.40.or.ksmod(im).eq.41) then
c                                 routines that invoke fluid EoS, as currently
c                                 configured these will set ins/isp arrays only
c                                 once. therefore some parameters and indices
c                                 can be saved in simple arrays
         call setsol (im,wham)

      end if
c                                 set independent species names and counters for output
c                                 special cases first:
      if (ksmod(im).eq.0.or.ksmod(im).eq.40.or.ksmod(im).eq.41) then
         spct(im) = jsp

         do i = 1, jsp
            spnams(i,im) = specie(ins(i))
         end do

      else

         spct(im) = nstot(im)
c                                 independent disordered species
         do i = 1, lstot(im)
            spnams(i,im) = mname(iorig(knsp(i,im)))
         end do
c                                 ordered species
         do i = lstot(im)+1, nstot(im)
            spnams(i,im) = mname(iorig(ndep(im)+i))
         end do

         if (ksmod(im).eq.29.or.ksmod(im).eq.32) then
c                                 BCC Fe-Si Lacaze and Sundman (29)
c                                 BCC Fe-Cr Andersson and Sundman (32)
            spct(im) = 4
            spnams(3,im) = 'osp'
            spnams(4,im) = 'asp'

         end if

      end if
c                                 by default assume simplicial models are bounded, this
c                                 can be overridden (unbd) by the unbounded_composition solution 
c                                 model keyword.
      if ((.not.unbd.and.lstot(im).eq.nstot(im).and.
     *     lstot(im).eq.mstot(im)).or.
     *    (.not.unbd.and.lorder(im).and..not.equimo(im))) then

          boundd(im) = .true.

      else 

          boundd(im) = .false.

      end if
c                                 make transformation matrices, p' is 
c                                 the nstot-1 independent p variables.
c                                 y2x
      call makayx (im)
c                                 p'2z 
      call makapz (im)
c                                 y2z, ayz must be called after makapz 
      call makayz (im)
c                                 y2c
      call makayc (im)
c                                 p'2c
      call makapc (im)
c                                 set derivatives for minfrc
      call setder (im,tname)

      end

      subroutine y2p0 (id)
c-----------------------------------------------------------------------
c y2p0 converts the y array of disordered dependent and independent
c species abundance to the pa array of the independent (ordered and
c disordered) species. pa is copied to p0a and converted to pp by 
c the call to makepp.

c for simplicial solutions the y and pa arrays are identical.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,k,l

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
c                                 initialize ordered species
      pa(lstot(id)+1:nstot(id)) = 0d0

      do k = 1, nstot(id)
c                                 initialize the independent species
c                                 other than the ordered species
         if (k.le.lstot(id)) pa(k) = y(knsp(k,id))
c                                 convert the dependent species to
c                                 idependent species
         do l = 1, ndep(id)
            pa(k) = pa(k) + y2pg(l,k,id) * y(knsp(lstot(id)+l,id))
         end do

      end do
c                                 convert the ordered species to 
c                                 the stoichiometric equivalent 
c                                 amounts of disordered species.
      call makepp (id)

      end

      subroutine makepp (id)
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,k,l,ind

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
c                                 order-disorder, need initial speciation
c                                 fully disordered for static compositions
c                                 but may be partially ordered for dynamic
c                                 compositions. 
      p0a(1:nstot(id)) = pa(1:nstot(id))
      pp(1:nstot(id)) = pa(1:nstot(id))

      do k = 1, nord(id)
         do l = 1, nrct(k,id)
            ind = ideps(l,k,id)
            pp(ind) = pp(ind) - dydy(ind,k,id) * pp(lstot(id)+k)
         end do
      end do
c                                 zero ordered pp's
      pp(lstot(id) + 1:nstot(id)) = 0d0

      end

      subroutine specis (g,id)
c----------------------------------------------------------------------
c subroutine to compute speciation of a solution with initial speciation
c p0a. the speciated composition is returned in array pa.
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      logical minfx, error

      double precision g, oldg, oldp(m14), g0, gordp0

      external gordp0

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      g0 = gordp0 (id)
      g = g0

      minfx = .false.

      if (iopt(37).lt.0) then
c                                 minfxc has been made the default solver:
         call minfxc (g,id,minfx)

      else if (nord(id).eq.1) then
c                                 as most models are single species, use 
c                                 special routines to avoid overhead.
         if (equimo(id)) then
c                                 initialize limit expressions
            call p0limt (id)

            call speci1 (g,id,1)

         else
c                                 assumes non-equimolar speciation is not
c                                 subject to site fraction restrictions
            call gpmlt1 (g,1,id,error)

         end if

      else

         if (equimo(id)) then
c                                 initialize limit expressions
            call p0limt (id)

            call speci2 (g,id,minfx)

         else

            call gpmelt (g,id,minfx)

         end if

      end if
c                                 if order_check (lopt(62)) or a routine
c                                 has set minfx = true, try bail out solution
      if (minfx.and.iopt(37).ne.5.or.lopt(62)) then
c                                 order_check option
         oldp(1:nstot(id)) = pa(1:nstot(id))
c                                 this is necessary for pinc0
         pa(1:nstot(id)) = p0a(1:nstot(id))

         oldg = g

         call minfxc (g,id,.false.)

         if (oldg-g.lt.-nopt(53)) then 
c                                 restore the old result
               g = oldg
               pa(1:nstot(id)) = oldp(1:nstot(id))

c              if (.not.equimo(id)) then 

c                  write (*,*) 'not worked',pa(13:14)

c               end if

c        else if (.not.equimo(id)) then 

c           write (*,*) 'worked',pa(13:14)

         end if

       end if
c                                 compare g for initial and eq speciation
c                                 return the lower
      if (g.gt.g0) then 
         g = g0
         pa = p0a
      end if

      end

      double precision function gordp0 (id)
c----------------------------------------------------------------------
c check that the current proportions of an o/d model are more stable
c than the initial (usually disordered) proportions, swap if not. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id

      double precision g, omega, gex

      external omega, gex

      double precision enth
      common/ cxt35 /enth(j3)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      g =  gex(id,p0a) - t*omega(id,p0a)

      do i = 1, nord(id)
         g = g + p0a(lstot(id)+i)*enth(i)
      end do

      gordp0 = g

      end

      subroutine oenth (id)
c----------------------------------------------------------------------
c subroutine to the enthalpy of ordering for speciation models
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id

      double precision enth
      common/ cxt35 /enth(j3)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      do i = 1, nord(id)
         enth(i) = deph(1,i,id) + t * deph(2,i,id) + p * deph(3,i,id)
      end do

      end

      subroutine gderiv (id,g,dg,minfx,error)
c----------------------------------------------------------------------
c subroutine to compute the g of a solution (id) and it's 1st and 2nd
c derivatives with respect to the fractions of the nord(id) ordered
c species. the formulation assumes atomic site fractions are linear
c functions of the ordered species fractions (p's) and that the
c excess function is second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical minfx, error

      integer i,k,l,i1,i2,i3,id,ipvt(j3)

      double precision g,dg(*),t,s,ds(j3),d2s(j3,j3),d2g(j3,j3)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision p,tk,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,tk,xc,u1,u2,tr,pr,r,ps

      double precision enth
      common/ cxt35 /enth(j3)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------
c                                 initialize
      g = 0d0
      norder = nord(id)
      dg(1:norder) = 0d0
      d2g(1:norder,1:norder) = 0d0

      if (lexces(id)) then

         do i = 1, jterm(id)
c                                 rather than write expensive
c                                 general code only the 2nd and 3rd
c                                 order cases are considered here:
            if (rko(i,id).eq.2) then

               i1 = jsub(1,i,id)
               i2 = jsub(2,i,id)

               g = g + w(i) * pa(i1) * pa(i2)

               do k = 1, norder

                  if (.not.pin(k)) cycle

                  dg(k) = dg(k) + w(i) * (pa(i1)*dydy(i2,k,id)
     *                                  + pa(i2)*dydy(i1,k,id))

                  do l = k, norder

                     d2g(l,k) = d2g(l,k) + w(i) * dppp(l,k,i,id)

                  end do

               end do

            else if (rko(i,id).eq.3) then

               i1 = jsub(1,i,id)
               i2 = jsub(2,i,id)
               i3 = jsub(3,i,id)

               g = g + w(i) * pa(i1) * pa(i2) * pa(i2)

               do k = 1, norder

                  if (.not.pin(k)) cycle

                  dg(k) = dg(k) + w(i) * (
     *                                     pa(i1)*pa(i2)*dydy(i3,k,id)
     *                                   + pa(i1)*pa(i3)*dydy(i2,k,id)
     *                                   + pa(i2)*pa(i3)*dydy(i1,k,id) )

                  do l = k, norder

                     d2g(l,k) = d2g(l,k) + w(i) * (
     *                            pa(i1)*(dydy(i2,l,id)*dydy(i3,k,id) +
     *                                    dydy(i2,k,id)*dydy(i3,l,id))
     *                          + pa(i2)*(dydy(i1,l,id)*dydy(i3,k,id) +
     *                                    dydy(i1,k,id)*dydy(i3,l,id))
     *                          + pa(i3)*(dydy(i1,l,id)*dydy(i2,k,id) +
     *                                    dydy(i1,k,id)*dydy(i2,l,id)) )

                  end do

               end do

            else

               call errdbg ('o > 3 gderiv')

            end if

         end do

         if (llaar(id)) then
c                                 so far this is unessecary because
c                                 t does not vary in h&p ordering models.
            t = 0d0
c                                 h&p van laar
            do i = 1, nstot(id)
               t = t + alpha(i)* pa(i)
            end do
c                                 coming out of this loop g, dg, and
c                                 d2g  are not the complete functions
c                                 because of the "tphi" term in the
c                                 van laar.

            do k = 1, norder

               if (.not.pin(k)) cycle
c                                 convert dg and d2g to the full derivative
               dg(k) = (dg(k) - g*dt(k)/t)/t
               do l = k, norder
                  d2g(l,k) = (d2g(l,k) - 2d0*dt(k)*dg(k))/t
               end do
            end do
c                                 and the full excess energy
            g = g/t

         end if

      end if
c                                 get the delta configurational entropy and derivatives
      call sderiv (id,s,ds,d2s)

      do k = 1, norder

         g = g + enth(k) * pa(lstot(id)+k)

         if (.not.pin(k)) cycle
c                                 dg is the negative of the differential of g
c                                 with respect to the kth species.
         dg(k) = -(enth(k) + dg(k) - tk*ds(k))

         do l = k, norder
            d2g(l,k) = d2g(l,k) - tk*d2s(l,k)
         end do

      end do
c                                 determininats, to check for a saddle point
c      if (norder.eq.2) then
c         detg = d2g(1,1)*d2g(2,2)-d2g(2,1)**2
c      else
c         detg = d2g(1,1)*(d2g(2,2)*d2g(3,3)-d2g(3,2)**2)
c     *        - d2g(2,2)*d2g(1,3)**2
c     *        + 2d0*d2g(2,1)*d2g(3,2)*d2g(3,1)-d2g(2,1)**2*d2g(3,3)
c      end if

      g = g - tk*s
c                                 if minfx just return with the gradient
      if (minfx) then

         dg(1:norder) = -dg(1:norder)

         return

      end if
c                          compute the newton-raphson increments:
      do k = 1, norder

         if (pin(k)) then
c                          flesh out the hessian
            do l = 1, k-1
               d2g(l,k) = d2g(k,l)
            end do

         else

            dg(k) = 1d0
            d2g(k,k) = 1d0

            do l = 1, norder
               if (l.eq.k) cycle
               d2g(l,k) = 0d0
               d2g(k,l) = 0d0
            end do

         end if
      end do
c                                 get the newton-raphson increments:
c                                 this is a general factorization routine, should
c                                 exploit that d2g is symmetric.
      call factor (d2g,j3,norder,ipvt,error)
c                                 solve for the increments by back-substitution,
c                                 this routine is also not efficient and should
c                                 be re written.
      if (.not.error) call subst (d2g,j3,ipvt,norder,dg,error)
c                                 substitute replaces the values of dg with the
c                                 newton-raphson increments for the ordered species
c                                 compositions.
      end

      subroutine sderiv (id,s,dsy,dsyy)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a
c solution with respect to the proportion of a dependent species.
c if maxs, then sderiv is being called by gsol4/minfxc for negentropy
c minimization/p2yx inversion the jacobian is unnecessary and no corrections
c are made for endmember configurational entropy.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,id

      double precision zt,qdzdy,s,dsy(*),dsyy(j3,*),q,zl,
     *                 z(m11,m10),s0,ztemp,zlnz
c                                 working arrays
      double precision zz, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------
      s = 0d0
      dsy(1:nord(id)) = 0d0
      dsyy(1:nord(id),1:nord(id)) = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         s0 = zt
c                                 get site fractions
         do j = 1, zsp(id,i)

            ztemp = dcoef(0,j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               ztemp = ztemp + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do

            call ckzlnz (ztemp,s0)

            zt = zt + ztemp

            z(j,i) = ztemp

         end do

         ztemp = 1d0 - zt

         call ckzlnz (ztemp,s0)

         z(j,i) = ztemp
         s = s - zmult(id,i) * s0

      end do
c                                 evaluate derivatives:
      do i = 1, msite(id)

         q = zmult(id,i)

         do j = 1, zsp(id,i) + 1

            zl = z(j,i)

            if (zl.gt.0d0) then 
               zlnz = 1d0 + dlog(zl)
            else
               zl = nopt(50)
               zlnz = nopt(54)
            end if

            do k = 1, nord(id)
c                                 skip species not in the model
               if (.not.pin(k)) cycle
c                                 sdzdp is (dz(i,j)/dp(k))
               qdzdy = q * sdzdp(k,j,i,id)

               if (qdzdy.eq.0d0) cycle
c                                 the first derivative is
               dsy(k) = dsy(k) - qdzdy * zlnz
c                                 and the jacobians are
               do l = k, nord(id)

                  if (.not.pin(l)) cycle

                  dsyy(l,k) = dsyy(l,k) - qdzdy * sdzdp(l,j,i,id) / zl

               end do

            end do

         end do

      end do

      if (maxs) then

         s = -s
         dsy(1:nord(id)) = -dsy(1:nord(id))
         return

      end if
c                                 endmember corrections
      do i = 1, nstot(id)

         s = s - pa(i) * scoef(i,id)

         do k = 1, nord(id)
           dsy(k) = dsy(k) - dydy(i,k,id) * scoef(i,id)
         end do

      end do

      end

      subroutine speci0 (g,h,w,n,fac,c0,f)
c----------------------------------------------------------------------
c subroutine to solve speciation of a 0-d solution with 1 ordering parameter
c by halving. assumes an ordered species in which A is on 1 site and B is on
c n sites, and a disordered state in which A and B are distributed over all
c n+1 sites.

c    h   - is the enthalpy of complete disordering
c    w   - is the interaction energy
c    fac - is an empirical correction to the entropy, supposedly accounting for SRO.
c    g   - is the change in G for the stable speciation relative to the ordered state.
c    y   - is the fraction of the ordered species

c                                                  JADC, Aug 29, 2010.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision g,h,w,sign,dy,odg,ndg,n,f,y,dgdy,rt,c0,c1,c2,fac

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
c                                 check ordered state
      y = 1d0 - nopt(50)
      rt = r*t*fac

      odg = dgdy(h,w,n,f,y,rt)
c                                 if dgdy < 0 must be fully ordered
c                                 (not really, non-zero w could make
c                                 a zero at intermediate y).
      if (odg.lt.0d0) then
c                                 fully ordered (y=1)
         y = 1d0

      else
c                                 initialize at halfway point
         dy = -0.5d0
c                                 iteration loop:
         do

            y = y + dy
            if (y.le.0d0) y = nopt(50)

            ndg = dgdy(h,w,n,f,y,rt)

            sign = odg*ndg

            if (sign.lt.0d0) then
c                                 crossed the zero, flip the search
               odg = ndg
               dy = -dy/2d0

            else if (dabs(dy/(1d0+y)).lt.nopt(50)) then
c                                 refined to tolerance
               exit

            else if (y.le.nopt(50)) then
c                                 fully disordered, y=0, c1 = c2
               y = 0d0
               exit

            end if

         end do

      end if

      c1 = (n+y)/c0

      if (c1.lt.nopt(56).and.c1.gt.nopt(50)) then
         g = rt*n*(c1*dlog(c1)+(1d0-c1)*dlog(1d0-c1))
      else
         g = 0d0
      end if

      c2 = (1d0-y)*n/c0

      if (c2.lt.nopt(56).and.c2.gt.nopt(50))
     *   g = g + rt*(c2*dlog(c2) + (1d0-c2)*dlog(1d0-c2))

      g = g + (1d0-y)*( w*y + h)

      end

      double precision function dgdy (h,w,n,f,y,rt)
c----------------------------------------------------------------------
c function to compute dg/dy for subroutine speci0
c----------------------------------------------------------------------
      implicit none

      double precision h,w,n,f,y,rt

      dgdy = (1d0-2d0*y)*w - h
     *       - rt*f*dlog(n*(1d0-y)**2/(n+y)/(1d0+n*y))

      end

      subroutine speci1 (g,id,k)
c----------------------------------------------------------------------
c subroutine to speciation of a solution with a single ordering parameter
c composition is returned in array pa.

c    k  - the ordered species
c    id - the solution.
c    g  - the change in G for the stable speciation relative to a mechanical
c         mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id, jd, k, itic, ind(m14), nr

      logical error, done, maxok, minok, usemax

      double precision g, pmax, pmin, dp, gord, dy(m14), gold, xdp,
     *                    gmax, gmin, wt

      external gord

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c----------------------------------------------------------------------
c                                 number of reactants to form ordered species k
      nr = nrct(k,id)

      do i = 1, nr
c                                 dependent disordered species
         ind(i) = ideps(i,k,id)
c                                 stoichiometric coefficients
         dy(i) = dydy(ind(i),k,id)

      end do 

      jd = lstot(id) + k

      error = .false.
c                                 starting point
      call plimit (pmin,pmax,k,id)
c                                 necessary?
      pin(k) = .true.
c                                 a composition for which no O/D 
c                                 is possible
      if (pmax-pmin.lt.nopt(50)) return
c                                 to avoid singularity set the initial
c                                 composition to the max - nopt(50), at this
c                                 condition the first derivative < 0,
c                                 and the second derivative > 0 (otherwise
c                                 the root must lie at p > pmax - nopt(50).
      pmax = pmax - nopt(50)
      pmin = pmin + nopt(50)
c                                 test if concave, convex, or mixed
c                                 at pmax
      call pincs (pmax - p0a(jd),dy,ind,jd,nr)
      call gderi1 (k,id,dp,gmax)

      if (dp.lt.0d0) then 
         maxok = .true.
      else
         maxok = .false.
      end if 
c                                 at pmin
      call pincs (pmin - p0a(jd),dy,ind,jd,nr)
      call gderi1 (k,id,dp,gmin)

      if (dp.gt.0d0) then 
         minok = .true.
      else 
         minok = .false.
      end if
c                                 decide where to start:
      usemax = .false.

      if (maxok.and.minok) then 
c                                 convex at edges
         if (gmax.le.gmin) usemax = .true.

      else if (maxok) then
c                                 convex at max
          usemax = .true.

      else if (.not.minok) then

         error = .true.

      end if

      if (.not.error) then 
c                                 initialize
         wt = 0.01

         if (usemax) then
            dp = (1d0-wt)*pmax + wt*pmin - p0a(jd)
         else
            dp = (1d0-wt)*pmin + wt*pmax - p0a(jd)
         end if
c                                 set starting point
         call pincs (dp,dy,ind,jd,nr)
c                                 iteration counter
         itic = 0
         gold = 0
         xdp = 0d0
c                                 newton raphson iteration
         do

            call gderi1 (k,id,dp,g)

            call pcheck (pa(jd),pmin,pmax,dp,done)
c                                 done means the search hit a limit
c                                 or dp < tolerance.
            if (done.or.dabs((gold-g)/(1d0+dabs(g))).lt.nopt(50)) then

               goodc(1) = goodc(1) + 1d0
               goodc(2) = goodc(2) + dfloat(itic)
c                                 use the last increment
               call pincs (pa(jd)-p0a(jd),dy,ind,jd,nr)

               exit

            else if (dp.eq.xdp) then 

               write (*,*) 'wroink! oscillating?',g-gold,id,itic

            else
c                                 apply the increment
               call pincs (pa(jd)-p0a(jd),dy,ind,jd,nr)

               if (itic.gt.iopt(21)) then
c                                 failed to converge. exit
c                 write (*,*) 'wroink2! failed, ',
c    *                        'increase speciation_max_it?',g-gold,id
                  error = .true.
                  badc(1) = badc(1) + 1d0
                  goodc(2) = goodc(2) + dfloat(itic)

                  exit

               end if

               xdp = dp
               gold = g
               itic = itic + 1

            end if

         end do

      end if
c                                 didn't converge or couldn't
c                                 find a starting point, set
c                                 ordered speciation, specis will
c                                 compare this the disordered case.
      if (error) then
c                                 concave, bail
         if (gmax.le.gmin) then 
c                                 ordered
            g = gmax
            call pincs (pmax-p0a(jd),dy,ind,jd,nr)

         else
c                                 anti-ordered
            g = gmin
            call pincs (pmin-p0a(jd),dy,ind,jd,nr)

         end if

      end if

      end

      subroutine speci2 (g,id,minfx)
c----------------------------------------------------------------------
c subroutine to multiple speciation of a solution with disordered composition
c p0a. the speciated composition is returned in array pa.
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c    minfx is true if the speciation cannot be solved by speci2
c    error is true if speci2 does not converge
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error, minfx

      integer i,k,id,lord,itic

      double precision g,dp(j3),tdp,gold,xtdp,scp(k5),scptot

      logical pin
      common/ cyt2 /pin(j3)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c---------------------------------------------------------------------
      if (idegen.gt.1000.and.nord(id).gt.1.and.icase(id).ne.0) then
c                                 compositional degeneracy can have the consequence
c                                 that the order parameters are dependent. get the
c                                 composition
         call getscp (scp,scptot,id,1)
c                                 look for degeneracy, this may not work 
c                                 for non-elemental components:
         do i = 1, idegen
            do k = 1, nord(id)
               if (endc(id,lstot(id)+k,idg(i)).ne.0d0) then
                  minfx = .true.
                  return
               end if
            end do
         end do
      end if
c                                 get initial p values
      if (refine) then 
         call nopinc (id,lord)
      else 
         call pinc0 (id,lord)
      end if 
c                                 lord is the number of possible species
      if (lord.lt.nord(id).and.icase(id).ne.0) then
c                                 most likely the model had degenerated
c                                 to ordering across two sites, but because
c                                 the ordered species are made of the same
c                                 endmembers (icase(id)=1), all the ordered
c                                 species are necessary to describe the ordering.
         minfx = .true.

      else if (lord.eq.1) then

         do i = 1, nord(id)
            if (.not.pin(i)) cycle
            call speci1 (g,id,i)


c           call opeci1 (oldg,id,i)

c           if (dabs(oldg-g).gt.1d-6) then 
c              write (*,*) 'oink',id,i,oldg-g,oldg,g
c           end if 

            exit
         end do

      else if (lord.gt.1) then

        if (icase(id).eq.2) then
c                                 compositional degeneracy can have the consequence
c                                 that the order parameters are dependent. get the
c                                 composition
            call getscp (scp,scptot,id,1)
c                                 look for degeneracy, this may not work 
c                                 for non-elemental components:
            do i = 1, idegen
               do k = 1, nord(id)
                  if (endc(id,lstot(id)+k,idg(i)).ne.0d0) then
                     minfx = .true.
                     return
                  end if
               end do
            end do
         end if
c                                 check if an odered species contains the
c                                 degenerate component:
         itic = 0
         gold = 0d0
         xtdp = 0d0
         minfx = .false.

         do

            call gderiv (id,g,dp,.false.,error)

            if (error) then
               badc(1) = badc(1) + 1d0
               exit
            end if

            tdp = 0d0

            do k = 1, nord(id)

               if (.not.pin(k)) cycle

               call pinc (dp(k),k,id,minfx)

               if (dp(k).eq.0d0) then
c                                 search has hit a constraint this should be 
c                                 bad news for newton-raphson, decide what to
c                                 to on the basis of iopt(37) - minfxc_solver:
                  if (iopt(37).eq.0) then
c                                 don't flag as a bad result and continue 
c                                 search
                     minfx = .false.

                  else if (iopt(37).eq.1) then
c                                  just continue, minfx set T by pinc.
                  else if (iopt(37).eq.2) then

                     if (icase(id).eq.0) then 
                        pin(k) = .false.
                     else 
                        return
                     end if

                  else if (iopt(37).eq.3) then

                     pin(k) = .false.

                  else if (iopt(37).ge.4) then

                     if (icase(id).eq.0) pin(k) = .false.
                     minfx = .false.

                  end if

               end if

               tdp = tdp + dabs(dp(k))

            end do

            if ((tdp.lt.nopt(50).or.
     *          dabs((gold-g)/(1d0+dabs(g))).lt.nopt(50))
     *         .and.itic.gt.1) then

c              if (tdp.lt.nopt(52).and.dabs((gold-g)/g).gt.nopt(53))
c    *            then 
c                 write (*,*) 'oink2',gold-g,g,itic,id
c              end if

               goodc(1) = goodc(1) + 1d0
               goodc(2) = goodc(2) + dfloat(itic)
               exit

            else if (itic.gt.5.and.gold.lt.g) then

               minfx = .true.
               exit

            else if (itic.gt.iopt(21)) then 

c              write (*,*) 'div2 ',gold-g,id,itic,g,tdp,tdp-xtdp
               minfx = .true. 
               exit

            else if (itic.gt.5.and.tdp.eq.xtdp) then 

               minfx = .true. 
c              write (*,*) 'wroink67 ',dp(1:lord),id,g
               exit

            end if

            itic = itic + 1

            xtdp = tdp

            gold = g

         end do

      end if

      end

      subroutine pincs (dp,dy,ind,jd,nr)
c----------------------------------------------------------------------
c subroutine to increment the proportions of endmembers involved in a
c predefined ordering reaction (called only by speci1, see pinc and speci2
c for the general case).

c this routine replicates dpinc, but the p's are computed from p0 and
c uses simple arrays.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jd, nr, ind(*)

      double precision dp, dy(*)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      pa(jd) = p0a(jd) + dp

      do i = 1, nr
         pa(ind(i)) = p0a(ind(i)) + dy(i)*dp
      end do

      end

      subroutine pinc (dp,k,id,minfx)
c----------------------------------------------------------------------
c subroutine to increment the k'th species of solution id, if the increment
c violates a stoichiometric limit, it's set to half it's maximum value.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id, jd

      logical minfx

      double precision dp,pmx,pmn

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
c                                 given dp check if it violates
c                                 stoichiometric constraints
      jd = lstot(id) + k

      call plimit (pmn,pmx,k,id)

      if (pa(jd)+dp.gt.pmx) then
         dp = pmx - pa(jd)
      else if (pa(jd)+dp.lt.pmn) then
         dp = pmn - pa(jd)
      end if
c                                 if dp is zero a constraint as been hit
c                                 and newton-raphson is ill-advised.
      if (pa(jd).eq.pmx.or.pa(jd).eq.pmn) then 
         minfx = .true.
      end if
c                                 adjust the composition by the increment
      call dpinc (dp,k,id,jd)

      end

      subroutine dpinc (dp,k,id,jd)
c----------------------------------------------------------------------
c subroutine to increment the k'th species of solution id.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,k,id,jd

      double precision dp

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c----------------------------------------------------------------------
c                                 adjust the composition by the increment
      do i = 1, nrct(k,id)

         pa(ideps(i,k,id)) = pa(ideps(i,k,id))
     *                     + dydy(ideps(i,k,id),k,id)*dp

      end do

      pa(jd) = pa(jd) + dp

      end

      subroutine pinc0 (id,lord)
c----------------------------------------------------------------------
c subroutine set initial species concentrations to half their
c stoichiometric limit. this requires that pa = p0a on entry!
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical degpin

      external degpin

      integer i,k,id,jd,lord,iout

      double precision dp,pmn,pmx,dpp(j3),dinc,tinc
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------

      lord = 0

      if (icase(id).eq.1) then
c                                 case 1: fully correlated
         dinc = 0.5d0/dfloat(nord(id))
         tinc = dinc

         do k = 1, nord(id)

            call plimit (pmn,pmx,k,id)

            if (pmn.ge.pmx.or.pmx-pmn.lt.nopt(50).or.degpin(k,id)) then
               pin(k) = .false.
               cycle
            else
               pin(k) = .true.
               lord = lord + 1
            end if

            jd = lstot(id) + k

            dp = pmn + (pmx - pmn) * tinc - pa(jd)
c                                 adjust the composition by the first increment
            call dpinc (dp,k,id,jd)

            tinc = tinc + dinc

         end do

      else if (icase(id).eq.2.or.icase(id).eq.0) then

         if (icase(id).eq.2) then
c                                 case 2: positive partial correlation
c                                         iterate 4 times for increments
            iout = 5
         else
c                                 case 0: no correlation/iteration
            iout = 1
         end if

         do i = 1, iout

            do k = 1, nord(id)

               call plimit (pmn,pmx,k,id)

               if (i.eq.1) then

                  if (pmn.ge.pmx.or.pmx-pmn.lt.nopt(50)
     *                                         .or.degpin(k,id)) then
                     pin(k) = .false.
                     cycle
                  else
                     pin(k) = .true.
                     lord = lord + 1
                 end if

               end if
c                                 adjust the composition by the first increment
               jd = lstot(id) + k
               dp = pmx - pa(jd)
               pa(jd) = pa(jd) + dp
               dpp(k) = pa(jd) - p0a(jd)

            end do
c                                 no species possible
            if (lord.eq.0) return

         end do
c                                 back off from maximum for final assignements
         do k = 1, nord(id)

            if (.not.pin(k)) cycle

            jd = lstot(id) + k
            pa(jd) = p0a(jd)

            dp = dpp(k)*0.9d0
c                                 adjust the composition by the first increment
            call dpinc (dp,k,id,jd)

         end do

      else if (nord(id).eq.1) then
c                                 only one order parameter, as currently programmed
c                                 this will never be called.
         call plimit (pmn,pmx,1,id)

         if (pmn.ge.pmx) then
            pin(1) = .false.
         else

            pin(1) = .true.
            lord = 1
            jd = lstot(id) + 1

            dp = pmn + (pmx - pmn) * 0.9d0 - pa(jd)
c                                 adjust the composition by the first increment
            call dpinc (dp,k,id,jd)

         end if

      else
c                                 unanticipated case?
         call error (999,p0a(1),i,
     *              'unanticipated correlation between ordered species')

      end if
c                                 check for degenerate compositions
!      if (lord.gt.0) then
!
!         iout = 0
!
!         do i = 1, lstot(id)
!            if (p0a(i).eq.0d0) then
!               iout = iout + 1
!               ibad(iout) = i
!            end if
!         end do
!c                                 the indices of the present components are igood(1..in)
!         if (iout.gt.0) then
!            do k = 1, nord(id)
!               if (pin(k)) then
!c                                 check that the ordered species are in the subcomposition
!                 do j = 1, nrct(k,id)
!                     do i = 1, iout
!                        if (ideps(j,k,id).eq.ibad(i)) then
!                           write (*,*) 'dbug'
!c                          lord = 0
!                           return
!                        end if
!                     end do
!                  end do
!               end if
!            end do
!         end if
!      end if

      end

      logical function degpin (k,id)
c----------------------------------------------------------------------
c check if ordered species k contains a component that the system is 
c degneratue in
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, k, id
c----------------------------------------------------------------------
      degpin = .false.

      do i = 1, idegen
         if (endc(id,lstot(id)+k,idg(i)).ne.0d0) then
            degpin = .true.
            return
         end if
      end do

      end

      subroutine nopinc (id,lord)
c----------------------------------------------------------------------
c subroutine to set lord during refinement
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical degpin

      external degpin

      integer k,id,lord

      double precision pmn,pmx

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------

      lord = 0

         do k = 1, nord(id)

            call plimit (pmn,pmx,k,id)

            if (pmn.ge.pmx.or.pmx-pmn.lt.nopt(50).or.degpin(k,id)) then
               pin(k) = .false.
               cycle
            else
               pin(k) = .true.
               lord = lord + 1
            end if

         end do

      end

      subroutine pcheck (x,xmin,xmax,dx,quit)
c-----------------------------------------------------------------------
c subroutine to increment x for a 1-d root search between xmin and xmax.
c modified to redefine limits according to gradient, Dec 20, 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical quit

      double precision x, xmin, xmax, dx, xt
c-----------------------------------------------------------------------
      quit = .false.

      xt = x + dx

      if (xt.eq.xmin.or.xt.eq.xmax) then
c                                 hit the limit, don't set x to
c                                 the limit to save revaluating x
c                                 dependent variables.
        quit = .true.

        return

      end if

      if (dx.lt.0d0) then
c                                 narrow limit to avoid oscillating
         if (x.lt.xmax) xmax = x
c                                 increment would make x < xmin
c                                 revise the increment
         if (xt.lt.xmin) dx = (xmin - x)/2d0

      else if (dx.gt.0d0) then
c                                 narrow limit to avoid oscillating
         if (x.gt.xmin) xmin = x
c                                 increment would make x > xmax
c                                 revise the increment
         if (xt.gt.xmax) dx = (xmax - x)/2d0

      end if

      x = x + dx
c                                 check if dx has dropped below
c                                 function precision
      if (dabs(dx/(1d0+dabs(x))).lt.nopt(50)) quit = .true.

      end

      subroutine gderi1 (k,id,dg,g)
c----------------------------------------------------------------------
c subroutine computes the newton raphson increment dg from the 1st and 2nd
c derivatives of the g of solution (id) with respect to the concentration
c of the kth ordered species.

c the formulation assumes:

c  1) the speciation reaction is equimolar (see gpder1 for non-equimolar
c     case.

c  2) atomic site fractions are linear functions of the ordered species
c     concentrations (p's).

c  3) the excess function is second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,k,i1,i2,id

      double precision g,dg,d2g,t,s,ds,d2s
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision enth
      common/ cxt35 /enth(j3)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 initialize, d2gx has been set in setw
      g = 0d0

      dg = g
      d2g = d2gx(k,k)

      if (lexces(id)) then

         do i = 1, jterm(id)
c                                 assuming regular terms
           i1 = jsub(1,i,id)
           i2 = jsub(2,i,id)

           g = g + w(i) * pa(i1) * pa(i2)
           dg = dg + w(i) * (pa(i1)*dydy(i2,k,id)
     *                     + pa(i2)*dydy(i1,k,id))

         end do
c                                 get derivative of excess function
         if (llaar(id)) then
c                                 for h&p van laar, this is unnecessary because
c                                 t is constant.
            t = 0d0
c                                 h&p van laar
            do i = 1, nstot(id)
               t = t + alpha(i)* pa(i)
            end do
c                                 coming out of this loop g, dg, and
c                                 d2g  are not the complete functions
c                                 because of the "tphi" term in the
c                                 van laar.

c                                 convert dg and d2g to the full derivative
            g = g/t
            dg = (dg - g*dt(k))/t
            d2g = (d2g - 2d0*dt(k)*dg)/t

         end if

      end if
c                                 get the configurational entropy derivatives
      call sderi1 (k,id,s,ds,d2s)
c                                 in case speci1 is being used for a degnerated
c                                 o/d problem add the enthalpy looping over all
c                                 ordered endmembers
      do i1 = 1, nord(id)
         g = g + pa(lstot(id)+i1)*enth(i1)
      end do

      g = g - v(2)*s
      dg  = dg + enth(k)  - v(2)*ds
      d2g = d2g - v(2)*d2s
c                                 dg becomes the newton raphson increment
      dg = -dg/d2g

      end

      subroutine sderi1 (l,id,s,ds,d2s)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a
c solution with respect to the proportion of the lth ordered species.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,id

      double precision zt,dzdy,dzy,dzyy,zl,ds,d2s,lnz,s,sy
c                                 working arrays
      double precision zz, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)
c----------------------------------------------------------------------
      s = 0d0
      ds = 0d0
      d2s = 0d0

      do i = 1, msite(id)

         sy = 0d0
         dzy = 0d0
         dzyy = 0d0

         zt = 0d0

         do j = 1, zsp(id,i)

            zl = dcoef(0,j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               zl = zl + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do
c                                 sdzdp is (dz(i,j)/dp(l))
            dzdy = sdzdp(l,j,i,id)

            if (zl.lt.nopt(50)) zl = nopt(50)

            zt = zt + zl

            lnz = dlog(zl)
c                                 the entropy
            sy = sy + zl * lnz
c                                 the first derivative
            dzy = dzy - dzdy * (1d0 + lnz)
c                                 the second
            dzyy = dzyy  - dzdy**2 / zl

         end do
c                                 add the contibution from the zsp(id,i)+1th
c                                 species:
         zl = 1d0 - zt

         dzdy = sdzdp(l,j,i,id)

         if (zl.lt.nopt(50)) zl = nopt(50)

         lnz = dlog(zl)
c                                 the entropy
         sy = sy + zl * lnz
c                                 the first derivative is
         dzy = dzy - dzdy * (1d0 + lnz)
c                                 and the second is
         dzyy = dzyy  - dzdy**2 / zl

         s = s - zmult(id,i)*sy
         ds = ds + zmult(id,i)*dzy
         d2s = d2s + zmult(id,i)*dzyy

      end do
c                                 for models with disordered
c                                 endmembers, correct first derivative for the
c                                 change in endmember configurational entropy
      do i = 1, nstot(id)
         s = s - pa(i)*scoef(i,id)
         ds = ds - dydy(i,l,id)*scoef(i,id)
      end do

      end

      subroutine p0limt (id)
c----------------------------------------------------------------------
c subroutine to compute the sums of the p0 terms in ordered species
c limit expressions for solution id.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)
c----------------------------------------------------------------------
      do k = 1, nord(id)
c                                 for ordered species k
         do i = 1, ln(k,id)
c                                 for limit i
            tsum(i,k) = l0c(1,i,k,id)

            do j = 1, lt(i,k,id)
c                                 for term j
               tsum(i,k) = tsum(i,k) + lc(j,i,k,id)*p0a(lid(j,i,k,id))

            end do

         end do

      end do

      end

      subroutine plimit (pmn,pmx,k,id)
c----------------------------------------------------------------------
c subroutine to compute minimum and maximum concentration of ordered
c species k in solution id from site fraction constraints, assumes the
c p0 terms have been accumulated in tsum (routine p0limt)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id

      double precision pmn,pmx,mini
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)
c----------------------------------------------------------------------
      pmx = 1d99
      pmn = -1d99

      do i = 1, ln(k,id)

         mini =  tsum(i,k)

         do j = 1, jt(i,k,id)

            mini = mini + jc(j,i,k,id)*pa(jid(j,i,k,id))

         end do

         if (mini.gt.pmn) pmn = mini
         if (l0c(2,i,k,id)+mini.lt.pmx) pmx = mini + l0c(2,i,k,id)

      end do

      end

      subroutine readlm (idim,tname)
c---------------------------------------------------------------------
c readlm - reads stoichiometric limits on ordered species concentrations
c as of 6.8.7 this is a dummy to read old (6.8.6-) solution model files.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer inds(k7),ict,idim,ier

      double precision coeffs(k7)

      character begin*5, tag*3, tname*10

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      call readcd (n9,ier,.true.)

      write (begin,'(5a)') chars(1:5)

      if (begin.eq.'begin') then

         do
c                                 read the limit equations for the
c                                 amounts of the ordered endmembers
            call readz (coeffs,inds,ict,idim,tname,tag)

            if (tag.eq.'end') exit

         end do

      else

         backspace (n9)

      end if

      end

      subroutine input9 (first)
c-----------------------------------------------------------------------
c given a list of solution phase names (fname(h9)) input9 searches a
c data file (on unit n9) for the relevant data and subdivides the
c solutions into pseudo-compounds.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, im, id, idsol, ixct, gcind, irjct, infnd, ifnd

      logical first, chksol, wham, ok, found

      character sname(h9)*10, new*3, tn1*6, tn2*22, rjct(h9)*10, 
     *          nfnd(h9)*10

      double precision zt

      external chksol

      character prject*100,tfname*100
      common/ cst228 /prject,tfname

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      character mname*8
      common/ cst18a /mname(m4)

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer ikp
      common/ cst61 /ikp(k1)

      integer iam
      common/ cst4 /iam

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq
c-----------------------------------------------------------------------
c                                 initialize counters
      ixct = 0
c                                 global compositional index counter
      gcind = 0
c                                 initialize model counter
      im = 0
c                                 rejected model counter
      irjct = 0
c                                 a flag to check if more than one solution model
c                                 references an internal molecular EoS.
      wham = .false.
c                                 pointer to aq solution model
      idaq = 0
c                                 no request for solutions
      if (io9.eq.1.or.isoct.eq.0) then

         isoct = 0
c                                 identify the fluid for aqrxdo
         call aqidst

         return

      end if
c                                 open pseudocompund list file
      if ((iam.eq.1.or.outprt).and.lopt(10).and.
     *    (iam.eq.1.or.iam.eq.15)) then

         call mertxt (tfname,prject,'_pseudocompound_list.txt',0)
         open (n8,file=tfname)

      end if
c                                 format test line
      read (n9,'(a)') new
c                                 check version compatability
      if (.not.chksol(new)) call error (3,zt,im,new)

      if (iam.lt.3.or.iam.eq.15)
     *           write (*,'(80(''-''),/,a,/)') 'Solution model summary:'

      do while (im.lt.isoct)
c                                 -------------------------------------
c                                 read the solution model
         call rmodel (tn1,tn2)
c                                 istot is zero, if eof:
         if (istot.eq.0) exit
c                                 -------------------------------------
c                                 check the solution model:
         call cmodel (im,idsol,first,found)

         if (jstot.eq.1.and.jsmod.eq.39.and.lopt(32)) then
c                                  lagged aqueous speciaton with a pure water solvent.
         else
c                                  normal solution.
            if (jstot.lt.2) then
               if (found) then
                  irjct = irjct + 1
                  rjct(irjct) = tname
               end if 
               cycle
            end if
c                                 -------------------------------------
c                                 reformulate the model so that it has
c                                 no missing endmembers:
            if (jstot.lt.istot) call reform (im,first)

            if (istot.lt.2) then
               irjct = irjct + 1
               rjct(irjct) = tname
               cycle
            end if

         end if
c                                 -------------------------------------
c                                 make various book keeping arrays (y2p,
c                                 jmsol, dydz, .....)
         call nmodel
c                                 check that the name has not already been found, i.e.,
c                                 that the name is duplicated in the solution model file
         if (first) then

            do i = 1, im - 1
               if (tname.eq.sname(i)) call error (75,0d0,i,tname)
            end do

         end if

         if (jsmod.eq.39.and.ifct.gt.0) then
c                                 check that a GFSM model is not being 
c                                 used together with a saturated fluid 
c                                 constraint
               write (*,1060) tname, tname

               call wrnstp

         end if 
c                                 save solution name
         sname(im) = tname
c                                 abbreviation
         aname(im) = tn1
c                                 long name
         lname(im) = tn2
c                                 save found solutions in global solution
c                                 model arrays
         call gmodel (im,wham)
c                                 generate pseudocompound compositions.
c                                 subdiv returns the total
c                                 number of pseudocompounds (ipcps) and
c                                 array y, of which element y(h,i,j) is
c                                 the site fraction of the jth species on
c                                 the ith site of the hth pseudocompound.
         if (iam.lt.3.or.iam.eq.15) then
c                                 vertex/meemum need static pseudocompounds
            do i = 1, kstot

               id = kdsol(knsp(i,im))
               if (iend(knsp(i,im)).eq.0) ikp(id) = im

            end do

            ophct = iphct

            if (outprt.and.lopt(10)) then 
c                                 write solution model name/endmembers for pseudocompound list file:
               if (lrecip(im)) then
                  j = mstot(im)
               else 
                  j = lstot(im)
               end if

               write (n8,1060) tname,(mname(iorig(knsp(i,im))), i= 1, j)

            end if

            if (.not.refine.or.iam.eq.15) then 
c                                 subdiv discretizes the composition of the 
c                                 solution and stores the data (soload)
               call subdiv (im,gcind)

               if (iphct-ophct.gt.0) then
c                                 write pseudocompound count
                  write (*,1100) iphct-ophct, tname
c                                 indicate site_check_override and refine endmembers
                  if (bdx(im)) write (*,1080) tname
                  if (.not.nrf(im).and..not.lopt(39)) 
     *               write (*,1090) tname

               end if

               jend(im,2) = iphct - ophct

            end if

         end if
c                               read next solution
      end do
c                               make lists of found/not-found solutions
      infnd = 0
      ifnd = 0

      do i = 1, isoct

         ok = .false.
c                                 check if fname was included:
         do j = 1, im
            if (fname(i).eq.sname(j)) then 
               ok = .true.
               ifnd = ifnd + 1
               solptr(ifnd) = j
               exit
            end if
         end do

         if (ok) cycle
c                                  check if fname was rejected:
         do j = 1, irjct
            if (fname(i).eq.rjct(j)) then 
               ok = .true.
               exit 
            end if
         end do

         if (ok) cycle
c                                  add to not found list:
         infnd = infnd + 1
         nfnd(infnd) = fname(i)

      end do

      if (iam.lt.3.or.iam.eq.15) then
c                                  total pseudocompound count:
         if (.not.refine.or.iam.eq.15) write (*,1110) iphct - ipoint
c                                  list of found solutions
         if (im.gt.0) then
            write (*,'(/,a,/)') 'Summary of included solution models:'
            write (*,'(8(a,1x))') (sname(i),i= 1, im)
         else
            write (*,'(/,a,/)') 'No solution models included!'
         end if

         if (irjct.gt.0) then
            write (*,'(/,a,/)') 'Summary of rejected solution models '//
     *                             '(see warnings above for reasons):'
         
            write (*,'(8(a,1x))') (rjct(i),i= 1, irjct)
         end if

         if (infnd.gt.0) then
            write (*,'(/,a,/)') 
     *             'Requested solution models that were not found:'
            write (*,'(8(a,1x))') (nfnd(i),i= 1, infnd)
         end if
c                               flush for paralyzer's piped output
         flush (6)
c                               scan for "killed endmembers"
         do i = 1, ipoint
c                               reset ikp
            if (ikp(i).lt.0) ikp(i) = 0
         end do

         if (io3.eq.0.and.outprt.and.(iam.eq.1.or.iam.eq.15)
     *                                     .and.isoct.ne.im) then

            write (n3,1020)
            write (n3,1010) (fname(i), i = 1, isoct)
            if (im.gt.0) then
               write (n3,1000)
               write (n3,1010) (sname(i), i = 1, im)
            else if (outprt) then
               write (n3,1040)
            end if

         end if

      end if

      do i = 1, im
         fname(i) = sname(i)
      end do

      isoct = im
c                              identify the fluid for aqrxdo
      call aqidst
c                              close pseudocompound list
      if ((iam.eq.1.or.outprt).and.lopt(10)) close (n8)
c                              close solution model file
      close (n9)

      write (*,'(80(''-''))')

      first = .false.

1000  format (/,'the following solution models will be considered:',/)
1010  format (7(2x,a10))
1020  format (/,'Of the requested solution models:',/)
1030  format ('**warning ver535** ',a,' is a generic fluid solution ',
     *        'model (GFSM) the presence',/,'of which is inconsistent ',
     *        'with saturated phase constraints if the saturated phase',
     *      /,'is a fluid. Possible courses of action are:',//,4x,
     *        '1) remove ',a,' and restart.',/,4x,
     *        '2) remove the phase saturation constraint and restart.',/
     *    ,4x,'3) ignore this warning and continue execution.',//,
     *        'Continue (Y/N)?')
1040  format (/,'no models will be considered.',/)
1060  format (/,'Solution: ',a,/,12x,'Endmember fractions:',
     *        /,12x,20(a,1x))
1080  format (9x,'reject_bad_compositions is on for ',a)
1090  format (9x,'refine_endmembers is on for ',a)
1100  format (i8,' pseudocompounds generated for: ',a)
1110  format (/,'Total number of pseudocompounds:',i8)
1120  format (/,'1 - Although the bulk composition of pseudocompounds'
     *        ,' for this solution is fixed,',/,' the proportions of'
     *        ,' its endmembers may vary due to respeciation.',/)
1130  format (/,'2 - Proportions output here may sum to <1 ',
     *          'because the ordered species',/,'may have non-zero ',
     *          'initial proportions.',/)
      end

      subroutine err41 (tag)
c---------------------------------------------------------------------
c if an entry will exceed dimension 'tag' and write apporpriate
c diagnostic on error.
c---------------------------------------------------------------------
      implicit none

      character tname*10, tag*(*)

      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c---------------------------------------------------------------------
c                                 error diagnostic
         if (refine) then
            call error (41,0d0,1,tag)
         else
            call error (41,0d0,0,tag)
         end if

      end

      subroutine satsrt
c---------------------------------------------------------------------
c routine to sort pseudocompounds consisting entirely of saturated
c components.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j,idc

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      do j = isat, 1, -1
         idc = icp + j
         if (cp(idc,iphct).ne.0d0) then
            isct(j) = isct(j) + 1
            if (isct(j).gt.h6) call error (17,cp(1,1),h6,'SATSRT')
            if (iphct.gt.k1) call error (72,cp(1,1),k1,
     *                                  'SATSRT increase parameter k1')
            ids(j,isct(j)) = iphct
            exit
         end if
      end do

      end

      subroutine soload (im,bad)
c--------------------------------------------------------------------------
c soload - loads/requires solution properties:

c   jend(h9,m14+2)  - h9 is the maximum number of solutions
c                   k12 is the maximum number of endmembers pers
c                   solution plus two.
c   jend(i,1)     - OBSOLETE! is the number of endmembers in solution i.
c   jend(i,2)     - is the number of pseudocompounds of solution i.
c   jend(i,3-3+j) - are the indices of the j endmembers in solution i.
c   sxs(k13)      - contains the mole fractions of the endmembers
c                   in the pseudocompounds.
c   ikp(i)        - the index of the solution corresponding to pseudocompound i.
c   exces(j,i)    - the excess function of pseudocompound i, accounts for
c                   excess properties and configurational entropy as a function
c                   of pressure and temperature:

c                       gexces(i) = exces(1) + exces(2)*T + exces(3)*P
c--------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character znm(3,2)*2, pnm(3)*2

      double precision zpr,omega,scp(k5)

      logical bad

      integer im, h, i, j, l, index, i228, oim

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      double precision wgl, wkl, vlar
      common/ cxt2r /wgl(m3,m1,h9),wkl(m16,m17,m18,h9),vlar(m3,m4,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ikp
      common/ cst61 /ikp(k1)

      double precision pa, p0a, zp, w, y, z, wl, pp
      common/ cxt7 /y(m4),zp(m4),pa(m4),p0a(m4),z(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision exces
      common/ cst304 /exces(m3,k1)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer iam
      common/ cst4 /iam

      save i228,oim
      data i228,oim/0,0/
c----------------------------------------------------------------------
c                                 reject special case:
c                                 ternary coh fluids above the CH4-CO join
      if (ksmod(im).eq.41.and.pa(1).ge.r13+pa(2)) then 
         bad = .true.
         return
      end if

      bad = .false.

      ikp(iphct) = im
c                                 -------------------------------------
c                                encode a name, this is archaic and only relevant 
c                                for CONVEX which is unlikely to be effective for
c                                multi-polytope composition models. 
      if (istg(im,1).gt.1) then
c                                make character nums for standard cases
c                                this is only to avoid run-time errors
c                                during debugging.
         do i = 1, istg(im,1)

            do j = 1, 2

               h = idint(1d2*z(1,i,j))

               if (h.eq.100.or.h.lt.0) then
                  znm(i,j) = '**'
               else
                  write (znm(i,j),'(i2)') h
               end if

            end do

         end do

      end if

      do j = 1, 3

         h = idint(1d2*pa(j))

         if (h.ge.1d2.or.h.lt.0d0) then
            pnm(j) = '**'
         else
            write (pnm(j),'(i2)') h
         end if

c use mname array to flag retained absent endmembers

      end do

      if (istg(im,1).eq.2.and.mstot(im).eq.4) then
c                                special case 1, bin-bin reciprocal solution
         write (names(iphct),1020) tname, znm(1,1),znm(2,1)

      else if (istg(im,1).eq.2.and.mstot(im).eq.6.and.ispg(im,1,1).eq.3)
     *        then
c                                special case 2, tern-bin reciprocal solution
         write (names(iphct),1060) tname, znm(1,1),znm(1,2),znm(2,1)

      else if (istg(im,1).eq.2.and.mstot(im).eq.6.and.ispg(im,1,1).eq.2)
     *        then
c                                special case 3, bin-tern reciprocal solution
         write (names(iphct),1060) tname, znm(1,1),znm(2,1),znm(2,2)

      else if (istg(im,1).eq.2.and.mstot(im).eq.9) then
c                                special case 4, tern-tern reciprocal solution
         write (names(iphct),1060) znm(1,1),znm(1,2),znm(2,1),znm(2,2)

      else if (mstot(im).eq.2) then
c                                binary solutions
         if (pa(1).gt.0.9999d0) then
            write (names(iphct),'(a3,a)') names(jend(im,3)),'_100*'
         else if (pa(1).ge.0.98d0) then
            write (names(iphct),'(a2,a,f5.2)')
     *             names(jend(im,3)),'_',1d2*pa(1)
         else if (pa(1).lt.1d-6) then
            write (names(iphct),'(a3,a)') names(jend(im,3)),'_0*'
         else if (pa(1).lt.0.02d0) then
            write (names(iphct),'(a2,a,f5.3)')
     *             names(jend(im,3)),'_',1d2*pa(1)
         else
            write (names(iphct),1070) names(jend(im,3)),1d2*pa(1)
         end if

      else if (mstot(im).eq.3) then
c                                ternary solutions
         write (names(iphct),1060) (names(jend(im,2+j)),
     *                             pnm(j), j = 1, 2)
      else if (mstot(im).eq.4) then
c                                quaternary solutions
         write (names(iphct),1060) tname, (pnm(j), j = 1, 3)

      else
c                                all the rest:
         if (iphct.lt.1000000) then
            write (names(iphct),1080) tname, iphct
         else if (iphct.lt.10000000) then
            write (names(iphct),1100) tname, iphct
         else
            write (names(iphct),'(i8)') iphct
         end if

      end if
c                                 get blanks out of name:
      if (mstot(im).lt.4) then
         call unblnk (names(iphct))
      else
         call reblnk (names(iphct))
      end if
c                                 -------------------------------------
      if ((iam.eq.1.or.outprt).and.lopt(10)) then

         if (lrecip(im)) then
            h = mstot(im)
         else 
            h = lstot(im)
         end if
c                                 write composition name to pseudocompound list file
          write (n8,1050) names(iphct),(pa(j), j = 1, h)

      end if
c                                 bulk composition stuff
      rkwak = .true.

      call getscp (scp,ctot(iphct),im,1)

      do l = 1, icomp

         if (scp(l).gt.-nopt(50).and.scp(l).lt.nopt(50)) then
            scp(l) = 0d0
         else if (scp(l).lt.0d0.and.im.ne.i228) then
            i228 = im
            call warn (228,scp(l),l,tname)
         end if

      end do
c                                 check if the phase consists
c                                 entirely of saturated components:
      if (ctot(iphct).lt.nopt(50)) then

         if (im.ne.oim) call warn (55,scp(1),l,tname)

         bad = .true.
         oim = im

         return

      end if
c                                 load the static composition matrix
      if (iam.eq.1.or.iam.eq.2) then 
c                                 MEEMUM/VERTEX
         do j = 1, icp
            a(j,iphct-jiinc) = scp(j)/ctot(iphct)
         end do

      else if (iam.eq.15) then 
c                                 CONVEX
         do j = 1, icomp
            a(j,iphct) = scp(j)
         end do

      end if
c                                 -------------------------------------
c                                 this section loads excess, configurational and
c                                 dqf contributions for static compositions
c                                 into exces; the messiness is probably 
c                                 not worth the effort

      exces(1:m3,iphct) = 0d0

      if (.not.lorder(im)) then 
c                                 configurational negentropy:
         if (msite(im).ne.0) exces(2,iphct) = -omega(im,pa)
c                                 load excess terms, if not Laar or ordered:
         if (extyp(im).eq.0) then

            do i = 1, jterm(im)

               zpr = 1d0

               do j = 1, rko(i,im)
                  zpr = zpr * pa(jsub(j,i,im))
               end do

               do j = 1, m3
                  exces(j,iphct) = exces(j,iphct) + zpr * wgl(j,i,im)
               end do

            end do

         else if (extyp(im).eq.1) then
c                                 redlich kister; expand polynomial
c                                 G Helffrich, 4/16
            do i = 1, jterm(im)
               do j = 1, rko(i,im)

                  zpr = pa(jsub(1,i,im))*pa(jsub(2,i,im))
     *                * (pa(jsub(1,i,im)) - pa(jsub(2,i,im)))**(j-1)
c                                 Sloppy fix for the linear P term
c                                 being wkl(6), i.e., 'wP', wkl(3..5)
c                                 wP0, wP1, wP2, are for brosh's murnaghan
c                                 and are always identified as a special 
c                                 case. This needs to be fixed, i.e., RK
c                                 should be identified as a special case 
c                                 like LAAR
                  do l = 1, 3

                     if (l.lt.3) then 
                        h = l
                     else 
                        h = 6
                     end if

                     exces(l,iphct) = exces(l,iphct)
     *                                + zpr * wkl(l,j,i,im)

                  end do

               end do
            end do
         end if
      end if
c                              dqf corrections are also be saved in the
c                              exces array this implies that speciation
c                              does not effect the amount of the dqf'd
c                              endmembers.
      do i = 1, jdqf(im)
c                              index points to the endmember in the full
c                              model:
         index = jndq(i,im)

         do j = 1, m3
            exces(j,iphct) = exces(j,iphct) + pp(index)*dqfg(j,i,im)
         end do

      end do

1020  format (a2,a2,'_',a2)
1050  format (a,2x,20(1x,f6.3,2x))
1060  format (a2,a2,a2,a2)
1070  format (a3,'_',f4.1)
1080  format (a2,i6)
1100  format (a1,i7)

      end

      subroutine reload (file)
c----------------------------------------------------------------------
c load the saved exploratory stage compositions into the static array
c for the auto-refine stage.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad, file

      integer id, i, j, ntot, ids, tmp, tco, ttot

      character sname(h9)*10, tag*11

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer ikp
      common/ cst61 /ikp(k1)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct,istart
      common/ cst111 /jphct,istart

      integer is
      double precision a,b,c
      common/ cst313 /a(k5*k1),b(k5),c(k1),is(k1+k5)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      if (file) then 
c                                 manual auto-refine, read the static 
c                                 compositions from the arf file
         read (n10,*) i
         read (n10,'(7(a,1x))') sname(1:i)
         read (n10,*) jend(1:i,2)
         if (i.ne.isoct) call error (63,y(1),i,'RELOAD/isoct')

         tcct = 0

         do i = 1, isoct

            if (sname(i).ne.fname(i)) 
     *         call error (63,y(1),i,'RELOAD/sname')

            tmp = jend(i,2)*tstot(i)

            read (n10,*) txco(tcct+1:tcct+tmp)

            tcct = tcct + tmp

         end do

         tcct = 0
         tpct = 0

         do i = 1, isoct

            ttot = tstot(i)

            do j = 1, jend(i,2)

               tpct = tpct + 1
               itxp(tpct) = tcct
               tcct = tcct + ttot

            end do

         end do

         iphct = ipoint + tpct

      else
c                                 automatic, read the data from memory
         if (.not.refine) then
c                                 in the first loop:
c                                 load stable static compositions to 
c                                 the list of dynamic compositions
            do i = ipoint + 1, iphct

               if (ststbl(i).or.lopt(30)) then
                  ids = ikp(i)
                  call setxyp (ids,i,bad)
                  if (bad) cycle
                  call savdyn (ids)
               end if

            end do

         else if (lopt(55)) then
c                                 re-refine option:
c                                 load the whole damn thing
            zcoct = 0
            tmp = 0 
c                                 compare the initial 
c                                 ipoint + 1: iphct 
c                                 compositions to the 
c                                 iphct - ipoint + 1 : tpct
c                                 new compositions
            stpct = iphct - ipoint + 1
c                                 at this point the static compositions
c                                 are in the first iphct-ipoint positions
c                                 in txco and the newly added compositions
c                                 are at iphct-ipoint:tpct

c                                 sort through the static compositions
c                                 to see if any (or all should be added)
c                                 to the future static array:
            do i = 1, isoct

               ntot = nstot(i)
               ttot = tstot(i)

               do j = 1, jend(i,2)

                  tmp = tmp + 1

                  if (ststbl(ipoint+tmp).or.lopt(30)) then

                     pa(1:ntot) = txco(itxp(tmp)+1:itxp(tmp)+ntot)
c                                 only for pp comparison
                     if (lorder(i)) call makepp (i)

                     call savdyn (i)

                     zcoct = zcoct + ttot
                     if (tcct+ttot.gt.m25) call errdbg ('increase m25')

                  end if
               end do
            end do
c                                 now shift the compositions for the future
c                                 static array to occupy the first 1:tpct
c                                 positions of txco and reset jend
            jend(1:isoct,2) = 0
            tmp = 0
            tco = 0
            zcoct = 0
c                                 iphct + 1 cause any old compounds
c                                 to be user have been added...
            do i = iphct + 1, tpct
c                                 solution model pointer
               ids = dkp(i)
c                                 curent position
               tco = itxp(i)
c                                 shift composition
               ttot = tstot(ids)

               txco(zcoct+1:zcoct+ttot) = txco(tco+1:tco+ttot)
c                                 counters, pointers:
               tmp = tmp + 1
               dkp(tmp) = ids
               itxp(tmp) = zcoct
               zcoct = zcoct + ttot
c                                 static counter
               jend(ids,2) = jend(ids,2) + 1

            end do

            tpct = tmp

         end if

         id = 0
         zcoct = 0

         do i = 1, isoct

            jend(i,2) = 0

            ttot = tstot(i)
c                                 for each solution cycle 
c                                 through the unsorted compositions
            do j = 1, tpct

               if (i.ne.dkp(j)) cycle
c                                 its a composition of solution i
               id = id + 1
               jend(i,2) = jend(i,2) + 1
c                                 load temporarily into the static compound 
c                                 a array
               is(id) = zcoct
               a(zcoct+1:zcoct+ttot) = txco(itxp(j)+1:itxp(j)+ttot)

               zcoct = zcoct + ttot

            end do

         end do
c                                 sort results 
         zcoct = 0
         id = 0
c                                 copy the sorted results back into txco
         do i = 1, isoct

            ttot = tstot(i)

            do j = 1, jend(i,2)

               id = id + 1

               txco(zcoct+1:zcoct+ttot) =  a(is(id)+1:is(id)+ttot)
               itxp(id) = zcoct
               zcoct = zcoct + ttot

            end do

         end do

      end if
c                                 reset iphct and reload static
      iphct = ipoint
      id = 0

      if (refine.and.lopt(55)) then
         tag = 'Re-refine  '
      else
         tag = 'Exploratory'
      end if

      write (*,'(80(''-''),/,a,'' stage generated:'',/)') tag

      do i = 1, isoct

         ntot = nstot(i)
c                                 set tname for soload diagnostics
         tname = fname(i)

         write (*,1100) jend(i,2), tname
         id = id + jend(i,2)

         do j = 1, jend(i,2)

            iphct = iphct + 1
            tmp = itxp(iphct - ipoint)
            dkp(iphct - ipoint) = i

            pa(1:ntot) = txco(tmp + 1:tmp + ntot)

            call makepp (i)

            call soload (i,bad)

         end do

      end do
c                                 reset counters, cold start, etc
      call initlp

      stpct = tpct + 1

      write (*,1110) tpct
      write (*,'(80(''-''))')

1100  format (i8,' compositions for: ',a)
1110  format (/,'Total number of compositions:',i8)

      end

      subroutine chkpa (ids)
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ids

      double precision sum

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c---------------------------------------------------------------------
      sum = 0d0

      do i = 1, nstot(ids)
         sum = sum + pa(i)
      end do

      if (dabs(sum-1d0).gt.nopt(50)) then 
         write (*,*) 'wowonka ',sum
      end if

      end

      subroutine initlp 
c--------------------------------------------------------------------
c initialize arrays and constants for lp minimization of static
c compositions.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision bl,bu
      common/ cstbup /bl(k1+k5),bu(k1+k5)

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct,istart
      common/ cst111 /jphct,istart

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer ipoint,kphct,imyn
      common/ cst60  /ipoint,kphct,imyn

      integer tphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct
c-----------------------------------------------------------------------
c                                 locate last point in dynamic/static lp arrays
      jpoint = ipoint - jiinc
      jphct = iphct - jiinc

      ctotal = 0d0

      do i = 1, icp
         ctotal = ctotal + cblk(i)
      end do 
c                                 composition constraint, normalized for reasons
c                                 of stupidity
      do i = 1, icp
         b(i) = cblk(i)/ctotal
      end do 
c                                 static/dynamic composition arrays for solutions
c                                 are loaded in soload/reload/loadgx. stoichiometric
c                                 compounds/endmembers loaded here:
      do i = 1, jpoint
         id = i + jiinc
c                                 jkp indicates which phase a dynamic composition 
c                                 is associated with, < 0 a compound, > 0 a solution
         jkp(i) = -id
         hkp(i) = 0
c                                 load all compounds into the static and dynamic 
c                                 arrays
         a(1:icp,i) = cp(1:icp,id)/ctot(id)
         cp2(1:icp,i) = a(1:icp,i)

      end do
c                                 stability flag for static compositions
      ststbl = .false.
c                                 cold start istart = 0
      istart = 0

      bl(1:jphct) = 0d0
      bu(1:jphct) = 1d0

      bl(jphct+1:jphct+icp) = b(1:icp)
      bu(jphct+1:jphct+icp) = b(1:icp)

      end

      double precision function gkomab (id,jd,vdp)
c---------------------------------------------------------------------
c evaluate g for iron according to the EoS of Komabayashi & Fei (JGR,2010)
c id points to a phase of iron
c jd points to its parameters in thermo
c vdp is the vdp integral for all phases except HCP
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,jd

      double precision  g,vdp

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      g = thermo(4,jd) + t*(thermo(5,jd) + thermo(6,jd)*dlog(t)
     *  + t*(thermo(7,jd) + t*thermo(8,jd))) + thermo(9,jd)/t

      if (id.eq.600) then
c                                 BCC iron
         if (t.gt.1811d0) then
            g = -25383.581d0 + t*(299.31255d0 - 46d0*dlog(t))
     *                       + 2.29603d31*t**(-9)
         end if

      else if (id.eq.601) then
c                                 FCC iron
           g = g - 2476.28 * dsqrt(t)

      else if (id.eq.602) then
c                                 HCP iron
           g = g - 2476.28 * dsqrt(t)
c                                 vdp from daewaele EoS
      else if (id.eq.603) then
c                                 liquid iron, destabilize at T < 1811
      end if

      gkomab = g + vdp

      end

      double precision function glacaz (id)
c---------------------------------------------------------------------
c evaluate various CALPHAD g functions
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision hserfe, hsersi, crbcc, hserc, fefcc

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c---------------------------------------------------------------------
      if (id.eq.610) then
c                                 Fe-bcc
         glacaz = hserfe(t)

      else if (id.eq.611) then
c                                 Si-bcc
         glacaz = 0.47D5 - 0.225D2*t + hsersi(t)

      else if (id.eq.612) then
c                                 Fe-fcc
         glacaz = fefcc(t)

      else if (id.eq.613) then
c                                 Si-fcc
         glacaz = 0.51d5 - 0.218D2*t + hsersi(t)

      else if (id.eq.614) then
c                                 Fe-liq
         if (t.lt.1811d0) then
            glacaz = 0.1204017d5 - 0.655843d1 * t - 0.36751551d-20 *
     *               t ** 7 + hserfe(t)
         else
            glacaz = -0.108397d5 + 0.291302d3 * t - 0.46d2 * t *
     *                dlog(t)
         end if

      else if (id.eq.615) then
c                                 Si-liq
         if (t.lt.1687d0) then
            glacaz = 0.506964D5 - 0.300994D2*t + 0.209307d-20*t**7
     *             + hsersi(t)
         else
            glacaz = 0.49828D5 - 0.295591D2*t + 0.420369d31/t**9 +
     *               hsersi(t)

         end if

      else if (id.eq.616) then
c                                 Fe2Si
         glacaz = -0.237522D5 - 0.354D1*t + 0.67D0*hserfe(t)
     *          + 0.33D0 * hsersi(t)

      else if (id.eq.617) then
c                                 Fe5si3
         glacaz = -0.30143D5 + 0.27D0*t + 0.625D0*hserfe(t)
     *          + 0.375D0 * hsersi(t)

      else if (id.eq.618) then
c                                 FeSi
         glacaz = -0.363806D5 + 0.222D1*t + hserfe(t)/2D0
     *           + hsersi(t)/2D0

      else if (id.eq.619) then
c                                 FeSi2
         glacaz = -0.27383D5 + 0.348D1*t + 0.33D0*hserfe(t)
     *           + 0.67D0*hsersi(t)

      else if (id.eq.620) then
c                                 Fe3Si7
         glacaz = -0.19649D5 - 0.92D0*t + 0.3D0 * hserfe(t)
     *           + 0.7D0 * hsersi(t)

      else if (id.eq.621) then
c                                 Si-diamond
         glacaz = hsersi(t)

      else if (id.eq.622) then
c                                 FeC-BCC
         glacaz = hserfe(t) + 0.269943d6 + 0.587857d3 * t
     *            - 0.729d2 * t * dlog(t) - 0.14169d-2 * t ** 2
     *            + 0.76878d7 / t - 0.7929d9 /t ** 2 +
     *            0.360d11 / t ** 3

      else if (id.eq.623) then
c                                 SiC-BCC
         glacaz = 0.47D5 - 0.225D2*t + hsersi(t) + 0.269944677d6 +
     *            0.436523d3 * t - 0.729d2 * t * dlog(t) - 0.14169d-2 *
     *            t ** 2 + 0.76878d7 / t - 0.7929d9 / t ** 2 +
     *            0.36d11 / t ** 3

      else if (id.eq.624) then
c                                 FeC-FCC
         if (t.lt.1811d0) then

            glacaz = 0.58376159d5 + 0.163135d3 * t - 0.2545D2 *
     *               t * dlog(t) + 0.1677d-3 * t ** 2 + 0.256260d7
     *               / t - 0.2643d9 / t ** 2 + 0.12D11 / t ** 3 +
     *               hserfe(t)

         else

            glacaz = 0.32740293d5 + 0.45510556d3 * t - 0.703d2 * t
     *               * dlog(t) - 0.4723d-3 * t ** 2 +
     *               0.25626d7 / t - 0.2643D9 / t ** 2 + 0.12D11 /
     *               t ** 3 + 0.278854D32 / t ** 9

         end if

      else if (id.eq.625) then
c                                 SiC-FCC
         glacaz = hsersi(t) - 0.37879d5 + 0.20943d3 * t - 0.243d2 * t
     *          * dlog(t) - 0.4723d-3 * t ** 2 + 0.25626d7 / t -
     *          0.2643d9 / t ** 2 + 0.120d11 / t ** 3

         else if (id.eq.626) then
c                                 C-LIQ
         glacaz = 117369d0 - 24.63*t + hserc(t)

      else if (id.eq.627) then
c                                 C-GPH
         glacaz = -0.17368441d5 + 0.17037d3 * t - 0.243d2 * t * dlog(t)
     *            - 0.4723d-3 * t ** 2 + 0.25626d7 / t - 0.2643d9 /
     *            t ** 2 + 0.12d11 / t ** 3

      else if (id.eq.628) then
c                                 SiC
         if (t.lt.700d0) then

            glacaz = -0.85572264D5 + 0.1732005D3 * t - 0.25856D2 * t *
     *                dlog(t) - 0.2107D-1 * t ** 2 + 0.32153D-5 * t ** 3
     *                + 0.438415D6 / t

         else if (t.gt.700d0.and.t.lt.2100) then

            glacaz = -0.95145902d5 + 0.300346d3 * t - 0.45093d2 * t *
     *                dlog(t) - 0.367d-2 * t ** 2 + 0.22d-6 * t ** 3 +
     *                0.1341065d7 / t

         else

            glacaz = -0.105007971d6 + 0.360309d3 * t - 0.53073d2 * t *
     *                dlog(t) - 0.74525d-3 * t ** 2 + 0.173167d-7 *
     *                t ** 3 + 0.3693345d7 / t

         end if

      else if (id.eq.629) then
c                                 Cementite
         glacaz = -0.10745d5 + 0.70604d3 * t - 0.1206d3 * t * dlog(t)


      else if (id.eq.630) then
c                                 Fe8Si2C
         glacaz = -0.210043d5 + 0.506d0 * t + 0.91d-1 * (-0.17368441d5
     *            + 0.17037d3 * t - 0.243d2 * t * dlog(t) - 0.4723d-3
     *            * t ** 2 + 0.25626d7 / t - 0.2643d9 / t ** 2 +
     *            0.12d11 / t ** 3) + 0.727d0 * hserfe(t) +
     *            0.182d0 * hsersi(t)

      else if (id.eq.631) then
c                                 C diam
         glacaz = -0.16359441D5 + 0.17561D3 * t - 0.2431D2 * t *
     *             dlog(t) - 0.4723D-3 * t ** 2 + 0.2698D7 / t -
     *             0.261D9 / t ** 2 + 0.111D11 / t ** 3

      else if (id.eq.632) then
c                                 Cr-BCC
           glacaz = crbcc(t)

      else if (id.eq.633) then
c                                 Cr-FCC
           glacaz = crbcc(t) + 7284d0 + 0.163d0*t

      else if (id.eq.634) then
c                                 Cr_LIQ
         if (t.lt.2180d0) then

            glacaz = crbcc(t) + 0.2433593D5 - 0.1142D2 * t +
     *               0.237615D-20 * t ** 7

         else

             glacaz = -0.16459d5 + 0.335618d3 * t -
     *                0.50d2 * t * dlog(t)

         end if

      else if (id.eq.635) then
c                                 FeCr-sig after Hertzman, 1980; used by Nastia
c          glacaz = (183802.0638d0 + 26d0*hserfe(t) +
c    *               4d0*crbcc(t))/30d0
c                                 FeCr-sig after Andersson & Sundman 1987; used by George
         glacaz = (8d0*fefcc(t) + 4d0*crbcc(t) + 18d0*hserfe(t) +
     *            117 300d0 - 95.96d0*t) / 30d0

      else if (id.eq.636) then
c                                 CrFe-sig after Hertzman, 1980; used by Nastia
c        glacaz =(-22624.05686d0 + 20d0*crbcc(t)
c    *            + 10d0*hserfe(t))/30d0
c                                 CrFe-sig after Andersson & Sundman 1987; used by George
         glacaz = (8d0*fefcc(t) + 22d0*crbcc(t)
     *            + 92 300d0 - 95.96d0*t) / 30d0

      else if (id.eq.637) then
c                                 Fe7C3 after Djurkovic et al., 2011
         glacaz =-0.2345062954D5 + 0.1761006488D4 * t -
     *            0.2975999679D3 * t * dlog(t) - 0.3148668241D-3 *
     *            t ** 2 + 0.1708400854D7 / t - 0.1762000881D9 /
     *            t ** 2 + 0.8000004D10 / t ** 3

      end if

      end

      subroutine chopit (ycum,wt,jst,jsp,lsite,lpoly,ids,jump,extra)
c---------------------------------------------------------------------
c subroutine to do cartesian or transform subdivision of species
c jst+1 through jsp on site k of solution ids. ycum is the smallest
c fraction possible (i.e., if the minimum bound for some species
c is > 0). the fractions are loaded into simp.

c extra - save space for an extra coordinate in simp (ksmod 20 or 9)
c wt   - factor to modify default resolution (1, except ksmod 9)
c jump  - offset for storing coordinates in simp (ksmod 9)
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer mres

      parameter (mres=12000)

      logical extra

      integer mode, ind(ms1), iy(ms1), jsp, lsite, indx, iexit,
     *        ieyit, i, j, k, ids, ic, jst, jump, lpoly

      double precision y(ms1,mres), ycum, ymax, dy, ync,
     *                 x, unstch, strtch, delt, wt, nlin

      external unstch, strtch

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
c----------------------------------------------------------------------
      if (.not.extra) then
c                                 chopit always generates jsp coordinates
         ic = jsp
      else
c                                 but in the case of charge balance save
c                                 space for the dependent coordinate.
         ic = jsp + 1
      end if

      nlin = 0d0

      do i = 1, jsp

         k = jst + i
c                                 electrolyte model skip limit on ns'th
c                                 species
         if (ksmod(ids).eq.20.and.k.eq.ns) k = k + 1
c                                 generate coordinates for i'th component
         iy(i) = 1
         y(i,1) = pxmn(lpoly,lsite,k)

         ync = pxnc(lpoly,lsite,k)/wt

         if (ync.gt.0.5d0) ync = 0.5d0

         if (ync.eq.0d0) cycle

         mode = imdg(k,lsite,lpoly,ids)
c                                 avoid impossible compositions 'cause a min > 0
         if (i.gt.1) then

            ycum = ycum + pxmn(lpoly,lsite,k-1)
c                                 1-ycum is the smallest fraction possible
            if (ycum.gt.nopt(55)) then
c                                 inconsistent limits
               write (*,'(/,a,/)') '#########BOOM WACKA BOOM###########'
c              write (*,*) ycum,ids,ksmod(ids),lsite,k,i,mode
c              write (*,*) (pxmx(lpoly,1,j),j=1,jsp)
c              write (*,*) (pxmn(lpoly,1,j),j=1,jsp)
c              write (*,*) (pxnc(lpoly,1,j),j=1,jsp)
c              write (*,*) (iy(j),j=1,jsp)
c              call warn (999,ycum,jsp,'cartes')

               cycle

            else
c                                 the smallest fraction possible is lt
c                                 than xmax
               ymax = pxmx(lpoly,lsite,k)

            end if
         else
            ymax = pxmx(lpoly,lsite,k)
         end if
c                                 two means of extracting y-range, cartesian
c                                 imod = 0 and transformation imod = 1
         if (mode.eq.0) then
c                                 cartesian
            delt = nopt(5)

            do

               iy(i) = iy(i) + 1
               if (iy(i).gt.mres) call error (50,ync,mres,fname(ids))

               y(i,iy(i)) = y(i,iy(i)-1) + ync

               if (dabs(y(i,iy(i))-ymax).lt.delt.or.
     *             y(i,iy(i)).gt.ymax) then
                  y(i,iy(i)) = ymax
                  exit
               end if

            end do

         else
c                                 see 6.8.0 for old multiple interval conformal transformation
c                                 y is the non-linear cartesian coordinate
c                                 x is the linear conformal coordinate.
            nlin = nlin + 1d0

            call setstc (ids,lpoly,lsite,k)

            delt = xmno(ids,1,lsite,k)
            if (delt.gt.nopt(5)) delt = nopt(5)

            x = unstch (pxmn(lpoly,lsite,k))

            do

               iy(i) = iy(i) + 1
               if (iy(i).gt.mres) call error (50,ync,mres,fname(ids))

               x = x + ync
               y(i,iy(i)) = strtch (x)

               if (dabs(y(i,iy(i))-ymax).le.delt.or.
     *             y(i,iy(i)).gt.ymax) then

                  y(i,iy(i)) = ymax
                  exit

               end if

            end do

         end if

      end do
c                                 the first coordinate
      npairs = 1

      do i = 1, jsp
         ind(i) = 1
         simp(jump+i) = y(i,1)
      end do
c                                 now make the array index run over all
c                                 values increasing the last index fastest
      iexit = 0
      ieyit = 0
      dy = 0d0

      do while (iexit.eq.0)
c                                 figure out which index to increment
         do i = jsp, 1, -1

            if (ind(i).lt.iy(i).and.ieyit.eq.0) then
c                                 this is the one to increment
               ind(i) = ind(i) + 1
               indx = i
               exit

            else if (i.gt.1) then
c                                 saturated the index
               ind(i) = 1
               ieyit = 0

            else
c                                 saturated first index, done.
                return

            end if

         end do
c                                 ok now we have the indices, check
c                                 the composition
         ycum = 0d0
         dy = 0d0

         do i = 1, jsp
            ycum = ycum + y(i,ind(i))
         end do
c                                 until 2/17/19 this was > 1.
         if (ycum.gt.nopt(55)) then
c                                 no matter what this is the last point
c                                 to be generated for ind(indx), set ieyit
c                                 to change indx
            ieyit = 1
c                                 but here is where it gets messy:
            if (indx.eq.1) then
c                                 we're at the first point, and already
c                                 over the top
               iexit = 1
               cycle

            else if (y(indx,ind(indx)) - y(indx,ind(indx)-1)
     *               - ycum + 1d0    .gt. delt ) then
c                                 the excess (ycum-1) is less then the
c                                 amount the variable was previously incremented
c                                 so it's possible to back off the composition
c                                 to zero excess. this will always be true for
c                                 cartesian transformations, it might not be true
c                                 for conformal stretching (i.e., increments could
c                                 grow in the direction of the nodal index).
               dy =  1d0 - ycum

            else
c                                 must have just hit on the last increment or
c                                 conformal.
               cycle

            end if

         end if

         npairs = npairs + 1
         j = jump + (npairs-1)*ic

         if (j+jsp.gt.k13) then
            write (*,*) 'k13, k1 = ',k13,k1
            write (*,*) 'k21, k18, k20, k24, k25'
            write (*,*) k21, k18, k20, k24, k25
            call error (180,nlin,lsite,fname(ids))
         end if

         do i = 1, jsp
            simp(j+i) = y(i,ind(i))
         end do

         simp(j+indx) = simp(j+indx) + dy

      end do

      end

      subroutine cartaq (ids)
c---------------------------------------------------------------------
c subroutine to cartesian or transform subdivision on a single site
c solution with charge balance. called by subdiv.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, m, n, ids, qpairs, np0

      double precision ycum, sum, q, ratio

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c----------------------------------------------------------------------
c                                 could save only nqs-2 coordinates, but
c                                 will try nqs-1 first.
      ycum = 0d0

      if (ns1.eq.0) then
c                                 only solvent, test for no solvent has
c                                 already been made in reform
         do j = 1, nqs1
c           prism(j) = 0d0
         end do

         npairs = 1

      else
c                                 subdivision of neutral ns+nn-1 species
         call chopit (ycum,1d0,0,ns1,1,1,ids,0,.true.)

         do i = 1, npairs

            k = (i-1)*sn
            l = (i-1)*nqs1

            sum = 0d0

            do j = 1, ns1
c              prism(l+j) = simp(k+j)
               sum = sum + simp(k+j)
            end do

c           if (nq.gt.0) prism(l+nqs1) = sum

         end do

      end if

      ntot = npairs

      if (nq.ne.0) then

         np0 = npairs
c                                 do the nq-1 species independently
         ycum = 0d0

         call chopit (ycum,1d0,sn,nq1,1,1,ids,0,.true.)
c                                 at this point simp contains all
c                                 possible compositions of the nq-1 species,
c                                 use charge balance to get the nqth species
         qpairs = 1

         do i = 1, npairs

            q = 0d0
            sum = 0d0

            k = (i-1)*nq
            l = (qpairs - 1)*nq
            m = 2 + sn

            do j = 1, nq1
               q = q + thermo(6,jend(ids,m+j))*simp(k+j)
               sum = sum + simp(k+j)
               simp(l+j) = simp(k+j)
            end do
c                                 charge ratio
            ratio = q/thermo(6,jend(ids,m+j))
c                                 the net charge has the same sign as the nqth
c                                 species or its amount violates closure, reject:
            if (ratio.gt.0d0.or.sum-ratio.ge.1e0) cycle
c                                 the amount of the species determined by charge balance
            simp(l+nq) = -ratio

            qpairs = qpairs + 1

         end do

         qpairs = qpairs - 1
c                                 for every charged species composition
c                                 load all neutral compositions that don't
c                                 violate closure:
         do i = 1, qpairs
c                                 get the sum of the charged species
            k = (i-1)*nq

            sum = 0d0

            do j = 1, nq
               sum = sum + simp(k+j)
            end do
c                                 now assemble full compositions:
            do j = 1, np0

               l = (j-1)*nqs1
c                                 test for closure
c              if (prism(l+nqs1)+sum.ge.1d0) cycle
c                                 acceptable composition
               m = ntot * nqs1
               if (m+nqs1.gt.k24) call err41 ('K24')

               ntot = ntot + 1
c                                 load neutral part
               do n = 1, ns1
c                 prism(m+n) = prism(l+n)
               end do
c                                 load charged part
               do n = 1, nq
c                 prism(m+ns1+n) = simp(k+n)
               end do

            end do

         end do
c                                  zero the charged species coordinates
c                                  for the first np0 neutral compositions
         do i = 1, np0

            l = (i-1)*nqs1

            do j = sn, nqs1
c              prism(l+j) = 0d0
            end do

         end do

      end if

      end

      subroutine setsol (ids, wham)
c-----------------------------------------------------------------------
c load species indices, charges, etc for aqueous model (ksmod = 20) into
c simple arrays
c-----------------------------------------------------------------------
      implicit none

      logical wham

      include 'perplex_parameters.h'

      integer i, ids

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      character specie*4
      integer ins, isp
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)
c----------------------------------------------------------------------
      if (wham) then
c                                 an internal molecular eos has already
c                                 been invoked
         call error (72,rt,i,'only one solution model may invoke an '//
     *                       'internal molecular fluid EoS.')

      else

         wham = .true.

      end if
c                                 load endmember indices into a simple array
      do i = 1, mstot(ids)
         jnd(i) = jend(ids,2+i)
      end do
c                                 model uses an internal speciation routine
c                                 set the ins array and isp pointer
      if (ksmod(ids).eq.0) then
c                                 fluid eos specified via ifug
         call setins (ifug)

      else if (ksmod(ids).eq.20) then
c                                 electrolyte
         do i = 1, nqs
            q2(i) = thermo(6,jnd(i))**2
         end do

         isp = ns

         do i = 1, isp
            ins(i) = jspec(ids,i)
         end do

         na1 = 1

      else if (ksmod(ids).eq.39) then
c                                 hybrid molecular
         isp = mstot(ids)
c                                 set ns in case of aqrxdo or aqlagd
         ns = isp
         sn1 = ns + 1
         nsa = ns + aqct
         na1 = nsa + 1
         na2 = nsa + 2
         na3 = nsa + 3
         nat = nsa + 9

         do i = 1, isp
            ins(i) = jspec(ids,i)
         end do

      else if (ksmod(ids).eq.40) then
c                                 MRK silicate vapor (40), EoS code 26
         call setins (26)

      else if (ksmod(ids).eq.41) then
c                                 MRK COH fluid (41), EoS code 27
         call setins (27)

      end if

      end

      subroutine slvnt1 (gsolv)
c-----------------------------------------------------------------------
c computes solvent p-t-composition dependent properties: dielectric cst (eps),
c molar mass, debye-hueckel (adh), and gHKF function.
c
c assumes pure species volumes and fugacity coefficients have been calculated.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, k

      double precision gsolv, vsolv, cdh, ysum, ysolv(nsp), hyvol

      double precision gfunc, ghybrid, gcpd

      external gfunc, ghybrid, gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      double precision vol
      common/ cst26 /vol

      double precision yf,gmrk,v
      common/ cstcoh /yf(nsp),gmrk(nsp),v(nsp)

      double precision gh,dvhy,gmrk0
      common/ csthyb /gh(nsp),dvhy(nsp),gmrk0(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)

      save cdh
      data cdh/-42182668.74d0/
c----------------------------------------------------------------------
      gsolv  = 0d0
      msol   = 0d0
      ysum   = 0d0

      do i = 1, ns
c                                 solvent mass, kg/mol compound
         msol   = msol + pa(i) * fwt(jnd(i))
c                                 g mech mix term for solvent:
         gsolv  = gsolv + aqg(i) * pa(i)

         ysum = ysum + pa(i)

      end do
c                                 compute normalized fractions and correct solvent
c                                 gibbs energy for ideal solute concentration
      do i = 1, ns

         ysolv(i) = pa(i)/ysum

      end do
c                                 compute and add in solvent activities by calling
c                                 mrkmix ghybrid gets pmvs for hybrid volume calculation
      gsolv = gsolv + ysum*( ghybrid (ysolv) + rt*dlog(ysum) )
c                                 compute pmvs required for dielectric constant
      hyvol = 0d0

      do k = 1, ns

         i = ins(k)
c                                 hybrid pmv
         vhyb(i) = dvhy(i) + v(i)
c                                 hybrid total volume
         hyvol = hyvol + yf(i) * vhyb(i)

      end do

      do k = 1, ns

         i = ins(k)
c                                 volume fractions
         vf(i) = yf(i) * vhyb(i) / hyvol

      end do

      vsolv = ysum * hyvol
c                                 get dielectric cst based on hybrid volumetric properties
      call geteps (epsln)
c                                 set reference dielectric cst, quick fix for mistake
c                                 of having made it composition dependent,
      epsln0 = 78.47d0
c                                 Debye-Hueckel factor, A[cgs] = -q^3*sqrt(NA)/(4*Pi*k^(3/2))
c                                 *(msolg/(10*vsoljbar))^(1/2)/(epsilon*T)^(3/2) for ln(gamma) = +A*....
c                                 A = cdh*(msolkg/(vsoljbar))^(1/2)/(epsilon*T)^(3/2)
      adh = cdh * dsqrt(1d1*msol/vsolv/(epsln*t)**3)
c                                 shock et al 1992 g function (cgs solvent density),
c                                 used by hkf
      gf = gfunc (msol*1d3/vsolv)

      end

      subroutine slvnt2 (gsolv)
c-----------------------------------------------------------------------
c computes: debye-hueckel (adh) and solute contribution of the fluid
c gibbs energy. assumes the molar speciation is stored in pa. called
c only for solution model type 20.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k

      double precision gsolv, mo(m4), lng0, is

      double precision gcpd, aqact

      external gcpd, aqact

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh
c----------------------------------------------------------------------
      is = 0d0
c                                 molalities and ionic strength
      do k = sn1, nqs
c                                 ln molality of solutes
         mo(k) = pa(k)/msol
         is = is + q2(k) * mo(k)

      end do
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
c                                 Davies extension.
      lng0 = dlog(aqact(is/2d0))
c                                 add in the solute gibbs energies
      do k = sn1, nqs

         if (pa(k).le.0d0) cycle

         gsolv = gsolv + pa(k) *
     *           (gcpd(jnd(k),.true.) + rt*(dlog(mo(k)) + lng0*q2(k)))

      end do

      end

      subroutine solut0 (id)
c-----------------------------------------------------------------------
c computes solute endmember g's for aqueous solutions, used only for
c output purposes by routine calpr0.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id

      double precision gso(nsp), gcpd

      external gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)
c----------------------------------------------------------------------
c                                 solvent properties
      call slvnt3 (gso,.true.,.true.,id)
c                                 solute gibbs energies
      do i = sn1, nqs

         aqg(i) = gcpd(jnd(i),.true.)

      end do

      end

      subroutine geteps (epsln)
c-----------------------------------------------------------------------
c computes dielectric constant for the molecular species of a solvent
c assumes species volumes or partial molar volumes in cohhyb have been
c initialized.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j

      double precision po(nsp,11), trt, rho, epsln, eps

      double precision epsh2o

      external epsh2o

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)
c----------------------------------------------------------------------
c                                non-polar, A_mu = C_alpha = 0:

c                                Harvey & Lemmon provide additional data for
c                                ethylene and long-chain hydrocarbons.

c                                Eq 5 of H&L 2005:

c                                P/rho (cm3/mol) = A + A_mu/T + B*rho + C*rho^D

c                                is simplified by dropping the second term the
c                                coefficients are a f(T) viz A = a0 + a1*(T/Tr - 1)...

c                                po(i,1:8) - a0, a1, A_mu, b0, b1, c0, c1, D

c                                polar Harvey & Mountain 2017:

c                                rho *(C_alpha + g(rho,T)*A_mu(T)/T)

c                                expanded as [see dielectric_harvey.mws]:

c                                rho*(po(i,3)+po(i,2)*(po(i,1)*exp(po(i,4)*T^po(i,5))*
c                                (1-exp(po(i,6)*rho^po(i,7)))+1)*(po(i,8)+po(i,9)*
c                                exp(po(i,10)*rho^po(i,11)))^2/T)

      data ((po(i,j),j=1,11),i=1,nsp)/
c                                1 - H2O
     *     11*0d0,
c                                2 - CO2, H&L 2005
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0,
     *     1.55d0, 3*0d0,
c                                3 - CO approximated by O2
     *     3.9578d0, 6.5d-3, 0d0, 0.575d0, 1.028d0, -8.96d0, -5.15d0,
     *     1.5d0, 3*0d0,
c                                4 - CH4, H&L 2005
     *     6.5443d0, 1.33d-2, 0d0,8.4578d0, 3.7196d0, -352.97d0,
     *     -100.65d0, 2d0, 3*0d0,
c                                5 - H2, H&L 2005
     *     2.0306d0, 5.6d-3, 0d0, 0.181d0, 0.021d0, -7.4d0, 0d0, 2d0,
     *     3*0d0,
c                                6 - H2S, H&M 2017
     *     1.18d0, 5829.059676d0, 9.232464738d0, -.1213537391d-1, .9d0,
     *     -453374.7482d0, 3.5d0, 1.241d0, -.241d0, -16.61833221d0,
     *     .5d0,
c                                7 - O2, H&L 2005
     *     3.9578d0, 6.5d-3, 0d0, 0.575d0, 1.028d0, -8.96d0, -5.15d0,
     *     1.5d0, 3*0d0,
c                                8 - SO2, H&M 2017
     *     2.516d0, 16242.2847d0, 10.31715322d0, -.225289526d-2, .98d0,
     *     -44.03397284d0, 1.2d0, 1.335d0, .335d0, -16.19171204d0,
     *     .75d0,
c                                9 - COS approximated by CO2
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0,
     *     1.55d0, 3*0d0,
c                                10 - N2, H&L 2005
     *     4.3872d0, 2.26d-3, 0d0, 2.206d0, 1.135d0, -169d0, -35.83d0,
     *     2.1d0, 3*0d0,
c                                11 - NH3 and 12-15 Si-O high T species
     *     55*0d0,
c                                16 - Ethane, H&L 2005
     *     11.1552d0, 0.0112d0, 0d0, 36.759d0, 23.639d0, -808.03d0,
     *     -378.84d0, 1.75d0, 3*0d0,
c                                17 - dilutant
     *      11*0d0/

      save po
c----------------------------------------------------------------------
      trt = t/273.16d0 - 1d0

      epsln = 0d0

      do i = 1, ns - 1

         j = ins(i)
c                                 rho = 1/vcm3, v(j) is initialized by lnfpur
c                                 in gcpd and is in cm3/mol
         rho = 1d0/vhyb(j)
c                                 Eq 5 of H&L 2005 for  polarization/rho
c                                 for polar species need to add po(j,3)
         if (po(j,3).eq.0d0) then
c                                 invert clausius-mosotti relation for dielectric
c                                 constant (Eq 1) non-polar molecules
            eps = po(j,1) + po(j,2)*trt + (po(j,4) + po(j,5)*trt)*rho
     *                + (po(j,6) + po(j,7)*trt)*rho**po(j,8)

            eps = (2d0*eps*rho + 1d0) / (1d0 - rho*eps)

         else
c                                 polarized species H&L 2017
            eps = rho*(po(j,3) + po(j,2)
     *                         * (po(j,1)*dexp(po(j,4)*t**po(j,5))
     *                         * (1d0-dexp(po(j,6)*rho**po(j,7))) + 1d0)
     *                         * (po(j,8) + po(j,9)*
     *                            dexp(po(j,10)*rho**po(j,11)))**2/t)
c                                 invert Kirkwood relation
            eps = 2.25d0*eps + 0.25d0
     *            + dsqrt((5.0625d0*eps + 1.125d0)*eps+ .5625d0)

         end if
c                                  Looyenga mixing rule justified by Mountain &
c                                  Harvey 2015, modified here to use volume fraction
c                                  computed from partial molar volumes.
         epsln  = epsln  + vf(j) * eps**r13

      end do
c                                  add in water (ns^th species)
      epsln  = epsln + vf(ins(i)) * epsh2o (vhyb(ins(i))/1d1)**r13

      epsln  = epsln**3

      end

      subroutine slvntg (gso,mu)
c-----------------------------------------------------------------------
c given chemical potentials compute solvent species partial molar gibbs
c energies for aqrxdo/aqrxlg
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j

      double precision gso(nsp), mu(k8)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)
c-----------------------------------------------------------------------

         do i = 1, ns

            gso(i) = 0d0

            do j = 1, kbulk

               if (isnan(mu(j))) cycle

               gso(i) = gso(i) + mu(j)*cp(j,jnd(i))

            end do

         end do

      end

      subroutine aqrxdo (jd,lu)
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation

c   jd - is the pointer to the solvent in the local assemblage
c   lu - is the LUN for output, also signals type of output, console
c        vs tab (<0).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, ind(l9), lu, jd, badct

      logical bad, output

      character text*200

      double precision mo(l9), blk(k5), dn, smot, err, is, gamm0, totm,
     *                 g0(l9), lnkw, gso(nsp), ph0, tmass, tsmas, tsmol,
     *                 smol(k5), errkw, smo, gcpd, posox

      external gcpd

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      character cname*5
      common/ csta4  /cname(k5)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision fwt
      common/ cst338 /fwt(k10)

      double precision atwt
      common/ cst45 /atwt(k0)

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)

      save badct
      data badct/0/
c----------------------------------------------------------------------
      if (.not.mus) then

         call warn (99,0d0,0,
     *       'no chemical potentials, cannot back-calculate solute '//
     *       'speciation')

         return

      end if

      output = .true.

      if (jdaq.ne.0.and.lopt(32)) then
c                                 a multi species solvent is present:
         do i = 1, ns
c                                 doing lagged speciation, set the
c                                 solvent mole fractions to the true
c                                 values
            ysp(i,jd) = caq(jd,i)

         end do

      end if
c                                 set feos = .true. because can't be
c                                 sure that the last solvent calculation was
c                                 at the present p-t condition.
      call slvnt3 (gso,.true.,.true.,jd)
c                                 slvnt3 doing lots needless calculations
      call slvntg (gso,mu)
c                                 iterate on speciation
      call aqsolv (g0,gso,mo,mu,is,gamm0,lnkw,bad)
c                                 back calculated bulk composition
      if (bad) then

         if (badct.lt.11) then

            badct = badct + 1

            call warn (99,0d0,0,'AQRXDO did not converge on solute '//
     *                          'speciation')

            if (badct.eq.10) call warn (49,0d0,99,'AQRXDO')

         end if

         do i = 1, iprop
            prop(i) = nopt(7)
         end do

         return

      end if
c                                 compute charge balance error
      err = 0d0

      do i = 1, ichg
         err = err + q(jchg(i)) * mo(jchg(i))
      end do
c                                 neutral pH
      ph0 = -lnkw/2d0/2.302585d0

      do i = 1, kbulk
         blk(i) = 0d0
      end do
c                                 total molality
      smot = 0d0
c                                 compute mole fractions, total moles first
      do i = 1, ns
c                                 moles/kg-solvent
         do j = 1, kbulk
            blk(j) = blk(j) + pa(i)*cp(j,jnd(i))/msol
         end do

        smot = smot + pa(i)/msol

      end do

      smo = 0d0

      do i = 1, aqct
c                                  total solute molality
         smo = smo + mo(i)
         ind(i) = i

         do j = 1, kbulk
            blk(j) = blk(j) + mo(i)*aqcp(j,i)
         end do

      end do
c                                  net charge
      posox = 0d0

      if (oxchg) then
c                                check on charge imbalance
         do j = 1, kbulk
            posox = posox + cox(j)*blk(j)
         end do

      end if

      smot = smot + smo
c                                 bulk fluid composition
      tmass = 0d0
      totm = 0d0

      do i = 1, kbulk
         totm = totm + blk(i)
         tmass = tmass + atwt(i)*blk(i)
      end do
c                                error in log10(K_w)
      errkw = (-lnkw + dlog(mo(ihy)*mo(ioh)*gamm0**2
     *             ))/2.302585d0

      if (output.and.lu.lt.0) then
c                                 WERAMI props on a grid
         do i = 1, kbulk
c                                 bulk composition
            if (iopt(2).eq.0) then
c                                 molar
               prop(i) = blk(i)/totm*1d2

            else
c                                 mass
               prop(i) = blk(i)*atwt(i)/tmass*1d2

            end if

         end do

         k = kbulk

         do i = 1, ns
c                                 solvent speciation
            k = k + 1

            if (lopt(26)) then
c                                 mole fraction
               prop(k) = pa(i)/msol/smot
            else
c                                 molality
               prop(k) = pa(i)/msol
            end if

         end do

         do i = 1, aqct
c                                 solute speciation
            k = k + 1

            if (lopt(27)) then
c                                 molality
               prop(k) = mo(i)
            else
c                                 mole fraction
               prop(k) = mo(i)/smot
            end if

         end do
c                                  other properties
         prop(k+2) = -dlog10(mo(ihy)*gamm0)
         prop(k+1) = prop(k+2) - ph0
         prop(k+3) = dabs(errkw)
         prop(k+4) = epsln
         prop(k+5) = is
         prop(k+6) = smo

      else if (output) then
c                                 WERAMI/MEEMUM console output
         call rankem (mo,ind,aqct,iopt(32))

         if (jdaq.eq.20.or.lopt(32).and.jdaq.eq.39) then

            write (lu,1000)

         else

            write (lu,1005)

         end if

         write (text,1050) -dlog10(mo(ihy)*gamm0),dabs(errkw),
     *                     -dlog10(mo(ihy)*gamm0)-ph0,is,gamm0
         call deblnk (text)
         write (lu,'(400a)') chars(1:length)
         write (text,1070) epsln,msol*1d3,smo
         call deblnk (text)
         write (lu,'(400a)') chars(1:length)
         write (text,1075) tmass/smot, posox/smot
         call deblnk (text)
         write (lu,'(400a)') chars(1:length)

         if (jdaq.eq.20.or.lopt(32).and.jdaq.eq.39) then

            write (lu,1100)

         else

            write (lu,1040)

         end if

         do i = 1, iopt(32)

            k = ind(i)

            if (mo(k).eq.0d0) cycle

            if (jdaq.eq.20) then
c                                 compare to forward speciation:
c                                 check if the species is the solution model
               l = 0

               do j = sn1, nqs
                  if (jnd(j)-aqst.eq.k) then
                     l = j
                     exit
                  end if
               end do

               if (l.ne.0) then

                  write (lu,1080) aqnam(k),int(thermo(6,k+aqst)),
     *                        mo(k),mo(k)/smot,ysp(l,jd),
     *                        int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                        int(g0(k))
               else

                  write (lu,1090) aqnam(k),int(thermo(6,k+aqst)),
     *                        mo(k),mo(k)/smot,
     *                        int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                        int(g0(k))
               end if

            else if (lopt(32)) then
c                                 compare to lagged speciation:
                  write (lu,1080) aqnam(k),int(thermo(6,k+aqst)),
     *                        mo(k),mo(k)/smot,caq(jd,k+ns)/caq(jd,na2),
     *                        int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                        int(g0(k))

            else
c                                 only back-calculated result:
               write (lu,1010) aqnam(k),int(thermo(6,k+aqst)),
     *                         mo(k),mo(k)/smot,
     *                         int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                         int(g0(k))

            end if
         end do

         if (lopt(32)) then
c                               lagged
            write (lu,1140)

            do i = 1, ns

               write (lu,1150) names(jnd(i)), pa(i)/msol, 
     *                         pa(i)/msol/smot,
     *                         ysp(i,jd), vf(ins(i)), vhyb(ins(i)),
     *                         int(gso(i)), int(gcpd(jnd(i),.true.))

            end do

         else

            if (jdaq.eq.20) then
c                               direct
               write (lu,1020)
            else
c                               back-calculated, indicate normalized
               write (lu,1160)

            end if

            do i = 1, ns

               write (lu,1030) names(jnd(i)), ysp(i,jd)/msol, ysp(i,jd),
     *                         vf(ins(i)), vhyb(ins(i)),
     *                         int(gso(i)), int(gcpd(jnd(i),.true.))

            end do

         end if
c                                 footnotes:
         if (lopt(32).or.jdaq.eq.20) then

            write (lu,1170)

         else

            write (lu,1060)

         end if
c                                 bulk chemistry:
         if (jdaq.eq.20.or.lopt(32)) then

            tsmol = 0d0
            tsmas = 0d0

            do i = 1, kbulk
               smol(i) = 0d0
            end do

            if (lopt(32)) then
c                                 lagged model
               do i = 1, nsa

                  do j = 1, kbulk

                     if (i.lt.sn1) then
                        dn = pa(i)/caq(jd,na3) * cp(j,jnd(i))
                     else
                        dn = caq(jd,i) * aqcp(j,i-ns)
                     end if

                     smol(j) = smol(j) + dn

                     tsmol = tsmol + dn
                     tsmas = tsmas + dn*atwt(j)

                  end do

               end do

            else
c                                 forward model
               do i = 1, nqs

                  k = jnd(i)

                  do j = 1, kbulk

                     if (i.lt.sn1) then
                        dn = ysp(i,jd) * cp(j,k)
                     else
                        dn = ysp(i,jd) * aqcp(j,k-aqst)
                     end if

                     smol(j) = smol(j) + dn

                     tsmol = tsmol + dn
                     tsmas = tsmas + dn*atwt(j)

                  end do

               end do

            end if

            write (lu,1130)

            do i = 1, kbulk
               write (lu,1110) cname(i),blk(i)/totm*1d2,
     *                                  smol(i)/tsmol*1d2,
     *                                  blk(i)*atwt(i)/tmass*1d2,
     *                                  smol(i)*atwt(i)/tsmas*1d2
            end do

         else

            write (lu,1120)

            do i = 1, kbulk
               write (lu,1110) cname(i),blk(i)/totm*1d2,
     *                                  blk(i)*atwt(i)/tmass*1d2
            end do

         end if

      end if

1000  format (/,'Back-calculated vs optimized solute speciation:',/)
1005  format (/,'Simple back-calculated solute speciation:',/)
1010  format (a8,4x,i2,3x,g12.6,3x,g12.6,5x,i8,5(2x,g12.6))
1020  format (/,'Solvent endmember properties:',//,
     *        9x,'molality',3x,'mol_fraction',2x,'vol_fraction#',
     *        1x,'v,cm3/mol*',2x,'g,J/mol*',3x,'g0,J/mol***')
1030  format (a8,2x,f7.4,5x,f7.5,8x,f7.5,5x,f7.3,3x,i8,4x,i8)
1040  format (/,'Solute endmember properties:',//,10x,'charge',3x,
     *       'molality',5x,'mol_fraction',6x,'g,J/mol*',5x,'g0,J/mol**')
1050  format ('pH = ',f6.3,'+/-',f5.3,
     *        '; Delta_pH = ',f6.3,'; ionic_strength = ',
     *        g10.4,'; gamma/q^2 = ', g10.4)
1060  format (/,'*partial molar, **molal ref. state, ***molar ref. ',
     *        'state.',/)
1070  format ('dielectric constant = ', g10.4,
     *        '; solvent molar mass, g/mol = ',f8.4,
     *        '; solute molality = ',g10.4)
1075  format ('total molar mass, g/mol-species = ',f8.4,
     *        '; ref_chg, 1/mol-species = ',f8.4)
1080  format (a8,4x,i2,3x,g12.6,3x,g12.6,3x,g12.6,5x,i8,5(2x,g12.6))
1090  format (a8,4x,i2,3x,g12.6,3x,g12.6,20x,i8,5(2x,g12.6))
1100  format (/,'Solute endmember properties:',/,
     *        48x,'optimized',/,10x,'charge',3x,
     *       'molality',5x,'mol_fraction',3x,'mol_fraction',6x,
     *       'g,J/mol*',5x,'g0,J/mol**')
1110  format (1x,a8,2x,4(g12.6,3x))
1120  format (/,'Simple back-calculated fluid bulk composition:',//,
     *        13x,'mol %',11x,'wt %')
1130  format (/,'Back-calculated vs optimized fluid bulk composition:',
     *        //,26x,'optimized',21x,'optimized',/,
     *        13x,'mol %',10x,'mol %',10x,'wt %',11x,'wt %')
1140  format (/,'Solvent endmember properties:',//,
     *        9x,'molality',3x,'mol_fraction',2x,'opt_mol_frac',
     *        2x,'vol_fraction# v,cm3/mol*',2x,'g,J/mol*',
     *        3x,'g0,J/mol***')
1150  format (a8,2x,f7.4,5x,2(f7.5,8x),f7.5,5x,f7.3,3x,i8,4x,i8)
1160  format (/,'Normalized solvent endmember properties:',//,
     *        9x,'molality',3x,'mol_fraction',2x,'vol_fraction',
     *        1x,'v,cm3/mol*',2x,'g,J/mol*',3x,'g0,J/mol***')
1170  format (/,'#normalized, *partial molar, **molal ref. state, ',
     *        '***molar ref. state.',/)
      end

      subroutine rankem (a,ind,right,n)
c-----------------------------------------------------------------------
c rank the n largest values of array a(left:right) in array ind. assumes ind has
c been initialized (left:right).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, imax, ind(*), left, right, n

      double precision a(*), amax
c----------------------------------------------------------------------
      left = 1

      do

         amax = -1d99

         do i = left, right
            if (a(ind(i)).gt.amax) then
               imax = i
               amax = a(ind(i))
            end if
         end do

         j = ind(left)
         ind(left) = ind(imax)
         ind(imax) = j

         left = left + 1

         if (left.eq.n) return

      end do

      end

      double precision function solve (c,q,x,jchg,ichg,bad)
c-----------------------------------------------------------------------
c function solve for hydronium molality (x) from charge balance for aqrxdo.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jchg(*), ichg, it, i, j

      logical bad

      double precision c(*), q(*), x, y, z, f, df, dx
c-----------------------------------------------------------------------

      it = 0

      do

         f = 0d0
         df = 1d0
         it = it + 1

         do i = 1, ichg

            j = jchg(i)

            y = x**(q(j)) * c(j)
            z = y*q(j)/x

            f = f + y
            df = df + z

         end do

         dx = -f/df

         x = x + dx

         if (x.le.1d-50.or.x.gt.1d3.or.it.gt.iopt(21)) then
            bad = .true.
            exit
         else if (dabs(dx)/(1d0+x).lt.nopt(50)) then
            bad = .false.
            exit
         end if

      end do

      solve = x

      end

      subroutine aqidst
c-----------------------------------------------------------------------
c identify the aqueous phase for aqrxdo/aqlagd
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k

      logical lagged

      double precision tot

      character name*100

      character prject*100, tfname*100
      common/ cst228 /prject,tfname

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      integer eos
      common/ cst303 /eos(k10)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      character specie*4
      integer ins, isp
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
c                                 set option flags if necessary
      if (lopt(25).or.lopt(32)) then

         if (ifct.gt.0.and.(iff(1).ne.0.or.iff(2).ne.0)) then 

            call warn (99,0d0,0,'aq_output and aq_lagged_speciation'
     *            //'cannot be used with saturated phase components'//
     *              'and have been disabled (AQIDST)')

             aqct = 0
             iopt(32) = 0
             lopt(25) = .false.
             lopt(32) = .false. 

             return

         end if
c                                reset iopt(32) [# aq species output]
         if (iopt(32).gt.aqct) iopt(32) = aqct

      else

         aqct = 0
         iopt(32) = 0
         return

      end if

      jdaq = 0
      lagged = .false.
c                                 look among solutions:
      do i = 1, isoct

         if (ksmod(i).eq.20.or.ksmod(i).eq.39) then

            idaq = i
            jdaq = ksmod(i)

            if (lopt(32)) then

               do j = 1, ns
c                                 set quack flag so the pure endmembers
c                                 won't be speciated by aqlagd
                  quack(jnd(j)) = .true.

               end do
c                                 identify non-solvent components
               isolc = 0

               do j = 1, icp

                  tot = 0d0

                  do k = 1, ns
                     tot = tot + cp(j,jnd(k))
                  end do

                  if (tot.gt.0d0) cycle

                  isolc = isolc + 1
                  solc(isolc) = j

               end do

               lagged =.true.

            end if

          end if

      end do

      if (jdaq.eq.0) then
c                                 no solution model found:
c                                 turn off lagged speciation just to be sure
         lopt(32) = .false.

        if (.not.lopt(25)) aqct = 0

c                                 else look for H2O
         do i = 1, ipoint

            if (eos(i).eq.101) then
               idaq = -i
               jnd(1) = i
c                                 set solvent/species pointers on the
c                                 off chance they will be used
               ns = 1
               ins(1) = 1
               isp = 1
               return
            end if
         end do

      end if

      if (lagged) then

         if (.not.lopt(39).and.nrf(idaq)) then
c                                lagged speciation, set
c                                refine_endmembers to true.
         write (*,'(/,a)') '**error ver099** aq_lagged_speciation is T,'
     *                      //' but refine_endmembers is F (AQIDST).'
         write (*,'(a)') 'Set refine_endmembers in either '//
     *                   fname(idaq)//' or perplex_option.dat'
         call errpau

         end if 

      end if
c                                open a bad point file for lagged and
c                                back-calculated speciation calculations
      if (lagged.and.iam.le.2) then

         if (iam.eq.1) then
            call mertxt (name,prject,'.pts',0)
         else
            call mertxt (name,prject,'_MEEMUM.pts',0)
         end if
         open (n13,file=name)
c                                meemum/vertex
       else if (.not.lagged.and.iam.eq.3.and.lopt(25)) then
c                                werami back-calc
          call mertxt (name,prject,'_WERAMI.pts',0)
          open (n13,file=name)

      end if

      end

      subroutine slvnt0 (gsolv,vsolv)
c-----------------------------------------------------------------------
c sets solvent p-t dependent properties for pure water: dielectric cst (eps),
c molar mass, debye-hueckel (adh), and gHKF function.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision cdh, gsolv, vsolv

      double precision gfunc, gcpd, epsh2o, lnfpur

      external gfunc, gcpd, epsh2o, lnfpur

      double precision vol
      common/ cst26 /vol

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      double precision yf,gmrk,v
      common/ cstcoh /yf(nsp),gmrk(nsp),v(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)

      integer iam
      common/ cst4 /iam

      save cdh
      data cdh/-5661800.47810d0/
c-----------------------------------------------------------------------
      if (iam.ne.5) then
c                                 not frendly
         gsolv  = gcpd(jnd(1),.false.)

      else
c                                 frendly
         gsolv = lnfpur(101)

      end if
c                                 calling gcpd will get the molar volume (cm3/mol)
c                                 in cst26 as vol
      vsolv = vol
      msol   = 18.01528d-3
      epsln0 = 78.47d0

      epsln  = epsh2o (vol/1d1)
c                                 Debye-Hueckel factor, A[cgs] = -q^3*sqrt(NA)/(4*Pi*k^(3/2))
c                                 *(NH2O/(10*vh2o))^(1/2)/(epsilon*T)^(3/2) for ln(gamma) = +A*....
c                                 for reasons of stupidity this is set up for v in j/bar
      adh = cdh/dsqrt(vol/1d1*(epsln*t)**3)
c                                 shock et al 1992 g function (cgs solvent density),
c                                 used by hkf
      gf = gfunc (msol*1d3/vol)

      yf(1) = 1d0
      vf(1) = 1d0

      end

      subroutine outtit
c-----------------------------------------------------------------------
c outtit writes title information and a brief description of the
c chemical system for each calculation requested.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j

      character*162 title
      common/ csta8 /title(4)

      integer iasmbl
      common/ cst27  /iasmbl(j9)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character cname*5
      common/ csta4  /cname(k5)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ixct,ifact
      common/ cst37 /ixct,ifact

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iam
      common/ cst4 /iam

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c----------------------------------------------------------------------
      write (n3,1000)
c                          title:
      write (n3,1190) title(1)
c                          data base
      write (n3,1210) dname
c                          fluid
      if (ifct.gt.0.or.gflu) call rfluid (2)
c                          independent potentials:
      write (n3,1070) (vname(jv(j)), j = 1, ipot)
c                          saturated phase components:
      if (ifct.gt.0) then
         j = icp + isat
         write (n3,1200) (cname(j+i), i = 1, ifct)
      end if
c                          saturated components:
      if (isat.gt.0) then
         j = icp + isat
         write (n3,1180) (cname(i), i = icp1, j)
      end if
c                          unconstrained components
      write (n3,1080) (cname(i), i = 1, icp)
c                          phases
      if (iam.eq.15) then
c                          CONVEX
         if (icp.gt.3) then
            write (n3,1150) (cname(i), i = 1, icp)
            do i = istct, iphct
               write (n3,'(3x,a,12(1x,f6.3,1x))')
     *            names(i), (a(j,i)/ctot(i), j = 1, icp)
            end do
         else if (icp.eq.3) then
            write (n3,1090) (cname(i), i = 2, 3)
            write (n3,'(3(1x,a,1x,f6.3,1x,f6.3,5x))')
     *         (names(i), a(2,i)/ctot(i), a(3,i)/ctot(i),
     *                                             i = istct, iphct)
         else if (icp.eq.2) then
            write (n3,1040) cname(2)
            write (n3,'(4(2x,a,1x,f6.3))')
     *         (names(i), a(2,i)/ctot(i), i = istct, iphct)
         else if (icp.eq.1) then
            write (n3,1130)
            write (n3,'(7(1x,a,1x))') (names(i), i = istct, iphct)
         end if

      else 
c                          MEEMUM/VERTEX
         if (icp.gt.3) then
            write (n3,1150) (cname(i), i = 1, icp)
            do i = istct, ipoint
               write (n3,'(3x,a,12(1x,f6.3,1x))')
     *            names(i), (cp(j,i)/ctot(i), j = 1, icp)
            end do
         else if (icp.eq.3) then
            write (n3,1090) (cname(i), i = 2, 3)
            write (n3,'(3(1x,a,1x,f5.3,1x,f6.3,5x))')
     *         (names(i), cp(2,i)/ctot(i), cp(3,i)/ctot(i),
     *                                             i = istct, ipoint)
         else if (icp.eq.2) then
            write (n3,1040) cname(2)
            write (n3,'(4(2x,a,1x,f6.3))')
     *         (names(i), cp(2,i)/ctot(i), i = istct, ipoint)
         else if (icp.eq.1) then
            write (n3,1130)
            write (n3,'(7(1x,a,1x))') (names(i), i = istct, ipoint)
         end if

      end if 
c                          saturation composant phases
      if (isat.ne.0) write (n3,'(/,a,/)')
     *          'Phases on saturation and buffering surfaces:'

      do i = 1, isat

         write (n3,'(/,3a,/)') ' for component ', cname(i+icp),':'
         write (n3,'(7(1x,a,1x))') (names(ids(i,j)), j = 1, isct(i))

      end do
c                          excluded phases
      if (ixct.ne.0) then
         write (n3,1140)
         write (n3,'(7(1x,a,1x))') (exname(i), i = 1, ixct)
      end if
c                          solution models
      if (isoct.ne.0) then 
         write (n3,1140)
         write (n3,'(6(1x,a,1x))') (fname(i), i = 1, isoct)
      end if

      write (n3,1000)

1000  format (/,80('-'),/)
1030  format (/,'Solution models considered:',/)
1040  format (/,'Phases and (projected) mol fraction ',a,':',/)
1070  format (/,'Independently constrained potentials:',//,3x,8(a,1x))
1080  format (/,'Components with unconstrained potent'
     *       ,'ials:',//,3x,10(a5,3x))
1090  format (/,'Phases and (projected) composition with respect to '
     *         ,a5,' and ',a5,':',/)
1130  format (/,'Phases:',/)
1140  format (/,'Excluded phases:',/)
1150  format (/,'Phases and (projected) compositions:',//,
     *        11x,12(1x,a5,2x),/)
1180  format (/,'Saturated or buffered components:',//,3x,7(a,3x))
1190  format (/,'Problem title: ',a,/)
1200  format (/,'Saturated phase components:',//,3x,5(a,3x))
1210  format ('Thermodynamic data base from: ',a)
      end

      recursive double precision function gphase (id)
c-----------------------------------------------------------------------
c gphase computes the gibbs free energy of static compound id.
c gphase does not assume that the properties of a pseudocompound endmembers
c are known (unlike gall) it is thus less efficient than gall.

c used only (i think) for mixed variable and schreinemaker's projections.
c alloy solution models are commented out. these must be reinstated for the
c aforementioned diagram types.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer id, ids

      double precision gzero, gerk, gph, gproj, gcpd, gfesi, gex,
     *                gfecr1, gfesic, gfes, gexces, gmchpr, gmech0

      external gzero, gerk, gproj, gcpd, gfesi, gfecr1, gfesic, gfes, 
     *         gex, gexces, gmchpr, gmech0

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer ikp
      common/ cst61 /ikp(k1)

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c----------------------------------------------------------------------
      ids = ikp(id)

      if (id.le.ipoint) then
c                                 a pure compound
         gph = gcpd (id,.true.)

      else if (lorder(ids)) then
c                                 get composition
         call setxyp (ids,id,bad)
c                                 compute margules coefficients
         call setw (ids)
c                                 evaluate enthalpies of ordering
         call oenth (ids)

         if (.not.noder(ids)) then

            call specis (gph,ids)

         else
c                                 no derivatives
            call minfxc (gph,ids,.false.)

         end if
c                                 add dqf corrections and mech
         gph = gph + gmchpr(ids) + gexces (id)

      else if (ksmod(ids).eq.0) then

         call setxyp (ids,id,bad)
c                              get the excess and/or ideal mixing effect
c                              and/or dqf corrections:
         call fexces (id,gph)

         gph = gph + gmech0(ids)

      else if (ksmod(ids).eq.40) then
c                                 si-o mrk fluid
         call setxyp (ids,id,bad)

         gph = gmech0(ids) + gerk (pa)

      else if (ksmod(ids).ge.29.and.ksmod(ids).le.32) then
c                                 nastia's models:
          call setxyp (ids,id,bad)

          if (ksmod(ids).eq.29) then
c                                 BCC Fe-Si Lacaze and Sundman
             gph = gfesi(pa(1), gproj (jend(ids,3)),
     *                         gproj (jend(ids,4)) )

          else if (ksmod(ids).eq.32) then
c                                 BCC Fe-Cr Andersson and Sundman
             gph = gfecr1(pa(1), gproj (jend(ids,3)),
     *                          gproj (jend(ids,4)) )

          else
c                                 Nastia's version of BCC/FCC Fe-Si-C Lacaze and Sundman
c                                 this model has to be called ahead of the standard models
c                                 because it sets lrecip(id) = true.
            gph =  gfesic (pa(1),pa(3),pa(4),
     *                     gproj (jend(ids,3)), gproj (jend(ids,4)),
     *                     gproj (jend(ids,5)), gproj (jend(ids,6)),
     *                     ksmod(ids))

          end if

      else if (ksmod(ids).eq.42) then

         call setxyp (ids,id,bad)
c                                 Fe-S fluid (Saxena & Eriksson 2015)
         gph = gfes (pa(2), gproj (jend(ids,3)), gproj (jend(ids,4)) )

      else
c                                 normal models (configurational
c                                 entropy fixed, excess function
c                                 linear in p-t) and special models
c                                 with normal gmech term
         call setxyp (ids,id,bad)

         if (ksmod(ids).eq.41) then
c                                 ternary coh fluid deltag
            call rkcoh6 (pa(2),pa(1),gph)

         else if (ksmod(ids).eq.26) then

            call hcneos (gph,pa(1),pa(2),pa(3))

         else
c                                 get the excess and/or ideal mixing effect
c                                 and/or dqf corrections:
             gph = gexces (id)

         end if
c                                 add gmech
         gph = gph + gmchpr(ids)
c                                 for van laar get fancier excess function
         if (llaar(ids)) then

            call setw (ids)

            gph = gph + gex (ids,pa)

         end if

      end if

      gphase = gph

      end

      subroutine factr1 (n,ier)
c-----------------------------------------------------------------------
c factr1 is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the dimension of the matrix a.
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c         ier- a flag, zero if a is of rank = n, and 1 if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,n,ip1,istr,ier

      double precision temp,ratio,tmax,rmax

      integer ipvt
      double precision a,d,x(k5)
      common/ cst301 /a(k5,k5),d(k5),ipvt(k5)
c-----------------------------------------------------------------------
      ier = 0
c                            initialize ipvt,d
      do i = 1, n
         ipvt(i) = i
         rmax = 0d0
         do j = 1,n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.1d-5) goto 9000
         x(i) = rmax
      end do
c                            begin decomposition:
      do i = 1, n - 1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/x(i)
         istr = i
         ip1 = i + 1

         do j = ip1, n
            tmax = dabs(a(j,i))/x(j)
            if (tmax.le.rmax) cycle
            rmax = tmax
            istr = j
         end do

         if (dabs(rmax).lt.1d-5) goto 9000
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then
            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = x(istr)
            x(istr) = x(i)
            x(i) = temp
            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do
         end if
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1,n
            a(j,i) = a(j,i)/a(i,i)
            ratio = a(j,i)
            do k = ip1, n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do
         end do

      end do

      if (dabs(a(n,n)).lt.1d-5) ier = 1

      return
c                           algoritmic singularity.
9000  ier = 1

      end

      subroutine uproj
c----------------------------------------------------------------------
c subroutine uproj computes the potentials of saturated phase components
c and saturated components.

c the energies of saturated components are projected through
c saturated volatile components.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,ict,ll,i1,id

      double precision uss(h6), fo2, gph, u, gphase

      external gphase

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k10)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      fo2 = 0d0
c                                 compute the chemical potentials of
c                                 saturated phase components.
      if (ifct.gt.0) call ufluid (fo2)

      do i = 1, isat
c                                 determine stable saturated composants
c                                 and the corresponding chemical potentials
         ict = isct(i)

         ll = icp+i

         do j = 1, ict

            k = ids(i,j)
            gph = gphase (k)

            if (ifct.gt.0) then
               do l = 1, 2
c                                 legendre transform for saturated phase
c                                 component potentials
                  if (iff(l).ne.0) gph = gph - cp(iff(l),k)*uf(l)
               end do
            end if

            uss(j) = gph

            if (i.gt.1) then
c                                 if multiple component saturation constraints
c                                 apply saturation hierarchy legendre transform:
               i1 = i-1
               do l = 1, i1
                  uss(j) = uss(j)-cp(icp+l,k)*mu(icp+l)
               end do
            end if

            g(k) = uss(j)
            uss(j) = uss(j)/cp(ll,k)
         end do
c                                 if O2, check if fo2 has been
c                                 determined by a fluid phase routine,
c                                 if so, add the transform:
         if (io2.eq.i) then
            do j = 1, ict
               uss(j) = uss(j) + r*t*fo2
            end do
         end if
c                           now find stable "composant":

         u = uss(1)

         id = 1

         if (ict.ne.1) then
            do j = 2, ict
               if (uss(j).gt.u) cycle
               id = j
               u = uss(j)
            end do
         end if
c                               save the id of the stable composant.
         idss(i) = ids(i,id)
c                               and its chemical potential.
         mu(icp+i) = u
c                               in case a phase in the component
c                               saturation space is an endmember of
c                               a solution transform the endmember G's:
         do j = 1, ict
            k = ids(i,j)
            g(k) = g(k) - cp(icp+i,k)*u
         end do

      end do

      end

      subroutine gall
c-----------------------------------------------------------------------
c subroutine gall computes molar free energies of all static compounds.
c for non-equimolar o/d returns the p0 normalized g.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer i, j, k, id

      double precision gval, dg, g0(m14)

      double precision gex, gfesi, gfesic, gerk, gproj, ghybrid, gzero,
     *                 gfecr1, gcpd, gfes, gmech, gexces

      external gerk, gzero, gex, gfesi, gfesic, gproj, ghybrid, gexces,
     *         gcpd, gfes, gmech

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision g
      common/ cst2 /g(k1)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision p,t,xco2,mu1,mu2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,mu1,mu2,tr,pr,r,ps

      double precision mmu
      common/ cst39 /mmu(i6)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
c-----------------------------------------------------------------------
c                                 compute the chemical potential
c                                 of the projected components. 
c                                 saturated component potentials are loaded
c                                 into mu(icp+1:icp+isat), currently no provision
c                                 for saturated phase components, this would be
c                                 necessary to allow back-calculated electrolytic
c                                 fluid speciation with a saturated fluid phase
      call uproj
c                                 load mobile components into mu(icp+isat+1:kbulk)
      do j = 1, jmct

         mu(icp+isat+j) = mmu(j)

      end do
c                                 first do the endmembers:
c                                 changed start index to 1 from kphct + 1, to
c                                 allow projection through endmembers in constrained
c                                 component space... feb 3, 2019 this seems to be
c                                 necessary for mobile components, but if so, why did
c                                 ah2o calculations work before??
      do id = 1, ipoint

         g(id) = gproj (id)

      end do
c                                 now do solutions:
      do i = 1, isoct

         if (lname(i).eq.'liquid'.and.t.lt.nopt(20)) then
c                                 a liquid below T_melt option threshold
            do j = 1, jend(i,2)

               g(id) = 1d6

               id = id + 1

            end do

         else if (lorder(i)) then
c                                 initialize margules, enthalpy of
c                                 ordering, internal dqfs (last for minfxc)
            call ingsol (i)
c                                 only for minfxc
            call ingend (i)

            do j = 1, jend(i,2)

               call setxyp (i,id,bad)

               if (.not.noder(i)) then

                  call specis (dg,i)

               else

                  call minfxc (dg,i,.false.)

               end if
c                                 for non-equimolar o/d gexces and gmech
c                                 are computed for the p0 mass, this is 
c                                 ok, because the g will be normalized by 
c                                 the static ctot value.
               g(id) = gexces(id) + dg + gmech(i)

               id = id + 1

            end do

         else if (.not.llaar(i).and.simple(i)) then
c                                 it's normal margules or ideal:
            do j = 1, jend(i,2)

               call setxyp (i,id,bad)
c                                 gexces returns with excess energy, dqf,
c                                 and configurational entropy terms for 
c                                 simple models.
               g(id) = gexces(id) + gmech(i)

               id = id + 1

            end do

         else if (ksmod(i).eq.0) then
c                                 it's a fluid compound, the way
c                                 things are now it must have two
c                                 components.
            do j = 1, lstot(i)
               g0(j) = gzero(jend(i,2+j))
            end do

            do j = 1, jend(i,2)

               call setxyp (i,id,bad)

               call fexces (id,gval)

               g(id) = g0(1) * pa(1) + g0(2) * pa(2) + gval

               id = id + 1

            end do

         else if (llaar(i)) then
c                                 compute margules coefficients
            call setw (i)
c                                 because the hp van laar may have p-t
c                                 dependent volumes, the full expression
c                                 must be evaluated here by gex:
            do j = 1, jend(i,2)

               call setxyp (i,id,bad)
c                                 for simple vanlaar gexces returns dqf
c                                 and configurational entropy terms
               g(id) = gexces(id) + gex(i,pa) + gmech(i)

               id = id + 1

            end do

         else if (ksmod(i).eq.20) then
c                                 electrolytic solution, assumes:
c                                 1) molal electrolyte standard state
c                                 for solutes.
c                                 2) water is the last species
c                                 solvent species Gibbs energies:
c                                 solvent Gibbs energies
            rt = r*t

            do k = 1, ns
               aqg(k) = g(jnd(k))
            end do
c                                 compute compound properties
            do j = 1, jend(i,2)
c                                 get the composition
               call setxyp (i,id,bad)
c                                 solvent properties
               call slvnt1 (g(id))
c                                 add in solute properties
               call slvnt2 (g(id))

               id = id + 1

            end do

         else if (ksmod(i).eq.26) then
c                                 H2O-CO2-Salt:
            do j = 1, jend(i,2)

               call setxyp (i,id,bad)

               call hcneos (g(id),pa(1),pa(2),pa(3))

               g(id) = g(id) + gmech(i)

               id = id + 1

            end do

         else if (ksmod(i).eq.39) then
c                                 generic hybrid EoS
            do j = 1, jend(i,2)
c                                 load composition array and pointers
               call setxyp (i,id,bad)
c                                 compute and add in activities
               g(id) = ghybrid(pa) + gmech(i)

               id = id + 1

            end do

         else if (ksmod(i).ge.29.and.ksmod(i).le.32) then
c                                 nastia's models:
            do j = 1, jend(i,2)
c                                 load composition array and pointers
               call setxyp (i,id,bad)

               if (ksmod(i).eq.29) then
c                                 BCC Fe-Si Lacaze and Sundman
                  g(id) = gfesi(pa(1),g(jend(i,3)),g(jend(i,4)))
               else if (ksmod(i).eq.32) then
c                                 BCC Fe-Cr Andersson and Sundman
                  g(id) = gfecr1 (pa(1),g(jend(i,3)),g(jend(i,4)))
               else

                  g(id) = gfesic (pa(1),pa(3),pa(4),
     *                            g(jend(i,3)),g(jend(i,4)),
     *                            g(jend(i,5)),g(jend(i,6)),ksmod(i))
               end if

               id = id + 1

            end do

         else if (ksmod(i).eq.41) then

            do j = 1, jend(i,2)
c                                 hybrid MRK ternary COH fluid
               call setxyp (i,id,bad)

               call rkcoh6 (pa(2),pa(1),g(id))

               g(id) = g(id) + gmech(i)

               id = id + 1

            end do

         else if (ksmod(i).eq.40) then

            do j = 1, jend(i,2)
c                                 MRK silicate vapor
               call setxyp (i,id,bad)

               g(id) = gmech (i) + gerk (pa)

               id = id + 1

            end do

         else if (ksmod(i).eq.42) then
c                                 Fe-S fluid (Saxena & Eriksson 2015)
            do j = 1, jend(i,2)

               call setxyp (i,id,bad)

               g(id) = gfes (1d0-pa(1),g(jend(i,3)),g(jend(i,4)))

               id = id + 1

            end do

         end if

      end do

      end

      subroutine ufluid (fo2)
c----------------------------------------------------------------------
c subroutine ufluid computes the potential of the components
c of a saturated fluid phase. if the mole fraction of a component is les
c less than 1.d-38 the chemical potential is set to -9.9d09.
c ufluid may call one of three molecular fluid equations of state, or
c alternatively users may supply their own routines, however,
c the routines currently in use return the log of a components fugacity
c which is then added to the reference state potential computed by the
c function gphase.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision xf(2), fo2, fs2, gcpd, gzero

      external gcpd, gzero

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)

      integer ifct,idfl
      common/ cst208 /ifct,idfl
c-----------------------------------------------------------------------
c                           compute the chemical potentials of
c                           fluid components in fluid saturated
c                           systems.
      call cfluid (fo2,fs2)

      if (idfl.ne.0) then

         uf(idfl) = gcpd (idfl,.false.) + r * t * f(idfl)

      else

         xf(1) = 1d0 - xco2
         xf(2) = xco2

         do i = 1, 2

            if (iff(i).ne.0) then

               if (xf(i).lt.1d-38) then

                  uf(i) = -1d10

               else

                  uf(i) = gzero (i) + r * t * f(i)

               end if

            end if

         end do

      end if

      end

      double precision function gzero (id)
c----------------------------------------------------------------------
c gzero computes the 1 bar reference pressure free energy of a compound
c identified by the argument 'id'. no fugacity terms are added for real
c gas species (cf, gcpd).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, j

      double precision g, vdp

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision mmu
      common/ cst39 /mmu(i6)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer eos
      common/ cst303 /eos(k10)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)
c----------------------------------------------------------------------

      g = thermo(1,id)
c                                 -sdt
     *      + t * (thermo(2,id) - thermo(4,id) * dlog(t)
     *      - t * (thermo(5,id) + (thermo(7,id) - thermo(24,id)*t) * t))
     *      - (thermo(6,id) + thermo(10,id) / t) / t
     *      + thermo(8,id) * dsqrt(t) + thermo(9,id)*dlog(t)

      do j = 1, jmct
c                                 -ndu
         g = g - vnumu(j,id) * mmu(j)
      end do
c                                 transitions
      vdp = 0d0

      if (lct(id).ne.0) call mtrans (g,vdp,id)

      gzero = g

      end

      subroutine lambw (dg,ld)
c---------------------------------------------------------------------
c calculate the energy of an order-disorder transition using the
c Bragg-Williams model (Holland and Powell, '96), 0-d speciation.
c    input: ld - pointer to the phase in therlm
c   output: dg - energy change of ordering
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ld

      double precision dg,h,w

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
c                                 enthalpy of complete disordering
      h = therlm(1,1,ld) + therlm(2,1,ld)*p
c                                 interaction energy
      w = therlm(3,1,ld) + therlm(4,1,ld)*p

      call speci0 (dg,h,w,therlm(5,1,ld),therlm(6,1,ld),
     *                    therlm(7,1,ld),therlm(8,1,ld))

      end

      subroutine outlim
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option. write arf file for convex
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer ii, i, j, k, ibad1, ibad2, ibad3, igood, ipop

      logical bad1, bad2, good

      double precision num, mnsum, mxsum
c                                 -------------------------------------
      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer iam
      common/ cst4 /iam

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)
c----------------------------------------------------------------------
      ibad1 = 0
      ibad2 = 0
      ibad3 = 0
      igood = 0

      if (lopt(11)) rewind (n11)

      if (isoct.eq.0) goto 99

      bad1 = .false.
      bad2 = .false.
      good = .false.

      do i = 1, isoct

         if (.not.stable(i)) then
            bad1 = .true.
            ibad1 = ibad1 + 1
         else
            good = .true.
            igood = igood + 1
         end if

         if (limit(i)) then
            bad2 = .true.
            ibad2 = ibad2 + 1
         end if

      end do

      if (.not.refine.and.iam.eq.15) then 
         rewind (n10)
         write (n10,*) ibad1,0,igood
      end if
c                                 write solutions present that are
c                                 not stable
      if (bad1) then

         write (*,1000)
         if (lopt(11)) write (n11,1000)

         do i = 1, isoct
            if (.not.stable(i)) then
               write (*,'(5x,a)') fname(i)
               if (.not.refine.and.iam.ne.1) write (n10,'(a)') fname(i)
               if (lopt(11)) write (n11,'(5x,a)') fname(i)
            end if
         end do
      end if

      if (.not.good) goto 99
c                                 write solutions that are on an internal
c                                 limit
      if (bad2) then

         if (icopt.gt.3) then
c                                 adaptive minimization
c                                 solutions whose limits could be relaxed
            write (*,1080)
            if (lopt(11)) write (n11,1080)

         else
c                                 non-adaptive minimization,
c                                 solutions on internal limits
            write (*,1010)
            if (lopt(11)) write (n11,1010)

         end if

         do i = 1, isoct
            if (limit(i)) then
               write (*,'(5x,a)') fname(i)
               if (lopt(11)) write (n11,'(5x,a)') fname(i)
            end if
         end do

         if (refine) then
            write (*,1091)
            if (lopt(11)) write (n11,1091)
         else
            write (*,1090)
            if (lopt(11)) write (n11,1090)
         end if

      end if

      do i = 1, isoct

         if (.not.stable(i)) cycle

         ipop = pop1(i)

         if (ipop.gt.1) then
c                                 composite composition space, check
c                                 for 0-wt subcompositions
            mnsum = 0d0
            mxsum = 0d0

            do k = 1, ndim(1,ipop,i)
               if (xlo(k,1,ipop,i).gt.xhi(k,1,ipop,i)) then 
                  xlo(k,1,ipop,i) = 0d0
                  xhi(k,1,ipop,i) = 0d0
               else
                  mnsum = mnsum + xlo(k,1,ipop,i)
                  mxsum = mxsum + xhi(k,1,ipop,i)
               end if
            end do

            if (xlo(k,1,ipop,i).gt.xhi(k,1,ipop,i)) then
               xlo(k,1,ipop,i) = 1d0 - mxsum
               xhi(k,1,ipop,i) = 1d0 - mnsum
            end if

         end if

         if (.not.refine.and.iam.eq.15) then

            write (n10,'(a)') fname(i)

            do ii = 1, ipop
               do j = 1, istg(i,ii)
                  do k = 1, ndim(j,ii,i)
                     write (n10,*) xlo(k,j,ii,i),xhi(k,j,ii,i)
                  end do
               end do
            end do

         end if
c                                 special case (1 component solution).
         if (ksmod(i).eq.39.and.ndim(1,1,i).eq.0) cycle

         call limprt (6,i)
         if (lopt(11)) call limprt (n11,i)

      end do

99    if (goodc(1)+badc(1).gt.0d0) then

         num = badc(1)/(badc(1)+goodc(1))*1d2
         write (*,1120) num, badc(1) + goodc(1)
         if (num.gt.1d-1) call warn (53,num,i,'OUTLIM')
         write (*,1140) goodc(2)/(badc(1)+goodc(1))

      end if

      if (iam.eq.15) close (n10)
      if (lopt(11)) close (n11)

1000  format (/,'The following solutions were input,'
     *         ,' but are not stable:',/)
1010  format (/,'**warning ver991** The following solutions have ',
     *          'compositions at an internal limit (i.e., 0<x<1):',/)
1080  format (/,'**warning ver991** The compositions of the following',
     *        ' solutions reached internal',/,
     *        'limits that were automatically relaxed:',/)
1090  format (/,'If the restrictions are unintentional, then relax ',
     *          'the corresponding limits',/,'in the solution model ',
     *          'file and restart the calculation.',/)
1091  format (/,'Restriction during the auto-refine stage is usually ',
     *          'unimportant. If desired, confirm',/,'by ',
     *          'comparing the ranges ',
     *          'below to those in the *.arf file.',//,'NOTE: ',
     *          'unintentional restrictions encountered during the ',
     *          'exploratory stage may be',/,'problematic, refer to ',
     *          'the *_auto_refine.txt file ',
     *          'for the exploratory stage warnings.',/)
1110  format (4x,a,t35,i2)
1120  format (/,'The failure rate during speciation (order-disorder) ',
     *        'calculations is ',f7.3,'%',/,'out of a total of ',f12.0,
     *        ' calculations.',/)
1140  format (/,'Average number of iterations per speciation ',
     *          'calculation:',f5.1,/)
      end

      subroutine outarf
c----------------------------------------------------------------------
c subroutine output arf file for vertex.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, tmp
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c----------------------------------------------------------------------
      rewind (n10)

      if (.not.refine.or.lopt(55)) then
c                                 load the former dynamic compositions
c                                 into the static arrays
         call reload (.false.)
c                                 output to arf
         write (n10,*) isoct
         write (n10,'(7(a,1x))') fname(1:isoct)
         write (n10,*) jend(1:isoct,2)

         tcct = 0

         do i = 1, isoct

            tmp = jend(i,2)*tstot(i)
            write (n10,*) txco(tcct+1:tcct+tmp)
            tcct = tcct + tmp

         end do

      end if

      close (n10)

      end


      subroutine limprt (lun,i)
c----------------------------------------------------------------------
c subroutine to print limits to LUN for outlim.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ii, i, j, k, ipop, lun

      character char8*8
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)
c----------------------------------------------------------------------
      ipop = pop1(i)

      if (istg(i,1).eq.1.and.ipop.eq.1.and.ksmod(i).ne.688) then
c                                 single site solution
         write (lun,1020) fname(i)

         if (ksmod(i).eq.20) then
c                                 charge balance model:
            do j = 1, ispg(i,1,1)

               if (j.lt.ns) then
                  char8 = names(jnd(j))
               else if (j.eq.ns) then
                  cycle
               else
                  char8 = aqnam(jnd(j)  - aqst)
               end if

               write (lun,1030) char8, xlo(j,1,1,i), xhi(j,1,1,i)

            end do

         else

            do j = 1, ispg(i,1,1)

               write (lun,1030) names(jend(i,2+j)),
     *                              xlo(j,1,1,i),xhi(j,1,1,i)

            end do

         end if

      else if (ipop.eq.1) then 
c                                 single polytope
         if (istg(i,1).gt.1) then
            write (lun,1040) 'prismatic model: '//fname(i)
         else 
            write (lun,1020) fname(i)
         end if

         do j = 1, istg(i,1)

            if (istg(i,1).gt.1) write (lun,1050) j,' '

            if (ispg(i,1,j).eq.1) then

               write (lun,1060)

            else

               do k = 1, ispg(i,1,j)

                  if (ksmod(i).ne.688) then 
                     write (lun,1070) k,xlo(k,j,1,i),xhi(k,j,1,i)
                  else 
                     write (lun,1080) poname(i,1,j,k),xlo(k,j,1,i),
     *                                                xhi(k,j,1,i)
                  end if

               end do

            end if

         end do

      else
c                                 composite polytope
         write (lun,1160) 'composite composition space model: '
     *                    //fname(i)
c                                  polytope weights:
         do ii = 1, poly(i)

            write (lun,1170) poname(i,pop1(i),1,ii), xlo(ii,1,ipop,i), 
     *                                     xhi(ii,1,ipop,i)

         end do
c                                  individual polytope compositions
         do ii = 1, poly(i)

c                                  polytope
               write (lun,'(/,a)') ' '//poname(i,pop1(i),1,ii)//
     *                             ' Subcomposition:'

               do j = 1, istg(i,ii)

                  if (ispg(i,ii,j).eq.1) then
                     write (lun,'(/,a)') ' Subcomposition '//
     *                     poname(i,pop1(i),1,ii)//' is 0-dimensional'
                     cycle
                  end if 

                  write (lun,1050) j,' '

                  do k = 1, ispg(i,ii,j)
                     write (lun,1080) poname(i,ii,j,k),xlo(k,j,ii,i),
     *                                                 xhi(k,j,ii,i)
                  end do

               end do

         end do

      end if
c                                 warn about bad p2yx inversions
      if (badinv(i,1).gt.0) write (lun,1010) float(badinv(i,1)) /
     *                      float(badinv(i,1)+badinv(i,2))*1d2,fname(i)

1010  format (/,'**warning ver204** ',f5.1,'% of the p2y inversions fo',
     *       'r ',a,' failed. A high',/,'failure rate may indicate th',
     *       'at the compositional polyhedron for the model does',/,
     *       'not span all possible model compositions.',/)
1020  format (/,'Compositions for simplicial model: ',a,//,5x,
     *          '              Minimum         Maximum')
1030  format (5x,a8,4x,g12.5,4x,g12.5)
1040  format (/,'Compositions for ',a)
1050  format (/,'  Simplex ',i1,a,/,5x,
     *               '              Minimum         Maximum')
1060  format (8x,'Dummy site generated by model reformulation',/)
1070  format (8x,i2,7x,g12.5,4x,g12.5)
1080  format (5x,a10,2x,g12.5,4x,g12.5)
1160  format (/,a,//,
     *          ' Subcomposition   Minimum         Maximum')
1170  format (4x,a,3x,g12.5,4x,g12.5)

      end

      subroutine meelim (x,i,ii,j,k)
c----------------------------------------------------------------------
c subroutine to write unnatural limit warnings for meemum
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer ii, i, j, k

      double precision x

      character char8*8
c                                 global variables:
c                                 working arrays
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m14+2)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)
c----------------------------------------------------------------------
         if (poly(i).eq.1.and.istg(i,1).eq.1) then
c                                 single site solution or orphan
            if (ksmod(i).eq.20) then

               if (k.lt.ns) then
                  char8 = names(jnd(k))
               else
                  char8 = aqnam(jnd(k)  - aqst)
               end if

            else
               char8 = names(jend(i,2+k))
            end if

            write (*,1010) char8, x, fname(i), xmng(i,ii,j,k),
     *                     xmxg(i,ii,j,k)

         else

            write (*,1020) ii, j, k, x, fname(i), xmng(i,ii,j,k),
     *                     xmxg(i,ii,j,k)

         end if

         if (refine) then
            write (*,1000) 'the *.arf file and restart MEEMUM.'
         else
            write (*,1000) 'the solution model file and restart MEEMUM.'
         end if

1000  format ('then relax the limit in ',a,/)
1010  format (/,'**warning ver991** X(',a,') = ',f6.4,' of'
     *       ,' solution ',a,' exceeds its current',/,'limits (XMIN = ',
     *  f6.4,', XMAX = ',f6.4,') if this restriction is unintentional,')
1020  format (/,'**warning ver991** X(',i1,i1,i1,') = ',f6.4,' of ',
     *       'solution ',a,' exceeds its',/,'current limits (XMIN = ',
     *  f6.4,', XMAX = ',f6.4,') if this restriction is unintentional,')

      end

      subroutine err993 (ids,ii,j,k,pos)
c----------------------------------------------------------------------
c subroutine to write ver993 error diagnostic.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical pos

      integer ii, ids, j, k

      character char8*8, incre*8, upper*5
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------

      if (pos) then 
         upper = 'upper'
         incre = 'increase'
         y(1) = xmxg(ids,ii,j,k) + xncg(ids,ii,j,k)/2d0
         if (y(1).gt.1d0) y(1) = 1d0
      else 
         upper = 'lower'
         incre = 'decrease'
         y(1) = xmng(ids,ii,j,k) - xncg(ids,ii,j,k)/2d0
         if (y(1).lt.0d0) y(1) = 0d0
      end if

      if (istg(ids,1).eq.1.and.pop1(ids).eq.1) then
c                                 single site solution
         if (ksmod(ids).eq.20) then
c                                 charge balance model:
            if (j.lt.ns) then
               char8 = names(jnd(k))
            else
               char8 = aqnam(jnd(k)  - aqst)
            end if

         else

            char8 = names(jend(ids,2+k))

         end if

         write (*,1020) fname(ids), char8
         write (*,1030) xmnh(ids,ii,j,k), xmxh(ids,ii,j,k), x(ii,j,k)

      else if (pop1(ids).eq.1) then 
c                                 single polytope
         if (ksmod(ids).ne.688) then
            write (*,1040) fname(ids),j,k
            write (*,1030) xmnh(ids,ii,j,k), xmxh(ids,ii,j,k), x(ii,j,k)
            write (*,1000)
         else
            write (*,1080) fname(ids), poname(ids,ii,j,k)
            write (*,1035) poname(ids,ii,j,k), xmnh(ids,ii,j,k), 
     *                     xmxh(ids,ii,j,k), x(ii,j,k)
         end if

      else if (ii.lt.pop1(ids)) then 
c                                 composite polytope
         write (*,1050) fname(ids), poname(ids,ii,j,k), 
     *                              poname(ids,pop1(ids),1,ii)
         write (*,1035) poname(ids,ii,j,k), xmnh(ids,ii,j,k), 
     *                  xmxh(ids,ii,j,k), x(ii,j,k)

      else

         write (*,1060) fname(ids), poname(ids,pop1(ids),1,1)
         write (*,1035) poname(ids,pop1(ids),1,1),xmnh(ids,ii,j,k),
     *                  xmxh(ids,ii,j,k), x(ii,j,k)

      end if

      write (*,1070) 'www.perplex.ethz.ch/perplex/faq/warning_'//
     *               'ver991_relax_solution_model_limits.txt'

1000  format (/,'*NOTE: if this solution model has been reformulated '
     *       ,'because of missing endmembers',/,'the variable indices ',
     *        'may not correspond to the indices in the solution model',
     *        ' file.')
1020  format (/,'**warning ver993** the composition of solution: '
     *       ,a,/'is beyond the subdivision range limits for '
     *       ,'endmember: ',a)
1030  format ('the original range was: ',
     *       g12.5,' - ',g12.5,/,'the current** value is: ',g12.5)
1035  format ('the original range of ',a,' was: ',
     *       g12.5,' - ',g12.5,/,'its current** value is: ',g12.5)
1040  format (/,'**warning ver993** the composition of solution: '
     *       ,a,/'is beyond the subdivision range limits for '
     *       ,'composition X(',i1,',',i2,')*.')
1050  format (/,'**warning ver993** the composition of solution '
     *       ,a,' is beyond',/,'the subdivision range of'
     *       ,' composition variable ',a,' of the ',a
     *       ,' subcomposition.')
1060  format (/,'**warning ver993** the composition of solution: '
     *       ,a,/'is beyond the subdivision range limits for '
     *       ,'subcomposition: ',a)
1070  format (/,'refer to: ',//,a,//'for additional information.',/)
1080  format (/,'**warning ver993** the composition of solution '
     *       ,a,' is beyond',/,'the subdivision range of'
     *       ,' composition variable: ',a)
      end


      subroutine gaqlgd (gtot,blk,totm,smo,id,bad,recalc)
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation
c configured to be called from resub with output to the (molar normalized)
c arrays g2/cp2, if recalc is true then aqlagd is being used to recover
c speciation and arrays g2/cp2 are not loaded.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id

      logical bad, recalc, lmus, feos

      double precision mo(l9), blk(*), gamm0, totm, g0(l9), lmu(k8),
     *                 tmu(k8),is, gso(nsp), lnkw, gtot, smo, err,
     *                 slvmo(nsp), solmol, negox, posox

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      character cname*5
      common/ csta4  /cname(k5)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision fwt
      common/ cst338 /fwt(k10)

      double precision atwt
      common/ cst45 /atwt(k0)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      double precision gh,dvhy,g0mrk
      common/ csthyb /gh(nsp),dvhy(nsp),g0mrk(nsp)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      double precision p,t,xco2,mmu,tr,pr,r,ps
      common/ cst5 /p,t,xco2,mmu(2),tr,pr,r,ps

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      logical abort
      common/ cstabo /abort
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      save lmu, lmus
c----------------------------------------------------------------------
      if ((.not.mus .and..not.recalc).or.
     *    (.not.lmus.and.recalc)) then

         lmus = .false.
         bad = .true.
         return

      else
c                                 load dependent chemical potentials
         if (recalc) then
c                                 use lagged chemical potentials
            tmu(1:kbulk) = lmu(1:kbulk)
c                                 set flag so slvnt3 evaluates pure
c                                 fluid eos (not clear if this is necessary).
            feos = .true.

         else

            lmus = .true.

            do i = 1, kbulk

               lmu(i) = mu(i)
               tmu(i) = mu(i)

               if (cblk(i).eq.0d0.and..not.lopt(36).and.i.le.jbulk) then
c                                 check that the solvent does not contain
c                                 the absent component
                  do j = 1, ns

                     if (pa(j).gt.0d0.and.cp(i,jnd(j)).gt.0d0) then

                        bad = .true.
                        return

                     end if

                  end do

               end if

               feos = .false.

            end do

         end if

         call slvnt3 (gso,.false.,feos,id)

         if (epsln.lt.nopt(34).or.abort) then
c                                 eos is predicting vapor phase
c                                 solvent densities
            bad = .true.
            return

         end if

         bad = .false.

      end if
c                                 iterate on speciation
      call aqsolv (g0,gso,mo,tmu,is,gamm0,lnkw,bad)
c                                 in new_solver.f needs to be debugged?
c     call aqsol2 (g0,gso,mo,tmu,is,gamm0,lnkw,bad)

      if (bad) return
c                                 back calculated bulk composition
      blk(1:kbulk) = 0d0
      smo = 0d0
      gtot = 0d0
      err = 0d0
c                                 everything on a molal basis
c                                 first the solutes
      do i = 1, aqct

         if (mo(i).eq.0d0) cycle
c                                 charge balance error
         err = err + q(i) * mo(i)
c                                 total g
         gtot = gtot + mo(i) * (g0(i) + rt*dlog(mo(i)*gamm0**q2(i)))
c                                 total molality
         smo = smo + mo(i)
c                                 accumulate component moles
         do j = 1, kbulk
            blk(j) = blk(j) + mo(i)*aqcp(j,i)
         end do

      end do

      solmol = smo
c                                 for the solvent mole fractions
c                                 need to accumulate total
c                                 molality first
      do i = 1, ns
c                                 solvent molality:
         slvmo(i) = yf(ins(i))/msol
c                                 total molality
         smo = smo + slvmo(i)
c                                 moles/kg-solvent
         do j = 1, kbulk
            blk(j) = blk(j) + slvmo(i)*cp(j,jnd(i))
         end do

      end do

      do i = 1, ns
c                                 solvent bulk mole fraction:
         if (recalc) caq(id,i) = slvmo(i)/smo
         if (slvmo(i).le.0d0) cycle
         gtot = gtot + slvmo(i) * (gso(i) + rt*dlog(slvmo(i)/smo))

      end do
c                                 bulk fluid composition
      totm = 0d0

      err = dabs(err)*1d1

      if (lopt(36).and.oxchg) then
c                                check on charge imbalance
         posox = 0d0
         negox = 0d0

         do j = 1, kbulk
            if (cox(j).gt.0) then
               posox = posox + cox(j)*blk(j)
            else
               i = j
               negox = negox + cox(j)*blk(j)
            end if
         end do

         blk(i) = blk(i) - (posox+negox)/cox(i)

      end if

      do j = 1, kbulk
c                                zero bulk compositions below chg balance error
         if (blk(j).lt.err) blk(j) = 0d0
c                                totm is the total number of moles of themodynamic components
c                                components in a solution of smo moles of
c                                species
         if (j.gt.icp) cycle

         totm = totm + blk(j)

      end do

      if (recalc) then
c                                 stuff needed for output:
         do i = 1, aqct
            caq(id,ns+i) = mo(i)
         end do
c                                 ionic strength
         caq(id,na1) = is
c                                 total molality
         caq(id,na2) = smo
c                                 solvent mass
         caq(id,na3) = msol
c                                 error in log10(Kw)
         caq(id,na3+1) = (-lnkw + dlog(mo(ihy)*mo(ioh)*gamm0**2
     *                    ))/2.302585d0
c                                  pH
         caq(id,na3+2) = -dlog10(mo(ihy)*gamm0)
c                                  Delta_pH
         caq(id,na3+3) = caq(id,na3+2) + lnkw/4.605170d0
c                                  solute molality
         caq(id,na3+4) = solmol
c                                  net charge
         posox = 0d0

         if (oxchg) then
c                                check on charge imbalance
            do j = 1, kbulk
               posox = posox + cox(j)*blk(j)
            end do

         end if

         caq(id,na3+5) = posox/smo
c                                  dielectric cst
         caq(id,nat) = epsln

      else
c                                 stuff need for optimization:
c                                 legendre transform for saturated/mobile components
         do j = icp+1, kbulk
            gtot = gtot - blk(j) * mu(j)
         end do

      end if

      end


      subroutine slvnt3 (gso,whysp,feos,id)
c-----------------------------------------------------------------------
c for back and lagged speciation calculations get solvent properties.
c   if whysp -> take the fluid fractions from ysp(:,id) and renormalize.
c   if feos -> call fluid eos (i.e., pure g's haven't been calculated).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical whysp, feos

      integer i, k, id

      double precision gso(nsp), dum, gcpd, ysum

      external gcpd

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision yf,gmrk,v
      common/ cstcoh /yf(nsp),gmrk(nsp),v(nsp)

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision gh,dvhy,gmrk0
      common/ csthyb /gh(nsp),dvhy(nsp),gmrk0(nsp)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)
c----------------------------------------------------------------------
      rt  = r*t

      if (ns.gt.1) then
c                                 a multi species solvent is present:
         if (whysp) then

            ysum = 0d0

            do i = 1, ns
               pa(i) = ysp(i,id)
               ysum = ysum + pa(i)
            end do
c                                 renormalize
            pa(1:ns) = pa(1:ns)/ysum

         end if

         if (feos) then

            do k = 1, ns
c                                 call to gcpd will put the pure species molar
c                                 volumes arrays vhyb0 and vmrk0 (cxt38) and
c                                 the mrk pure fugacity coeff in gp.
               aqg(k) = gcpd(jnd(k),.false.)

            end do

         else

            do k = 1, ns
c                                 use previously computed g
               aqg(k) = g(jnd(k))
c                                 unproject for mobile and saturated components:
               do i = icp+1, kbulk
                  aqg(k) = aqg(k) + cp(i,jnd(k))*mu(i)
               end do

            end do

         end if
c                                  dum is just a dummy.
         call slvnt1 (dum)

         do k = 1, ns
c                                 add in the solvent species activity coefficients
c                                 under the assumption that these are independent of the
c                                 solute speciation the partial g of a solvent species will be
c                                 g(i) = gs0(i) + RT ln x(i).
            i = ins(k)

            gso(k) = aqg(k) + rt * dlog(gmrk(i)/gmrk0(i))

         end do

      else
c                                  solvent is pure water
         pa(1) = 1d0
         ysp(1,id) = 1d0

         call slvnt0 (gso(1),dum)

      end if

      end

      subroutine aqsolv (g0,gso,mo,mu,is,gamm0,lnkw,bad)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, it, jt, iexp, iwarn

      logical bad, kill, switch

      double precision c(l9), mo(*), mu(*), dg, xis, dn, xdix,
     *                 d(l9), is, gamm0, g0(*), lnkw, dix,
     *                 gso(*), xdn, qb(l9), dnmax

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
            else if (dabs(dg/rt).gt.nopt(57)) then
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
               c(i) = d(i)

            else
c                                  neutral species assumed to be ideal, molality is
               mo(i) = dg

            end if

         end do
c                                  initialize iteration loop
         lnkw = (gso(ns) - g0(ioh))/rt

         if (c(ion).eq.0d0) then
c                                 no hydrogen or no oxygen
            bad = .true.
            return

         end if

         mo(ion) = dexp(lnkw/2d0)
         gamm0 = 1d0
         is = 0d0

         xdn = 1d0
         iexp = 1
         it = 0
         jt = 0
         xdix = 1d99
         dnmax = 0.5d0
         bad = .false.
         switch = .true.
c                                  iteration loop for ionic strength
         do
c                                  solve charge balance for ion
            mo(ion) = solve(c,qr,mo(ion),jchg,ichg,bad)

            xis = is

            if (bad) then
               kill = .false.
               exit
            end if
c                                  back calculate charged species molalities
c                                  and ionic strength
            is = 0d0

            do i = 1, ichg
               j = jchg(i)
               mo(j) = c(j) / q(j) * mo(ion)**qr(j)
               is = is + q2(j) * mo(j)
            end do

            is = is / 2d0

            dn = is - xis

            if (dnmax.gt.zero) then 

               if (dn/xdn.lt.0d0.and.switch) then
c                                 switching directions
                  dnmax = dabs(dn/10d0)
                  switch = .false.

               else if (dn/xdn.gt.0d0.and..not.switch) then

                  dnmax = dabs(dn/10d0)
                  switch = .true.

               end if

            end if

            if (dabs(dn).gt.dnmax) then 

               is = xis + dabs(dn)/dn*dnmax

            end if

c           if (dabs(dn).gt.1d0/2d0**iexp) then
c              dn = dn/dabs(dn)/2d0**iexp
c              if (dn*xdn.lt.0d0) iexp = iexp + 1
c           end if
c                                 DH law activity coefficient factor (gamma = gamm0^(q^2))
            gamm0 = aqact(is)
            if (gamm0.lt.nopt(50)) gamm0 = nopt(50)
c                                 check for convergence
            dix = dabs(xis-is)/(1d0+is)

            if (dix.lt.nopt(50)) then
c                                 converged
c              call aqsol2 (g0,gso,mo,mu,is,gamm0,lnkw,bad)

               return

            else if (it.gt.iopt(21)) then

               if (xdix.gt.dix.and.jt.lt.10) then
c                                 try again?
                  it = 0
                  jt = jt + 1
                  xdix = dix

               else
c                                 diverging
                  kill = .true.
                  bad = kill

c                 call aqsol2 (g0,gso,mo,mu,is,gamm0,lnkw,bad)

                  exit

               end if

            end if

            xdn = dn

            it = it + 1
c                                 update coefficients
            do i = 1, ichg
               j = jchg(i)
               c(j) = d(j)*gamm0**qb(j)
            end do

         end do
c                                 switch to the backup ion
         if (lopt(44)) then
            ion = ioh
         else
            ion = ihy
         end if

      end do
c                                 failure is the only path here
      if (kill.and.iwarn.lt.11) then

         call warn (64,is,it,' ')

         call prtptx

         if (iwarn.eq.10) call warn (49,0d0,64,'AQSOLV')

         iwarn = iwarn + 1

      end if

      end


      double precision function aqact (is)
c-----------------------------------------------------------------------
c compute Debye-Hueckel-type activity coefficient factor
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision is

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh
c----------------------------------------------------------------------
      aqact = dexp(adh*dsqrt(is)/(1d0 + dsqrt(is)) + 0.2d0*is)

      end

      subroutine dumper (iclos,id,tkp,ukp,amt,lambda)
c----------------------------------------------------------------------
c dump phase data from yclos routines:
c     iclos = 1 - static, 2 - dynamic
c     hkp refinement point
c     lkp solution model
c     amt - amount
c     lambda - lambda
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, iclos, id, tkp, ukp

      double precision amt, lambda

      character name*14

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------
      call getnam (name,ukp)

      if (iclos.eq.1) then
c                                 static
         write (*,1000) id,tkp,ukp,name,amt,lambda,c(id),
     *                  (a(i,id),i=1,jbulk)
      else
c                                 dynamic
         write (*,1000) id,tkp,ukp,name,amt,lambda,g2(id),
     *                  (cp2(i,id),i=1,jbulk)
      end if

1000  format (i7,1x,i3,1x,i4,1x,a,20(g14.6,1x))

      end

      double precision function ginc (dt,dp,id)
c-----------------------------------------------------------------------
      implicit none

      double precision dt,dp,gee,gsol,gfrnd

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 put NaN check to prevent a NaN increment
c                                 from setting p/t to a NaN, this should only
c                                 happen if a phase has no volumetric and/or
c                                 caloric EoS. JADC 9/18/2016.
      if (isnan(dp)) dp = 0d0
      if (isnan(dt)) dt = 0d0

      p = p + dp
      t = t + dt

      if (iam.eq.5) then
c                                 frendly
         gee = gfrnd(-id)

      else
c                                 meemum/werami
         gee = gsol(id)

      end if

      p = p - dp
      t = t - dt

      ginc = gee

      end

      double precision function gsol (id)
c-----------------------------------------------------------------------
c gsol computes the total (excess+ideal) free energy of solution
c for a solution identified by index ids and composition pa(m4) input
c from cxt7

c gsol assumes the endmember g's have not been calculated by gall and is
c      called by WERAMI and MEEMUM via ginc for output purposes.

c gsol1 is identical to gsol but can only been called after gall and is
c      called by VERTEX and MEEMUM via minfrc.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, k, id, ntot

      double precision omega, gfluid, gzero, aqact, gmchpt, gmech0,
     *                 gex, gfesi, gcpd, gerk, gfecr1, ghybrid,
     *                 gfes, gfesic, g, gso(nsp), gamm0, gdqf

      external gphase, omega, gfluid, gzero, gex, gfesic, gdqf,
     *         gfesi, gerk, gfecr1, ghybrid, gcpd, aqact, gfes,
     *         gmchpt, gmech0

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer jnd
      double precision aqg,q2s,rt
      common/ cxt2 /aqg(m4),q2s(m4),rt,jnd(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision yf,gh,v
      common/ cstcoh /yf(nsp),gh(nsp),v(nsp)

      integer jd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,jd

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)
c----------------------------------------------------------------------
      if (id.lt.0) then

         gsol = gcpd (-id,.true.)

      else

         g = 0d0
c                                 initialize p-t dependent coefficients
         call ingsol (id)

         if (specil(id)) then

            g =  gfesic (pa(1),pa(3),pa(4),
     *           gcpd(jend(id,3),.true.),gcpd(jend(id,4),.true.),
     *           gcpd(jend(id,5),.true.),gcpd(jend(id,6),.true.),
     *           ksmod(id))

         else if (lorder(id)) then
c                                 -------------------------------------
c                                 as gsol may be called many times for the same 
c                                 bulk composition, the pa array will change, reset
c                                 to the p0a values
            ntot = nstot(id)

            pa(1:ntot) = p0a(1:ntot)

            if (.not.noder(id)) then 
c                                 get the speciation, excess and entropy effects.
               call specis (g,id)

            else 
c                                 derivative free speciation
               call minfxc (g,id,.false.)

            end if
c                                 for dnu ~=0 this g is normalized to the 
c                                 p0a mass.
            g = g + gmchpt(id) + gdqf(id)

         else if (lrecip(id).or.simple(id)) then
c                                 -------------------------------------
c                                 macroscopic reciprocal solution w/o order-disorder
            g = gmchpt(id) + gdqf(id) - t*omega(id,pa) + gex(id,pa)

         else if (ksmod(id).eq.20) then
c                                 electrolytic solution
c                                 -------------------------------------
c                                 compute solvent mass and gibbs energy:
            rt = r*t

            do k = 1, ns
c                                 solvent species gibbs energy and volumes
               if (pa(k).le.0d0) cycle
               aqg(k) = gcpd(jnd(k),.true.)

            end do
c                                 solvent properties
            call slvnt1 (g)
c                                 add in solute properties
            call slvnt2 (g)

         else if (ksmod(id).eq.26) then
c                                 ------------------------------------
c                                 andreas salt model
            call hcneos (g,pa(1),pa(2),pa(3))

            g = g + gmchpt (id)

         else if (ksmod(id).eq.29) then
c                                 -------------------------------------
c                                 BCC Fe-Si Lacaze and Sundman
            g = gfesi(pa(1), gcpd (jend(id,3),.true.),
     *                       gcpd (jend(id,4),.true.) )

         else if (ksmod(id).eq.32) then
c                                 -------------------------------------
c                                 BCC Fe-Cr Andersson and Sundman
            g =  gfecr1(pa(1), gcpd (jend(id,3),.true.),
     *                         gcpd (jend(id,4),.true.) )

         else if (ksmod(id).eq.39) then
c                                 -------------------------------------
c                                 generic hybrid EoS
            if (lopt(32).and.caq(jd,na1).ne.0d0) then
c                                 lagged speciation
               call slvnt3 (gso,.false.,.true.,id)
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
               gamm0 = aqact(caq(jd,na1))
c                                 solvent species (caq => mole fraction)
               do k = 1, ns
                  if (caq(jd,k).eq.0d0) cycle
                  g = g + caq(jd,k) * (gso(k) + rt * dlog(caq(jd,k)))
               end do
c                                 solute species (caq => molality)
               do k = sn1, nsa

                  i = k - ns

                  if (caq(jd,k).eq.0d0) cycle

                  g = g + caq(jd,k)/caq(jd,na2) * (gcpd(aqst+i,.false.)
     *                  + rt*dlog(caq(jd,k)*gamm0**q2(i)))

               end do

            else
c                                 mech + activities
               g = gmchpt (id) + ghybrid (pa)

            end if

         else if (ksmod(id).eq.41) then
c                                 hybrid MRK ternary COH fluid
            call rkcoh6 (pa(2),pa(1),g)

            g = g + gmchpt (id)

         else if (ksmod(id).eq.40) then
c                                 MRK silicate vapor
            g = gmech0 (id) + gerk (pa)

         else if (ksmod(id).eq.42) then
c                                 ------------------------------------
c                                 Fe-S fluid (Saxena & Eriksson 2015)
            g = gfes (pa(2), gcpd (jend(id,3),.true.),
     *                       gcpd (jend(id,4),.true.) )

         else if (ksmod(id).eq.0) then
c                                 ------------------------------------
c                                 internal fluid eos. hardwired to special
c                                 component choices.
c                                 don't know whether it's a speciation routine
c                                 so set the fluid species fractions just in case,
c                                 this is only necessay for species output by
c                                 WERAMI/MEEMUM, these will be reset if it actually
c                                 is a speciation routine.
            yf(2) = pa(1)
            yf(1) = 1d0 - yf(2)
c
            g = gmech0(id) + gfluid (yf(2))

         else

            write (*,*) 'what the **** am i doing here?'
            stop

         end if

         gsol = g

      end if

      end

      double precision function gfrnd (id)
c-----------------------------------------------------------------------
c function to get g's for frendly. differs from gphase in that it checks
c for special components O2, H2O, CO2. sloppy but who cares?
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision gee, fo2, fs2, gcpd

      external gcpd

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer eos
      common/ cst303 /eos(k10)

      integer ifct,idfl
      common/ cst208 /ifct,idfl
c-----------------------------------------------------------------------

      gee = gcpd (id,.false.) + r * t * dlog(act(id))

      if (ifct.gt.0.and.eos(id).lt.100) then
c                                 this is a quick fix that will
c                                 call the fluid routine way more
c                                 than necessary.
         call cfluid (fo2,fs2)

         if (id.eq.idf(3)) then

            gee = gee + r*t*fo2

         else if (id.eq.idf(1)) then

            gee = gee + r*t*fh2o

         else if (id.eq.idf(2)) then

            gee = gee + r*t*fco2

         end if

      end if

      gfrnd = gee

      end

      subroutine setstc (ids,h,i,j)
c----------------------------------------------------------------------
c set stretching transformation for the site i, species j, of solution ids,
c it seems doubtful this is worth the effort, functions strtch and unstch
c could just use the stch array directly.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids,h,i,j

      double precision stch
      common/ cst47 /stch(h9,h4,mst,msp,4)

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm
c----------------------------------------------------------------------
      bm1  = stch(ids,h,i,j,1)
      bp1  = stch(ids,h,i,j,2)
      bpm  = stch(ids,h,i,j,3)
      lbpm = stch(ids,h,i,j,4)

      end

      double precision function stinc (x,dy,ids,h,i,j)
c----------------------------------------------------------------------
c stinc increments stretched x by the conformal increment dy for solution
c ids, site i, species j. this is used to set limits for the conformal
c y in terms of the stretched x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids,h,i,j

      double precision x, y, dy, strtch, unstch

      external strtch, unstch
c----------------------------------------------------------------------
c                                 set stretching parameters
      call setstc (ids,h,i,j)
c                                 unstretch x and decrement/increment it
c                                 by +/- one conformal increament, then
c                                 restretch.
      y =  unstch (x) + dy

      if (y.gt.1d0) then
         y = 1d0
      else if (y.lt.0d0) then
         y = 0d0
      end if

      stinc = strtch ( y )

      end

      subroutine gpmlt1 (g,k,id,error)
c----------------------------------------------------------------------
c subroutine to speciation of the green et al (JMG, 2016) melt model. this
c model is a special case because the model has a single ordering parameter, which
c green et al take as the fraction of the ordered species (an). this formulation is
c unfortunate because p(an) is not orthogonal to the disordered speciation
c (p0, because the moles of the species is not constant with changing speciation).
c here the model is recast as g(p0,q) where q is the number of moles of an that can be
c formed given p0.

c    id - identifies the solution.
c    g  - change in G for the stable speciation relative to a mechanical
c          mixture of the endmembers.
c    pc is the mass normalization factor, sum(p0*ctot)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, k, id, itic

      logical error, done

      double precision g, qmax, qmin, q, q0, dqq, rqmax, rqmin, gold

      double precision omega, gex
      external omega, gex

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      double precision enth
      common/ cxt35 /enth(j3)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c----------------------------------------------------------------------
      error = .false.
c                                 rqmax the maximum amount of the
c                                 ordered species that can be formed
c                                 from the fully disordered species
c                                 fractions

c                                 this solver DOES NOT account for the
c                                 antiordered state! is there one? i donut
c                                 think so. There is, but it's nonsense, JADC 11/22
      rqmax = 1d0
      rqmin = 0d0 

      do i = 1, nrct(k,id)
c                                 this is probably ok for HP melt models
c                                 as the endmember fractions are generally
c                                 related to a site fraction
         if (dydy(ideps(i,k,id),k,id).lt.0d0) then

            if (-p0a(ideps(i,k,id))/dydy(ideps(i,k,id),k,id).lt.rqmax)
     *              rqmax = -p0a(ideps(i,k,id))/dydy(ideps(i,k,id),k,id)

         else 

            if (-p0a(ideps(i,k,id))/dydy(ideps(i,k,id),k,id).gt.rqmin)
     *              rqmin = -p0a(ideps(i,k,id))/dydy(ideps(i,k,id),k,id)

         end if

      end do

      q0 = p0a(nstot(id))
      rqmin = q0 + rqmin
      rqmax = q0 + rqmax
c                                 to avoid singularity set the initial
c                                 composition to the max - nopt(50), at this
c                                 condition the first derivative < 0,
c                                 and the second derivative > 0 (otherwise
c                                 the root must lie at p > pmax - nopt(50).
      if (rqmax.gt.nopt(50)) then

         pin(k) = .true.
         qmax = rqmax - nopt(50)
         qmin = rqmin + nopt(50)
c                                 the p's are computed in gpderi
         call gpder1 (k,id,qmax-q0,dqq,g,.false.)

         if (dqq.lt.0d0) then
c                                 at the maximum concentration, the
c                                 first derivative is positive, if
c                                 the second is also > 0 then we're
c                                 business
            q = qmax

         else
c                                 try the min
            call gpder1 (k,id,qmin-q0,dqq,g,.false.)

            if (dqq.gt.0d0) then
c                                 ok
               q = qmin

            else
c                                 no search from either limit possible
c                                 set error .true. to compare g at the
c                                 limits.
               error = .true.
               goto 90

            end if
         end if
c                                 increment and check p
         call pcheck (q,qmin,qmax,dqq,done)
c                                 iteration counter to escape
c                                 infinite loops
         itic = 0

         gold = g
c                                 newton raphson iteration
         do

            call gpder1 (k,id,q-q0,dqq,g,.false.)

            call pcheck (q,qmin,qmax,dqq,done)
c                                 done is just a flag to quit
            if (done.or.dabs((gold-g)/(1d0+dabs(g))).lt.nopt(50)) then

c              if (done.and.dabs((gold-g)/g).gt.nopt(53)) then 
c                 write (*,*) 'oink3',gold-g,g,itic,id
c              end if

               goodc(1) = goodc(1) + 1d0
               goodc(2) = goodc(2) + dfloat(itic)
c                                 in principle the p's could be incremented
c                                 here and g evaluated for the last update.
               return

            else

               gold = g

            end if

            itic = itic + 1

            if (itic.gt.iopt(21)) then
c                                 fails to converge.
               error = .true.
               badc(1) = badc(1) + 1d0
               goodc(2) = goodc(2) + dfloat(itic)
               exit

            end if

         end do

      else
c                                 speciation is stoichiometrically frustrated
         g = -t*omega(id,p0a) + gex(id,p0a)
         return

      end if

90    if (error) then
c                                 didn't converge or couldn't find a good
c                                 starting point compute the fully ordered
c                                 g, specis will compare this to the
c                                 disordered g and take the lowest:
         do i = 1, nstot(id)
            pa(i) = (p0a(i) + dydy(i,k,id)*rqmax)/(1d0 +dnu(k,id)*rqmax)
         end do

         g = (pa(nstot(id))*enth(k) - t*omega(id,pa) + gex(id,pa)) *
     *       (1d0 + dnu(k,id)*rqmax)

      end if

      end

      subroutine qlim (dqmin,dqmax,lord,id)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, k, id, jd, lord

      double precision dqmax(*), dqmin(*)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      lord = 0
c                                 rqmax the maximum amount of the
c                                 ordered species that can be formed
c                                 from the fully disordered species
c                                 fractions
      do k = 1, nord(id)

         dqmax(k) = 1d0

         do i = 1, nrct(k,id)
c                                 this is probably ok for HP melt models
c                                 as the endmember fractions are generally
c                                 related to a site fraction
            jd = ideps(i,k,id)

            if (dydy(jd,k,id).gt.0d0) cycle

            if (-p0a(jd)/dydy(jd,k,id).lt.dqmax(k))
     *                           dqmax(k) = -p0a(jd)/dydy(jd,k,id)

         end do

         dqmin(k) = -p0a(lstot(id)+k) + nopt(50)
         dqmax(k) = dqmax(k) - nopt(50)

         if (dqmax(k)-dqmin(k).gt.nopt(50)) then
            pin(k) = .true.
            lord = lord + 1
         else
            pin(k) = .false.
         end if

      end do

      end

      subroutine gpmelt (g,id,minfx)
c----------------------------------------------------------------------
c subroutine to non-equilimolar speciation order-disorder. this
c model is a special case because the model has a single ordering parameter, which
c green et al take as the fraction of the ordered species (an). this formulation is
c unfortunate because p(an) is not orthogonal to the disordered speciation
c (p0, because the moles of the species is not constant with changing speciation).
c here the model is recast as g(p0,q) where q is the number of moles of an that can be
c formed given p0.

c    id - identifies the solution.
c    g  - change in G for the stable speciation relative to a mechanical
c          mixture of the endmembers.
c    pc is the mass normalization factor, sum(p0*ctot)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id, itic, lord

      logical error, minfx

      double precision g, dqmax(j3), dqmin(j3), dqq(j3), ddq(j3), gold,
     *                 tdp, xtdp

      double precision omega, gex
      external omega, gex

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      double precision enth
      common/ cxt35 /enth(j3)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c----------------------------------------------------------------------
      error = .false.
      minfx = .false.
c                                 count free order parameters and get limits
      call qlim (dqmin,dqmax,lord,id)
c                                 if lord = 1, switch to 1d solver
      if (lord.eq.1) then
         do k = 1, nord(id)
            if (pin(k)) then
               call gpmlt1 (g,k,id,error)
               return
            end if
         end do
      end if

c                                 set initial q values
      if (refine) then
c                                 use known speciation
         do k = 1, nord(id)
            dqq(k) = 0d0
         end do 

      else
c                                 assume parameters are independent
c                                 set each to 0.9*qmax as in speci2
         do k = 1, nord(id)
            dqq(k) = 0.9d0 * (dqmax(k)-dqmin(k))
         end do

      end if

      if (lord.gt.0) then

         itic = 0
         gold = 0d0
         xtdp = 0d0
         minfx = .false.

         do

            call gpderi (id,dqq,g,ddq,.false.,error)

            if (error) then
               badc(1) = badc(1) + 1d0
               exit
            end if

            tdp = 0d0
c                                 increment q's
            do k = 1, nord(id)

               if (.not.pin(k)) cycle

               if (dqq(k)+ddq(k).gt.dqmax(k)) then
                  ddq(k) = dqmax(k) - dqq(k)
                  dqq(k) = dqmax(k)
               else if (dqq(k)+ddq(k).lt.dqmin(k)) then
                  ddq(k) = dqmin(k) - ddq(k)
                  dqq(k) = dqmin(k)
               else
                  dqq(k) = dqq(k) + ddq(k)
               end if

               tdp = tdp + dabs(ddq(k))

            end do
c                                 check for convergence
            if ((tdp.lt.nopt(50).or.
     *          dabs((gold-g)/(1d0+dabs(g))).lt.nopt(50))
     *         .and.itic.gt.1) then

c              if (tdp.lt.nopt(52).and.dabs((gold-g)/g).gt.nopt(53))
c    *            then 
c                 write (*,*) 'oink2',gold-g,g,itic,id
c              end if

               goodc(1) = goodc(1) + 1d0
               goodc(2) = goodc(2) + dfloat(itic)
               exit

            else if (itic.gt.25.and.gold.lt.g) then

               minfx = .true.
               exit

            else if (itic.gt.iopt(21)) then 

c              write (*,*) 'div2 ',gold-g,id,itic,g,tdp,tdp-xtdp
               minfx = .true. 
               exit

            else if (itic.gt.5.and.tdp.eq.xtdp) then 

               minfx = .true. 
c              write (*,*) 'wroink67 ',dp(1:lord),id,g
               exit

            end if

            itic = itic + 1

            xtdp = tdp

            gold = g

         end do

      end if

      if (error) then
c                                 speciation is stoichiometrically frustrated,
c                                 didn't converge, or couldn't find a good
c                                 starting point compute the fully ordered
c                                 g, specis will compare this to the
c                                 disordered g and take the lowest:
         dqq(1:nord(id)) = dqmax(1:nord(id))
         call gpderi (id,dqq,g,ddq,.false.,error)
         error = .true.

      else if (lord.eq.0) then

         g = 1d99

      end if 

      end

      subroutine gpderi (id,q,g,dg,minfxc,error)
c----------------------------------------------------------------------
c subroutine to compute the newton-raphson increment (dp) in the ordering
c parameter from the 1st and 2nd derivatives of the g of a solution with
c non-equimolar. 

c   minfxc - input, indicates calling program specis/minfxc 
c   id - is the index of the solution model.
c   q - input, the current order parameters
c   g  - output, the gibbs energy of the solution at q
c   dg - output: if minfxc - the derivatives dg/dq, else the newton-raphson increments
c   error - output, true if couldn't solve for the increments

c assumptions:

c   2nd order excess function.
c   temkin s evaluation assumes no disordered endmembers.
c   no independent endmember is involved in more than one speciation reaction.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, m, i1, i2, id, ipvt(j3)

      logical minfxc, error

      double precision g, dg(*), d2g(j3,j3), s, ds(j3), d2s(j3,j3), 
     *                 q(*), dp(m14,j3), d2p(m14,j3,j3), d2gn(j3,j3), 
     *                 theta, dtheta(j3), d2thet(j3,j3), nn, dnn, norm,
     *                 nt, dnt(j3), d2nt(j3,j3), dz(j3),
     *                 d2z(j3,j3), lnz1, zlnz, dzlnz(j3), d2zlnz(j3,j3),
     *                 z, n(m11), dn(m11,j3), d2n(m11,j3,j3),
     *                 t, dt(j3), d2t(j3,j3), dnorm(j3)
c                                 working arrays
      double precision zz, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)
c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision alpha,dt0
      common/ cyt0  /alpha(m4),dt0(j3)

      double precision enth
      common/ cxt35 /enth(j3)

      logical pin
      common/ cyt2 /pin(j3)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 initialize, assume pa is initialized 
c                                 to p0a
      norder = nord(id)

      g   = 0d0
      dg(1:norder) = 0d0
      d2g(1:norder,1:norder) = 0d0

      s = 0d0
      ds(1:norder) = 0d0
      d2s(1:norder,1:norder) = 0d0

      norm = 1d0
      dnorm(1:norder) = dnu(1:norder,id)
c                                 the total number of moles after disordering
      do k = 1, norder
         norm = norm + dnorm(k) * q(k)
      end do

      theta = 1d0/norm
c                                 derivatives of theta with respect to q
      do k = 1, norder
         dtheta(k) = -dnorm(k)*theta/norm
         do j = 1, norder 
            d2thet(k,j) = -2d0*dtheta(k)*dnorm(j)/norm
         end do
      end do

      do i = 1, nstot(id)
c                                 compute the unnormalized species fractions, it is assumed that
c                                 each disordered species is related to only one 
c                                 ordered species, use of a pointer would eliminate
c                                 this loop
         nn = p0a(i)
         dnn = 0d0

         do k = 1, norder

            dnn = dydy(i,k,id)
            nn = nn + dnn * q(k)
            dp(i,k) = nn*dtheta(k) + dnn*theta
            if (dnn.ne.0d0) exit

         end do

         pa(i) = nn*theta

         do k = 1, norder
            do j = 1, norder
c                              d2n is always zero
               d2p(i,k,j) =   dnn * (dtheta(k) + dtheta(j))
     *                      +  nn * d2thet(k,j)
            end do
         end do
      end do

      if (llaar(id)) then
c                                 h&p van laar, initialize
         t = 0d0
         dt(1:norder) = 0d0
         d2t(1:norder,1:norder) = 0d0
c                                 t-derivatives
         do i = 1, nstot(id)

            t = t + alpha(i)* pa(i)

            do k = 1, norder

               dt(k) = dt(k) + alpha(i)* dp(i,k)

               do j = 1, norder
                  d2t(k,j) = d2t(k,j) + alpha(i)* d2p(i,k,j)
               end do
            end do
         end do
c                                 excess terms
         do i = 1, jterm(id)

            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)

            g = g + w(i) * pa(i1) * pa(i2)

            do k = 1, norder

               dg(k) = dg(k) 
     *                 + w(i) * (pa(i1)*dp(i2,k) + pa(i2)*dp(i1,k))

               do j = 1, norder

                  d2g(k,j) = d2g(k,j) + w(i) * 
     *                   (pa(i1)*d2p(i2,k,j) + pa(i2)*d2p(i1,k,j) 
     *                                       + dp(i1,j)*dp(i2,k) 
     *                                       + dp(i2,j)*dp(i1,k))
               end do
            end do
         end do
c                                 laar size normalization
         g = g/t

         do k = 1, norder
c                                 use ds, d2s as temporary storage for normalized
c                                 derivatives
c                                 WAS THIS IN ERROR BEFORE?
            ds(k) = (dg(k) - g*dt(k))/t

            do j = 1, norder

              d2s(k,j) = (d2g(k,j) - g*d2t(k,j)
     *                   + ((2d0*g*dt(k) - dg(k))*dt(j) - dg(j)*dt(k))/t
     *                   ) / t
            end do
         end do
c                                 reassign dg, d2g
         dg(1:norder) = ds(1:norder)
         d2g(1:norder,1:norder) = d2s(1:norder,1:norder)

      else

         do i = 1, jterm(id)
c                                 excess g assuming regular terms
            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)

            g = g + w(i) * pa(i1) * pa(i2)

            do k = 1, norder

               dg(k) = dg(k) + w(i) *
     *                   (pa(i1)*dp(i2,k) + pa(i2)*dp(i1,k))

               do j = 1, norder

                  d2g(k,j) = d2g(k,j) + w(i) * ( pa(i1)*d2p(i2,k,j)
     *                                + pa(i2)*d2p(i1,k,j) 
     *                                + dp(i1,j)*dp(i2,k) 
     *                                + dp(i2,j)*dp(i1,k) )
               end do
           end do
         end do
      end if
c                                 get the configurational entropy derivatives
      do i = 1, msite(id)

         nt = 0d0
         dnt(1:norder) = 0d0
         d2nt(1:norder,1:norder) = 0d0
         zlnz = 0d0
         dzlnz(1:norder) = 0d0
         d2zlnz(1:norder,1:norder) = 0d0

         if (zmult(id,i).eq.0d0) then
c                                 temkin
            do j = 1, zsp(id,i)

               n(j) = dcoef(0,j,i,id)
               dn(j,1:norder) = 0d0
               d2n(j,1:norder,1:norder) = 0d0

               do k = 1, lterm(j,i,id)
c                                 n(j) is molar site population
                  n(j) = n(j) + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))

                  do l = 1, norder

                     dn(j,l) = dn(j,l) 
     *                         + dcoef(k,j,i,id) * dp(ksub(k,j,i,id),l)

                     do m = 1, norder
                        d2n(j,l,m) = d2n(j,l,m) + dcoef(k,j,i,id) 
     *                               * d2p(ksub(k,j,i,id),l,m)
                     end do
                  end do
               end do

               nt = nt + n(j)

               do l = 1, norder
                  dnt(l) = dnt(l) + dn(j,l)
                  do m = 1, norder
                     d2nt(l,m) = d2nt(l,m) + d2n(j,l,m)
                  end do
               end do
            end do

            if (nt.gt.nopt(50)) then
c                                 site has non-zero multiplicity
               do j = 1, zsp(id,i)

                  z = n(j) / nt
c                                 zlnz is accumulated z*ln(z), lnz1 is 1 + ln(z)
                  call ckdzlz (z,zlnz,lnz1)
c                                 accumulate 1st derivative
                  do k = 1, norder

                     dz(k) = (dn(j,k) - z*dnt(k)) / nt
                     dzlnz(k) = dzlnz(k) + dz(k) * lnz1

                     do l = 1, norder

                        d2z(k,l) = (d2n(j,k,l) 
     *                             + (2d0 * z * dnt(k) * dnt(l)
     *                                         - dn(j,k) * dnt(l) 
     *                                         - dn(j,l) * dnt(k)) / nt
     *                                         - z * d2nt(k,l) ) / nt
                     end do
                  end do
c                                 accumulate 2nd derivative
                  do k = 1, norder

                     do l = 1, norder

                        d2zlnz(k,l) = d2zlnz(k,l) 
     *                                + d2z(k,l) * lnz1
     *                                + dz(k) * dz(l) / z
                     end do
                  end do
               end do
c                                 entropy units
               s = s - nt * zlnz

               do k = 1, norder

                  ds(k) = ds(k) - nt * dzlnz(k) - zlnz * dnt(k)

                  do l = k, norder

                     d2s(k,l) = d2s(k,l)
     *                          - nt * d2zlnz(k,l) - dnt(l) * dzlnz(k)
     *                          - zlnz * d2nt(k,l) - dzlnz(l) * dnt(k)
                  end do
               end do
            end if

         else
c                                 non-temkin
c                                 here nt is zt, dnt is dz, d2nt is d2z
            do j = 1, zsp(id,i)

               z = dcoef(0,j,i,id)
               dz(1:norder) = 0d0
               d2z(1:norder,1:norder) = 0d0
c                                 for each term:
               do k = 1, lterm(j,i,id)

                  z = z + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
c                                 accumulate first and second derivatives
                  do l = 1, norder
                     dz(l) = dz(l) 
     *                       + dcoef(k,j,i,id) * dp(ksub(k,j,i,id),l)
                     do m = 1, norder
                        d2z(l,m) = d2z(l,m) 
     *                       + dcoef(k,j,i,id) * d2p(ksub(k,j,i,id),l,m)
                     end do
                  end do
               end do

               call ckdzlz (z,zlnz,lnz1)
c                                 why bother? this is gonna be 1
               nt = nt + z

               do l = 1, norder

                  dzlnz(l) = dzlnz(l) + dz(l) * lnz1
c                                 the nth species fraction is evaluated by
c                                 difference, so this sum(dz(l),l=1..n-1) is non zero
                  dnt(l) = dnt(l) + dz(l)

                  do m = 1, norder

                     d2nt(l,m) = d2nt(l,m) + d2z(l,m)
c                                  dlnz1 = dz(l)/z
                     d2zlnz(l,m) = d2zlnz(l,m) 
     *                             + d2z(l,m) * lnz1
     *                             + dz(l) * dz(m) / z
                  end do
               end do
            end do
c                                 add the nth species:
            z = 1d0 - nt

            call ckdzlz (z,zlnz,lnz1)

            s = s - zmult(id,i)*zlnz/r

            do k = 1, norder
c                                  the contribution of the nth species
c                                  to dzlnz(k) is -dnt(k) * lnz1
               ds(k) = ds(k) - zmult(id,i)/r*(dzlnz(k) - dnt(k) * lnz1)

               do l = k, norder
c                                  dlnz1 = dz(l)/z
                  d2s(k,l) = d2s(k,l) - zmult(id,i)/r *
     *                       (d2zlnz(k,l) - d2nt(k,l) * lnz1
     *                                    - dz(l)/z * dnt(k))
               end do
            end do

         end if

      end do
c                                 add in entropic term to g
      g = g - r*v(2)*s
c                                 and add in enthalpic terms
      do k = 1, norder

         g      = g      + enth(k) * pa(lstot(id)+k)
         dg(k)  = dg(k)  + enth(k) * dp(lstot(id)+k,k) - r*v(2)*ds(k)

         do l = k, norder
            d2g(k,l) = d2g(k,l) + enth(k) * d2p(lstot(id)+k,k,l) 
     *                          - r*v(2)*d2s(k,l)
         end do
      end do
c                                 normalize
      do k = 1, norder

         dg(k)  = g * dnorm(k) + norm * dg(k)

         do l = k, norder
            d2gn(k,l) =   dg(l) * dnorm(k)
     *                  + dnorm(l) * dg(k) + norm * d2g(k,l)
         end do
      end do

      g = g * norm

      if (minfxc) return
c                          compute the newton-raphson increments:
      do k = 1, norder

         if (pin(k)) then

            dg(k) = -dg(k)
c                          flesh out the hessian
            do l = 1, k-1
               d2g(k,l) = d2g(l,k)
            end do

         else

            dg(k) = 1d0
            d2g(k,k) = 1d0

            do l = 1, norder
               if (l.eq.k) cycle
               d2g(l,k) = 0d0
               d2g(k,l) = 0d0
            end do

         end if
      end do
c                                 get the newton-raphson increments:
c                                 this is a general factorization routine, should
c                                 exploit that d2g is symmetric.
      call factor (d2g,j3,norder,ipvt,error)
c                                 solve for the increments by back-substitution,
c                                 this routine is also not efficient and should
c                                 be re written.
      if (.not.error) call subst (d2g,j3,ipvt,norder,dg,error)
c                                 substitute replaces the values of dg with the
c                                 newton-raphson increments for the ordered species
c                                 compositions.
c     dg(1) = -dg(1)/d2gn(1,1)

      end

      subroutine gpder1 (l,id,q,dg,g,minfxc)
c----------------------------------------------------------------------
c subroutine to compute the newton-raphson increment (dg) in the ordering
c parameter from the 1st and 2nd derivatives of the g of a temkin model
c with one ordering parameter. id is the index of the solution model.

c temkin s evaluation assumes no disordered endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, i1, i2, id

      logical minfxc

      double precision g, dg, d2g, s, ds, d2s, q, pnorm, pnorm2, 
     *                 d2p(m11), dng, gnorm, dgnorm, nt, dnt, d2nt, dz,
     *                 d2z, lnz1, zlnz, dzlnz, d2zlnz, nu, dp(m11),
     *                 z, n(m11), dn(m11), d2n(m11), t, dt, d2t
c                                 working arrays
      double precision zz, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)
c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision alpha,dt0
      common/ cyt0  /alpha(m4),dt0(j3)

      double precision enth
      common/ cxt35 /enth(j3)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 initialize
      g   = 0d0
      dg  = 0d0
      d2g = 0d0

      s = 0d0
      ds = 0d0
      d2s = 0d0

      gnorm  = 1d0 + dnu(l,id) * q
      dgnorm = dnu(l,id)
      pnorm  = 1d0/gnorm
      pnorm2 = 2d0*pnorm
c                                 the difficulty in this model is the
c                                 non-equimolar speciation reaction, this
c                                 causes the number of moles of the components
c                                 in a mole of solution to change as a function
c                                 of the order parameter even if composition is
c                                 held constant.

c                                 to keep the number of moles of the components
c                                 in the solution constant the gibbs energy
c                                 is multiplied by gnorm = 1 + q*sum(nu(i)), where
c                                 the nu(i) are the stoichiometric coefficients of
c                                 the endmembers in the ordering reaction (it being
c                                 assumed that nu(jd) = 1 and p0(jd) = 0). this gives
c                                 the solutions g when it has the same amounts of the
c                                 components as in the disordered limit (p = p0). the
c                                 amounts of the species (p) for a partially or completely
c                                 disordered state are p(i) = (p0(i) + nu(i))*q/gnorm.
c                                 q is the molar amount of the ordered species formed
c                                 by the ordering reaction from the amounts of the
c                                 reactant species in the disordered limit.

c                                 for the green et al melt model sum(nu(i)) for the
c                                 reaction wo + als = an is -1, therefore
c                                 gnorm = (1 - q) and pnorm = 1/(gnorm)
      do i = 1, nstot(id)
c                                 calculate pa, dp(i)/dq, d2p(i)/dq.
         nu = dydy(i,l,id)
         pa(i) = (p0a(i) + nu*q) * pnorm
         dp(i) = (nu - pa(i)*dnu(l,id)) * pnorm
         d2p(i) = dp(i) * pnorm2

      end do

      if (llaar(id)) then

         t = 0d0
         dt = 0d0
         d2t = 0d0
c                                 h&p van laar
         do i = 1, nstot(id)
            t = t + alpha(i)* pa(i)
            dt = dt + alpha(i)* dp(i)
            d2t = d2t + alpha(i)* d2p(i)
         end do

         do i = 1, jterm(id)

            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)

            g = g + w(i) * pa(i1) * pa(i2)
            dg = dg + w(i) * (pa(i1)*dp(i2) + pa(i2)*dp(i1))
            d2g = d2g + w(i) * (pa(i1)*d2p(i2) + pa(i2)*d2p(i1) 
     *                                         + 2d0*dp(i1)*dp(i2) )
         end do

         g = g/t
         dg =  dg - g*dt
         d2g = (d2g - 2d0*dt/t*dg - g*d2t)/t

      else

         do i = 1, jterm(id)
c                                 excess g assuming regular terms
            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)

            g = g + w(i) * pa(i1) * pa(i2)
            dg = dg + w(i) * (pa(i1)*dp(i2) + pa(i2)*dp(i1))
            d2g = d2g + w(i) * (      d2p(i1) * pa(i2)
     *                           + 2d0*dp(i2) * dp(i1)
     *                           +    d2p(i2) * pa(i1) )

         end do

      end if
c                                 get the configurational entropy derivatives
      do i = 1, msite(id)

         nt = 0d0
         dnt = 0d0
         d2nt = 0d0
         zlnz = 0d0
         dzlnz = 0d0
         d2zlnz = 0d0

         if (zmult(id,i).eq.0d0) then
c                                 temkin
            do j = 1, zsp(id,i)

               n(j) = dcoef(0,j,i,id)
               dn(j) = 0d0
               d2n(j) = 0d0

               do k = 1, lterm(j,i,id)
c                                 n(j) is molar site population
                  n(j) = n(j) + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
                  dn(j) = dn(j) + dcoef(k,j,i,id) * dp(ksub(k,j,i,id))
                  d2n(j) = d2n(j) + dcoef(k,j,i,id) *d2p(ksub(k,j,i,id))

               end do

               nt = nt + n(j)
               dnt = dnt + dn(j)
               d2nt = d2nt + d2n(j)

            end do

            if (nt.gt.nopt(50)) then
c                                 site has non-zero multiplicity
               do j = 1, zsp(id,i)

                  z = n(j)/nt
                  dz = (dn(j) - z*dnt)/nt
                  d2z = (2d0*dnt*(z*dnt-dn(j)) + nt*d2n(j) - n(j)*d2nt)
     *                  /nt**2

                  call ckdzlz (z,zlnz,lnz1)

                  dzlnz = dzlnz + dz * lnz1
                  d2zlnz = d2zlnz + d2z * lnz1 + dz**2/z

               end do
c                                 entropy units
               s = s - nt * zlnz
               ds = ds - nt * dzlnz - zlnz * dnt
               d2s = d2s - d2nt * zlnz - 2d0*dnt*dzlnz - d2zlnz*nt

            end if

         else
c                                 non-temkin
c                                 here nt is zt, dnt is dz, d2nt is d2z
            do j = 1, zsp(id,i)

               z = dcoef(0,j,i,id)
               dz = 0d0
               d2z = 0d0
c                                 for each term:
               do k = 1, lterm(j,i,id)
                  z = z + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
                  dz = dz + dcoef(k,j,i,id) * dp(ksub(k,j,i,id))
                  d2z = d2z + dcoef(k,j,i,id) * d2p(ksub(k,j,i,id))
               end do

               call ckdzlz (z,zlnz,lnz1)

               dzlnz = dzlnz + dz * lnz1
               d2zlnz = d2zlnz + d2z * lnz1 + dz**2 / z

               nt = nt + z
               dnt = dnt + dz
               d2nt = d2nt + d2z

            end do
c                                 add the contibution from the last species:
            z = 1d0 - nt

            call ckdzlz (z,zlnz,lnz1)

            s = s - zmult(id,i)*zlnz/r
            ds = ds - zmult(id,i)*(dzlnz - dnt * lnz1)/r
            d2s = d2s - zmult(id,i)*
     *                     (d2zlnz - d2nt * lnz1 + dnt**2 / z)/r

         end if

      end do

      g   = g   + enth(l)*pa(nstot(id))  - r*v(2)*s
      dg  = dg  + enth(l)*dp(nstot(id))  - r*v(2)*ds
c                                 the normalized g derivative
      dng  = g * dgnorm + gnorm * dg
c                                 the normalized g:
      g   = g * gnorm

      if (minfxc) then
         dg = dng
         return
      end if
c                                 and second derivative
      d2g = gnorm * (d2g + enth(l)*d2p(nstot(id)) - r*v(2)*d2s)
     *       + 2d0 * dg * dgnorm
c                                 dg becomes the newton-raphson increment:
      dg = -dng/d2g

      end

      subroutine zmake (z,i,l,ids)
c----------------------------------------------------------------------
c subroutine to the site fraction of ksp+1 th species on site i of
c endmember l
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision y(m4),zt,z

      integer i,j,k,l,ids

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
         y(1:nstot(ids)) = 0d0

         y(l) = 1d0

         zt = 0d0
c                                 get site fractions
         do j = 1, zsp(ids,i)

            z = dcoef(0,j,i,ids)
c                                 for each term:
            do k = 1, lterm(j,i,ids)
               z = z + dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))
            end do

            zt = zt + z

         end do

         z = 1d0 - zt

      end


      subroutine rmoden
c---------------------------------------------------------------------
c rmoden - reads 688 solution models from LUN n9.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical eor

      character tag*3, key*22, val*3, values*80, strg80*80,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40

      integer nreact, i, j, k, l, m, ier, idim

      double precision coeffs(k7), rnums(m4), enth(3), sum

      integer ijk(mst),inds(k7),ict

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer nsub,nterm
      double precision acoef
      common/ cst107 /acoef(m10,m11,0:m0),
     *     nterm(m10,m11),nsub(m10,m11,m0)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      character mname*8
      common/ cst18a /mname(m4)

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------
c                                read number of sub-polytopes
      call readda (rnums,1,tname)
      poly(h0) = idint(rnums(1))
c                                composite compositional simplex
      isimp(poly(h0)+1) = 1
      ivert(poly(h0)+1,1) = poly(h0)

      if (poly(h0).gt.h4) call error (1,rnums(1),poly(h0),
     *    'h4 (maximum number of subcompositions for solution model: '
     *     //tname//')')
c                                read subdivision ranges for the polytopes
      if (poly(h0).gt.1) then 
         call redsub (poly(h0)+1,tname)
      else
         poname(h0,1,1,1) = tname
      end if 
c                                initialize total number of polyyope vertices
      istot = 0 
c                                read data for each polytope
      do i = 1, poly(h0)
c                                number of simplices
         call readda (rnums,1,tname)

         isimp(i) = idint(rnums(1))

         if (isimp(i).gt.mst) call error (1,rnums(1),isimp(i),
     *      'mst (maximum number of simplices in a prism for '//
     *      'solution model: '//tname//')')

c                                number of vertices on each simplex:
         call readda (rnums,isimp(i),tname)

         ipvert(i) = 1

         pvptr(i,1) = istot + 1

         do j = 1, isimp(i)
            ivert(i,j) = idint(rnums(j))
c                                number of vertices in the sub-polytope
            ipvert(i) = ipvert(i)*ivert(i,j)
         end do

         pvptr(i,2) = istot + ipvert(i)

         if (istot.gt.m4) call error (1,rnums(1),istot,
     *                    'm4 (maximum number of endmembers)')

c                                read vertex names into mname
         call readn (istot,ipvert(i),tname)
c                               read subdivision data for each sub-polytope
         call redsub (i,tname)
c                               create pointer from the endmember l to its
c                               polytope vertex
         do j = 2, isimp(i)
            ijk(j) = 1
         end do

         ijk(1) = 0

         do l = istot+1, istot + ipvert(i)

            do m = 1, isimp(i)

               if (ijk(m).lt.ivert(i,m)) then

                  ijk(m) = ijk(m) + 1
c                                increment only one index per endmember
                  do j = 1, isimp(i)
                     jmsol(l,j) = ijk(j)
                  end do

                  if (ijk(m).eq.ivert(i,m).and.
     *                                     l.lt.istot+ipvert(i)) then
c                                saturated counter, increment first
c                                unsaturated counter and reset all
c                                lower site counters to first index,
                     do j = 2, isimp(i)

                        if (ijk(j).lt.ivert(i,j)) then

                           ijk(j) = ijk(j) + 1

                           do k = 2, j-1
                              ijk(k) = 1
                           end do

                           ijk(1) = 0

                           exit 

                        end if

                     end do

                  end if

                  exit

               end if

            end do

         end do

         istot = istot + ipvert(i)

      end do
c                               look for and read optional data this may be, 
c                               sequentially:
c                                  1) begin_ordered_endmembers
c                                  2) begin_dependent_endmembers

      norder = 0
      mdep = 0
      idim = istot

      do

         call redcd1 (n9,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'begin_ordered_endmembe') then

            do
c                               on input nreact = -1 signals ordering reaction
               nreact = -1

               call readr (coeffs,enth,inds,idim,nreact,tname,eor)

               if (eor) then

                  exit

               else

                  ordmod = .true.
                  norder = norder + 1

                  if (istot+norder.gt.m4) call error (1,rnums(1),
     *                istot+norder,'m4 (maximum number of endmembers)')

                  if (norder.gt.j3) call error (5,rnums(1),norder,tname)

               end if

               do j = 1, 3
                  denth(norder,j) = enth(j)
               end do

               nr(norder) = nreact

               do j = 1, nreact
                  depnu(j,norder) = coeffs(j+1)
                  iddeps(j,norder) = inds(j+1)
               end do

            end do

         else if (key.eq.'begin_dependent_endmem') then

            do 
c                               nreact is returned by readr
               nreact = 0

               call readr (coeffs,enth,inds,idim,nreact,tname,eor)

               if (eor) then

                  exit

               else

                  depmod = .true.
                  recip = .true.
                  mdep = mdep + 1
                  if (mdep.gt.m15) call error (1,enth(1),mdep,'m15')

               end if

               jdep(mdep) = inds(1)
               ndph(mdep) = nreact - 1

               if (ndph(mdep).gt.j4) 
     *            call error (1,enth(1),ndph(mdep),'j4')

               sum = 0d0

               do j = 1, ndph(mdep)
                  sum = sum + coeffs(j+1)
                  nu(mdep,j) = coeffs(j+1)
                  idep(mdep,j) = inds(j+1)
               end do

              if (dabs(sum-1d0).gt.nopt(50)) then 
                  write (*,'(/,a,g12.6,/)') 'coefficient sum = ', sum
                  call error (72,sum,i,'dependent endmember '//
     *                 mname(inds(1))//' definition coefficients do not'
     *                 //' sum to 1, solution model: '//tname)
               end if

            end do

         else 
c                                 done, must be at excess function
            backspace (n9)
            exit

         end if

      end do
c                                 read excess function
      call readx (idim,tname)
c                                 expansion for S(configurational)
      call readda (rnums,1,tname)
c                                 total number of sites (including non-mixing)
      msite(h0) = idint(rnums(1))

      if (msite(h0).gt.m10) call error (1,0d0,i,'m10')
c                                 for each site
      do i = 1, msite(h0)
c                                 688: read site name
         call redcd0 (n9,ier,key,values,strg80)
         znames(h0,i,0) = key
c                                 # of species, effective and true site multiplicty.
         call readda (rnums,3,tname)
c                                 effective multiplicity
         zmult(h0,i) = rnums(2)
c                                 true multiplicity
         tzmult(h0,i) = rnums(3)
c                                 number of species 
         zsp1(h0,i) = idint(rnums(1))
c                                 number of independent site fractions
         if (zmult(h0,i).gt.0d0) then
            zsp(h0,i) = zsp1(h0,i) - 1
         else
            zsp(h0,i) = zsp1(h0,i)
         end if

         if (zsp1(h0,i).gt.m11) call error (1,0d0,zsp(h0,i),'m11')
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, zsp1(h0,i)
c                                 read expression for site
c                                 fraction of species j on
c                                 site i.
            call readz (coeffs,inds,ict,idim,tname,tag)

            znames(h0,i,j) = tag

            acoef(i,j,0) = coeffs(1)
            nterm(i,j) = ict - 1
            if (nterm(i,j).gt.m0) call error (33,0d0,m0,tname)
c                                 for each term:
            do k = 2, ict
c                                 all terms 1 order type, this
c                                 saved for compatability with
c                                 old versions:
               acoef(i,j,k-1)  = coeffs(k)
               nsub(i,j,k-1)   = inds(k)

            end do

         end do

      end do
c                                 read suffix used to complete structural formula
c                                 on output. 
      call redcd0 (n9,ier,key,values,strg80)
      zuffix(h0) = key
c                                 initialize endmember flags
      do i = 1, istot
         iend(i) = 0
      end do 
c                              look for van laar and/or dqf parameters
c                              endmember flags or the end_of_model tag
      call readop (idim,istot-mdep,tname)
c                                 save original indices, need this for
c                                 melt models etc that have species specific
c                                 equations of state.
      do i = 1, istot + norder
         iorig(i) = i
      end do

      end

      subroutine readef (idim,tname)
c----------------------------------------------------------------------
c readef - read solution model endmembers to be flagged so that they 
c are not identified by the solution model on output, assumes
c data on one line of less than 240 characters, the expected format is

c        name

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ier, match, idim, index

      character name*8, eod*3, tname*10

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot
c----------------------------------------------------------------------

      do

         call readcd (n9,ier,.true.)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') chars(1:3)

         if (eod.eq.'end') exit

         i = 1

         call readnm (i,index,com,ier,name)
         if (ier.ne.0) goto 90

         index = match (idim,ier,name)
         if (ier.ne.0) goto 90

         iend(index) = 1

      end do

      return

90    write (*,1000) tname,chars(1:com)
      write (*,1001)

      call errpau

1000  format ('**error ver200** READEF bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a,/)
1001  format (/,'usually this error is caused by a mispelled ',
     *          'endmember name.',/)

      end

      subroutine reforn (im,first)
c---------------------------------------------------------------------
c reforn - counts the number of species that can be respresented for a
c solution given the present endmembers.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, killed, nokill

      integer kill, ikill, jkill, kill1, i, j, kosp(mst,msp), kill2,
     *        k, l, im, ii, jpoly, jsimp, jvct, ksimp(mst)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot
c----------------------------------------------------------------------
      dedpol(1:poly(h0)) = .false.
c                                the increment from the polytope vertex
c                                to the endmember index
      do ii = 1, poly(h0)

         killed = .false.

         do

            do i = 1, isimp(ii)

               do j = 1, ivert(ii,i)
                  kosp(i,j) = 0
               end do

            end do

            nokill = .true.
c                                 for the sum(ivert) species compute
c                                 the difference between the number
c                                 of endmembers that do not have and
c                                 do have each species.

            do k = pvptr(ii,1), pvptr(ii,2)

               if (kdsol(k).ne.0) then 

                  do j = 1, isimp(ii)
                     kosp(j,jmsol(k,j)) = kosp(j,jmsol(k,j)) - 1
                  end do

               else 

                  nokill = .false.
                  killed = .true.

                  do j = 1, isimp(ii)
                     kosp(j,jmsol(k,j)) = kosp(j,jmsol(k,j)) + 1
                  end do

               end if 

            end do
c                                 no endmembers to kill on polytope ii
            if (nokill) exit
c                                 find the species that is missing
c                                 the most from the model
            kill  = -99
            kill1 = 0

            do i = 1, isimp(ii)

               do j = 1, ivert(ii,i)

                  if (kosp(i,j).gt.kill) then

                     ikill = i
                     jkill = j
                     kill = kosp(i,j)

                  else if (kosp(i,j).eq.kill) then
c                                 a tie:
c                                 of the two choices take the one that 
c                                 will kill the most dependent endmembers
                     kill2 = 0

                     do k = pvptr(ii,1), pvptr(ii,2)
c                                 count the number that will be killed
                        if (jmsol(k,ikill).eq.jkill.and.kdsol(k).eq.-2)
     *                     kill2 = kill2 + 1

                     end do

                     if (kill2.gt.kill1) then
c                                 this is more than before (kill1)
                        kill1 = kill2
                        ikill = i
                        jkill = j

                     else if (i.ne.ikill) then

                        do k = 1, isimp(ii)

                           ksimp(k) = 99

                           do l = 1, ivert(ii,k)
                              if (kosp(k,l).lt.ksimp(k)) 
     *                                         ksimp(k) = kosp(k,l)
                           end do
                        end do

                        if (ksimp(i).lt.ksimp(ikill)) then 
                           ikill = i
                           jkill = j
                        end if

                     end if
                  end if
               end do
            end do
c                                 kill the species jkill on site ikill
c                                 and reformulate the model (this is
c                                 inefficient, but who cares). kill02
c                                 does not clean the composition space,
c                                 this is done afterwards by repoly 
c                                 after the final set of endmembers 
c                                 has been identified.
            call kill02 (ii,ikill,jkill)

            if (ipvert(ii).eq.0) then

               dedpol(ii) = .true.

               if (first.and.isimp(ii).gt.1) call warn (100,0d0,101,
     *             'eliminated subcomposition '
     *             //poname(h0,poly(h0)+1,1,ii)/
     *             /'during reformulation of '//tname//
     *             ' due to missing endmembers.')
               exit

            end if

         end do

         if (ipvert(ii).gt.0.and.killed.and.first.and.isimp(ii).gt.1) 
     *      call warn (100,0d0,102,
     *          'reformulated subcomposition '
     *          //poname(h0,poly(h0)+1,1,ii)/
     *          /' of '//tname//' due to missing endmembers.')
c                                 next polytope
      end do
c                                 clean up the model by eliminating empty/
c                                 redundant polytopes and/or simplices:
      if (istot.lt.2) then
c                                 failed, rejected too many endmembers
         im = im - 1
         if (first) call warn (25,wg(1,1),jstot,tname)
         jstot = 0

      else 
c                                 check if the composition space includes
c                                 redundant polytopes and/or simplices:
         jvct  = 0 
         jpoly = 0
c                                 first eliminate dead polytopes
         do ii = 1, poly(h0)

            if (dedpol(ii)) cycle

            jpoly = jpoly + 1

            pvptr(jpoly,1) = pvptr(ii,1)
            pvptr(jpoly,2) = pvptr(ii,2)
c                                  shift the species indices
            do i = 1, ipvert(ii)

               do j = 1, isimp(ii)
                  jmsol(pvptr(jpoly,1)-1+i,j) = jmsol(pvptr(ii,1)-1+i,j)
               end do 

            end do
c                                 composition space vertex counter
            jvct = jvct + ipvert(ii)

            if (ii.lt.poly(h0)) then 
c                                 shift the subdivision ranges
c                                 for the composition space down
               pxmn(poly(h0)+1,1,jpoly) = pxmn(poly(h0)+1,1,ii)
               pxmx(poly(h0)+1,1,jpoly) = pxmx(poly(h0)+1,1,ii)
               pxnc(poly(h0)+1,1,jpoly) = pxnc(poly(h0)+1,1,ii)
               pimd(poly(h0)+1,1,jpoly) = pimd(poly(h0)+1,1,ii)
               poname(h0,poly(h0)+1,1,jpoly) = 
     *                            poname(h0,poly(h0)+1,1,ii)

            end if
c                                 shift all polytopes down
            ipvert(jpoly) = ipvert(ii)
            isimp(jpoly) = isimp(ii)

            do j = 1, isimp(ii)

               ivert(jpoly,j) = ivert(ii,j)

               do k = 1, ivert(ii,j) - 1
                  pxmn(jpoly,j,k) = pxmn(ii,j,k)
                  pxmx(jpoly,j,k) = pxmx(ii,j,k)
                  pxnc(jpoly,j,k) = pxnc(ii,j,k)
                  pimd(jpoly,j,k) = pimd(ii,j,k)
                  poname(h0,jpoly,j,k) = poname(h0,ii,j,k)
               end do

               poname(h0,jpoly,j,k) = poname(h0,ii,j,k)

            end do

         end do

         if (jpoly.eq.0) then 

            im = im - 1
            if (first) call warn (25,wg(1,1),jstot,tname)
            jstot = 0
            return 

         end if

         j = 0

         do ii = 1, poly(h0)

            if (dedpol(ii)) cycle

            j = j + 1
c                                shift composition space subdivision ranges left
            pxmn(jpoly+1,1,j) = pxmn(poly(h0)+1,1,ii)
            pxmx(jpoly+1,1,j) = pxmx(poly(h0)+1,1,ii)
            pxnc(jpoly+1,1,j) = pxnc(poly(h0)+1,1,ii)
            pimd(jpoly+1,1,j) = pimd(poly(h0)+1,1,ii)
            poname(h0,jpoly+1,1,j) = poname(h0,poly(h0)+1,1,ii)

         end do

         poly(h0) = j
c                                 ---------------------------------------------
c                                 eliminate redundant simplicies from polytopes
         recip = .false.

         do ii = 1, poly(h0)

            if (isimp(ii).gt.1) then

               jsimp = 0

               do j = 1, isimp(ii)

                  if (ivert(ii,j).eq.1) cycle

                  jsimp = jsimp + 1

                  do i = pvptr(ii,1), pvptr(ii,2)

                     jmsol(i,jsimp) = jmsol(i,j)

                  end do

                  ivert(ii,jsimp) = ivert(ii,j)

                  do k = 1, ivert(ii,j) - 1
                     pxmn(ii,jsimp,k) = pxmn(ii,j,k)
                     pxmx(ii,jsimp,k) = pxmx(ii,j,k)
                     pxnc(ii,jsimp,k) = pxnc(ii,j,k)
                     pimd(ii,jsimp,k) = pimd(ii,j,k)
                     poname(h0,ii,jsimp,k) = poname(h0,ii,j,k)
                  end do

                  poname(h0,ii,jsimp,k) = poname(h0,ii,j,k)

               end do

               isimp(ii) = jsimp

               if (jsimp.gt.1) recip = .true.

            end if

         end do

         istot = jvct

      end if
c                                 check if polytope model can be 
c                                 reduced to a simplex?
      end

      subroutine kill02 (pkill,ikill,jkill)
c---------------------------------------------------------------------
c killsp - eliminates species jkill from simplex ikill of polytope pkill
c in a solution model and reformulates the model accordingly
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical skip, bad

      integer jsp,jtic,morder,pkill,ii,ivct,
     *        i,j,ikill,jkill,kill,kdep,kdqf,ktic,jold,
     *        i2ni(m4),kwas(m4),k,l,itic,ijkill(m4),
     *        j2oj(msp),j2nj(msp),i2oi(m4)
c                                 dqf variables
      integer indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c                                 local input variables
      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer nsub,nterm
      double precision acoef
      common/ cst107 /acoef(m10,m11,0:m0),
     *                nterm(m10,m11),nsub(m10,m11,m0)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)
c----------------------------------------------------------------------

      do i = 1, isimp(pkill)
c                                 nothing happens
         if (i.ne.ikill) cycle
c                                 on a simplex where a vertex will be eliminated
         jsp = ivert(pkill,i) - 1
c                                 make old-to-new and new-to-old pointers for the
c                                 vertices
         jtic = 0

         do j = 1, ivert(pkill,i)

            if (j.ne.jkill) then
               jtic = jtic + 1
c                              pointer from new j to old j
               j2oj(jtic) = j
c                              pointer from old j to new j
               j2nj(j) = jtic
            end if

         end do
c                              reset vertex counter
         ivert(pkill,i) = jsp

         if (ivert(pkill,i).gt.1) then
c                              shift subdivision ranges
            do j = 1, ivert(pkill,i) - 1
               pxmn(pkill,i,j) = pxmn(pkill,i,j2oj(j))
               pxmx(pkill,i,j) = pxmx(pkill,i,j2oj(j))
               pxnc(pkill,i,j) = pxnc(pkill,i,j2oj(j))
               pimd(pkill,i,j) = pimd(pkill,i,j2oj(j))
               poname(h0,pkill,i,j) = poname(h0,pkill,i,j2oj(j))
            end do

            poname(h0,pkill,i,j) = poname(h0,pkill,i,j2oj(j))

         end if

         exit

      end do
c                                the endmembers to be eliminated are in the range
c                                pvptr(pkill,1):pvptr(pkill,2)
      kdep = 0

      if (depmod) then
c                                create an array which gives the
c                                original locations of the dependent
c                                endmembers, need this to be able to
c                                reorder the y2p array:
         do i = 1, istot
            if (kdsol(i).eq.-2) then
               kdep = kdep + 1
            end if
         end do

      end if

      do i = pvptr(pkill,1), pvptr(pkill,2)
c                                 kill endmembers with the species
c                                 to be deleted:
         if (jmsol(i,ikill).eq.jkill) kdsol(i) = -3
      end do
c                                 reset pvptr values
      ivct = 0

      do ii = 1, poly(h0)

         ipvert(ii) = 1

         do j = 1, isimp(ii)
c                                number of vertices in the sub-polytope
            ipvert(ii) = ipvert(ii)*ivert(ii,j)
         end do

         pvptr(ii,1) = ivct + 1
         ivct = ivct + ipvert(ii)
         pvptr(ii,2) = ivct

      end do
c                                 check the ordered species
      morder = 0
c                                 first check if the ordered endmember
c                                 may be stable
      do k = 1, norder
c                                 check if a missing constituent
         bad = .false.

         do j = 1, nr(k)
            if (kdsol(iddeps(j,k)).eq.-3) then
               bad = .true.
               exit
            end if
         end do

         if (bad) then
c                                 add species to the kill list
            kdsol(istot+k) = -3

         else

            morder = morder + 1
            kdsol(istot+k) = -1
            kwas(morder) = k

         end if

      end do

      itic = 0
      jtic = 0
      ktic = 0
      kill = 0

      do i = 1, istot + norder

         if (kdsol(i).ge.-2) then
c                                 replacement for istot (itic)
            itic = itic + 1
c                                 total vertex count
            if (i.le.istot) ktic = ktic + 1
c                                 pointers from new to old endmember index (i2oi)
            i2oi(itic) = i
c                                 pointers from new to old endmember index (i2ni)
            i2ni(i) = itic
c                                 pointer to original species index
            iorig(itic) = iorig(i)
c                                 number of missing endmembers (jtic)
            if (kdsol(i).eq.0) jtic = jtic + 1
c                                 reset the kdsol array
            kdsol(itic) = kdsol(i)

         else
c                                 kill records the killed endmembers
            kill = kill + 1
            ijkill(kill) = i

         end if

      end do
c                                 reset total and present counters
      istot = ktic

      jstot = ktic - jtic

      do i = 1, itic

         if (i.ge.pvptr(pkill,1).and.i.le.pvptr(pkill,2)) then
c                                 the endmember is on a polytope where something 
c                                 was killed:
            do j = 1, isimp(pkill)
               if (j.eq.ikill) then
                  jmsol(i,j) = j2nj(jmsol(i2oi(i),j))
               else
                  jmsol(i,j) = jmsol(i2oi(i),j)
               end if
            end do

         else
c                                 shift the jmsol indices for endmembers on 
c                                 polytopes where nothing happened
            do j = 1, mst
               jmsol(i,j) = jmsol(i2oi(i),j)
            end do

         end if

      end do
c                                --------------------------------------
c                                excess terms:
      itic = 0

      do i = 1, iterm
c                                check for forbidden terms (i.e., terms
c                                with a missing endmember
         skip = .false.
c                                 macroscopic formulation
         do j = 1, kill
c                                 check if subscript points to a killed
c                                 endmember
            do k = 1, rkord(i)
               if (isub(i,k).eq.ijkill(j)) then
                  skip = .true.
                  exit
               end if
            end do

            if (skip) exit

         end do

         if (skip) cycle
c                               the term is acceptable
         itic = itic + 1

         rkord(itic) = rkord(i)

         isub(itic,1:rkord(i)) = i2ni(isub(i,1:rkord(i)))

         if (xtyp.eq.0) then
c                                save the coefficient
            do j = 1, m3
               wg(itic,j) = wg(i,j)
            end do

         else
c                                 redlich kistler
            do j = 1, rkord(itic)
               do k = 1, m16
                  wk(k,j,itic) = wk(k,j,i)
               end do
            end do

         end if

      end do
c                                reset counter
      iterm = itic
c                                --------------------------------------
c                                van laar volume functions
      if (laar) then
         do i = 1, istot + morder
            do j = 1, m3
               vlaar(j,i) = vlaar(j,i2oi(i))
            end do
         end do
      end if
c                                 --------------------------------------
c                                 dqf corrections, this is sloppy since
c                                 uses istot instead of kstot
      if (idqf.gt.0) then

         kdqf = 0
c                                 check if a retained species has a dqf
c                                 correction
         do j = 1, idqf
c                                 the itoi index must be in the inner loop
c                                 in case the values of indq are not sequential
            do i = 1, istot
               if (indq(j).eq.i2oi(i)) then
c                                 found a dqf'd endmember
                  kdqf = kdqf + 1
                  indq(kdqf) = i
                  do k = 1, m3
                     dqf(k,kdqf) = dqf(k,j)
                  end do
                  exit
               end if
            end do

            if (kdqf.eq.idqf) exit

         end do

         idqf = kdqf

      end if
c                                 --------------------------------------
c                                 configurational entropy model

c                                 site fractions as a function of endmember fractions
      do i = 1, msite(h0)
c                                 for each species, read function to define
c                                 the site fraction of the species and eliminate
c                                 killed species

c                                 species counter is incremented in advance
c                                 and must be decremented before saving the
c                                 final value:
         jtic = 1

         do j = 1, zsp1(h0,i)

            ktic = 0
c                                 for each term:
            do k = 1, nterm(i,j)

               dead = .false.
               do l = 1, kill
                  if (nsub(i,j,k).eq.ijkill(l)) then
                     dead = .true.
                     exit
                  end if
               end do

               if (.not.dead) then
c                                 the term has survived (and therefore
c                                 also the species):
                  ktic = ktic + 1
c                                 but my dear peanut brained friend, do
c                                 not forget to move the pointer:
                  nsub(i,jtic,ktic) = i2ni(nsub(i,j,k))
                  acoef(i,jtic,ktic) = acoef(i,j,k)
               end if
            end do
c                                 ktic is the number of terms representing
c                                 the jth species.
            if (ktic.gt.0) then
c                                 increment the species counter
               znames(h0,i,jtic) = znames(h0,i,j)
               nterm(i,jtic) = ktic
               acoef(i,jtic,0) = acoef(i,j,0)
               jtic = jtic + 1

            end if

         end do

         zsp1(h0,i) = jtic - 1

         if (zmult(h0,i).gt.0d0) then
c                                 non-temkin site
            zsp(h0,i) = zsp1(h0,i) - 1
         else 
c                                 temkin
            zsp(h0,i) = zsp1(h0,i)
         end if

      end do
c                                 ---------------------------------------
c                                 ordered species:
      if (ordmod) then

         norder = morder

         if (morder.eq.0) then
c                                 there are no ordered species left
            ordmod = .false.

         else
c                                 shift the ordered species pointers
c                                 and data to eliminate kill ordered
c                                 species.
            do j = 1, morder

               jold = kwas(j)

               do i = 1, 3
                  denth(j,i) = denth(jold,i)
               end do

               nr(j) = nr(jold)

               do i = 1, nr(j)
                  iddeps(i,j) = i2ni(iddeps(i,jold))
                  depnu(i,j) = depnu(i,jold)
               end do

            end do

         end if

      end if
c                                 --------------------------------------
c                                 dependent endmember properties, the
      if (depmod) then
c                                 dependent endmembers have been reordered
c                                 in redep, but are still expressed in
c                                 terms of the old indices, so reset the
c                                 indices:
         do i = 1, mdep
            jdep(i) = i2ni(jdep(i))
            do j = 1, ndph(i)
               idep(i,j) = i2ni(idep(i,j))
            end do
         end do

      end if

      end

      subroutine subdiv (ids,gcind)
c---------------------------------------------------------------------
c static subdivision and data storage

c ids   - points to the solution/subdivision for the static case
c phct  - static pseudocompound counter (now same as iphct)
c gcind - global simplicial composition counter

c both ids and kds are necessary for dynamic, kds is not used for the
c static case.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical simpl

      integer i, j, ii, ids, ncomp, nind(h4), pos, nc, gcind,
     *        stind(h4), ipop1

      double precision twt

      integer ntot,npairs
      common/ cst86 /ntot,npairs
c---------------------------------------------------------------------
      if (ksmod(ids).eq.20) then
c                                 subdivision with charge balance
         call cartaq (ids)
c                                 assign to y()?
         return

      end if

      scoct = 0
c                                 if pop1(ids) = 1, then the solution has a simple 
c                                 polytopic composition space.

c                                 the starting position - 1 of the simplicial
c                                 coordinates in xco/zco
      stind(pop1(ids)) = icoct

      pwt(pop1(ids)) = 1d0

      if (pop1(ids).gt.1) then
c                                 composite polytopic composition space,
c                                 do the simplical composition space to
c                                 generates npair of coordinates used to 
c                                 weight each polytope. subpol loads the 
c                                 coordinates in xco/zco and the simplex
c                                 indices (unnecessarily) in the local 
c                                 sco array.
         call subpol (1d0,ids,pop1(ids))
c                                 the number of subdivisions will be:
         ipop1 = npairs

         simpl = .false.

      else

         ipop1 = 1

         simpl = .true.

      end if

      pos = stind(pop1(ids))

      do i = 1, ipop1

         if (simpl) then 

            scoct = 0

         else
c                                 reset the simplicial coordinate counter so as 
c                                 not to over-write the pop1 coordinates
            scoct = ipop1
c                                 get the polytope weights
            twt = 0d0

            nind(pop1(ids)) = i

            do j = 1, ndim(1,pop1(ids),ids)

               pos = pos + 1

               pwt(j) = xco(pos)
               twt = twt + pwt(j)

            end do

            pwt(j) = 1d0 - twt

         end if
c                                 initialize the total number of polytopic
c                                 compositions
         nc = 1
c                                 do the subdivisions for each polytope
         do ii = 1, poly(ids)

            if (pwt(ii).le.nopt(50)) then
               pwt(ii) = 0d0
               npol(ii) = 0
               cycle
            end if
c                                 the starting position of the simplicial
c                                 compositions for polytope ii
            stind(ii) = scoct

            call subpol (pwt(ii),ids,ii)
c                                  the number of simplicial compositions
c                                  generated for polytope ii
            npol(ii) = ntot
            nc = nc * ntot

         end do
c                                  generate all permutations of the polytopic
c                                  compositions at constant wt, initialization:
         ncomp = 1

         do ii = 1, poly(ids)
            nind(ii) = 1
         end do

         call setind (ids,stind,nind,gcind)
c                                  now generate all permutations of the polytopic 
c                                  compositions:
         do 

            ncomp = ncomp + 1

            if (ncomp.gt.nc) exit 
c                                  figure out the index to be incremented
            do j = 1, poly(ids)

               if (nind(j).lt.npol(j)) then

                  nind(j) = nind(j) + 1

                  exit

               else 

                  nind(j) = 1

               end if

            end do
c                                 save the indexes
            call setind (ids,stind,nind,gcind)

         end do 

      end do

      end

      subroutine subpol (wt,ids,ii)
c---------------------------------------------------------------------
c subpol does jth subdivision of polytope ii of solution ids.

c ids   - points to the solution/subdivision for, respectively, the 
c         static and dynamic cases
c wt    - the effective resolution will be res/wt.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, h, ids, nt, isite, nind(mst), ii

      double precision wt

      integer ntot,npairs
      common/ cst86 /ntot,npairs
c---------------------------------------------------------------------
      ntot = 1

      isite = istg(ids,ii)
c                                 subdivide each simplex of the polytope
      do i = 1, isite
c                                 starting position of the compositional coordinates
c                                 for simplex i
         spx(ii,i) = icoct
c                                 cartes loads the simplicial coordinates into
c                                 array simp
         call cartes (wt,i,ii,ids)
c                                 copy these into the static or dynamic array
         do h = 1, npairs*ndim(i,ii,ids)

            icoct = icoct + 1
c
            if (icoct.gt.k18) call err41 ('K18')
            xco(icoct) = simp(h)

         end do
c                                 the number of compositions in the simplex
         snp(i) = npairs
c                                 number of compositions in the polytope
         ntot = ntot * npairs

      end do

      nt = 1

      do i = 1, isite
c                                 initialize the indices
         nind(i) = 1
         scoct = scoct + 1
         if (scoct.gt.k13) call err41 ('K13')

         sco(scoct) = 1

      end do
c                                 generate all compositons in the polytope
      do

         nt = nt + 1

         if (nt.gt.ntot) exit 
c                                 figure out which index to increment
         do i = 1, isite

            if (nind(i).lt.snp(i)) then

               nind(i) = nind(i) + 1

               exit

            else 

               nind(i) = 1

            end if

         end do

         do i = 1, isite

            scoct = scoct + 1
            if (scoct.gt.k13) call err41 ('K13')

            sco(scoct) = nind(i)

         end do 

      end do

      end

      subroutine setexs (ids,id)
c-----------------------------------------------------------------------
c recover the static polytopic composition id of solution ids.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, id, ii, i, k, pos, ipop, jpos

      double precision sum

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
c                                 get the simplicial composition indices:
      ipop = pop1(ids)
c                                 static:
         if (ipop.gt.1) then 
c                                 composite composition space, load 
c                                 weights
            sum = 0d0
            pos = jcox(icox(id))

            do ii = 1, ndim(1,ipop,ids)

               pwt(ii) = xco(pos+ii)
               x(ipop,1,ii) = pwt(ii)
               sum = sum + pwt(ii)

            end do 

            jpos = icox(id) + 1

            if (sum.lt.1d0) then
               pwt(ii) = 1d0 - sum
            else 
               pwt(ii) = 0d0
            end if

            x(ipop,1,ii) = pwt(ii)

         else

            jpos = icox(id)
            pwt(1) = 1d0

         end if

         do ii = 1, poly(ids)
c                                  cycle on 0-wt polytopes for static composiions
c                                  because don't need absolute positions
            if (pwt(ii).eq.0d0) cycle
c                                 recover the polytope compositions
            do i = 1, istg(ids,ii)
c                                 skip 0-d simplices
               if (ndim(i,ii,ids).eq.0) then 
                  x(ii,1,1) = 1d0
                  cycle
               end if

               sum = 0d0
               pos = jcox(jpos)

               do k = 1, ndim(i,ii,ids)

                  x(ii,i,k) = xco(pos+k)
                  sum = sum + xco(pos+k)

               end do

               if (sum.lt.1d0) then
                  x(ii,i,k) = 1d0 - sum
               else 
                  x(ii,i,k) = 0d0
               end if

               jpos = jpos + 1

            end do

         end do

      end

      subroutine setind (ids,stind,nind,gcind)
c-----------------------------------------------------------------------
c after a call to subpol, setind loads the local simplicial indices into
c the static global index arrays and sets the local composition
c arrays. 
c  ii    - is the polytope index
c  stind(ii) - locates the starting position of the simplicial indices
c              for polytope ii in sco
c  nind(ii) - indicates the polytopic composition to be used to 
c             generate the bulk composition.
c  cind  - is the starting local composition index - 1
c  gcind - is the starting global composition index - 1 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer ii, i, ids, gcind, stind(h4), nind(h4), pos, ipop

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c-----------------------------------------------------------------------
      iphct = iphct + 1
      ipop = pop1(ids)
c                                 load simplicial composition indices
c                                 into static arrays:
         if (iphct.gt.k1) call err41 ('K1 [SOLOAD/SETIND]')

         icox(iphct) = gcind + 1

         if (ipop.gt.1) then 
c                                 composite space, save location of 
c                                 polytopic wts
           gcind = gcind + 1
           if (gcind.gt.k24) call err41 ('K24 [SOLOAD/SETIND]')
           jcox(gcind) = spx(ipop,1) + (nind(ipop)-1)*ndim(1,ipop,ids)

         end if
c                                 save location of each set of simplicial
c                                 coordinates in each polytope
         do ii = 1, poly(ids)

            if (pwt(ii).le.0d0) cycle

            pos = stind(ii) + (nind(ii)-1)*istg(ids,ii)

            do i = 1, istg(ids,ii)
c                                 skip 0-d simplices
               if (ndim(i,ii,ids).eq.0) cycle

               gcind = gcind + 1
               if (gcind.gt.k24) call err41 ('K24 [SOLOAD/SETIND]')
               jcox(gcind) = spx(ii,i) 
     *                       +  (sco(pos+i) - 1) * ndim(i,ii,ids)
            end do

         end do

         call setxyp (ids,iphct,bad)

         if (.not.bad) call soload (ids,bad)

         if (bad) then
            gcind = icox(iphct) - 1
            iphct = iphct - 1
         end if

      end

      subroutine setxyp (ids,id,bad)
c-----------------------------------------------------------------------
c load compositional coordinates from the static xco or dynamic zco
c arrays into simple compositional arrays for the of compound 
c id of solution ids.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad, zbad

      external zbad

      integer ids, id, tmp

      double precision zsite(m10,m11)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
      bad = .false.

      if (refine.and.iam.ne.15) then
c                                 auto-refine in vertex
         tmp = itxp(id-ipoint)
         pa(1:nstot(ids)) = txco(tmp + 1:tmp + nstot(ids))
         call makepp (ids)

         return

      end if
c                                 get the polytopic compositions:
      call setexs (ids,id)
c                                 convert to 1-d polytopic compositions, the bad
c                                 test is unnecessary for static compositions once
c                                 they have been loaded by soload. this could be
c                                 eliminated to save time.
      call xtoy (ids,bad)
c                                 xtoy returns bad if the composition is of a 
c                                 optionally non-refineable endmember, otherwise
c                                 xtoy sets the y's for the composite polytopic
c                                 composition.
      if (bad) return

      if (bdx(ids)) then
c                                 as ridiculous as this may seem, if ~sck, then this
c                                 is a relict equipartition model, BUT because people 
c                                 prefer the result that they get with site-checking do
c                                 a site check here and reject compositions with negative 
c                                 site fractions:
        bad = zbad(pa,ids,zsite,fname(ids),.false.,fname(ids))
        if (bad) return

      end if 
c                                 convert the y's into p0a/pp/pa arrays indexed
c                                 only by independent endmembers.
      call y2p0 (ids)

      end

      subroutine xtoy (ids,bad)
c----------------------------------------------------------------------
c subroutine to convert composite polytopic solution compositions (x)
c to geometric endmember fractions (y) for solution model ids. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, ii, k, l, m

      logical bad

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd
c----------------------------------------------------------------------

      bad  = .false.

         do ii = 1, poly(ids)

            k = 0

            if (pwt(ii).lt.nopt(50)) then 

               do l = pvert(ids,ii,1), pvert(ids,ii,2)

                  y(l) = 0d0

               end do

               cycle

            end if 

            do l = pvert(ids,ii,1), pvert(ids,ii,2)

               y(l) = 1d0

               do m = 1, istg(ids,ii)
                  y(l) = y(l)*x(ii,m,kmsol(ids,l,m))
               end do

               if (y(l).gt.nopt(56)) then
                  k = l
                  exit
               end if

            end do

            if (k.ne.0) then
c                                 reject pure independent endmember compositions.
               if (ldsol(k,ids).gt.0.and.nrf(ids)
     *                              .and.pwt(ii).gt.nopt(56)) then

                  bad = .true.

                  return

               end if

               y(k) = 1d0

               do l = pvert(ids,ii,1), pvert(ids,ii,2)

                  if (l.eq.k) cycle

                  y(l) = 0d0

               end do

            end if

            do l = pvert(ids,ii,1), pvert(ids,ii,2)

               y(l) = y(l) * pwt(ii)

            end do

         end do

      end


      subroutine reset (phct,gcind)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer phct, gcind

      gcind = icoz(phct) - 1
      phct = phct - 1

      end 


      subroutine getxcp (xc,ntot,ids)
c-----------------------------------------------------------------------
c getxcp gets the normalized chemical composition of solution ids from 
c its endmember fractions. called by minfxc for o/d models. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, ids

      double precision xc(*), ntot

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      xc(1:icomp) = 0d0
      ntot = 0d0

      do i = 1, nstot(ids)

         do j = 1, icomp 
            xc(j) = xc(j) + pa(i) * endc(ids,i,j)/endt(ids,i)
         end do 

         ntot = ntot + pa(i) * endt(ids,i)

      end do

      end


      subroutine getscp (scp,scptot,ids,jd)
c-----------------------------------------------------------------------
c getscp gets the bulk chemical composition of solution ids from the composition
c of its endmembers. the composition of the solution in terms of its endmembers
c must be set by a prior call to setxyp.

c jd is a pointer that is used only for lagged speciation. After a call to
c savrpc in meemum/vertex it points to the the composition in the cp2 array.
c However, prior to a call to savrpc it getscp cannot be used!
c For werami it points to the composition in the local caq array.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, jd, ids

      double precision scp(*), scptot, xx

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer iam
      common/ cst4 /iam

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd
c-----------------------------------------------------------------------

      scp(1:icomp) = 0d0

      if (lopt(32).and.ksmod(ids).eq.39) then

         if (rkwak) then
c                                  called by soload or after gaqlgd failure
            do i = 1, ns
               do j = 1, icomp 
                  scp(j) = scp(j) + pa(i) * cp(j,jnd(i))
               end do 
            end do

         else if (iam.eq.1.or.iam.eq.2) then 
c                                  meemum, vertex => during dynamic optimization
            do j = 1, icomp
               scp(j) = cp2(j,jd)*c2tot(jd)
            end do

         else
c                                  werami
            if (rkwak.or.caq(jd,na1).eq.0d0) then
c                                  pure solvent
               do i = 1, ns
                  do j = 1, icomp 
                     scp(j) = scp(j) + pa(i) * cp(j,jnd(i))
                  end do 
               end do

            else 
c                                  impure solvent
               do i = 1, ns
                  do j = 1, icomp 
                     scp(j) = scp(j) + caq(jd,i) * cp(j,jnd(i))
                  end do 
               end do

               do i = sn1, nsa

                  k = i - ns
c                                 convert molality to mole fraction (xx)
                  xx = caq(jd,i)/caq(jd,na2)

                  do j = 1, icomp
                     scp(j) = scp(j) + xx * aqcp(j,k)
                  end do  

               end do

           end if

         end if

      else if (ksmod(ids).eq.20) then 
c                                 electrolyte:
c                                 solute species  
         do i = sn1, nqs
            do j = 1, icomp
               scp(j) = scp(j) + pa(i) * aqcp(j,jnd(i) - aqst)
            end do
         end do 
c                                 solvent species 
         do i = 1, ns 
            do j = 1, icomp
               scp(j) = scp(j) + pa(i) * cp(j,jnd(i))
            end do
         end do

      else
c                                 normal solutions:
c                                 pp coordinates used to compute the composition
         do i = 1, lstot(ids)
            do j = 1, icomp
               scp(j) = scp(j) + pp(i) * endc(ids,i,j)
            end do
         end do

      end if

      scptot = 0d0
c                                 note normalization is to the total amount of
c                                 thermodynamic components.
      do i = 1, icp

         if (dabs(scp(i)).lt.nopt(50)) scp(i) = 0d0
         scptot = scptot + scp(i)

      end do

      end

      subroutine input2 (first)
c----------------------------------------------------------------------
c input2 reads the thermodynamic data file for most perplex programs, 
c a (the?) notable exception being frendly that calls the parallel 
c routine jnput2.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character mnames(k16*k17)*8

      double precision twt(k5),tsel(k5),tcox(k5),cst
 
      integer i, j, k, l, im, ict, ifer,inames, jphct, imak(k16), iox

      logical eof, good, first, tpro(k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      character cname*5
      common/ csta4 /cname(k5)

      character name*8
      common/ csta6 /name

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ic
      common/ cst42 /ic(k0)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision atwt
      common/ cst45 /atwt(k0) 

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer make
      common / cst335 /make(k10)

      integer eos
      common/ cst303 /eos(k10)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer iam
      common/ cst4 /iam

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec
c-----------------------------------------------------------------------
c                               initialization for each data set
c                               for k10 endmembers
      do i = 1, k10
         make(i) = 0 
         names(i) = ' '
      end do
c                               for k1 phases:
      do i = 1, k1
         ikp(i) = 0
      end do 
c                               other counters and flags:
      do i = 1, h5
         isct(i) = 0
      end do 
c                               counters for bounds
      iphct = 0
      lamin = 0 
      idsin = 0 
      idfl = 0
c                               flag for GFSM endmembers in the saturated component space
      ltemp1 = .false.
c                               read data base header, do component
c                               transformations, read make definitions.
      call topn2 (0)
c                               general input data for main program

c                               reorder thermodynamic components
c                               if the saturated phase components are 
c                               present
      if (lopt(7)) then

         do k = 1, ispec 
                             
            do i = 1, icp

               if (cname(i).eq.cmpnt(idspe(k))) then 

                  if (i.eq.k) exit 

                  cname(i) = cname(k)

                  do j = 1, 3
                     cst = dblk(j,i)
                     dblk(j,i) = dblk(j,k) 
                     dblk(j,k) = cst
                  end do 

                  cname(k) = cmpnt(idspe(k))

                  exit            

               end if 

            end do 

         end do 

      end if  
c                              load the old cbulk array
      if (ifct.gt.0) iphct = 2
c                               identify nonzero components.
c                               initialize icout(i) = 0
      do i = 1, icmpn
         icout(i) = 0
      end do

      do i = 1, icomp

         im = 0

         do j = 1, icmpn

            if (cname(i).eq.cmpnt(j)) then 

               twt(i) = atwt(j)
               tsel(i) = sel(j)
               tcox(i) = cox(j)
               tpro(i) = dispro(j)

               ic(i) = j
               icout(j) = 1

               do k = 1, ispec
                  if (j.eq.idspe(k)) then 
                     iff(k) = i
                     idfl = idfl + 1
                  end if 
               end do 
 
               im = 1

            end if 

         end do 
c                               write error message if a component
c                               was not found:
         if (im.eq.0) then 
            write (*,1230) cname(i), (cmpnt(k), k = 1, icmpn)
            write (*,1240)
            stop
         end if 
 
      end do 
c                                 this segment is to check if
c                                 a possible saturated phase component
c                                 has been made a mobile component,
c                                 if there is also a saturated phase
c                                 component idfl is the identity of the
c                                 mobile component otherwise idfl = 0.
      if (ifct.eq.1.and.idfl.eq.2) then

         do i = 1, ispec
            if (zname.ne.cmpnt(idspe(i))) cycle 
            idfl = i
            exit 
         end do 

      else 
         idfl = 0
      end if
c                                 load atwts, sel in updated order
      do i = 1, icomp
         atwt(i) = twt(i)
         sel(i)  = tsel(i)
         cox(i)  = tcox(i)
         dispro(i) = tpro(i)
         if (cox(i).lt.0d0) iox = i 
      end do 
c                                 convert weight to molar amounts
      if (jbulk.ne.0) then 

         if (iwt.eq.1) then 
            do i = 1, jbulk
               do j = 1, 3
                  dblk(j,i) = dblk(j,i)/atwt(i)
               end do 
            end do 
         end if

         call iniblk

      end if 
c                                 get composition vectors for entities
c                                 defined by a make definition:
      call makecp (inames,mnames,first)
c                                 loop to read reference phase data for
c                                 activity/fugacity variables
      ict = 0 

      if (ifact.gt.0) then
c                                 rewind and read 'til end of header
         call eohead (n2)

         good = .false.

         do

            call getphi (name,.false.,eof)

            if (eof) then 

               write (*,1000) (afname(i),i=1,jmct)
               write (*,1010)
               call errpau

            end if 
c                                 now look for a match with the 
c                                 reference phase names
            do i = 1, jmct

               if (name.eq.afname(i)) then 
c                                 got a match, count
                  iphct = iphct + 1

                  ict = ict + 1

                  idaf(i) = iphct
c                                 store thermodynamic parameters:
                  call loadit (iphct,.false.,.true.)
c                                 zero the component
c                 vnumu(i,iphct) = 0d0

                  if (imaf(i).eq.2) then 
c                                 if some cretin chooses fugacity, prevent
c                                 gphase from calling the EoS.   
                     eos(iphct) = ieos 

                  else if (lopt(7)) then 
c                                 check for special component names
c                                 this is necessary because loadit 
c                                 will not set isfp if ifct > 0.
                     do k = 1, ispec
                        if (name.ne.cmpnt(idspe(k))) cycle
                        eos(iphct) = 100 + k 
                        exit 
                     end do 
 
                  end if 
c                                 blank the name, this has two purposes,
c                                 it prevents problems if an entry is 
c                                 replicated in the data file, and flags
c                                 tagged entries 
                  afname(i) = ' '

                  if (ict.eq.ifact) good = .true.

                  exit 

               end if 

            end do 

            if (good) exit 

         end do 

      end if 
c                                 begin first read loop for data on
c                                 saturated components.
      if (isat.eq.0.and.ifct.eq.0) goto 40
c                                 read 'til end of header
      call eohead (n2)
c                                 loop to read real saturated
c                                 entities:
      ifer = 0

      do 

         call getphi (name,.false.,eof)

         if (eof) exit
 
         call chkphi (0,name,good)

         if (good) call sattst (ifer,.false.,good)

      end do 
c                                 loop to load made saturated entities
      do i = 1, nmak

         if (.not.mksat(i)) cycle
c                                 load make data 
         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                 redundant check:
         call chkphi (2,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)
c                                 set eos flag
         ieos = meos(i)

         call sattst (ifer,.true.,good)

         if (.not.good) call error (57,comp(1),iphct,name)

         if (good) then 
            make(iphct) = i
c                                 pointer used for iemod.
            imak(i) = iphct
         end if 

      end do 
c                                 check that there is data for
c                                 every fluid component.
      if (ifct.gt.0.and.ifer.ne.ifct) call error (36,r,i,'INPUT2')
c                                 check that there is one phase
c                                 for each saturation constraint
40    do i = 1, isat
         if (isct(i).lt.1) call error (15,r,i,cname(icp+i))
      end do

      if (isat.gt.0.and.first.and.(iam.lt.4.or.iam.eq.15)) then

         write (*,'(/,80(''-'')/,a)') 
     *         'Summary of saturated-component entities:'

         do i = 1, isat

            write (*,*) ' '
            write (*,1040) (cname(icp+j),j=1, i)

            do k = 1, isct(i), 6
               l = k + 5
               if (l.gt.isct(i)) l = isct(i)
               write (*,1050) (names(ids(i,j)), j = k, l)
            end do 
         end do

         if (iam.eq.15.and.isoct.gt.0) write (*,'(/,a)')
     *         '*solutions may also have compositions'
     *       //' consisting entirely of saturated components'

         write (*,'(80(''-''))')

      end if 
c                                 save endmembers that consist entirely 
c                                 of saturated phase or mobile components:
      kphct = iphct 

      if (ifct+jmct.gt.0) then 

         call eohead (n2)

         do 

            call getphi (name,.false.,eof)

            if (eof) exit

            call chkphi (4,name,good)

            if (.not.good) cycle 
c                                 reject phases already in the list
            do i = 1, kphct
               if (names(i).eq.name) then
                  good = .false.
                  exit
               end if 
            end do 

            if (.not.good) cycle             
c                                 matched a name
            iphct = iphct + 1
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)

         end do

      end if 
c                                 -------------------------------------
c                                 real entities in the thermodynamic 
c                                 composition space:
      istct = iphct + 1
c                                 increment between iphct and jphct counters
      jiinc = istct - 1
c                                 read till end of header
      call eohead (n2)
c                                 loop to load normal thermodynamic data:
      do  
    
         call getphi (name,.false.,eof)

         if (eof) exit 
c                                 check if valid phase:
         call chkphi (1,name,good)

         if (good) then 
c                                 acceptable data, count the phase:
            iphct = iphct + 1
c                                 for normalized composition:
            ctot(iphct) = tot
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)

         end if

      end do 
c                                 -------------------------------------
c                                 made entities (as opposed to the required
c                                 data read later):
      do i = 1, nmak

         if (mksat(i)) cycle
c                                 load make data 
         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                 redundant check, but makes ctot.
         call chkphi (3,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)

         iphct = iphct + 1
         ctot(iphct) = tot
c                                 set ieos flag to that of the first
c                                 real entity in the make definition
         ieos = meos(i)

         call loadit (iphct,.true.,.true.)

         make(iphct) = i
c                                 pointer used for iemod.
         imak(i) = iphct

      end do 
c                                 load thermodynamic data for make definitions and
c                                 solute species, at this point iphct points to the 
c                                 last real entity, save this value and restore it later.
      jphct = iphct
c                                 -------------------------------------
c                                 make definition data: this
c                                 data is saved in the arrays thermo
c                                 and cp by loadit, but are not counted,
c                                 i.e., the counters ipoint and iphct
c                                 are reset. soload will then load the
c                                 cp array over the values loaded here,
c                                 but thermo should not be affected. gmake
c                                 then gets the data using the array 
c                                 mkind. the names array will also be 
c                                 overwritten.
      call eohead (n2)

      do 

         call getphi (name,.true.,eof)

         if (eof) exit

         do i = 1, inames

            if (name.ne.mnames(i)) cycle
c                                 matched a name
            iphct = iphct + 1
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.false.)

         end do

      end do 

      do i = 1, nmak
c                                remake pointer array for makes 
         do j = 1, mknum(i)
            do k = jphct + 1, iphct
               if (names(k).ne.mknam(i,j)) cycle
               mkind(i,j) = k
            end do
         end do 
      end do  
c                                 -------------------------------------
c                                 aqueous species, thermo data, as is the
c                                 case for make data is loaded in thermo;
c                                 names and composition loaded into 
c                                 aqnam and aqcp.
      aqst = iphct 
c
      call eohead (n2)
c                                 loop to load solute data:
      do  
    
         call getphi (name,.true.,eof)

         if (eof) exit
c                                 skip non-solute standard state data
         if (ieos.ne.15.and.ieos.ne.16) cycle
c                                 check if valid species:
         call chkphi (5,name,good)
c                                 check for oxidation state of aqueous
c                                 data if aq_oxides is set:
c         if (good.and.lopt(36).and.oxchg) then

c            qchg = thermo(6,k10)

c            if (qchg.eq.0d0.and.comp(ic(iox)).ne.0d0.or.
c     *          qchg-cox(iox)*comp(ic(iox)).ne.0d0) then 

c               call warn (100,r,102,
c     *              name//' has been rejected; to retain '//name//
c     *              ' set aq_oxide_components to false.')

c               good = .false.

c            end if

c         end if 

         if (good) then 
c                                 acceptable data, count the phase:
            iphct = iphct + 1
c                                 for normalized composition, probably
c                                 con't need this, but could be used to
c                                 save molar wt or something like that:
            ctot(iphct) = tot
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)

         end if 

      end do
c                                write summary and checks
      if (aqct.gt.0) then 

         if (lopt(25).and.(ihy.eq.0.or.ioh.eq.0)) then 
            call warn (99,0d0,0,'missing H+ or OH- species, '//
     *                          'aq_output set = F (INPUT2)')
         end if 

         ichg = 0
         
         do i = 1, aqct 

            q(i) = thermo(6, aqst + i)
            q2(i) = q(i)**2

            if (q(i).ne.0d0) then 
               ichg = ichg + 1
               jchg(ichg) = i 
            end if

         end do 

         if (first.and.iam.lt.3) then 
            write (*,1020)
            do i = 1, aqct, 6
               k = i + 5
               if (k.gt.aqct) k = aqct
               write (*,1030) (aqnam(j),int(thermo(6,j+aqst)),j=i,k)
            end do 
            write (*,'(//)')
         end if

      else if (lopt(32).or.lopt(25)) then 

         if (first.and.iam.lt.4) 
     *       call warn (99,0d0,0,' no data for aqueous species, '
     *                 //'aq_output and aq_lagged_speciation disabled.')

         lopt(32) = .false.
         lopt(25) = .false.
  
      end if 
c                                reset ipoint counter, but do not 
c                                reset iphct, because the compositions
c                                of the make phases are necessary for
c                                chemical potential variables.
c                                really? then why was it reset here?
      iphct = jphct
      ipoint = jphct

      do i = 1, nmak
c                                make an iemod flag for made
c                                endmembers:
         pmod(h9) = .true.
         smod(h9) = .true.

         do j = 1, mknum(i)

            if (iemod(mkind(i,j)).eq.0) then 
               pmod(h9) = .false.
               smod(h9) = .false.
               exit
            else if (iemod(mkind(i,j)).eq.1) then
               pmod(h9) = .false.
            else if (iemod(mkind(i,j)).eq.3) then
               smod(h9) = .false.
            end if 

         end do 

         if (pmod(h9).and.smod(h9)) then 
            iemod(imak(i)) = 2
         else if (pmod(h9)) then 
            iemod(imak(i)) = 3
         else if (smod(h9)) then 
            iemod(imak(i)) = 1
         else 
            iemod(imak(i)) = 0
         end if

      end do
c                                 if saturated phase components, 
c                                 then check to make sure no GFSM
c                                 endmembers are in use as saturated 
c                                 components or in the thermodynamic
c                                 composition space:
      if (ifct.gt.0) then

         do i = ifct+1, iphct

            if (eos(i).gt.100.and.eos(i).lt.200) then 
c                                 got one
               write (*,1060) names(i), names(i)

               call wrnstp

            end if

         end do

      end if

      if (ltemp1) then

         do i = istct, iphct

            if (eos(i).gt.100.and.eos(i).lt.200) then 
c                                 got one
               write (*,1070) names(i), names(i)

               call wrnstp

            end if

         end do

      end if

1000  format ('**error ver007** at least one of the reference ',
     *        'endmembers:',/,5(a,1x))
1010  format ('needed to define an independent fugacity/activity ',
     *    'variable is missing',/,'most likely the endmember has ',
     *    'been rejected, if so then set',/,'the auto_exclude ',
     *    'option to FALSE.',/)
1020  format (/,'Summary of aqueous solute species:',//,
     *        6('name     chg   ')) 
1030  format (6(a,2x,i2,3x))
1040  format (2x,'for: ',5(a,1x))
1050  format (4x,6(a,2x))
1060  format (/,'**warning ver533** ',a,' is a molecular fluid species '
     *       ,'the presence of which is ',/,'inconsistent with satura',
     *        'ted phase component constraints if the saturated phase',
     *      /,'is a fluid. Possible courses of action are:',//,4x,
     *        '1) exclude ',a,' and restart.',/,4x,
     *        '2) remove the phase saturation constraint and restart.',/
     *    ,4x,'3) ignore this warning and continue execution.',//,
     *        'Continue (Y/N)?')
1070  format (/,'**warning ver534** ',a,' is a molecular fluid species '
     *       ,'the presence of which is ',/,'inconsistent with existe',
     *        'nce of molecular fluid species in the saturated ',/,
     *        'component composition space.',//,
     *        'Possible courses of action are:',//,4x,
     *        '1) exclude ',a,' and restart.',/,4x,
     *        '2) remove the component saturation constraint and ',
     *        'restart.',/,4x,
     *        '3) ignore this warning and continue execution.',//)
1230  format ('**error ver013** ',a,' is an incorrect component'
     *       ,' name, valid names are:',/,12(1x,a))
1240  format ('check for upper/lower case matches or extra blanks',/)

      close (n2)

      end

      subroutine setvr0 (i,j)
c--------------------------------------------------------------------
c setvr1 computes nodal variables for node ij, three cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont
c----------------------------------------------------------------------

      if (icont.eq.1) then 

         v(iv1) = vmin(iv1) + (i-1)*dv(iv1)
         v(iv2) = vmin(iv2) + (j-1)*dv(iv2)
         call incdp0

      else if (icont.eq.2) then 

         v(iv1) = vmin(iv1) + (j-1)*dv(iv1)
         call incdep (iv1)

         cx(1) =  (i-1)*dvr(1)
         call setblk 

      else 

         cx(1) = (i-1) * dvr(1)
         cx(2) = (j-1) * dvr(2)
         call setblk

      end if 

      end

      subroutine setvar 
c--------------------------------------------------------------------
c setvar initializes the variables for gridded minimization, three
c cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision rloopy,rloopx

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc 

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      if (iam.eq.3) then 
c                                 WERAMI (3), PSSECT (7):
c                                 jinc will only be ~= 1 only for 
c                                 2d intermediate grid results
         rloopy = dfloat((loopy-1)/jinc)
         rloopx = dfloat((loopx-1)/jinc)

      else

         rloopy = dfloat(loopy-1)
         rloopx = dfloat(loopx-1)

      end if 
c                                 for 1d calculations
      if (loopx.eq.1.or.loopx.eq.0) rloopx = rloopy

      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do

      call incdp0

      if (icopt.eq.7.and.fileio) then 
c                                using nodal coordinate system
         dvr(1) = 1d0

      else if (icopt.eq.9.or.icopt.eq.11) then 
c                                using non-thermodynamic coordinate frame
         dvr(1) = (vmx(1) - vmn(1))/rloopx
         dvr(2) = (vmx(2) - vmn(2))/rloopy

      else if (icopt.eq.12) then 

         dvr(1) = nopt(36)
         dvr(2) = 1
         loopx = iopt(36)
         rloopx = dfloat(loopx)

      else if (icont.eq.1) then 
c                                v(iv1) on x, v(iv2) on y
         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopx
         dvr(1) = dv(iv1)

         dv(iv2) = (vmax(iv2) - vmin(iv2))/rloopy
         dvr(2) = dv(iv2)

      else if (icont.eq.2) then 
c                               composition is on x, v(iv1) on y
         dvr(1) = 1d0/rloopx
         cx(1) = 0d0

         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopy
         dvr(2) = dv(iv1)

      else 
c                                compositions on both axes
         dvr(1) = 1d0/rloopx
         dvr(2) = 1d0/rloopy 
         cx(1) = 0d0
         cx(2) = 0d0

      end if 
c                                set the bulk composition:
      call iniblk

      end 

      subroutine iniblk
c---------------------------------------------------------------------
c iniblk initializes the bulk composition (1st comp read by input2)
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont
c----------------------------------------------------------------------

      do i = 1, jbulk
         cblk(i) = dblk(1,i)
      end do

      end

      subroutine inipot 
c--------------------------------------------------------------------
c inipot initializes the independent potential variables to their 
c minimum values
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------
c                                 initialize potentials
      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do 
c                                 set dependent potential, if it exists
      call incdp0

      end

      subroutine bplinp (err)
c-----------------------------------------------------------------------
c read the b-plot file that contains the information on the assemblages
c stable at each grid node
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical err

      integer jxco, kxco, i, j, ids, ier
c                                 -------------------------------------
c                                 global variables
c                                 global assemblage data
      integer icog,jcog
      common/ cxt17 /icog(k2),jcog(k2)

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      double precision bg
      common/ cxt19 /bg(k5,k2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision amu
      common/ cst48 /amu(k8,k2)

      integer iam
      common/ cst4 /iam

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd
c----------------------------------------------------------------------
c                                 assemblage counter
      ibulk = 0
c                                 pointer to solution compositional coordinates
      jxco = 0 
      kxco = 0

      err = .false. 

      do 

         ibulk = ibulk + 1

         if (ibulk.gt.k2) call error (183,0d0,k2,'BLINP')

         read (n5,*,end=99) icog(ibulk),jcog(ibulk),iap(ibulk)

         ias = iap(ibulk)
c                                if ias = 0, probably reading 
c                                an inconsistent blk file in unsplt
         if (ias.le.0) then 
            ier = 1
            exit 
         end if 
c                                phase molar amounts
         read (n5,*,iostat=ier) (bg(i,ibulk),i=1,iavar(3,ias))
         if (ier.ne.0) goto 99

         icox(ibulk) = jxco

         do i = 1, iavar(1,ias)

            ids = idasls(i,ias)     

            kxco = jxco + nstot(ids) 
            jxco = jxco + 1

            if (kxco.gt.k18) call error (61,0d0,k18,'BPLINP')

            read (n5,*,iostat=ier) (xco(j), j = jxco, kxco)

            if (ier.ne.0) goto 99

            if (lopt(32).and.ksmod(ids).eq.39) then 
c                                lagged speciation

               jxco = kxco + 1
               kxco = kxco + nat

               if (kxco.gt.k18) call error (61,0d0,k18,'BPLINP')

               read (n5,*,iostat=ier) (xco(j), j = jxco, kxco)
               if (ier.ne.0) goto 99

            end if  
         
            jxco = kxco

         end do 

         jxco = kxco  
c                                 read mu's
         read (n5,*,iostat=ier) (amu(i,ibulk), i = 1, kbulk)

         if (ier.ne.0) then 
c                                 if error on read most probably its
c                                 because of NaN's for the chemical 
c                                 potentials
            amu(1:kbulk,ibulk) = nopt(7)

            ier = 0

         end if

      end do

99    ibulk = ibulk - 1

      if (ier.ne.0) err = .true.

      end


      subroutine plinp (err)
c---------------------------------------------------------------------- 
c plinp - subroutine to read assemblage info for gridded min calculations.
c if icopt = 7 and fileio also reads nodal coordinates.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, jst, irep, kd, jend, ier, iend

      logical kount, err

      character text*(lchar)

      integer igrd
      common/ cst311/igrd(l7,l7)

      integer iam
      common/ cst4 /iam

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vip
      common/ cst28 /vip(l2,k2)

      character*100 cfname
      common/ cst227 /cfname

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      integer idstab,nstab,istab
      common/ cst34 /idstab(i11),nstab(i11),istab

      integer idsol,nrep,nph
      common/ cst38/idsol(k5,k3),nrep(k5,k3),nph(k3)

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)
c----------------------------------------------------------------------
      err = .false.

      if (iam.eq.7.and.plopt(3)) then
c                                 open assemblage list for PSSECT
         call mertxt (tfname,prject,'_assemblages.txt',0)
         open (n8, file = tfname, status = 'unknown', iostat = ier)

         write (*,'(a,a)') 'Assemblage list will be written to file: ',
     *                     tfname

         if (ier.ne.0) then 
            write (*,*) 'error cannot open: ',tfname
            write (*,*) 'file is probably open in an editor'
            call errpau
         end if

      end if
c                                 top of plot file
      read (n4,*,iostat=ier) loopx, loopy, jinc
c                                 check if the file was generated by unsplt
c                                 if so set unsplt flag for sample_on_grid
      if (jinc.eq.-1) then

         jinc = 1
         jlev = grid(3,2)
         lopt(47) = .true.

      else
         lopt(47) = .false.

      end if

      if (ier.ne.0) goto 99
c                                 prior to 6.8.5 vertex did not write 
c                                 the final value of jinc to the plot 
c                                 file, reset it here for back-compatibility
      if (loopx.eq.1.or.loopy.eq.1) jinc = 1
c                                 decompress the grid data
      do i = 1, loopx, jinc
         jst = 1
         do while (jst.le.loopy)
            read (n4,*,iostat=ier) irep, kd
            if (ier.ne.0) goto 99
            if (kd.eq.0) write (*,*) 'bad un at i, j',i,j
            jend = jst + irep 
            do j = jst, jend
               if (j.gt.l7) call error (2,nopt(1),j,
     *                      'coordinates (routine PLINP), increase L7')
               igrd(i,j) = kd
            end do 
            jst = jend + 1
         end do 
      end do 
c                                 read assemblages
      read (n4,*,iostat=ier) iasct
      if (ier.ne.0) goto 99
c                                 global stable phase counter
      istab = 0 
c                                 min/max number of phases in an assemblage
      piopt(1) = 100
      piopt(2) = 0

      do i = 1, iasct

         read (n4,*,iostat=ier) iavar(1,i),iavar(2,i),iavar(3,i)
         if (ier.ne.0) goto 99

         if (iavar(3,i).lt.piopt(1)) piopt(1) = iavar(3,i)
         if (iavar(3,i).gt.piopt(2)) piopt(2) = iavar(3,i)

         read (n4,*,iostat=ier) (idasls(j,i), j = 1, iavar(3,i))
         if (ier.ne.0) goto 99

         if (iam.eq.7.and.plopt(3)) then
            call psbtxt (i, text, iend)
            write (n8,'(i3,a,a)') i,' - ',text
         end if
c                                 make a cumulative list of stable phases
c                                 first get the number of occurrences of 
c                                 each phase in the assemblage
         nph(i) = 0

         do j = 1, k5
            idsol(j,i) = 0
            nrep(j,i) = 0
         end do 

         do j = 1, iavar(3,i)
c                                 loop over all phases
            kount = .true.

            if (j.le.iavar(1,i)) then 
c                                 a solution phase
               do k = 1, nph(i)

                  if (idsol(k,i).eq.idasls(j,i)) then 
c                                 the phase has already been found
c                                 in the assemblage, count the replicate
                     kount = .false.
                     nrep(k,i) = nrep(k,i) + 1
                     exit

                  end if

               end do
 
            end if

            if (kount) then
c                                  the phase as not yet been found in 
c                                  the assemblage.
               nph(i) = nph(i) + 1
               idsol(nph(i),i) = idasls(j,i)
               nrep(nph(i),i) = 1

            end if

         end do
c                                 at this point nph(i) is the number of 
c                                 unique phases, nrep(i) is the number of
c                                 time it is repeated.

c                                 make a global array in which each id 
c                                 occurs only once

c                                 next compare to the existing list
         do k = 1, nph(i)  

            kount = .true.

            do j = 1, istab

               if (idsol(k,i).eq.idstab(j)) then

                  if (nrep(k,i).gt.nstab(j)) nstab(j) = nrep(k,i)
                  kount = .false.
                  exit

               end if

            end do 

            if (kount) then

               istab = istab + 1
               if (istab.gt.k10) call error (999,0d0,istab,'ISTAB ')
               nstab(istab) = nrep(k,i)
               idstab(istab) = idsol(k,i)

            end if

         end do 

      end do
c                                 close n8 for assemblage list (plopt(3), PSSECT)
      close (n8) 
c                                 make the "null" assemblage
      iap(k2) = k3
      iavar(1,k3) = 0
      iavar(2,k3) = 0 
      iavar(3,k3) = 0 

      if (icopt.eq.7.and.fileio) then 
c                                 if coodinates from a file, read
c                                 coordinate file.
         open (n8,file=cfname,status='old',iostat=ier)
         if (ier.ne.0) call error (6,vip(1,1),i,cfname)
         if (loopy.gt.k2) call error (1,vip(1,1),loopy,'k2')
         do j = 1, loopy
            read (n8,*,iostat=ier) (vip(i,j), i = 1, ipot)
            if (ier.ne.0) then 
               write (*,1000) cfname
               stop
            end if 
         end do 
         close (n8)

      end if 

99    if (ier.ne.0) err = .true.

1000  format (/,'**error ver635** Coordinate file ',a,/,
     *       'is inconsistent with plot file, re-run VERTEX.',/)

      end

      subroutine redsub (jpoly,tname)
c----------------------------------------------------------------------
c subroutine to read polytope/subdivision ranges for polytope jpoly
c of a 688 solution model.
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer ier, jpoly, j, k

      character key*22, values*80, strg*80, tname*10

      double precision rnums(m4)

      character mname*8
      common/ cst18a /mname(m4)
c----------------------------------------------------------------------
      if ((poly(h0).gt.1.and.ivert(jpoly,isimp(jpoly)).gt.1)
     *                                     .or.isimp(jpoly).gt.1) then

         ier = 0
c                              reading a composite model or a prism, 
c                              a name is associated with each subdivision range
         do j = 1, isimp(jpoly)

            do k = 1, ivert(jpoly,j) - 1

               call redcd0 (n9,ier,key,values,strg)

               if (ier.ne.0) exit 

               poname(h0,jpoly,j,k) = key
               read (values,*,iostat=ier) pxmn(jpoly,j,k), 
     *              pxmx(jpoly,j,k), pxnc(jpoly,j,k), pimd(jpoly,j,k)

               if (ier.ne.0) exit

            end do

            if (ier.ne.0) exit

            call redcd0 (n9,ier,key,values,strg)

            poname(h0,jpoly,j,k) = key

         end do

         if (ier.ne.0) call error (99,0d0,k,'error while reading polyto'
     *                     //'pe/subdivision data for solution '//tname)

      else
c                              reading a simplicial model no names are
c                              read with the subdivision scheme, 
c                              so why the j-loop?
         do j = 1, isimp(jpoly)
            do k = 1, ivert(jpoly,j) - 1
               call readda (rnums,4,tname)

               poname(h0,jpoly,j,k) = 'X_'//mname(k)
               pxmn(jpoly,j,k) = rnums(1)
               pxmx(jpoly,j,k) = rnums(2)
               pxnc(jpoly,j,k) = rnums(3)
               pimd(jpoly,j,k) = idint(rnums(4))

            end do

            poname(h0,jpoly,j,k) = 'X_'//mname(k)

         end do

      end if

      end 

      subroutine endpa (ld,jd,ids)
c----------------------------------------------------------------------
c generate compositional coordinates for endmember jd of solution ids
c during adaptive optimization, ld is the associated refinement point.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ids, jd, ld

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer tphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct
c----------------------------------------------------------------------
c                                 set refinement point index
      hkp(jd) = ld

      pa(1:nstot(ids)) = 0d0
c                                 locate the endmember in the solution
      do i = 1, lstot(ids)
         if (jend(ids,2+i).eq.jd) then
            pa(i) = 1d0
            exit
         end if 
      end do

      call makepp (ids)

      end


      subroutine makapz (id)
c----------------------------------------------------------------------
c subroutine to construct the apz matrix for the p' to independent z 
c limits, where p' is the first nstot-1 elements of p.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, i, j, k, l, m, nvar, zct

      double precision dbz

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c---------------------------------------------------------
      nvar = nstot(id) - 1
c                                 to be counted:
      nz(id) = 0
c                                 for each site
      do i = 1, msite(id)
c                                 for each species
         if (zmult(id,i).eq.0d0.or.ksmod(id).eq.688) then
            zct = zsp1(id,i)
         else 
            zct = zsp(id,i)
         end if

         do j = 1, zct
c                                 initial az, bz
            nz(id) = nz(id) + 1

            dbz = 0d0
c                                 both Temkin and non-Temkin have
c                                 Az*p >= 0 constraints:
            apz(id,nz(id),1:nstot(id)) = 0d0

            do k = 1, lterm(j,i,id)

               m = ksub(k,j,i,id)
               if (m.ne.nstot(id).or..not.equimo(id)) then

                  apz(id,nz(id),m) = apz(id,nz(id),m) 
     *                               + dcoef(k,j,i,id)

               else 
c                                 decompose p(ntot) into 1 - p(1) - ...- p(nvar)
                  dbz = dcoef(k,j,i,id)

                  do l = 1, nvar

                     apz(id,nz(id),l) = apz(id,nz(id),l) - dbz

                  end do

               end if 

            end do

            zl(id,nz(id)) = - dbz
c                                 non-Temkin have the Az*p <= 1 constraint
            if (zmult(id,i).ne.0d0) then

               zu(id,nz(id)) =  1d0 - dbz

            else 

               zu(id,nz(id)) = 1d20

            end if

         end do

         if (zmult(id,i).gt.0d0.and.ksmod(id).ne.688) then 
c                                 a pre-688 model, need to make
c                                 a limit expression for the missing
c                                 site fraction

c                                 pointer to the previous constraints for
c                                 the site
            k = nz(id) - zsp(id,i) + 1

            nz(id) = nz(id) + 1
            dbz = 1d0
            apz(id,nz(id),1:nvar) = 0d0

            do j = 1, zsp(id,i)
               dbz = dbz + zl(id,k)
               do l = 1, nvar
                  apz(id,nz(id),l) = apz(id,nz(id),l) - apz(id,k,l)
               end do
               k = k + 1
            end do

            zl(id,nz(id)) =  -dbz
            zu(id,nz(id)) =  1d0 - dbz

         end if

      end do

      end

      subroutine sety2x (id)
c----------------------------------------------------------------------
c subroutine to convert independent disordered y to subcomposition
c x's, assumes y's are normalized. if the model is a CSS
c then all fraction in 0 wt prisms are zeroed and the result scanned for 
c negative fractions, if these are found bad is set to true.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ii, i, j, k, l, m, id

      double precision xt

      integer mx
      double precision ayx
      common/ csty2x /ayx(h9,h4,mst*msp,m4),mx(h9,h4)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
c                                  test
      do ii = 1, poly(id)
c                                  get the polytope weight
         if (pop1(id).eq.1) then 

            pwt(ii) = 1d0

         else

            pwt(ii) = 0d0 

            do k = pvert(id,ii,1), pvert(id,ii,2)
               pwt(ii) = pwt(ii) + y(k)
            end do

            if (dabs(pwt(ii)).lt.nopt(50)) then

               pwt(ii) = 0d0

               do k = pvert(id,ii,1), pvert(id,ii,2)
                  y(k) = 0d0
               end do

            else if (pwt(ii).gt.nopt(55)) then 

               pwt(ii) = 1d0

            end if

         end if

         l = 1
         m = 1

         do i = 1, mx(id,ii)

            xt = 0d0
c                                 the algebra
            do k = pvert(id,ii,1), pvert(id,ii,2)

               j = k - pvert(id,ii,1) + 1

               xt = xt + ayx(id,ii,i,j)*y(k)

            end do

            if (pwt(ii).gt.nopt(50)) xt = xt / pwt(ii)

            if (xt.lt.nopt(50)) then 
               xt = 0d0
            else if (xt.gt.nopt(56)) then 
               xt = 1d0
            end if

            x(ii,l,m) = xt

            m = m + 1

            if (m.gt.ispg(id,ii,l)) then 
               l = l + 1
               m = 1
            end if

         end do

      end do
c                                 set the prism weights, just in case
c                                 they are used, which i doubt.
      if (pop1(id).gt.1) then 
         do k = 1, poly(id)
            x(pop1(id),1,k) = pwt(k)
         end do
      end if

      end


      subroutine p2zind (p,z,l,ids)
c----------------------------------------------------------------------
c subroutine to compute independent site fractions (or molar amounts 
c for temkin sites) and load them sequentially into the 1d array z
c
c     l - the total number of indpendent site fractions
c non-temkin models:
c     zsp - number of independent site fractions for site i (zsp1 - 1)
c     z(l) - molar site fraction of species j on site i
c temkin models:
c     zsp - number of species (zsp1) for site i
c     z(l) - molar amount of species j on site i
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision p(*), z(*)

      integer i,j,k,l,ids

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
c                                 for each site
      l = 0

      do i = 1, msite(ids)
c                                 get site fractions
         do j = 1, zsp(ids,i)

            l = l + 1

            z(l) = dcoef(0,j,i,ids)
c                                 for each term:
            do k = 1, lterm(j,i,ids)

               z(l) = z(l) + dcoef(k,j,i,ids) * p(ksub(k,j,i,ids))

            end do

         end do

      end do

      end

      subroutine p2zall (y,z,ldz,ids)
c----------------------------------------------------------------------
c subroutine to compute all site fractions (or molar amounts for temkin
c sites) computed from equations
c
c non-temkin models:
c     zsp1 - number of site fractions
c     zsp - number of independent site fractions zsp-1
c     z(i,j) - molar site fraction of species j on site i
c temkin models:
c     zsp1 - number of species
c     zsp - zsp1
c     z(i,j) - molar amount of species j on site i
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,ldz,ids

      double precision y(*), zt, z(ldz,*)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
c                                 for each site
      do i = 1, msite(ids)

         zt = 0d0

         if (zmult(ids,i).ne.0d0.and.ksmod(ids).ne.688) then
c                                 get site fractions
            do j = 1, zsp(ids,i)

               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) +
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do

               zt = zt + z(i,j)

            end do

            z(i,j) = 1d0 - zt

         else if (zsp1(ids,i).gt.1) then
c                                 temkin or 688 model format, all species fractions are available
            do j = 1, zsp1(ids,i)
c                                 molar site population
               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) + 
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do

            end do

         end if

      end do

      end

      subroutine setder (ids,tname)
c---------------------------------------------------------------------
c evaluate coefficients for diff(g,p') computed by getder as called by
c minfrc.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, ind, i, j, k, l, ntot, nvar

      character tname*(*), reason*20

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c                                 configurational entropy variables:
      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      if (ksmod(ids).eq.0.or.
     *    (ksmod(ids).ge.20.and.ksmod(ids).le.50)) then 

          deriv(ids) = .false.
          reason =  'special case'

      else if (extyp(ids).eq.1) then 

          deriv(ids) = .false.
          reason = 'redlich-kistler ex'

      else if (.not.equimo(ids)) then

          deriv(ids) = .false.
          reason = 'non-equimolar O/D'

      else

          deriv(ids) = .true.

      end if

      if (.not.deriv(ids)) then
         if (iam.lt.3) write (*,1000) tname, reason
         return
      end if

      ntot = nstot(ids)
      nvar = ntot - 1
c----------------------------------------------------------------------
c                                 configurational negentropy derivatives:
c                                 for each site
      do i = 1, msite(ids)
c                                 site fraction derivatives
         dzdp(1:zsp1(ids,i),i,1:nvar,ids) = 0d0
c                                 for each species
         do j = 1, zsp(ids,i)
c                                 for each term:
            do k = 1, lterm(j,i,ids)
c                                 endmember index
               ind = ksub(k,j,i,ids)

               if (ind.le.nvar) then 
c                                 coefficient of endmembers in dz/dp
                  dzdp(j,i,ind,ids) = dzdp(j,i,ind,ids) 
     *                                + dcoef(k,j,i,ids)
               else
c                                  the eliminated endmember contributes to 
c                                  to all remaining endmembers
                  do l = 1, nvar
                     dzdp(j,i,l,ids) = dzdp(j,i,l,ids) 
     *                                 - dcoef(k,j,i,ids)
                  end do

               end if

            end do

         end do
c                                need dzdp for j = zsp + 1 species
c                                this would not be necessary for
c                                688 format if the zsp1(ids,i) 
c                                counter were used above.
         do l = 1, nvar
            do k = 1, zsp(ids,i)
               dzdp(j,i,l,ids) = dzdp(j,i,l,ids) - dzdp(k,i,l,ids) 
            end do
         end do

         if (zmult(ids,i).ne.0d0) then
c                                scale the derivatives by r*multiplicity
            do j = 1, zsp(ids,i) + 1
               dzdp(j,i,1:nvar,ids) = dzdp(j,i,1:nvar,ids)*zmult(ids,i)
            end do

         end if

      end do
c----------------------------------------------------------------------
c                                 endmember configurational derivatives
      do l = 1, nvar
c                                 these are negentropy derivatives
         ds0dp(l,ids) = scoef(l,ids) - scoef(ntot,ids)

      end do
c----------------------------------------------------------------------
c                                 endmember configurational derivatives
      dcdp(1:icp,1:nvar,ids) = 0d0
c                                 solutions with no dependent endmembers:
c                                 pa coordinates used to compute the composition
      do i = 1, nvar
         do j = 1, icp
            dcdp(j,i,ids) = endc(ids,i,j) - endc(ids,ntot,j)
         end do
      end do
c----------------------------------------------------------------------
c                                 excess function derivatives
      do i = 1, jterm(ids)

         do j = 1, rko(i,ids)

            k = jsub(j,i,ids)

            do l = 1, nvar

              if (k.eq.l) then
                  dgex(l,j,i,ids) = 1d0
              else if (k.eq.ntot) then
                  dgex(l,j,i,ids) = -1d0
               else 
                  dgex(l,j,i,ids) = 0d0
               end if

            end do 

         end do

      end do

1000  format (/,'**warning ver212** No MINFRC derivatives for: ',a,/,
     *          'Reason: ',a,/)

      end

      subroutine getder (g,dgdp,ids)
c---------------------------------------------------------------------
c compute g(p') and diff(g,p') for minfrc.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, l, ntot, nvar

      double precision g, dgdp(*), gx, dgxdp(m14)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision p,t,xco2,mu1,mu2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,mu1,mu2,tr,pr,r,ps
c----------------------------------------------------------------------
      ntot = nstot(ids)
      nvar = ntot - 1
c                                 get the bulk composition needed for 
c                                 leveling.
      call getscp (rcp,rsum,rids,rids)

      g = 0d0
      dgdp(1:nvar) = 0d0
c                                 configurational negentropy and derivatives
      call p2sds (g,dgdp,nvar,ids)
c                                 correct derivatives for mechanical configurational 
c                                 negentropy, multiply by t to convert to configurational
c                                 gibbs energy
      do l = 1, ntot
         g = g + pa(l)*scoef(l,ids)
         if (l.gt.nvar) exit
         dgdp(l) = t*(dgdp(l) + ds0dp(l,ids))
      end do
c                                 compute excess gibbs energy and it's
c                                 derivatives are added to dgdp
      call p2gdg (gx,dgxdp,nvar,ntot,ids)
c                                 at this point dsdp is really -dsdp (entropy units).
c                                 add in excess as well
      g = t*g + gx
c                                 add mechanical mix and derivatives
      do l = 1, ntot
         g = g + pa(l) * gend(l)
         if (l.gt.nvar) exit
         dgdp(l) = dgdp(l) + dgxdp(l) + gend(l) - gend(ntot)
      end do

      end

      subroutine ingend (id)
c-----------------------------------------------------------------------
c make a generic endmember g0' array for solution id, used only by gsol2
c for minfrc. g0' = g0(p,t) + dqf(p,t)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,i,k,l,ind

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision g
      common/ cst2 /g(k1)

      double precision enth
      common/ cxt35 /enth(j3)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c----------------------------------------------------------------------
      do i = 1, lstot(id)
         gend(i) = g(jend(id,2+i))
      end do
c                                 add in the dqf's, these can only
c                                 be for the first lstot endmembers:
      do i = 1, jdqf(id)
         gend(iq(i)) = gend(iq(i)) + dq(i)
      end do
c                                 make the ordered endmembers:
      do k = 1, nord(id)

         i = lstot(id) + k

         gend(i) = enth(k)

         do l = 1, nrct(k,id)

            ind = ideps(l,k,id)

            gend(i) = gend(i) - dydy(ind,k,id) * gend(ind)

         end do

      end do

      end

      subroutine p2gdg (g,dgdp,nvar,ntot,ids)
c----------------------------------------------------------------------
c subroutine to compute the excess energy and its derivatives with
c respect to the p' endmember fractions
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision dgdp(*), g, t, tcum

      integer i, j, l, k, nvar, ntot, ids

      double precision z, pa, p0a, x, w, yy, wl, pp
      common/ cxt7 /yy(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)
c                                 local alpha
      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c----------------------------------------------------------------------
      g = 0d0
      dgdp(1:nvar) = 0d0

      do i = 1, jterm(ids)

         t = 1d0

         do j = 1, rko(i,ids)
            t = t * pa(jsub(j,i,ids))
         end do

         g = g + w(i) * t

         do l = 1, nvar

            tcum = 0d0

            do j = 1, rko(i,ids)

               t = 1d0

               do k = 1, rko(i,ids)

                  if (k.eq.j) then

                     t = t * dgex(l,k,i,ids)

                  else 

                     t = t * pa(jsub(k,i,ids))

                  end if

                  if (t.eq.0d0) exit

               end do

               tcum = tcum + t

            end do

            dgdp(l) = dgdp(l) + w(i) * tcum

         end do

      end do

      if (llaar(ids)) then

         t = 0d0 

         do i = 1, ntot

            t = t + alpha(i)*pa(i)

         end do

         g = g / t

         do i = 1, nvar

            dgdp(i) = (dgdp(i) - g * (alpha(i)-alpha(ntot))) / t

         end do

      end if

      end

      subroutine p2sds (s,dsdp,nvar,ids)
c----------------------------------------------------------------------
c subroutine to configurational negentropy  and its derivatives with
c respect to the endmember fractions p
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision zt, z(m11), dsdp(*), dzlnz, s, zlnz, lnz

      integer i, j, k, l, ids, nvar

      double precision p,t,xco2,mu1,mu2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,mu1,mu2,tr,pr,r,ps

      double precision zz, pa, p0a, x, w, yy, wl, pp
      common/ cxt7 /yy(m4),zz(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
c                                 for each site
      do i = 1, msite(ids)

         zt = 0d0
         zlnz = 0d0

         if (zmult(ids,i).ne.0d0) then
c                                 non-temkin:
c                                 get site fractions
            do j = 1, zsp(ids,i)

               z(j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)
                  z(j) = z(j) + dcoef(k,j,i,ids) * pa(ksub(k,j,i,ids))
               end do

               zt = zt + z(j)

               call ckzlnz (z(j),zlnz)

               lnz = 1d0 + dlog(z(j))

               do l = 1, nvar

                  dsdp(l) = dsdp(l) + dzdp(j,i,l,ids) * lnz

               end do

            end do

            zt = 1d0 - zt

            call ckzlnz (zt,zlnz)
c                                 site negentropy
            zlnz = zmult(ids,i) * zlnz
c                                 for non-temkin sites dzdp is already scaled by q*R
            lnz = 1d0 + dlog(zt)

            do l = 1, nvar

               dsdp(l) = dsdp(l) + dzdp(j,i,l,ids) * lnz

            end do

         else 
c                                 temkin
            do j = 1, zsp(ids,i)
c                                 molar site population
               z(j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)
                  z(j) = z(j) + dcoef(k,j,i,ids) * pa(ksub(k,j,i,ids))
               end do
c                                 zt is the multiplicity here
               zt = zt + z(j)

            end do
c                                 site doesn't exist if zt = 0
            if (zt.lt.nopt(50)) cycle
c                                 convert molar amounts to fractions
            z(1:zsp(ids,i)) = z(1:zsp(ids,i)) / zt

            do j = 1, zsp(ids,i)
c                                 ckzlnz sets z < nopt(50) = nopt(50)
               call ckzlnz (z(j),zlnz)

            end do
c                                 site negentropy
            zlnz = r * zt * zlnz
c                                 derivatives
            do l = 1, nvar

               dzlnz = 0d0
c                                 for each species
               do j = 1, zsp(ids,i)

                  dzlnz = dzlnz +  dzdp(j,i,l,ids) * dlog(z(j))

               end do

               dsdp(l) = dsdp(l) + r*dzlnz

            end do

         end if

         s = s + zlnz

      end do

      end

      subroutine ckdzlz (z,zlnz,dlnz)
c----------------------------------------------------------------------
c subroutine to test and accumulate z*ln(z) and set 1 + ln(z)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision z, zlnz, dlnz
c----------------------------------------------------------------------
      if (z.gt.1d0) then

         z = 1d0

      else if (z.lt.nopt(50)) then

         z = nopt(50)

      end if

      dlnz = dlog(z) + 1d0
      zlnz = zlnz + z * dlog(z)

      end 

      subroutine ckzlnz (z,zlnz)
c----------------------------------------------------------------------
c subroutine to test/reset site fraction value z and accumulate z*ln(z)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision z, zlnz
c----------------------------------------------------------------------
      if (z.gt.1d0) then
         z = 1d0
      else if (z.lt.nopt(50)) then
         z = nopt(50)
      end if

      zlnz = zlnz + z * dlog(z)

      end


      subroutine makayz (id)
c----------------------------------------------------------------------
c subroutine to make the ayz matrix for ayz*y = z, z is the independent
c subset of the site fractions.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, id, nz

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision ayz
      common/ csty2z /ayz(h9,m20,m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c----------------------------------------------------------------------
      ayz(id,:,:) = 0d0
c                                 independent endmembers:
      do k = 1, lstot(id)
c                                 endmember is in column knsp(k,id)
         pa(:) = 0d0
         pa(k) = 1d0
         call p2zind (pa,z,nz,id)
         ayz(id,1:nz,knsp(k,id)) = z(1:nz)

      end do
c                                 dependent endmembers:
      do k = 1, ndep(id)
c                                 the index of the dependent endmember
c                                 in y is
         j = knsp(lstot(id)+k,id)

         do l = 1, ndph(k)
c                                 the dependent disordered endmember decomposes
c                                 to independent disordered endmember idep(k,l):
c                                 this is insanely inefficient, but who cares?
            pa(:) = 0d0
            pa(iy2p(idep(k,l))) = 1d0
            call p2zind (pa,z,nz,id)

            do i = 1, nz
               ayz(id,i,j) = ayz(id,i,j) + nu(k,l)*z(i)
            end do

         end do

      end do

      end

      subroutine makayc (id)
c----------------------------------------------------------------------
c subroutine to make the ayc matrix for ayc*y = c
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, id

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision ayc
      common/ csty2c /ayc(h9,k5,m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)
c----------------------------------------------------------------------
      ayc(id,:,:) = 0d0
c                                 independent endmembers:
      do k = 1, lstot(id)
c                                 endmember is in column j
         j = knsp(k,id)

         do i = 1, icomp
            ayc(id,i,j) = endc(id,k,i)
         end do

      end do
c                                 dependent endmembers:
      do k = 1, ndep(id)
c                                 the index of the dependent endmember
c                                 in y is
         j = knsp(lstot(id)+k,id)

         do l = 1, ndph(k)
c                                 the dependent endmember decomposes
c                                 to independent endmember idep(k,l):
            do i = 1, icomp
               ayc(id,i,j) = ayc(id,i,j) + 
     *                       nu(k,l) * endc(id,iy2p(idep(k,l)),i)
            end do

         end do

      end do

      end

      subroutine makapc (id)
c----------------------------------------------------------------------
c subroutine to make the ap'*c matrix for ap'*c = c, where p' is the 
c the independent subset of p. c is the normalized composition.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------

      do j = 1, nstot(id)

         do i = 1, icomp
            apc(id,i,j) = endc(id,j,i)
c           apc(id,i,j) = endc(id,j,i) / endt(id,j)
         end do

      end do
c                                   eliminate p(nstot) as 1 - sum(p,nstot-1)
c                                   this method costs an extra column dimension in
c                                   p2c, but what the heck... it can be used 
c                                   to form the constraint b vector
      do j = 1, nstot(id) - 1
         do i = 1, icomp + 1
            apc(id,i,j) = apc(id,i,j) - apc(id,i,nstot(id))
         end do
      end do

      end

      subroutine setord (im)
c---------------------------------------------------------------------
c set global order/disorder models parameters, call by gmodel for 
c solution model im. 
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical zbad

      integer im, i, j, ind, id, k, l,itic, ii, imatch, 
     *        il, ik, kk, jp1

      double precision dzt, dx, delta, c0(0:20), c1(0:20), zij

      external zbad

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision p,t,xco2,mu1,mu2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,mu1,mu2,tr,pr,r,ps

      integer nsub,nterm
      double precision acoef
      common/ cst107 /acoef(m10,m11,0:m0),
     *                nterm(m10,m11),nsub(m10,m11,m0)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m14+2)

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)

      character specie*4
      integer jsp, ins
      common/ cxt33 /jsp,ins(nsp),specie(nsp)

      character mname*8
      common/ cst18a /mname(m4)

      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c----------------------------------------------------------------------
      equimo(im) = .true.
c                                 models with speciation:
      do j = 1, norder

         do i = 1, 3
            deph(i,j,im) = denth(j,i)
         end do

         nrct(j,im) = nr(j)

         do i = 1, nr(j)
            ideps(i,j,im) = iy2p(iddeps(i,j))
         end do

      end do
c                                 classify multiple species models according
c                                 to whether the disordered reactants are
c                                 partially or completely correlated, assume
c                                 anti-correlation is not a possible case.
      icase(im) = 0

      if (norder.gt.1) then

         imatch = 0

         do j = 1, nr(1)

            id = ideps(j,1,im)

            do i = 1, nr(2)
               if (id.eq.ideps(i,2,im)) then
                  imatch = imatch + 1
                  exit
               end if
            end do

         end do

         if (imatch.eq.1) then
c                                 if match = 1 one species didn't match
c                                 assume partial correlation
            icase(im) = 2
         else if (imatch.ge.2) then
            icase(im) = 1
         end if

      end if
c                                first create derivatives of endmember
c                                fractions with respect to the ordered
c                                species:
      do j = 1, norder

         do i = 1, nstot(im)
            dydy(i,j,im) = 0d0
         end do
c                                derivative of the ordered species with
c                                respect to itself:
         dydy(kstot+j,j,im) = 1d0
c                                each ordered species decomposes to
c                                two disordered species iddeps(1-2,j)
c                                depnu is the stoichiometric coefficient
c                                of the disordered species in the ordered
c                                species.

c                                derivatives of the consituent species
c                                with respect to the ordered species
         dnu(j,im) = 1d0

         do i = 1, nr(j)
            dydy(ideps(i,j,im),j,im) = dydy(ideps(i,j,im),j,im)
     *                                  - depnu(i,j)
            dnu(j,im) = dnu(j,im) + dydy(ideps(i,j,im),j,im)
         end do
c                                dnu =~ 0 => speciation reaction is not equimolar
         if (dabs(dnu(j,im)).lt.nopt(50)) then
            dnu(j,im) = 0d0
         else 
            equimo(im) = .false.
         end if

      end do
c                                evaluate the second derivative of each
c                                pi*pj term in the excess function with
c                                respect to kth species
      do i = 1, iterm
         do j = 1, norder
            do k = 1, norder

                  dppp(k,j,i,im) =  dydy(jsub(1,i,im),k,im)
     *                             *dydy(jsub(2,i,im),j,im)
     *                           +
     *                              dydy(jsub(2,i,im),k,im)
     *                             *dydy(jsub(1,i,im),j,im)

            end do
         end do
      end do
c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      do i = 1, msite(im)
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         if (zsp(im,i)+1.gt.m11) call error (1,dx,zsp(im,i)+1,'m11')

         do k = 1, norder
            sdzdp(k,zsp(im,i)+1,i,im) = 0d0
         end do

         do j = 1, zsp(im,i)
c                                 # of terms in the
c                                 site fraction function and a0.
            do l = 1, norder
               sdzdp(l,j,i,im) = 0d0
            end do
c                                 for each term:
            do k = 1, lterm(j,i,im)
c                                 endmember indexes
               ind = ksub(k,j,i,im)
c                                 get derivatives, of species fractions
c                                 with respect to ordered species
               do l = 1, norder
                  itic = 0
                  do ii = 1, nr(l)
                     if (ind.eq.ideps(ii,l,im)) then
                        sdzdp(l,j,i,im) = sdzdp(l,j,i,im)
     *                  + dydy(ideps(ii,l,im),l,im)*dcoef(k,j,i,im)
                        itic = itic + 1
c                                 high order terms not allowed
                        if (itic.gt.1) call error (999,r,801,'SETORD')
                     end if
                  end do
c                                 the derivative of a term with the
c                                 ordered species.
                  if (ind.eq.kstot+l)
     *               sdzdp(l,j,i,im) = sdzdp(l,j,i,im)
     *                               + dcoef(k,j,i,im)

               end do
            end do
         end do
      end do
c                                 multiply each dzdp by qmult (R*site
c                                 multiplicity) to reduce operation
c                                 count in evaluating derivatives.
      do k = 1, norder
         do i = 1, msite(h0)

            dzt = 0d0

            do j = 1, zsp(im,i)
               if (dabs(sdzdp(k,j,i,im)).lt.nopt(50)) 
     *                                       sdzdp(k,j,i,im) = 0d0
               dzt = dzt + sdzdp(k,j,i,im)
            end do

            if (dabs(dzt).lt.nopt(50)) dzt = 0d0
            sdzdp(k,j,i,im) = -dzt

         end do
      end do
c                                 ----------------------------------------------
c                                 derive z2p limit expressions for O/D models
      do k = 1, nord(im)
c                                 number of limits for ordered species k
         ln(k,im) = 0
      end do

      if (ksmod(im).ne.688) then 

         if (nord(im).gt.1) call error (72,c0(0),i,
     *            'solution '//tname//': multiple order parameters '//
     *            'only allowed in 688 format models')

c                                 all z expressions may be necessary to
c                                 formulate limits, make the ksp'th + 1
c                                 species expression by differnce
         do i = 1, msite(im)
c                                 qmult = 0, temkin, all expressions are
c                                 available
            if (zmult(im,i).eq.0d0) cycle

            jp1 = zsp(im,i) + 1
c                                 initialize the term counter
            lterm(jp1,i,im) = 0
c                                 cycle through the endmembers to work out
c                                 if it has a non zero fraction
            do l = 1, nstot(im)

               call zmake (zij,i,l,im)

               if (dabs(zij).lt.nopt(50)) cycle

               lterm(jp1,i,im) = lterm(jp1,i,im) + 1
               dcoef(lterm(jp1,i,im),jp1,i,im) = zij
               ksub(lterm(jp1,i,im),jp1,i,im) = l

            end do
         end do
      end if 

      do i = 1, msite(im)

         if (zmult(im,i).eq.0d0) then
            jp1 = 0
         else
            jp1 = 1
         end if

         do j = 1, zsp(im,i) + jp1

            do k = 1, nord(im)

               c0 = 0d0
               c1 = 0d0

               if (dcoef(0,j,i,im).ne.0d0) call error (72,c0(0),i,
     *            'solution '//tname//': constants not allowed in '//
     *            'O/D model configurational entropy site fraction '//
     *            'expressions')

               do l = 1, lterm(j,i,im)

                  il = ksub(l,j,i,im)

                  if (il.le.lstot(im)) then

                     c0(il) = c0(il) + dcoef(l,j,i,im)

                     do ik = 1, nord(im)

                        kk = lstot(im) + ik
c                                  coefficient on p0
                        c0(kk) = c0(kk) - dydy(il,ik,im)
     *                                  * dcoef(l,j,i,im)
c                                  coefficient on p
                        c1(kk) = c1(kk) + dydy(il,ik,im)
     *                                  * dcoef(l,j,i,im)

                     end do

                  else

                     c1(il) = c1(il) + dcoef(l,j,i,im)

                  end if
c                                 at this point the c's are the coefficients for a
c                                 z(p,p0), below they are rearranged to get p(kk) = f(p0,p[~kk],z[0,1])
c                                 in other words the loop on mord(im) is superfluous, but what the heck...
                  kk = lstot(im)+k

                  if (l.eq.lterm(j,i,im).and.dabs(c1(kk)).gt.nopt(50))
     *                                                              then

                     do ik = 0, nstot(im)
                        c0(ik) = -c0(ik)/c1(kk)
                     end do
c                                the constant for the p(k) limit when z(j) = 1
                     c1(0) = c0(0) + 1d0/c1(kk)

                     do ik = lstot(im) + 1, nstot(im)
                        if (ik.eq.kk) cycle
                        c1(ik) = -c1(ik)/c1(kk)
                     end do

                     c1(kk) = 0d0

                     if  (c1(0).gt.c0(0)) then
c                                 z = 1 is the upper limit, the constant is c0(0) and
                        delta = c1(0)-c0(0)

                     else
c                                 z = 0 is the upper limit, the constant is c1(0) and
                        delta = c0(0)-c1(0)
                        c0(0) = c1(0)

                     end if
c                                 -----------------------------------------------------
c                                 found a limit, increment limit counter
                     ln(k,im) = ln(k,im) + 1
c                                 initialize p0 term counter for limit
                     lt(ln(k,im),k,im) = 0
c                                 load the p0 coefficients, if simplicial p0 > lstot(im) = 0
                     do ik = 1, nstot(im)

                        if (dabs(c0(ik)).lt.nopt(50)) cycle
c                                 increment term counter:
                        lt(ln(k,im),k,im) = lt(ln(k,im),k,im) + 1
c                                 save the coefficient and index:
                        lc(lt(ln(k,im),k,im),ln(k,im),k,im) = c0(ik)
                        lid(lt(ln(k,im),k,im),ln(k,im),k,im) = ik

                     end do
c                                 initialize p term counter for limit
                     jt(ln(k,im),k,im) = 0
c                                 load the p coefficients
                     do ik = lstot(im) + 1, nstot(im)

                        if (dabs(c0(ik)).lt.nopt(50)) cycle
c                                 increment term counter:
                        jt(ln(k,im),k,im) = jt(ln(k,im),k,im) + 1
c                                 save the coefficient and index:
                        jc(jt(ln(k,im),k,im),ln(k,im),k,im) = c1(ik)
                        jid(jt(ln(k,im),k,im),ln(k,im),k,im) = ik

                     end do
c                                load the constant and delta:
                     l0c(1,ln(k,im),k,im) = c0(0)
                     l0c(2,ln(k,im),k,im) = delta

                  end if

               end do

            end do

         end do

      end do

      end

      subroutine makayx (id)
c----------------------------------------------------------------------
c subroutine to construct the ayx matrices for the y to x conversion
c of each subcomposition.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ii, i, j, k, l, m, id

      double precision xt

      integer mx
      double precision ayx
      common/ csty2x /ayx(h9,h4,mst*msp,m4),mx(h9,h4)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      do ii = 1, pop1(id)
c                                 for subcomposition ii 
c                                 the Ax matrix in Ax*y = x
c                                 m = sum( ispg(1:istg) )
c                                 n = pvert(2) - pvert(1) + 1

         mx(id,ii) = 0

         do i = 1, istg(id,ii)
            mx(id,ii) = mx(id,ii) + ispg(id,ii,i)
         end do 

         do i = 1, ncoor(id)
            do j = 1, pvert(id,ii,2) - pvert(id,ii,2) + 1
               ayx(id,ii,i,j) = 0d0 
            end do 
         end do

         do k = pvert(id,ii,1), pvert(id,ii,2)
c                                 for each endmember in the subcomposition
c                                 load the column of ax:
            j = k - pvert(id,ii,1) + 1

            i = 0

            do l = 1, istg(id,ii)

               do m = 1, ispg(id,ii,l)
c                                 kmsol indicates the species on 
c                                 the site m of endmember k
                  if (kmsol(id,k,l).eq.m) then

                     ayx(id,ii,i+m,j) = 1d0

                     exit 

                  end if

               end do

               i = i + ispg(id,ii,l)

            end do 

         end do 

      end do
c                                  test
      do ii = 1, poly(id)
c                                  get the polytope weight
         if (pop1(id).eq.1) then 

            pwt(ii) = 1d0

         else

            pwt(ii) = 0d0 

            do k = pvert(id,ii,1), pvert(id,ii,2)
               pwt(ii) = pwt(ii) + y(k)
            end do

         end if

         l = 1
         m = 1

         do i = 1, mx(id,ii)

            xt = 0d0
c                                 the algebra
            do k = pvert(id,ii,1), pvert(id,ii,2)

               j = k - pvert(id,ii,1) + 1

               xt = xt + ayx(id,ii,i,j)*y(k)

            end do

c           write (*,1000) ii,l,m,xt,x(ii,l,m)

            m = m + 1

            if (m.gt.ispg(id,ii,l)) then 
               l = l + 1
               m = 1
            end if

         end do

      end do

1000  format (3(i2,1x),3(3x,g14.6))

      end

      subroutine interm (finish,err)
c-----------------------------------------------------------------------
c if finish (only vertex) close plt/blk and delete interim results else 
c if ~finish open/read plt/blk files for PSSECT, UNSPLT, and WERAMI.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character yes*1, text*3, name*100

      integer ier, jnd(12,2), i, j, ind1, ind2

      logical err, finish, inter

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      if (finish) then

         close (n4)
         close (n5)

         if (iopt(34).eq.1) then 
c                                 delete interim results
            call mertxt (tfname,prject,'.irf',0)
            open (1000, file = tfname, status = 'old', iostat = ier)
            if (ier.ne.0) return

            do

               read (1000,*,iostat=ier) i,j
c                                 file is in use or end of irf file
               if (ier.ne.0) exit 
c                                 make the root
               write (text,'(a,i1,i1)') '_',i,j
               call mertxt (name,prject,text,0)

               call mertxt (tfname,name,'.plt',0)
               open (1001, file = tfname, status = 'old', iostat = ier)
               if (ier.ne.0) exit 
               close (1001, status = 'delete')

               call mertxt (tfname,name,'.blk',0)
               open (1001, file = tfname, status = 'old', iostat = ier)
               if (ier.ne.0) exit 
               close (1001, status = 'delete')

            end do

            close (1000, status = 'delete')

         end if 

         return

      end if 

      if (iopt(34).ne.2.or.icopt.ne.5.or.iam.eq.14) then 
c                                 for all calculations other than 2d gridded 
c                                 min OR if interim_results (iopt(34)) < 2
c                                 try to open final plt and blk files
         name = prject

         call redplt (name,err)

         if (err) then

            if (iam.eq.14) then

               return

            else if (icopt.ne.5.or.iopt(34).eq.0) then 

               call error (72,nopt(1),i,'missing/corrupt plt/blk files '
     *                     //'VERTEX may still be running or the files'
     *                     //' are locked by another program')

            else 

                call warn (99,nopt(1),i,'error occurred while attemptin'
     *          //'g to read final plt/blk files; looking for interim '
     *          //'results...')

            end if

         else 

            return

         end if

      end if
c                                 the only paths here are
c                                 1) iopt(34) = 2 => man.
c                                 2) iopt(34) = 1 => auto and no final results.
c                                 3) iopt(34) > 0 => off
      inter = .false.
c                                 reset err 
      err = .false.
c                                 only hope is interim results:
      call mertxt (tfname,prject,'.irf',0)
      open (1000, file = tfname, status = 'old', iostat = ier)

      if (ier.ne.0) then 

         if (iopt(34).eq.1) then 
c                                  end of the line
            call error (72,nopt(1),i,'no IRF file: interim '//
     *                               'results are not available')
         else 
c                                  maybe the user deleted the irf file
            call warn (99,nopt(1),i,'no IRF file: interim '//
     *                              'results are not available')

            i = 0 

         end if

      else 
c                                 make a list of the result files
         i = 1
c                                 make a list of the result files
         do

            read (1000,*,iostat=ier) jnd(i,1),jnd(i,2)

            if (ier.ne.0) then

               if (i.eq.1) then

                  call error (72,nopt(1),i,'empty IRF file: interim '//
     *                                     'results are not available')

               else

                  i = i - 1
                  exit

               end if

            end if

            i = i + 1

         end do

      end if 

      if (iopt(34).eq.1) then 

         if (i.eq.0) then 
            write (*,'(a)') 'VERTEX has not completed the calculation '
     *                    //'and no interim results are available.'

            stop

         end if 
c                                 interim_results is auto, and the final results
c                                 are not available, find/use last interim result:
         write (*,'(a,/,a)') 'VERTEX has not completed the calculation'
     *         //', continue with the','latest interim result (Y/N)?'

         if (refine.and.jnd(i,1).eq.0) write (*,'(2(/,a))')
     *      'WARNING: VERTEX is currently in, or was interrupted '//
     *      'during, the auto-refine stage, but the','latest interim '//
     *      'result is from the exploratory stage, the result may be '//
     *      'inconsistent or unreadable.'

         read (*,'(a)') yes

         if (yes.ne.'y'.and.yes.ne.'Y') then 
            stop
         else if (refine.and.jnd(i,1).eq.0) then 
c                                 try reading solutions without refine data
            write (*,'(/,3(a,/))')
     *            'If an error follows change T to F in the TOF file '//
     *            'and restart PSSECT.'

         end if

         inter = .true.

         write (text,'(a,i1,i1)') '_',jnd(i,1),jnd(i,2)
         call mertxt (name,prject,text,0)

      else
c                                 if here must be auto and an irf file exists
         if (i.gt.0) then 

            write (*,'(a)') 'Do you want to plot/analyze interim '//
     *                        'results (Y/N)?'
            read (*,'(a)') yes

            if (yes.eq.'y'.or.yes.eq.'Y') then
c                                 use intermediate results
               write (*,'(/,a,/)') 'Choose from the following interim'//
     *                             ' results [default is the last]:'

               do j = 1, i 

                  if (jnd(j,1).eq.0) then

                     write (*,'(4x,i1,a,i1)') j,
     *                      ' - exploratory stage, grid level ',jnd(j,2)
                  else

                     write (*,'(4x,i1,a,i1)') j,
     *                      ' - auto-refine stage, grid level ',jnd(j,2)

                  end if

               end do

               call rdnumb (nopt(1),0d0,i,i,.false.)
               write (*,'(/)')

               ind1 = jnd(i,1)
               ind2 = jnd(i,2)

               if (refine.and.ind1.eq.0) write (*,'(3(a,/))')
     *            'WARNING: VERTEX is in, or has completed, the '//
     *            'auto-refine stage, interim results ',
     *            'from the exploratory stage may be '//
     *            'inconsistent or unreadable.','if VERTEX has been '//
     *            'terminated and the next message is **error ver072'//
     *            '**, then edit T to F in the TOF file'

               write (text,'(a,i1,i1)') '_',ind1, ind2
               call mertxt (name,prject,text,0)
               inter = .true.

            else 

               name = prject

            end if

         else
c                                 use final results
            name = prject

         end if

      end if

      call redplt (name,err)

      if (err) then
         if (inter) then 
            call error (72,nopt(1),i,'corrupt interim results, '//
     *                             'use auto-refine stage results.')
         else
            call error (72,nopt(1),i,'missing/corrupt plt/blk files '
     *                     //'VERTEX may still be running or the files'
     *                     //' are locked by another program')
         end if
      end if

      end

      subroutine redplt (name,err)
c-----------------------------------------------------------------------
c open/read plt/blk files for PSSECT and WERAMI.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character name*100

      integer ier

      logical err

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      err = .false.
c                                 open the plot file
      call mertxt (tfname,name,'.plt',0)
      open (n4, file = tfname, iostat = ier, status = 'old')
      if (ier.ne.0) then 
         err = .true.
         return
      end if 
c                                 open assemblage file
      call mertxt (tfname,name,'.blk',0)
      open (n5, file = tfname, iostat = ier, status = 'old')
      if (ier.ne.0) then 
         err = .true.
         return
      end if 
c                                 read grid data:
      call plinp (err)
      if (err) return
c                                 read bulk composition data:
      call bplinp (err)

      end

      double precision function gfesi (y,g1,g2)
c-----------------------------------------------------------------------
c gfesi returns the Gibbs free energy for BCC FeSi alloy after
c Lacaze & Sundman 1990. See FeSiBCC.mws.

c    y   - the bulk Fe mole fraction
c    g01 - free energy of Bcc Fe, without Gmag
c    g02 - free energy of Bcc Si
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical done

      integer itic

      double precision g1, g2, y, x, w0, w1, w2, rt, dg, xmin,
     *                 d2g, gord, xmax, dx, gfesi0, g0, g12, gmag

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save w1, w2, gord
      data w1, w2, gord/-11544d0, 3890d0, -10475.64d0/
c----------------------------------------------------------------------

c      g1 = g1p - gmag(1d0)

      if (y.le.nopt(50).or.y.ge.nopt(56)) then
c                                 endmember compositions, no order possible
         gfesi = y*g1 + (1d0-y)*g2 + gmag(y)
         return
      end if

c!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!

      w0  = -27809d0 + 11.62d0 * t
      gord = -10475.64d0*2 + ( (g1 + g2)/2d0 + w0)
      rt  = r*t
      g12 = 2d0*(gord - w0) - g1 - g2
c                                 max concentration of ordered species
      if (y.gt.0.5d0) then
         xmax = 1d0
c                                 the true xmin (commented) allows for
c                                 anti-ordering, but because the model
c                                 is symmetric, i up xmin to y
c        xmin = 2d0*(y-.5d0)
      else
         xmax = 2d0*y
c        xmin = 0d0
      end if

      xmax = xmax - nopt(50)
      xmin = y + nopt(50)
      x = xmax
c                                 get 1st and 2nd derivatives
      call dgfesi (dg,d2g,y,x,g12,rt)

      done = .false.
c                                 find starting point for newton-raphson
c                                 search
      if (dg.gt.0d0.and.d2g.gt.0d0) then
c                                 the max order concentration is a
c                                 good starting point
         dx = -dg/d2g

      else if (dg.lt.0d0) then
c                                 the max order is a minimum
         x = y
         done = .true.

      else
c                                 try the max disordered concentration
         x = xmin

         call dgfesi (dg,d2g,y,x,g12,rt)

         if (dg.lt.0d0.and.d2g.gt.0d0) then
c                                 ok
            dx = -dg/d2g

         else
c                                 full disordered
            done = .true.

         end if

      end if
c                                 iteration loop
      if (.not.done) then
c                                 increment and check bounds
         call pcheck (x,xmin,xmax,dx,done)
c                                 iteration counter
         itic = 0

         do

            call dgfesi (dg,d2g,y,x,g12,rt)

            dx = -dg/d2g

            call pcheck (x,xmin,xmax,dx,done)

            if (done) then

               exit

            else

               itic = itic + 1
               if (itic.gt.iopt(21)) exit

            end if

         end do

      end if

      gfesi = gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)
c                                 order check, compare to the
c                                 max order g
      g0 = gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)
      if (gfesi.gt.g0) gfesi = g0
c                                 min order g
      g0 = gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)
      if (gfesi.gt.g0) gfesi = g0
c                                 add magnetic component
      gfesi = gfesi + gmag(y)

      end

      double precision function gfes (y,g1,g2)
c-----------------------------------------------------------------------
c gfes returns the Gibbs free energy for Fe-S fluid after
c Saxena & Eriksson 2015.

c coded by ecrg Dec 2017 with cribbing from the Fe-Si models

c    y   - the bulk S mole fraction
c    g1 - free energy of S liq
c    g2 - free energy of Fe liq

c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical done

      integer itic

      double precision g1, g2, y, x, g00, g01, g02, g04, g10, g20, g30,
     *                 rt, xmin, xmax, dg, d2g, dx, gfes0, g0

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      if (y.le.nopt(50).or.y.ge.nopt(56)) then
c                          endmember compositions, no order possible
         gfes =  y*g2 + (1d0 - y)*g1
         return
      end if


      g00 = -104888.1d0 + 3.3884608d-1*t + 9.489d-2*p
     *                                      + 3.4769476d-5*t*p
c                       or   + 1.7687597d-1*p - 8.5431919d-6*t*p in green2.dat
      g01 = -8626.2578d0
      g02 = 72954.295d0 - 26.1780d0*t
      g04 = 25106d0
      g10 = 35043.323d0 - 9.880908d0*t - 5.1303766d-1*p
     *                                      - 2.5038372d-7*t*p
      g20 = -23972.273d0
      g30 = 30436.822d0

      rt  = r*t

c                          max/min concentrations of ordered species.
c                          for y=b/(a+b) and a-b formation limited
c                          by b,
c                          xmax = (2 y Zab Zba)/
c                                   (Zaa Zab - y Zaa Zab - y Zaa Zba + 2 y Zab Zba);
c                          the case below is for ZFeFe = ZSS = 6;  ZFeS = ZSFe = 2
      xmin = nopt(50)

      if (y.lt.0.5d0) then
         xmax = 2.d0*y/(3d0 - 4d0*y) - nopt(50)
      else
         xmax = 2d0*(1d0-y)/(3d0 - 4d0*(1d0-y)) - nopt(50)
      end if

c                                 get 1st and 2nd derivatives
      x = xmax

      call dgfes (dg,d2g,y,x,rt,g00,g01,g02,g04,g10,g20,g30)

      done = .false.
c                                 find starting point for newton-raphson
c                                 search
      if (dg.gt.0d0) then
c                                 max ordered concentration
         dx = -dg/d2g

      else
c                                 most disordered concentration
         x = xmin

         call dgfes (dg,d2g,y,x,rt,g00,g01,g02,g04,g10,g20,g30)

         if (d2g.gt.0d0) then
c                                 sanity check
            dx = -dg/d2g

         else
c                                 full disordered - shouldn't be in here
            done = .true.

         end if

      end if
c                                 iteration loop
      if (.not.done) then
c                                 increment and check bounds
         call pcheck (x,xmin,xmax,dx,done)
c                                 iteration counter
         itic = 0

         do

            call dgfes (dg,d2g,y,x,rt,g00,g01,g02,g04,g10,g20,g30)

            dx = -dg/d2g

            call pcheck (x,xmin,xmax,dx,done)

            if (done) then

               exit

            else

               itic = itic + 1
               if (itic.gt.iopt(21)) exit

            end if

         end do

      end if

      gfes = gfes0 (y,x,g1,g2,rt,g00,g01,g02,g04,g10,g20,g30)
c                                 compare to
c                                 max order g
      g0 = gfes0 (y,xmax,g1,g2,rt,g00,g01,g02,g04,g10,g20,g30)
      if (gfes.gt.g0) gfes = g0
c                                 min order g
      g0 = gfes0 (y,xmin,g1,g2,rt,g00,g01,g02,g04,g10,g20,g30)
      if (gfes.gt.g0) gfes = g0

      end

      subroutine mtrans (gval,vdp,id)
c----------------------------------------------------------------------
c mtrans sorts through and evaluates transition functions
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision gval, dg, vdp, tc, b, pee, gmags

      external gmags

      integer eos
      common/ cst303 /eos(k10)

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
         if (ltyp(id).eq.1) then
c                                 ubc-type transitions
            call lamubc (p,t,dg,lmda(id),lct(id))
            gval = gval + dg

         else if (ltyp(id).eq.2) then
c                                 standard transitions
            call lamhel (p,t,gval,vdp,lmda(id),lct(id))

         else if (ltyp(id).eq.3) then
c                                 supcrt q/coe lambda transition
            call lamqtz (p,t,gval,lmda(id),id)

         else if (ltyp(id).eq.4) then

            if (eos(id).ne.8.and.eos(id).ne.9) then
c                                 putnis landau model as implemented incorrectly
c                                 in hp98 (ds5)
               call lamla0 (dg,vdp,lmda(id))

            else
c                                 putnis landau model as implemented correctly
c                                 in hp11 (ds6)
               call lamla1 (dg,vdp,lmda(id))

            end if

            gval = gval + dg

         else if (ltyp(id).eq.5) then
c                                 holland and powell bragg-williams model
            call lambw (dg,lmda(id))
            gval = gval + dg

         else if (ltyp(id).eq.7) then
c                                 George's Hillert & Jarl magnetic transition model
            if (lct(id).gt.1) write(0,*)'**>1 type = 7 trans.!?'
            tc = therlm(1,1,lmda(id))
            b = therlm(2,1,lmda(id))
            pee = therlm(3,1,lmda(id))
            gval = gval + gmags (tc,b,pee)

         else

            call errdbg ('no such transition model')

         end if

      end

      subroutine savdyn (ids)
c----------------------------------------------------------------------
c subroutine to save exploratory stage dynamic compositions for use
c as static compositions during auto-refine, pa loaded by sollim
c  ids - pointer to solution model
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical rplica, isend

      integer ids

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      external rplica, isend
c----------------------------------------------------------------------
      if (refine.and..not.lopt(55)) return
c                                 currently all calls to savdyn make replicate 
c                                 test with hardwired tolerance nopt(35)
      if (rplica(ids)) return
c                                 deleted degeneracy test.
      if (isend(ids)) return

      tpct = tpct + 1

      if (tpct.gt.m24) call errdbg ('increase m24')
      if (tcct+nstot(ids).gt.m25) call errdbg ('increase m25')
c                                 solution pointer
      dkp(tpct) = ids
c                                 save the composition
      txco(tcct+1:tcct+nstot(ids)) = pa(1:nstot(ids))

      if (lorder(ids))
     *       txco(tcct+nstot(ids)+1:tcct+tstot(ids)) = 
     *       pp(1:lstot(ids))
c                                 save the starting position - 1
      itxp(tpct) = tcct
c                                 increment the counter
      tcct = tcct + tstot(ids)

      end 

      logical function isend (id)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, i, j

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c-----------------------------------------------------------------------
      i = 0

      do j = 1, nstot(id)
         if (dabs(pa(j)).gt.zero) then 
            i = i + 1
            if (i.gt.1) then 
               isend = .false.
               return
            end if
         end if
      end do

      isend = .true.

      end


      logical function rplica (id)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, i, j, tmp, ltot, ntot

      double precision tol, diff

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c-----------------------------------------------------------------------
      ltot = lstot(id)
      ntot = nstot(id)

      tol = nopt(35)

      if (.not.lorder(id).and.ntot.ne.ltot) call errdbg ('oink')
c                                o/d models use the pp array, which is 
c                                not normalized for non-equimolar o/d, do
c                                the normalization here
      if (.not.equimo(id)) then

         diff = 0d0

         do j = 1, ltot
            diff = diff + pp(j)
         end do

         do j = 1, ltot
            pp(j) = pp(j)/diff
         end do

      end if

      do i = stpct, tpct

         if (dkp(i).ne.id) cycle

         diff = 0d0

         if (lorder(id)) then

            tmp = itxp(i) + ntot

            do j = 1, ltot
               diff = diff + dabs(pp(j) - txco(tmp+j))
            end do

         else

            tmp = itxp(i)

            do j = 1, ltot
               diff = diff + dabs(pa(j) - txco(tmp+j))
            end do

         end if

         if (diff.gt.tol) cycle

         rplica = .true.
         return

      end do

      rplica = .false.

      end
