c-----------------------------------------------------------------------
 
c FLIB - fluid phase subroutines common to FRENDLY, VERTEX, COHSRK,
c        RK, and BUILD.

c Unless otherwise noted, the subroutines herein were written by
c J. A. D. Connolly.
 
c Please do not distribute this source.
      
c-----------------------------------------------------------------------
      subroutine cfluid (fo2,fs2)
c-----------------------------------------------------------------------
c subroutine cfluid call fluid equations of state depending on
c the value of ifug. The GCOH eos return ln(fo2) as fo2. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision fo2,fs2,yo2,yc,dg

      double precision f
      common/ cst11 /f(3)

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      if (xco2.gt.1d0) then 
         xco2 = 1d0
      else if (xco2.lt.0d0) then 
         xco2 = 0d0
      end if 

      if (ifug.eq.0) then
         call mrk
      else if (ifug.eq.1) then 
         call hsmrk
      else if (ifug.eq.2) then 
         call qrkmrk
      else if (ifug.eq.5) then 
         call hprk
      else if (ifug.eq.8) then  
         call cohfo2 (fo2)
      else if (ifug.eq.10) then
         call gcohx6 (fo2)
      else if (ifug.eq.12) then 
         call cohsgr (fo2,fs2)
      else if (ifug.eq.13) then 
         call hh2ork (fo2,.false.)
      else if (ifug.eq.14) then 
         call pshp 
      else if (ifug.eq.15) then 
         call hh2ork (fo2,.true.)
      else if (ifug.eq.16) then 
         call homrk (fo2)
      else if (ifug.eq.17) then 
         call hosrk5 (fo2,fs2)
      else if (ifug.eq.19.or.ifug.eq.20) then 
         call xoxsrk (fo2,fs2)
      else if (ifug.eq.24) then 
         call cohngr (fo2)
      else if (ifug.eq.25) then 
         call waddah
      else if (ifug.eq.26) then 
c         call rkboth
         call idsi5 
      else if (ifug.eq.27) then 
c                                 coming in from fluids xco2 is no/(no+nh)
c                                 fs2 = nc/(no+nh+nc), convert to c-o2-h2
c                                 bulk coordinates
         yo2 = xco2 * (1d0 - fs2) / (1d0 + fs2)
         yc  = 2d0 * fs2 / (1d0 + fs2)
         call rkcoh6 (yo2,yc,dg)
      else 
         call error (11,xco2,ifug,'EoS (routine CFLUID)') 
      end if 

      end

      subroutine rfluid (irk)
c---------------------------------------------------------------------
c irk = 1 - write/read prompt for fluid equations of state
c irk = 2 - write fluid equation of state for outtit to unit n3
c irk = 3 - write fluid equation of state to console
c---------------------------------------------------------------------
      implicit none
   
      include 'perplex_parameters.h'

      integer nrk,i,irk,ier

      parameter (nrk=27)
   
      character rkname(0:nrk)*63, y*1

      character vname*8, xname*8
      common / csta2 /xname(k5),vname(l2)

      double precision buf
      common/ cst112 /buf(5)

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer iam
      common/ cst4 /iam

      save rkname 

      data (rkname(i), i = 0, 11)/
c 0
     *'X(CO2) H2O-CO2 Modified Redlich-Kwong (MRK) DeSantis et al 74',
c 1
     *'X(CO2) H2O-CO2 HSMRK Kerrick & Jacobs 81',
c 2
     *'X(CO2) H2O-CO2 MRK hybrid-EoS*',
c 3
     *'Disabled Eos',
c 4
     *'Disabled Eos',
c 5
     *'X(CO2) H2O-CO2 CORK Holland & Powell 91, 98',
c 6-7
     *'Disabled Eos',
     *'Disabled Eos',
c 8
     *'f(O2/CO2) C-buffered COH MRK hybrid-EoS*',
c 9
     *'Disabled Eos',
c 10
     *'X(O) C-buffered COH MRK hybrid-EoS Connolly & Cesare 93*',
c 11
     *'Disabled Eos'/

      data (rkname(i), i = 12, nrk)/
c 12
     *'X(O)-f(S2) C-buffered COHS MRK hybrid-EoS Connolly & Cesare 93*',
c 13
     *'X(H2) H2O-H2 MRK hybrid-EoS*',
c 14
     *'X(CO2) H2O-CO2 Pitzer & Sterner 94; Holland & Powell mixing 03',
c 15
     *'X(H2) H2O-H2 low T MRK hybrid-EoS*',
c 16
     *'X(O) H-O MRK hybrid-EoS*',
c 17
     *'X(O)-f(S2) H-O-S MRK hybrid-EoS*',
c 18
     *'Disabled Eos',
c 19
     *'X(O)-X(S) C-buffered COHS MRK hybrid-EoS Connolly & Cesare 93*',
c 20
     *'X(O)-X(C) COHS MRK hybrid-EoS Connolly & Cesare 93*',
c 21-23
     *'Disabled Eos',
     *'Disabled Eos',
     *'Disabled EoS',
c 24
     *'f(O2/CO2)-N/C C-buffered COHN MRK hybrid-EoS*',
c 25
     *'X(CO2)-X(NaCl) H2O-CO2-NaCl Aranovich et al 10',
     *'X(O) O-Si MRK Connolly 16',
     *'X(O)-X(C) C-O-H MRK hybrid-EoS*'/
c---------------------------------------------------------------------
      if (irk.eq.2) then

         write (n3,1060) rkname(ifug)
         call hybout (ifug,n3)
         return

      else if (irk.eq.3) then

         write (*,1060) rkname(ifug)
         call hybout (ifug,6)
         return
      
      end if     

      elag = 0d0
      ibuf = 1
      dlnfo2 = 0d0

10    write (*,1000)

      do i = 0, nrk
         if (i.eq.4.or.i.eq.6.or.i.eq.9.or.i.eq.18.or.i.eq.21.or.    
     *       i.eq.3.or.i.eq.22.or.i.eq.23.or.i.eq.7.or.i.eq.11) cycle 

         write (*,1070) i,rkname(i)

      end do 
c                                 write hybrid eos blurb
      call hybout (-1,6)

      read (*,*,iostat=ier) ifug
      if (ifug.gt.nrk.or.ifug.eq.4.or.ifug.eq.6.or.ifug.eq.7.or.
     *    ifug.eq.18.or.ifug.eq.21.or.ifug.eq.22.or.ifug.eq.23.or.
     *    ifug.lt.0.or.ifug.eq.3.or.ifug.eq.9) ier = 1

      call rerror (ier,*10)

      if (iam.ne.11.and.ifug.eq.27) then 
         write (*,56)
         goto 10 
      end if 
c                                 warn about variable mass-speciation 
c                                 routines
      if (iam.eq.5.or.iam.eq.11) then 
         if (ifug.eq.5.or.ifug.eq.14) then
            write (*,1191)
         else if (ifug.eq.26.or.ifug.eq.27.or.
     *          ifug.eq.8.or.ifug.eq.10.or.ifug.eq.11.or.
     *          ifug.eq.12.or.ifug.eq.16.or.ifug.eq.17.or.ifug.eq.19.or.
     *          ifug.eq.20.or.ifug.eq.24) then 
            write (*,1190)
         end if 
      end if 

      if (ifug.eq.12.or.ifug.eq.17.or.ifug.eq.20) then
c                                 COHS & HOS equations of state
c                                 get sulfur fugacity constraint:
         vname(3) = 'X(O)'
12       write (*,1090) 
         read (*,*,iostat=ier) ibuf
         if (ibuf.gt.3.or.ibuf.lt.1) ier = 1
         call rerror (ier,*12)

         if (ibuf.eq.2) then 
c                                 if ibuf = 2 dlnfo2 is the fe/s
c                                 of pyrrhotite
13          write (*,1100) 
            read (*,*,iostat=ier) dlnfo2
            call rerror (ier,*13)
         else if (ibuf.eq.3) then
c                                 if ibuf = 2 dlnfo2 is the fs2
14          write (*,1110)
            read (*,*,iostat=ier) dlnfo2
            call rerror (ier,*14)
            dlnfo2 = 2.302585093d0 * dlnfo2

         end if 

      else if (ifug.eq.10.or.ifug.eq.12.or. 
     *         ifug.gt.16.and.ifug.lt.21) then

         if (ifug.ne.18) vname(3) = 'X(O)'

      else if (ifug.eq.13.or.ifug.eq.15) then
 
         vname(3) = 'X(H2)'

      else if (ifug.eq.26) then
 
         vname(3) = 'X(Si)'

      else if (ifug.eq.8.or.ifug.eq.24) then
c                                 chosen COH speciation option
c                                 check that XCO2 isn't a independent
c                                 variable:
         if (iv(1).eq.3.or.iv(2).eq.3) then 
            call warn (172,dlnfo2,ifug,'RFLUID')
            goto 10
         end if
c                                 change default buffer
20       write (*,1020)
         read (*,1030) y
         ibuf = 2
         dlnfo2 = 0d0
         if (y.eq.'y'.or.y.eq.'Y') then 

            write (*,1010) 
            read (*,*,iostat=ier) ibuf
            call rerror (ier,*20)

            if (ibuf.gt.5.or.ibuf.lt.1) then
               call warn (173,dlnfo2,ifug,'RFLUID')
               goto 20
            end if 
c                                 ibuf = 5, define own buffer
            if (ibuf.eq.5) then
45             write (*,1180)  
               read (*,*,iostat=ier) buf
               call rerror (ier,*45)
               goto 30
c                                 ibuf = 3, constant fo2.
            else if (ibuf.eq.3) then
40             write (*,1140) 
               read (*,*,iostat=ier) dlnfo2
               call rerror (ier,*40)
               dlnfo2 = 2.302585093d0 * dlnfo2
               goto 50
            end if 
c                                 ibuf 2 or 1, permit del(fo2)
30          write (*,1040) 
            read (*,1030) y
            if (y.eq.'y'.or.y.eq.'Y') then 
               write (*,1050) 
               read (*,*,iostat=ier) dlnfo2
               call rerror (ier,*30)
               dlnfo2 = 2.302585093d0 * dlnfo2
            end if 
         end if
      end if 
c                                 for graphite EoS's allow
c                                 reduced gph activity:
c                                 this could be done for all
c                                 but here we only allow simple
c                                 EoS's because the X(C),S/C and N/C
c                                 routines use the variable elag to
c                                 store these ratios.
50    if (ifug.ge.8.and.ifug.le.12) then
c                                 if build, then ask whether to 
c                                 set H2 and O2 as the independent
c                                 fugacities (hu = 1), otherwise they 
c                                 are H2O and CO2 (hu = 0).
         hu = 0 

         if (iam.eq.4) then 
            write (*,1170)
            read (*,1030) y 
            if (y.eq.'y'.or.y.eq.'Y') hu = 1
         end if 

         write (*,1120) 
         read (*,1030) y

         if (y.eq.'y'.or.y.eq.'Y') then
            write (*,1130) 
            read (*,*,iostat=ier) elag
            call rerror (ier,*50)
            elag = dlog (elag)
         end if

      else if (ifug.eq.19) then 
c                                 for XO-XS EoS get S/C
         write (*,1150) 
         read (*,*,iostat=ier) elag
         call rerror (ier,*50)

      else if (ifug.eq.20) then
c                                 for XO-XC EoS get X(C)
         write (*,1160) 
         read (*,*,iostat=ier) elag
         call rerror (ier,*50)
    
      else if (ifug.eq.25) then 
c                                 special conditions for H2O-CO2-NaCl EoS
         write (*,1200)
         read (*,*,iostat=ier) ibuf
         call rerror (ier,*50)
         if (ibuf.lt.1.or.ibuf.gt.2) then
            write (*,*) 'invalid choice, try again ...'
            goto 50 
         end if 
c                                get the salt content (elag):
         if (ibuf.eq.1) then 
            write (*,1210) 'weight'
         else if (ibuf.eq.2) then 
            write (*,1210) 'molar '
         end if 

         read (*,*,iostat=ier) elag
         call rerror (ier,*50)
         if (elag.gt.1d0.or.elag.lt.0d0) then
            write (*,*) 'salt fraction must be > 0 and < 1, try again'
            goto 50 
         end if 

         vname(3) = 'Y(CO2)*'

      end if 

1000  format (/,'Select fluid equation of state:',/)
1010  format (/,'Select buffer: ',//,
     *          ' 1 - aQFM, 298-1200K',/,
     *          ' 2 - Maximum H2O content, 523-1273K, .5-30kbar',/,
     *          ' 3 - user specified f(O2)',/,
     *          ' 4 - aQ-Ru-Cc-Tn-Gph',/,
     *          ' 5 - ln(f(O2))= a + (b + c*p)/t + d/t**2 + e/t**3 ',/)
1020  format (/,'Modify default buffer (max H2O) (Y/N)? ')
1030  format (a)
1040  format (/,'Modify calculated fO2 by a constant (Y/N)?',/)
1050  format (/,'Enter constant in units of log10(fO2):',/)
1060  format (/,'Fluid equation of state: ',a)
1070  format (2x,i2,' - ',a)
1090  format (/,'Choose a sulfur buffer:',//
     *         ,'  1 - Pyrite + Pyrrhotite',
     *        /,'  2 - Pyrrhotite',/,'  3 - f(S2)',/)
1100  format (/,'Enter atomic Fe/S of pyrrhotite:',/)
1110  format (/,'Enter log10[f(S2)]:',/)
1120  format (/,'Reduce graphite activity (Y/N)?',/)
1130  format (/,'Enter activity of graphite:',/)
1140  format (/,'Enter log10[f(O2)]:',/)
1150  format (/,'Enter X(S), i.e., {n(S)/[n(S) + n(C)]}:',/)
1160  format (/,'Enter X(C), i.e., {n(C)/[n(S)+n(C)+n(O)+n(H)]}:',/)
1170  format (/,'Compute f(H2) & f(O2) as the dependent fugacities',
     *        /,'(do not unless you project through carbon) (Y/N)?',/)
1180  format ('Enter a-e :',/)
1191  format (/,'**warning ver119** volume for this EoS',
     *       ' computed by FRENDLY or FLUIDS',/,'is estimated by',
     *       ' finite difference over a 1 bar interval',/)
1190  format (/,'**warning ver119** partial molar volumes for this EoS',
     *       ' are estimated by finite',/,
     *       'difference under the assumption that fluid',
     *       ' speciation does not change over',/,'the 1 bar',
     *       ' finite difference interval. The estimates'
     *       ' should be considered',/,'with caution.',/)
1200  format (/,'For this EoS Y(CO2)* is defined as:',/,
     *          '  Y(CO2)* = n(CO2)/[n(H2O)+n(CO2)]',/,
     *          'i.e., Y(CO2)* may vary from 0 -> 1 ',
     *          'regardless of salt content',//,
     *          'Choose how salt content is to be specified:',/,
     *          ' 1 - weight fraction',/,
     *          ' 2 - mole fraction',/)
1210  format (/,'Enter ',a,' salt fraction (0->1) in the fluid:',/)
56    format (/,'**warning ver056** C-O-H hybrid cannot be specified ',
     *          'for saturated components',/,'or a saturated phase. To',
     *          'use this EoS see:',/,
     *          'perplex.ethz.ch/perplex/faq/calculations_with_un',
     *          'buffered_COH_fluids.txt',/)

      end 

      subroutine hybout (icheck,nout)
c----------------------------------------------------------------------
c write hybrid fluid equations of state for outtit to unit nout
c if icheck is >= 0 it identifies an eos and output is only
c written if icheck corresponds to a hybrid EoS
c----------------------------------------------------------------------
      implicit none
   
      include 'perplex_parameters.h'

      integer i, j, nhyb, icheck, nout

      parameter (nhyb=7)
   
      character hyname(0:nhyb)*32

      character specie*4
      integer ins, isp
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      save hyname 

      data (hyname(i), i = 0, 7)/
     *  'MRK DeSantis et al 1974',
     *  'HSMRK Kerrick & Jacobs 1981',
     *  'CORK Holland & Powell 1998',
     *  'BRMRK Bottinga & Richet 1981',
     *  'PSEoS Pitzer & Sterner 1994',
     *  'Haar et al 1982',
     *  'Zhang & Duan 2005',
     *  'Zhang & Duan 2009'/
c----------------------------------------------------------------------

      if (icheck.lt.0) then 

         write (nout,1000)

         do j = 1, 3
            if (j.eq.1) then 
               write (nout,1020) specie(j),hyname(iopt(25))
            else if (j.eq.2) then 
               write (nout,1020) specie(j),hyname(iopt(26))
            else if (j.eq.3) then 
               write (nout,1020) specie(j+1),hyname(iopt(27))
            end if 
         end do 

      else if (icheck.eq.2.or.icheck.eq.8.or.icheck.eq.10.or.
     *         icheck.eq.12.or.icheck.eq.13.or.icheck.eq.27.or.
     *        (icheck.ge.15.and.icheck.le.20)) then 

         write (nout,1010) 

         do i = 1, isp
            j = ins(i)
            if (j.eq.1) then 
               write (nout,1020) specie(j),hyname(iopt(25))
            else if (j.eq.2) then 
               write (nout,1020) specie(j),hyname(iopt(26))
            else if (j.eq.4) then 
               write (nout,1020) specie(j),hyname(iopt(27))
            end if 

         end do

      end if

1000  format (/,'*Hybrid EoS use the following pure species EoS, ',
     *      'to change these associations',/,'modify the hybrid_EoS ',
     *      'keywords in the perplex_option file:',/)
1010  format (/,'*This hybrid EoS uses the following pure species EoS,',
     *      ' to change these associations',/,'modify the hybrid_EoS ',
     *      'keywords in the perplex_option file:',/)
1020  format (7x,a,' - ',a)

      end 

      subroutine setins (jfug)
c-----------------------------------------------------------------------
c set species for internal eos choices and, if appropriate, variable name.

c    jfug - eos pointer
c    isp  - species counter
c    ins  - pointer to species

c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jfug, i

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)
c----------------------------------------------------------------------- 
      if (jfug.le.5.or.jfug.eq.14.or.jfug.eq.25) then 
c                                 xco2 EoS's
            vname(3) = 'X(CO2)  '

            isp = 2

            ins(1) = 1
            ins(2) = 2

      else if (jfug.ge.8.and.jfug.le.12.or.jfug.eq.19.or.
     *         jfug.eq.20.or.jfug.eq.24.or.jfug.eq.27) then
c                                 standard COHS species  
         if (jfug.eq.8.or.jfug.eq.24) then
            vname(3) = 'log(fO2)'
         else 
            vname(3) = 'X(O)    '
         end if 

         isp = 5

         do i = 1, 6
            ins(i) = i
         end do
 
         if (jfug.eq.10) then 

            isp = 6 
            ins(6) = 16

         else if (jfug.eq.19.or.jfug.eq.20) then

            isp = 8
            ins(7) = 8
            ins(8) = 9
  
         else if (jfug.ge.12.and.jfug.le.18) then

            isp = 9
            ins(7) = 7
            ins(8) = 8
            ins(9) = 9

         else if (jfug.eq.24) then

            isp = 7
            ins(6) = 10
            ins(7) = 11

         else if (jfug.eq.27) then
c                                 C-O-H free
            vname(4) = 'Y(C)    '

            isp = 6
            ins(6) = 7

         end if

      else if (jfug.eq.13.or.jfug.eq.15) then 

         vname(3) = 'X(H2)   '

         isp = 2
         ins(1) = 1
         ins(2) = 5

      else if (jfug.eq.16) then

         vname(3) = 'X(O)    '

         isp = 3
         ins(1) = 1
         ins(2) = 5
         ins(3) = 7

      else if (jfug.eq.17) then

         vname(3) = 'X(O)    '

         isp = 5
         ins(1) = 1
         ins(2) = 5
         ins(3) = 6
         ins(4) = 7
         ins(5) = 8

      else if (jfug.eq.26) then
c                                 silica vapor 
         vname(3) = 'X(Si)   '

         isp = 5
         ins(1) = 14
         ins(2) = 13
         ins(3) = 12
         ins(4) = 7
         ins(5) = 15

      else 

         call error (74,1d0,jfug,vname(3))

      end if 

      end 

      subroutine cohsgr (fo2,fs2)
c----------------------------------------------------------------------
c program to calculate graphite saturated C-H-O-S speciation as 
c a function of XO using an MRK/HSMRK hybrid.

c Species are CO2 CH4 CO H2 H2O H2S O2 SO2 COS. The latter 3 species
c are only significant for reduced graphite activities. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ins(9),jns(3),nit,ier

      logical lmt

      double precision fo2,fs2,oh2o,c1,c2,c3,c4,c5,c6,c7,
     *                 ek1,ek2,ek3,ek4,ek5,ek6,ek7

      integer ibuf,hu,hv,hw,hx
      double precision rat,elag,gz,gy,gx
      common/ cst100 /rat,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      double precision vol
      common/ cst26 /vol

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision eqk
      common / csteqk /eqk(nsp)

      save ins, jns
      data ins, jns/1,2,3,4,5,6,7,8,9,1,2,4/
c----------------------------------------------------------------------
      nit = 0
      oh2o = 2d0
c                                check for in bounds composition
      call xcheck (xo,lmt)
c                                fs2 = 1/2 ln (fs2)
      call setfs2 (fs2)
c                                compute equilibrium constants in csteqk
      call seteqk (ins,9,elag)
c                                compute pure mrk fluid properties
      call mrkpur (ins,9)
c                                compute hybrid pure fluid props
      call hybeos (jns,3)

      call zeroys

      if (lmt) return

      y(5) = 0.00001d0

      c3 = dexp (eqk(4)) * p
      c1 = dexp (eqk(2) - 2d0*eqk(3)) * p 
      c2 = dexp (eqk(1) - eqk(3)) * p
      c4 = dexp (eqk(6) + fs2)
      c5 = dexp (eqk(9) + fs2)
      c6 = dexp (eqk(8) - 2d0*eqk(3) + fs2) * p
      c7 = dexp (-2d0*eqk(3)) * p
c                                outer iteration loop: 
10    ek1 = c1 * g(3)**2/g(2) 
      ek2 = c2 * g(3) * g(5)/g(1)
      ek3 = c3 * g(5)**2/g(4) 
      ek4 = c4 * g(5)/g(6)
      ek5 = c5 * g(3)/g(9)
      ek6 = c6 * g(3)**2/g(8)
      ek7 = c7 * g(3)**2/g(7)
c                                 solve for xh2, xco
      call evlxh1 (ek1,ek2,ek3,ek4,ek5,ek6,ek7,xo,y(5),y(3),ier)

      if (ier.ne.0) call warn (501,xo,ier,'COHSGR')

      y(4) = ek3 * y(5)**2 
      y(1) = ek2 * y(5) * y(3)
      y(2) = ek1 * y(3)**2
      y(6) = ek4 * y(5)
      y(8) = ek6 * y(3)**2
      y(9) = ek5 * y(3)
      y(7) = ek7 * y(3)**2

      nit = nit + 1

      if (nit.gt.iopt(21)) call warn (502,xo,ier,'COHSGR')

      if (dabs(y(1)-oh2o).lt.nopt(50)) goto 90

      oh2o = y(1)

      call mrkhyb (ins, jns, 9, 3, 1)

      goto 10 

90    do i = 1, 3 
         vol = vol + y(jns(i))*vh(jns(i))
      end do 

      goto (91), hu

      fh2o = dlog(g(1)*p*y(1))
      fco2 = dlog(g(2)*p*y(2))
      fo2 = 2d0 * (dlog(g(3)*p*y(3)) - eqk(3))

      goto 99

91    fh2o = dlog(g(5)*p*y(5))
      fco2 = 2d0 * (dlog(g(3)*p*y(3)) - eqk(3))

99    end

      subroutine evlxh1 (ek1,ek2,ek3,ek4,ek5,ek6,ek7,xo,xh2,yco,ier)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier,it
      double precision ek1,ek2,ek3,ek4,ek5,ek6,ek7,xo,xh2,yco,f0,e1,e2,
     *                 e3,e4,e5,e6,e7,e8,e9,e0,r0,t2,t10,t11,t15,c1,g,dg
c----------------------------------------------------------------------
      it = 0 
      ier = 0 

      f0 = 2d0*(ek7 + ek6 + ek1)
      e0 = 2d0*xo
      e1 = 1d0/f0
      e2 = 1d0 + ek5**2 + 2d0*(ek5 + f0)
      e3 = 2d0*ek2*(1d0 + ek5) - 2d0*f0*(ek4 + 1d0)
      e4 = ek2**2 - 2d0*ek3*f0
      e5 = e0 + e0*ek4 
      e6 = 4d0*xo*ek3
      e7 = xo - ek5 - 1d0 + xo*ek5
      e8 = f0 * (xo - 1d0)
      e9 = ek2*(3d0* xo - 1d0)

10    r0 = xh2
      t2 = xh2**2
      t10 = e2 + e3*xh2 + e4*t2

      if (t10.lt.0d0) then
c                                 if t10 < 0, bad guess
c                                 for xh2, find roots:
         c1 = dsqrt (e3**2 - 4d0*e4*e2)
         xh2 = 0.9d0*(-c1 - e3/2d0/ e4)
         r0 = xh2 
         t2 = xh2**2
         t10 = e2 + e3*xh2 + e4*t2

      end if 

      t10 = dsqrt (t10)

      t11 = t10 - 1d0 - ek2*xh2 - ek5
      yco = e1*t11
      g = e5*xh2 + e6*t2 + (e7 + e8*yco + e9*xh2)*yco
      t15 = (e3+2d0*e4*xh2)/2d0/t10 - ek2
      dg = e5 + 2d0*e6*xh2 + e1*t15*(e9*xh2 + e7)
     *        + t11*(2d0*e8*e1**2*t15 + e9*e1)

      xh2 = r0 - g/dg
      if (xh2.lt.0d0) xh2 = r0/2d0
c                                 converged:
      if (dabs(xh2-r0).lt.nopt(50)) goto 999

      it = it + 1

      if (it.gt.iopt(21)) then 

         ier = 2

         goto 99

      end if 

      goto 10 

999   yco = e1*(dsqrt(e2 + (e3 + e4*xh2)*xh2) - 1d0 - ek2*xh2 - ek5)

99    end


      subroutine evlxh2 (c1,c2,c3,c5,c6,xo,xs,xh2,yco,xh2o,ier)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier,jt,it

      double precision c1,c2,c3,c5,c6,c7,xo,xs,xh2,yco,xh2o,d1,d2,d3,
     *                 d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,
     *                 d17,d18,d19,d20,d21,d22,d23,d24,d25,f10,r1,g1,
     *                 g2,r0,e1,e2,e3,e4,t14,t31,t37,t39,t43,t45,g,
     *                 t27,t49,t55,t61,t65,t67,t71,t73,t74,t80,t88,t98,
     *                 dg,t10,t11,t25,f2,f3,f4,f7,f8,f9,f11,f12,t4,f13,
     *                 f1,t32,f5,f6,f,t13,t38,t47,df
c----------------------------------------------------------------------
      jt = 0 
      ier = 0 
      c7 = 2d0*c5
      d2 = c5 - xs * c5
      d3 = c3 - xs * c3
      d4 = c1 - 2d0*xs*c1 - xs
      d6 = 2d0 * c2
      d7 = xs * c2
      d8 = xs * c6
      d9 = 4d0 * c6
      d10 = 2d0 * c3
      d11 = 12d0 * c6
      d12 = 5d0 * d8
      d25 = d10 + 2d0 
      f10 = 3d0*d2

100   r1 = xh2o

      it = 0 
      ier = 0 

      d1 = xh2o**2
      d5 = c7 * d1
      d13 = d7 * d1
      d14 = 2d0*d4*xh2o
      d15 = d2*d1
      d16 = d3*xh2o
      d17 = d4*xh2o  
      d18 = -3d0* d16
      d19 = -4d0*c5*d1
      d20 = 3d0*xh2o
      d21 = d9/xh2o
      d22 = d11*d3
      d23 = d11/xh2o
      d24 = d15*xh2o

      g1 = -6d0*c2*d1*d3
      g2 = d6*xh2o
c                                 solve g for xh2, guessed xh2o
10    r0 = xh2
c                                 evaluate g:
      e1 = xh2**2
      e2 = e1*xh2
      e3 = e2*xh2
      e4 = e3*xh2
      t14 = d24 + d16*e2
      t31 = d17*e1 - xh2*d13 - d8*e4
      t39 = -t14/t31
      t37 = g2*t39/xh2
      t43 = d5/e1
      t45 = c1*t39

      g = (t37 + t39 + t43 + xh2o + t45)/
     *    (t37 + t39 + t43 + d20  + t45 + 2d0*xh2  
     *         - d21*t14/t31*e2 + d10*xh2) - xo
c                                 evaluate dg:
      t27 = g1*xh2/t31
      t37 = t31**2
      t49 = d14*xh2 - d13 - d12*e3
      t55 = g2*t14/t37/xh2*t49
      t61 = -g2*t39/e1
      t65 = d18*e1/t31
      t67 = t14/t37*t49
      t71 = d19/e2
      t73 = c1*t65
      t74 = c1*t67
      t80 = -g2*t14/t31/xh2
      t88 = c1*t39
      t98 = t80 + t39 + t43 + d20 + t88 + 2d0*xh2 
     *          + d21*t39*e2 + d10*xh2

      dg = (t27 + t55 + t61 + t65 + t67 + t71 + t73 + t74)/t98
     *     - (t80 + t39 + t43 + xh2o + t88)/t98/t98
     *      *(t27 + t55 + t61 + t65 + t67 + t71 + t73 + t74 
     *        - d22*e4/t31 + d21*e2*t67 + d23*t39*e1 + d25)

      xh2 = r0 - g/dg
      if (xh2.lt.0d0) xh2 = r0/2d0
c                                 converged:
      if (dabs(xh2-r0).lt.nopt(50)) goto 40

      it = it + 1

      if (it.gt.iopt(21)) then 

         ier = 2

         goto 40

      end if 

      goto 10 

40    it = 0 
c                                 use xh2 to refine guesssed xh2o
      t10 = xh2**2
      t11 = xh2**3
      t25 = t10*t11
      f2 = 2d0*d7*xh2
      f3 = c2/xh2
      f4 = xh2 - 1d0 + c3*xh2
      f7 = c2*xh2
      f8 = d3*t11
      f9 = d7*xh2
      f11 = c5/t10
      f12 = d4*t10

30    r0 = xh2o
c                                 evaluate f:
      t4 = xh2o**2
      f13 = d2*t4
      f1 = f12*xh2o - t4*f9 - d8*t25
      t14 = f13*xh2o + xh2o*f8
      t32 = t14/f1
      f5 = c6*t11
      f6 = f5/xh2o

      f = -t32 - f3*xh2o*t32 - f5*t32/xh2o - c1*t32 
     *    + t4*f11 + xh2o + f4
c                                 evaluate df:
      t13 = f10*t4 + f8
      t31 = t13/f1
      t37 = d2*t4*xh2o + xh2o*f8
      t38 = f1**2
      t45 = d4*t10 - f2*xh2o
      t47 = t37/t38*t45
      t49 = -f3*f1

      df = t47 - t31 + t37*t49 + xh2o*t13*t49 
     *     + f7*xh2o*t47 - f6*t31
     *     + f6*t47 + f5*t37/f1/t4
     *     - c1*t31 + c1*t47 + c7*xh2o/t10 + 1d0

      xh2o = r0 - f/df

      if (xh2o.lt.0d0) then
         xh2o = r0/2d0
      else if (xh2o.ge.1d0) then
         xh2o = r0 + (1d0 - r0)/2d0
      end if         
c                                 converged:
      if (dabs(xh2o-r0).lt.nopt(50)) goto 99

      it = it + 1

      if (it.gt.iopt(21)) then 

         ier = 2

         goto 99

      end if 

      goto 30 

99    t4 = xh2o**2

      yco = -(d2*t4*xh2o + d3*xh2o*t11)/
     *       (d4*t10*xh2o - d7*xh2*t4 - d8*t25)
c                                 is new xh2o same as guess?:
      if (dabs(xh2o-r1).lt.nopt(50)) goto 9999

      jt = jt + 1

      if (jt.gt.100) then 

         ier = 2

         goto 9999

      end if 

      goto 100 

9999  end

      subroutine evlxh3 (c1,c2,c3,c5,c6,xo,xc,xh2,yco,xh2o,ier)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier,jt,it

      double precision c1,c2,c3,c5,c6,c7,xo,xh2,yco,xh2o,r1,g1,
     *                 g2,r0,t39,t43,t45,t27,t49,t61,t65,t71,
     *                 t10,t11,f2,f3,f4,f7,f8,t4,
     *                 f1,t32,f5,f6,f,t47,df,g0,g3,g4,g5,g6,g7,
     *                 g8,g9,g10,g11,g12,g13,c11,c15,c16,c17,c18,c19,
     *                 c20,c25,c26,c27,t1,t2,c8,c9,c10,c12,c14,c21,c22,
     *                 c23,t5,t21,c13,t16,t17,t19,t33,c,t12,t42,t101,
     *                 xc,t46,t58,t60,t62,t79,t81,t87,t91,t93,t95,t99,
     *                 t107,t111,t115,dc,t3,t59,t63,t76,t7
c----------------------------------------------------------------------

      jt = 0 
      ier = 0 
      g0 = 3d0*c5
      g1 = 2d0*c5
      g2 = xo*g1
      g3 = 3d0*xo
      g4 = 2d0*xo
      g5 = g4*c3
      g6 = 2d0*c2
      g7 = g4*c2
      g8 = xo*c1
      g9 = 4d0*xo*c6
      g10 = 6d0*xo
      g11 = 3d0*g4
      g12 = 3d0*g5
      g13 = 6d0*c5
      c11 = 3d0*c3
      c15 = c11 + 2d0
      c16 = 2d0*g3
      c17 = 3d0*g2
      c18 = 3d0*g1
      c19 = 2d0*g6
      c20 = 2d0*g7
      c25 = 2d0*c1
      c26 = 2d0*g8
      c27 = 5d0*g9

100   r1 = xh2o
      t1 = xh2o**2
      t2 = t1*xh2o
 
      c7 = g1*t2
      c8 = g2*t2
      c9 = g0*t1
      c10 = 3d0*xh2o
      c12 = c2*xh2o
      c14 = g13*t1
      c21 = c6/xh2o
      c22 = g7*t1
      c23 = g6*t1

      it = 0 
      ier = 0 
c                                 solve c for xh2, guessed xh2o
10    r0 = xh2
c                                 evaluate c:
      t4 = xh2**2
      t5 = t4*t1
      t21 = t4**2
      t10 = t4*xh2
      c13 = t10/xh2o
      t11 = t10*xh2o
      t16 = c7 + t5 - c8 - g3*t5 - g4*t11 - g5*t11
      t17 = t1*xh2
      t19 = t4*xh2o
      t32 = (g6-g7)*t17+t19+(c1-xo-g8)*t19-g9*t21*xh2
      t33 = t16/t32
      t39 = c12*t33/xh2
      t45 = c21*t33*t10
      t47 = c1*t33

      c = (-t33 - t39 - t45 - t47)/
     *    (-3d0*t39 - 2d0*t33 + c9/t4 + c10 - 3d0*t47 + 2d0*xh2 
     *     -5d0*t45 + c11*xh2) - xc
c                                 evaluate dc:
      t12 = 2d0*t17 - g10*t17 - g11*t19 - g12*t19
      t27 = t12/t32
      t42 = c7 + t5 - c8 -g3*t5 - g4*t11 - g5*t11
      t43 = t32**2
      t101 = t42/t32
      t46 = xh2*xh2o
      t58 = c23 + 2d0*t46 + c25*t46  -c22 - g4*t46 - c26*t46 - c27*t21
      t60 = t42/t43*t58
      t62 = t32*xh2
      t65 = c12*t12/t62
      t71 = c12/xh2*t60
      t76 = c12*t101/t4
      t79 = c6*c13/t32
      t81 = t12*t79
      t87 = c6*c13*t60
      t91 = c21*t101*t4
      t93 = c1*t27
      t95 = c1*t60
      t99 = c12*t42/t62
      t107 = c1*t101
      t111 = t42*t79
      t115 = -3d0*(t99 + t107) + 2d0*(xh2 - t101) + c9/t4 + c10  
     *       -5d0*t111 + c11*xh2

      dc = (-t27+t60-t65+t71+t76-t81+t87-3*t91-t93+t95)
     *     /t115 - (-t101 - t99 - t111 - t107)/t115**2        
     *     *(+ 3d0*(t71 - t65  + t76 - t93 + t95) + 2d0*(t60 - t27)
     *       -c14/t10  + 5d0*(t87- t81) -15d0*t91 + c15)

      xh2 = r0 - c/dc

      if (xh2.lt.0d0) xh2 = r0/2d0
c                                 converged:
      if (dabs(xh2-r0).lt.nopt(50)) goto 40

      it = it + 1

      if (it.gt.iopt(21)) then 

         ier = 2

         goto 40

      end if 

      goto 10 

40    it = 0 
c                                 use xh2 to refine guesssed xh2o
      t4 = xh2**2 
      t10 = t4*xh2
      t27 = t4**2
      f1 = c5/t4
      f2 = c3*xh2 + xh2 -1d0
      f3 = g9*xh2
      f4 = xo*t4
      f5 = c2/xh2
      f6 = + t4 + c1*t4 - f4 - g8*t4
      f7 = g1/t4
      f8 = c6*t10

30    r0 = xh2o
c                                 evaluate f:
      t1 = xh2o**2
      t2 = t1*xh2o
      t17 = t1*xh2
      t5 = t4*t1
      t11 = t10*xh2o
      t19 = t4*xh2o

      t7 = g1*t2 + t5 - g2*t2 - g3*t5 - g4*t11 - g5*t11
      t3 = g6*t17 + t19 + c1*t19 - g7*t17 - xo*t19 - g8*t19 - t27*f3
      t33 = t7/t3

      f = -t33 - f5*xh2o*t33 - f8*t33/xh2o - c1*t33
     *    + f1*t1 + xh2o + f2
c                                 evaluate df:
      t16 = c18*t1 + 2d0*t19- c17*t1- c16*t19 - g4*t10 - g5*t10
      t32 = t16/t3
      t49 = xh2*xh2o
      t59 = c19*t49 - c20*t49 - f6 
      t61 = t7/t3**2*t59
      t63 = c2*t3/xh2

      df = -t32 + t61 - t7*t63 - xh2o*t16*t63 + f5*xh2o*t61
     *     - f8*t32/xh2o + f8/xh2o*t61 + f8*t7/t3/t1 - c1*t32 
     *     + c1*t61 + f7*xh2o + 1d0

      xh2o = r0 - f/df

      if (xh2o.lt.0d0) then
         xh2o = r0/2d0
      else if (xh2o.ge.1d0) then
         xh2o = r0 + (1d0 - r0)/2d0
      end if         
c                                 converged:
      if (dabs(xh2o-r0).lt.nopt(50)) goto 99

      it = it + 1

      if (it.gt.iopt(21)) then 

         ier = 2

         goto 99

      end if 

      goto 30 

99    yco = -((g1-g2)*xh2o**2+((1d0-g3)*xh2o-(g4+g5)*xh2)*xh2**2)*xh2o/
     *       ((((g6-g7)*xh2o+(1d0+c1-xo-g8)*xh2)*xh2o-g9*xh2**4)*xh2)
c                                 is new xh2o same as guess?:
      if (dabs(xh2o-r1).lt.nopt(50)) goto 9999

      jt = jt + 1

      if (it.gt.iopt(21)) then 

         ier = 2

         goto 9999

      end if 

      goto 100 

9999  end



      subroutine evalg (k1,k2,k3,xt,xh,g,dg,sign)
c----------------------------------------------------------------------
      implicit none

      double precision k1,k2,k3,xt,xh,g,dg,sign,u1,k2x,xk2x,t8,t9,t15,
     *                 t20,t21,t24,t32,t34,t36,t40,t41,t43,t53,t54,t66,
     *                 t68,k22
c----------------------------------------------------------------------

      u1 = xt*k1
      k2x = k2*xh
      xk2x = xt*k2x
      k22 = k2**2
      t8 = xh**2
      t9 = k22*t8
      t15 = xt**2
      t20 = k3*t8
      t21 = k1*t20
      t24 = k1*xh

      t32 = dsqrt  (4d0*(t9 - xk2x) 
     *        + t15*(6d0*k2x + 9d0*t9 + 1d0 - 32d0*t21 - 16d0*t24) 
     *        + xt*(32d0*t21 + 16d0*t24 - 12d0*t9))

      t34 = - 2d0*k2x + 3d0*xk2x + xt - sign * t32
      t36 = t34/u1
      t40 = xt*k2
      t41 = t34**2
      t43 = k22*xh
      t53 = k3*xh
      t54 = k1*t53
      t66 = - 2d0*k2 + 3d0*t40 - (4d0*t43 
     *       + xt*(32d0*t54-12d0*t43) - 2d0*t40  
     *       + t15*(9d0*t43 + 3d0*k2 - 32d0*t54 - 8d0*k1)  
     *       + 8d0*u1 )/t32

      t68 = t66/u1

      g = -k2x*t36/4d0 + t41/k1/t15/16d0 + t20 - t36/4d0 + xh - 1d0 

      dg = (-k2/u1*t34 - k2x*t68 + t34/k1/t15*t66/2d0)/4d0
     *      + 2d0*t53 - t68/4d0 + 1d0

      end

      subroutine fo2buf (fo2)
c----------------------------------------------------------------------
c this routine returns the ln(fO2) of buffers as a function of P-T.
c corrections for graphite activity and delta fo2 are made from 
c cst100.

c ibuf - 1 - fit to aQFM from Holland and Powell in the range
c            298-1200 K, probably not valid above aQ/bQ trans.
c            should be refit to revised hp data.

c ibuf - 2 - fit of fo2 to obtain maximum H2O content.

c ibuf - 3 - fo2 input by user

c ibuf - 4 - fit of rutile-titantite-a-quartz-calcite-graphite 
c            holland & powell 1993. 
c            expansivity and compressibility ignored.

c ibuf - 5 - user buffer function
c----------------------------------------------------------------------
      implicit none

      double precision t2,t3,p2,fo2,lp,lt

      double precision p,t,xo,uc,u2,tr,pr,r,ps
      common / cst5  /p,t,xo,uc,u2,tr,pr,r,ps

      double precision ab,bb,cb,db,eb
      common/ cst112 /ab,bb,cb,db,eb

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      double precision ao,bo,co,do,eo,fo,go,ho,io,jo,ko,lo,mo,no,oo,po,
     *                 qo,ro,so,to

      save ao,bo,co,do,eo,fo,go,ho,io,jo,ko,lo,mo,no,oo,po,qo,ro,so,to

      data to,ao,bo,co,do,eo,fo,go,ho,io,jo,ko,lo,mo,no,oo,po,qo,ro,so/
     *-804.2316d0  , -.1652445d0  , -.5376252d-02, -4037433d0    ,
     *-.2091203d-06, -.4638105d-08, 0.3753368d-04, -.3853404d-02,
     *-.5442896d-08, 0.6484263d-13, -121.6754d0  ,  127.5998d0  ,
     *-.1486220d0  , -164866.6d0  , -.1863209d-05, 0.9622612d0  ,
     *  2.097447d0 , -.9838123d-03, 0.3077560d-02, 0.7829503d-03/

      t2 = t*t
      t3 = t2 * t
      p2 = p*p

      if (ibuf.eq.1) then 

         fo2 = 13.5029012d0 + (-46704.69695d0 +.2190281453d0*p)/t
     *         -6145687.892d0/t2 + 754294046.5d0/t3

      else if (ibuf.eq.2) then

         lp = dlog(p)
         lt = dlog(t)

         fo2  = to + t*(ao + do*p + t*(fo + ho*t) + (po+qo*t)/p +
     *               ro*lp) 
     *          + p*(bo + p*(eo + io*p) + so*lt)
     *          + p/t*(jo/t + no*p + oo)
     *          + ko*lt + lo*lp + co/t2 + go*dsqrt(p*t) + mo/p2

      else if (ibuf.eq.3) then

         fo2 = dlnfo2
         return

      else if (ibuf.eq.4) then

         fo2 = 16.8582d0 + (0.2131248d0*p - 53946.36d0)/t 
     *                    - 767509.6d0/t2 + 0.9371923d0/t3

      else if (ibuf.eq.5) then
 
         fo2 = ab + (bb + cb*p)/t + db/t2 + eb/t3

      else

         call error (28,r,ibuf,'FO2BUF') 
       
      end if 

      fo2 = fo2 + dlnfo2
  
      end

      subroutine setfs2 (fs2)
c-----------------------------------------------------------------------
c set fs2 according to buffer choice (ibuf)
c-----------------------------------------------------------------------
      implicit none

      double precision fs2,xf

      integer ibuf,hu,hv,hw,hx
      double precision rat,elag,gz,gy,gx
      common/ cst100 /rat,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps
c                                get sulfur fugacity according to
c                                the value of ibuf:
      if (ibuf.eq.1) then 
c                                get po-py 1/2 ln sulfur fugacity:
c                                from simon poulson:
         fs2 = .005388049d0*t + 10.24535d0 - 15035.91d0/t 
     *                        + 0.03453878d0/t*p

      else if (ibuf.eq.2) then
c                                get suflur fugacity from fe/s ratio
c                                of po, expression derived from 
c                                toulmin & barton 1964, with p
c                                correction from Craig & Scott 1982.
c                                (this only takes into account V(S2),
c                                hope it's right.
         xf = rat/( rat + 1d0)
         fs2 = 197.6309d0*xf + 45.2458d0*dsqrt(1d0- 1.9962d0*xf) 
     *         - 94.33691d0
     *         + (80624.79d0 + 0.2273782d0*p - 197630.9d0*xf)/t
      
      else 

         fs2 = rat/2d0

      end if

      end 

      subroutine xoxsrk (fo2,fs2)
c----------------------------------------------------------------------
c program to calculate C-H-O-S speciation as a function of XO and XS
c or XC using an MRK/HSMRK hybrid. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(8),jns(3),nit,ier,i

      logical lmt 

      double precision fo2,fs2,oh2o,c1,c2,c3,c5,c6,
     *                 ek1,ek2,ek3,ek5,ek6

      integer ibuf,hu,hv,hw,hx
      double precision rat,xs,ag,gy,gx
      common/ cst100 /rat,xs,ag,gy,gx,ibuf,hu,hv,hw,hx

      double precision vol
      common/ cst26 /vol

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      save ins, jns
      data ins, jns/ 1,2,3,4,5,6,8,9,1,2,4/
c----------------------------------------------------------------------
      nit = 0
      oh2o = 2d0
c                                this fs2 = 1/2 ln (fs2),
      call setfs2 (fs2)
c                                 check for in bounds composition
      call xcheck (xo,lmt)
c                                 compute equilibrium constants in csteqk

c                                 who knows how this works, i guess if 
c                                 ifug = 19 graphite saturated and if 
c                                 ifug = 20 it isn't (but then what happens
c                                 to sulfur?)
      call seteqk (ins,8,0d0)
c                                 compute pure mrk fluid properties
      call mrkpur (ins,8)
c                                 compute hybrid pure fluid props
      call hybeos (jns,3)

      call zeroys

      if (lmt) return

      y(5) = 0.00001d0
      y(1) = 0.1d0
c
      c1 = dexp (eqk(9) + fs2) 
      c2 = dexp (eqk(2) - eqk(3) - eqk(1)) 
      c3 = dexp (eqk(6) + fs2)
      c5 = dexp (eqk(8) + 4d0*(eqk(3)-eqk(2)) + 2d0*eqk(1) + fs2)/p
      c6 = dexp (eqk(4) + eqk(1) - eqk(3)) * p * p
c                                outer iteration loop: 
10    ek1 = c1 * g(3)/g(9) 
      ek2 = c2 * g(3) * g(1)/g(5)/g(2)
      ek3 = c3 * g(5)/g(6)     
      ek5 = c5 * g(2)**4 * g(5)**2/g(3)**4/g(8)/g(1)/g(1)
      ek6 = c6 * g(5)**3 * g(3)/g(1)/g(4)
c                                 solve for xh2, yco
      if (ifug.eq.19) then
         call evlxh2 (ek1,ek2,ek3,ek5,ek6,xo,xs,y(5),y(3),y(1),ier)
      else 
         call evlxh3 (ek1,ek2,ek3,ek5,ek6,xo,xs,y(5),y(3),y(1),ier)
      end if 

      if (ier.ne.0) call warn (501,y(5),ier,'XOXSRK')

      y(4) = ek6 * y(3) * y(5)**3/y(1)
      y(2) = ek2 * y(1) * y(3)/y(5)
      y(6) = ek3 * y(5)
      y(8) = ek5 * y(1)**2/y(5)**2
      y(9) = ek1 * y(3)

      nit = nit + 1

      if (nit.gt.iopt(21)) then
         call warn (175,y(5),ier,'XOXSRK')  
         goto 90
      end if     

      if (dabs(y(1)-oh2o).lt.nopt(50)) goto 90

      oh2o = y(1)

      call mrkhyb (ins, jns, 8, 3, 1)

      goto 10 

90    fh2o = dlog(g(1)*p*y(1))
      fco2 = dlog(g(2)*p*y(2))

      fo2 = 2d0 * (fh2o - dlog(g(5)*p*y(5)) - eqk(1))
c                                 compute graphite activity:
      ag = fco2 - fo2 - eqk(2)

      do i = 1, 3 
         vol = vol + y(jns(i))*vh(jns(i))
      end do 

      end


      subroutine hosrk5 (fo2,fs2)
c-----------------------------------------------------------------------
c program to calculate H-O-S speciation as a function of XO using
c an MRK/HSMRK hybrid. Species are H2 O2 H2O H2S SO2.

c this routine is probably not accurate cause it solves on yh2o. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision fo2,fs2,ek3,xom,xop,xos,c0,c1,c2,c3,
     *                 c4,c5,c6,c7,a,b,c,d,xl,xi,h,dh,oy5

      integer ins(5), jns(1), j, i

      logical lmt

      double precision vol
      common/ cst26 /vol

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save ins, jns
      data ins, jns/ 1,5,6,7,8,1/
c----------------------------------------------------------------------
c                                check if xo is <1, >0,
c                                reset if necessary
      call xcheck (xo,lmt) 
c                                 compute equilibrium constants in csteqk
      call seteqk (ins,5,-1d0)
c                                this fs2 = 1/2 ln (fs2),
      call setfs2 (fs2)
c                                 compute pure mrk fluid properties
      call mrkpur (ins,5)
c                                 compute hybrid pure fluid props
      call hybeos (jns,1)

      call zeroys

      if (lmt) return
c                               
      ek3 = dexp(eqk(1))
c                                 get first guess:
      if (xo.lt.r13) then 
         if (xo.gt.r13-nopt(50)) xo = r13 - nopt(50)
         xl = 2d0 *  xo/(1d0 - xo)
      else if (xo.ge.r13) then
         if (xo.lt.r13+nopt(50)) xo = r13 + nopt(50)
         xl = 2d0 * (1d0 - xo)/(1d0 + xo)
      end if

      xom = xo - 1d0
      xop = xo + 1d0
      xos = xo * xo

      c1 = dexp(eqk(6) + fs2)
      c3 = dexp(eqk(8) + fs2)
      c5 = 1d0/p/ek3/ek3
      a = -8d0*xo*xom**3
      b = -4d0*(3d0*xos+1d0)*xom**2 
      c0 = 2d0*xom*(-xop * (3d0*xo*xom + 2d0))
      c7 = 8d0*xom*c5
      d = -xom**2 * xop**2
      oy5 = 2d0 
c                                 outer iteration loop:
      do 30 j = 1, iopt(21)

         c2 = g(5)/g(6)
         c4 = g(7)/g(8)
         c6 = g(1)**2/g(5)**2/g(7)
         c =  c0 + c7*c6*(1d0+c1*c2)**2*(1d0+c3*c4)
         y(1) = xl
         xi = xl

         do i = 1, iopt(21)
c                                 inner iteration loop:
            h = a + (b + (c + d * y(1)) * y(1)) * y(1)
    
            dh = b + (2d0 * c  + 3d0 * d * y(1)) * y(1)
            y(1) = xi - h/dh

            y(5)  = -0.5d0*(xo*y(1)+y(1)+2d0*xo-2d0)/(1d0+c1*c2)
            y(6) = c1*c2*y(5)
            y(7) = c5*c6*(y(1)*y(1))/(y(5)*y(5))
            y(8) = c3*c4*y(7)
 
            if (dabs(xi-y(1)).lt.nopt(50)) goto 20
            if (y(1).ge.1d0) y(1) = xi + (1d0-xi)/2d0
            xi = y(1)

         end do

         call warn (176,y(1),i,'HOSRK5')
         stop 

20       y(5)  = -0.5d0*(xo*y(1)+y(1)+2d0*xo-2d0)/(1d0+c1*c2)
         y(6) = c1*c2*y(5)
         y(7) = c5*c6*(y(1)*y(1))/(y(5)*y(5))
         y(8) = c3*c4*y(7)

         if (dabs(oy5-y(5)).lt.nopt(50)) goto 40

         call mrkhyb (ins, jns, 5, 1, 1)

         oy5 = y(5) 

30    continue

      call warn (176,y(1),j,'HOSRK5')
      stop 

40    f(1) = dlog(g(5)*p*y(5)) 

      vol = vol + y(1) * vh(1)

      if (y(7).lt.y(5)) then
         fo2 = 2d0 * (dlog(g(1)*p*y(1)) - f(1) - dlog(ek3))
      else
         fo2 = dlog (g(7) * p * y(7))
      end if

      f(2) = fo2

      end 

      subroutine homrk (fo2)
c-----------------------------------------------------------------------
c program to calculate H-O speciation as a function of XO using
c an MRK/HSMRK hybrid. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision fo2,c1,c10,c11,c12,c13,c14,xl,xi,x12

      integer ins(3), jns(1), i, j

      logical lmt 

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision vol
      common/ cst26 /vol

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
      save / cst59 /

      save ins,jns
      data ins, jns/1,5,7,1/
c----------------------------------------------------------------------
c                                check if xo is <1, >0,
c                                reset if necessary
      call xcheck (xo,lmt)  
c                                 compute equilibrium constants in csteqk
      call seteqk (ins,1,-1d0)
c                                 compute pure mrk fluid properties
      call mrkpur (ins,3)
c                                 compute hybrid pure fluid props
      call hybeos (jns,1)

      call zeroys

      if (lmt) return

      c1 = 1d0/dsqrt(p)/eqk(1)

      if (xo.lt.r13) then 
c                                 this trap avoids singularity, a better approach
c                                 would be to acknowledge that xh2 and xo2 are not
c                                 independent at xo=1/3
         if (xo.gt.r13-nopt(50)) xo = r13-nopt(50)
         xl = 2d0 *  xo/(1d0 - xo)

      else if (xo.ge.r13) then

         if (xo.lt.r13+nopt(50)) xo = r13+nopt(50)
         xl = 2d0 * (1d0 - xo)/(1d0 + xo)

      end if 

      c10 = (xo - 1d0)/2d0
      c11 = 1d0 - xo
      c12 = 1d0 + c10
c                                 outer iteration loop:
      do 30 j = 1, iopt(21)
c                                 
         c13 = c1 * g(1)/g(5)/dsqrt(g(7))
         c14 = c10 * c13/2d0
         xi = xl

         do i = 1, iopt(21)
c                                 inner iteration loop:
            y(7) = xo + c10 * y(1)

            if (y(7).gt.nopt(50)) then  

               x12 = dsqrt (y(7))

               y(1) = xi + 
     *               (c11 - c12 * y(1) - c13 * y(1)/x12) /
     *               (c12 + c13 * x12 + c14 * y(1)/x12)

            else 

               y(1) = 2d0 *  xo/c11 

            end if 

            if (dabs(xi-y(1)).lt.nopt(50)) goto 20
            if (y(1).ge.1d0) y(1) = xi + (1d0-xi)/2d0
            xi = y(1)

         end do

         call warn (176,y(1),i,'HOMRK')

         goto 99

20       if (y(7).lt.0d0) y(7) = 0d0
         y(5) = 1d0 - y(7) - y(1)
         if (j.gt.1.and.dabs(xl-y(1)).lt.nopt(50)) goto 40

         call mrkhyb (ins, jns, 3, 1, 1)

         xl = y(1)

30    continue

      call warn (176,y(1),j,'HOMRK')
      goto 99 

40    f(1) = dlog(g(1)*p*y(1))

      vol = vol + y(1) * vh(1) 
      
      if (y(7).lt.y(5)) then
         fo2 = 2d0 * (f(1) - dlog(g(5)*p*y(5))- eqk(1))
      else
         fo2 = dlog (g(7) * p * y(7))
      end if

      f(2) = fo2

      return 

99    f(1) = dlog(1d12*p)
      f(2) = f(1)

      end 


      subroutine hh2ork (fo2,low)
c----------------------------------------------------------------------
c program to calculate fh2o, fo2 for H2-H2O mixtures using hsmrk/mrk
c hybrid EoS
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical low, lmt 

      integer ins(2), jns(1)

      double precision fo2

      double precision fh2o,fh2,funk
      common/ cst11 /fh2o,fh2,funk

      double precision p,t,xv,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xv,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      save ins, jns
      data ins, jns/1, 5, 1/
c----------------------------------------------------------------------
      call zeroys

      y(5) = xv
c                                check if xh2 is <1, >0,
c                                reset if necessary.
      call xcheck (y(5),lmt)

      y(1) = 1d0 - y(5)
c                                get pure species fugacities
      call hybeos (jns,1)
c                                 compute equilibrium constants in csteqk
      call seteqk (ins,1,-1d0)

      if (low) then 
c                                get mrk fugacities:
         call lomrk (ins, 2)

      else 
c                                get mrk fugacities:
         call mrkmix (ins, 2, 1)

      end if 
                    
      g(1) = gh(1) * g(1)

      fh2o = dlog(g(1)*p*y(1))

      fh2 = dlog(g(5)*p*y(5))
  
      fo2 = 2d0 * (fh2o - fh2 - eqk(1))

      vol = vol + y(1) * vh(1)
 
      end

      subroutine mrkhyb (ins,jns,isp,jsp,imix)
c---------------------------------------------------------------------
c routine to make mrk mixture calculations and convert to hybridized
c fugacity coefficients.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ins(*),isp,jns(*),jsp,imix
 
      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)
c----------------------------------------------------------------------
      call mrkmix (ins, isp, imix)

      do i = 1, jsp 
         g(jns(i)) = gh(jns(i)) * g(jns(i))
      end do 

      end 


      subroutine cohfo2 (fo2)
c----------------------------------------------------------------------
c subroutine to compute H2O and CO2 fugacities in a COH fluid
c consistent with a specifed graphite activity and oxygen fugacity
c specified by ibuf in routine fo2buf. Uses MRK activities with
c pure fluid fugacities as specified by the hybrid_EoS option.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ins(5),jns(3),nit

      double precision fo2,kh2o,kco2,kco,kch4,qa,qb,oh2o

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      save ins, jns
      data ins, jns/1,2,3,4,5,1,2,4/
c----------------------------------------------------------------------
      nit = 0

      call fo2buf (fo2)

      call seteqk (ins,5,elag)
c                                 compute pure mrk fluid properties
      call mrkpur (ins,5)
c                                 compute hybrid pure fluid props
      call hybeos (jns,3)

      call zeroys

      kco2 = dexp(eqk(2) + fo2)/p
      kco  = dexp(eqk(3) + fo2/2d0)/p
c                                check for graphite saturation:
      y(2) = kco2/g(2)
      y(3) = kco/g(3)

      if (y(2)+y(3).ge.1d0) then 
 
         write (*,1000) fo2,p,t
         fco2 = dlog(p*g(2))
         y(2) = 1d0
         y(3) = 0d0
         return

      end if  
c                                evaluate other k's:
      kh2o = dexp(eqk(1) + fo2/2d0)
      kch4 = dexp(eqk(4))*p
c                                solve for x's
      oh2o = 2d0

      do

         y(2) =  kco2/g(2)
         y(3) = kco/g(3)
c                                make quadratic 0 = qa * xh2**2 
c                                 + qb * xh2 + qc
         qb = kh2o * g(5)/g(1) + 1d0
         qa = kch4 * g(5)**2/g(4)

         y(5) = (-qb + dsqrt(qb**2 - 4d0*qa*(y(3) + y(2) -1d0)))/2d0/qa
         y(4) = kch4 * g(5)**2 * y(5)**2/g(4)
         y(1) = kh2o * g(5) * y(5)/g(1)
 
         nit = nit + 1

         if (nit.gt.iopt(21)) then 
            call warn (176,y(1),nit,'COHFO2')
            if (y(2)+y(3).gt.0.9999d0) then
               y(2) = 1d0
               y(1) = 1d-20
               call mrkpur (ins, 5)
               exit 
            else 
               stop
            end if 
            
         end if 

         if (dabs(y(1)-oh2o).lt.nopt(50)) exit

         oh2o = y(1)

         call mrkhyb (ins, jns, 5, 3, 1)

      end do

      do i = 1, 3 
         vol = vol + y(jns(i))*vh(jns(i))
      end do 


      xc = y(2) 

      if (hu.ne.1) then 

         fh2o = dlog(g(1)*p*y(1))
         fco2 = dlog(g(2)*p*y(2))

      else 

         fh2o = dlog(g(5)*p*y(5))
         fco2 = fo2

      end if 

1000  format ('**warning ver222** routine COHFO2, specified lnfO2 (',
     *        g12.6,')',/,'is inconsistent with graphite saturation',
     *        ' at P(bar)=',g12.6,' T(K)=',g12.6,/,'XCO2=1 assumed.') 

      end 

      double precision function hsmrkf (v,i)
c---------------------------------------------------------------------
c subprogram to get ln(fugacity) of pure Kerrick and Jacobs (1981) HSMRK 
c species i given an initial volume estimate v:

c   i = 1 - H2O
c   i = 2 - CO2
c   i > 2 - CH4

c on return v is the HSMRK volume of species i
c---------------------------------------------------------------------
      implicit none

      logical bad

      integer i

      double precision b,bw,bc,bm,rtt,t2,c,d,e,yz,fugp,rr,v,t12

      external fugp

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      save bw, bc, bm, rr
      data bw, bc, bm, rr /29d0, 58d0, 60d0, 83.144126d0/
c----------------------------------------------------------------------

      t12 = dsqrt(t)
      rtt = rr*dsqrt(t**3)
      t2 = t*t

      if (i.eq.1) then 
         b = bw
         c = 290.78d6-0.30276d6*t+0.00014774d6*t2
         d = -8374d6+19.437d6*t-0.008148d6*t2
         e = 76600d6-133.9d6*t+0.1071d6*t2
      else if (i.eq.2) then 
         b = bc 
         c = 28.31d6+0.10721d6*t-0.00000881d6*t2
         d = 9380d6-8.53d6*t+0.001189d6*t2
         e = -368654d6+715.9d6*t+0.1534d6*t2
      else 
         b = bm 
         c = 13.403d6 + 9.28d4 * t + 2.7d0 * t2
         d = 5.216d9 - 6.8d6 * t + 3.28d3 * t2
         e = -2.3322d11 + 6.738d8 * t + 3.179d5 * t2
      end if

      call nurap (bm,c,d,e,yz,v,t12,rr,bad)

      if (bad) then 
         hsmrkf = dlog(1d12*p)
      else 
         hsmrkf = dlog(p) + fugp (rtt,bm,yz,c,d,e,v)
      end if

      end 

      subroutine nurap (b,c,d,e,yz,vi,t12,r,bad)
c----------------------------------------------------------------------
c newton-raphson iteration to solve for hsmrk volume, my
c idea here was to evaluate all the constants outside of 
c the iteration loop, which might make sense if a lot of
c iterations were necessary, as it is i'm unsure. someone
c should do a comparison with my old nurap routine (xnurap).
c                                 jadc, july 96.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad 

      integer k

      double precision b,c,d,e,yz,vi,s1,s2,s3,p0,b2,q0,q1,q2,q3,q4,q5,
     *                 q6,q7,q8,q9,p1,p2,p3,p4,p5,p6,p7,p8,cor,t12,r

      double precision p,t,xc,u1,u2,tr,pr,r0,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r0,ps

      double precision vol
      common/ cst26 /vol
c----------------------------------------------------------------------
      s1 = r*t*t12
      s2 = b*s1
      s3 = t12*p*b
      p0 = -256d0*s1
      b2 = b*b
      q0 = 256d0*t12*p
      q1 = 256d0*(s3 - s1)
      q2 = (-160d0*s3 - 512d0*s1)*b + 256d0*c
      q3 = (-80d0*s3 + p0)*b2 + 256d0*d
      q4 = ((65d0*s3 + 8d0*s1)*b - 160d0*c)*b2 + 256d0*e
      q5 = -b2*(((14d0*s3 - 15d0*s1)*b - 80d0*c)*b + 160d0*d)
      q6 = b2*((((s3 + 6d0*s1)*b - 15d0*c)*b + 80d0*d)*b - 160d0*e)
      q7 = b**3*(((-s1*b + c)*b - 15d0*d)*b + 80d0*e)
      q8 = b**4*(-15d0*e + d*b)
      q9 = e*b**5
      p1 = 512d0*c-768d0*s2
      p2 = (-832d0*s2-256d0*c)*b+768d0*d
      p3 = ((-368d0*s2 - 64d0*c)*b - 256d0*d)*b + 1024d0*e 
      p4 = -b*(((33d0*s2 - 64d0*c)*b + 224d0*d)*b + 256d0*e)
      p5 = 2d0*b2*(b*((s2-c)*7d0*b + 72d0*d) - 192d0*e)
      p6 = -b**3*(b*((s2-c)*b + 29d0*d) - 224d0*e)
      p7 = 2d0*b**4*(d*b-22d0*e)
      p8 = 3d0*q9

      k = 0
      bad = .false.

      do 

         cor = ((((((((((q0*vi+q1)*vi+q2)*vi+q3)*vi+q4)*vi+q5)*vi
     *         +q6)*vi+q7)*vi+q8)*vi+q9)*vi)/((((((((p0*vi+p1)*vi
     *         +p2)*vi+p3)*vi+p4)*vi+p5)*vi+p6)*vi+p7)*vi+p8)
         vi = vi + cor

         if (dabs(cor/vi).lt.nopt(50)) then
            exit
         else if (vi.lt.0d0) then 
            bad = .true.
            exit 
         else 
            k = k + 1
            if (k.le.iopt(21)) cycle
            bad = .true.
            exit 
         end if

      end do 
 
      yz = vi * p/r/t

      vol = vi
 
      end

      subroutine mrkmix (ins, isp, iavg)
c-----------------------------------------------------------------------
c subroutine to calculate the log fugacities and volume of mixed
c species fluids using the RK/MRK EoS. 

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp    - the number of species to be calculated.
c        p,t    - input from cst5, p(bars), t(K)
c        iavg   - a flag indicating whether the a cross term is to 
c                 be computed as the geometric mean (iavg = 1), the
c                 arithmetic mean (iavg = 2), or the harmonic mean; 
c                 provided the cross term is not specified explicitly, 
c                 as in the deSantis et al. 1974 formulation for water
c                 and CO2.  

c output (to common cstcoh):

c        g(i)    - fugacity coefficient of ith species
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(*),i,j,k,l,iroot,isp,ineg,ipos,iavg,jrt

      logical max
 
      double precision f(nsp),aj2(nsp),ev(3),c1,c2,ax,dv,ch,bx,aij,pdv,
     *                 c3,vmin,vmax,d1,d2,d3,d4,d5,d6,d7,rt,dsqrtt,r,
     *                 vpb, vmb, yl

      double precision p,t,xco2,u1,u2,tr,pr,rbar,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,rbar,ps

      double precision fg
      common/ cst11 /fg(3) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision vol
      common/ cst26 /vol

      double precision a, b
      common/ rkab /a(nsp),b(nsp)

      double precision vrt
      integer irt
      logical sroot
      common/ rkroot /vrt,irt,sroot

      integer iroots
      logical switch, rkmin, min
      common/ rkdivs /iroots,switch,rkmin,min
 
      save r
      data r, max /83.1441d0, .false./
c---------------------------------------------------------------------- 
      dsqrtt = dsqrt(t)
      rt = r*t

      call rkparm (ins,isp)

      bx = 0d0
      aij = 0d0

      do k = 1, isp

         i = ins(k)
c DEBUG DEBUG 691 RATSO
         if (y(i).lt.0d0) y(i) = 0d0

         aj2(i) = 0d0
         bx = bx + b(i)*y(i)

      end do 
 
      ch = dexp(-11.218d0 + (6032d0 + (-2782000d0 + 4.708d8/t)/t)/t) * 
     *          6912.824964d0 *t*t*dsqrtt + 79267647d0

      do k = 1, isp

         i = ins(k)

         do l = 1, isp

            j = ins(l)
 
            if (i.eq.1.and.j.eq.2.or.i.eq.2.and.j.eq.1) then
               aij = aij + y(i)*y(j)*ch/2d0
               aj2(i) = aj2(i) + y(j)*ch
            else 

               if (i.eq.14.and.j.eq.15.or.i.eq.15.and.j.eq.14) then
c                                 special mixing rule
                  ax = 2d0/(1d0/a(i) + 1d0/a(j))
               else if (iavg.eq.1) then 
c                                 geometric mean mixing rule
                  ax = dsqrt(a(i)*a(j))
               else if (iavg.eq.2) then 
c                                 arithmetic mean mixing rule
                  ax = (a(i) + a(j))/2d0
               else 
c                                 harmonic mean mixing rule
                  ax = 2d0/(1d0/a(i) + 1d0/a(j))
               end if 

               aij = aij + y(i)*y(j)*ax
               aj2(i) = aj2(i) + 2d0*y(j)*ax

            end if

         end do 

      end do 
c                                 solve for mixture molar volume
      c1 = -rt/p
      c3 = -aij*bx/p/dsqrtt
      c2 = c1*bx + aij/dsqrtt/p - bx*bx

      call roots3 (c1,c2,c3,ev,vmin,vmax,iroots,ineg,ipos)

      if (sroot) then 
c                                use characteristics of previous solution 
c                                to choose the, potentially metastable, root
         if (irt.eq.3.and.iroots.eq.3.and.ineg.eq.0.and.vmin.gt.bx) then

            if (max) then 
               vol = vmax
            else 
               vol = vmin
            end if 

         else if (iroots.eq.3.and.irt.eq.3) then

            vol = vmax

         else 
c                                 different number of roots, minimize
c                                 the difference
            dv = 1d99

            do i = 1, iroots 
               if (ev(i).lt.0d0) cycle
               if (dabs(ev(i) - vrt).lt.dv) then
                  jrt = i
                  dv = dabs(ev(i) - vrt)
               end if 
            end do

            if (dv.eq.1d99) then 
               write (*,*) 'rats'
            else 
               vol = ev(jrt)
            end if 

         end if 

      else if (iroots.eq.3.and.ineg.eq.0.and.vmin.gt.bx) then
c                                choose the root with lowest gibbs energy
c                                by evaluating p*delta(v) - int(pdv)
         pdv = p*(vmax-vmin) - 
     *         dlog((vmax-bx)/(vmin-bx))  * rt - 
     *         dlog((vmax+bx)/(bx+vmin)*vmin/vmax) * aij/bx/dsqrtt

         if (pdv.gt.0d0) then
            vol = vmin
            max = .false.
         else 
            vol = vmax
            max = .true.
         end if 

      else if (iroots.eq.3) then

         vol = vmax

      else 

         vol = ev(ipos)

      end if

      if (.not.sroot) then 
c                                 for finite difference property 
c                                 computations, save the root
         irt = iroots
         vrt = vol

      end if

      if (vol.lt.0d0) then
         write (*,*) 'rats II'
         return
      end if 

      iroot = iroots
c                                 compute fugacities and pmvs:
      vpb = vol + bx
      vmb = vol - bx
      d1 = rt*dsqrtt*bx
      d2 = dlog(vpb/vol)/d1
      d3 = aij*d2/bx - aij/vpb/d1 + 1d0/vmb 
      d6 = dlog(rt/vmb)
c                                 pmv constants
      d4 = vmb**2/vpb/(rt*dsqrtt)/vol
      d5 = aij*d4*(1d0/vol+1d0/vpb) - 1d0
      d7 = -d4*aij/vpb

      do i = 1, isp

         l = ins(i)
         
         if (y(l).gt.nopt(50)) then 
            yl = y(l)
         else
            yl = nopt(50)
         end if

         f(l) = dlog(yl) + b(l)*d3 - aj2(l)*d2 + d6
         g(l) = dexp(f(l))/p/yl
c                                 pmv
         v(l) = (d4*aj2(l) - b(l) - vmb + d7*b(l)) / d5

         if (l.lt.3) fg(l) = f(l)

      end do 
  
      end

      subroutine crkh2o (pbar,t,vol,fh2o)
c-----------------------------------------------------------------------
c compute ln(f[h2o], bar) and volume h2o (cm3/mol) from CORK EoS Holland 
c & Powell CMP 109:265-273. Input pbar - pressure (bars); tk - temp (K).
c                                 J.A.D. Connolly, 1992.
c-----------------------------------------------------------------------
      implicit none

      double precision x(3),b,r,p0,rt,rtp,t12,psat,t,pbar,vol,
     *                 fh2o,a1,a2,a3,xmin,xmax,cc,gam,dp,c,d,e,vmin,
     *                 vmax,gam2,a,p

      integer iroots, i, ineg, ipos

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
      save / cst59 /

      save b,r,p0
      data b,r,p0 /1.465d0,8.314d-3,2d0/
c----------------------------------------------------------------------

      p = pbar/1d3
      rt = r*t
      rtp = rt/p
      t12 = dsqrt(t)

      if (t.lt.695d0) then
         psat = -13.627d-3 + t**2*(0.729395d-6 - 0.234622d-8*t 
     *          + 0.483607d-14*t**3)

         if (p.lt.psat.and.t.lt.673d0) then
c                          this is a(gas)
             a = 16138.87d0 
     *           - t*(69.66291d0 - t*(0.1161905d0 - 0.68133d-4*t))
         else
            if (t.lt.673d0) then 
               a = -1449.009d0
     *             + t*(12.70068d0 - t*(0.02208648d0 - 0.13183d-4*t))
            else
               a = 1036.975d0 
     *             + t*(0.5306079d0 - t*(0.7394203d-3 - 0.17791d-6*t))
            end if
         end if
      else 
          psat = 0d0 
          a = 1036.975d0 
     *        + t*(0.5306079d0 - t*(0.7394203d-3 - 0.17791d-6*t))
      end if

      a1 = -rtp
      a2 = a/t12/p - b*(rtp+b)
      a3 = -a*b/t12/p

      call roots3 (a1,a2,a3,x,xmin,xmax,iroots,ineg,ipos) 

      if (iroots.eq.1) then
         vol = x(1)
      else 

         if (p.lt.psat) then
            vol = xmax
         else if (t.lt.7d2.and.xmin.gt.0d0) then
c                                 this seems hokey, buy there must be
c                                 a reason for it.
            vol = xmin
         else 
c                                 cork has 3 roots sometimes at
c                                 high P, this is a fitting artifact
c                                 and in my experience the only real
c                                 root is positive.
            do i = 1, 3
               if (x(i).gt.0d0) then
                  vol = x(i)
                  goto 2
               end if
            end do 
         end if
      end if 

2     cc = a/b/rt/t12

      gam = vol/rtp - 1d0 - dlog((vol-b)/rtp) - cc*dlog(1d0+b/vol)

      if (p.gt.p0) then
c                            new tjbh cork, 97:
         dp = p - p0
         c = 1.9853d-3*dp
         d = -8.909d-2*dsqrt(dp)
         e = 8.0331d-2*dp**0.25d0

         vol = vol + c + d + e

         gam = gam + dp*(c/2d0 + d*r23 + 0.8d0*e)/rt

      end if
c                              add check to keep T > 273
      if (t.lt.695d0.and.p.gt.psat.and.t.gt.273d0) then
c                              tinkham version, Sept. 12, 2002
         p = psat
         rtp = rt/p
c
         a1 = -rtp
         a2 = a/t12/p - b*(rtp+b)
         a3 = -a*b/t12/p

         call roots3 (a1,a2,a3,x,vmin,vmax,iroots,ineg,ipos)
c                                 gamma liq 
         gam2 = vmin/rtp-1d0-dlog((vmin-b)/rtp)-cc*dlog(1d0+b/vmin)
c
         if (t.lt.673d0) then 
c                                 gas phase has a different a, recalculate
c                                 vmax 
             a = 16138.87d0  - t*(69.66291d0  
     *                       - t*(0.1161905d0 - 0.68133d-4*t))
c                                 prior to Jun 11, 2004, cc was not
c                                 updated for new a value, pj gorman.
             cc = a/b/rt/t12

             a1 = -rtp
             a2 = a/t12/p - b*(rtp+b)
             a3 = -a*b/t12/p        
             call roots3 (a1,a2,a3,x,vmin,vmax,iroots,ineg,ipos)  
c                                 could check that there really are three real roots
         end if 

         gam = vmax/rtp-1d0-dlog((vmax-b)/rtp)-cc*dlog(1d0+b/vmax) 
     *         - gam2 + gam

      end if 
c                                 convert volume from j/bar to cm3/mol 
      vol = vol * 1d1

      fh2o = gam + dlog(pbar)

      end

      subroutine crkco2 (pbar,t,vol,fco2)
c-----------------------------------------------------------------------
c compute ln(f[co2], bar) and volume co2 (cm3/mol) from CORK EoS Holland 
c & Powell CMP 109:265-273. Input pbar - pressure (bars); tk - temp (K).
c                                 J.A.D. Connolly, 1992.
c-----------------------------------------------------------------------
      implicit none

      integer i, iroots, ineg, ipos

      double precision x(3),pbar,t,vol,fco2,b,r,p0,p,dp,rt,rtp,t12,a,
     *                 cc,a1,a2,a3,xmin,xmax, c, d

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
      save / cst59 /

      save b,r,p0

      data b,r,p0/3.057d0,8.314d-3,5d0/
c----------------------------------------------------------------------
      p = pbar/1d3
      rt = r*t
      rtp = rt/p
      t12 = dsqrt(t)

c      a = 741.2 - t*(0.10891 + 3.4203d-4*t)
c                                 a from roger:
      a = 659.8 + 0.21078 * t - 6.3976d-4 * t * t
c
      a1 = -rtp
      a2 = a/t12/p - b*(rtp+b)
      a3 = -a*b/t12/p 
c                                 get volume:
      call roots3 (a1,a2,a3,x,xmin,xmax,iroots,ineg,ipos) 

      if (iroots.eq.1) then
         vol = x(1)
      else         
c                                 cork has 3 roots sometimes at
c                                 high P or T, this is a fitting 
c                                 artifact and in my experience 
c                                 the only real root is positive.
         do i = 1, 3
            if (x(i).gt.0d0) then
               vol = x(i)
               goto 2
            end if
         end do

         call error (999,xmax,iroots,'CRKCO2')

      end if 

2     cc = a/b/rt/t12

      fco2 = dlog(pbar) + vol/rtp - 1d0 
     *     - dlog((vol-b)/rtp) - cc*dlog(1d0+b/vol)

      if (p.gt.p0) then
c                                 add virial component:
c                                 coefficients from tjbh '95
          dp = p - p0
c                       
c         d = 5.40776d-3 - 1.59046d-6*t
c         c = -1.78198d-1 + 2.45317d-5*t

c         fco2 = fco2 + (2d0*c*dp**(1.5d0)/3d0
c     *               + d/2d0*dp*dp)/rt

          c = 1.5d0 * (0.1967099672d-2 - 14.28899046d0/t)
          d = 2d0 * (0.3252201107d0/t - 0.9564950686d-4)

          vol = vol + c + d

          fco2 = fco2 + dp*(c*r23*dsqrt(dp) + d/2d0*dp)

      end if
c                                 convert volume from j/bar to cm3/mol 
      vol = vol*1d1

      end


      subroutine roots3 (a1,a2,a3,x,vmin,vmax,iroots,ineg,ipos)
c----------------------------------------------------------------------
c returns real roots (in x) of: x**3 + a1*x**2 + a2*x + a3
c     iroots - the number of real roots
c     ineg - the number of negative real roots, if iroots = 3
c     ipos - the index of the real positive root in x
c----------------------------------------------------------------------
      implicit none 

      integer iroots, i, ineg, ipos 

      double precision x(*),a1,a2,a3,a4,a5,a6,rr,qq,vmin,vmax,dif,phi,
     *                 a7,dphi,v
c----------------------------------------------------------------------
      qq = (a1**2 - 3d0 * a2) / 9d0
      rr = (a1 * (2d0 * a1**2 - 9d0*a2) + 27d0 * a3) / 54d0
      a5  = a1 / 3d0

      dif = qq**3 - rr**2

      if (dif.ge.0d0) then

         if (dif.gt.0d0) then 
            phi = dacos( rr / qq**(1.5d0) )
         else
            phi = 0d0
         end if 

         a4  = -2d0 * dsqrt(qq)
         a6  = phi / 3d0

         dphi = 0d0
         vmin = 1d9
         vmax = -1d9
         ineg = 0 

         do i = 1, 3
            v = a4 * dcos(a6 + dphi) - a5
            if (v.gt.vmax ) vmax = v
            if (v.lt.vmin ) vmin = v
            if (v.le.0d0  ) then 
               ineg = ineg + 1
            else 
               ipos = i
            end if 
            x(i) = v
            dphi = dphi + 2.094395102497915d0
         end do 

         iroots = 3

      else

         a7 = ( dsqrt(-dif) + dabs(rr) )**(1d0/3d0)
         x(1) = -rr / dabs(rr) * ( a7 + qq/a7 ) - a5
         iroots = 1
         ineg = 0 
         ipos = 1

      end if 

      end

      subroutine lomrk (ins, isp)
c-----------------------------------------------------------------------
c a semi-fudged (see comments below) MRK to calculate the fugacities
c of mixed species fluids in the vicinity of the water critical point. 

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp     - the number of species to be calculated.
c        p,t     - input from cst5, p(bars), t(K)

c output (to common cstcoh):

c        g(i)    - fugacity coefficient of ith species
c        v(i)    - volume of the ith species

c species indices:

c         1 = H2O
c         2 = CO2
c         3 = CO
c         4 = CH4 
c         5 = H2
c         6 = H2S
c         7 = O2
c         8 = SO2
c         9 = COS
c        10 = N2
c        11 = NH3
c        12 = O
c        13 = SiO
c        14 = SiO2
c        15 = Si  
c        16 = C2H6
c        17 = HF
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(*), i, j, k, l, isp, iroots, ineg, ipos
 
      double precision f(nsp),aj2(nsp),ev(3),t2,rt,d2,
     *                 dsqrtt,ch,bx,aij,c1,c2,c3,aij12,vmin,vmax,vol,
     *                 d1,d3,d4,d5,d6,r

      double precision p,t,xco2,u1,u2,tr,pr,rbar,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,rbar,ps

      double precision fg
      common/ cst11 /fg(3) 
 
      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision a, b
      common/ rkab /a(nsp),b(nsp)

      save r 
      data r /83.1441d0/
c----------------------------------------------------------------------
      t2 = t*t
      dsqrtt = dsqrt(t)
      rt = r*t
c
      call rkparm (ins,isp)
c                             this is a cheap trick the a-fit 
c                             was derived using b=15, then b
c                             was adjusted to lower the H2-H2O
c                             critical T, however the fugacities
c                             seem reasonable so who knows?
      a(1) =  .3930568949d9 - .1273025840d7 * t + 2049.978752 * t2
     *       -1.122350458 * t2 * t
c                             compute dispersion term for co2
      a(2) =  0.9293554d8 - 0.8213073d5*t + 0.2129d2*t2
 
      ch = dexp(-11.218d0 + 6032d0/t - 2782d3/t2 + 4.708d8/t2/t) * 
     *          6912.824964d0 * t2 * dsqrtt + 79267647d0
c                             composition dependent mrk-terms
      bx = 0d0
      aij = 0d0
 
      do k = 1, isp
         i = ins(k)
         aj2(i) = 0d0
         bx = bx + b(i)*y(i)
      end do 
 
      do k = 1, isp
         i = ins(k)
         do l = 1, isp
            j = ins(l)
 
            if (i.eq.1.and.j.eq.2.or.i.eq.2.and.j.eq.1) then
               aij = aij + y(i)*y(j)*ch/2d0
               aj2(i) = aj2(i) + y(j)*ch
            else
               aij12 = y(j)*dsqrt(a(i)*a(j))
               aij = aij + y(i)*aij12
               aj2(i) = aj2(i) + 2d0*aij12
            end if
         end do 
      end do 
c                                 solve for molar volume of mixture
      c1 = -rt/p
      c3 = -aij*bx/p/dsqrtt
      c2 = c1*bx + aij/dsqrtt/p - bx*bx

      call roots3 (c1,c2,c3,ev,vmin,vmax,iroots,ineg,ipos)

      if (iroots.eq.3) then
         vol = vmax
      else
         vol = ev(1)
      end if
c                                 segment to compute the
c                                 fugacities.
      d2 = dlog((vol+bx)/vol)
      d1 = rt*dsqrtt*bx
      d3 = d2 -bx/(bx+vol)
      d4 = vol -bx
      d5 = aij*d3/d1/bx
      d6 = dlog(rt/d4)
 
      do 70 i = 1, isp

         l = ins(i)
 
         if (y(l).gt.0d0) goto 80
         f(l) = 0d0
         g(l) = 1d0
         goto 75
80       f(l) = dlog(y(l))+b(l)/d4-aj2(l)/d1*d2+b(l)*d5+d6
         g(l) = dexp(f(l))/p/y(l)
75       if (l.lt.3) fg(l) = f(l)
70    continue
  
      end

      subroutine mrkpur (ins, isp)
c-----------------------------------------------------------------------
c subroutine to calculate the log fugacities and volume of single
c species fluids using the RK/MRK EoS. 

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp     - the number of species to be calculated.
c        p,t     - input from cst5, p(bars), t(K)

c output (to common cstcoh):

c        g(i)    - fugacity coefficient of ith species
c        v(i)    - volume of the ith species

c species indices:

c         1 = H2O
c         2 = CO2
c         3 = CO
c         4 = CH4 
c         5 = H2
c         6 = H2S
c         7 = O2
c         8 = SO2
c         9 = COS
c        10 = N2
c        11 = NH3
c        12 = O
c        13 = SiO
c        14 = SiO2
c        15 = Si  
c        16 = C2H6
c        17 = HF
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision f(nsp),ev(3),dsqrtt,rt,
     *                 d1,d2,d4,bx,v1,v2,aij,c1,c2,c3,pdv,r

      integer ins(*), isp, k, iroots, i, ineg, ipos 
 
      double precision p,t,xco2,u1,u2,tr,pr,rbar,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,rbar,ps

      double precision fg
      common/ cst11 /fg(3) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision vol
      common/ cst26 /vol

      double precision a, b
      common/ rkab /a(nsp),b(nsp)
 
      save r
                             
      data r /83.1441/
c----------------------------------------------------------------------
      dsqrtt = dsqrt(t)
      rt = r*t

      call rkparm (ins,isp)

      do k = 1, isp

         i = ins(k)

         aij = a(i)
         bx = b(i)

         c1 = -rt/p
         c3 = -aij*bx/p/dsqrtt
         c2 = c1*bx + aij/dsqrtt/p - bx*bx
c                                 v2 is the max root, v1 the min
         call roots3 (c1,c2,c3,ev,v1,v2,iroots,ineg,ipos)

         if (iroots.eq.3.and.ineg.eq.0.and.v1.gt.bx) then
c                                choose the root with lowest gibbs energy
c                                by evaluating p*delta(v) - int(pdv)
            pdv = p*(v2-v1) - 
     *            dlog((v2-bx)/(v1-bx))  * rt - 
     *            dlog((v2+bx)/(bx+v1)*v1/v2) * aij/bx/dsqrtt

           if (pdv.gt.0d0) then
              vol = v1
           else 
              vol = v2
           end if 

         else if (iroots.eq.3) then

            vol = v2

         else
c                                 assume only one positive real root,
c                                 could check for no positive or 2 positive
            vol = ev(ipos)

         end if
c                                 compute fugacities.
         d1 = vol + bx
         d2 = dlog(d1/vol)
         d4 = vol - bx

         v(i) = vol
         f(i) = bx/d4 - (1d0/d1 + d2/bx)*aij/rt/dsqrtt + dlog(rt/d4)

         if (i.lt.3) fg(i) = f(i)
         g(i) = dexp(f(i))/p 

      end do 

      end

      function fug (rtt,cij,dij,eij,x1,x2,b,yz,c,d,e,b1,c1,d1,e1)
c---------------------------------------------------------------------
c streamlined JADC July 96.
c---------------------------------------------------------------------
      implicit none

      double precision rtt,cij,dij,eij,x1,x2,b,yz,c,d,e,b1,c1,d1,e1,y,
     *                 vpb,dvbv,y1,dvbvb,vi,vi2,fug
 
      double precision v
      common/ cst26 /v
c----------------------------------------------------------------------
      y = b/4d0/v
      vpb = v + b
      dvbv = dlog(vpb/v)
      y1 = 1d0 - y
      dvbvb = dvbv/b
      vi = 1d0/v
      vi2 = .5d0/v/v

      fug =  ((4d0-3d0*y)*y + b1/b*(2d0-y)*y*2d0/y1)/y1/y1
     *     + (c*b1*(dvbvb - 1d0/vpb)
     *     - (c1*x1+cij*x2)*2d0*dvbv
     *     + (2d0*(x1*d1+dij*x2)+d)*(dvbvb - vi)
     *     + d*b1*((vi+2d0/b)/vpb - 2d0*dvbvb/b)
     *     + 2d0*(e1*x1+eij*x2+e)*((vi - dvbvb)/b - vi2) 
     *     + e*b1*((vi2 - (1.5d0/v+3d0/b)/b)/vpb + 3d0*dvbvb/b/b)
     *       )/rtt/b-dlog(yz)
 
      end
 
      double precision function fugp (rtt,b,yz,c,d,e,v)
c----------------------------------------------------------------------
c streamlined JADC July 96.
c----------------------------------------------------------------------
      implicit none

      double precision y,y1,vpb,dvbv,rtt,b,yz,c,d,e,v
c----------------------------------------------------------------------
  
      y = b/4d0/v
      vpb = v + b
      dvbv = dlog(vpb/v)
      y1 = 1d0 - y

      fugp = ((4d0-3d0*y)*y + (2d0-y)*2d0*y/y1)/y1/y1
     *     + (d*(dvbv/b + (4d0*y + 2d0)/vpb - 3d0/v)
     *      - c*(dvbv + b/vpb)
     *      + e*((4d0/b-2d0/v)/v - dvbv/b/b 
     *         + ((2d0*y - 1.5d0)/v - 3d0/b)/vpb)
     *        )/rtt/b - dlog(yz)

      end

 
      subroutine brmrk (v2,fco2)
c-----------------------------------------------------------------------
c brmrk routine to calculate molar volume and fugacity of co2 from
c bottinga and richet 1981, ajs v 281.
c                                 j. connolly 1988.
c
c  for the unitiated, input and output is done through common blocks
c  cst5, cst11, and cst26, the relevant variables are:
c
c  for input:
c               p    = pressure in bars
c               t    = temperature in kelvins
c               r    = gas constant
c               xco2 = mole fraction of co2 in fluid phase
c  for output:
c               v    = molar volume (cm3/mole) at p and t
c               fco2 = the natural log of the co2 fugacity
c-----------------------------------------------------------------------
      implicit none

      double precision vdpdv,v1,v2,f1,f2,fco2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 external function to calculate
c                                 the product v(dp/dv)
      external vdpdv
c------------------------------------------------------------------------
c                                 call brvol to calculate volume
c                                 of co2, 1st 1 bar, use ideal gas guess
      v1 = 1d1*r*t/p
      call brvol (1d0,t,v1)
c
      call brvol (p,t,v2)
c                                 call qromb to evaluate integral
c                                 of (dp/dv)dv
      if (v2.ge.180d0) then
c                                 if b&r equation is continuous
c                                 through limits v1 to v do 
c                                 integration in one chunk.
         call qromb (vdpdv,v1,v2,fco2)
c
      else if (v2.gt.47.22d0) then
c                                 else do the integration in parts
         call qromb (vdpdv,v1,180d0,f1)
         call qromb (vdpdv,180d0,v2,f2)
         fco2 = f1 + f2
c
      else
c
         call qromb (vdpdv,v1,180d0,f1)
         call qromb (vdpdv,180d0,47.22d0,f2)
         call qromb (vdpdv,47.22d0,v2,fco2)
         fco2 = fco2 + f1 + f2
 
      end if

      fco2 = fco2/(1d1*r*t)
 
      end
 
      subroutine brvol (p,t,vol)
c-----------------------------------------------------------------------
c brvol function to calculate molar volume of co2 by newton-raphson
c iteration of eq 3 of bottinga and richet 1981, ajs v 281.
 
c input:  p = pressure in bars, t = temperature in kelvins, v = molar
c         volume of co2 (cm3) estimated by the mrk equation of state.
 
c output: brvol = molar volume of co2 in (cm3/mole)
 
c                                 j. connolly 1988.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer itic

      double precision dv,a1,a2,a3,rbar,rt,t12,b1,b2,a,b,vp,bp,fv,p,
     *                 dpdv,t,ap,vol

      double precision v
      common/ cst26 /v
 
      save rbar, dv, a1, a2, a3
 
      data rbar,dv/83.143d0,0.001d0/,a1,a2,a3/6.566d7,7.276d7,37.3d0/
c---------------------------------------------------------------------
      rt = rbar * t
      t12 = dsqrt(t)
      itic = 0
c                                 brvol appoximates the derivative dp/dv
c                                 by finite difference with dv (in cm3)
      dv = 0.00005d0
c                                 choose parameters for b&r's equation
10    if (v.le.47.22d0) then
         b1 = 1.856669d0
         b2 = 0.0637935d0
      else if (v.lt.180d0) then
         b1 = 11.707864d0
         b2 = 0.363955d0
      else
         b1 = 7.352629d0
         b2 = 0.241413d0
      end if
c                                 at v:
      b = (dlog(v/a3)+b1)/b2
      a = a1*( (a3/v)**3 - (a3/v)**6 ) + a2
c                                 and at vp:
      vp = v + dv
c
      bp = (dlog(vp/a3)+b1)/b2
      ap = a1*( (a3/vp)**3 - (a3/vp)**6 ) + a2
c
      fv = rt/(v-b) - a/(v*(v+b)*t12) - p
      dpdv = (fv - (rt/(vp-bp) - ap/(vp*(vp+bp)*t12) - p))/dv
c                                 refine estimate:
      v = v + fv/dpdv
c                                 if the change from the last estimate
c                                 is less than 0.001 cm3 leave.
      if (dabs(fv/dpdv).lt.1d-3) goto 99
      itic = itic + 1
      if (itic.gt.50) goto 999
 
      goto 10
 
999   call warn (176,ap,itic,'BRVOL')
      stop
99    vol = v
      end
 
      function vdpdv (v)
c-----------------------------------------------------------------------
c dpdv function to calculate the finite difference approximation of
c the product v(dp/dt) for the equation of bottinga and richet 1981,
c ajs v 281, eq 3.
c input: v - molar volume of co2 in j/bar
c        t - temperature in k (common block cst5)
c output:vdpdv - the approximation of v(dp/dv)

c                                 j. connolly 1988.
c-----------------------------------------------------------------------
      implicit none

      double precision rbar,dv,a1,a2,a3,t12,b1,b2,v,bp,ap,vdpdv,rt,a,b,
     *                 vp,pp

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps
 
      save rbar, dv, a1, a2, a3
 
      data rbar,dv/83.143d0,1d-3/a1,a2,a3/6.566d7,7.276d7,37.3d0/
c---------------------------------------------------------------------
 
      rt = rbar * t
      t12 = dsqrt(t)

      if (v.le.47.22d0) then
         b1 = 1.856669d0
         b2 = 0.0637935d0
      else if (v.lt.180d0) then
         b1 = 11.707864d0
         b2 = 0.363955d0
      else
         b1 = 7.352629d0
         b2 = 0.241413d0
      end if

      b = (dlog(v/a3)+b1)/b2
      a = a1*( (a3/v)**3 - (a3/v)**6 ) + a2
 
      vp = v + dv
 
      bp = (dlog(vp/a3)+b1)/b2
      ap = a1*( (a3/vp)**3 - (a3/vp)**6 ) + a2
 
      pp = rt/(v-b) - a/(v*(v+b)*t12)
      vdpdv = -v*(pp - (rt/(vp-bp) - ap/(vp*(vp+bp)*t12) ))/dv
 
      end
 
      subroutine simps (func,a,b,dx,ss)
c------------------------------------------------------------------
c subroutine to numerically integrate the function func over the
c limits a to b, in n increments by simpsons method, the integral
c is returned as ss. from conte and deboor 1981.
c                               jadc, aug 1990.
c------------------------------------------------------------------
      implicit none
 
      integer i,n

      double precision func,b,a,dx,ss,x,h,hover2,half
 
      external func
 
      n = idint (dabs(b-a)/dx)
      if (n.lt.100) n = 100
 
      h = (b-a)/dfloat(n)
      hover2 = h/2d0
      ss = 0d0
 
      half = func (a + hover2)
 
      do i = 1, n-1
         x = a + h * dfloat(i)
         ss = ss + func(x)
         half = half + func (x + hover2)
      end do 
 
      ss = (h/6d0) * (func(a) + 4d0*half + 2d0*ss + func(b))
 
      end
 
      subroutine qromb (func,a,b,ss)
c-----------------------------------------------------------------------
c subroutine to numerically integrate a function func over the limits
c a to b by romberg's method, the value of the integral is returned as
c ss.

c eps is the relative error acceptable for the estimate.
c jmax is the number of iterations permitted.
c k is the minimum number of iterations(?)
c
c from numerical recipes by press et al. 1988.
c                                 jadc
c-----------------------------------------------------------------------
      implicit none

      integer j,jmax,jmaxp,k

      double precision func,eps,a,b,ss,dss

      parameter (eps=1d-8, jmax = 20, jmaxp = jmax+1, k = 5)

      double precision s(jmaxp), h(jmaxp)
 
      external func
  
      h(1) = 1d0
      do j = 1, jmax
         call trapzd (func,a,b,s(j),j)
         if (j.ge.k) then
            call polint(h,s,j,0d0,ss,dss)
            if (dabs(dss).lt.eps*dabs(ss)) return
         end if
         s(j+1)= s(j)
         h(j+1)= h(j)/4d0
      end do 

      write (*,*) '**error ver410** didnt converge in qromb'
      stop
      end
c
      subroutine polint (h,s,k,x,ss,dss)
c-----------------------------------------------------------------------
c subroutine to extrapolate a polynomial of degree k - 1 to obtain its
c value ss at x, where dss is an estimate of the error in ss.
c the array s contains the values of the function at the corresponding
c values of the array h.
c
c from numerical recipes by press et al. 1988.
c                                 jadc
c-----------------------------------------------------------------------
      implicit none

      integer i,k,ns,m

      double precision s(k),h(k),c(40),d(40),x,ss,dss,dif,dift,ho,hp,w,
     *                 den 

      if (k.gt.40) then
         write (*,*) '**error ver409** ugabugga polint k=',k
         stop
      end if

      ns = 1
      dif = dabs(x-h(1))
      do i = 1, k
         dift = dabs(x-h(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         end if
         c(i)=s(i)
         d(i)=s(i)
      end do 

      ss = s(ns)
      ns = ns - 1
      do m = 1, k-1
         do i = 1, k-m
            ho = h(i) - x
            hp = h(i+m) - x
            w = c(i+1) - d(i)
            den = ho - hp
            if (den.eq.0d0) then
               write (*,*) '**error ver498** polint'
               stop
            end if
            den = w/den
            d(i) = hp * den
            c(i) = ho * den
         end do 
         if (2*ns.lt.k-m) then
            dss = c(ns+1)
         else
            dss = d(ns)
            ns = ns - 1
         end if
         ss = ss + dss

      end do 

      end

      subroutine trapzd (func,a,b,s,j)
c-----------------------------------------------------------------------
c subroutine to compute the j'th stage of refinement of an extended
c trapezoidal rule, returned as s, of the function func between a and b

c from numerical recipes by press et al. 1988.
c                                 jadc
c-----------------------------------------------------------------------
      implicit none

      integer i,j,it

      double precision func,a,b,s,tnm,del,x,sum
 
      external func
 
      if (j.eq.1) then
         s = (b-a)*(func(a)+func(b))/2d0
      else  
         it = j
         tnm = it
         del = (b-a)/tnm
         x = a + del/2d0
         sum = 0d0
 
         do i = 1, it
            sum = sum + func(x)
            x = x + del
         end do 
         s = (s+(b-a)*sum/tnm)/2d0
      end if
 
      end
 
      subroutine hprk
c-----------------------------------------------------------------------
c hprk routine to compute h2o-co2 fugacities from the holland and
c powell (cmp 1991)
c                                 j. connolly 1992.
 
c  for the unitiated input and output is done through common blocks
c  cst5, cst11, and cst26, the relevant variables are:
 
c  for input:
c               pbars= pressure in bars
c               t    = temperature in kelvins
c               r    = gas constant
c               xco2 = mole fraction of co2 in fluid phase
c  for output:
c               v    = molar volume (cm3/mole) at p and t
c               fco2 = the natural log of the co2 fugacity
c-----------------------------------------------------------------------
      implicit none

      double precision rt,p,vco2,vol,xh2o,wco2,wh2o,gco2,gh2o
 
      double precision pbars,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /pbars,t,xco2,u1,u2,tr,pr,r,ps

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk
c----------------------------------------------------------------------
      rt = r*t/1d3
 
      p = pbars/1d3
 
      if (xco2.eq.1d0) then
 
         call crkco2 (pbars,t,vco2,fco2) 
         fh2o = dlog(1d12*p)

      else if (xco2.eq.0d0) then
 
         call crkh2o (pbars,t,vol,fh2o) 
         fco2 = dlog(1d12*p)

      else

         call crkco2 (pbars,t,vco2,fco2) 
         call crkh2o (pbars,t,vol,fh2o) 
  
         xh2o = 1d0 - xco2
 
         wco2 = (13.2d0 - .290d0 * dsqrt(t)) * p**0.25d0
         wh2o = (7.0d0  - 0.15d0 * dsqrt(t)) * p**0.25d0
 
         gco2 = xh2o*xh2o*(wco2+2d0*xco2*(wh2o-wco2))/rt
         gh2o = xco2*xco2*(wh2o+2d0*xh2o*(wco2-wh2o))/rt
 
         fco2 = fco2 + gco2 + dlog(xco2)
         fh2o = fh2o + gh2o + dlog(xh2o)
 
      end if

      end

      subroutine mrk
c---------------------------------------------------------------------
c  routine to calculate properties of H2O-CO2 mixtures using the
c  MRK EoS.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(1), jns(2)

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk
 
      save jns 
      data jns/ 1, 2/
c----------------------------------------------------------------------
      if (xc.eq.1d0) then

         ins(1) = 2
         call mrkpur (ins, 1)
         fh2o = dlog(1d12*p)
         return
        
      else if (xc.eq.0d0) then

         ins(1) = 1
         call mrkpur (ins, 1)
         fco2 = dlog(1d12*p)
         return 

      end if

      call zeroys

      y(2) = xc
      y(1) = 1d0 - xc

      call mrkmix (jns, 2, 1)

      end
 
      subroutine hsmrk
c---------------------------------------------------------------------
c  main program and subprograms written by gary k. jacobs june 1980
c  modified for therm april 1982      jadc
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(1)
 
      double precision bc,bw,xw,t12,rtt,t2,cc,dc,ec,cw,dw,ew,bm,cij,
     *       dij,eij,xc2,xw2,xwc2,cm,dm,em,yzm,fug,rr,hsmrkf

      external fug, hsmrkf

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)
 
      save bw, bc, rr
      data bw, bc, rr/ 29d0, 58d0, 83.144126d0/
c----------------------------------------------------------------------
      if (xc.eq.1d0) then

         ins(1) = 2
         call mrkpur (ins,1)
         fco2 = hsmrkf (v(2),2)
         fh2o = dlog(1d12*p)
         return 
         
      else if (xc.eq.0d0) then

         ins(1) = 1
         call mrkpur (ins,1)
         fh2o = hsmrkf (v(1),1)
         fco2 = dlog(1d12*p)
         return

      end if
 
      xw = 1d0 - xc

      t12 = dsqrt(t)
      rtt = rr*dsqrt(t**3)
      t2 = t*t
 
      cc = 28.31d6+0.10721d6*t-0.00000881d6*t2
      dc = 9380d6-8.53d6*t+0.001189d6*t2
      ec = -368654d6+715.9d6*t+0.1534d6*t2
      cw = 290.78d6-0.30276d6*t+0.00014774d6*t2
      dw = -8374d6+19.437d6*t-0.008148d6*t2
      ew = 76600d6-133.9d6*t+0.1071d6*t2
 
      bm = bc*xc+bw*xw
 
      cij = cc * cw
      dij = dc * dw
      eij = ec * ew
      if (dij.lt.0d0.or.eij.lt.0d0.or.cij.lt.0d0) then
         write (*,1000) p,t
         eij = 0d0
         dij = 0d0
         cij = 0d0
      else
         cij = dsqrt(cij)
         dij = dsqrt(dij)
         eij = dsqrt(eij)
      end if
 
      xc2 = xc*xc
      xw2 = xw*xw
      xwc2 = 2d0*xc*xw
      cm = cc*xc2+cw*xw2+xwc2*cij
      dm = dc*xc2+dw*xw2+xwc2*dij
      em = ec*xc2+ew*xw2+xwc2*eij
c                                 get initial volume estimate for
c                                 newrap from mrk, implicit in newrap
 
c                                 solve for hsmrk volume
      call newrap (bm,cm,dm,em,yzm,t12,rr)
c                                 calculate hsmrk (log) fugacities:
      fco2 = dlog(xc*p)+ fug(rtt,cij,dij,eij,xc,xw,bm,yzm,cm,
     *                      dm,em,bc,cc,dc,ec)
      fh2o = dlog(xw*p)+ fug(rtt,cij,dij,eij,xw,xc,bm,yzm,cm,
     *                      dm,em,bw,cw,dw,ew)
 
1000  format ('**warning ver678** p,t (',g9.3,1x,g9.3,
     *        ') conditions are out of range for HSMRK',/,
     *        'your results may be incorrect.')
 
      end
 
      subroutine newrap (b,c,d,e,yz,t12,r)
c----------------------------------------------------------------------
      implicit none

      integer k

      double precision x,y,bi,bi2,vi2,vi3,x3,y3,pn,pa1,d1,d3,df,yz,
     *                 diff,v,b,c,d,e,t12,r

      double precision p,t,xc,u1,u2,tr,pr,r0,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r0,ps

      double precision vi
      common/ cst26 /vi
c----------------------------------------------------------------------
c                                 call mrk to get initial volume
      call mrk
 
      do k = 1, 50
 
      y = b/4d0/vi
      x = 1d0 - y
 
      bi = vi + b
      bi2 = bi * bi
      vi2 = vi * vi
      vi3 = vi2 * vi
      x3 = x**3
      y3 = y**3
 
      pn = 1d0 + y + y*y - y3
      pa1 = c + d/vi + e/vi2
 
      d1 = -.75d0 * b/vi3/x/x3  - 1d0/vi2/x3
      d3 = (-b/4d0/vi2 - 2d0 * b * b/16d0/vi3
     *     + 0.046875d0 * b**3/vi/vi3 )/ vi/x3
 
      df = ( pn * d1 + d3 )*r*t
     *      - (pa1 * (-1d0/vi/bi2 - 1d0/vi2/bi)
     *      + ( -d/vi2 - 2d0 * e/vi3)/ vi/bi )/t12
 
      v = vi - (pn/vi/x3 * r * t - pa1/t12/vi/bi - p)/df
      diff = dabs(v-vi)
      vi = v
 
      if (diff.lt.0.01d0) exit

      end do 
 
      yz = vi * p/83.14d0/t
 
      end
 
      subroutine qrkmrk
c---------------------------------------------------------------------
c qrkmrk calculates properties of a binary H2O-CO2 fluids using MRK
c acivities and pure fluid routines specified in hybeos. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jns(2)

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xc,u1,u2,tr,pr,rcal,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,rcal,ps

      double precision f
      common/ cst11 /f(3)

      double precision vol
      common/ cst26 /vol

      save jns
      data jns/1, 2/
c----------------------------------------------------------------------- 
c                                 compute pure mrk fluid properties
      call mrkpur (jns,2)
c                                 compute hybrid pure fluid props
      call hybeos (jns,2)

      if (xc.eq.1d0) then 
         f(1) = dlog(1d12*p)
         return
      else if (xc.eq.0d0) then 
         f(2) = dlog(1d12*p)
         return
      end if 

      call zeroys

      y(2) = xc
      y(1) = 1d0 - xc

      call mrkhyb (jns, jns, 2, 2, 1)

      vol = vol + y(1)*vh(1) + y(2)*vh(2)

      end
 
      subroutine haar (vh2o,fh2o)
c-----------------------------------------------------------------------
c     this code is modified from a source file obtained from
c     Bern University and written by Christian DeCapitani at ubc.
 
c     for a documented and more general program see
c     post 1984 steam tables.
 
c haar routine to calculate molar volume and fugacity of h2o from the
c equation of haar et al (steam tables).
 
c                                 jadc 1989.
 
c  for the unitiated input and output is done through common blocks
c  cst5, cst11, and cst26, the relevant variables are:
 
c  for input:
c               p    = pressure in bars
c               t    = temperature in kelvins
c               rref = gas constant
c               xco2 = mole fraction of co2 in fluid phase
c  for output:
c               vh2o = molar volume (cm3/mole) at p and t
c               fh2o = the natural log of the h2o fugacity
c-----------------------------------------------------------------------
      implicit none

      double precision r,t0,amh2o,b,bb,ps,rhn,rh,rh2,y,er,y3,aly,rt,
     *                 bety,f1,f2,pr,dpr,s,del,rhoi2,tau,abc,q10,qm,
     *                 x,dp,dr,gh2o,aid,gid,aa, psat2, vh2o, fh2o

      external psat2

      integer i,loo,nlow,nhigh

      double precision p,t,v,tref,pref,rref,v4
      common/ cst5 /p,t,v(3),tref,pref,rref,v4
 
      integer ki(40), li(40)

      double precision taui(0:6), ermi(0:9), gi(40)
      double precision rhoi(37:40),ttti(37:40),alpi(37:40),beti(37:40)
 
c -----gi are in (bar cc / g)  =  10 * (j / g)
 
      save ki, li, gi, rhoi, ttti, alpi, beti, r, t0, amh2o
 
      data gi /-.53062968529023d4, .22744901424408d5, .78779333020687d4,
     1     -.69830527374994d3, .17863832875422d6, -.39514731563338d6,
     2     .33803884280753d6, -.13855050202703d6, -.25637436613260d7,
     3     .48212575981415d7, -.34183016969660d7, .12223156417448d7,
     4     .11797433655832d8, -.21734810110373d8, .10829952168620d8,
     5     -.25441998064049d7, -.31377774947767d8, .52911910757704d8,
     6     -.13802577177877d8, -.25109914369001d7, .46561826115608d8,
     7     -.72752773275387d8, .41774246148294d7, .14016358244614d8,
     8     -.31555231392127d8, .47929666384584d8, .40912664781209d7,
     9     -.13626369388386d8, .69625220862664d7, -.10834900096447d8,
     *     -.22722827401688d7, .38365486000660d7, .68833257944332d5,
     1     .21757245522644d6, -.26627944829770d5, -.70730418082074d6,
     2     -.225d1, -1.68d1, .055d1, -93.0d1/
 
      data ki /4*1, 4*2, 4*3, 4*4, 4*5, 4*6, 4*7, 4*9, 2*3, 1, 5, 3*2,
     1     4/
 
      data li /1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4,
     1     6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 0, 3*3, 0, 2, 0, 0/
 
      data rhoi /0.319d0, 0.310d0, 0.310d0, 1.55d0/
      data ttti /640d0, 640d0, 641.6d0, 270d0/
      data alpi /34d0, 40d0, 30d0, 1050d0/
      data beti /2d4, 2d4, 4d4, 25d0/
      data r, t0, amh2o/4.6152d0,647.073d0,18.0152d0/
c----------------------------------------------------------------------
      rt = r * t
      nlow = 40
      nhigh = 20
      if (t .lt. 449.35d0) nhigh = 40
 
c -----the values (t/t0)**i are stored in the array taui(i)
 
      taui(0) = 1d0
      taui(1) = t/t0
 
      do i = 2, 6
         taui(i) = taui(i - 1) * taui(1)
      end do 
 
      b = -0.3540782d0 * dlog(taui(1)) + 0.7478629d0 +0.007159876d0 /
     1    taui(3) - 0.003528426d0/taui(5)
      bb = 1.1278334d0 - 0.5944001d0/taui(1) - 5.010996d0/taui(2) +
     1     0.63684256d0/taui(4)
c
      ps = 220.55d0

      if (t .le. 647.25d0) then 
         ps = psat2(t)
         if (p.gt.ps) vh2o = 18d0
      end if 

      rhn = amh2o/vh2o 
 
c -----find the true(?) rh(t,p)
c -----note: pr = pressure corresponding to guessed rh
c            dpr = (dp/drh)
c            the values (1-exp(-rh))**i are stored in the array ermi(i)
 
      do loo = 1, 100
        rh = rhn
        if (rh .le. 0d0) rh = 1d-8
        if (rh .gt. 1.9d0) rh = 1.9d0
        rh2 = rh** 2
        y = rh * b/4d0
        er = dexp(-rh)
        y3 = (1d0 - y)** 3
        aly = 11d0 * y
        bety = 44.33333333333333d0 * y**2
        f1 = (1d0 + aly + bety)/y3
        f2 = 4d0 * y * (bb/b - 3.5d0)
        ermi(0) = 1d0
        ermi(1) = 1d0 - er
        do i = 2, 9
           ermi(i) = ermi(i - 1) * ermi(1)
        end do 
  
        pr = 0d0
        dpr = 0d0
        do i = 1, 36
          s = gi(i)/taui(li(i)) * ermi(ki(i) - 1)
          pr = pr + s
          dpr = dpr + (2d0 + rh*(ki(i)*er - 1d0)/ermi(1)) * s
        end do 

        do i = nlow, nhigh
          del = rh/rhoi(i) - 1d0
          rhoi2 = rhoi(i) * rhoi(i)
          tau = t/ttti(i) - 1d0
          abc = -alpi(i) * del** ki(i) - beti(i) * tau** 2
          if (abc .gt. - 1d2) then
            q10 = gi(i) * del** li(i) * dexp(abc)
          else
            q10 = 0d0
          end if
          qm = li(i)/del - ki(i) * alpi(i) * del** (ki(i) - 1)
          s = q10 * qm * rh2/rhoi(i)
          pr = pr + s
          dpr = dpr + s * (2d0/rh + qm/rhoi(i)) - rh2/rhoi2 * q10 *
     1    (li(i)/del/del + ki(i)*(ki(i) - 1)*alpi(i)*del**(ki(i) - 2))
        end do 
        pr = rh * (rh*er*pr + rt*(f1 + f2))
        dpr = rh * er * dpr + rt * ((1d0 + 2d0*aly + 3d0*bety)/y3 
     *        + 3d0*y*f1/(1d0 - y) + 2d0*f2)
 
        if (dpr .le. 0d0) then
          if (p .le. ps) then
            rhn = rhn * 0.95d0
          else
            rhn = rhn * 1.05d0
          end if
        else
          if (dpr .lt. 1d-2) dpr = 1d-2
          x = (p - pr)/dpr
          if (dabs(x) .gt. 0.1d0) x = 1d-1 * x/dabs(x)
          rhn = rh + x
        end if
        dp = dabs(1d0 - pr/p)
        dr = dabs(1d0 - rhn/rh)
        if (dp .lt. 5d-2 .and. dr .lt. 5d-2) go to 60
      end do 
   60 rh = rhn
 
      y = rh * b/4d0
      x = 1d0 - y
      er = dexp(-rh)
      ermi(0) = 1d0 
      ermi(1) = 1d0 - er
 
      do i = 2, 9
         ermi(i) = ermi(i - 1) * ermi(1)
      end do 
 
c -----calculate base function
 
      aa = rt * (-dlog(x) - 43.33333333333333d0/x 
     *           + 28.16666666666667d0/x/x 
     *           + 4d0*y*(bb/b - 3.5d0) + 15.16666666666667d0
     *           + dlog(rh*rt/1.01325d0) )
 
c -----calculate residual function
 
      do i = 1, 36
         aa = aa + gi(i)/ki(i)/taui(li(i)) * ermi(ki(i))
      end do 
 
      do i = nlow, nhigh
        del = rh/rhoi(i) - 1d0
        tau = t/ttti(i) - 1d0
        abc = -alpi(i) * del** ki(i) - beti(i) * tau**2
        if (abc .gt. - 1d2) then
          aa = aa + gi(i) * del** li(i) * dexp(abc)
        else
          aa = aa
        end if
      end do 
 
      call aideal (t/1d2,rt,aid)
 
      aa = aa + aid
 
      gid = (aid  * amh2o * 1d-1) + rref * t
      gh2o = (aa + p/rh) * amh2o * 1d-1
      fh2o = (gh2o - gid)/rref/t
 
c                           to correct to 1 bar 298 gf +
c                           sliding scale in t ideal gas
c                           gref should be 182051 j / mol
c
      vh2o = amh2o/rh
 
      end

c----------------------------------------------------------------
      double precision function psat2 (t)

      implicit none
 
      double precision t, ps, a(8), w, wsq, v, ff

      integer i
 
      save a
 
      data a /-7.8889166d0, 2.5514255d0, -6.716169d0, 33.239495d0,
     1     -105.38479d0, 174.35319d0, -148.39348d0, 48.631602d0/
 
      if (t .le. 314.00d0) then
        ps = dexp(6.3573118d0 - 8858.843d0/t + 607.56335d0/(t**0.6d0))
      else
         v = t/647.25d0 
         w = dabs(1d0 - v)
         wsq = dsqrt(w)
         ff = 0d0 
         do i = 1, 8
            ff = ff + a(i) * w
            w = w * wsq
         end do 
         ps = 220.93d0 * dexp(ff/v)
      end if

      psat2 = ps 
 
      end
c----------------------------------------------------------------
      subroutine aideal (tr,rt,aid)
 
c subroutine to compute the helmoltz free energy of ideal steam.
 
      implicit none

      integer i
 
      double precision ci(18),tr,rt,aid,w
 
      save ci
 
      data ci /.19730271018d2, .209662681977d2, -.483429455355d0,
     1         .605743189245d1, 22.56023885d0, -9.87532442d0,
     2         -.43135538513d1, .458155781d0, -.47754901883d-1,
     3         .41238460633d-2, -.27929052852d-3, .14481695261d-4,
     4         -.56473658748d-6, .16200446d-7, -.3303822796d-9,
     5         .451916067368d-11,-.370734122708d-13,.137546068238d-15/
 
      w = tr** (-3)
 
      aid = 1d0 + (ci(1)/tr + ci(2)) * dlog(tr)
 
      do i = 3, 18
         aid = aid + ci(i) * w
         w = w * tr
      end do 
 
      aid = -rt * aid
 
      end

      subroutine cohngr (fo2)
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in graphite
c saturated COHN fluid at oxygen fugacity specified 
c by ibuf in routine fo2buf.

c fo2  - log of the oxygen fugacity (from subroutine fo2buf)
c elag - the natural log of graphite activity
c p    - pressure, bar
c t    - temperature, K
c gz   - molar N/C ratio

c derivation and data sources in maple work sheet num_cohn.mws

c                                 JADC 2/04

c fo2 changed to ln(fo2), JADC 11/17/10
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ins(7), nit, i

      double precision fo2,t2,t3,x2,x3,d6,d36,d67,df,xt,d678x,x,tx,
     *                 rad,eq9,dxnh3,deq9,dxh2o,sign,c1,c2,c3,c4,c5

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      save ins
      data ins/1,2,3,4,5,10,11/
c----------------------------------------------------------------------
      t2 = t * t
      t3 = t2 * t
      x = gz

      call fo2buf (fo2)
c                                evaluate lnk's and correct for pressure, carbon 
c                                activity and oxygen fugacity
      c1 = p*dexp(-0.1386241656D2+(0.1230903706D5+0.6372383931D-1*p)/t
     *            -0.8793147005D6/t2+0.7754138439D8/t3+elag)
      c2 = dexp(0.4078341613D-1+(0.47681676177D5+0.6372383931D-1*p)/t
     *         -0.1346621904D6/t2+0.1701579431D8/t3 + elag + fo2)/p
      c3 = dexp(0.1032730663D2+(0.140627396777D5+0.6372383931D-1*p)/t
     *         -0.3712371571D6/t2+0.5351536595D8/t3 + elag + fo2/2d0)/p
      c4 = dexp(-0.7028214449D1+0.3060734044D5/t-0.4750344632D6/t2+
     *          0.5087984255D8/t3 + fo2/2d0)
      c5 = dexp(0.2527543051D8/t3-0.4017985659D6/t2+0.7323735697D4/t
     *         -0.1439146998D2)*p**2
c                                get pure species fugacities
      call mrkpur (ins, 7)
c                                check for graphite saturation:
      y(ins(2)) = c2/g(ins(2))
      y(ins(3)) = c3/g(ins(3))

      if (y(ins(2))+y(ins(3)).ge.1d0) then 
 
         write (*,1000) fo2,p,t

         if (hu.eq.0) then 

            fco2 = dlog(g(ins(2))*p*y(ins(2)))

         else
c                                 return f(h2) in fh2o, f(O2) in fco2
c                                 for carbon projection. 
            fco2 = fo2

         end if

         call zeroys

         y(ins(2)) = 1d0

         return 

      end if 
c                                since ynh3 is quadratic in yh2o we
c                                have two roots to choose from, differentiated
c                                by the sign on the radical "rad", the 
c                                first root seems to be relevant for 
c                                high f(O2).
      sign = -1d0 

      do i = 1, 2 
c                                guess xh2o, this might not be a bad guess for
c                                oxidixed fluids, but for reduced?
         y(ins(1)) = 1d0 - y(ins(3)) - y(ins(2)) 
         nit = 0 

         do 
c                                 the problem in this loop is that
c                                 the equation has 4 roots, when we
c                                 have the right root for xnh3 then
c                                 there is only one real root between
c                                 zero and one, if we have the wrong
c                                 root for xnh3 we have two real roots
c                                 between zero and one, but these give
c                                 negative values for nh3, and then there
c                                 are the imaginary roots that i don't 
c                                 even want to think of.
            x2    = y(ins(1))*y(ins(1)) 
            x3    = y(ins(1)) * x2  

            d6    = c4*g(ins(5))/g(ins(1))
            d36   = c1/g(ins(4))/c4**2*g(ins(1))**2
            d67   = c4**3/g(ins(1))**3/c5*g(ins(7))**2/g(ins(6))
            xt    = y(ins(2)) + y(ins(3))
            df    = (1d0+d6)/d6
            d678x = d67*8d0*x 
            rad   = y(ins(1))*(x3+d678x*(d36*x2+xt))
            if (rad.lt.0d0) then 
               sign = -sign
               exit
            end if 
c                                quadratic root for nh3:
            rad   = sign*dsqrt(rad)
            y(ins(7))  = (-x2+rad)*y(ins(1))/4d0/d67
c                                a negative value for xnh3 is ominous
c                                but does it really mean the end of this
c                                root?  
            if (y(ins(7)).lt.0d0) then 
               sign = -sign
               exit
            end if 
       
            eq9   = 1d0-(d36*y(ins(1))+df)*y(ins(1))-xt
     *                   -y(ins(7))*(1d0+y(ins(7))*d67/x3)
            dxnh3 = ((-3d0*y(ins(1))+1d0/rad*
     *                ((4d0*y(ins(1))+3d0*d678x*d36)*x2
     *                +d678x*xt)/2d0)*y(ins(1))+rad)/d67/4d0
            deq9  = -2d0*d36*y(ins(1))-df-dxnh3
     *            + (3d0*y(ins(7))/y(ins(1))-2d0*dxnh3)*d67*y(ins(7))/x3

            dxh2o = eq9/deq9

            y(ins(2))  = c2/g(ins(2))
            y(ins(3))   = c3/g(ins(3))
            y(ins(5))  = y(ins(1))/d6
            y(ins(4)) = d36*x2
            y(ins(6))  = y(ins(7))**2*d67/x3

            tx = xt + y(ins(5)) + y(ins(1)) + y(ins(4)) + y(ins(7)) 
     *              + y(ins(6)) - 1d0

            nit = nit + 1

            if (nit.gt.iopt(21)) then
c                                 not converging to much
c                                 of anything, try other root
               write (*,1000) t,p

               sign = -sign

               exit

            else if (dabs(dxh2o).lt.nopt(50).and.
     *               dabs(tx).lt.nopt(50)) then 
c                                 seems to have converged
               if (y(ins(1)).gt.1d0.or.y(ins(1)).lt.0d0.or.
     *             y(ins(7)).gt.1d0.or.y(ins(7)).lt.0d0) then
c                                 but a bad root
                  sign = -sign
                  exit
               else
c                                 everything seems ok
                  goto 90
               end if 
            end if 
c                                 get new gamma's
            call mrkmix (ins, 7, 1)
            y(ins(1)) = y(ins(1)) - dxh2o

         end do 
      end do 
c                                 if we get here, no good solution
      write (*,*) 'fd'
      stop

90    if (hu.eq.0) then
 
         fh2o = dlog(g(ins(1))*p*y(ins(1)))
         fco2 = dlog(g(ins(2))*p*y(ins(2)))

      else
c                                 return f(h2) in fh2o, f(O2) in fco2
c                                 for carbon projection. 
         fh2o = dlog(g(ins(5))*p*y(ins(5)))
         fco2 = fo2

      end if 

1000  format (/,'**warning ver222** routine COHNGR, specified lnfO2 (',
     *        g12.6,')',/,'is inconsistent with graphite saturation',
     *       ' at P(bar)=',g12.6,' T(K)=',g12.6,/,'XCO2=1 assumed.',/) 
      end

      subroutine hcneos (gex,x3,x1,x2)
c----------------------------------------------------------------------
c subroutine to compute excess energy (gex) of a H2O-CO2-NaCl
c mixture. x3 = mol fraction NaCl; x2 = mol fraction CO2; 
c x1 = mol fraction H2O. This routine is for the true ternary
c fluid. see routine waddah for h2o-co2 with constant mole
c or weight fraction of salt.
c----------------------------------------------------------------------
      implicit none 

      double precision pbar,t,xco2,u1,u2,tr,pr,r,ps,alpha,
     *                 w2,w3,w4,w5,a1,gex,x1,x2,x3,fh2o,fco2,
     *                 p,v1,v2,rt,rid,sid,tid,uid,vid

      common/ cst5  /pbar,t,xco2,u1,u2,tr,pr,r,ps

c                                 these routines return volume cm3/mol
      call crkco2 (pbar,t,v2,fco2) 
      call crkh2o (pbar,t,v1,fh2o)

      rt = r*t

      p = pbar/1d3

      w2 = 906.12d0 - 57.277d0 * p
      w3 = 101788d0 - 2916d0 * p
      w4 = 38007d0  + 2445d0 * p
      w5 = -37371d0 + 916d0 * p

      alpha = dexp(4.04d0 - 0.1611d0*v1) - 134.2d0 * p/t

      if (alpha.lt.0d0) then 
         alpha = 0d0
      else if (alpha.gt.1d0) then 
         alpha = 1d0
      end if 
         
      a1 = 1d0 + alpha

      sid = 0d0
      tid = sid
      uid = sid
      vid = sid

      if (x1.gt.1d-8) sid = x1*dlog(x1)
      if (x2.gt.1d-8) sid = sid + x2*dlog(x2)  
      if (x3.gt.1d-8) then
         sid = sid + x3*dlog(x3)
         rid = x3/(x1+x3)
         vid = x3*(    a1 * dlog(a1/(1d0+alpha*rid))
     *            + alpha * dlog(rid))
     *         -x1*dlog(1d0+alpha*rid)
      end if 

      if (x2+x3.gt.1d-8) uid = (x2*w3+x3*w4)/(x2+x3)
      if (x2+x1.gt.1d-8) tid = 202046.4d0*(x1+x2)/(v1*x1+v2*x2)

      gex = rt*(sid+vid)+x2*(x1*tid+x3*(uid+x1*w5))+x1*x3*w2

      end 

      subroutine waddah 
c----------------------------------------------------------------------
c subroutine to compute fh2o and fco2 for a H2O-CO2-NaCl mixture. 
c using cork (holland and powell '98) for pure h2o and co2 properties
c and leonid aranovich's van-laar model for activities.

c  x3 = mol fraction NaCl
c  x2 = mol fraction CO2 
c  x1 = mol fraction H2O 

c  elag - the salt fraction
c  ibuf - if 1, salt fraction is mass fraction
c         else if 2, salt fraction is mole fraction

c  xco2 - nco2/(nco2+nh2o)
c  fco2 - natural log of co2 fugacity
c  fh2o - natural log of h2o fugacity
c                                 JADC, 4/27/04.
c----------------------------------------------------------------------
      implicit none 

      double precision pbar,t,xco2,u1,u2,tr,pr,r,ps,alpha,w1,
     *                 w2,w3,w4,w5,x1,x2,x3,t0,t1,p,v1,v2,rt,xt,
     *                 wh2o,wco2,wnacl,nh2o,nco2,nnacl,ntot

      common/ cst5  /pbar,t,xco2,u1,u2,tr,pr,r,ps

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx
c                                 molar weights
      save wh2o,wco2,wnacl
      data wh2o,wco2,wnacl/18.016,44.01,58.446/
c                                 convert xco2 and elag to
c                                 mole appropriate mole fractions:
      if (ibuf.eq.1) then
c                                 from weight salt:
         if (xco2.eq.1d0) then 
            nh2o = 0d0
            nco2 = 1d0 
            nnacl = -elag*wco2/wnacl/(elag-1d0)
         else if (xco2.eq.0d0) then 
            nh2o = 1d0
            nco2 = 0d0 
            nnacl = -elag*wh2o/wnacl/(elag-1d0)
         else 
            nh2o = (elag-1d0)*(xco2-1d0)/(xco2*(wco2-wh2o)+wh2o)
            nco2 = xco2*nh2o/(1d0-xco2)
            nnacl = elag/wnacl
         end if 

         ntot = nh2o + nco2 + nnacl
         x1 = nh2o/ntot
         x2 = nco2/ntot
         x3 = 1d0 - x1 - x2

      else if (ibuf.eq.2) then 
c                                 from salt mole fraction:
         x3 = elag
         xt = 1d0 - x3
         x2 = xco2*xt
         x1 = 1d0 - x2 - x3

      else
c                                 program or input error?
         call error (999,t,ibuf,'WADDAH')
      end if
c                                 supposedly these routines 
c                                 return volume in J/bar
      call crkco2 (pbar,t,v2,fco2) 
      call crkh2o (pbar,t,v1,fh2o)
       
      if (x1.eq.1d0.or.x2.eq.1d0.or.x3.eq.1d0) return

      rt = r*t

      p = pbar/1d3

      w1 = 202046.4d0
      w2 = 906.12d0 - 57.277d0 * p
      w3 = 101788d0 - 2916d0 * p
      w4 = 38007d0  +  2445d0 * p
      w5 = -37371d0 + 916d0 * p

      alpha = dexp(4.04d0 - 0.1611d0*v1) - 134.2d0 * p/t
c                                 modified 6/6/05 to restict
c                                 alpha to physical values. JADC
      if (alpha.lt.0d0) then 
         alpha = 0d0
      else if (alpha.gt.1d0) then 
         alpha = 1d0
      end if 

      t0 = (v1*x1+v2*x2)**2
      t1 = x2 + x3
c                                 this is ln(fH2O)
      if (x1.ne.0d0) then 
         fh2o = fh2o + ( w2*x3*t1 - w5*x2*(x1-x2-x3)*x3
     *        - x2*x3*(w3*x2+w4*x3)/t1
     *        + w1*x2*(v1*x1**2*x3+v2*x2*(x1+x2+x1*x3))/t0)/rt
     *        + dlog(x1*(x1+x3)/(1d0+x3*alpha))
      else 
         fh2o = dlog(p*1d4)
      end if 
c                                 this is ln(fCO2)
      if (x2.ne.0d0) then 
         fco2 = fco2 + ( w5*x1*x3*(x1-x2+x3) - w2*x1*x3 
     *        + w1*x1*(v2*x2**2*x3+v1*x1*(x1+x2+x2*x3))/t0
     *        + x3/t1**2*(w4*x3*(-x2**2+x1*x3+x3**2)
     *        + w3*x2*(2d0*x3*t1+x1*(t1+x3))))/rt + dlog(x2)
      else 
         fco2 = dlog(p*1d4)
      end if 

      end 

      subroutine dimond (agph)
c-----------------------------------------------------------------------
c dimond tests if p-t conditions are in the diamond stability field
c if they are it computes the activity of graphite needed to represent
c diamond stability for C-O-H fluid speciation routines. 


c  elag - natural log of the activity of the stable C polymorph
c         requested by the user (usually 0).
        
c  agph - natural log of graphite activity, this is adjusted
c         to account for deviations imposed by the user via
c         the input variable elag.

c polynomial functions fit by S. Poli, 4/27/04
c                                 JADC, 4/27/04
c-----------------------------------------------------------------------
      implicit none 

      double precision ptrans,agph

      integer ibuf,hu,hv,hw,hx      
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

c                                 get transition pressure
c                                 (should really use the G
c                                 function for consistency):

      ptrans = 5284.165053d0 + (33.21515773d0 - .002106330992d0*t)*t

      if (p.lt.ptrans) then 
         agph = elag
      else 
c                                 compute corrected graphite activity
         agph = elag + 0.008423508384179629d0 
     *        + (4.693008650307614d-11 * p - 3.850380793502567d-5)*p
     *        + (1.4126916053951515d-3 + 1.393226795939807d-8*p 
     *        - 5.887505938975768d-7 * t)*t

      end if 

      end 

      subroutine pseos (v,f,jam)
c-----------------------------------------------------------------------
c compute ln(f[iam], bar) and volume (cm3/mol) for pure fluids from Pitzer 
c & Sterner JCP '94. 

c      p   - pressure (bars)
c      t   - temp (K)
c      iam - 1 for H2O, 2 for CO2.

c                                 J.A.D. Connolly, Nov 26, 2010.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision r,prt,rt,f,c1,c2,c3,c4,c5,c6,c7,c8,c9,c0,v,a1,
     *                 a2,a3,c12,c20,c33,c34,c36,c44,c46,c55,c56,c66,
     *                 c64,c53,c42,e1,e2,dv,t2

      integer it, jam, iwarn

      double precision p,t,xco2,u1,u2,tr,pr,rc,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,rc,ps

      save r, iwarn
      data r, iwarn/83.14d0, 0/
c----------------------------------------------------------------------
      t2 = t*t
c                                 temperature dependent coefficients
      if (jam.eq.1) then 
c                                 h2o 
         c1 = 0.24657688d6/t + 0.51359951d2
         c2 = 0.58638965/t -0.28646939d-2 + 0.31375577d-4*t
         c3 = -0.62783840d1/t + 0.14791599d-1 
     *                      + t*(0.35779579d-3 + 0.15432925d-7*t)
         c4 = -0.42719875 - 0.16325155d-4*t
         c5 = 0.56654978d4/t - 0.16580167d2 + 0.76560762d-1*t
         c6 = 0.10917883
         c7 = ((0.38878656d13/t**2 - 0.13494878d9)/t + 0.30916564d6)/t 
     *                     + 0.75591105d1
         c8 = -0.65537898d5/t + 0.18810675d3
         c9 = ((-0.14182435d14/t**2 + 0.18165390e9)/t - 0.19769068e6)/t 
     *                     - 0.23530318e2
         c0 = 0.92093375d5/t + 0.12246777d3
c                                 CORK volume guess and backup fugacity
         call crkh2o (p,t,v,f) 

      else if (jam.eq.2) then 
c                                 co2 
         c1 = 0.18261340d7/t + 0.79224365d2
         c2 = 0.66560660d-4 + 0.57152798d-5*t + 0.30222363d-9*t2
         c3 = 0.59957845d-2 + 0.71669631d-4*t + 0.62416103d-8*t2
         c4 = -0.13270279d1/t - 0.15210731d0 
     *                         + 0.53654244d-3*t - 0.71115142d-7*t2
         c5 = 0.12456776/t + 0.49045367d1 
     *                         + 0.98220560d-2*t + 0.55962121d-5*t2
         c6 = 0.75522299d0
         c7 = ((-0.39344644d12/t2 +  0.90918237d8)/t + 0.42776716d6)/t 
     *         - 0.22347856d2
         c8 = 0.40282608d3/t + 0.11971627d3
         c9 = (0.22995650d8/t -0.78971817d5)/t -0.63376456d2
         c0 = 0.95029765d5/t + 0.18038071d2
c                                 CORK volume guess and backup fugacity
         call crkco2 (p,t,v,f) 

      else 

         call error (11,xco2,jam,'species (routine pseos)')

      end if 

      c12 = 12d0*c5
      c20 = 20d0*c6
      c46 = 6d0*c4
      c53 = 3d0*c5
      c42 = 2d0*c4
      c64 = 4d0*c6
      c33 = 2d0*c3*c3
      c34 = 8d0*c3*c4
      c36 = -16d0*c3*c6 - c12*c42
      c44 = 8d0*c4*c4 + c3*c12
      c55 = -32d0*c4*c6 - 18d0*c5*c5
      c56 = -c12*c64
      c66 = 32d0*c6*c6

      rt = r*t
      prt = p/rt
      it = 0
c                                 iteration loop for volume
      do

         a1 = c2 + (c3 + (c4 + (c5 + c6/v)/v)/v)/v
         a2 = a1*a1
         a3 = a2*a1
         e1 = c7*dexp(-c8/v)
         e2 = c9*dexp(-c0/v)
c                                 refine guess by dp/(dp/dv)
         dv = 
c                                 (p - p(v))/rt
     *         (prt - (1d0 + (c1 + e1 + e2)/v 
     *          - (c3 + (c42 + (c53 + c64/v)/v)/v)/v/a2)/v)/
c                                 dpdv/rt
     *         ((-1d0 + (2d0*(c3/a2-c1-e1-e2) + (c8*e1 + c0*e2 
     *          + (c46*a1 - c33)/a3 + (a1*c12 - c34 
     *          + (c20*a1 - c44 + (c36 + (c55 + (c56 
     *          - c66/v)/v)/v)/v)/v)/v/a3)/v)/v)/v/v)

         if (dv.lt.0d0.and.v+dv.lt.0d0) then 

            v = v*0.8d0

         else 

            v = v + dv

         end if 

         if (dabs(dv/v).lt.nopt(50)) then
c                                 converged, compute ln(fugacity)
            f = c1/v+1d0/a1-1d0/c2-(e1-c7)/c8-(e2-c9)/c0
     *          + dlog(rt/v) + p*v/rt - 1d0

            exit
          
         else if (v.lt.0d0.or.it.gt.iopt(21)) then
c                                 will use cork fugacities
            iwarn = iwarn + 1

            if (iwarn.le.50) then 
               write (*,1000) p,t,v
               if (iwarn.eq.50) call warn (49,p,93,'PSEOS')
            end if 

            exit 

         end if 

         it = it + 1

      end do 

1000  format (/,'**warning ver093** PSEoS did not converge at:',
     *        3(1x,g12.6))

      end 

      subroutine pshp
c-----------------------------------------------------------------------
c pshp routine to compute h2o-co2 fugacities from the pitzer & sterner
c (1994) eos, using excess properties from holland & powell (2003).
c                                 JADC, Nov 28 2010.
 
c  for the unitiated input and output is done through common blocks
c  cst5, cst11, and cst26, the relevant variables are:
 
c  for input:
c               pbars= pressure in bars
c               t    = temperature in kelvins
c               r    = gas constant
c               xco2 = mole fraction of co2 in fluid phase
c  for output:
c               fh2o = the natural log of the co2 fugacity
c               fco2 = the natural log of the co2 fugacity
c-----------------------------------------------------------------------
      implicit none

      double precision vh2o,vco2,w,xh2o,whc
 
      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk
                  
      save whc
      data whc/13000d0/ 
c
c 252.586d0 
c----------------------------------------------------------------------
      if (xco2.eq.1d0) then
 
         call pseos (w,fco2,2) 
         fh2o = dlog(1d12*p)

      else if (xco2.eq.0d0) then
 
         call pseos (w,fh2o,1) 
         fco2 = dlog(1d12*p)

      else

         call pseos (vco2,fco2,2) 
         call pseos (vh2o,fh2o,1) 

         xh2o = 1d0 - xco2
c                                 whc = 2*ahc/R, where ahc is the h2o-co2
c                                 interaction parameter estimated by H&P
c                                 2003, 10.5 kJ^2/kbar, the volumes coming
c                                 out of are cm3 = 0.1 J/bar therefore 
c                                 ahc' * vcm3 = (ahc * 0.1) * vcm3  = ahc * vjbar
c                                 ahc' = 10500 * 0.1 => whc = .2*ahc/8.314 = 252.586
         w = whc/t/(xh2o*vh2o+xco2*vco2)**2
c                                 ln(gamma) = w*v*x^2
c                                 verified again, JADC, 11/2021
         fco2 = fco2 + w*vco2*xh2o**2 + dlog(xco2)
         fh2o = fh2o + w*vh2o*xco2**2 + dlog(xh2o)

      end if

      end

      subroutine rksi4 (bad,iavg)
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 4 species silica vapor

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet si-o_rk4_R=R.mws
c                                 JADC 7/12
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ins(4), isp, nit, i1, i2, i3, i4, iroots, ineg, ipos, 
     *        i, icon, iavg, iwarn

      logical bad

      double precision c1,c2,rat,rp1,rm1,r2m1,oldy,a0,a1,a2,
     *                 vmin,vmax,x(3),lnk2,lnk3

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      save ins, isp, i1, i2, i3, i4, iwarn
      data isp, ins, i1, i2, i3, i4, iwarn/4, 14, 13, 12, 7,  
     *                                        14, 13, 12, 7, 0/
c----------------------------------------------------------------------
c                                 rat = nsi/no = xc/(1-xc) 
      rat = xc/(1d0-xc)
c                                 evaluate K's and correct for pressure
c                                 c1 = exp(lnK_1)*p => 2 O = O2, HSC K
      c1 = dexp((-0.9214495D6/t + 0.6234471D5)/t - 0.1631235D2) * p
c                                 c2 = exp(lnK_2)/p => SiO2 = SiO + O, HSC K 
c     c2 = dexp((-0.1133204D7/t - 0.5491882D5)/t + 0.1710990D2) / p
c                                  k_2 from shornikov enthalpy
      lnk2 = (-1.133204d+06/t - 5.491882d+04)/t + 1.710990d+01
      c2 = dexp(lnk2)/p
c                                 some inner loop constants
      rp1 = rat + 1d0
      rm1 = rat - 1d0
      r2m1 = 2d0*rat - 1d0
c                                 get pure species fugacities
      call mrkpur (ins, isp)

      nit = 0 
      oldy = 0d0

      icon = i1 
 
      do 
c                                 solve (yo^3 + a1*yo^2 + a2*yo + a3) for yo
         a0 = g(i4)*c2*g(i1)*rm1/c1/g(i3)**3/g(i2)
         a1 = g(i4)*(r2m1/g(i3)**2 + c2*g(i1)/g(i2)/g(i3)**3)/c1
         a2 = (c2*g(i1)*g(i3)/g(i2)*rp1 - rm1*g(i4)/c1)/g(i3)**2

         call roots3 (a2,a1,a0,x,vmin,vmax,iroots,ineg,ipos)

         do i = 1, iroots
c                                 find an acceptable root, if any... 
            if (x(i).le.0d0.or.x(i).ge.1d0) cycle  
c                                 monatomic O                 
            y(i3) = x(i)
c                                 back calculate remaining fractions:
c                                 K1 => O2: 
            y(i4) = c1/g(i4)*(y(i3)*g(i3))**2
c                                 mass balance => sio2:
            y(i1) = c2*g(i1)*(1d0 - y(i3) - y(i4)) / 
     *                (g(i2)*y(i3)*g(3) + c2*g(i1))

            if (y(i1).lt.0d0) then

               if (dabs(y(i1)).lt.nopt(50)) then
                  y(i1) = 0d0
               else 
                  cycle
               end if 
                     
            else if (y(i1).gt.0.5d0) then 

               icon = i1 

            end if 

            y(i2) = 1d0 - y(i1) - y(i3) - y(i4)

            if (y(i2).lt.0d0) then

               if (dabs(y(i2)).lt.nopt(50)) then
                  y(i2) = 0d0
               else 
                  cycle
               end if 

            else if (y(i2).gt.0.5d0) then
 
               icon = i2 

            end if

            bad = .false. 

            exit 
 
         end do 

         if (dabs(oldy-y(icon)).lt.nopt(50)) exit  
c                                 get new gamma's
         call mrkmix (ins, isp, iavg)

         oldy = y(icon)

         nit = nit + 1

         if (nit.lt.iopt(21)) cycle
 
         bad = .true.
         exit 

      end do 

      if (bad) then 

         if (nit.gt.iopt(21).and.iwarn.lt.100) then

             write (*,'(a,2(g12.6,1x))') 
     *            'ugga rksi4 not converging T,P:',t,p

         else if (iwarn.lt.100) then

             write (*,'(a,5(g12.6,1x))') 
     *            'ugga rksi4 not valid solution T,P:',t,p,x
             
         end if 

         iwarn = iwarn + 1
    
         if (iwarn.eq.100) call warn (49,t,0,'RKSI4')

         call setbad (fh2o)

         return 

      end if 
c                                lnK_3 => SiO = Si + O, shornikov H_SiO 
      lnk3 = (1.906315d6/t - 1.005993d5)/t + 1.664069d1
      fh2o = dlog(p*g(i3)*y(i3))

      if (y(i2).gt.0d0) then 

         fco2 = lnk3 + dlog(g(i2)*y(i2)/g(i3)/y(i3)) 

      else       
c                                  K_2 => SiO2 = SiO + O, HSC K 
c     c2 = dexp((-0.1133204D7/t - 0.5491882D5)/t + 0.1710990D2) / p
c                                  k_2 from shornikov enthalpy
         lnk2 = (-1.133204d+06/t - 5.491882d+04)/t + 1.710990d+01
         fco2 = lnk3 + lnk2 + dlog(g(i1)*y(i1)/(g(i3)*y(i3))**2) 

      end if 

      end

      subroutine rksi5 (c1, c2, c3, bad)
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 5 species silica vapor

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet Si-O_rk5_R=R_speciation.mws
c                                 JADC 6/12
c the Si-O2_rk4_v1_speciation.mws sheet for rksi4a may be more rational.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer iavg, ins(5), isp, nit, i1, i2, i3, i4, i5, i

      logical bad, both

      double precision c1,c2,c3,rat,rp1,rm1,nsi,no,oymin,nymin,oy(nsp),
     *                 r2p1,r2m1,dquart,oymax,nymax,by(nsp), bg(nsp),
     *                 bvol, gold, gnew

      external dquart

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3 

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer iroots
      logical switch, rkmin, min
      common/ rkdivs /iroots,switch,rkmin,min

      double precision vol
      common/ cst26 /vol

      double precision vrt
      integer irt
      logical sroot
      common/ rkroot /vrt,irt,sroot

      save isp, ins, i1, i2, i3, i4, i5, both, iavg
      data isp, ins, i1, i2, i3, i4, i5, both, iavg
     *                                      /5, 14, 13, 12, 7, 15, 
     *                                          14, 13, 12, 7, 15, 
     *                                          .false., 1/
c----------------------------------------------------------------------
c                                 get pure species fugacities
      call mrkpur (ins, isp)
c                                 zero species in case of degenerate compositions
      do i = 1, isp
         y(ins(i)) = 0d0
      end do 
c                                 some inner loop constants
c                                 rat = nsi/no = xc/(1-xc) 
      rat = xc/(1d0-xc)
      rp1    = rat + 1d0
      rm1    = rat - 1d0
      r2p1    = 2d0*rat + 1d0
      r2m1    = 2d0*rat - 1d0

      switch = .false.
      both = .false.
c                                switch root entry point
10    nit = 0 
      oymin = 1d0
      oymax = 0d0 
      bad = .false. 
c                                 initial guess:
      if (rat.ge.0.5.and.rat.le.1d0) then 
c                                 SiO-SiO2 mix 
         y(14) = (1d0 - rat)/rat
         y(13) = 1d0 - y(14)
         y(15) = nopt(50)
         y(7)  = nopt(50)/4d0
         y(12) = nopt(50)/2d0

      else
c                                 Si-SiO mix
         y(14) = nopt(50)
         y(13) = 1d0/rat
         y(15) = 1d0 - 1/rat
         y(7)  = nopt(50)
         y(12) = nopt(50)

      end if 

      call mrkmix (ins, isp, iavg)    

      do 
c                                 solve mass balance: yo^4 + a3*yo^3 + a2*yo^2 + a1*yo + a0) 
c                                 for yo
         a0 = -c2 * c3 / c1 * g(i1) * g(i4) / g(i5) / g(i3) ** 4 
         a1 = c2 * g(i1) * g(i4) * (rm1 * g(i3) / g(i2) + c3 * rp1 
     *            / g(i5)) / g(i3) ** 4 / c1
         a2 = (c2 * g(i1) * g(i3) / g(i5) * c3 * r2p1 + g(i4) * 
     *           (r2m1 * g(i3) + c2 * g(i1) / g(i2)) / c1) / g(i3) ** 3

         a3 = (c2 * g(i1) * g(i3) / g(i2) * rp1 - rm1 * g(i4) / c1) 
     *           / g(i3) ** 2
c                                 monatomic O     
         call newton (dquart,1d0,0d0,1d-12,y(i3),bad)

         if (bad) then 
           
            exit 

         else if (y(i3).eq.0d0) then 
c                                 may not find the root if switch on 
c                                 first iteration       
            y(i3) = nopt(50)

         else if (isnan(y(i3)).or.y(i3).le.0d0.or.y(i3).eq.nopt(50)) 
     *                                                              then

            bad = .true.
            exit 

         end if 
c                                 back calculate remaining fractions:
c                                 K1 => O2: 
         y(i4) = c1/g(i4)*(y(i3)*g(i3))**2
c                                 mass balance => sio: this might be singular
c                                 at R = 1/2?
         y(i2)  = g(i5)*y(i3)*g(i3)*((2d0 - y(i3))*rat 
     *             - 1d0 + y(i3) + y(i4)) / rat 
     *             / (g(i5)*y(i3)*g(i3) + 2d0*c3*g(i2))
c
         if (y(i2).le.0d0) then 

            if (dabs(y(i2)).lt.nopt(50)) then 

               y(i2) = 0d0

            else 

               bad = .true.
               exit 

               do i = 1, isp
                  y(ins(i)) = 0d0
               end do 
c                                 two 4 species bailouts
c                                 4a has Si-O2-SiO-SiO2
c                                 4 has SiO2-SiO-O-O2,
c                                 neither is appropriate for reduced 
c                                 liquids.
               call rksi4a (c1,c2,iavg,bad)

               if (bad) then
 
                  do i = 1, isp
                     y(ins(i)) = 0d0
                  end do 

                  call rksi4 (bad,iavg)
     
               end if 

               return

            end if 

         end if 
c                                  K3 => Si
         y(i5) = c3/g(i5)/y(i3)/g(i3)*y(i2)*g(i2)
c                                 closure => sio2: 
         y(i1) = 1d0 - y(i2) - y(i3) - y(i4) - y(i5)

         if (y(i1).lt.0d0) then

            if (dabs(y(i1)).lt.nopt(50)) then 
               y(i1) = 0d0
            else 
               bad = .true.
               exit
            end if  

         end if 

         nymin = 1d0
         nymax = 0d0 
         no = 0d0      
c                                 renormalize, this helps! 
         do i = 1, isp
            no = no + y(ins(i))
         end do

         do i = 1, isp
            y(ins(i)) = y(ins(i))/no
         end do 

         do i = 1, isp

            if (y(ins(i)).gt.nymax) nymax = y(ins(i))
            if (y(ins(i)).lt.nymin.and.
     *          y(ins(i)).gt.0d0) nymin = y(ins(i))

         end do

         nsi = y(14) + y(13) + y(15)
         no  = 2d0*(y(7)+y(14)) + y(13) + y(12) 

         if (dabs(nymax-oymax).lt.nopt(50).and.
     *        dabs(nymin-oymin).lt.nopt(50).and.
     *        dabs(xc-nsi/(nsi+no)).lt.nopt(50) ) then

            exit 

         else if (nit.gt.iopt(21)) then  
 
            bad = .true.
            exit
 
         end if 
c                                 make change less violent
c                                 this does help!
         if (nit.gt.1) then 

            do i = 1, isp

               y(ins(i)) = (0.5*y(ins(i))+0.5*oy(ins(i)))

            end do  

         end if 
c                                 get new gamma's
            call mrkmix (ins, isp, iavg)           
c                                could converge to the wrong speciation
            if (sroot.and.iroots.eq.3) then 
               both = .true.
            else if ((.not.both).and.iroots.eq.3) then
               both = .true.
               min = rkmin
            end if

         nit = nit + 1
         oymin = nymin
         oymax = nymax
         do i = 1, isp
            oy(ins(i)) = y(ins(i))
         end do 

      end do 

      if (both) then

         if (switch) then
c                                on the second solution
            gnew = 1d9
            gold = 1d9
            if (y(i3).gt.1d-20.and.y(i5).gt.1d-20)  
     *         gnew =  (1d0-xc) * dlog(g(i3)*y(i3)) 
     *                     + xc * dlog(g(i5)*y(i5))
            if (by(i3).gt.1d-20.and.by(i5).gt.1d-20) 
     *         gold =  (1d0-xc) * dlog(bg(i3)*by(i3)) 
     *                     + xc * dlog(bg(i5)*by(i5))

            if (gold.lt.gnew) then 
c                                swap to the old solution
               do i = 1, isp
                  y(ins(i)) = by(ins(i))
                  g(ins(i)) = bg(ins(i))
               end do 

               vol = bvol
c                                for sroot, 1st solution is good
               switch = .false.
   
            end if             

         else 
c                                save old solution
            do i = 1, isp
               by(ins(i)) = y(ins(i))
               oy(ins(i)) = 0d0 
               bg(ins(i)) = g(ins(i))
            end do 

            bvol = vol

            switch = .true.

            call mrkpur (ins, isp)

            goto 10

         end if  

      end if 

      end

      subroutine rkboth 
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 5 species silica vapor

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet Si-O_rk5_R=R_speciation.mws
c                                 JADC 6/12
c the Si-O2_rk4_v1_speciation.mws sheet for rksi4a may be more rational.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical bad

      integer ins(5), isp, i1, i2, i3, i4, i5, i, itic, ibad, igood

      double precision oy(nsp), og(nsp), xfo, xfs, xv, lnk1, 
     *                 n, dg, g1, g2, rat, c1, c2, c3, lnk2, lnk3

c      double precision go2, go, gsi, gsio, gsio2, gzero
c      external gzero

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision vol
      common/ cst26 /vol

      save isp, ins, i1, i2, i3, i4, i5, itic
      data isp, ins, i1, i2, i3, i4, i5, itic 
     *                                  /5, 14, 13, 12, 7, 15, 
     *                                      14, 13, 12, 7, 15, 0/
c----------------------------------------------------------------------
c                                  low T hack
      if (t.lt.2400d-3.and.v(14).lt.1d2.and.
     *     xc.gt.0.326.and.xc.lt.0.340) then
c                                 conditional to destabilize the 
c                                 multispecies fluid at low P in favor
c                                 of SiO2(g)

c                                 a more elegant solution would be to 
c                                 look at the actual solution and reject
c                                 if it's the dense phase at T < ~2400.
         fh2o = dlog(1d12*p)
         fco2 = dlog(1d12*p)
         return

      end if 

      if (xc.eq.0d0.or.xc.eq.1d0) then 
c                                 degenerate cases
c                                 get pure species fugacities, formerly
c                                 in these cases xc was set to tol or 1 - tol
c                                 probably to stabilize optimization? or maybe for
c                                 chemical potentials.
         call mrkpur (ins, isp)
c                                 zero species in case of degenerate compositions
         do i = 1, isp
            y(ins(i)) = 0d0
         end do 

         if (xc.eq.1d0) then 
c                                 pure Si
            fh2o = dlog(1d8*p)
            fco2 = dlog(p*g(i5))
            y(i5) = 1d0 

         else 
c                                 pure O
c                                 evaluate K's and correct for pressure
c                                 c1 = exp(lnK_1)*p => 2 O = O2, HSC K
c            go =  gzero (2)
c            go2 = gzero (ipoint+2)

            lnk1 = (-0.9214495D6/t + 0.6234471D5)/t - 0.1631235D2
c            lnk1 = (2d0*go - go2)/r/t

            c1 = dexp(lnk1) * p

            if (c1.gt.1d0/nopt(50)) then 
c                                assume pure O2
               fh2o = (dlog(p*g(i4)) - lnk1)/2d0
               fco2 = dlog(1d12*p)
               y(i4) = 1d0 

            else 

               call rko2 (c1,1)

            end if 
         
         end if 

      else 

         rat = xc/(1d0-xc)
c                                 to make this compatible with the normal perplex
c                                 codes use lnk expressions, for the lv-version use
c                                 gzero.

c                                 evaluate K's and correct for pressure
c         go =  gzero (2)
c         go2 = gzero (ipoint+2)
c         gsio  = gzero (ipoint+1)
c         gsio2 = gzero (ipoint+3)
c         gsi = gzero (1)
c         lnk1 = (2d0*go - go2)/r/t      
c         lnk2 = (gsio2 - go - gsio)/r/t
c         lnk3 = (gsio - go - gsi)/r/t
c                                  k_1 2 O = O2, HSC K
         lnk1 = (-0.9214495D6/t + 0.6234471D5)/t - 0.1631235D2
c                                  k_2 from shornikov enthalpy, SiO2 = SiO + O
c         lnk2 = (-1.133204d+06/t - 5.491882d+04)/t + 1.710990d+01
c                                  HSC   
         lnk2 = (-0.1133204D7/t - 0.5491882D5)/t + 0.1710990D2 
c                                  SiO = Si + O, shornikov H_SiO 
c         lnk3 = (1.906315d6/t - 1.005993d5)/t + 1.664069d1
c                                  HSC
         lnk3 = (4.187931d5/t - 9.909023d4)/t + 1.627557d1
c                                 c1 = exp(lnK_1)*p => 2 O = O2
c                                 c3 = exp(lnK_3)/p => SiO = Si + O
c                                 c2 = exp(lnK_2)/p => SiO2 = SiO + O
         c1 = dexp(lnk1) * p
         c2 = dexp(lnk2) / p
         c3 = dexp(lnk3) / p
c                                 non-degenerate compositions:
         if (rat.ge.0.5d0) then 
c                                 use a liquid like starting guess
            call rksi5 (c1, c2, c3, bad)

            if (bad) then 

               ibad = ibad + 1 

               call setbad (fh2o)

            else 

               itic = itic + 1
               igood = igood + 1

               fh2o = dlog(p*g(i3)*y(i3)) 
               if (y(i5).ne.0d0) then 
                  fco2 = dlog(p*g(i5)*y(i5))
               else if (y(i2).ne.0d0) then 
                  fco2 = lnk3 + dlog(g(i2)*y(i2)/g(i3)/y(i3))
               else if (y(i1).ne.0d0) then 
                  fco2 = lnk2 + lnk3 + 
     *                   dlog(g(i1)*y(i1)/p/(g(i3)*y(i3))**2)
               else
                  write (*,*) 'wugga rksi5 ',t,p,xc,y
               end if 

            end if 

            g1 = fco2*xc + fh2o*(1d0-xc)

            do i = 1, isp
               oy(ins(i)) = y(ins(i))
               og(ins(i)) = g(ins(i))
            end do 

            xfs = fco2
            xfo = fh2o
            xv = vol

         else

            g1 = 1d99

         end if 

         call rksi5a (c1, c2, c3, bad)

         if (bad) then 

            ibad = ibad + 1 

            call setbad (fh2o)

         else 

            itic = itic + 1
            igood = igood + 1

            fh2o = dlog(p*g(i3)*y(i3)) 
            if (y(i5).ne.0d0) then 
               fco2 = dlog(p*g(i5)*y(i5))
            else if (y(i2).ne.0d0) then 
               fco2 = lnk3 + dlog(g(i2)*y(i2)/g(i3)/y(i3))
            else if (y(i1).ne.0d0) then 
               fco2 = lnk2 + lnk3 + 
     *                dlog(g(i1)*y(i1)/p/(g(i3)*y(i3))**2)
            else
               write (*,*) 'wugga rksi5 ',t,p,xc,y
            end if 

         end if 

         g2 = fco2*xc + fh2o*(1d0-xc)

         dg = g1 - g2

         if (dg.lt.1d-5) then 

            do i = 1, isp
               y(ins(i)) = oy(ins(i))
               g(ins(i)) = og(ins(i))
            end do 

            fco2 = xfs 
            fh2o = xfo
            vol = xv 

         end if 

      end if 
c                                 convert to g-atom vol
      n = 3d0*y(14) + 2d0*y(13) + 2d0*y(7) + y(12) + y(15)
c                                 convert to "molar amounts"         
      vol = vol/n/1d1

      end 

      subroutine rksi5a (c1,c2,c3,bad)
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 5 species silica vapor

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)

c derivation and data sources in maple work sheet Si-O_rk5_R=R_speciation.mws
c                                 JADC 6/12
c the Si-O2_rk4_v1_speciation.mws sheet for rksi4a may be more rational.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer iavg, ins(5), isp, nit, i1, i2, i3, i4, i5, i

      logical bad, both

      double precision c1,c2,c3,rat,rp1,rm1,nsi,no,oymin,nymin,oy(nsp),
     *                 r2p1,r2m1,dquart,oymax,nymax,by(nsp), bg(nsp), 
     *                 bvol, gold, gnew

      external dquart

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer iroots
      logical switch, rkmin, min
      common/ rkdivs /iroots,switch,rkmin,min

      double precision vol
      common/ cst26 /vol

      double precision vrt
      integer irt
      logical sroot
      common/ rkroot /vrt,irt,sroot

      save isp, ins, i1, i2, i3, i4, i5, both, iavg
      data isp, ins, i1, i2, i3, i4, i5, both, iavg
     *                                      /5, 14, 13, 12, 7, 15, 
     *                                          14, 13, 12, 7, 15, 
     *                                          .false., 1/
c----------------------------------------------------------------------
c                                 get pure species fugacities
      call mrkpur (ins, isp)
c                                 zero species in case of degenerate compositions
      do i = 1, isp
         y(ins(i)) = 0d0
      end do 
c                                 some inner loop constants
c                                 rat = nsi/no = xc/(1-xc) 
      rat = xc/(1d0-xc)

      rp1    = rat + 1d0
      rm1    = rat - 1d0
      r2p1    = 2d0*rat + 1d0
      r2m1    = 2d0*rat - 1d0

      switch = .false.
      both = .false.
c                                switch root entry point
10    nit = 0 
      oymin = 1d0
      oymax = 0d0 
      bad = .false. 
c                                 solve ideal gas speciation
      a0 = -c2 * c3 / c1  
      a1 = c2 * (rm1 + c3 * rp1) / c1
      a2 = (c2 * c3 * r2p1 + (r2m1 + c2) / c1) 
      a3 = (c2 * rp1 - rm1 / c1) 
c                                 monatomic O     
      call newton (dquart,1d0,0d0,1d-12,y(i3),bad)

      if (bad) then
         write (*,*) 'paused yo',p,t,xc
         call errpau
      end if
c                                 back calculate remaining fractions:
c                                 K1 => O2: 
      y(i4) = c1*y(i3)**2
c                                 mass balance => sio: this might be singular
c                                 at R = 1/2?
      y(i2)  = y(i3)*((2d0 - y(i3))*rat - 1d0 + y(i3) + y(i4)) / rat 
     *             / (y(i3) + 2d0*c3)

      if (y(i2).lt.0d0) y(i2) = nopt(50)
c                                  K3 => Si
      y(i5) = c3/y(i3)*y(i2)
c                                 closure => sio2: 
      y(i1) = 1d0 - y(i2) - y(i3) - y(i4) - y(i5)

      if (y(i1).lt.0d0) y(i1) = nopt(50)


      call mrkmix (ins, isp, iavg)    

      do 
c                                 solve mass balance: yo^4 + a3*yo^3 + a2*yo^2 + a1*yo + a0) 
c                                 for yo
         a0 = -c2 * c3 / c1 * g(i1) * g(i4) / g(i5) / g(i3) ** 4 
         a1 = c2 * g(i1) * g(i4) * (rm1 * g(i3) / g(i2) + c3 * rp1 
     *            / g(i5)) / g(i3) ** 4 / c1
         a2 = (c2 * g(i1) * g(i3) / g(i5) * c3 * r2p1 + g(i4) * 
     *           (r2m1 * g(i3) + c2 * g(i1) / g(i2)) / c1) / g(i3) ** 3

         a3 = (c2 * g(i1) * g(i3) / g(i2) * rp1 - rm1 * g(i4) / c1) 
     *           / g(i3) ** 2
c                                 monatomic O     
         call newton (dquart,1d0,0d0,1d-12,y(i3),bad)

         if (bad) then 
           
            exit 

         else if (y(i3).eq.0d0) then 
c                                 may not find the root if switch on 
c                                 first iteration       
            y(i3) = nopt(50)

         else if (isnan(y(i3)).or.y(i3).le.0d0.or.y(i3).eq.nopt(50)) 
     *                                                             then

            bad = .true.
            exit 

         end if 
c                                 back calculate remaining fractions:
c                                 K1 => O2: 
         y(i4) = c1/g(i4)*(y(i3)*g(i3))**2
c                                 mass balance => sio: this might be singular
c                                 at R = 1/2?
         y(i2)  = g(i5)*y(i3)*g(i3)*((2d0 - y(i3))*rat 
     *             - 1d0 + y(i3) + y(i4)) / rat 
     *             / (g(i5)*y(i3)*g(i3) + 2d0*c3*g(i2))
c
         if (y(i2).le.0d0) then 

            if (dabs(y(i2)).lt.nopt(50)) then 

               y(i2) = 0d0

            else 

               bad = .true.
               exit 

               do i = 1, isp
                  y(ins(i)) = 0d0
               end do 
c                                 two 4 species bailouts
c                                 4a has Si-O2-SiO-SiO2
c                                 4 has SiO2-SiO-O-O2,
c                                 neither is appropriate for reduced 
c                                 liquids.
               call rksi4a (c1,c2,iavg,bad)

               if (bad) then
 
                  do i = 1, isp
                     y(ins(i)) = 0d0
                  end do 

                  call rksi4 (bad,iavg)
     
               end if 

               return

            end if 

         end if 
c                                  K3 => Si
         y(i5) = c3/g(i5)/y(i3)/g(i3)*y(i2)*g(i2)
c                                 closure => sio2: 
         y(i1) = 1d0 - y(i2) - y(i3) - y(i4) - y(i5)

         if (y(i1).lt.0d0) then

            if (dabs(y(i1)).lt.nopt(50)) then 
               y(i1) = 0d0
            else 
               bad = .true.
               exit
            end if  

         end if 

         nymin = 1d0
         nymax = 0d0 
         no = 0d0      
c                                 renormalize, this helps! 
         do i = 1, isp
            no = no + y(ins(i))
         end do

         do i = 1, isp
            y(ins(i)) = y(ins(i))/no
         end do 

         do i = 1, isp

            if (y(ins(i)).gt.nymax) nymax = y(ins(i))
            if (y(ins(i)).lt.nymin.and.
     *          y(ins(i)).gt.0d0) nymin = y(ins(i))

         end do

         nsi = y(14) + y(13) + y(15)
         no  = 2d0*(y(7)+y(14)) + y(13) + y(12) 

         if (dabs(nymax-oymax).lt.nopt(50).and.
     *        dabs(nymin-oymin).lt.nopt(50).and.
     *        dabs(xc-nsi/(nsi+no)).lt.nopt(50) ) then 

            exit 

         else if (nit.gt.iopt(21)) then  
 
            bad = .true.
            exit
 
         end if 
c                                 make change less violent
c                                 this does help!
         if (nit.gt.1) then 

            do i = 1, isp

               y(ins(i)) = (y(ins(i))+oy(ins(i)))/2d0

            end do  

         end if 
c                                 get new gamma's
            call mrkmix (ins, isp, iavg)           
c                                could converge to the wrong speciation
            if (sroot.and.iroots.eq.3) then 
               both = .true.
            else if ((.not.both).and.iroots.eq.3) then
               both = .true.
               min = rkmin
            end if

         nit = nit + 1
         oymin = nymin
         oymax = nymax

         do i = 1, isp
            oy(ins(i)) = y(ins(i))
         end do 

      end do 

      if (both) then

         if (switch) then
c                                on the second solution
            gnew = 1d9
            gold = 1d9
            if (y(i3).gt.1d-20.and.y(i5).gt.1d-20)  
     *         gnew =  (1d0-xc) * dlog(g(i3)*y(i3)) 
     *                     + xc * dlog(g(i5)*y(i5))
            if (by(i3).gt.1d-20.and.by(i5).gt.1d-20) 
     *         gold =  (1d0-xc) * dlog(bg(i3)*by(i3)) 
     *                     + xc * dlog(bg(i5)*by(i5))

            if (gold.lt.gnew) then 
c                                swap to the old solution
               do i = 1, isp
                  y(ins(i)) = by(ins(i))
                  g(ins(i)) = bg(ins(i))
               end do 

               vol = bvol
c                                for sroot, 1st solution is good
               switch = .false.
   
            end if             

         else 
c                                save old solution
            do i = 1, isp
               by(ins(i)) = y(ins(i))
               oy(ins(i)) = 0d0 
               bg(ins(i)) = g(ins(i))
            end do 

            bvol = vol

            switch = .true.

            call mrkpur (ins, isp)

            goto 10

         end if  

      end if 

      end 

      subroutine idsi5 
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 5 species ideal silica vapor

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet Si-O_rk5_R=R_speciation.mws
c                                 JADC 6/12
c the Si-O2_rk4_v1_speciation.mws sheet for rksi4a may be more rational.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ins(5), isp, i1, i2, i3, i4, i5, i, itic, igood, ibad

      logical bad

      double precision c1,c2,c3,rat,rp1,rm1,
     *                 r2p1,r2m1,lnk1,lnk2,lnk3,dquart

      external dquart 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      integer iam
      common/ cst4 /iam

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3

      save isp, ins, i1, i2, i3, i4, i5, itic, igood, ibad 
      data isp, ins, i1, i2, i3, i4, i5, itic, igood, ibad 
     *                                      /5, 14, 13, 12, 7, 15, 
     *                                          14, 13, 12, 7, 15, 3*0/
c----------------------------------------------------------------------
c                                 zero species in case of degenerate compositions
      do i = 1, isp
         y(ins(i)) = 0d0
         g(ins(i)) = 1d0
      end do 
c                                 degenerate compositions:
      if (xc.eq.1d0) then 
c                                 pure Si
         fh2o = dlog(1d8*p)
         fco2 = dlog(p*g(i5))
         y(i5) = 1d0 
         
         return

      end if 
c                                 evaluate K's and correct for pressure
c                                 c1 = exp(lnK_1)*p => 2 O = O2, HSC K
      lnk1 = (-0.9214495D6/t + 0.6234471D5)/t - 0.1631235D2
      c1 = dexp(lnk1) * p

      if (xc.eq.0d0) xc = nopt(50) 
c                                 c2 = exp(lnK_2)/p => SiO2 = SiO + O, HSC K 
c     c2 = dexp((-0.1133204D7/t - 0.5491882D5)/t + 0.1710990D2) / p
c                                  k_2 from shornikov enthalpy
      lnk2 = (-1.133204d+06/t - 5.491882d+04)/t + 1.710990d+01
      c2 = dexp(lnk2)/p
c                                 c3 = exp(lnK_3)/p => SiO = Si + O, shornikov H_SiO 
      lnk3 = (1.906315d6/t - 1.005993d5)/t + 1.664069d1
      c3 = dexp(lnk3)/p
c                                 some inner loop constants
c                                 rat = nsi/no = xc/(1-xc) 
      rat = xc/(1d0-xc)
c                                 set singular compositions if within 
c                                 speciation tolerance.
      if (dabs(rat-0.5d0).lt.nopt(50)) then 
         rat = 0.5d0
      else if (dabs(rat-1d0).lt.nopt(50)) then 
         rat = 1d0
      end if 

      rp1    = rat + 1d0
      rm1    = rat - 1d0
      r2p1    = 2d0*rat + 1d0
      r2m1    = 2d0*rat - 1d0
c
      bad = .false. 
c                                 solve mass balance: yo^4 + a3*yo^3 + a2*yo^2 + a1*yo + a0) 
c                                 for yo
         a0 = -c2 * c3 / c1  
         a1 = c2 * (rm1 + c3 * rp1) / c1
         a2 = (c2 * c3 * r2p1 + (r2m1 + c2) / c1) 
         a3 = (c2 * rp1 - rm1  / c1) 
c                                 monatomic O     
         call newton (dquart,1d0,0d0,1d-12,y(i3),bad)

         if (isnan(y(i3)).or.y(i3).le.0d0.or.y(i3).eq.nopt(50)) then

            bad = .true.

         end if 
c                                 back calculate remaining fractions:
c                                 K1 => O2: 
         y(i4) = c1/g(i4)*(y(i3)*g(i3))**2
c                                 mass balance => sio: this might be singular
c                                 at R = 1/2?
         y(i2)  = g(i5)*y(i3)*g(i3)*((2d0 - y(i3))*rat 
     *             - 1d0 + y(i3) + y(i4)) / rat 
     *             / (g(i5)*y(i3)*g(i3) + 2d0*c3*g(i2))
c                                  K3 => Si
         y(i5) = c3/g(i5)/y(i3)/g(i3)*y(i2)*g(i2)
c                                 closure => sio2: 
         y(i1) = 1d0 - y(i2) - y(i3) - y(i4) - y(i5)

         if (y(i1).lt.0d0) then

            if (dabs(y(i1)).lt.nopt(50)) then 
               y(i1) = 0d0
            else 
               bad = .true.
            end if  

         end if 



      if (bad) then 

         ibad = ibad + 1 

         call setbad (fh2o)
         return 

      end if 

      fh2o = dlog(p*g(i3)*y(i3)) 
      if (y(i5).ne.0d0) then 
         fco2 = dlog(p*g(i5)*y(i5))
      else if (y(i2).ne.0d0) then 
         fco2 = lnk3 + dlog(g(i2)*y(i2)/g(i3)/y(i3))
      else if (y(i1).ne.0d0) then 
         fco2 = lnk2 + lnk3 + dlog(g(i1)*y(i1)/p/(g(i3)*y(i3))**2)
      else
         write (*,*) 'wugga rksi5 ',t,p,xc,y
      end if 

c     if (itic.gt.200000) then 
c        itic = 0 
c        write (*,*) 'good,bad:',igood,ibad
c     end if 

c      if (nit.gt.20) write (*,*) 'rk5 long nit',nit
      end

      subroutine rkparm (ins, isp)
c-----------------------------------------------------------------------
c subroutine to return standard rk a and b terms. 

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp     - the number of species to be calculated.
c        t       - input from cst5, p(bars), t(K)
c        cm3/mol volumetric units (R = 83.1441)

c species indices:

c         1 = H2O
c         2 = CO2
c         3 = CO
c         4 = CH4 
c         5 = H2
c         6 = H2S
c         7 = O2
c         8 = SO2
c         9 = COS
c        10 = N2
c        11 = NH3
c        12 = O
c        13 = SiO
c        14 = SiO2
c        15 = Si  
c        16 = C2H6
c        17 = HF
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision brk(nsp), ark(nsp)

      integer ins(*), isp, i, k

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision a, b
      common/ rkab /a(nsp),b(nsp)

      save ark, brk
                             
      data 
     *     ark /0d0, 0d0, 16.98d6, 32.154811d6, 2.821d5, 89d6,
     *          174026d2, 133.1d6, 130d6, 136d5 , 631d5,
c             O, SiO, SiO2, Si, unk:
c             asio = aco/aco2 * asio2 => 1/16.722 => 424441664.6d0
c bp values for si, sio2
c     *          174026d2, 424441664.6d0,  7097834092d0, 2348660042d0, 
c max a values for si, sio2
c    *          174026d2, 424441664.6d0,  7500726468d0, 3767833334d0, 
     *          174026d2, 424441664.6d0,  7373939618d0, 3767833334d0, 
     *          98774720.4d0, 7284049.7d0/,
     *     brk /14.6,  29.7,  27.38, 29.681, 15.7699, 29.94,
     *          22.07,  37.4,  43d0,  23.42,   18.84, 
c             O, SiO, SiO2, Si, unk:
c             bsio = bco/bco2*bsio2 => 1/4.83
c             bo  ~ bo2
c using actual bco/bco2:
     *          22.07, 23.81d0, 25.83798814d0, 10.35788774, 45.139d0,
     *          17.93096733d0/
c rhocrit o2 - 1449 o - 724. sio - 1851 sio2 - 2325 si - 2711
c using bco/bco2 = 1/4.83 (true critical b's):
c    *          22.07, 5.148732998d0, 25.83798814d0, 10.35788774, 0d0/
c who knows what this was
c     *          22.07, 5.148732998d0, 24.85507977, 10.35765058, 0d0/
c----------------------------------------------------------------------
c      if (first) then 
c         write (*,*) 'enter fac'
c         read (*,*) fac
c         if (fac.eq.0d0) fac=16.77
c         first = .false.
c      end if 

      do k = 1, isp

         i = ins(k)
         b(i) = brk(i)

         if (i.eq.1) then 

            if (t.gt.3d2) then 
c                                 MRK dispersion term for H2O
               a(1) = .1452535403d8 
     *                          + t*(306893.3587d0 + t*(-307.9995871d0
     *                          + t*(.09226256008d0 -.2930106337d-5*t)))
            else 
c                                 this is a gimmick to make the MRK predict
c                                 liquid water at the reference 298 K, these 
c                                 parameters give a molar voume of 18.0686 cm3/mol
c                                 at 298.15 and would give a boiling point ~370 K.
c                                 see mrk_water.mws. 
              b(1) = 16d0
              a(1) = 127354240d0

            end if 

         else if (i.eq.2) then 
c                                 MRK dispersion term for CO2
            a(2) =  92935540d0 + t*(-82130.73d0 + 21.29d0*t)

c         else if (i.eq.5) then 
c                                 MRK fit to NIS table for H2 at 10 kbar
c                                 400-1000 K.
c            b(5) = 12.81508162d0
c            a(5) = 0.391950132949994654D8 
c     *           + t * (-0.881231157499978144D5) 
c     *           + t**2 * 0.890185987380923081D2 
c     *           + t**3 * (-0.286881183333320412D-1)


         else if (i.eq.14) then 
c                                 MRK dispersion term for SiO2, from
c                                 sio2_mp_fit3.mws 
            a(14) = (-0.370631646315402D9 -88784.52
     *            + 0.710713269453173D8*dlog(t)
     *               -0.468778070702675D7/t + 
     *               (0.194790021605110D4*dsqrt(t) -0.110935131465938D6
     *               -0.120230245951606D2 * t) * t)*1d2
     *            +  32300.*(t-1999.) + 14.25*(t-1999.)**2

         else if (i.eq.15) then 
c                                 MRK dispersion term for Si 
            a(15) = (.131596431388077d7 - ((.380259023635694d-1*t 
     *      + .124090483523393d4)*t + .170392520137105d7)*dsqrt(t)
     *      + .151371320806448d6/dsqrt(t) + .427563259532326d7*dlog(t) 
     *      + (.108181901455347d2*t + 0.711400073165747d5)*t
     *      + 17737.22d0 
     *      -50.5d0*(t-1687.) - 2.04d-2*(t-1687.)**2
     *      )*1d2

         else 

            a(i) = ark(i)
 
         end if

         if (a(i).lt.0d0) a(i) = 1d0

      end do 

      a(13) = ark(14)/20d0
      b(13) = brk(13)

      end 

      double precision function dquart (x)
c----------------------------------------------------------------------
c function quart gives the value of a 4th order polynomial for x divided
c by it's derivative. coefficients if common block coeffs
c                                                  JADC, July 26, 2012. 
c----------------------------------------------------------------------
      implicit none

      double precision x,f,fx

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3 
c----------------------------------------------------------------------
      f  = a0 + x*(a1 + x*(a2 + x*(a3 + x)))

      fx = a1 + x*(2d0*a2 + x*(3d0*a3 + 4d0*x))

      if (fx.ne.0d0) then 
         dquart = -f/fx
      else
         dquart = 0d0
      end if 

      end 

      double precision function d32 (x)
c----------------------------------------------------------------------
c function d32 gives the value of (x + a2)*x + (a3*x  + a1)*x^(1/2) + a0
c by it's derivative. coefficients if common block coeffs
c                                                  JADC, Aug 31, 2012. 
c----------------------------------------------------------------------
      implicit none

      double precision x,f,fx,x12

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3 
c----------------------------------------------------------------------
      if (x.eq.0d0) then 
         d32 = 0d0
         return
      end if 

      x12 = dsqrt(x)       

      f  = (x + a2)*x + (a3*x  + a1)*x12 + a0

      fx = 2d0*x + a2 + (3d0*a3*x12 + a1/x12)/2d0

      d32 = -f/fx

      end 

      subroutine halver (func,max,min,tol,x)
c----------------------------------------------------------------------
c subroutine to halver locates a zero of func between min and max within 
c tolerance tol
c                                                  JADC, July 26, 2012. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision func

      external func

      double precision sign,max,min,dx,xdx,x,tol 
c----------------------------------------------------------------------
      x = min 
      dx = (max - min)/1d1
      xdx = x + dx

      do 

         sign = func(x)*func(xdx)

         if (sign.gt.0d0) then 

            if (xdx.lt.max) then 

               x = xdx
               xdx = x + dx

               if (xdx.gt.max) then 
c                                 hit upper limit
                  xdx = max 
                  dx = max - x
               end if 

            else if (xdx.eq.max) then
c                                 no root found
               x = -1d0
               exit 

            end if 

         else if (dx.gt.tol) then 
c                                 crossed the root between
c                                 x and dx
            dx = dx/2d0
            xdx = x + dx            
 
         else

            exit 

         end if 

      end do 

      end 


      subroutine newton (func,max,min,tol,x,bad)
c----------------------------------------------------------------------
c subroutine to locate a root of a function between min and max with 
c relative tolerance tol, func returns the value of -f(x)/diff(f,x)
c
c                                                  JADC, June 26, 2012. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad 

      double precision func

      external func

      integer nit 

      double precision max,min,dx,oldx,x,tol
c----------------------------------------------------------------------
      x = min 
      dx = func(x)

      if (isnan(dx).or.dx.le.0d0) then 
c                                 grad points to an out of range solution
         x = max
         dx = func(x)

         if (isnan(dx).or.dx.ge.0d0) then 

            bad = .true.
            return

         else if (x+dx.le.0d0) then 

            dx = -x/2d0

         end if 

      else if (x+dx.ge.1d0) then 

         dx = (1d0-x)/2d0

      end if 

      nit = 0 

      do 

         oldx = x
         x = x + dx
          
         if (dabs(x-oldx)/x.lt.tol) exit

         if (nit.gt.iopt(21)) then 

            bad = .true.
            exit

         end if 

         nit = nit + 1
         dx = func(x)

         if (dx.lt.0d0.and.x+dx.le.0d0) then 
            dx = -x/2d0
         else if (x+dx.ge.1d0) then 
            dx = (1d0-x)/2d0
         end if 

      end do 

      end 

      subroutine rksi30 
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 3 species silica liquid
c (a bail out routine for rksi5, assumes failure of rksi5 occurs for near
c stoichiometric liquid)

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet Si-O_rk3_R=R_speciation.mws
c                                 JADC 7/12
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ins(3), nit, i1, i2, i3, icon, ineg, ipos, iroots, i

      logical bad

      double precision c1,rat,rp1,rm1,vmin,vmax,x(3),oldy,lnk2,lnk3

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      integer iam
      common/ cst4 /iam

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3

      save ins 
      data ins, i1, i2, i3/14, 12, 15, 14, 12, 15/
c----------------------------------------------------------------------
c                                 rat = nsi/no = xc/(1-xc) 
      rat = xc/(1d0-xc)
c                                 evaluate K's and correct for pressure
c                                 SiO2 = SiO + O, HSC K 
c                                 k_2 from shornikov enthalpy
      lnk2 = (-1.133204d+06/t - 5.491882d+04)/t + 1.710990d+01
c                                 lnK_3 => SiO = Si + O, shornikov H_SiO 
      lnk3 = (1.906315d6/t - 1.005993d5)/t + 1.664069d1
      c1 =  dexp(lnk2+lnk3)/p**2
c                                 some inner loop constants
      rm1    = rat - 1d0
      rp1    = -c1 * (rat + 1d0) / rm1
      a2     = -(2d0*rat - 1d0) / rm1
      rm1    = c1/rm1
c                                 get pure species fugacities
      call mrkpur (ins, 5)

      nit = 0 
      oldy = 0d0
c                                 choose species for convergence test
c                                 SiO2, may blow at high T?
      icon = i2

      do 
c                                 solve (yo^3 + a2*yo^2 + a1*yo + a0) for yo
c                                 a2 is stoichiometric (defined above) and = 0 for rat = 1/2
         a0 = g(i1) * rm1 / g(i2)** 2 / g(i3) 
         a1 = g(i1) * rp1 / g(i2)** 2 / g(i3) 

         call roots3 (a2,a1,a0,x,vmin,vmax,iroots,ineg,ipos)

            bad = .true.

            do i = 1, iroots
c                                 find an acceptable root, if any... 
               if (isnan(x(i)).or.x(i).gt.1d0.or.x(i).lt.0d0) cycle  
c                                 monatomic O                 
               y(i2) = x(i)
c                                 back calculate remaining fractions:
c                                 Si: 
               y(i3) = (1d0 - y(i2)) / 
     *                 ((y(i2)*g(i2))** 2 * g(i3) / c1 / g(i1) + 1d0)

               if (isnan(y(i3)).or.y(i3).lt.0d0.or.y(i3).gt.1d0) cycle  
c                                 closure => sio2: 
               y(i1) = 1d0 - y(i2) - y(i3)

               if (y(i1).lt.0d0) cycle

               bad = .false.

               exit 

            end do 
        
            if (bad .or. dabs(oldy-y(icon)).lt.nopt(50)) exit
c                                 get new gamma's
            call mrkmix (ins, 5, 1)

            oldy = y(icon)

            nit = nit + 1

            if (nit.gt.iopt(21)) then 

               write (*,*) 'wug'

            end if 

            if (nit.lt.iopt(21)) cycle
 
            exit 

         end do

      if (bad) then 

          if (nit.gt.iopt(21)) then

             write (*,'(a,2(g12.6,1x))') 
     *            'ugga rksi30 not converging T,P:',t,p,xc

          else 

             write (*,'(a,5(g12.6,1x))') 
     *            'ugga rksi30 not valid solution T,P:',t,p,xc

          end if 

          call setbad (fh2o)

      else 

         fh2o = dlog(p*g(i2)*y(i2))
         fco2 = dlog(p*g(i3)*y(i3))

      end if 

      end

      subroutine mrkhen (ins,isp,ir,iavg)
c-----------------------------------------------------------------------
c subprogram to compute the henryian fugacity coefficient of mrk species 
c in solvent species ir. 

c assumes prior call to rkparm and mrkpur (for v(ir))

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp    -  the number of species to be calculated.
c        ir     -  solvent pointer
c        iavg   -  specifies the averaging method for the a cross term

c output (to common cstcoh):

c        g(i)    - fugacity coefficient of ith species

c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(*),i,k,isp,ir,iavg
 
      double precision rt,e0,e1,e2,r,ax
 
      double precision p,t,xco2,u1,u2,tr,pr,rbar,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,rbar,ps

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision a, b
      common/ rkab /a(nsp),b(nsp)

      double precision vol
      common/ cst26 /vol
 
      save r
                             
      data r /83.1441/
c---------------------------------------------------------------------- 
      rt = r*t
      e1 = b(ir)*rt*dsqrt(t)
      e0 = dlog( 1d0 +b(ir)/v(ir) )/e1
      e2 = 1d0/(v(ir)-b(ir))

      do k = 1, isp

         i = ins(k)
c                                 skip solvent
         if (i.eq.ir) cycle 


         if (i.eq.14.and.ir.eq.15.or.i.eq.15.and.ir.eq.14) then
c                                 arithmetic mean mixing rule
            ax = 2d0/(1d0/a(ir) + 1d0/a(i))
         else if (iavg.eq.1) then 
c                                 geometric mean mixing rule
            ax = dsqrt(a(ir)*a(i))
         else if (iavg.eq.2) then 
c                                 arithmetic mean mixing rule
            ax = (a(ir) + a(i))/2d0
         else 
c                                 harmonic mean mixing rule
            ax = 2d0/(1d0/a(ir) + 1d0/a(i))
         end if 

         g(i) = dexp (
     *          b(i) * (a(ir) * (e0/b(ir) - 1d0 / (v(ir) + b(ir)) / e1)
     *          + e2) - 2d0 * ax * e0 + dlog(rt*e2/p) )

      end do 
   
      end

      subroutine rksi3  
c----------------------------------------------------------------------
c subroutine to compute henryian fugacities of O and S above an SiO2
c fluid (bailout routine for rksi5)

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet Si-O_rk3_R=R_speciation.mws
c                                 JADC 7/12
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ins(3), isp, i1, i2, i3, ineg, ipos, iroots, i

      logical bad

      double precision c1,rat,rp1,rm1,vmin,vmax,x(3),lnk2,lnk3

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      integer iam
      common/ cst4 /iam

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3 

      save ins, isp 
      data isp, ins, i1, i2, i3/3, 14, 12, 15, 14, 12, 15/
c----------------------------------------------------------------------
c                                 rat = nsi/no = xc/(1-xc) 
      rat = xc/(1d0-xc)
c                                 evaluate K's and correct for pressure
c                                 SiO2 = SiO + O, HSC K 
c                                 k_2 from shornikov enthalpy
      lnk2 = (-1.133204d+06/t - 5.491882d+04)/t + 1.710990d+01
c                                 lnK_3 => SiO = Si + O, shornikov H_SiO 
      lnk3 = (1.906315d6/t - 1.005993d5)/t + 1.664069d1
      c1 =  dexp(lnk2+lnk3)/p**2

      rm1    = rat - 1d0
      rp1    = rat + 1d0
c                                 pure species sio2 fugacity
      call mrkpur (ins, 1)

      call mrkhen (ins, isp, i1, 2) 
c                                 solve (yo^3 + a2*yo^2 + a1*yo + a0) for yo
c                                 a2 is stoichiometric (defined above) and = 0 for rat = 1/2
      a0 = c1* g(i1) / g(i2)** 2 / g(i3) / rm1
      a1 = -a0 * rp1
      a2 = (1d0 - 2d0*rat) / rm1

      call roots3 (a2,a1,a0,x,vmin,vmax,iroots,ineg,ipos)

      bad = .true.

      do i = 1, iroots
c                                 find an acceptable root, if any... 
         if (isnan(x(i)).or.x(i).gt.1d0.or.x(i).le.0d0) cycle  
c                                 monatomic O                 
         y(i2) = x(i)
c                                 back calculate remaining fractions:
c                                 Si: 
         y(i3) = (1d0 - y(i2)) / 
     *           ((y(i2)*g(i2))** 2 * g(i3) / c1 / g(i1) + 1d0)

         if (isnan(y(i3)).or.y(i3).le.0d0.or.y(i3).gt.1d0) cycle  
c                                 closure => sio2: 
         y(i1) = 1d0 - y(i2) - y(i3)

         if (y(i1).le.0d0) cycle

         bad = .false.

         exit 

      end do 
   

      if (bad) then 


         write (*,'(a,5(g12.6,1x))') 
     *            'ugga wugga not valid solution T,P:',t,p,xc

         call setbad (fh2o)

      else 

         fh2o = dlog(p*g(i2)*y(i2))
         fco2 = dlog(p*g(i3)*y(i3))

      end if 

      end

      subroutine rko2 (c1,iavg) 
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 2 species oxygen
c call by rksi5

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet Si-O_rk5_R=R_speciation.mws
c                                 JADC 6/12
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ins(2), nit, i3, i4, iavg

      double precision c1,oldy,a0,a1

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      save ins
      data ins, i3, i4 /12, 7, 12, 7/
c----------------------------------------------------------------------
      nit = 0 
      oldy = 0d0

      call zeroys
c                                iterate for non-ideality
      do 

         a0 = 2d0*c1*g(i3)**2
         a1 = dsqrt(g(i4)*(g(i4)+2d0*a0))
         y(i3) = (a1-g(i4))/a0
         if (y(i3).gt.1d0.or.y(i3).lt.0d0) y(i3) = -(a1+g(i4))/a0
         y(i4) = 1d0 - y(i3)
         if (dabs(oldy-y(i3)).lt.nopt(50)) exit  
c                                 get new gamma's
         call mrkmix (ins, 2, iavg)

         oldy = y(i3)
         nit = nit + 1

         if (nit.lt.iopt(21)) cycle 
         write (*,*) 'ugga wugga not converging on pure O'
         exit 

      end do 

      fco2 = dlog(1d12*p)
      fh2o = dlog(p*g(i3)*y(i3))

      end

      subroutine rksi4a (c1,c2,iavg,bad)
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in SiO2-SiO-Si-O2 silica 
c vapor, bailout routine for rksi5.

c p    - pressure, bar
c t    - temperature, K
c xc   - bulk Si/(Si+O) (molar)
c fh2o - ln(fO)
c fco2 - ln(fSi)

c derivation and data sources in maple work sheet Si-O2_rk4_v1_speciation.mws 
c                                 JADC 9/12
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer iavg, ins(4), isp, nit, i1, i2, i3, i4, i, ir, 
     *        ibad, igood, imed, itic  

      logical bad, henry

      double precision c1,c2,c12,rat,rp1,rm1,nsi,no,oymin,nymin,oy(nsp),
     *                 r2m1,r2p1,r2,a4,lnk1,d32,oymax,nymax

      external d32 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      integer iam
      common/ cst4 /iam

      double precision a0,a1,a2,a3 
      common/ coeffs /a0,a1,a2,a3 

      save isp, ins, i1, i2, i3, i4, ibad, igood, imed, itic 
      data isp, ins, i1, i2, i3, i4, ibad, igood, imed, itic/
     *                                   4, 14, 13, 15, 7,
     *                                      14, 13, 15, 7, 4*0/
c----------------------------------------------------------------------
c                                 get pure species fugacities
      call mrkpur (ins, isp)
c                                 some inner loop constants
c                                 rat = nsi/no = xc/(1-xc) 
      rat   = xc/(1d0-xc)
      r2    = 2d0 * rat 
      rm1   = rat - 1d0
      rp1   = rat + 1d0 
      r2m1  = r2 - 1d0
      r2p1  = r2 + 1d0
      c12   = dsqrt(c1*c2)
c
      nit = 0  
      bad = .false. 
      oymin = 1d0
      oymax = 0d0 

      do 
c                                 solve (yo2^2 + a3*yo2^(3/2) + a2*yo2 + a1*yo2^(1/2) + a0) for yo2
         a4 = g(i1) / dsqrt(g(i4)/c1) / g(i2)
         a0 = -g(i1) / g(i4) * c12 / g(i3)
         a1 = rm1 * a4
         a2 = -r2p1 * a0  + r2m1
         a3 = rp1 * a4
  
         call newton (d32,1d0,0d0,1d-12,y(i4),bad)
c                                 may not find the root if switch on 
c                                 first iteration
         if (bad) then 

            exit 

         else if (y(i4).eq.0d0) then 
       
            y(i4) = nopt(50)

         else if (isnan(y(i4)).or.y(i4).le.0d0.or.y(i4).eq.nopt(50)) 
     *                                                             then

            bad = .true.
            exit 

         end if 
c                                 back calculate remaining fractions, use
c                                 mass balance and equilibrium constants to assure > 0.
c                                 sio: mass balance is singular for y(sio) at R = 1
c                                 and singular for y(sio2) at R = 1/2. 
c                                 use closure instead:
         y(i2) = dsqrt(y(i4))*(1d0-y(i4))/(y(i4)*dsqrt(g(i4)/c1) 
     *   * g(i2)/g(i1) + dsqrt(y(i4)) + dsqrt(c2/g(i4))*g(i2)/g(i3))
c                                 sio2: equilibrium
         y(i1) = dsqrt(y(i4) * g(i4) / c1) * g(i2) * y(i2) / g(i1)
c                                 si: equilibrium
         y(i3) = dsqrt(c2/y(i4)/g(i4)) * g(i2)*y(i2)/g(i3)
c                                 calculate molar totals to test for
c                                 convergence
         henry = .false.
         nymin = 1d0
         nymax = 0d0 
         no = 0d0      

         do i = 1, isp
            no = no + y(ins(i))
            if (y(ins(i)).lt.0d0.or.no.gt.2d0) then 
               write (*,*) 'wock'
            end if 
         end do

         do i = 1, isp
            y(ins(i)) = y(ins(i))/no
         end do 

         do i = 1, isp

            if (y(ins(i)).gt.nymax) nymax = y(ins(i))
            if (y(ins(i)).lt.nymin.and.
     *          y(ins(i)).gt.0d0) nymin = y(ins(i))

            if (y(ins(i)).gt.1d0-nopt(50).and.y(ins(i)).le.1d0) then 
               ir = i
               henry = .true.
               exit 
            end if 

         end do

         nsi = y(14) + y(13) + y(15)
         no  = 2d0*(y(7)+y(14)) + y(13)  

         do i = 1, isp

            if (y(ins(i)).gt.nymax) nymax = y(ins(i))
            if (y(ins(i)).lt.nymin.and.
     *          y(ins(i)).gt.0d0) nymin = y(ins(i))

         end do

         if (dabs(nymax-oymax).lt.nopt(50).and.
     *        dabs(nymin-oymin).lt.nopt(50).and.
     *        dabs(xc-nsi/(nsi+no)).lt.nopt(50).and.
     *        dabs(nsi+y(7)-1d0).lt.nopt(50) ) then

            igood = igood + 1

            exit 

         else if (nit.gt.400.and.
     *        dabs(nymax-oymax)/nymax.lt.1d-3.and.
     *        dabs(nymin-oymin)/nymin.lt.1d0.and.
     *        dabs(xc-nsi/(nsi+no)).lt.nopt(50).and.
     *        dabs(nsi+y(7)-1d0).lt.nopt(50) ) then
          
            imed = imed + 1

            exit 

         else if (nit.gt.iopt(21)) then 

            bad = .true.
            exit           

         end if 
c                                 make change less violent
c                                 this does help!
         if (nit.gt.1.and.
     *       dabs(nymax-oymax)/nymax.gt.1d-3.or.
     *       dabs(nymin-oymin)/nymin.gt.1d0) then 

            do i = 1, isp

               y(ins(i)) = (0.5*y(ins(i))+0.5*oy(ins(i)))

            end do  

         end if 
c                                 get new gamma's
         if (v(14).lt.0*1d2.and.xc.gt.0.326.and.xc.lt.0.340) then

            fh2o = dlog(1d12*p)
            fco2 = dlog(1d12*p)
            return

         else 

            call mrkmix (ins, isp, iavg)           

         end if 

         nit = nit + 1
         oymin = nymin
         oymax = nymax
         do i = 1, isp
            oy(ins(i)) = y(ins(i))
         end do 

      end do 

      itic = itic + 1

      if (bad) then 

         ibad = ibad + 1

c         if (nit.ge.60.and.iwarn.lt.100) then

c            write (*,'(a,2(g12.6,1x))') 
c     *            'ugga wugga rk4a not converging T,P:',t,p,xc

c         else if (iwarn.lt.1e6) then 

c            write (*,'(a,5(g12.6,1x))') 
c     *            'ugga wugga rk4a not valid solution T,P:',t,p,xc

c         end if 

c         iwarn = iwarn + 1 

c         if (iwarn.eq.100) call warn (49,t,0,'RKSI4a')

         call setbad (fh2o)

         return 

      end if 

      if (itic.gt.iopt(21)) then
         write (*,*) 'rk4a: igood,imed,ibad: ',igood,imed,ibad
         itic = 0
      end if 

c      if (nit.gt.20) write (*,*) 'rk4a long it:',nit

      fco2 = dlog(p*g(i3)*y(i3))
      fh2o = (dlog(p*g(i4)*y(i4)) - lnk1) / 2d0

      end

      double precision function gerk (x)
c----------------------------------------------------------------------
c subroutine to compute free energy of solution of an Si-O MRK fluid
c assumes
c x(1) = ysio2
c x(2) = ysio
c x(3) = yo
c x(4) = yo2
c x(5) = ysi
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer iavg, ins(5), isp, i

      double precision x(*)

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol

      save isp, ins, iavg
      data isp, ins, iavg /5, 14, 13, 12, 7, 15,  1/
c----------------------------------------------------------------------
c                                 load composition
      do i = 1, isp
         y(ins(i)) = x(i)
      end do 

      call mrkmix (ins, isp, iavg)    

      gerk = 0d0

      do i = 1, isp
         if (x(i).eq.0d0) cycle
         gerk = gerk + x(i)*dlog(g(ins(i))*p*x(i))
      end do 

      gerk = r*t*gerk
c                                 convert to j/bar from cm3, only for lv version
      vol = vol/1d1

      end

      double precision function ghybrid (x)
c----------------------------------------------------------------------
c  generic hybrid EoS called for solution model type ksmod = 27.
c  reads composition from solution model 

c  10/28/16, JADC.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k

      double precision gtemp, x(*)

      character specie*4
      integer ins, isp
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)
c----------------------------------------------------------------------
      call zeroys
c                                 load the solution model y's into the
c                                 fluid species y array
      do k = 1, isp
         y(ins(k)) = x(k)
      end do 
c                                 calculate hybrid fugacity coefficients
      call mrkmix (ins, isp, 1)
c                                 compute hybrid activity effect
      gtemp = 0d0

      do k = 1, isp
         if (x(k).le.0d0) cycle
         gtemp = gtemp + x(k) * dlog(x(k)*g(ins(k))/g0(ins(k)))
      end do
  
      ghybrid = r*t*gtemp

      end

      subroutine setbad (fo2)
c----------------------------------------------------------------------
c set standard fugacities to "bad" values if fluid routine speciation
c routine fails. 5/17, JADC
c----------------------------------------------------------------------
      implicit none

      double precision fo2

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol

      double precision f
      common/ cst11 /f(3)
c----------------------------------------------------------------------
      vol = 0d0
      f(1) = dlog(1e4*p)
      f(2) = f(1)
c                                 it's not at all clear to me that f(3)
c                                 is ever used, it ought to be deleted
c                                 to reduce confusion. 
      f(3) = f(1)
      fo2 = f(1)
      end 

      subroutine gcohx6 (fo2)
c----------------------------------------------------------------------
c  program to calculate GCOH fluid properties as a function of XO 
c  see Connolly (1995) and/or coh_speciation_with_ethane.mws for 
c  details. Uses MRK activities with pure fluid fugacities as specified 
c  by the hybrid_EoS option.

c  replaces hocgra and hocmrk.

c  10/28/16, JADC.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(6),jns(3),nit,i

      logical bad, lmt

      double precision oy5,fo2,ytot,t4y3,t3y3,t2y5,t4y5,det,x,dy5,dy3,
     *       c1,c2,c3,c4,t1,t2,t3,t4,m,dm3,dm5,c,dc3,dc5,nh,rat,y5,y3

      double precision dinc
      external dinc

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
      save / cst59 /

      double precision eqk
      common / csteqk /eqk(nsp)

      save ins, jns
      data ins, jns/1,2,3,4,5,16,1,2,4/
c----------------------------------------------------------------------
      nit = 0
      oy5 = 0d0
      bad = .false.
c                                 check for in bounds composition
      call xcheck (xo,lmt)
c                                 compute equilibrium constants, returned
c                                 in csteqk
      call seteqk (ins,6,elag)
c                                 compute pure mrk fluid properties
      call mrkpur (ins,6)
c                                 compute hybrid pure fluid props
      call hybeos (jns,3)

      call zeroys

      if (lmt) then 
         call setbad (fo2)
         return
      end if

      c1 = dexp (eqk(4)) * p
      c2 = dexp (2d0*eqk(16) - 3d0*eqk(4)) * p
      c3 = dexp (eqk(2) - 2d0*eqk(3)) * p 
      c4 = dexp (eqk(1) - eqk(3)) * p
c                                 initial guess, assume near binary
      x = 1d0 + xo 
      rat = xo/(1d0-xo) 

      if (dabs(xo-r13).lt.nopt(5)) then 
         y3 = 1d0/dsqrt(c4)
         y5 = y3 
      else if (xo.gt.r13) then
         y3 = dsqrt((c3*x*(3d0*xo - 1d0)))/c3/x
         y5 = 2d0 * (1d0 - xo)/c4/y3/x
      else 
         y5 = dsqrt((c1*x*(1d0 - 3d0*xo)))/c1/x
         y3 = 4d0*xo/c4/y5/x
      end if

      do
c                                 iteration loop
         t1 = c1 / g(ins(4)) * g(ins(5)) ** 2
         t2 = dsqrt(c2 * (t1 * g(ins(4))) ** 3) / g(ins(6))
         t3 = c3 / g(ins(2)) * g(ins(3)) ** 2
         t4 = c4 / g(ins(1)) * g(ins(5)) * g(ins(3))

         t4y3 = t4*y3
         t3y3 = t3*y3
         t2y5 = t2*y5
         t4y5 = t4*y5

         nh = ((3d0*t2y5 + 2d0*t1)*y5 + t4y3 + 1d0)*y5
c                                 no/nh
         x = (t3*y3 + (1d0 + t4y5)/2d0)*y3/nh
c                                 mass balance eq
         m = rat - x
c                                 diff(m,y5)
         dm5 = (x*((9d0*t2y5 + 4d0*t1)*y5 + t4y3 + 1d0) - t4y3/2d0)/nh
c                                 diff(m,y3)
         dm3 = (x*t4y5 - 2d0*t3y3 - 0.5d0 - t4y5/2d0)/nh
c                                 closure eq
         c = ((-t2y5 - t1)*y5 - t4y3 - 1d0)*y5 + 1d0 - (t3y3 + 1d0)*y3 
c                                 diff(c,y5)
         dc5 = (-3d0*t2y5 - 2d0*t1)*y5 - t4y3 - 1d0
c                                 diff(c,y3)
         dc3 = -2d0*t3y3 - t4y5 - 1d0
c                                 newton-raphson increments
         det = dc5 * dm3 - dc3 * dm5
         dy5 = -(c * dm3 - dc3 * m) / det
         dy3 =  (c * dm5 - dc5 * m) / det
c                                 add the increment
         y5 = dinc(y5,dy5)
         y3 = dinc(y3,dy3)
c                                 back calculate y's
         y(ins(5)) = y5
         y(ins(3)) = y3
         y(ins(4)) = t1 * y5**2 
         y(ins(6)) = t2 * y5**3
         y(ins(2)) = t3 * y3**2
         y(ins(1)) = t4 * y5 * y3
c                                 compute ytot
         ytot = 0d0

         do i = 1, 6 
            ytot = ytot + y(ins(i))
         end do 
c                                 renormalize
         do i = 1, 6
            y(ins(i)) = y(ins(i)) / ytot
         end do 
c                                 check for convergence, could 
c                                 do better than this
         if (dabs(y5-oy5).lt.nopt(50)) exit
c                                 check if iteration count exceeded
         if (nit.gt.iopt(21)) then 
            call warn (176,y5,nit,'GCOHX6')
            bad = .true.
            exit
         end if
c                                 calculate new fugacity coefficients
         call mrkhyb (ins, jns, 6, 3, 1)

         oy5 = y5

         nit = nit + 1

         y5 = y(ins(5))
         y3 = y(ins(3))

      end do 

      if (bad) then

         call setbad (fo2)

      else 

         if (hu.eq.0) then 
c                                 normal fugacities
            fh2o = dlog(g(ins(1))*p*y(ins(1)))
            fco2 = dlog(g(ins(2))*p*y(ins(2)))
            fo2 = 2d0*(dlog(g(ins(3))*p*y3) - eqk(3))

         else 
c                                 projecting through graphite
            fh2o = dlog(g(ins(5))*p*y5)
            fco2 = 2d0*(dlog(g(ins(3))*p*y3) - eqk(3))

         end if

      end if 

      do i = 1, 3
         vol = vol + y(jns(i))*vh(jns(i))
      end do 

      end

      double precision function dinc (y,dy)
c----------------------------------------------------------------------
c function to increment 0 < y < 1 by dy, if the increment would cause
c y to exceed its bounds, the increment is taken to be half the 
c interval to the bound. 10/2016, JADC.
c----------------------------------------------------------------------
      implicit none

      double precision y, dy
c----------------------------------------------------------------------
      if (y+dy.ge.1d0) then
         dinc = 0.5d0 + y / 2d0
      else if (y+dy.le.0d0) then 
         dinc = y / 2d0
      else 
         dinc = y + dy
      end if 

      end 

      subroutine seteqk (ins,isp,ac)
c----------------------------------------------------------------------
c program to compute standard (from the elements convention) ln equilibrium 
c constants for 1 mole of the nsp fluid species. 

c ac - ln(a[graphite])
c 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision dg,t2,t3,ac

      integer i,j,isp,ins(*)

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

c----------------------------------------------------------------------
      t2 = t * t
      t3 = t2 * t

      if (ac.ge.0d0) then 
c                                correct activity of graphite
c                                for diamond stability if necessary:
         call dimond (dg)
c                                graphite pressure effect and activity
c                                corrections:
         dg = dg + ac + p*( 1.8042d-06 + (0.058345d0 - 8.42d-08*p)/t ) 

      else 

         dg = 0d0

      end if 
 
      do i = 1, isp 

         j = ins(i)

         if (j.eq.1) then 
c                                h2o/robie
            eqk(j) = -7.028214449d0 + 30607.34044d0/t  
     *               - 475034.4632d0/t2 + 50879842.55d0/t3
         else if (j.eq.2) then 
c                                co2/robie
            eqk(j) = .04078341613d0 + 47681.676177d0/t 
     *              - 134662.1904d0/t2 + 17015794.31d0/t3 + dg
         else if (j.eq.3) then 
c                                co/robie
            eqk(j) = 10.32730663d0  + 14062.7396777d0/t
     *              - 371237.1571d0/t2 + 53515365.95d0/t3 + dg
         else if (j.eq.4) then 
c                                ch4/robie
            eqk(j) = -13.86241656d0 + 12309.03706d0/t
     *             - 879314.7005d0/t2 + .7754138439d8/t3 + dg
         else if (j.eq.6) then 
c                                h2s/ohmoto and kerrick
            eqk(j) = 10115.3d0/t - 0.791d0 * dlog (t) + 0.30164d0
         else if (j.eq.8) then 
c                                so2/ohmoto and kerrick
            eqk(j) = 43585.63147d0/t - 8.710679055d0
         else if (j.eq.9) then 
c                                cos/ohmoto and kerrick
            eqk(j) = 10893.52964d0/t - 9.988613730d0
         else if (j.eq.16) then 
c                                c2h6/HSC, 10/2016
            eqk(j) =  4.09702552d7/t3 - 8.01186095d5/t2 + 1.39350247d4/t
     *                      - 2.64306669d1 + 2d0 * dg
         end if 

      end do 

      end 

      subroutine dimon1 (agph)
c-----------------------------------------------------------------------
c dimon1 tests if p-t conditions are in the diamond stability field
c if they are it computes the activity of graphite needed to represent
c diamond stability for C-O-H fluid speciation routines. 

c dimon1 should replace dimond which automatically adds a graphite
c activity correction

c  agph - natural log of graphite activity to be added to graphite
c        free energy to equate g(diamond) = g(graphite). 

c polynomial functions fit by S. Poli, 4/27/04
c                                 JADC, 4/27/04
c-----------------------------------------------------------------------
      implicit none 

      double precision ptrans,agph

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

c                                 get transition pressure
c                                 (should really use the G
c                                 function for consistency):

      ptrans = 5284.165053d0 + (33.21515773d0 - .002106330992d0*t)*t

      if (p.gt.ptrans) then 
c                                 compute corrected graphite activity
         agph = 0.008423508384179629d0 
     *        + (4.693008650307614d-11 * p - 3.850380793502567d-5)*p
     *        + (1.4126916053951515d-3 + 1.393226795939807d-8*p 
     *        - 5.887505938975768d-7 * t)*t

      end if 

      end 

      subroutine zeroys
c---------------------------------------------------------------------
c zero complete fluid speciation array to avoid potential conflicts.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)
c---------------------------------------------------------------------
      do i = 1, nsp
         y(i) = 0d0
      end do

      end 

      subroutine xcheck (x,lmt)
c----------------------------------------------------------------------
c x check for speciation routines that can't handle xo = 0, 1. this
c causes too many possibilities for error in vertex/meemum so modified
c merely to avoid x out of range. feb 7, 2017. jadc,
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x

      logical lmt
c----------------------------------------------------------------------
c                                check if xo is <1, >0,
c                                reset if necessary
      if (x.lt.nopt(50)) then
         x = nopt(50)
c        lmt = .true.
      else if (x.gt.nopt(56)) then
         x = nopt(56)
c        lmt = .true.
c     else
c        lmt = .false.
      end if 

      lmt = .false.

      end 

      subroutine rkcoh6 (yo2,yc,deltag)
c----------------------------------------------------------------------
c subroutine to compute speciation and fugacites in 6 species COH fluid

c p    - pressure, bar
c t    - temperature, K
c yo2   - bulk fraction O2
c yc    - bulk fraction C

c derivation and data sources in maple work sheet rkcoh6.mws
c                                 JADC 1/15
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical bad, quit, limito, limitc

      integer ins(6), jns(3), i, j, k, itct, icase

      double precision a(3,4),yc,yo2,x0,xo3,deltag,yh2,xc,xo,
     *                 dg(nsp),eg(nsp),oy(7),s4,s5,gtot,ogtot,
     *                 xc5,xc3,xo1, x1,x2,x3,x4,x5,x6,x7,x8,x9,t4,t5,
     *                 d1,d2,d11,d22,d12,det,dp(2),b(3),dgs1,dgs2,s

      double precision vol
      common/ cst26 /vol

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision muc,muo2,muh2
      common / cst11 /muc,muo2,muh2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save t4,t5,ins,jns
      data t4,t5,ins,jns/1d-3,1d-3,1,2,4,3,5,7,1,2,4/
c----------------------------------------------------------------------
      bad = .false.

      yh2 = 1d0 - yo2 - yc

20    if (yc.ge.1d0/3d0+yo2.or.bad
     *                     .or.yh2.lt.nopt(50).or.yh2.ge.nopt(56)
     *                     .or.yo2.lt.nopt(50).or.yo2.ge.nopt(56)
     *                     .or.yc .lt.nopt(50).or.yc .ge.nopt(56)
     *                            ) then
c                                 above the ch4-co join
         deltag = 1d5

         return

      end if 
c                                 for reasons of stupidity convert 
c                                 c-o2-h2 coordinates to c-o-h
      xo = 2d0*yo2/(2d0-yc)
c                                 do not allow degenerate compositions:
      call xcheck (xo,limito)

      xc = yc/(2d0-yc)

      call xcheck (xc,limitc)
c                                 compute equilibrium constants in csteqk
      call seteqk (ins,6,0d0)
c                                 compute pure mrk fluid properties
      call mrkpur (ins,6)
c                                 compute hybrid pure fluid props
      call hybeos (jns,3)
c                                 for safety, zero the species array
      call zeroys

      if (limito.or.limitc) return
c                                 gh => gcork/gmrk, the g(i)'s are the fugacity coefficients
c                                 of the pure gas according to the various eos's (CORK, HSMRK, MRK)
c                                 could replace CORK with pitzer and sterner.
c                                 ---------------------------------------------
c                                 free energies for go2=gh2=gc=0
c                                 g0h2o/R/T - ln phi0
      dg(1) = (-eqk(1) - dlog(g(5)*dsqrt(g(7)*p)))
c                                 gco2/R/T 
      dg(2) = (-eqk(2) - dlog(g(7)))
c                                 gco/R/T 
      dg(3) = (-eqk(3)  + dlog(p/g(7))/2d0)
c                                 gch4/R/T 
      dg(4) = (-eqk(4) - dlog(g(5)*g(5)*p))
c                                 
      dg(5)  = - dlog(g(5))
      dg(7)  = - dlog(g(7))
c                                 ---------------------------------------------
      xc5 = 5d0*xc
      xc3 = 3d0*xc
      xo1 = xo + 1d0
      xo3 = 3d0*xo
      x1 = -xo + 1d0 + xc
      x2 = xc - xo1
      x4 = xc5 - xo3 - 1d0
      x5 = xo - 1d0 + xc
      x7 = xc + xo3 + 1d0
      x8 = xc3 + xo - 1d0
      x9 = xc - xo
c                                 choose independent species
c                                 first assume co2-h2o-ch4 ternary
      y(1) = (xc5 - xo1)/x2
      y(2) = (1d0 - xo3 - xc5)/x2/2d0
      y(4) = 1d0 - y(1) - y(2)
      y(7) = 0d0 

      do i = 1, 7
         oy(i) = y(i)
      end do 

      s4 = t4
      s5 = t5 

30    if (y(1).ge.0d0.and.y(2).ge.0d0.and.y(4).ge.0d0) then
c                                 inside co2-h2o-ch4 ternary
         ins(1) = 1
         ins(2) = 2
         ins(3) = 4
         ins(4) = 3
         ins(5) = 5

         icase = 1

         y(ins(1)) = -( (2d0*s4 + 3d0*s5 - 5d0)*xc 
     *                + (1d0 - s5)*xo + 1d0 - s5) / x2
         y(ins(2)) = (( (s4 + 3d0*s5 - 5d0)*xc 
     *                + (s5 + 3d0*s4 - 3d0)*xo 
     *                - s5 + s4 + 1d0) / x2) / 2d0
c                                 stoichiometric coefficients:
         x3 = -1d0 + 3d0*x9
         x6 = xo + xc3 - 1d0

         a(1,1) = -x1/x3
         a(1,2) = -x1/x3
         a(1,3) = -2d0*x9/x3
         a(2,1) = -x6/x3
         a(2,2) = 2d0*(xo1 - xc3)/x3
         a(2,3) = -2d0*xo/x3
         a(3,1) = x7/x3
         a(3,2) = 4d0*xc/x3
         a(3,3) = x4/x3

         b(1) = 1.5d0*a(2, 1) + 3d0*a(1,1) + 1.5d0 + a(3,1)
         b(2) = 2d0 + 1.5d0*a(2,2) + 3d0*a(1,2) + a(3,2)
         b(3) = 1.5d0*a(2,3) + 3d0*a(1,3) + a(3,3)

      else if (y(1).ge.0d0.and.y(2).lt.0d0.and.y(4).ge.0d0) then
c                                 inside h2-h2o-ch4 ternary
         ins(1) = 1
         ins(2) = 5
         ins(3) = 4
         ins(4) = 3
         ins(5) = 2

         icase = 2

         y(ins(1)) = -(  (6d0*s5 + 3d0*s4)*xc 
     *                 + (2d0 - 3d0*s4 - 2d0*s5)*xo - 2d0*s5 - s4) / x8
         y(ins(2)) = -(  (s4 - 2d0*s5 - 5d0)*xc 
     *                 + (3d0*s4 + 2d0*s5 - 3d0)*xo 
     *                 + s4 + 2d0*s5 + 1d0) / x8

         x0 = 4d0*xc

         a(1,1) = -x1/x0
         a(1,2) = -x1/x0
         a(1,3) = -x5/x0
         a(2,1) = -x2/xc/2d0
         a(2,2) = (xo1 - xc3)/xc/2d0
         a(2,3) = (xc5 -xo1)/xc/2d0
         a(3,1) = -x7/x0
         a(3,2) = (3d0*x9 - 1d0)/x0
         a(3,3) = -x4/x0

         b(1) = 3d0*a(1,1) + 1.5d0*a(2,1) + 2d0*a(3,1) + 1.5d0
         b(2) = 3d0*a(1,2) + 1.5d0*a(2,2) + 2d0*a(3,2) + 1d0
         b(3) = 3d0*a(1,3) + 1.5d0*a(2,3) + 2d0*a(3,3)

      else if (y(1).lt.0d0.and.y(2).ge.0d0.and.y(4).ge.0d0) then
c                                 inside co-co2-ch4 ternary
         ins(1) = 3
         ins(2) = 2
         ins(3) = 4
         ins(4) = 1
         ins(5) = 5

         icase = 3

         y(ins(1)) = (-s4 - 3d0*s5 + 5d0 + xo1/xc*(s4 - 1d0 + s5))/2d0
         y(ins(2)) = (-s4 + 3d0*s5 - 5d0 
     *                    + (xo3 + 1d0)/xc*(1d0 - s4 - s5))/4d0
c                                 stoichiometric coefficients:
         a(1,1) = x1/x8
         a(1,2) = x1/x8
         a(1,3) = -2d0*xc/x8
         a(2,1) = (1d0 - 3d0*x9)/x8
         a(2,2) = 2d0*(xo1 - xc3)/x8
         a(2,3) = -2d0*xo/x8
         a(3,1) = -x7/x8
         a(3,2) = 2d0*x2/x8
         a(3,3) = (xc5 + xo3 - 1d0)/x8

         b(1) = 1.5d0 + 3d0*a(1,1) + 1.5d0*a(2,1) + a(3,1)
         b(2) = 2d0   + 3d0*a(1,2) + 1.5d0*a(2,2) + a(3,2)
         b(3) = 3d0*a(1,3) + 1.5d0*a(2,3) + a(3,3)

      else
c                                 inside h2o-co2-o2 simplex
         ins(1) = 1
         ins(2) = 2
         ins(3) = 7
         ins(4) = 3
         ins(5) = 5

         icase = 4

         y(4) = 0d0 

         y(ins(1)) = -((s4 + 2d0*s5 - 2d0)*xc + (s4 - 2d0)*xo 
     *                 + 2d0 - s4 - 2d0*s5) / x2
         y(ins(2)) = ((s4 + 2d0*s5 -4d0) * xc + s4 * xo1) / x2

         a(1,1) = -x1/2d0
         a(1,2) = -x1/2d0
         a(1,3) = -x9
         a(2,1) = xc
         a(2,2) = xc - 1d0
         a(2,3) = 2d0*xc
         a(3,1) = (-xo1 - xc)/2d0
         a(3,2) = -x5/2d0
         a(3,3) = -x5

         b(1) = 1.5d0*a(2,1) + a(1,1) + 1.5d0 + a(3,1)
         b(2) = 2d0 + 1.5d0*a(2,2) + a(1,2) + a(3,2)
         b(3) = 1.5d0*a(2,3) + a(1,3) + a(3,3)

      end if
c                                 check if ok
      do i = 1, 2

         if (y(ins(i)).lt.0d0) then 
            s4 = s4/2d0
            s5 = s4

            do j = 1, 7
               y(j) = oy(j)
            end do 

            if (s5.gt.1d-10) goto 30
            bad = .true.
            goto 20 

         end if 

      end do 
c                                 pre-loop initialization
      do i = 1, 3

         j = i + 2

         y(ins(j)) = a(i,1)*y(ins(1)) + a(i,2)*y(ins(2)) + a(i,3)

      end do 

      gtot = 0d0
      ogtot = -99d0
      dp(1) = 0d0
      dp(2) = 0d0 

      call mrkhyb (ins, jns, 5, 3, 1)

      itct = 0
c                                 iteration loop
      do  

         do i = 1, 5
c                                 save old values for convergence test
            oy(i) = y(ins(i))

         end do 
c                                 update independent fractions
10       do i = 1, 2

            if (y(ins(i)) + dp(i) .le. 0d0) then 
               y(ins(i)) = y(ins(i))/2d0
            else if (y(ins(i)) + dp(i) .ge. 1d0) then 
               y(ins(i)) = y(ins(i)) + (1d0 - y(ins(i)))/2d0
            else 
               y(ins(i)) = y(ins(i)) + dp(i)
            end if 
c                                 check convergence
            if (dabs(y(ins(i))-oy(i)).gt.nopt(50)) then
               quit = .false.
            end if 

         end do 
c                                 update dependent fractions
         do i = 1, 3

            j = i + 2

            y(ins(j)) = a(i,1)*y(ins(1)) + a(i,2)*y(ins(2)) + a(i,3)

            if (y(ins(j)).le.0d0.or.y(ins(j)).ge.1d0) then 

               do k = 1, 2
                  dp(k) = dp(k)/2d0
                  if (dabs(dp(k)).lt.1d-30) then 
                     bad = .true.
                     goto 20
                  end if 
                  y(ins(k)) = oy(k)
               end do 

               quit = .true.

               goto 10 

            end if
c                                 check convergence
            if (dabs(y(ins(j))-oy(j)).gt.nopt(50)) then
              quit = .false.
            end if 

         end do 

         if (itct.gt.1.and.dabs((gtot-ogtot)/gtot).lt.nopt(50)) then
            quit = .true.
            do j = 1, 5

               if (dabs(y(ins(j))-oy(j)).gt.nopt(50)) then
               quit = .false.
               end if 

             end do      

            if (quit) exit

            if (itct.gt.iopt(21)) then 
               bad = .true.
               goto 20 
            end if 

         else if (itct.gt.iopt(21)) then 
            bad = .true.
            goto 20 
         end if 
c                                 get activity coefficients
         call mrkhyb (ins, jns, 5, 3, 1)
c                                 gi/T/T + ln P + ln gamma
         ogtot = gtot
         gtot = 0d0 

         do i = 1, 5

            eg(i) = dg(ins(i)) + dlog(g(ins(i)))
            gtot = gtot + y(ins(i))*(eg(i) + dlog(y(ins(i))))

         end do 
c                                 get newton-raphson increments
c                                 first derivatives of Gf
         d1  =   (eg(1) + dlog(y(ins(1))) + 1d0) 
     *         + (eg(3) + dlog(y(ins(3))) + 1d0)*a(1,1)
     *         + (eg(4) + dlog(y(ins(4))) + 1d0)*a(2,1)
     *         + (eg(5) + dlog(y(ins(5))) + 1d0)*a(3,1)
         d2  = + (eg(2) + dlog(y(ins(2))) + 1d0) 
     *         + (eg(3) + dlog(y(ins(3))) + 1d0)*a(1,2)  
     *         + (eg(4) + dlog(y(ins(4))) + 1d0)*a(2,2)  
     *         + (eg(5) + dlog(y(ins(5))) + 1d0)*a(3,2) 
c                                 second derivatives 
         d11 = 1d0/y(ins(1)) + a(1,1)*a(1,1)/y(ins(3)) 
     *                       + a(2,1)*a(2,1)/y(ins(4)) 
     *                       + a(3,1)*a(3,1)/y(ins(5))
         d22 = 1d0/y(ins(2)) + a(1,2)*a(1,2)/y(ins(3)) 
     *                       + a(2,2)*a(2,2)/y(ins(4)) 
     *                       + a(3,2)*a(3,2)/y(ins(5))
         d12 =                 a(1,1)*a(1,2)/y(ins(3))
     *                       + a(2,1)*a(2,2)/y(ins(4)) 
     *                       + a(3,1)*a(3,2)/y(ins(5))

         s = b(1)*y(ins(1)) + b(2)*y(ins(2)) + b(3)
c                                 complete negative first derivatives
         dgs1  = (gtot*b(1)/s - d1)/s
         dgs2  = (gtot*b(2)/s - d2)/s
c                                 complete second derivatives
         d11 = (d11 + 2d0*dgs1*b(1))/s 
         d22 = (d22 + 2d0*dgs2*b(2))/s 
         d12 = (((2d0*gtot*b(2)/s - d2)*b(1) - d1*b(2))/s + d12 )/s
c                                 determinant of the jacobian
         det = d11*d22 - d12*d12

         if (det.eq.0d0) then 
            write (*,*) 'dead'
            bad = .true.
            goto 20 
         end if 
c                                 newton-raphson increments for the ordered species
c                                 compositions. 
         dp(1) = (dgs1*d22 - dgs2*d12)/det 
         dp(2) = (dgs2*d11 - dgs1*d12)/det 

         itct = itct + 1

      end do 

      deltag = gtot*r*t/s

      end

      subroutine hybeos (jns, jsp)
c---------------------------------------------------------------------
c set up routine for hybrid fluid EoS calculations. computes the 
c (unecessay?) delta volumes and pure fluid fugacity coefficient 
c ratios used to convert the mrk fugacities to hybrid fugacities.

c the choice of the pure fluid eos are specified by the perplex_option
c keywords hybrid_EoS_H2O (iopt(24)), hybrid_EoS_CO2 (iopt(25)), and
c hybrid_EoS_H2O (iopt(26)).

c the routine mrkpur must be called prior to hybeos to set initial 
c guesses for volume. 

c 11/2016 JADC

c modified to return fh2o and fco2 in special fugacity cst11. 5/2017.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,jns(*),jsp

      double precision hsmrkf

      external hsmrkf
 
      double precision gh,vh,g0
      common/ csthyb /gh(nsp),vh(nsp),g0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)

      double precision vol
      common/ cst26 /vol
c----------------------------------------------------------------------
      do i = 1, jsp

         j = jns(i)

         vh(j) = -v(j)
         gh(j) = g(j)

         if (j.eq.1) then 
c                                 water
            if (iopt(25).eq.0) then 
c                                 mrk, already called
               f(j) = dlog(p*g(j))

            else if (iopt(25).eq.1) then 
c                                 hsmrk
               f(j) = hsmrkf (v(j),j)

            else if (iopt(25).eq.2) then
c                                 cork
               call crkh2o (p,t,v(j),f(j))

            else if (iopt(25).eq.4) then
c                                 pseos, pitzer & sterner 1994
               call pseos (v(j),f(j),j)

            else if (iopt(25).eq.5) then
c                                 haar, haar et el. 1982
               call haar (v(j),f(j))

            else if (iopt(25).eq.6) then 
c                                 zhang & duan 2005
               call zhdh2o (v(j),f(j))

            else if (iopt(25).eq.7) then 
c                                 zhang & duan 2009
               call zd09pr (v(j),f(j),j)

            end if 
c                                 the fugacity coefficient of the pure gas
            g(j) = dexp(f(j))/p
c                                 vol is set in case of a call by fluids for
c                                 the pure species. 
            vol = v(j)

         else if (j.eq.2) then 
c                                CO2
            if (iopt(26).eq.0) then 
c                                 mrk, already called
               f(j) = dlog(p*g(j))

            else if (iopt(25).eq.1) then 
c                                 hsmrk
               f(j) = hsmrkf (v(j),j)

            else if (iopt(26).eq.2) then
c                                 cork
               call crkco2 (p,t,v(j),f(j))

            else if (iopt(26).eq.3) then
c                                 brmrk, bottinga & richet 1981
               call brmrk (v(j),f(j))

            else if (iopt(26).eq.4) then
c                                 pseos, pitzer & sterner 1994
               call pseos (v(j),f(j),j)

            else if (iopt(26).eq.7) then 
c                                 zhang & duan 2009
               call zd09pr (v(j),f(j),j)

            end if 
c                                 the fugacity coefficient of the pure gas
            g(j) = dexp(f(j))/p

            vol = v(j)

         else if (j.eq.4) then
c                                CH4
            if (iopt(27).eq.0) then 
c                                 mrk, already called, do nothing
            else if (iopt(27).eq.1) then 
c                                 methane hsmrk kerrick and jacobs 1981.
               g(j) = dexp(hsmrkf (v(j),j))/p

            else if (iopt(27).eq.7) then 
c                                 zhang & duan 2009
               call zd09pr (v(j),f(j),j)

            end if

         end if
c                                 the hybrid delta volume (hyb-mrk), it's 
c                                 doubtful this thing is really used, if it 
c                                 is it must be in fluids.
         vh(j) = vh(j) + v(j)
c                                 the hybrid/mrk pure fluid fugacity ratio
         gh(j) = g(j)/gh(j)

       end do

      end

      double precision function lnfpur (id)
c---------------------------------------------------------------------
c function to return log fugacity of pure fluid of the nsp standard
c perple_x gas species, id is 100 + the index of the standard species. 
c the fugacity is calculated by the EoS indicated by the hyb_EoS
c option or, if no such option exists, with the MRK EoS. 

c 11/2016 JADC
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j, ins(1), id

      double precision hsmrkf,ftemp

      external hsmrkf

      double precision vol
      common/ cst26 /vol

      double precision gh,dvhy,gmrk0
      common/ csthyb /gh(nsp),dvhy(nsp),gmrk0(nsp) 

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)
c----------------------------------------------------------------------
      j = id - 100
      ins(1) = j
c                                call mrkpur in all cases to, at least, 
c                                get a volume guess
      call mrkpur (ins,1)
c                                save the pure mrk molar volume
      vmrk0(j) = vol

      ftemp = dlog(g(j)*p)
c                                for the generic hybrid eos solution model
c                                save the pure mix_EoS fugacity
      gmrk0(j) = g(j)
c                                then check if hyb_EoS indicates another
c                                EoS:
      if (j.eq.1) then 
c                                 water
         if (iopt(25).eq.1) then 
c                                 hsmrk
            ftemp = hsmrkf (v(j),j)

         else if (iopt(25).eq.2) then
c                                 cork
            call crkh2o (p,t,v(j),ftemp)

         else if (iopt(25).eq.4) then
c                                 pseos, pitzer & sterner 1994
            call pseos (v(j),ftemp,j)

         else if (iopt(25).eq.5) then
c                                 haar, haar et el. 1982
            call haar (v(j),ftemp)

         else if (iopt(25).eq.6) then 
c                                 zhang & duan 2005
            call zhdh2o (v(j),ftemp)

         else if (iopt(25).eq.7) then 
c                                 zhang & duan 2009
            call zd09pr (v(j),ftemp,1)

         end if

      else if (ins(1).eq.2) then 
c                                CO2
         if (iopt(26).eq.1) then
c                                 hsmrk
            ftemp = hsmrkf (v(j),j)

         else if (iopt(26).eq.2) then
c                                 cork
            call crkco2 (p,t,v(j),ftemp)

         else if (iopt(26).eq.3) then
c                                 brmrk, bottinga & richet 1981
            call brmrk (v(j),ftemp)

         else if (iopt(26).eq.4) then
c                                 pseos, pitzer & sterner 1994
            call pseos (v(j),ftemp,j)

         else if (iopt(26).eq.7) then
c                                 zhang & duan 2009
            call zd09pr (v(j),ftemp,1)

         end if

      else if (j.eq.4) then
c                                CH4
         if (iopt(27).eq.1) then 
c                                 methane hsmrk kerrick and jacobs 191.
            ftemp = hsmrkf (v(j),j)

         else if (iopt(27).eq.7) then 
c                                 zhang & duan 2009
            call zd09pr (v(j),ftemp,1)

         end if

      end if 

      lnfpur = ftemp
c                                 set bulk v variable for frendly
c                                 dielectric constant calculation, 
c                                 hopefully this doesn't interfere 
c                                 with anything.
c                                           JADC, march 6 2018
      vol = v(j)
      dvhy(j) = vol - vmrk0(j)

      end



      subroutine zd09pr (vol,lnfug,i)
c----------------------------------------------------------------------
c Zhand & Duan 2009 pure fluid EoS for species as below.
c See: duan2009_CO2_volume_fugacity.mws

c         1 = H2O
c         2 = CO2
c         3 = CO
c         4 = CH4 
c         5 = H2
c         7 = O2
c        16 = C2H6
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer it, iwarn, ins(1), i

      double precision prt,b,c,d,e,f,ge,expg,gamm,vi,veq,dveq,lnfug,dv

      double precision vol, vmrk, fmrk, eps(16), sig3(16), et, et2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      save iwarn, eps, sig3
      data iwarn/0/

      data eps /510d0, 235d0, 105.6d0, 154d0, 31.2d0, 0d0, 124.5d0,
     *          8*0d0, 246.1d0/

      data sig3 / 23.887872d0, 54.439939d0, 49.027896d0, 50.28426837d0,
     *            25.153757d0, 0d0, 37.933056d0, 8*0d0, 82.312875d0/
c---------------------------------------------------------------------
c                                 MRK volume guess and backup fugacity
      ins(1) = i
      call mrkpur (ins, 1)

      vmrk = v(i)
      fmrk = dlog(p*g(i))

      vol = vmrk
c                                 zd pv units MPa, J/MPa (cm3).
      prt = p/1d1/r/t
      it = 0

      gamm = 6.123507682d0*sig3(i)**2
      et = eps(i)/t
      et2 = et * et

      b = (.5870171892d0 + (-5.314333643d0 -1.498847241d0*et)*et2) *
     *    sig3(i)
      c = (.5106889412d0 + (-2.431331151d0 + 8.294070444d0*et)*et2) * 
     *    sig3(i)**2
      d = (.4045789083d0 + (3.437865241d0 - 5.988792021d0*et)*et2) * 
     *    sig3(i)**4
      e = (-.7351354702d-1 + (.7017349038d0 - .2308963611d0*et)*et2) *
     *    sig3(i)**5
      f = 1.985438372d0*et2*et*sig3(i)**2
      ge = 16.60301885d0*et2*et*sig3(i)**4 
c                                 iteration loop for volume
      do 

         vi = 1d0/vol

         expg = dexp(-gamm*vi*vi)
c                                 p(v)/rt
         veq = -vi - b*vi**2 + (-f*expg-c)*vi**3 + (-ge*expg-d)*vi**5 
     *             - e*vi**6
c                                 diff(veq,v)
         dveq = -veq*vi + b*vi**3 + 2d0*(f*expg+c)*vi**4  
     *           + (-2d0*f*expg*gamm + 4d0*ge*expg + 4d0*d)*vi**6 
     *           + 5d0*e*vi**7 - 2d0*ge*expg*gamm*vi**8

         dv = -(prt + veq)/dveq

         if (dv.lt.0d0.and.vol+dv.lt.0d0) then 

            vol = vol*0.8d0

         else 

            vol = vol + dv

         end if 

         if (dabs(dv/vol).lt.nopt(50)) then

            exit
          
         else if (vol.lt.0d0.or.it.gt.iopt(21)) then
c                                 use cork fugacities
            iwarn = iwarn + 1

            if (iwarn.le.50) then 
               write (*,1000) p,t,vol
               if (iwarn.eq.50) call warn (49,p,93,'zh09pr')
            end if 

            exit 

         end if 

         it = it + 1

      end do

      expg = dexp(gamm/vol/vol)

      lnfug = dlog(r*t/vol/pr/1d-1) 
     *        + 0.5d0*(f+ge/gamm)*(1d0-1d0/expg)/gamm
     *        + (2d0*b + (1.5d0*c + (f-0.5d0*ge/gamm)/expg 
     *        + (1.25d0*d + ge/expg + 1.2d0*e/vol)/vol**2)/vol)/vol

c                                 convert volume from j/bar to cm3/mol
      vol = 1d1 * vol

1000  format (/,'**warning ver093** zd09pr did not converge at:',
     *        3(1x,g12.6))

      end 

      subroutine zhdh2o (v,lnfug)
c----------------------------------------------------------------------
c Zhang & Duan 2005 pure water volume and fugacity of water. 
c see duan_water_volume_fugacity.mws
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision prt,b,c,d,e,f,g,expg,gamm,v,vi,veq,dveq,lnfug,dv,
     *                 fcrk, vcrk

      integer it, iwarn

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      save iwarn
      data iwarn/0/
c---------------------------------------------------------------------
c                                 CORK volume guess and backup fugacity
      call crkh2o (p,t,vcrk,lnfug)

      v = vcrk /1d1
      fcrk = lnfug

      prt = p/r/t
      it = 0

      gamm = .3317993788d0

      b = 1.957197778 - 6821674.863d0/t**2 + 3047984261d0/t**3
      c = 3.531471196 + 9821873.173d0/t**2 - 7411448875d0/t**3
      d = 16.71639581 - 6007496.747d0/t**2 + .1540316803d11/t**3
      e = -4.611555959 + 11372008.36d0/t**2 - .136192675d11/t**3
      f = -2033.267066d0 / t
      g = -0.002765323035d0 * t
c                                 iteration loop for volume
      do 

         vi = 1d0/v

         expg = dexp(-gamm/v/v)
c                                 p(v)/rt
         veq = -vi - b*vi**2 + (-f*expg-c)*vi**3 + (-g*expg-d)*vi**5 
     *             - e*vi**6
c                                 diff(veq,v)
         dveq = -veq*vi + b*vi**3 + 2d0*(f*expg+c)*vi**4  
     *           + (-2d0*f*expg*gamm + 4d0*g*expg + 4d0*d)*vi**6 
     *           + 5d0*e*vi**7 - 2d0*g*expg*gamm*vi**8

         dv = -(prt + veq)/dveq

         if (dv.lt.0d0.and.v+dv.lt.0d0) then 

            v = v*0.8d0

         else 

            v = v + dv

         end if 

         if (dabs(dv/v).lt.nopt(50)) then

            exit
          
         else if (v.lt.0d0.or.it.gt.iopt(21)) then
c                                 use cork fugacities
            iwarn = iwarn + 1

            if (iwarn.le.50) then 
               write (*,1000) p,t,v
               if (iwarn.eq.50) call warn (49,p,93,'ZHDH2O')
            end if

            v = vcrk

            lnfug = fcrk

            return

         end if 

         it = it + 1

      end do

      expg = dexp(gamm/v/v)

      lnfug = dlog(r*t/v) + 0.5d0*(f+g/gamm)*(1d0-1d0/expg)/gamm
     *        +(2d0*b + (1.5d0*c + (f-0.5d0*g/gamm)/expg 
     *        + (1.25d0*d + g/expg + 1.2d0*e/v)/v**2)/v)/v
c                                 convert volume from j/bar to cm3/mol
      v = 1d1 * v

1000  format (/,'**warning ver093** ZHDH2O did not converge at:',
     *        3(1x,g12.6))

      end 