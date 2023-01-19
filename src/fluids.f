      program cohsrk
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *    cohsrk.6.1993     *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------

c COHSRK is a fortran program to call various subroutines to calculate  
c C-O-H-S fluid speciation as a function of X(O), X(CO2), X(S), X(C) 
c a(C), f(S2), f(O2) or O2-S2-buffer assemblage. 

c The variables X(O), X(S), and X(C) are defined as:

c  X(O) = n(O)/{n(O)+n(H)}
c  X(S) = n(S)/{n(S)+n(C)}
c  X(C) = n(C)/{n(S)+n(C)+n(O)+n(H)}

c where n(C), n(S), n(O) and n(H)  are the total number of moles of 
c carbon, sulfur, oxygen and hydrogen (as opposed to the amounts
c of the species) in the fluid.

c The user may choose from the following routines, identified by
c number:

c    0 - Modified Redlich-Kwong (MRK)
c    1 - Kerrick & Jacobs 1981 hard sphere MRK (HSMRK)
c    2 - Hybrid MRK/HSMRK
c    3 - Saxena & Fei 1987 pseudo-virial expansion
c    4 - Bottinga & Richet 1981 (CO2 RK)
c    5 - Holland & Powell 1990 (CORK)
c    6 - Hybrid Haar et al 1979/HSMRK (TRKMRK)
c    7 - Graphite buffered COH MRK fluid
c    8 - Graphite buffered COH hybrid-RK fluid
c    9 - Maximum X(H2O) GCOH fluid Cesare & Connolly 1993
c   10 - X(O) GCOH-fluid hybrid-MRK Connolly & Cesare 1993
c   11 - X(O) GCOH-fluid MRK Connolly and Cesare 1993
c   12 - X(O) GCOHS Connolly & Cesare 1993
c   13 - X(H) H2-H2O-hybrid
c   14 - hogbrd, don't use this if you dont know what it is.
c   15 - X(H) low T H2-H2O-hybrid
c   16 - X(O) H-O HSMRK/MRK hybrid
c   17 - X(O) H-O-S HSMRK/MRK hybrid
c   18 - Delany/HSMRK/MRK hybrid, for P > 10 kb
c   19 - X(O)-X(S) GCOHS Connolly & Cesare 1993
c   20 - X(O)-X(C) GCOHS Connolly & Cesare 1993
c   24 - f(O2)-N/C graphite saturated COHN MRK fluid

c Routines 0-3, 5-6, and 18 are for conventional P-T-X(CO2) 
c calculations, where X(CO2) is the mole fraction of CO2 in 
c a binary H2O-CO2 mixture. The Saxena & Fei routine was
c programmed by someone else, and I suspect it has been 
c entered incorrectly. 

c Routine 4 is for pure CO2 fluids.

c Routines 7 and 8 are for C-O-H fluids as a function of 
c P-T and f(O2) or f(O2)-buffer at specified 
c graphite activity (usually 1). I recommend choice 8.

c Routine 9 returns the fugacities of O2, H2O, and CO2
c for a H:O = 2 graphite saturated C-O-H fluid, i.e.,
c the fugacities computed with routine 10 at X(O) = 1/3
c and a(graphite) = 1.

c Routines 10 and 11 are for C-O-H fluids as a function 
c of X(O) at specified graphite activity. I recommend 10.

c Routines 12 and 17 are for C-O-H-S and H-O-S fluids as
c a function of P-T-X(O) conditions at specified sulfur
c fugacity or buffer, for C-O-H-S fluids graphite activity
c must also be specified. 

c Routines 13 and 15 are for binary H2/H2O mixtures 
c (effectively H-O fluids at X(O) < 1/3), routine 15 uses
c parameters optimized for the solvus region.

c Routine 16 is for H-O fluids as a function of P-T-X(O).

c Routines 2, 8, 12, 13, 14, 16 and 17 are a hybrids of the 
c HSMRK EoS (Kerrick & Jacobs) and the MRK (Redlich & Kwong/
c DeSantis et al/Holloway) EoS as described by Connolly & 
c Cesare, J Met Geol, 11:379-388, and for most purposes I would
c recommend these equations.

c Routine 19 calculates C-O-H-S graphite-UNDERSATURATED fluid 
c properties as a function of X(O) and the atomic S/C ratio 
c (expressed by the variable X(S)) at specified sulfur fugacity. 
c This routine often may not converge. The EoS is identical to
c that used in Routine 12.

c Routine 20 calculates C-O-H-S graphite-UNDERSATURATED fluid 
c properties as a function of X(O) and the atomic fraction of
c carbon in the fluid (X(C)) at specified sulfur fugacity.
c This routine can be made to calculate simple COH fluid 
c speciation as a function of bulk composition by setting 
c sulfur fugacity to such a low value that the concentration of 
c sulfur becomes negligible.
c This routine often may not converge. The EoS is identical to
c that used in Routine 12.

c Routine 24 calculates C-O-H-N graphite saturated fluid properties
c as a function of f(O2) and the fluid N/C ratio using a modified
c Redlich-Kwong EoS. The routine will fail in the limits x(H2O)->0
c and N/C->0.

c COHSRK is primarily intended to provide an example of how the
c various speciation routines in PERPLEX can be called. If you
c can follow the logic of the program then you shall probably be
c able to customize the program to use input, and generate output,
c more in line with your specific needs. I do have another program
c that generates a plot file with various species, concentrations,
c and fugacities plotted against X(O), and I shall be happy to 
c provide it to interested users.

c----------------------------------------------------------------------

c Limitations/warnings:

c The C-O-H-S and H-O-S routines do not consider the species C2H6,
c O2, and SO3 and pure S species, these will not be significant at
c the conditions realized in crustal metamorphic environments. 
c It is relatively easy to incorporate new species in the routines
c so if you suspect the species are important, please contact me
c and I will add the species.

c The HSMRK equation of state for water seems unreliable at P > 
c 15 Kbar, and hence all the hybrid EoS speciation routines
c as well. This has not bothered me because I am mostly concerned  
c with crustal metamorphic conditions; however, should you wish
c to make calculations at high pressure you need only replace the
c calls which get pure fluid fugacities with calls to the high
c pressure routine of your choice (e.g., DeSantis-Holloway MRK,
c Delany & Helgeson, Bottinga & Richet, etc. etc.).

c Functions for equilibrium constants are fit in the range 400-1400K.

c No warnings are issued if P-T-X(O) conditions are within a solvus.

c All X(O) routines fail at X(O) = 0 or 1, to avoid this problem
c the routines reset extreme values of X(O) to 1d-10 or 0.9999999999
c as appropriate.

c All X(O) routines may fail at X(O) = 1/3 if the fluid becomes an
c essentially mono-species H2O fluid, this is most likely to occur
c in simple H-O fluids.

c The C-O-H-S and H-O-S routines permit calculations for a fluid
c in equilibrium with pyrrhotite with a fixed atomic Fe/S ratio,
c this is not the same as N(FeS), the mole fraction of troilite
c relative to S2, as often reported for pyrrhotite analyses.

c Routines 19 and 20 may not converge, or may converge to 
c invalid roots, in the vicinity of graphite saturated and
c supersaturated conditions; and at very carbon under-saturated 
c conditions numerical slop may lead to errors on the order of
c 0.4 log units in calculated properties such as f(O2).

c In speciation routines, if the concentration of a species becomes 
c negligible (as determined by the numerical precision of the
c computer being used), its properties may be undefined or
c may have meaningless fluctuations.
c----------------------------------------------------------------------

c NOTES ON THIS SOURCE AND COMPILING:

c The MAIN program COHSRK is included in this file,
c all the subroutines called by COHSRK are in the file flib.f.
c The files may be compiled separately and the objects linked together,
c or you may concatenate the two files. 

c This program was put together from routines used in PERPLEX (see
c below) and for this reason it is not as compact or as flexible
c as might be desired. The size of the flib.f may cause problems
c for DOS compilers, these can be overcome by splitting the source
c into two or three blocks that can then be combined during linking.
c Alternatively, routines that are not considered may be useful
c can be eliminated if the calls to these routines from subprogram
c cfluid are also eliminated.

c Specifically I would recommend eliminating:

c      subroutine brmrk
c      subroutine simps 
c      subroutine qromb 
c      subroutine polint 
c      subroutine trapzd 
c      function brvol 
c      function vdpdv 
c      subroutine hosrk5 
c      subroutine cohfit 
c      subroutine haar
c      subroutine psat2
c      subroutine aideal 
c      subroutine trkmrk
c      subroutine saxfei
c      subroutine hprk
c      subroutine cohgra 
c      subroutine hh2ork 
c      subroutine lohork
c      subroutine lomrk 

c this requires elimination of the calls to: brmrk, cohfit,
c trkmrk, saxfei, hprk, cohgra, hh2ork, and lohork from cfluid.

c subroutines warn and error can also be eliminated if the 
c calls to these routines, from subroutines rfluid, brmrk,
c cohgra, cohsgr and cohhyb, are replaced by statements that write 
c an appropriate warning or error message. 

c----------------------------------------------------------------------

c These routines are incorporated in the PERPLEX programs for calculating
c phase equilibria and diagrams, if you would like a copy of these 
c programs, or if you have any problems with this program, please
c contact me at:

c                     James Connolly
c                     IGP-ETHZ
c                     CH-8092 Zuerich

c by e-mail at:

c                     james.connolly@erdw.ethz.ch

c or by telephone/fax at:

c                     0041-44-632-7804/0041-44-632-1088

c-----------------------------------------------------------------------
c I/O: most input and output is done through the common blocks below,
c the significance of the variables as named in the main program is:

c ifug   - number indexing the requested EoS.
c p      - pressure, bars.
c t      - temperature, kelvins.
c xo     - X(O) for multispecies routines, and X(CO2) or X(H2) for
c          binary routines.
c vol    - molar volume for all multispecies routines, and binary
c          routines 0, 1, 13, 15.
c fhc(1) - natural log (f(H2O)) for all routines.
c fhc(2)   - natural log of the species other than H2O (i.e., CO2 or H2)
c          in all binary routines.
c fo2    - natural log (f(O2)).
c fs2    - 1/2 natural log (f(S2)).

c the following variables are only used for multispecies calculations:

c ins(i) - pointers indicating the species are to be calculated.
c isp    - number of species to be calculated.
c nsp    - dimensioning for the maximum number of species. 
c ibuf   - pointer indicating method of calculating f(O2) for
c          routines 7 and 8, or f(S2) for routines 12 and 17.
c dlnfo2 - for routines 7 and 8 the displacement of the f(O2)
c          relative to a buffer (in log units) or the absolute
c          ln(f(O2)) (if ibuf = 3). For routines 12 and 17,
c          if ibuf = 2 then dlnfo2 is the atomic Fe/S of pyrrhotite,
c          if ibuf = 3 then dlnfo2 is 1/2 the natural log (f(S2)). 
c elag   - natural log of graphite activity.
c g(i)   - fugacity coefficient of the ith species.
c x(i)   - mole fraction of the ith species.

c the indices of eleven species presently defined are:

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

c O2 should be replaced by SO3.
c-----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      logical log, bad, fileio

      character y*1, n4name*100, title*100, tags(33)*14

      integer ier, igo, i, isp, j, k, l, kmax, kount, nel

      double precision nc, nh, no, ns, nn, nsi, tentoe, fo2, fs2, fh2,
     *                 ag, tot, totx, var(l2), f, prop(40), vdif,
     *                 vpar(nsp), xxs(nsp), xg(nsp)

      double precision fhc
      common / cst11 /fhc(3)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer inc,jpot
      common/ cst101 /inc(l2),jpot

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision p,t,xo,u
      common/ cst5 /p,t,xo,u(6)

      double precision xs,g,v
      common/ cstcoh /xs(nsp),g(nsp),v(nsp)

      character specie*4
      integer ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      double precision vol
      common/ cst26 /vol

      integer ibuf,hu,hv,hw,hx 
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      integer iam
      common/ cst4 /iam

      double precision eqk
      common / csteqk /eqk(nsp)

       save tentoe
       data tentoe/2.302585093d0/
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 11
      vname(1) = 'P(bar)'
      vname(2) = 'T(K)'
c                                 max number of indendent potentials, may be
c                                 increased depending on EoS choices. 
      ipot = 2
c                                 max number of element fractions to be output
      nel = 6
c                                 version info
      call vrsion (6)
c                                 read options
      n4name = 'perplex_option.dat'
      call redop1 (.false.,n4name)

      open (39,file='partial_molar_volumes.dat')

      do 
c                                 configure EoS
         igo = 0
         fs2 = -9999d0*tentoe/2d0
         elag = 0d0
c                                 get the users choice of EoS:   
         call rfluid (1)
c                                 for multispecies fluids set
c                                 up species indices and name
c                                 of independent variable
         call setins (ifug)

         write (*,'(/,a)') 'Tabulate properties (y/n)?'
         read (*,'(a)') y

         if (y.eq.'y'.or.y.eq.'Y') then 
c                                 tabulated properties
            if (ifug.le.3.or.ifug.eq.13..or.ifug.eq.14.or.
     *               ifug.eq.17.or.ifug.eq.25) then 
c                                 binary x EoS's
               ipot = 3

            else if (ifug.lt.9) then
c                                 fo2 EoS's
               if (ibuf.eq.3) then
                  ipot = 4
                  vname(4) = 'log(aC)'
               else 
                  ipot = 3    
                  vname(3) = 'log(aC)'
               end if
               
            else if (ifug.le.11) then
c                                 X(O)-a(C) EoS
               ipot = 4
               vname(4) = 'log(aC)'

            else if (ifug.eq.12) then 
c                                 X(O)-a(C)-f(S2) EoS
               ipot = 4
               vname(4) = 'log(aC)'

               if (ibuf.eq.3) then 
                  ipot = 5
                  vname(5) = 'log(fS2)'
               end if 

            else if (ifug.eq.15.or.ifug.eq.16) then 
c                                X(O) HO
               ipot = 3

            else if (ifug.eq.17) then 
c                                X(O)-fs2 HOS 
               ipot = 3

               if (ibuf.eq.3) then 
                  ipot = 4
                  vname(4) = 'log(fS2)'
               end if                   

            else if (ifug.eq.19) then 
c                                X(O)-X(S) COHS
               ipot = 4    
               vname(4) = 'X(S)    '

            else if (ifug.eq.20) then 
c                                X(O)-X(C) COH
               ipot = 4
               vname(4) = 'X(C)    '

               if (ibuf.eq.3) then 
                  ipot = 5
                  vname(5) = 'log(fS2)'
               end if 

            else if (ifug.eq.24) then 
c                                fo2-aC-N/C COHN
               vname(3) = 'log(aC) '
               if (ibuf.eq.3) then
                  ipot = 5 
                  vname(4) = 'log(fO2)'
                  vname(5) = 'N/C     '          
               else 
                  ipot = 4 
                  vname(4) = 'N/C     '
               end if 

            else if (ifug.eq.26) then 

               ipot = 3

            else if (ifug.eq.27) then 
c                                X(O)-X(C) COH
               ipot = 4
               vname(4) = 'X(C)    '

            end if 

            do 

               write (*,'(/,a)') 
     *               'How many independent variables, 1 or 2?'
               call rdnumb (p,0d0,jpot,1,.false.)

               if (jpot.eq.5) then 
c                                 hidden file p-t input option
                  jpot = 1
                  fileio = .true.
c                                 special file output
                  do 
                     write (*,*) 'Enter P-T file name '
                     read (*,*) title
                     open (666,file=title,status='old',iostat=ier)
                     if (ier.ne.0) cycle
                     exit
                  end do 

                  exit

               else

                  fileio = .false.

               end if 

               if (jpot.gt.0.and.jpot.lt.3) exit

            end do 
            
            if (jpot.eq.1) then 
               write (*,'(/,a)') 'Select the independent variable:'
            else 
               write (*,'(/,a)') 
     *               'Select the primary independent variable:'
            end if 

            do i = 1, ipot
               write (*,'(3x,i1,a,a)') i,' - ',vname(i)
               iv(i) = 0
            end do 

            do 
               read (*,*,iostat=ier) iv(1)
               if (ier.eq.0.and.iv(1).gt.0.and.iv(1).le.ipot) exit
               call rerr
            end do 

            if (jpot.eq.2) then 

               write (*,'(/,a)') 
     *               'Select the secondary independent variable:'
               do i = 1, ipot
                  if (i.eq.iv(1)) cycle 
                  write (*,'(3x,i1,a,a)') i,' - ',vname(i)
               end do 

               do 
                  read (*,*,iostat=ier) iv(2)
                  if (ier.eq.0.and.iv(2).ne.iv(1).and.
     *                iv(2).gt.0.and.iv(2).le.ipot) exit
                  call rerr
               end do 

            end if 

            write (*,'(/)')
c                                 get independent variable range and increments
            do i = 1, jpot

               write (*,'(a,a)') 
     *               'Enter minimum, maximum and increment for: ', 
     *               vname(iv(i))
               read (*,*,iostat=ier) vmin(iv(i)),vmax(iv(i)),dv(iv(i))

               var(iv(i)) = vmin(iv(i))

            end do 
c                                 get sectioning values for the remainder
            j = jpot

            do i = 1, ipot

               if (jpot.eq.1.and.iv(1).eq.i) then
                  cycle 
               else if (iv(1).eq.i.or.iv(2).eq.i) then 
                  cycle
               end if 

               j = j + 1
               iv(j) = i

               do 

                  write (*,'(a,a)') 'Enter the value for: ',
     *                               vname(i)
                  read (*,*,iostat=ier) var(i)
                  vmin(i) = var(i)
                  if (ier.eq.0) exit 
                  call rerr

               end do 
                              
            end do

            do i = 1, jpot
               inc(iv(i)) = 
     *            idint(dabs(vmax(iv(i))-vmin(iv(i)))/dv(iv(i))) + 1
            end do 

            do i = jpot+1, ipot
               inc(iv(i)) = 1
               dv(iv(i)) = 1
            end do 
c                                 query output file name
            write (*,'(/,a,a)') 'Enter a name for the output file ',
     *                          '[without the .tab suffix]:'
c                                 readrt loads the root into prject
            call readrt 
            call mertxt (n4name,prject,'.tab',0)
            open (n4,file=n4name)
c                                 query for title
            write (*,'(/,a)') 'Enter calculation title:'
            read (*,'(a)') title

            write (*,'(/,a,a)') 'Output logarithmic values for ',
     *           'species fractions and fugacities (y/n)?'
            read (*,'(a)') y

            if (y.eq.'y'.or.y.eq.'Y') then 
               log = .true.
            else 
               log = .false.
            end if 
c                                 write header
c                                 write version flag
            write (n4,'(a)') '|6.6.6'
            write (n4,'(a)') title 
            write (n4,*) jpot

            do i = 1, jpot
               write (n4,*) vname(iv(i))
               write (n4,*) vmin(iv(i))
               write (n4,*) dv(iv(i))
               write (n4,*) inc(iv(i))
            end do 
c                                 create column tags
            do i = 1, ipot
               tags(i) = vname(iv(i))
            end do 

            tags(ipot+1) = 'vol[cm3/mol]'
c                                 species fractions
            j = 0
            do i = ipot+2, ipot+isp+1

               j = j + 1
               if (log) then 
                  write (tags(i),'(a,a,a)') 'log[y(',specie(ins(j)),')]'
               else 
                  write (tags(i),'(a,a,a)') 'y(',specie(ins(j)),')'
               end if 
               call unblnk (tags(i))

            end do
c                                 atomic fractions
            tags(ipot+isp+2) = 'Y_C'
            tags(ipot+isp+3) = 'Y_O'
            tags(ipot+isp+4) = 'Y_H'
            tags(ipot+isp+5) = 'Y_N'
            tags(ipot+isp+6) = 'Y_S'
            tags(ipot+isp+7) = 'Y_Si'

            kount = ipot+isp+8
c                                  species fugacities
            j = 0
            do i = kount, kount+isp-1
               j = j + 1 
               if (log) then 
                  write (tags(i),'(a,a,a)') 'log[f(',specie(ins(j)),')]'
               else 
                  write (tags(i),'(a,a,a)') 'f(',specie(ins(j)),')'
               end if 
               call unblnk (tags(i))
            end do

            kount = kount + isp + 3
            tags(kount-3) = 'log[f(O2)]'
            tags(kount-2) = 'log[f(S2)]'
            tags(kount-1) = 'log[f(H2)]'
            tags(kount)   = 'log[a(C)] '

            write (n4,*) kount 
            write (n4,'(40(a14,1x))') (tags(i),i=1,kount)
c                                 terminal info on variables
            write (*,'(/,a,/)') 'Table columns will be:'
            write (*,'(5(a14,1x))') (tags(i),i=1,kount)
         
            if (jpot.eq.1) then 
               call plblrb (4)
            else
               call plblrb (1)
            end if 
c                                 computational loop, initialize p,t to be safe
            p = var(1)
            t = var(2) 

            if (ifug.le.5 .or.ifug.eq.13.or.ifug.eq.14.or.
     *          ifug.eq.15.or.ifug.eq.16.or.ifug.eq.25.or.
     *          ifug.eq.26) then 
c                                 xco2, xo, xh2 EoS's
                xo = var(3)

            else if (ifug.eq.8) then
c                                 fo2 EoS's
               if (ibuf.eq.3) then
                  dlnfo2 = var(3) * tentoe
                  elag = var(4) * tentoe
               else 
                  elag = var(3) * tentoe
               end if
               
            else if (ifug.eq.10) then
c                                 X(O)-a(C) EoS
               xo = var(3)
               elag = var(4) * tentoe

            else if (ifug.eq.12) then 
c                                 X(O)-a(C)-f(S2) EoS
               xo = var(3)
               elag = var(4) * tentoe
               if (ibuf.eq.3) dlnfo2 = var(5) * tentoe

            else if (ifug.eq.17) then 
c                                X(O)-fs2 HOS 
               xo = var(3)
               if (ibuf.eq.3) dlnfo2 = var(4) * tentoe               

            else if (ifug.eq.19) then 
c                                X(O)-X(S) COHS
               xo = var(3)
               gz = var(4)  

            else if (ifug.eq.20) then 
c                                X(O)-X(C) COH
               xo = var(3)
               gz = var(4)  
               if (ibuf.eq.3) dlnfo2 = var(5) * tentoe

            else if (ifug.eq.24) then 
c                                fo2-aC-N/C COHN
               elag = var(3) * tentoe
               if (ibuf.eq.3) then
                  dlnfo2 = var(4) * tentoe
                  gz = var(5)            
               else 
                  gz = var(4) 
               end if 

            else if (ifug.eq.27) then 
c                                 C-O-H XO-YC
               xo  = var(3)
               fs2 = var(4) 

            end if 

            do j = 1, inc(iv(2))

               var(iv(2)) = vmin(iv(2)) + dfloat(j-1)*dv(iv(2))

               do i = 1, inc(iv(1))

                  var(iv(1)) = vmin(iv(1)) + dfloat(i-1)*dv(iv(1))

                  do k = 1, jpot
c                                 assign values to local variables
                     if (iv(k).eq.1) then 
                        p = var(iv(k))
                     else if (iv(k).eq.2) then 
                        t = var(iv(k))
                     else if (iv(k).eq.3) then 

                        if (ifug.eq.8) then

                           if (ibuf.eq.3) then
                              dlnfo2 = var(iv(k)) * tentoe
                           else 
                              elag = var(iv(k)) * tentoe
                           end if

                        else if (ifug.eq.24) then 

                           elag = var(iv(k)) * tentoe

                        else 

                           xo = var(iv(k))

                        end if 

                     else if (iv(k).eq.4) then 

                        if (ifug.eq.7.or.ifug.eq.8) then
 
                           if (ibuf.eq.3) elag = var(iv(k)) * tentoe
               
                        else if (ifug.eq.10.or.ifug.eq.11.or.
     *                          ifug.eq.12.) then

                           elag = var(iv(k)) * tentoe

                        else if (ifug.eq.17) then 

                           if (ibuf.eq.3) dlnfo2 = var(iv(k)) * tentoe

                        else if (ifug.eq.19.or.ifug.eq.20) then 

                           gz = var(iv(k))  

                        else if (ifug.eq.24) then 

                           if (ibuf.eq.3) then
                              dlnfo2 = var(iv(k)) * tentoe        
                           else 
                              gz = var(iv(k)) 
                           end if 

                        else if (ifug.eq.27) then 

                           fs2 = var(iv(k)) 

                        else 

                           write (*,*) 'fugga wugga!'
                           call errpau

                        end if 

                     else if (iv(k).eq.5) then 

                        if ((ifug.eq.12.or.ifug.eq.20)
     *                                    .and.ibuf.eq.3) then 

                           dlnfo2 = var(iv(k)) * tentoe

                        else if (ifug.eq.24.and.ibuf.eq.3) then 

                           gz = var(iv(k))             

                        else 

                           write (*,*) 'fugga wugga!'
                           call errpau

                        end if
 
                     end if 

                  end do 

c                                 calculate properties
                  bad = .false.

                  if (fileio) read (666,*,iostat=ier) p, t

                  call cfluid (fo2, fs2)    
c                                 variables
                  do k = 1, ipot 
                     prop(k) = var(iv(k))
                  end do 
c                                 species fractions/fugacities
                  if (ifug.le.5 .or.ifug.eq.14.or.ifug.eq.21.or.
     *                ifug.eq.22.or.ifug.eq.25) then 
c                                 xco2 EoS's 
                     xs(ins(1)) = 1d0 - xo
                     xs(ins(2)) = xo

                     do k = 1, 2

                        f = dexp(fhc(k))

                        if (log) then 
                           prop(ipot+1+k) = dlog10(xs(ins(k)))
                           prop(ipot+isp+1+nel+k) = dlog10(f)
                        else 
                           prop(ipot+1+k) = xs(k)
                           prop(ipot+isp+1+nel+k) = f
                        end if

                        if (dabs(prop(ipot+isp+1+nel+k)).gt.1d99) 
     *                           prop(ipot+isp+1+nel+k) = nopt(7) 

                     end do

                  else 
c                                 assume multispecies fluids                
                     do k = 1, isp

                        f = xs(ins(k))*p*g(ins(k))

                        if (log.and.f.le.0d0) then  
                           prop(ipot+1+k) = nopt(7)
                           prop(ipot+isp+1+nel+k) = nopt(7)
                        else if (log) then 
                           prop(ipot+1+k) = dlog10(xs(ins(k)))
                           prop(ipot+isp+1+nel+k) = dlog10(f)
                        else 
                           prop(ipot+1+k) = xs(ins(k))
                           prop(ipot+isp+1+nel+k) = f
                        end if

                        if (dabs(prop(ipot+isp+1+nel+k)).gt.1d99) 
     *                           prop(ipot+isp+1+nel+k) = nopt(7) 

                     end do

                  end if 
c                                 atomic fractions 
                  call elmnts (nc,no,nh,nn,ns,nsi)

                  tot = 0d0 

                  do k = 1, isp
                     tot = tot + xs(ins(k))
                  end do 

                  if (dabs(tot-1d0).gt.nopt(50)+0.5d0) then
                     write (*,*) ' bad total :',tot-1d0,xo,elag
                     bad = .true. 
                  end if 

                  tot = ns + no + nh + nc + nn + nsi

                  prop(ipot+isp+2) = nc/tot
                  prop(ipot+isp+3) = no/tot
                  prop(ipot+isp+4) = nh/tot
c                                 ternary coordinates (XC = Y, XO = X)
c                 prop(ipot+isp+2) = nc/tot * 0.866025d0
c                 prop(ipot+isp+3) = (no + nc/2d0)/tot
c
                  prop(ipot+isp+5) = nn/tot
                  prop(ipot+isp+6) = ns/tot
                  prop(ipot+isp+7) = nsi/tot

                  prop(kount-3) = fo2/tentoe
                  prop(kount-2) = fs2/tentoe

                  if (ifug.eq.27) then 
c                                   ac, fo2, fh2 computed from major
c                                   species
                     prop(kount-3) = fhc(2)/tentoe
                     prop(kount-2) = nopt(7)
                     prop(kount-1) = fhc(3)/tentoe
                     prop(kount)   = fhc(1)/tentoe

                  else 

                     prop(kount-3) = fo2/tentoe
                     prop(kount-2) = fs2/tentoe
                     prop(kount-1) = nopt(7)
                     prop(kount)   = elag/tentoe

                  end if 

                  if (ifug.le.2.or.ifug.eq.13.or.ifug.eq.15) then
c                                 use analytic vol
                  else 
c                                 compute volume by finite difference
                     vdif = 0d0 

                     do k = 1, isp
                        vpar(k) = 0d0
                        xxs(ins(k)) = xs(ins(k))
                     end do 

                     p = var(1) + 0.5d0
                     f = 1d0 

                     do l = 1, 2

                        call cfluid (fo2,fs2)

                        do k = 1, isp
                           if (g(ins(k))*p*xxs(ins(k)).eq.0d0) cycle
                           vpar(k) = vpar(k) +  
     *                         83.14d0*t*f*dlog(g(ins(k))*p*xxs(ins(k)))
                           vdif = vdif + 
     *                       f*xxs(ins(k))*dlog(g(ins(k))*p*xxs(ins(k)))
                        end do 

                        f = -1d0
                        p = var(1) - 0.5

                     end do 

                     p = var(1)
c                                 use finite difference total volume only if non-hybrid EoS
                    if (ifug.eq.5.or.ifug.eq.14.or.ifug.eq.25) 
     *                 vol = 83.14d0*t*vdif

                     write (39,'(12(g14.6,1x))') p, t, (vpar(k),k=1,isp)

                  end if 

                  prop(ipot+1) = vol

                  if (bad) then 
                     do k = ipot+1, kount
                        prop(k) = nopt(7)
                     end do
                  end if 

                  write (n4,'(40(g14.7,1x))') (prop(k),k=1,kount)
     
               end do 

            end do 
    
            write (*,'(a,a)') 'The table has been written to file: ',
     *                        n4name
            close (n4)

         else  
c                                 properties at arbitrary conditions        
            write (*,'(/,a,/)') 'Enter a zero for pressure to quit.'

            do 
c                                 just a check
               vol = 0d0

               if ((igo.eq.0.or.ibuf.ne.3).and.ifug.eq.8) then 
c                                 get P-T conditions:
                  write (*,'(/,a)') 'Enter p(bar), T(K): '

                  xo = 1d0
                  igo = 1

                  do 
                     read (*,*,iostat=ier) p, t
                     if (ier.eq.0) exit
                     call rerr
                  end do
        
               else if (ifug.eq.19) then

                  write (*,'(/,a)') 'Enter p(bar), T(K), X(O), X(S):'

                  do 
                     read (*,*,iostat=ier) p, t, xo, elag
                     if (ier.eq.0) exit
                     call rerr
                  end do

               else if (ifug.eq.20) then

                  write (*,'(/,a)') 'Enter p(bar), T(K), X(O), X(C):'
                  do 
                     read (*,*,iostat=ier) p, t, xo, elag
                     if (ier.eq.0) exit
                     call rerr
                  end do

               else if (ifug.eq.24) then 
         
                  if (ibuf.ne.3) then
 
                     write (*,'(/,a)') 'Enter p(bar), T(K), molar N/C:'

                     do 
                        read (*,*,iostat=ier) p, t, gz
                        if (ier.eq.0) exit
                        call rerr
                     end do

                  else 

                     write (*,'(/,a)') 
     *                    'Enter p(bar), T(K), log10[f(O2)], molar N/C:'

                     do 
                        read (*,*,iostat=ier) p, t, dlnfo2, gz
                        if (ier.eq.0) exit
                        call rerr
                     end do

                     dlnfo2 = tentoe*dlnfo2

                  end if 

                  if (igo.eq.0) then 

                     igo = 1
 
                     write (*,'(/,a)') 'Enter log10[a(gph/dia)]:'

                     do 
                        read (*,*,iostat=ier) elag
                        if (ier.eq.0) exit
                        call rerr
                     end do

                  end if 

               else if (ifug.eq.27) then 
         
                     write (*,'(/,a)') 
     *                    'Enter p(bar), T(K), XO, YC:'

                    do 
                        read (*,*,iostat=ier) p, t, xo, fs2
                        if (ier.eq.0) exit
                        call rerr
                     end do

               else 
c                                  or get P-T-X/f conditions:
                  write (*,'(/,a,a,a)') 
     *                   'Enter p(bar), T(K), ',vname(3),':'

                  do 
                     read (*,*,iostat=ier) p, t, xo
                     if (ier.eq.0) exit
                     call rerr
                  end do

                  if (ifug.eq.7.or.ifug.eq.8.or.ifug.eq.24) 
     *               dlnfo2 = tentoe * xo

               end if 
c                                  quit if p = 0
               if (p.eq.0d0) exit
c                                  if sulfur dependent, get
c                                  fs2 if user hasn't opted for
c                                  a buffer:
               if (ibuf.eq.3.and.igo.ne.0.and.(ifug.eq.12.or.ifug.eq.17
     *             .or.ifug.eq.19.or.ifug.eq.20)) then 

                  write (*,'(/,a)') 'Enter log10[f(S2)]:'

                  do 
                     read (*,*,iostat=ier) dlnfo2
                     if (ier.eq.0) exit
                     call rerr
                  end do 

                  dlnfo2 = tentoe * dlnfo2

               end if 
c                                  call fluid routine:
               call cfluid (fo2, fs2)
               write (*,*) ' '
               call rfluid (3)
               write (*,1280) p,t
c                                  output results:
               igo = 1

               if (ifug.ne.27) then 
                  fhc(1) = dexp(fhc(1)) 
                  fhc(2) = dexp(fhc(2)) 
                  fo2 = fo2 / tentoe
               end if 

               if (ifug.le.2.or.ifug.eq.5.or.ifug.eq.14) then

                  write (*,1130) fhc(1), fhc(2)
c                                  finite difference estimate of volume:
                  vdif = 0d0
                  p = p + 0.5d0
                  f = 1d0 
                  do l = 1, 2
                     call cfluid (fo2,fs2)
                     vdif = vdif + f*((1d0-xo)*fhc(1) + xo*fhc(2))
                     f = -1d0
                     p = p - 1d0
                  end do 

                  write (*,1300) 83.14d0*t*vdif

               else if (ifug.eq.13.or.ifug.eq.15) then

                  write (*,1160) fhc(1),fhc(2),fo2
                  write (*,1300) vol

               else if (ifug.eq.25) then 
c                                 primitive output for Leonya's EoS
                  write (*,1130) fhc(1), fhc(2)

               else 

                  if (ifug.eq.16.or.ifug.eq.17) then
                     ag = 0d0
                  else if (ifug.eq.19.or.ifug.eq.20) then
                     ag = dexp (gz)
                  else if (ifug.ne.26) then 
                     ag = dexp (elag)
                  end if 
c                                  routine cfluid returns ln(fs2)/2
                  if (ifug.eq.26) then
                     write (*,1120) fhc(1)/tentoe,fhc(2)/tentoe
                  else if (ifug.eq.27) then 
c                                   ac, fo2, fh2 computed from major
c                                   species, and direct values:
                     if (xs(7).eq.0d0.or.xs(4).eq.0d0.or.
     *                                   xs(3).eq.0d0) then
                        fh2 = dlog(p*xs(5)*g(5))
                        fo2 = 2d0*(dlog(p*xs(1)*g(1)) - eqk(1) - fh2)
                        ag  = p*xs(2)*g(2)*dexp(-fo2 - eqk(2))
                     else if (xs(5).eq.0d0) then 
                        fo2 = dlog(p*xs(7)*g(7))
                        fh2 = dlog(p*xs(1)*g(1)) - eqk(1) - fo2/2d0
                        ag  = p*xs(2)*g(2)*dexp(-fo2 - eqk(2))
                     else if (xs(1).eq.0d0.or.xs(2).eq.0d0) then
                        fh2 = dlog(p*xs(5)*g(5))
                        fo2 = dlog(p*xs(7)*g(7))
                        ag  = p*xs(3)*g(3)*dexp(-fo2/2d0 - eqk(3))
                     end if 

                     write (*,1110) ag,fo2/tentoe,fh2/tentoe

                  else

                     write (*,1170) fo2, 2d0*fs2/tentoe, ag

                  end if 
c                                  compute atomic sums:
                  call elmnts (nc,no,nh,nn,ns,nsi) 
c                                  save old speciation and initialize for fd
                  totx = 0d0

                  do k = 1, isp
                     totx = totx + xs(ins(k))
                     vpar(k) = 0d0
                     xxs(ins(k)) = xs(ins(k))
                     xg(ins(k)) = g(ins(k))
                  end do 
c                                  volumes by finite difference
                  p = p + 0.5d0
                  f = 1d0 

                  do l = 1, 2

                     call cfluid (fo2,fs2)

                     do k = 1, isp

                        if (g(ins(k))*p*xs(ins(k)).eq.0d0) cycle

                        vpar(k) = vpar(k) +  
     *                       83.14d0*t*f*dlog(g(ins(k))*p*xs(ins(k)))

                     end do 

                     f = -1d0
                     p = p - 1d0

                  end do 
c                                  output speciation:
                  write (*,1230)

                  do j = 1, isp, 4
                     kmax = j + 3
                     if (kmax.gt.isp) kmax = isp
                     write (*,1180) (specie(ins(k)), k = j, kmax)
                     write (*,1190) ' x    ',
     *                            (xxs(ins(k)), k = j, kmax)
                     write (*,1190) ' f,bar',
     *                            (xg(ins(k))*p*xxs(ins(k)), k = j,kmax)
                     write (*,1190) ' v,cm3',
     *                             (vpar(k), k = j,kmax)
                     write (*,'(/)')

                  end do 

                  write (*,1370) totx
                  if (dabs(totx-1d0).gt.1d-3) call warn (177,p,i,y)
c                                  output bulk properties and V:
                  write (*,1240)

                  tot = ns + no + nh + nc + nn + nsi    

                  write (*,1250) nc/tot, nh/tot, no/tot, ns/tot, nn/tot,
     *                           nsi/tot
                  if (nh.ne.0d0)  write (*,1270) no/(no+nh)
                  if (ns.ne.0d0)  write (*,1310) ns/(ns+nc)
                  if (nc.ne.0d0)  write (*,1350) nc/(nc+no+nh+ns)
                  if (nsi.ne.0d0) write (*,1390) nsi/(no+nsi)
                  if (nn.ne.0d0)  write (*,1400) nn/nc
                  if (nh.ne.0d0.and.ns.ne.0d0) write (*,1290) ns/nh
                  if (ifug.eq.19.or.ifug.eq.20.and.ag.gt.1d0) 
     *               write (*,1330)

                  write (*,1420) vol

               end if 

            end do 

         end if 

         write (*,'(/,a)') 'More calculations (y/n)?'
         read (*,'(a)') y
         if (y.ne.'y'.and.y.ne.'Y') exit

      end do 
1110  format (/,10x,'a(gph/dia) = ',g12.5,
     *        /,10x,'log[f(O2)] = ',g12.5,
     *        /,10x,'log[f(H2)] = ',g12.5,/)
1120  format (/,10x,'log[f(O) ] = ',g12.5,
     *        /,10x,'log[f(Si)] = ',g12.5,/)
1130  format (/,10x,'f(H2O) = ',g12.5,/,10x,'f(CO2) = ',g12.5,/)
1140  format (/,10x,'f(',a,') = ',g12.5,/)
1150  format (/,10x,'f(H2O)     = ',g12.5,/
     *         ,10x,'f(CO2)     = ',g12.5,/
     *         ,10x,'log[f(O2)] = ',g12.5,/)
1160  format (/,10x,'f(H2O)     = ',g12.5,
     *        /,10x,'f(H2 )     = ',g12.5,
     *        /,10x,'log[f(O2)] = ',g12.5,/)
1170  format (/,10x,'log[f(O2)] = ',g12.5,
     *        /,10x,'log[f(S2)] = ',g12.5,
     *        /,10x,'a(gph/dia) = ',g12.5,/)
1180  format (10x,4(4x,a,5x))
1190  format (2x,a,2x,4(g12.5,1x))
1230  format (/,16x,'Speciation/Fugacities/Partial Molar Volumes',/)
1240  format (/,22x,'Atomic Proportions',//,
     *        10x,'C',12x,'H',12x,'O',12x,'S',12x,'N',12x,'Si')
1250  format (4x,6(g12.5,1x),/)
1270  format (5x,'Back-calculated X(O) = ',g16.9)
1280  format (/,10x,'p(bar)     = ',g12.5,/,10x,'T(K)       = ',g12.5)
1290  format (5x,'S/H = ',g12.5,/)
1300  format (10x,'V(cm3/mol) = ',g12.5,/)
1310  format (5x,'Back-calculated X(S) = ',g16.9)
1330  format (/,5x,'COMPOSITION IS SUPERSATURATED WITH ',
     *           'RESPECT TO CARBON!!',/)
1350  format (5x,'Back-calculated X(C) = ',g16.9)
1370  format (/,5x,'Sum of species fractions: ',f14.9,/)
1380  format (/,5x,'INVALID SPECIATION!!',/)
1390  format (5x,'Back-calculated X(Si) = ',g16.9)
1400  format (5x,'Back-calculated N/C  = ',g16.9)
1420  format (/,5x,'Molar Volume (cm3/mol) = ',g12.5,/)

      end 

      subroutine elmnts (nc,no,nh,nn,ns,nsi) 
c----------------------------------------------------------------------
c compute element totals per mol of fluid species, hardwired stoichiometries
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision nc,no,nh,nn,ns,nsi

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)
c----------------------------------------------------------------------
      ns = y(6) + y(8) + y(9) 
      no = y(1) + y(2)*2d0 + y(3) + y(7)*2d0 + y(12)
     *     + y(8)*2d0 + y(9) + y(14)*2d0 + y(13)
      nc = y(2) + y(3) + y(4) + y(9) + y(16)*2d0
      nh = (y(1) + y(5) + y(6))*2d0 + y(4)*4d0 
     *     + y(11)*3d0 + y(16)*6d0
      nn = 2d0*y(10) + y(11)
      nsi = y(13) + y(14) + y(15)

      end
