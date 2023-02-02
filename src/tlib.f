c----------------------------------------------------------------------

c TLIB - a library of subprograms called by the PERPLEX programs.

c Copyright (C) 1986-2023 James A D Connolly

c This file is part of Perple_X.

c Perple_X is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2, or (at your option)
c any later version.

c Perple_X is distributed in the hope that it will be useful,
c but without any warranty; without even the implied warranty of
c merchantability or fitness for a particular purpose.  See the
c GNU General Public License for more details.

c You should have received a copy of the GNU General Public License
c along with Perple_X (file license.txt). If not see
c <http://www.gnu.org/licenses/>.

c----------------------------------------------------------------------

      subroutine vrsion (n)
c----------------------------------------------------------------------
c a version stamp for each executable
c----------------------------------------------------------------------
      implicit none

      integer n

      write (n,'(/,a,//,a)') 
     *     'Perple_X release 7.0.5, February 2, 2023.',

     *     'Copyright (C) 1986-2023 James A D Connolly '//
     *     '<www.perplex.ethz.ch/copyright.html>.'

      end

      logical function chksol (new)
c----------------------------------------------------------------------
c check that the version flag in the solution model file is consistent.
c This creates headaches, but is used to prevent old versions of Perple_X 
c from crashing while reading a new solution model format
c----------------------------------------------------------------------
      implicit none

      character*3 new

      if (new.eq.'682'.or.new.eq.'683'.or.new.eq.'688'.or.
     *    new.eq.'685'.or.new.eq.'687') then

         call error (3,0d0,0,new)

      else if (new.eq.'008'.or.new.eq.'011'.or.new.eq.'670'.or.
     *         new.eq.'672'.or.new.eq.'673'.or.new.eq.'674'.or.
     *         new.eq.'675'.or.new.eq.'676'.or.new.eq.'678'.or.
     *         new.eq.'679'.or.new.eq.'689'.or.new.eq.'690'.or.
     *         new.eq.'691') then 

         chksol = .true.

      else 

         chksol = .false.

      end if 

      end 

      subroutine redop1 (output,opname)
c----------------------------------------------------------------------
c redop1 - redop1 looks for the perplex_option.dat file, if it finds
c the option file it scans the file for keywords and sets options
c accordingly. also sets machine precision dependent (wmach) and 
c fractional constants. also checks for goofball parameter choices.

c iam - indicates calling program 1 - vertex
c                                 2 - meemum
c                                 3 - werami
c                                 all other values no output

c special internal values for lopt, iopt, nopt

c               lop_28-30 
c               iop_28-30
c               nop_28-30
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, jer, i, loopx, loopy, ibeg, iend, ik1, ik21

      logical output

      character*3 key*22, val, nval1*12, nval2*12,
     *            nval3*12,opname*100,strg*40,strg1*40

      double precision r2

      double precision dnan

      external dnan

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9

      integer iam
      common/ cst4 /iam

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      double precision wmach
      common/ cstmch /wmach(10)

      integer itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *        maxact, mxfree, maxnz
      common/ ngg010 /itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *                maxact, mxfree, maxnz
c----------------------------------------------------------------------
c                                 periodic fractions
      r13 = 1d0/3d0
      r23 = 2d0/3d0
      r43 = 4d0/3d0
      r59 = 5d0/9d0
c                                 seismic speed conversion, this shouldn't be
c                                 hardwired.
      units = dsqrt(1d1)/1d1
c                                 -------------------------------------
c                                 loop to find machine precision (mainly
c                                 for nag)
      r1 = 1d-12
      r2 = 0d0

      do
         if (1d0+r1.eq.1d0) exit
         r2 = r1 
         r1 = r1/2d0
      end do 

      if (r2.eq.0d0) call error (72,r1,i,
     *                         'starting precision for r1 < zero')
c                                 wmach(1-2,5,9) do not have the 
c                                 standard BLAS values, additionally
c                                 BLAS routines a dimension
c                                 of 15 for wmach.
c                                 as currently set wmach(1,2) are

c                                 function precision 
      wmach(1) = r2**0.9d0
c                                 optimality tolerance.
      wmach(2) = r2**0.8d0
c                                 precision (epsmch)
      wmach(3) = r2
c                                 feasibility tolerance (often used as numeric zero)
      wmach(4) = dsqrt(r2)

      epspt3 = r2**0.33d0
      epspt5 = wmach(4)
      epspt8 = wmach(2)
      epspt9 = wmach(1)

      wmach(9) = dmax1(1d0/wmach(4),1d2)

      nopt(51) = wmach(2)
      nopt(52) = 2d0*wmach(3)
      nopt(53) = wmach(1)
c                                 infinite log + 1, for configurational entropy derivatives
      nopt(50) = wmach(1)

      nopt(55) = 1d0 + nopt(50)
      nopt(56) = 1d0 - nopt(50)
c                                 the biggest exponent
      nopt(57) = dlog(huge(0d0)**(0.99d0))
      wmach(5) = 1d0 + r2
c                                 largest number, sqrt
      wmach(7) = huge(0d0)**(0.99d0)
      wmach(8) = dsqrt(wmach(7))
c                                 smallest number, sqrt
      wmach(10) = 1d0/wmach(7)
      if (wmach(10).eq.0d0) 
     *    call errdbg ('1/wmach(7) = 0, reduce wmach(7)')
      wmach(6) = dsqrt(wmach(10))
c                                 solution composition zero and one
      zero = dsqrt(r2)
      r1 = 1d0 + zero
      one = 1d0 - zero
c                                 -------------------------------------
c                                 default option values:
c                                 reserved for temporary use:
      do i = 28, 30
         nopt(i) = 1d0
         lopt(i) = .false.
         iopt(i) = 0
      end do 
c                                 -------------------------------------
c                                 minimum replicate label distance
      nopt(4) = 0.025
c                                 speciation_factor
      nopt(5) = 1d-5
c                                 vrh_weighting keyword
      nopt(6) = 0.5d0
c                                 bad_number keyword
      nopt(7) = dnan()
c                                 zero_mode (<0 off)
      nopt(9) = 1d-6
c                                 set zero threshold for fractionation calculations
      if (icopt.eq.7.or.icopt.eq.9.or.icopt.eq.12) nopt(9) = zero
c                                 tolerance below which a component is considered to 
c                                 be zero during fractionation
      nopt(11) = 1d-6
c                                 quench temperature (K)
      nopt(12) = 0d0
c                                 initial_resolution, solvus_tolerance (output)
      if (iam.eq.15) then 
c                                 convex
         nopt(13) = 1d0/16d0
         nopt(8) = 1.5*nopt(13)
      else
c                                 meemum, vertex
         nopt(13) = 0.2d0
         nopt(8) = 5d-2
      end if
c                                 solvus_tolerance_II (computational)
      nopt(25) = 1d0/16d0
c                                 compositional resolution for conformal
c                                 subdivision
      nopt(14) = 2d-3
c                                 perturbation to eliminate pseudocompound
c                                 degeneracies
      nopt(15) = 5d-3
c                                 poisson ratio to be use for calculating
c                                 missing shear moduli, see iopt(15)
      nopt(16) = 0.35d0
c                                 increase in resolution for adaptive minimization 
      nopt(17) = 3d0
c                                 increase in resolution for composition and mixed variable diagram calculations
      nopt(18) = 1d1
c                                 increase in resolution for Schreinemakers diagram calculations   
      nopt(19) = 3d0 
c                                 T_melt cutoff 
      nopt(20) = 873d0
c                                 optimization_precision, absolute
      nopt(21) = 1d-4
c                                 finite_difference_p threshold for finite difference estimates
      nopt(26) = 1d4
c                                 finite_difference_p fraction for first order difference estimates
      nopt(27) = 1d-3
c                                 fd_expansion_factor is the factor by which finite difference
c                                 increments increase for higher order derivatives
      nopt(31) = 2d0
c                                 fractionation_lower_threshold
      nopt(32) = 0d0
c                                 fractionation_upper_threshold
      nopt(33) = 0d0
c                                 aq_vapor_epsilon
      nopt(34) = 1d0
c                                 replicate_threshold (savdyn)
      nopt(35) = 1d-2
c                                 rep_dynamic_threshold (savrpc)
      nopt(37) = 1d-3
c                                 scatter_increment
      nopt(48) = 1d-2
c                                 MINFRC_diff_increment
      nopt(49) = 1d-7
c                                 -------------------------------------
c                                 composition_phase
      iopt(2) = 0 
      valu(2) = 'mol'
c                                 porportions
      iopt(3) = 0 
      valu(3) = 'vol'
c                                 interpolation keyword
      iopt(4) = 2
      valu(4) = 'on '
c                                 autorefine, 2 - automatic, 1 - manual, 0 - no
      iopt(6) = 2
      valu(6) = 'aut'
c                                 LP_max_it (formerly L6), the max number of iterations
c                                 allowed during an LP optimization, theoretically
c                                 5*(phct+icp) iterations may be required, generally
c                                 the actual number of iterations is < 100. 
      iopt(7) = itmax2
c                                 subdivision_override, 0 - solution model choice
c                                 1 - all cartesian, 2 - all stretch
      iopt(13) = 0 
      valu(13)  = 'off'
c                                 poisson ratio switch fpr calculating
c                                 missing shear moduli
      iopt(16) = 1
      valu(15) = 'on '
c                                 assume linear boundaries within a cell during gridded minimization
      iopt(18) = 1
      valu(18) = 'on '
c                                 seismic data output for WERAMI/MEEMUM, 0 - none, 1 - some, 2 - all
      iopt(14) = 1
      valu(14) = 'som'
c                                 optimization_max_it, max number of iterations
      iopt(20) = 40
c                                 speciation_max_it - for speciation calculations
      iopt(21) = 100
c                                 aq_bad_results 
c                                       0 - err - treat as bad result (optimization error)
c                                       1 - 101 - cease iteration at solute component saturation
c                                       2 - 102
c                                       3 - 103
c                                      99 - ign - ignore 102/103 conditions and cease as in 101.
      iopt(22) = 0
      valu(5) = 'err'
c                                 for infiltration/fractionation calculations set default to 101
c     if (icopt.gt.6) then 
c        iopt(22) = 1
c        valu(5) = '101'
c     end if 
c                                 solution_names 0 - model, 1 - abbreviation, 2 - full
      iopt(24) = 0
      valu(22) = 'mod'
c                                 hyb_h2o - eos to be used for pure h2o, 0-2, 4-5, 6-7
      iopt(25) = 4
c                                 hyb_co2 - eos to be used for pure co2, 0-4, 6-7
      iopt(26) = 4
c                                 hyb_ch4 - eos to be used for pure ch4, 0-1, 6-7
      iopt(27) = 0

c     iopt(28-30)                 reserved as debug options iop_28 - iop_30

c                                 refinement_points_II renamed refinement_points
      iopt(31) = 5
c                                 maximum number of aqueous species
      iopt(32) = 20
c                                 interim_results, 1 - auto, 0 - off, 2 - man
      iopt(34) = 1
      valu(34) = 'aut'
c                                 minfxc_solver
c                                 -1 - use minfxc as default optimizer
c                                  0 - only use minfx when speci2 sets minfx, do not set minfx when on a constraint
c                                  1 - set minfx on any constraint, but allow speci2 to continue
c                                  2 - set minfx on any constraint, only continue for icase = 0
c                                  3 - set minfx on any constraint, continue for all cases.
      iopt(37) = 0
c                                 dynamic_LP_start
c                                  0 - cold
c                                  1 - warm
c                                  2 - hot
      iopt(38) = 1
      valu(38) = 'war'
c                                 static_LP_start
c                                  0 - cold
c                                  1 - warm
c                                  2 - hot
      iopt(39) = 1
      valu(39) = 'war'
c                                 keep_max
      iopt(52) = 20000
c                                 -------------------------------------
c                                 closed or open compositional space
      lopt(1) = .true.
c                                 lopt(2) set by WERAMI
c     lopt(2)
c                                 hard_limits for solution model refinement
      lopt(3) = .false.
      valu(16) = 'off'
c                                 Anderson-Gruneisen correction
      lopt(4) = .false.
c                                 auto_exclude 
      lopt(5) = .true.
c                                 melt_is_fluid
      lopt(6) = .true.
c                                 set locally
c     lopt(7) 
c                                 approx_alpha
      lopt(8) = .true.
c                                 automatic solvus tolerance
      lopt(9) = .true.
c                                 pseudocompound_file
      lopt(10) = .false.
c                                 auto_refine_file
      if (iam.eq.15) then 
         lopt(11) = .true.
      else 
         lopt(11) = .false.
      end if
c                                 option_list_files
      lopt(12) = .false.
c                                 automatic solvus_tolerance_II
      lopt(13) = .true.
c                                 logarithimic P
      lopt(14) = .false.
c                                 spreadsheet format
      lopt(15) = .true.
c                                 averaging scheme
      lopt(16) = .true.
      valu(19) = 'VRH'
c                                 use explicit bulk modulus when available
      lopt(17) = .true.
c                                 refine_bad_nodes -> not used
      lopt(18) = .true. 
c                                 pause_on_error
      lopt(19) = .true.
c                                 poisson_test (reject results with poisson < 0)
      lopt(20) = .false.
c                                 species_output
      lopt(21) = .true. 
c                                 composition_constant
      lopt(22) = .false.
c                                 composition_system (true = wt)
      lopt(23) = .true. 
      valu(21) = 'wt '
c                                 output endmember gibbs energies (werami/meemum)
      lopt(24) = .false. 
c                                 aq_output output back-calculated solute speciation
      lopt(25) = .true.
c                                 aq_solvent_composition (true = molar)
      lopt(26) = .true.
      valu(26) = 'y'
c                                 aq_solute_composition (true = molal)
      lopt(27) = .true.
      valu(27) = 'm'
c                                 lop_28, lop_29, lop_30 temporary options
c     lopt(28-30)
c                                 warning_ver637, werami interpolation
      lopt(31) = .true.
c                                 aq_lagged_speciation
      lopt(32) = .false.
c                                 output_iteration_g
      lopt(33) = .false.
c                                 output_iteration_details
      lopt(34) = .false.
c                                 fractionation threshold, set by nopt(32)/nopt(33)
      lopt(35) = .false.
c                                 aq_oxide_components
      lopt(36) = .false.
c                                 logarithmic X
      lopt(37) = .false.
c                                 non_linear_switch
      lopt(38) = .false.
c                                 refine_endmembers
      lopt(39) = .false.
c                                 automatic specification of refinement_points_II
      lopt(40) = .true.
c                                 absolute (amounts)
      lopt(41) = .false.
c                                 cumulative (amounts)
      lopt(42) = .false.
c                                 reject_negative_sites
      lopt(43) = .true.
c                                 aq_ion_H+ 
      lopt(44) = .true.
c                                 fancy_cumulative_modes
      lopt(45) = .false.
c                                 aq_solvent_solvus
      lopt(46) = .true.
c                                 sample_on_grid 
      lopt(48) = .true.
c                                 refinement_switch
      lopt(49) = .false.
c                                 seismic_data_file
      lopt(50) = .true.
c                                 structural formula options
      lopt(51) = .true.
c                                 keep_auto
      lopt(52) = .true.
c                                 scatter-points
      lopt(54) = .true.
c                                 re-refine
      lopt(55) = .false.
c                                 warning_ver017, relict equipartition -> now warn_interactive
      lopt(56) = .true.
c                                 intermediate_savrpc calls
      lopt(57) = .false.
c                                 intermediate_savdyn calls
      lopt(58) = .false.
c                                 keep_all_rpcs
      lopt(59) = .true.
c                                 warning_ver013, negative composition
      lopt(60) = .true.
c                                 timing
      lopt(61) = .true.
c                                 order_check
      lopt(62) = .false.
c                                 allow GFSM/disable saturated phase
      lopt(63) = .false.
c                                 override counter limits for (some) warnings
      lopt(64) = .false.
c                                 fluid_shear_modulus
      lopt(65) = .true.
c                                 compute_FD_increments for MINFRC
      lopt(66) = .false.
c                                 phi_d
      nopt(65) = 0.36
c                                 initialize mus flag lagged speciation
      mus = .false.
c                                 -------------------------------------
c                                 for gridded minimization:
c                                 # nodes in i direction
      grid(1,1) = 10 
      grid(1,2) = 40
c                                 # nodes in j direction
      grid(2,1) = 10 
      grid(2,2) = 40
c                                 # of levels
      grid(3,1) = 1
      grid(3,2) = 4
c                                 1d fractionation path
      grid(4,1) = 40 
      grid(4,2) = 150
c                                 -------------------------------------
c                                 for schreinemakers etc:
c                                 max variance 
      grid(5,1) = 1
      grid(5,2) = 99
c                                 default increment (relative)
      rid(1,1) = 0.1d0
      rid(1,2) = 0.025d0
c                                 reaction format
      ifull = 0 
      valu(7) = 'min'
c                                 reaction list 
      lopt(53) = .false.
      valu(9) = 'off'
c                                 console msgs
      imsg = 0
      valu(8) = 'on '
c                                 efficiency level
      isec = 3
c                                 short print
      io3p = 1
      valu(10) = 'on '
c                                 -------------------------------------
c                                 look for file
      open (n8, file = opname, iostat = jer, status = 'old')
c                                 if no option file (jer.ne.0) use defaults
      ier = jer
c                                 read cards to end of 
c                                 option file
      do while (ier.eq.0) 

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 ier ne 0 = eof
         if (ier.ne.0) exit 
c                                 cycle on default specification
         if (strg.eq.'default') cycle
c                                 if here we should have a keyword and
c                                 value
         if (key.eq.'composition'.or.key.eq.'composition_phase') then 
c                                 phase composition key
            if (val.eq.'wt') then
               iopt(2) = 1
               valu(2) = val 
            end if 

         else if (key.eq.'aq_species') then

            read (strg,*) iopt(32)

         else if (key.eq.'auto_exclude') then 

            if (val.ne.'T') lopt(5) = .false.

         else if (key.eq.'aq_output') then

            if (val.ne.'T') lopt(25) = .false.

         else if (key.eq.'aq_lagged_speciation') then 

            if (val.eq.'T') lopt(32) = .true.

         else if (key.eq.'aq_oxide_components') then 

            if (val.eq.'T') lopt(36) = .true.

         else if (key.eq.'aq_ion_H+') then 

             if (val.eq.'F') lopt(44) = .false.

         else if (key.eq.'fancy_cumulative_modes') then 

             if (val.eq.'T') lopt(45) = .true.

         else if (key.eq.'non_linear_switch') then 

            if (val.eq.'T') lopt(38) = .true.

         else if (key.eq.'output_iteration_G') then 

            if (val.eq.'T') lopt(33) = .true.

         else if (key.eq.'output_iteration_detai') then 

            if (val.eq.'T') lopt(34) = .true.

         else if (key.eq.'aq_bad_results') then 

            if (val.eq.'err') then 
c                                 abort on any hint of trouble
               iopt(22) = 0
            else if (val.eq.'101') then 
c                                 continue on solute undersaturation (unwise)
               iopt(22) = 1
            else if (val.eq.'102') then 
c                                 continue if pure solvent coexists with immiscible impure solvent
               iopt(22) = 2
            else if (val.eq.'103') then
c                                 abort if pure solvent is stable
               iopt(22) = 3
            else if (val.eq.'ign') then 
               iopt(22) = 99
            end if

            valu(5) = val

         else if (key.eq.'refine_endmembers') then 

            if (val.eq.'T') lopt(39) = .true.

         else if (key.eq.'absolute') then 

            if (val.eq.'T') lopt(41) = .true.

         else if (key.eq.'cumulative') then 

            if (val.eq.'T') lopt(42) = .true.

         else if (key.eq.'reject_negative_sites') then 

            if (val.eq.'F') lopt(43) = .false.

         else if (key.eq.'aq_solvent_composition') then

            if (val.ne.'y') then
               lopt(26) = .false.
               valu(26) = 'm'
            end if 

         else if (key.eq.'aq_solute_composition') then

            if (val.ne.'m') then
               lopt(27) = .false.
               valu(27) = 'y'
            end if 

         else if (key.eq.'hybrid_EoS_H2O') then

            read (strg,*) iopt(25)

            if (iopt(25).lt.0.or.iopt(25).gt.7.or.iopt(25).eq.3) then 
               write (*,1180) strg,key
               call errpau
            end if 

         else if (key.eq.'hybrid_EoS_CO2') then

            read (strg,*) iopt(26)

            if (iopt(26).lt.0.or.(iopt(26).gt.4.and.iopt(26).ne.7)) then 
               write (*,1180) strg,key
               call errpau
            end if

         else if (key.eq.'hybrid_EoS_CH4') then

            read (strg,*) iopt(27)

            if (iopt(27).lt.0.or.(iopt(27).gt.1.and.iopt(27).ne.7)) then 
               write (*,1180) strg,key
               call errpau
            end if 

         else if (key.eq.'proportions') then 
c                                 phase proportion key
c                                 volume is default
            if (val.eq.'wt') then
               iopt(3) = 1
            else if (val.eq.'mol') then 
               iopt(3) = 2 
            end if
 
            valu(3) = val

         else if (key.eq.'interpolation') then 
c                                 interpolation key
            if (val.ne.'on ') then
               iopt(4) = 0
               valu(4) = val
            else 
               read (nval1,*) iopt(4)
            end if 

         else if (key.eq.'bounds') then 
c                                 
            if (val.eq.'HS'.or.val.eq.'hs')  then
              lopt(16) = .false.
              valu(19) = 'HS '
            end if 

         else if (key.eq.'vrh/hs_weighting'.or.
     *            key.eq.'vrh_weighting') then 
c                                 vrh/hs weighting key
            read (strg,*) nopt(6)

         else if (key.eq.'bad_number') then
c                                 bad number key 
            if (val.eq.'NaN'.or.val.eq.'nan') then
               nopt(7) = dnan()
            else 
               read (strg,*) nopt(7)
            end if 

         else if (key.eq.'solution_names') then 

              if (val.eq.'abb') then 
                 iopt(24) = 1
              else if (val.eq.'ful') then 
                 iopt(24) = 2
              end if 

              valu(22) = val

         else if (key.eq.'solvus_tolerance') then 

            if (val.ne.'aut') then 
               lopt(9) = .false.
               read (strg,*) nopt(8)
            end if 

         else if (key.eq.'solvus_tolerance_II') then 

            if (val.ne.'aut') then 
               lopt(13) = .false.
               read (strg,*) nopt(25)
            end if 

         else if (key.eq.'speciation_factor') then 

         else if (key.eq.'replicate_threshold') then 

            read (strg,*) nopt(35)

         else if (key.eq.'rep_dynamic_threshold') then 

            read (strg,*) nopt(37)

         else if (key.eq.'MINFRC_diff_increment') then 

            read (strg,*) nopt(49)

         else if (key.eq.'speciation_precision') then 

            read (strg,*) nopt(5)

         else if (key.eq.'optimization_precision') then 

            read (strg,*) nopt(21)

         else if (key.eq.'LP_max_it') then 

            read (strg,*) iopt(7)
            itmax2 = iopt(7)

         else if (key.eq.'optimization_max_it') then 

            read (strg,*) iopt(20)

         else if (key.eq.'zero_bulk') then
c                                 zero_bulk key
            read (strg,*) nopt(11)

         else if (key.eq.'aq_vapor_epsilon') then
c                                 "vapor" threshold
            read (strg,*) nopt(34)

         else if (key.eq.'aq_max_molality') then

c             obsolete
         else if (key.eq.'Tisza_test') then

c             to be implemented

         else if (key.eq.'warning_ver637') then
c                                  allow for solvent immiscisibiliy
            if (val.eq.'F') lopt(31) = .false.

         else if (key.eq.'aq_solvent_solvus') then
c                                  allow for solvent immiscisibiliy
            if (val.eq.'F') lopt(46) = .false.

         else if (key.eq.'interim_results') then
c                                  output interim results (VERTEX/PSSECT/WERAMI)
            if (val.eq.'off') then 
               iopt(34) = 0
            else if (val.eq.'man') then 
               iopt(34) = 2
            end if

            valu(34) = val

         else if (key.eq.'minfxc_solver'.or.key.eq.'MINFXC_solver') then
c                                 mnfxc_solver
c                                 -1 - use minfxc as default optimizer
c                                  0 - only use minfx when speci2 sets minfx, do not set minfx when on a constraint
c                                  1 - set minfx on any constraint, but allow speci2 to continue
c                                  2 - set minfx on any constraint, only continue for icase = 0
c                                  3 - set minfx on any constraint, continue for all cases.
c                                  4 - 3 + but don't use minfxc
            read (strg,*) iopt(37)

         else if (key.eq.'scatter-points') then
c                                 generate points scattered about refinement
c                                 point compositions
            if (val.eq.'F') lopt(54) = .false.

         else if (key.eq.'scatter-increment') then
c                                 scatter point increment
            read (strg,*) nopt(48)

         else if (key.eq.'re-refine') then
c                                 allow re-refinement in VERTEX
            if (val.eq.'T') lopt(55) = .true.

         else if (key.eq.'warn_interactive') then
c                                 override interactive warnings with the bad choice.
            if (val.eq.'F') lopt(56) = .false.

         else if (key.eq.'warn_no_limit') then
c                                 override counter limits for (some) warnings
            if (val.eq.'T') lopt(64) = .true.

         else if (key.eq.'fluid_shear_modulus') then
c                                 compute shear modulus assuming textural eq
            if (val.eq.'F') lopt(65) = .false.
 
         else if (key.eq.'compute_FD_increments') then
c                                 compute finite difference increments for MINFRC
            if (val.eq.'T') lopt(66) = .true.

         else if (key.eq.'phi_d') then
c                                 disaggregation porosity for fluid_shear_modulus
            read (strg,*) nopt(65)

         else if (key.eq.'dynamic_LP_start') then
c                                 use cold starts for dynamic LP
            if (val.eq.'col') then 
               iopt(38) = 0
            else if (val.eq.'hot') then 
               iopt(38) = 2
            end if

            valu(38) = val

         else if (key.eq.'static_LP_start') then
c                                 use cold starts for dynamic LP
            if (val.eq.'col') then 
               iopt(39) = 0
            else if (val.eq.'hot') then 
               iopt(39) = 2
            end if

            valu(39) = val


         else if (key.eq.'timing') then
c                                 timing for VERTEX
            if (val.eq.'F') lopt(61) = .false.

         else if (key.eq.'order_check') then
c                                 check multi-species o/d by QP
            if (val.eq.'T') lopt(62) = .true.

         else if (key.eq.'GFSM') then
c                                 allow GFSM solution models/disable saturated phase
            if (val.eq.'T') lopt(63) = .true.

         else if (key.eq.'intermediate_savrpc') then
c                                 use all g's generated by minfrc
            if (val.eq.'T') lopt(57) = .true.

         else if (key.eq.'intermediate_savdyn') then
c                                 save all minfrc solutions for auto-refine stage
            if (val.eq.'T') lopt(58) = .true.

         else if (key.eq.'keep_all_rpcs') then
c                                 accumulate rpcs during iteration
            if (val.eq.'F') lopt(59) = .false.

         else if (key.eq.'sample_on_grid') then
c                                 sample on computational grid (WERAMI)
            if (val.eq.'F') lopt(48) = .false.

         else if (key.eq.'zero_mode') then
c                                 zero_mode key
            read (strg,*) nopt(9)

         else if (key.eq.'iteration'.or.
     *            key.eq.'resolution_factor') then
c                                 how fast resolution improves with iteration
            write (*,'(a)') 'iteration and resolution_factor are obso'//
     *                      'lete options, 6.9.1+'

         else if (key.eq.'initial_resolution') then
c                                 initial_resolution key
            read (strg1,'(40a)') chars(1:40)
            ibeg = 1
            call readfr (nopt(13),ibeg,iend,40,ier)
            if (ier.ne.0) call error (72,nopt(1),iopt(1),key//
     *                                 'has an invalid value.')
c                                 initial resolution
            if (nopt(13).ge.1d0.or.nopt(13).lt.0d0) then
               nopt(13) = 1d0/5d0
               if (iam.eq.15) nopt(13) = 1d0/16d0
               write (*,1050)
            end if

            if (iam.eq.15) then 
c                                  686+ read second value as arf value
               ibeg = iend + 1

               call readfr (nopt(17),ibeg,iend,40,ier)

               if (ier.ne.0) then
c                                  special backward compatibility msg
                  write (*,1010) opname
                  call errpau

               end if 

               if (nopt(17).ge.1d0.or.nopt(17).lt.0d0) then
                  nopt(17) = nopt(13)/3d0
                  write (*,1050)
               end if

               nopt(17) = nopt(13)/nopt(17)
               nopt(18) = nopt(17)
               nopt(19) = nopt(17)

             end if

         else if (key.eq.'final_resolution') then
c                                 final_resolution keys 
            write (*,'(a)') 'final_resolution is an obso'//
     *                      'lete option, 6.9.1+'

         else if (key.eq.'fd_expansion_factor') then 

            read (strg,*) nopt(31)

         else if (key.eq.'finite_difference_p') then 
c                                 p threshold
            read (strg,*) nopt(26)
c                                 p fraction
            read (nval1,*) nopt(27)

         else if (key.eq.'global_reach_increment') then
          
            write (*,'(a)') 'global_reach_increment is an obso'//
     *                      'lete option, 6.9.1+'

         else if (key.eq.'seismic_output') then 
c                                 seismic data output WERAMI/MEEMUM/FRENDLY
            valu(14) = val

            if (val.eq.'non') then 
               iopt(14) = 0
            else if (val.eq.'all') then
               iopt(14) = 2
            else
               valu(14) = 'som'
            end if

         else if (key.eq.'refinement_points_II'.or.
     *            key.eq.'refinement_points') then
c                                 refinement points
            if (val.ne.'aut') then 
               read (strg,*) iopt(31)
               lopt(40) = .false.
            end if

         else if (key.eq.'refinement_switch') then
c                                 during iteration allow metastable 
c                                 refinement points for stable solutions.
            if (val.ne.'F') lopt(49) = .true.

         else if (key.eq.'seismic_data_file') then
 
            if (val.ne.'T') lopt(50) = .false.

         else if (key.eq.'structural_formulae') then

            if (val.ne.'T') lopt(51) = .false. 

         else if (key.eq.'keep_auto') then

            if (val.ne.'T') lopt(52) = .false. 

         else if (key.eq.'keep_max') then 

            read (strg,*) iopt(52)

         else if (key.eq.'max_aq_species_out') then 
c                                 max number of aq species output for
c                                 back-calculated and lagged speciation
            read (strg,*) iopt(32)

         else if (key.eq.'reach_increment_switch') then 
c                                 reach_increment_switch
            write (*,'(a)') 'reach_increment_switch is an obso'//
     *                      'lete option, 6.9.1+'

         else if (key.eq.'stretch_factor') then
c                                 stretch_factor key = b - 1       
            read (strg,*) nopt(14)

         else if (key.eq.'subdivision_override') then 
c                                 subdivision override key
            valu(13) = val

            if (val.eq.'lin') then
               iopt(13) = 1
            else if (val.eq.'str') then
               iopt(13) = 2
            end if 

         else if (key.eq.'auto_refine') then
c                                 autorefine
            valu(6) = val

            if (val.eq.'off') then
               iopt(6) = 0
            else if (val.eq.'man') then
               iopt(6) = 1
            else if (val.eq.'aut') then
               iopt(6) = 2
            end if 

         else if (key.eq.'auto_refine_factor_I') then
c   
            read (strg,*) nopt(17)

         else if (key.eq.'auto_refine_factor_II') then
c   
            read (strg,*) nopt(18)

         else if (key.eq.'auto_refine_factor_III') then
c   
            read (strg,*) nopt(19)

         else if (key.eq.'speciation_max_it') then
c   
            read (strg,*) iopt(21)

         else if (key.eq.'x_nodes') then
c                                 number of x nodes at level 1 before autorefine
            read (strg,*) grid(1,1)
c                                 number of x nodes for autorefine
            read (nval1,*,iostat=ier) grid(1,2)

            if (ier.ne.0) call error (72,r1,i,'the x_nodes option '//
     *                               'requires two values')

         else if (key.eq.'y_nodes') then
c                                 number of y nodes at level 1
            read (strg,*) grid(2,1)
c                                 number of y nodes for autorefine
            read (nval1,*,iostat=ier) grid(2,2)

            if (ier.ne.0) call error (72,r1,i,'the y_nodes option '//
     *                               'requires two values')

         else if (key.eq.'grid_levels') then 
c                                 number of grid levels before autorefine
            read (strg,*) grid(3,1)
c                                 number of grid levels for autorefine
            read (nval1,*,iostat=ier) grid(3,2)

            if (ier.ne.0) call error (72,r1,i,'the grid_levels option'//
     *                               ' requires two values')

         else if (key.eq.'1d_path') then 
c                                 number of grid points for 1d path before autorefine
            read (strg,*) grid(4,1)
c                                 number of grid points for 1d path for autorefine
            read (nval1,*,iostat=ier) grid(4,2)

            if (ier.ne.0) call error (72,r1,i,'the 1d_path option '//
     *                               'requires two values')

         else if (key.eq.'variance') then 
c                                 max variance of traced equilibria before autorefine
            read (strg,*) grid(5,1)
c                                 max variance of traced equilibria for autorefine
            read (nval1,*,iostat=ier) grid(5,2)

            if (ier.ne.0) call error (72,r1,i,'the variance option '//
     *                               'requires two values')

         else if (key.eq.'increment') then 
c                                 default exploratory relative increment    
            read (strg,*) rid(1,1)  
c                                 default autorefine relative increment
            read (nval1,*,iostat=ier) rid(1,2)

            if (ier.ne.0) call error (72,r1,i,'the increment option '//
     *                               'requires two values')

         else if (key.eq.'reaction_format') then 

            valu(7) = val

            if (val.eq.'ful') then 
               ifull = 1
            else if (val.eq.'sto') then 
               ifull = 2
            else if (val.eq.'S+V') then 
               ifull = 3
            else if (val.eq.'eve') then
               ifull = 4
            end if 
  
         else if (key.eq.'console_messages') then 
            
            if (val.eq.'off') then 
               valu(8) = val
               imsg = 1
            end if 

         else if (key.eq.'reaction_list') then

            if (val.eq.'on ') then 
               valu(9) = val
               lopt(53) = .true.
            end if

         else if (key.eq.'efficiency') then 

            read (strg,*) isec

            if (isec.lt.1.or.isec.gt.5) isec = 3 

         else if (key.eq.'short_print') then 

            if (val.eq.'off') then 
               io3p = 0
               valu(10) = val
            end if 

         else if (key.eq.'dependent_potentials') then 

            write (*,'(a)') 'dependent potentials is an obso'//
     *                      'lete option, 6.9.1+'

         else if (key.eq.'hard_limits') then 

            if (val.eq.'on ') then 
               lopt(3) = .true.
               valu(16) = val
            end if

         else if (key.eq.'T_stop') then 
c                                 equilibrium cutoff T (K)    
            read (strg,*) nopt(12)

         else if (key.eq.'T_melt') then 
c                                 cutoff T (K) for melt endmember stability    
            read (strg,*) nopt(20)

         else if (key.eq.'fractionation_hi_limit') then 
c                                 upper fractionation threshold
            read (strg,*) nopt(33)

         else if (key.eq.'fractionation_lo_limit') then 
c                                 lower fractionation threshold
            read (strg,*) nopt(32)

         else if (key.eq.'order_check') then 
c                                 compare local and max disorder state for o-d models
c                                 obsolete.
         else if (key.eq.'linear_model') then   
c                                 assume linear boundaries within a cell during gridded minimization
            if (val.eq.'off') then 
               iopt(18) = 0
               valu(18) = 'off'
            end if 

         else if (key.eq.'closed_c_space') then
 
            if (val.ne.'T') lopt(1) = .false. 

         else if (key.eq.'pause_on_error') then
 
            if (val.ne.'T') lopt(19) = .false. 

         else if (key.eq.'poisson_test') then
 
            if (val.ne.'F') lopt(20) = .true. 

         else if (key.eq.'species_output'.or.
     *            key.eq.'output_species') then
 
            if (val.ne.'T') lopt(21) = .false. 

         else if (key.eq.'composition_constant') then
 
            if (val.ne.'F') lopt(22) = .true. 

         else if (key.eq.'composition_system') then
 
            if (val.ne.'wt') then
               lopt(23) = .false. 
               valu(21) = val
            end if

         else if (key.eq.'species_Gibbs_energies'.or.
     *            key.eq.'output_species_props') then
 
            if (val.ne.'F') lopt(24) = .true.

         else if (key.eq.'logarithmic_p') then
 
            if (val.eq.'T') lopt(14) = .true.

         else if (key.eq.'logarithmic_X') then
 
            if (val.eq.'T') lopt(37) = .true.

         else if (key.eq.'spreadsheet') then
 
            if (val.eq.'F') lopt(15) = .false.

         else if (key.eq.'Anderson-Gruneisen') then

            if (val.eq.'T') lopt(4) = .true.

         else if (key.eq.'approx_alpha') then

            if (val.eq.'F') lopt(8) = .false.

         else if (key.eq.'melt_is_fluid') then 

            if (val.eq.'F') lopt(6) = .false.

         else if (key.eq.'pc_perturbation') then
c                                 perturbation to eliminate pseudocompound degeneracies  
            read (strg,*) nopt(15)

         else if (key.eq.'pseudocompound_file') then

            if (val.eq.'T') lopt(10) = .true.

         else if (key.eq.'auto_refine_file') then

            if (val.eq.'F') then 
               lopt(11) = .false.
            else if (val.eq.'T') then 
               lopt(11) = .true.
            end if

         else if (key.eq.'option_list_files') then

            if (val.eq.'T') lopt(12) = .true.

         else if (key.eq.'explicit_bulk_modulus') then 

            if (val.eq.'F') lopt(17) = .false.

         else if (key.eq.'poisson_ratio') then 
c                                 handle missing shear moduli
            if (val.eq.'on ') then
               read (nval1,*) nopt(16)
            else if (val.eq.'off') then 
               valu(15) = val
               iopt(16) = 0
            else if (val.eq.'all') then 
               read (nval1,*) nopt(16)
               valu(15) = val
               iopt(16) = 2
            end if   
          
         else if (key.eq.'lop_28') then
c                                 reserved values for debugging, etc
            if (val.eq.'T') lopt(28) = .true.
         else if (key.eq.'lop_29') then
            if (val.eq.'T') lopt(29) = .true.
         else if (key.eq.'lop_30') then
            if (val.eq.'T') lopt(30) = .true.
         else if (key.eq.'iop_28') then
            read (strg,*) iopt(28) 
         else if (key.eq.'iop_29') then
            read (strg,*) iopt(29) 
         else if (key.eq.'iop_30') then
            read (strg,*) iopt(30) 
         else if (key.eq.'nop_28') then
            read (strg,*) nopt(28) 
         else if (key.eq.'nop_29') then
            read (strg,*) nopt(29) 
         else if (key.eq.'nop_30') then
            read (strg,*) nopt(30) 
         else if (key.ne.'|') then

            call error (72,nopt(1),iopt(1),key//' is not a valid Perpl'
     *                 //'e_X option file keyword and must be deleted '
     *                 //'or corrected.')

         end if

      end do

      close (n8)
c                                 -------------------------------------
c                                 set permanent parameters for lpsol
c                                 common blocks ngg010, ng011, ngg005
      call lpset 
c                                 -------------------------------------
c                                 set permanent parameters for nlpsol
c                                 common blocks ngg017, ng019.
      call nlpset
c                                 -------------------------------------
c                                 computation dependent options
c                                 -------------------------------------
c                                 automatic specification of metastable
c                                 refinement points
      if (lopt(40)) iopt(31) = icp + 2
c                                 write and optional file choices
      if (iam.ne.14) then 
         if (jer.ne.0) then 
            write (*,1120) opname
         else 
            write (*,1130) opname
         end if
      end if 

      if (iam.eq.1.or.iam.eq.15) then 
c                                 vertex only files:
         if (icopt.eq.1.or.icopt.eq.3) then 

            if (lopt(53)) then 
               call mertxt (tfname,prject,'_reaction_list.txt',0)
               open (n6,file=tfname)
            else 
               tfname = 'not requested'
            end if 

            write (*,1170) tfname

         end if 
c                                 auto refine summary
         if (iopt(6).ne.0.and.io9.eq.0) then
 
            if (lopt(11)) then 
               call mertxt (tfname,prject,'_auto_refine.txt',0)
            else 
               tfname = 'not requested'
            end if 

            write (*,1150) tfname

         end if

      end if

      if (iam.lt.3) then

         if (lopt(50)) then 
            call mertxt (tfname,prject,'_seismic_data.txt',0)
         else 
            tfname = 'not requested'
         end if 

         write (*,'(a)') 'Writing seismic data options to: '//tfname

      end if 
c                                 pseudocompound glossary
      if (io9.eq.0.and.(iam.lt.3.or.iam.eq.15)) then
 
         if (lopt(10)) then 
            call mertxt (tfname,prject,'_pseudocompound_glossary.txt',0)
         else 
            tfname = 'not requested'
         end if 

         write (*,1140) tfname

      end if 
c                                 computational options this is redundant
      if (iam.lt.4.or.iam.eq.15) then 
         if (lopt(12)) then 
            if (iam.eq.1) then           
               call mertxt (tfname,prject,'_VERTEX_options.txt',0)
            else if (iam.eq.2) then 
               call mertxt (tfname,prject,'_MEEMUM_options.txt',0)
            else if (iam.eq.3) then 
               call mertxt (tfname,prject,'_WERAMI_options.txt',0)
            else if (iam.eq.15) then 
               call mertxt (tfname,prject,'_CONVEX_options.txt',0)
            end if
         else 
            tfname = 'not requested'
         end if 
 
         write (*,1160) tfname

      end if 
c                                 -------------------------------------
c                                 dependent parameters and error traps:
c                                 fractionation theshold flag
      if (nopt(33).gt.nopt(32)) lopt(35) = .true. 

      if (iopt(31).gt.k5+2.or.iopt(31).lt.1) then 
         write (*,1090) icp + 2
         iopt(31) = icp + 2
      end if 
c                                 stretching parameters
      if (nopt(14).le.0d0) then 
         write (*,1060)
         nopt(14) = 2d-3
      end if 
c                                 grid parameters
      do i = 1, 2

         if (grid(3,i).le.0.or.grid(3,i).gt.l8) grid(3,i) = 4

         loopy = (grid(2,i)-1) * 2**(grid(3,i)-1) + 1
         loopx = (grid(1,i)-1) * 2**(grid(3,i)-1) + 1

         if (loopy.gt.l7) then 
            call warn (92,nopt(1),loopy,'y_node')
            grid(2,i) = (l7 - 1)/2**(grid(3,i)-1) + 1
         end if 

         if (loopx.gt.l7) then 
            call warn (92,nopt(1),loopx,'x_node')
            grid(1,i) = (l7 - 1)/2**(grid(3,i)-1) + 1
         end if 

         if (grid(4,i).gt.l7) then 
            call warn (92,nopt(1),loopx,'1dpath')
            grid(4,i) = l7 - 1
         end if  

         if (grid(5,i).lt.1) then 
            call warn (113,rid(1,i),grid(5,i),'VARIAN')
            grid(5,i) = 1
         end if 

         if (rid(1,i).lt.1d-2) then 
            call warn (114,rid(1,i),i,'INPUT1')
         end if 

      end do
c                                 --------------------------------------
c                                 program/computation specific settings
c                                 set autorefine factor
      if (icopt.eq.1) then
         nopt(17) = nopt(19)
      else if (icopt.le.3) then 
         nopt(17) = nopt(18)
      end if
c                                 effective initial resolution
      if (nopt(13).eq.0d0) then 
c                                 user wants to use solution model values, set
c                                 initial resolution to a representative value
          rid(3,1) = 0.1d0
      else 
          rid(3,1) = nopt(13)
      end if 

      rid(3,2) = rid(3,1)/nopt(17)
c                                 --------------------------------------
c                                 output
      if (((iam.eq.1.or.iam.eq.15).and.output).or.iam.eq.2
     *                                        .or.iam.eq.3
     *                                        .or.iam.eq.5) then
c                                 console
         call outopt (6)

         if (lopt(12).and.iam.ne.5) then 
c                                 file version, create the file name 
            if (iam.eq.1) then           
               call mertxt (tfname,prject,'_VERTEX_options.txt',0)
            else if (iam.eq.2) then 
               call mertxt (tfname,prject,'_MEEMUM_options.txt',0)
            else if (iam.eq.3) then 
               call mertxt (tfname,prject,'_WERAMI_options.txt',0)
            else if (iam.eq.15) then 
               call mertxt (tfname,prject,'_CONVEX_options.txt',0)
            end if 

            open (n8, file = tfname)
            call outopt (n8)
            close (n8) 

            write (*,1000) tfname

         end if 

      end if 
c                                 -------------------------------------
c                                 recalculate parameters:
c                                 proportionality constant for shear modulus
      nopt(1) = nopt(16)
      nopt(16) = 1.5d0*(1d0-2d0*nopt(16))/(1d0+nopt(16))
c                                 -------------------------------------
c                                 check dynamic memory allocation:
      if (iam.lt.3) then

c        write (*,*) k21, k19
c                                 recommended value for k21, to be absolutely 
c                                 safe multiply by 2.
         ik21 = iopt(52)*icp

         if (k21.le.ik21) then 
c                                 consequent value for k1
            ik1 = (memory - (1 + k32)*ik21)/(2*k31+k32+2)

            write (*,'(a,i8,2a,i7,a/,a,i7,a,/)') 
     *            '**warning ver072** parameter k21 (',k21,') is below',
     *            ' its recommended value (',ik21,')','decrease k1 to ',
     *            ik1,' or reduce keep_max to avoid potential problems'

         end if 

         if (k21.le.100) then

            ik1 = (memory - (1 + k32)*50000)/(2*k31+k32+2)

            write (*,'(a,i7,2a/,a,i7,a,/)') 
     *            '**error ver072** k21 (',k21,') is too far below the',
     *            ' minimum safe value (~10000)','decrease k1 to < ',
     *            ik1,' and recompile Perple_X'

            call errpau

         end if

      end if

1000  format ('Context specific options are echoed in: ',a)
1010  format (/,'ERROR: reading option file: ',a50,/,
     *          'most likely you are using a pre-6.8.6 option',
     *          ' file with only one',/,'initial_resolution value. Corr'
     *         ,'ect the error by adding a second value,',/,'typically',
     *          ' 1/3 the first (exploratory stage) value.')
1050  format (/,'Warning: initial_resolution values must be ',
     *         '< 1',/,'the keyword will be',
     *         ' assigned its default values.',/)
1060  format (/,'Warning: stretch_factor cannot be < 0 ',
     *          'the keyword will be assigned its default value',/)
1070  format (/,'Warning: initial_resolution values must be ',
     *         '< 1',/,'the keyword will be',
     *         ' assigned its default value (',i2,').',/)
1090  format (/,'Warning: refinement_points must be ',
     *         ' >1 and <k5+2',/,'refinement_points will be',
     *         ' assigned its default value [',i2,'].',/)
1120  format (/,'Warning: the Perple_X option file: ',a,/,
     *       'was not found, default option values will be used.',/) 
1130  format (/,'Reading Perple_X options from: ',a)
1140  format ('Writing pseudocompound glossary to: ',a)
1150  format ('Writing auto refine summary to: ',a)
1160  format ('Writing Perple_X option summary to: ',a)
1170  format ('Writing complete reaction list to: ',a)
1180  format (/,'Error: value ',a,' is invalid for Perple_X option ',
     *        'keyword ',a,/,'see www.perplex.ch/perplex_options.html ',
     *        'for a list of valid values',/)
      end 

      subroutine lpset
c----------------------------------------------------------------------
c                                 set permanent parameters for lpsol
c                                 common blocks ngg010, ng011, ngg005,
c                                 itmax2 may be reset in 
c---------------------------------------------------------------------
      implicit none

      integer itnfix, kdegen, ndegen, nfix
      double precision tolinc, tolx0
      common/ ngg005 /tolx0, tolinc, kdegen, ndegen, itnfix, nfix(2)

      integer itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *        maxact, mxfree, maxnz
      common/ ngg010 /itmax1, itmax2, kchk, kcycle, lcrash, lprob, 
     *                maxact, mxfree, maxnz

      double precision bigbnd, bigdx, bndlow, bndupp, tolact, tolfea, 
     *                 tolrnk
      common/ ngg011 /bigbnd, bigdx, bndlow, bndupp,
     *                tolact, tolfea, tolrnk
c----------------------------------------------------------------------
      itmax2 = 500
      kchk = 50
      kcycle = 10000
c                                 kdegen:  expand frequency
      kdegen = kcycle
      tolact = 1d-2
      bigbnd = 0.99999d20
      bigdx = bigbnd
c                                 tolinc: scaled increment to the current featol
      tolinc = 0.49d0/kcycle
c                                 tolx0: the minimum (scaled) feasibility tolerance
      tolx0 = 0.5d0

      end

      subroutine nlpset
c----------------------------------------------------------------------
c                                 set permanent parameters for nlpsol
c                                 common blocks ngg017, ngg019, ngg021
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision fac 

      double precision wmach
      common/ cstmch /wmach(10)

      double precision cdint, ctol, dxlim, epsrf, eta, fdint, ftol,
     *                 hcndbd
      common/ ngg021 /cdint, ctol, dxlim, epsrf, eta,
     *                fdint, ftol, hcndbd

      double precision epspt3, epspt5, epspt8, epspt9
      common/ ngg006 /epspt3, epspt5, epspt8, epspt9
c----------------------------------------------------------------------
      rhomax = 1d0/wmach(3)
      tolact = 1d-2
      bigbnd = 0.99999d20
      bigdx  = bigbnd
c                                 fac = 1d-2 was being used in 691:
      fac = 1d0
c                                 EPSRF, function precision
      epsrf = (wmach(3)*fac)**(0.9)
c                                 FTOL, optimality tolerance
      ftol = (wmach(3)*fac)**(0.8)
c                                 CTOL, feasibility tolerance
      ctol = epspt5
c                                 DXLIM, step limit < nopt(5) leads to bad results
c                                 relates to alfamax as (1+norm(x))*dxlim 
      dxlim = 0.05d0
c                                 ETA, linesearch tolerance, low values -> more accurate search 
c                                 -> more function calls, 0.05-.4 seem best
      eta = 0.225d0
c                                 FDINT, finite difference interval, forward.
      fdint = nopt(49)
c                                 CDINT, 2nd order forward finite difference interval
      cdint = fdint**(0.67d0)
c                                 objective function calls for minfrc
      count = 0
      rcount = 0
      lcount = 0

      end 

      subroutine outopt (n)
c----------------------------------------------------------------------
c outopt - for program IAM, outputs context specific options to LUN N,
c called by redop1  
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, siz, n

      character nval1*12, text(14)*1, numb*5, nval2*12

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      logical oned
      common/ cst82 /oned

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 version
      if (n.ne.6) call vrsion (n)
c                                 generic blurb
      if (iam.eq.1) then 
         write (n,1000) 'VERTEX'
      else if (iam.eq.2) then 
         write (n,1000) 'MEEMUM'
      else if (iam.eq.3) then 
         write (n,1000) 'WERAMI'
      else if (iam.eq.5) then 
         write (n,1000) 'FRENDLY'
      else if (iam.eq.15) then 
         write (n,1000) 'CONVEX'
      end if 

      if (iam.le.2.or.iam.eq.15) then 
c                                 VERTEX/MEEMUM:
c                                 solvus tolerance text
         if (lopt(9)) then 

           nval1 = 'aut    '

         else 

           call numtxt (nopt(8),text,siz)
           write (nval1,'(14a)') (text(i),i=1,siz)

         end if

         if (lopt(13)) then 

           nval2 = 'aut    '

         else 

           call numtxt (nopt(25),text,siz)
           write (nval2,'(14a)') (text(i),i=1,siz)

         end if 

         if (iam.eq.1.or.iam.eq.15) write (n,1015) valu(6), nopt(35),
     *    nopt(37), lopt(55), lopt(57), lopt(58), lopt(59)
c                                 only vertex:
c                                 context specific parameters:
         if (icopt.le.3.and.(iam.eq.1.or.iam.eq.15)) then 
c                                 non-adaptive calculations
c                                 reaction format and lists
            if (icopt.gt.0) write (n,1160) grid(5,1),grid(5,2),rid(1,1),
     *                    rid(1,2),isec,valu(7),valu(9),valu(8),valu(10)
         else 
c                                 adaptive optimization
            write (n,1180) nopt(49),iopt(37),iopt(20),nopt(21),
     *                     valu(38),valu(39),lopt(62),iopt(31),k5,
     *                     lopt(49),lopt(54),nopt(48),nval2,nopt(9)
c                                 gridding parameters
            if (iam.eq.1.and.icopt.eq.5.and.oned) then
c                                 1d multilevel grid
               write (n,1190) grid(2,1),grid(2,2),l7,
     *                  (grid(2,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(2,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(3,1),grid(3,2),l8

            else if (iam.eq.1.and.icopt.eq.5) then
c                                 2d multilevel grid
               write (n,1200) grid(1,1),grid(1,2),l7,
     *                  (grid(1,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(1,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(2,1),grid(2,2),l7,
     *                  (grid(2,1)-1) * 2**(grid(3,1)-1) + 1,
     *                  (grid(2,2)-1) * 2**(grid(3,2)-1) + 1,
     *                  grid(3,1),grid(3,2),l8,valu(18)

            else if (iam.eq.1.and.icopt.eq.7) then 
c                                 1d fractionation grid
                write (n,1210) grid(4,1),grid(4,2),l7

            end if 
c                                 closed or open composition space
            if (iam.eq.1.and.icont.gt.1) write (n,1220) lopt(1)

         end if 
c                                 generic subdivision parameters:
         if (iam.eq.15) then 
            write (n,1011) nopt(13),nopt(13)/nopt(17),numb,nopt(14),
     *                  lopt(38),valu(13),valu(16),lopt(39),nopt(15)
         else 
            write (n,1010) nopt(13),nopt(14),
     *                     lopt(38),valu(13),lopt(39)
         end if 
c                                 generic thermo parameters:
         write (n,1012) nval1,nopt(12),nopt(20),
     *                  lopt(8),lopt(4),nopt(5),iopt(21),lopt(63),
     *                  iopt(25),iopt(26),iopt(27),valu(5),
     *                  lopt(32),lopt(44),lopt(36),lopt(46),nopt(34)
c                                 for meemum add fd stuff
         if (iam.eq.2) write (n,1017) nopt(31),nopt(26),nopt(27)

         if (iam.eq.1.or.iam.eq.15) then 
c                                 vertex output options, dependent potentials
c                                 pause_on_error
            write (n,1013) lopt(19),lopt(61)
c                                 auto_exclude, warn_interactive, 
c                                 output_iteration_details, output_iteration_g
            write (n,1234) lopt(5),lopt(56),lopt(64),lopt(33),lopt(34)
c                                 logarithmic_p, bad_number, interim_results
            if (iam.eq.1) write (n,1014) lopt(14),lopt(37),nopt(7),
     *                                   valu(34)

         end if 

      end if

      if (iam.eq.3) then 
c                                 WERAMI input/output options
         write (n,1230) lopt(25),iopt(32),l9,valu(26),valu(27),lopt(15),
     *                  lopt(14),lopt(37),nopt(7),lopt(22),valu(2),
     *                  valu(21),valu(3),lopt(41),lopt(42),lopt(45),
     *                  valu(4),lopt(6),valu(22),lopt(51),lopt(21),
     *                  lopt(24),
     *                  valu(14),lopt(19),lopt(20),valu(34),lopt(48)
c                                 WERAMI info file options
         write (n,1241) lopt(12)       
c                                 WERAMI thermodynamic options
         write (n,1016) lopt(8),lopt(4),iopt(25),iopt(26),iopt(27)
         write (n,1017) nopt(31),nopt(26),nopt(27)

      else if (iam.eq.2) then 
c                                 MEEMUM input/output options
         write (n,1231) lopt(25),iopt(32),l9,valu(26),valu(27),
     *                  lopt(14),lopt(37),nopt(7),lopt(22),valu(2),
     *                  valu(21),valu(3),lopt(6),valu(22),lopt(51),
     *                  lopt(21),lopt(24),valu(14),lopt(19),
     *                  lopt(20),lopt(61)
c                                 auto_exclude, warn_interactive, etc
         write (n,1234) lopt(5),lopt(56),lopt(64),lopt(33),lopt(34)

      else if (iam.eq.5) then 
c                                 FRENDLY input/output options
         write (n,1232) lopt(15),lopt(37),lopt(14),nopt(7),lopt(6),
     *                  lopt(19),.false.

      end if 
c                                 seismic property options
      if (iam.eq.2.or.iam.eq.3.or.iam.eq.5) write (n,1233) lopt(50),
     *         valu(19),
     *         nopt(6),lopt(17),valu(15),nopt(16),valu(14),lopt(20),
     *         .false.,lopt(65),nopt(65)

      if (iam.eq.5) then 
c                                 FRENDLY thermo options
         write (n,1016) lopt(8),lopt(4),iopt(25),iopt(26),iopt(27)
         write (n,1017) nopt(31),nopt(26),nopt(27)

      end if 

      if (iam.le.2) then 
c                                 info file options
         write (n,1240) lopt(12),lopt(10)
         if (iam.eq.1.or.iam.eq.15) write (n,1250) lopt(11)
         if (iam.eq.1) write (n,'(4x,a,l1,9x,a)') 
     *                    'seismic_data_file       ',lopt(50),'[F] T;'//
     *                    ' echo seismic wavespeed options'

      end if 

      write (n,1020) 

1000  format (/,'Perple_X computational option settings for ',a,':',//,
     *      '    Keyword:               Value:     Permitted values ',
     *          '[default]:')

1010  format (/,2x,'Solution subdivision options:',//,
     *        4x,'initial_resolution:     ',f6.4,4x,
     *                                   '[1/5] 0->1; 0 => off',/,
     *        4x,'stretch_factor          ',f6.4,4x,'[2d-3] >0 ',/,
     *        4x,'non_linear_switch       ',l1,9x,'[F] T',/,
     *        4x,'subdivision_override    ',a3,7x,'[lin] off str',/,
     *        4x,'refine_endmembers       ',l1,9x,'[F] T')

1011  format (/,2x,'Solution subdivision options:',//,
     *        4x,'initial_resolution:    ',/,
     *        4x,'  exploratory stage    ',f6.4,5x,
     *           '0->1 [1/16], 0 => off',/,
     *        4x,'  auto-refine stage    ',f6.4,5x,
     *           '0->1 [',a,'], 0 => off',/,
     *        4x,'stretch_factor         ',f6.4,5x,'>0 [2d-3]',/,
     *        4x,'non_linear_switch       ',l1,9x,'[F] T',/,
     *        4x,'subdivision_override    ',a3,7x,'[off] lin str',/,
     *        4x,'hard_limits             ',a3,7x,'[off] on',/,
     *        4x,'refine_endmembers       ',l1,9x,'[F] T',/,
     *        4x,'pc_perturbation        ',f6.4,5x,'[5d-3]')
c                                 generic thermo options
1012  format (/,2x,'Thermodynamic options:',//,
     *        4x,'solvus_tolerance        ',a7,3x,          
     *           '[aut] or 0->1; aut = automatic, 0 => ',
     *           'p=c pseudocompounds, 1 => homogenize',/,
     *        4x,'T_stop (K)           ',f6.1,7x,'[0]',/,
     *        4x,'T_melt (K)             ',f6.1,5x,'[873]',/,
     *        4x,'approx_alpha            ',l1,9x,'[T] F',/,
     *        4x,'Anderson-Gruneisen      ',l1,9x,'[F] T',/,
     *        4x,'speciation_precision   ',g7.1E1,4x,
     *           '[1d-5] <1; absolute',/,
     *        4x,'speciation_max_it      ',i4,7x,'[100]',/,
     *        4x,'GFSM                    ',l1,9x,
     *           '[F] T GFSM/special_component toggle',/,
     *        4x,'hybrid_EoS_H2O          ',i1,9x,'[4] 0-2, 4-7',/,
     *        4x,'hybrid_EoS_CO2          ',i1,9x,'[4] 0-4, 7',/,
     *        4x,'hybrid_EoS_CH4          ',i1,9x,'[0] 0-1, 7',/,
     *        4x,'aq_bad_results          ',a3,7x,'[err] 101 102 103',
     *                                           ' ignore',/,
     *        4x,'aq_lagged_speciation    ',l1,9x,'[F] T',/,
     *        4x,'aq_ion_H+               ',l1,9x,'[T] F => use OH-',/,
     *        4x,'aq_oxide_components     ',l1,9x,'[F] T',/,
     *        4x,'aq_solvent_solvus       ',l1,9x,'[T] F',/,
     *        4x,'aq_vapor_epsilon        ',f3.1,7x,'[1.]')
1013  format (/,2x,'Input/Output options:',//,
     *        4x,'pause_on_error          ',l1,9x,'[T] F',/,
     *        4x,'timing                  ',l1,9x,'[T] F')
1014  format (4x,'logarithmic_p           ',l1,9x,'[F] T',/,
     *        4x,'logarithmic_X           ',l1,9x,'[F] T',/,
     *        4x,'bad_number          ',f7.1,7x,'[NaN]',/,
     *        4x,'interim_results         ',a3,7x,'[auto] off manual')
1015  format (/,2x,'Auto-refine options:',//,
     *        4x,'auto_refine             ',a3,7x,'[auto] manual off',/,
     *        4x,'replicate_threshold    ',g7.1E1,4x,
     *           '[1e-2]; static opt; <0 => no replica test',/,
     *        4x,'rep_dynamic_threshold  ',g7.1E1,4x,
     *           '[1d-3]; dynamic opt; <0 => no replica test',/,
     *        4x,'re-refine               ',l1,9x,'[F] T',/,
     *        4x,'intermediate_savrpc     ',l1,9x,'[F] T',/,
     *        4x,'intermediate_savdyn     ',l1,9x,'[F] T',/,
     *        4x,'keep_all_rpcs           ',l1,9x,'[T] F')
c                                 thermo options for frendly
1016  format (/,2x,'Thermodynamic options:',//,
     *        4x,'approx_alpha            ',l1,9x,'[T] F',/,
     *        4x,'Anderson-Gruneisen      ',l1,9x,'[F] T',/,
     *        4x,'hybrid_EoS_H2O          ',i4,6x,'[4] 0-2, 4-7',/,
     *        4x,'hybrid_EoS_CO2          ',i4,6x,'[4] 0-4, 7',/,
     *        4x,'hybrid_EoS_CH4          ',i4,6x,'[0] 0-1, 7')
1017  format (4x,'fd_expansion_factor     ',f3.1,7x,'[2] >0',/,
     *        4x,'finite_difference_p     ',g7.1,3x,'[1d4] >0; ',
     *           'fraction = ',g7.1,3x,'[1d-2]')
1020  format (/,'To change these options see: ',
     *        'www.perplex.ethz.ch/perplex_options.html',/)
1100  format (/,2x,'Adapative minimization will be done with: ',
     *        //,3x,i2,' iterations in the exploratory stage',/,
     *           3x,i2,' iterations in the autorefine stage')
1160  format (/,2x,'Schreinemakers and Mixed-variable diagram ',
     *           'options:',//,
     *        4x,'variance               ',i2,' /',i2,5x,
     *           '[1/99], >0; maximum true variance',/,
     *        4x,'increment           ',f5.3,'/',f5.3,3x,
     *           '[0.1/0.025], ',
     *           'default search/trace variable increment',/,
     *        4x,'efficiency               ',i1,8x,'[3] >0, <6',/,
     *        4x,'reaction_format         ',a3,7x,'[min] ',
     *           'full stoichiometry S+V everything',/,
     *        4x,'reaction_list           ',a3,7x,'[off] on',/,
     *        4x,'console_messages        ',a3,7x,'[on] off',/,
     *        4x,'short_print_file        ',a3,7x,'[on] off')

1180  format (/,2x,'Free energy minimization options:',//,
     *        4x,'MINFRC_diff_increment  ',g7.1E1,4x,
     *           '[1e-7] 1e-3 => 1e-9',/,
     *        4x,'MINFXC_solver           ',i2,8x,
     *           '[0] >= 0 - speci2, -1 - MINFXC',/,
     *        4x,'optimization_max_it     ',i2,8x,'[40] >1',/,
     *        4x,'optimization_precision ',g7.1E1,4x,
     *           '[1e-4], 1e-1 => 1e-6, absolute',/,
     *        4x,'dynamic_LP_start        ',a3,7x,'[warm] cold hot',/,
     *        4x,'static_LP_start         ',a3,7x,'[hot] cold warm',/,
     *        4x,'order_check             ',l1,9x,'[F] T',/,
     *        4x,'refinement_points       ',i2,8x,'[auto] 1->',i2,/,
     *        4x,'refinement_switch       ',l1,9x,'[T] F',/,
     *        4x,'scatter-points          ',l1,9x,'[T] F',/,
     *        4x,'scatter-increment      ',g7.1E1,4x,
     *           '[1e-2] 1e-2 => 1e-7',/,
     *        4x,'solvus_tolerance_II     ',a7,3x,'[0.2] 0->1 ',/,
     *        4x,'zero_mode              ',e7.1E1,4x,
     *           '[1e-6] 0->1; < 0 => off')
1190  format (/,2x,'1D grid options:',//,
     *        4x,'y_nodes               ',i3,' /',i3,4x,'[40/40] >0, '
     *          ,'<',i4,'; effective y-resolution ',i4,' /',i4,
     *           ' nodes',/
     *        4x,'grid_levels             ',i1,' /',i2,5x,'[1/4] >0, '
     *          ,'<',i2,/)
1200  format (/,2x,'2D grid options:',//,
     *        4x,'x_nodes                ',i3,' /',i3,3x,'[10/40] >0, '
     *          ,'<',i4,'; effective x-resolution ',i4,' /',i4
     *          ,' nodes',/
     *        4x,'y_nodes                ',i3,' /',i3,3x,'[10/40] >0, '
     *          ,'<',i4,'; effective y-resolution ',i4,' /',i4,
     *           ' nodes',/
     *        4x,'grid_levels             ',i1,' /',i2,5x,'[1/4] >0, '
     *          ,'<',i2,/,
     *        4x,'linear_model            ',a3,7x,'[on] off')
1210  format (/,2x,'Fractionation path options:',//,
     *        4x,'1d_path               ',i3,' /',i3,4x,
     *           '[20/150] >0, <',i4)
1220  format (/,2x,'Composition options:',//,
     *        4x,'closed_c_space          ',l1,9x,'[T] F')
1230  format (/,2x,'Input/Output options:',//,
     *        4x,'aqueous_output          ',l1,9x,'[F] T',/
     *        4x,'aqeuous_species        ',i3,8x,'[20] 0-',i3,/,
     *        4x,'aq_solvent_composition  ',a3,7x,
     *        '[y] m: y => mol fraction, m => molality',/,
     *        4x,'aq_solute_composition   ',a3,7x,
     *        'y [m]: y => mol fraction, m => molality',/,
     *        4x,'spreadsheet             ',l1,9x,'[F] T',/,
     *        4x,'logarithmic_p           ',l1,9x,'[F] T',/,
     *        4x,'logarithmic_X           ',l1,9x,'[F] T',/,
     *        4x,'bad_number         ',f7.1,8x,'[NaN]',/,
     *        4x,'composition_constant    ',l1,9x,'[F] T',/,
     *        4x,'composition_phase       ',a3,7x,'[mol] wt',/,
     *        4x,'composition_system      ',a3,7x,'[wt] mol',/,
     *        4x,'proportions             ',a3,7x,'[vol] wt mol',/,
     *        4x,'absolute                ',l1,9x,'[F] T',/,
     *        4x,'cumulative              ',l1,9x,'[F] T',/,
     *        4x,'fancy_cumulative_modes  ',l1,9x,'[F] T',/,
     *        4x,'interpolation           ',a3,7x,'[on] off ',/,
     *        4x,'melt_is_fluid           ',l1,9x,'[T] F',/,
     *        4x,'solution_names          ',a3,7x,'[model] ',
     *                                           'abbreviation full',/,
     *        4x,'structural_formulae     ',l1,9x,'[T] F',/,
     *        4x,'output_species          ',l1,9x,'[T] F',/,
     *        4x,'output_species_props    ',l1,9x,'[F] T',/,
     *        4x,'seismic_output          ',a3,7x,'[some] none all',/,
     *        4x,'pause_on_error          ',l1,9x,'[T] F',/,
     *        4x,'poisson_test            ',l1,9x,'[F] T',/,
     *        4x,'interim_results         ',a3,7x,'[auto] off manual',/,
     *        4x,'sample_on_grid          ',l1,9x,'[T] F')
1231  format (/,2x,'Input/Output options:',//,
     *        4x,'aq_output               ',l1,9x,'[T] F',/
     *        4x,'aq_species             ',i3,8x,'[20] 0-',i3,/,
     *        4x,'aq_solvent_composition  ',a3,7x,
     *        '[y] m: y => mol fraction, m => molality',/,
     *        4x,'aq_solute_composition   ',a3,7x,
     *        'y [m]: y => mol fraction, m => molality',/,
     *        4x,'logarithmic_p           ',l1,9x,'[F] T',/,
     *        4x,'logarithmic_X           ',l1,9x,'[F] T',/,
     *        4x,'bad_number         ',f7.1,8x,'[NaN]',/,
     *        4x,'composition_constant    ',l1,9x,'[F] T',/,
     *        4x,'composition_phase       ',a3,7x,'[mol] wt',/,
     *        4x,'composition_system      ',a3,7x,'[wt] mol',/,
     *        4x,'proportions             ',a3,7x,'[vol] wt mol',/,
     *        4x,'melt_is_fluid           ',l1,9x,'[T] F',/,
     *        4x,'solution_names          ',a3,7x,'[mod] abb ful',/,
     *        4x,'structural_formulae     ',l1,9x,'[T] F',/,
     *        4x,'species_output          ',l1,9x,'[T] F',/,
     *        4x,'endmember_Gs            ',l1,9x,'[F] T',/,
     *        4x,'seismic_output          ',a3,7x,'[some] none all',/,
     *        4x,'pause_on_error          ',l1,9x,'[T] F',/,
     *        4x,'poisson_test            ',l1,9x,'[F] T',/,
     *        4x,'timing                  ',l1,9x,'[F] T')
1232  format (/,2x,'Input/Output options:',//,
     *        4x,'spreadsheet             ',l1,9x,'[F] T',/,
     *        4x,'logarithmic_p           ',l1,9x,'[F] T',/,
     *        4x,'logarithmic_X           ',l1,9x,'[F] T',/,
     *        4x,'bad_number         ',f7.1,8x,'[NaN]',/,
     *        4x,'melt_is_fluid           ',l1,9x,'[T] F',/,
     *        4x,'pause_on_error          ',l1,9x,'[T] F',/,
     *        4x,'Tisza_test              ',l1,9x,'[F] T')
1233  format (/,2x,'Seismic wavespeed computational options:',//,
     *        4x,'seismic_data_file       ',l1,9x,'[F] T',/,
     *        4x,'bounds                  ',a3,7x,'[VRH] HS',/,
     *        4x,'vrh/hs_weighting        ',f3.1,7x,'[0.5] 0->1',/,
     *        4x,'explicit_bulk_modulus   ',l1,9x,'[T] F',/,
     *        4x,'poisson_ratio           ',a3,7x,'[on] all off; ',
     *        'Poisson ratio = ',f4.2,/,
     *        4x,'seismic_output          ',a3,7x,'[some] none all',/,
     *        4x,'poisson_test            ',l1,9x,'[F] T',/,
     *        4x,'Tisza_test              ',l1,9x,'[F] T',/,
     *        4x,'fluid_shear_modulus     ',l1,9x,'[T] F',/,
     *        4x,'phi_d                   ',f4.2,6x,'[0.36] 0->1')
1234  format (4x,'auto_exclude            ',l1,9x,'[T] F',/,
     *        4x,'warn_interactive        ',l1,9x,'[T] F',/,
     *        4x,'warn_no_limit           ',l1,9x,'[F] T',/,
     *        4x,'output_iteration_detai  ',l1,9x,'[F] T',/,
     *        4x,'output_iteration_g      ',l1,9x,'[F] T')
1240  format (/,2x,'Information file output options:',//,
     *        4x,'option_list_files       ',l1,9x,'[F] T; ',
     *           'echo computational options',/,
     *        4x,'pseudocompound_file     ',l1,9x,'[F] T; ',
     *           'echo static pseudocompound compositions')
1241  format (/,2x,'Information file output options:',//,
     *        4x,'option_list_files       ',l1,9x,'[F] T; ',
     *           'echo computational options')
1250  format (4x,'auto_refine_file        ',l1,9x,'[T] F; ',
     *           'echo auto-refine compositions')

      end

      subroutine redcd1 (lun,ier,key,val,nval1,nval2,nval3,strg,strg1)
c----------------------------------------------------------------------
c this routine seeks a card containing a keyword and as many as 
c six values (char variables nval, nval1...), the first 3 letters of nval
c is also returned in val, if the second value is longer than 12 characters
c it is also saved in the character variable strg*80
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer lun, ier, iscan, iscnlt, ibeg, iend, ist, lend

      character card*(lchar), key*22, val*3,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
 
      ier = 0 
      key = ' '

      do 

         read (lun,'(a)',iostat=ier) card
         if (ier.ne.0) return

         if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')

            if (ibeg.ge.com) cycle
c                                 for programs (actcor,ctransf) 
c                                 that echo data read the full card
            length = iscnlt (lchar,1,' ')

            exit 

         end if 

      end do 
c                                 find end of keyword 
      iend = ibeg + 1
      iend = iscan (iend,lchar,' ') - 1

c      if (iend-ibeg.gt.21) then
c         call warn (99,0d0,ier,'invalid keyword in '
c     *   //'REDCD1, keywords must be < 23 characters.')
c         ier = 1
c         return 
c      end if 

      if (iend-ibeg.gt.21) then
         lend = ibeg + 21
      else 
         lend = iend
      end if 
c                                 load chars into key
      write (key,'(22a)') chars(ibeg:lend)
c                                 initialize other values
      strg = ' '
      strg1 = ' '
      nval1 = '0'
      nval2 = '0'
      nval3 = '0'

      iend = iend + 1
c                                 now locate the value:
      ibeg = iscnlt (iend,com,' ')
c                                 now find trailing blank
      iend = iscan (ibeg,lchar,' ') 
c                                 return if just a keyword
      if (iend.gt.lchar) return
c                                 look if it contains a comment character
      ist = iscan (ibeg,iend,'|') 
      if (ist.lt.iend) iend = ist - 1
c                                 save longer versions (only on first value)
c                                 this is done in case it's long text or 
c                                 several numbers on certain options. 
      if (iend-ibeg.gt.39) iend = ibeg+39
      write (strg,'(40a)') chars(ibeg:iend)
      write (strg1,'(40a)') chars(ibeg:ibeg+39)
c                                 read value:
      if (ibeg+2.gt.iend) then 
         ist = iend
      else 
         ist = ibeg + 2
      end if 

      write (val,'(3a)') chars(ibeg:ist)
c                                 look for a second value
      ist = iscan (ibeg,lchar,' ')
      if (ist.gt.com) return

      ibeg = iscnlt (ist,com,' ')
      if (ibeg.gt.com) return 

      iend = iscan (ibeg,com,' ')
      if (iend-ibeg.gt.11) iend = ibeg + 11 
      write (nval1,'(12a)') chars(ibeg:iend)
c                                 look for a third value
      ist = iscan (ibeg,lchar,' ')
      if (ist.gt.com) return 

      ibeg = iscnlt (ist,com,' ')
      if (ibeg.gt.com) return 

      iend = iscan (ibeg,com,' ')
      if (iend-ibeg.gt.11) iend = ibeg + 11 
      write (nval2,'(12a)') chars(ibeg:iend) 
c                                 look for a fourth value
      ist = iscan (ibeg,lchar,' ')
      if (ist.gt.com) return

      ibeg = iscnlt (ist,com,' ')
      if (ibeg.gt.com) return

      iend = iscan (ibeg,com,' ')
      if (iend-ibeg.gt.11) iend = ibeg + 11 
      write (nval3,'(12a)') chars(ibeg:iend)

      end

      subroutine rdstrg (lun,nstrg,string,eof)
c----------------------------------------------------------------------
c rdstrg - read 240 column card images from unit lun until a non-blank
c (i.e., with data other than comments) record, then read up to three
c strings from the record. on output nstrg is the number of strings read
c from the record. 
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer lun, iscan, iscnlt, ibeg, iend, ier, nstrg, imax

      external iscan, iscnlt

      logical eof

      character card*(lchar), string(3)*8

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      eof = .false.

      do 
c                                 read cards till a non-blank
c                                 card.
         read (lun,'(a)',iostat=ier) card

         if (ier.ne.0) then
c                                 error on read = eof
            eof = .true.

            return

         else if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1

            if (com.eq.0) cycle 
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')

            exit 

         end if 

      end do 
c 
c                                 we have a non-blank card
      nstrg = 1

      do 
c                                 find the end of the string
         iend = iscan (ibeg,lchar,' ') - 1

         if (iend-ibeg.gt.7) then

            imax = ibeg + 7 

         else

            imax = iend
 
         end if 
c                                 load chars into string
         write (string(nstrg),'(8a)') chars(ibeg:imax)
c                                 find the next string
         ibeg = iscnlt (iend+1,com,' ')

         if (ibeg.gt.com.or.nstrg.eq.3) return
 
         nstrg = nstrg + 1

      end do 

      end

      subroutine rdnumb (numb,def,inumb,idef,reel)
c----------------------------------------------------------------------
c rdnumb - reads a line from terminal input for numeric input, if blank
c assigns default (def, idef); if non numeric prompts for new value.
c----------------------------------------------------------------------    
      implicit none

      integer inumb, idef, ier

      double precision numb, def

      logical defalt, reel

      character card*80
c----------------------------------------------------------------------
      defalt = .true.

      do 
c                                 read input
         read (*,'(a)',iostat=ier) card
c                                 user enters a blank
         if (ier.ne.0.or.card.eq.' ') exit
c                                 read data from card
         if (reel) then 
            read (card,*,iostat=ier) numb
         else 
            read (card,*,iostat=ier) inumb
         end if 

         if (ier.ne.0) then 
            call rerr 
         else
            defalt = .false.
            exit 
         end if 

      end do 

      if (defalt) then 
         if (reel) then 
            numb = def
         else 
            inumb = idef
         end if 
      end if 

      end

      subroutine eohead (n)
c----------------------------------------------------------------------
c eohead reads cards from n until an 'END ' or 'end ' is found in
c the first 4 columns
c----------------------------------------------------------------------
      implicit none

      integer n, ier

      character tag*4

      rewind n

      do 
         read (n,'(a)',iostat=ier) tag
         if (ier.ne.0) call error (37,1d0,n,'EOHEAD')
         if (tag.eq.'end'.or.tag.eq.'END') exit
      end do 

      end

      subroutine rerror (ier,*)
c---------------------------------------------------------------------
c rerror - routine to check for errors during list directed i/o
 
      implicit none

      integer ier
c---------------------------------------------------------------------
 
      if (ier.eq.0) then
         return
      else
         write (*,1000)
         ier = 0
         return 1
      end if
 
1000  format (/,'Your input is incorrect, probably you have specified ',
     *        'an invalid numerical value',/,'or you are using ',
     *        'a character where you should be using a number ',
     *        'or vice versa.',/,'try again...',/)
 
      end

      subroutine rerr 
c---------------------------------------------------------------------
c rerror - routine to write bad input message for interactive i/o
 
      implicit none
c---------------------------------------------------------------------

      write (*,1000)
 
1000  format (/,'Your input is incorrect, probably you are using ',
     *        'a character where',/,'you should be using a number ',
     *        'or vice versa, try again...',/)
 
      end

      subroutine errdbg (msg)
c---------------------------------------------------------------------
c print debugging error msg and pause.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character msg*(*)

      write (*,'(/,a,/)') msg

      call errpau

      end 

      subroutine errpau 
c---------------------------------------------------------------------
c pause on error option
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character a*1
      
      if (lopt(19)) then 
         write (*,'(/,a,/)') 'Press Enter to quit...' 
         read (*,'(a)') a
      end if 
      
      stop
      
      end 

      subroutine error (ier,realv,int,char)
c---------------------------------------------------------------------
c write error messages and terminate execution
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, int
 
      character char*(*), tag*3

      double precision realv

      integer iam
      common/ cst4 /iam

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      if (refine) then
         tag = '2nd'
      else 
         tag = '1st'
      end if

      if (ier.eq.1.or.ier.eq.2) then 
         write (*,1) char,int
      else if (ier.eq.3) then 
         write (*,3) char
      else if (ier.eq.4) then
         write (*,4) char 
      else if (ier.eq.5) then 
         write (*,5) int,char,j3
      else if (ier.eq.6) then 
         write (*,6) char
      else if (ier.eq.7) then 
         write (*,7) char
      else if (ier.eq.8) then 
         write (*,8) 
      else if (ier.eq.9) then 
         write (*,9) char
      else if (ier.eq.10) then 
         write (*,10) char, int
      else if (ier.eq.11) then 
         write (*,11) char,int
      else if (ier.eq.12) then 
         write (*,12) char
      else if (ier.eq.13) then
         write (*,13) h8
      else if (ier.eq.14) then
         write (*,14) char
      else if (ier.eq.15) then
         write (*,15) char
      else if (ier.eq.16) then
         write (*,16) h5
      else if (ier.eq.17) then
         write (*,17) int
      else if (ier.eq.18) then
         write (*,18) char
      else if (ier.eq.19) then
         write (*,19) char
      else if (ier.eq.20) then
         write (*,20) int, char
      else if (ier.eq.21) then
         write (*,21) char 
      else if (ier.eq.22) then
         write (*,22) int, char
      else if (ier.eq.23) then
         write (*,23) char
      else if (ier.eq.24) then
         write (*,24) int
      else if (ier.eq.25) then
         write (*,25) h9
      else if (ier.eq.26) then
         write (*,26) int, char
      else if (ier.eq.27) then
         write (*,27) char
      else if (ier.eq.28) then 
         write (*,28) int, char
      else if (ier.eq.29) then 
         write (*,29) int, char
      else if (ier.eq.32) then 
         write (*,32)
      else if (ier.eq.33) then 
         write (*,33) char, int
      else if (ier.eq.34) then
         write (*,34)
      else if (ier.eq.35) then
         write (*,35)
      else if (ier.eq.36) then
         write (*,36)
      else if (ier.eq.37) then
         write (*,37) int
      else if (ier.eq.38) then 
         write (*,38) 
      else if (ier.eq.39) then
         write (*,39) int
      else if (ier.eq.40) then
         write (*,40) int, char
      else if (ier.eq.41) then

         write (*,41)

         if (int.eq.0) then 
            write (*,410)
         else if (int.eq.1) then 
            write (*,411)
         else if (int.eq.2) then 
            write (*,*) ' are you kidding me? please report'
         end if

         write (*,414) char
         write (*,413)

      else if (ier.eq.42) then 
         write (*,42) char
      else if (ier.eq.43) then
         write (*,43) char
      else if (ier.eq.44) then 
         write (*,44) 
      else if (ier.eq.45) then 
         write (*,45) 
      else if (ier.eq.46) then 
         write (*,46) iopt(16), int, char
      else if (ier.eq.47) then 
         write (*,47) char
      else if (ier.eq.48) then 
         write (*,48) char,int
      else if (ier.eq.49) then 
         write (*,49) char,int
      else if (ier.eq.50) then 
         write (*,50) realv, char, int
      else if (ier.eq.51) then
         write (*,51) char
      else if (ier.eq.52) then 
         write (*,52) h9
      else if (ier.eq.53) then 
         write (*,53) 
      else if (ier.eq.54) then 
         write (*,54)
      else if (ier.eq.55) then 
         write (*,55) k16
      else if (ier.eq.56) then 
         write (*,56) k17
      else if (ier.eq.57) then 
         write (*,57) char, char, char
      else if (ier.eq.60) then 
         write (*,60) k22, char
      else if (ier.eq.61) then 
         write (*,61) k18, char
      else if (ier.eq.62) then
         write (*,62) char, int, realv
      else if (ier.eq.63) then 
         write (*,63) char
      else if (ier.eq.64) then
         write (*,64) char
      else if (ier.eq.65) then
         write (*,65) char
      else if (ier.eq.66) then
         write (*,66) 
      else if (ier.eq.67) then
         write (*,67) char
      else if (ier.eq.68) then
         write (*,68) char
      else if (ier.eq.69) then 
         write (*,69) char
      else if (ier.eq.72) then 
         write (*,72) char
      else if (ier.eq.73) then 
         write (*,73) char
      else if (ier.eq.74) then 
         write (*,74) int
      else if (ier.eq.75) then 
         write (*,75) char
      else if (ier.eq.76) then 
         write (*,76) char, char
      else if (ier.eq.78) then 
         write (*,78) char,char
      else if (ier.eq.89) then
         write (*,89) 
      else if (ier.eq.90) then
         write (*,90) l6
      else if (ier.eq.106) then
         write (*,106) char
      else if (ier.eq.107) then
         write (*,107) int
      else if (ier.eq.108) then
         write (*,108) int
      else if (ier.eq.109) then
         write (*,109) int
      else if (ier.eq.110) then
         write (*,110)
      else if (ier.eq.111) then
         write (*,111)
      else if (ier.eq.112) then
         write (*,112) char
      else if (ier.eq.116) then
         write (*,116)
      else if (ier.eq.117) then
         write (*,117)
      else if (ier.eq.118) then
         write (*,118)
      else if (ier.eq.120) then
         write (*,120) char
      else if (ier.eq.125) then 
         write (*,125) char, int
      else if (ier.eq.169) then
         write (*,169) int
      else if (ier.eq.180) then
c                                 static
            write (*,180) k13,int,tname,tag
            write (*,1803)
            write (*,1804)

         if (realv.gt.0d0) then 
c                                 non-linear schemes are in use, instruct
c                                 accordingly:
            if (.not.lopt(38).and.iam.ne.15) then
               write (*,1805) tname,'T (true).'
            else if (lopt(38).and.iam.eq.15) then
               write (*,1805) tname,'F (false).'
            end if

         end if

      else if (ier.eq.181) then
         write (*,181) int
      else if (ier.eq.182) then
         write (*,182) k2
      else if (ier.eq.183) then
         write (*,183) k2,char
      else if (ier.eq.197) then
         write (*,197) int, k5, char
      else if (ier.eq.200) then
         write (*,200)
      else if (ier.eq.204) then
         write (*,204) int
      else if (ier.eq.206) then
         write (*,206) int
      else if (ier.eq.207) then
         write (*,207) realv,char
      else if (ier.eq.208) then
         write (*,208) char
      else if (ier.eq.227) then
         write (*,227) char, int
      else
         write (*,999) ier, realv, int, char
      end if
 
      call errpau

1     format (/,'**error ver001** increase parameter ',a,' to ',i7,' in'
     *       ,' perplex_parameters.h and recompile Perple_X',/)
2     format (/,'**error ver002** solution model file versions ',
     *         '6.8.2-6.8.7 are inconsistent with',/,
     *         'this version of Perple_X. Update the solution ',
     *         'file.',/)
3     format (/,'**error ver003** the solution model file format (',a,
     *         ') is inconsistent with',/,
     *         'this version of Perple_X. Update the file and/or ',
     *         'Perple_X',/)
4     format (/,'**error ver004** you must use ',a,' to analyze this '
     *         ,'type of calculation.',/)
5     format (/,'**error ver005** too many ordered species (',i2,') in',
     *        ' solution model ',a,/,'increase dimension j3 (',i2,')',/)
6     format (/,'**error ver006** fractionation path coordinate file: '
     *          ,a,/,'does not exist.',/)
7     format (/,'**error ver007** reference phase ',a,' is not in the ',
     *          'thermodynamic data file.',/)
8     format (/,'**error ver008** the thermodynamic data file ',
     *          'is out of date, download',/,'the current ',
     *          'version from: www.perplex.ethz.ch',/)
9     format (/,'**error ver009** invalid tag (',a,') in the ',
     *          'thermodynamic data file.',/)
10    format (/,'**error ver010** the text string:',/,a,/,'is too long',
     *        ' to be processed.',/,'The maximum allowed length is ',i3,
     *        ' characters.',/)
11    format (/,'**error ver011** invalid ',a,' choice (',i3,')',a,/)
12    format (/,'**error ver012** file:',a,/,'is missing or formatted ',
     *        'incorrectly, run paralyzer or create/edit it manually',/)
13    format ('**error ver013** too many excluded phases, ',
     *        'increase dimension h8 (',i3,')')
14    format ('**error ver014** programming error, routine ',a)
15    format (/,'**error ver015** missing composant for: ',a,/)
16    format (/,'**error ver016** too many saturated components, ',
     *        'increase dimension h5 (',i2,')')
17    format (/,'**error ver017** too many composants for a saturation',
     *        ' constraint increase dimension h6 (',i3,')')
18    format (/,'**error ver018** ',a,' is defined as a saturated ',
     *        'phase component in the thermodynamic data file.')
19    format (/,'**error ver019** probable cause missing composant,',
     *        ' executing routine ',a)
20    format (/,'**error ver020** error ',i2,' reading solution model',
     *        ' file.',/,'   Reading model: ',a,' Check format.',/)
21    format (/,'**error ver021**error reading ',
     *        'header section of',/,'the thermodynamic data ',
     *        'file, last data read:',/,a,/,'Check formatting',/)
22    format (/,'**error ver022** too many divariant assemblages, ',
     *        'increase dimension j9 (',i8,') routine: ',a)
23    format (/,'**error ver023**error reading',
     *        ' thermodynamic data file.',/,'Last data read:',
     *        /,a,/,'Check formatting.',/)
24    format (/,'**error ver024** too many solution models in',
     *        ' solution model file',/,' increase parameter i9 (',
     *        i3,')')
25    format (/,'**error ver025** too many solution models ',
     *          'increase parameter h9 (',i3,')')
26    format (/,'**error ver026** the number of fixed components (',
     *        i2,') in ',a,/,' is >= the number of components ',/)
27    format (/,'**error ver027** Error reading the problem',
     *        ' definition file:',//,a,//,
     *        'Probable cause: You are using an ',
     *        'input file created by an out-of-date',/,
     *        '                version of BUILD, or you have',
     *        ' incorrectly edited the',/'                input file',
     *        ' created by BUILD',/)
28    format (/,'**error ver028** invalid buffer choice (',i3,') in',
     *          ' routine: ',a,/)
29    format (/,'**error ver029** unknown term type ',i6,' for',
     *          ' solution model: ',a,/)
32    format (/,'**error ver032** stability field calculations (',
     *          'option 2) are disabled in this version of PERPLEX',/)
33    format (/,'**error ver033** expression with too many terms in ',a
     *       ,/,'increase m0 or j6 to',i2,'.',/)
34    format (/,'**error ver034** vmax is lt vmin, check input.')
35    format (/,'**error ver035** dv is lt 0, check input.')
36    format (/,'**error ver036** missing composant for the saturated',
     *        ' phase,',/,'you have probably excluded either H2O or',
     *        ' CO2,',/,'or a composant is duplicated in the',
     *        ' thermodynamic data file',/)
37    format (/,'**error ver037** no end marker in header',/,
     *        'section of thermodynamic data file unit ',i2,/)
38    format (/,'**error ver038** you have configured a ',
     *       'problem with only one independent variable.',/,
     *       'This case cannot be handled by constrained minimization',
     *       ' use the unconstrained computational mode.'/)
39    format (/,'**error ver039** too many end-members, ',
     *        'increase dimension k12 (',i2,') Routine: ',a)
40    format (/,'**error ver040** too many compositional coordinates, ',
     *        'increase dimension k13 (',i7,')  Routine: ',a)
41    format (/,'**error ver041** too many static compositions this ',
     *          'error can be eliminated',/,'by one of the ',
     *          'following actions (best listed first):',/)
410   format (2x,'- increase exploratory stage initial_resolution',
     *           ' in perplex_option.dat',/,
     *        2x,'- restrict subdivision ranges of complex solutio',
     *           'ns in solution_model.dat')
411   format (2x,'- increase auto-refine stage initial_resolution ',
     *           'in perplex_option.dat')
413   format (2x,'- simplify the calculation, e.g., eliminate ',
     *           'components and/or simplify solution models')
414   format (2x,'- increase parameter ',a,' and recompile ',
     *           'Perple_X')
42    format (/,'**error ver042** cannot open file:',a,/,'check that it'
     *       ,' is not being used by another program',/)
43    format (/,'**error ver043** you cannot simultaneously treat: ',
     *          a,/,'as a thermodynamic solution and as a saturated',
     *          ' phase.',/)
44    format (/,'**error ver044** too many saturated phase components.'
     *        /)
45    format (/,'**error ver045** too many mobile components.'/)
46    format (/,'**error ver046** the first value of the iteration ',
     *          'keyword exceeds (',i2,') the value',/,'of MRES (',i3,
     *          ') specified in routine ',a,'. Either reduce the',/,
     *          'iteration keyword value or increase MRES.',/) 
47    format (/,'**error ver047** solution model ',a,' is incorrectly ',
     *        'formatted (van Laar).',/)
48    format (/,'**error ver048** too many terms in solution model ',a,
     *        ' increase parameter m1 (',i2,').',/)
49    format (/,'**error ver049** the order of solution model ',a,
     *        ' is too high, increase parameter m2 (',i2,').',/)
50    format (/,'**error ver050** requested resolution ',
     *          '(',f6.0,') for a component in solution:',a,/,
     *          'exceeds 1/MRES (MRES=',i5,') ',
     *          'reduce requested resolution or inrease',/,
     *          'MRES in routine CARTES',/)
51    format (/,'**error ver051** DUMMY1 could not find the auxilliary'
     *         ,' input file:',/,a,/,'required for open system model ',
     *          'computations (ICOPT=9).',/)
52    format (/,'**error ver052** too many solution models in your'
     *         ,' calculation',/,'reduce the number of models or ',
     *          'increase parameter H9 (',i2,').',/)
53    format (/,'**error ver053** phase fractionation calculations '
     *         ,'require >1 thermodynamic component.',/)
54    format (/,'**error ver054** unanticipated condition, probable ',
     *          'cause is incorrect ordering of',/,'endmembers in the',
     *          ' solution model, which leads to inconsistent site ',
     *          'occupancies',/)
55    format (/,'**error ver055** too many make definitions, delete '
     *         ,'unused definitions from the',/
     *         ,'thermodynamic data file or '
     *         ,'increase parameter K16 (',i3,') and recompile.',/)
56    format (/,'**error ver056** too many phases in a make definition'
     *         ,', increase parameter K17 (',i2,') and recompile.',/)
57    format (/,'**error ver057** failed on an accepted make definition'
     *         ,' for ',a,/,'routine INPUT2. Exclude ',a,' and restart',
     *          ' the calculation.',/,'If ',a,/,' is legitimate. please'
     *         ,' report this error.',/)
60    format (/,'**error ver060** too many coordinates generated by ',
     *        'refinement, increase dimension k22 (',i8,') routine: ',a)
61    format (/,'**error ver061** too many solution coordinates, ',
     *        'increase dimension k18 (',i8,') routine: ',a)
62    format (/,'**error ver062** solution model ',a,' specifies non-',
     *          'Cartesian subdivision (',i1,')',/,' and must be refor',
     *          'mulated for adapative minimization, but VERTEX cannot',
     *        /,' do the reformulation because the initial_reolution ',
     *          'keyword specified in',/,' perplex_option.dat (',f5.2,
     *          ') is invalid',/)
63    format (/,'**error ver063** inconsistent auto-refine data: ',a,/,
     *        'Suppress or reinitialize auto-refinement.',/) 
64    format (/,'**error ver064** PSVDRAW plots only ',
     *          'binary mixed-variable and ',/,
     *          'ternary composition diagrams (',a,').',/)
65    format (/,'**error ver065** dimensioning error (',a,').',/)
66    format (/,'**error ver066** invalid format, most probably this',
     *          ' result should be plotted with PSSECT.',/)
67    format (/,'**error ver067** file ',a,/,
     *        'is not formatted correctly for PSVDRAW.',/)
68    format (/,'**error ver068** solution model: ',a,
     *          ' is in a format that is no longer supported',/,
     *          'Use a more recent solution model file, e.g., copy ',
     *          'the current version from: ',//,
     *          'www.perplex.ethz.ch/datafiles/solution_model.dat',/)
69    format (/,'**error ver069** too many points (',a,'), increase ',
     *          'parameter L5',/)
72    format (/,'**error ver072** ',a,/)
73    format (/,'**error ver073** the thermodynamic data file has ', 
     *          'more than one entity named: ',a,/,'delete or rename ',
     *          'the entities, this error is often caused by make ',
     *          'definitions.')
74    format (/,'**error ver074** unrecognized EoS pointer (',i3, 
     *          ') in routine SETINS',/)
75    format (/,'**error ver075** more than one solution model is ',
     *          'named ',a,/,'delete or rename the replicate models in',
     *          ' the solution model file.',/)
76    format (/,'**error ver076** the ',a,' solution model was not ',
     *          'reformulated correctly.',//,
     *          'If this error was preceeded by warning ver114 for a ',
     *          'dependent endmember, then',/,'it is possible that ',
     *          'the endmember definition is missing from the model.',
     *           //,
     *          'Otherwise excluding either more or fewer ',a,
     *          ' endmembers may bypass/rectify this error.',/)
78    format (/,'**error ver078** ',a,' has dependent endmembers with ',
     *        'invalid site populations',/,'it cannot be used unless ',
     *        'it is corrected or the site_check_override keyword is',/,
     *        'specified at the end of the ',a,' model.',/)
89    format (/,'**error ver089** SMPLX programming error. Change ',
     *        'minimnization method.',/)
90    format (/,'**error ver090** SMPLX failed to converge within ', 
     *        i6,' iterations.',/,'Probable cause: the possible ',
     *        'phases do not span the systems composition',/,
     *        'To avoid this problem add phases or modify the bulk ',
     *        'composition.',/,'Alternatively, although less ',
     *        'probably, increasing parameter L6 in perplex_',
     *        'parameters.h',/,
     *        'and recompiling VERTEX permit SMPLEX to converge.',/)
106   format (/,'**error ver106** programming error in ',a)
107   format (/,'**error ver107** the assemblage you input is ',
     *        'metastable (ier=',i3,').')
108   format (/,'**error ver108** the assemblage you input is ',
     *        'degenerate (ier=',i3,').')
109   format (/,'**error ver109** the assemblage you input does not '
     *       ,'span the specified bulk composition (ier=',i3,').')
110   format (/,'**error ver110** you have requested a calculation ',
     *        'with the composition ',/,'of a saturated phase as a ',
     *        'variable, but you have not defined its composition')
111   format (/,'**error ver111** you have requested a calculation ',
     *        'with the composition ',/,'of a saturated phase as a ',
     *        'variable, but the phase has only one component')
112   format (/,'**error ver112** the maximum value of an independent '
     *       ,'variable',/,'is ',a,' to the minimum value')
116   format (/,'**error ver116** an independent variable, or at least'
     *       ,' its name, is undefined')
117   format (/,'**error ver117** vmax(iv(3)) ne vmin(iv(3) but no ',
     *        'sectioning variable v(iv(3)) is defined')
118   format (/,'**error ver118** the default increment of the ',
     *        'sectioning variable will result ',/,
     *        'in the generation of more ',
     *        'than 10 sections, to avoid this',/,' error increase ',
     *        'the increment or modify this test')
120   format (/,'**error ver120** file:',/,a,/,
     *        'could not be opened, check that it exists or that it is',
     *        ' not in use by another program.',/) 
125   format (/,'**error ver125** a site fraction is out of range for ',
     *          'solution: ',a,' endmember ',i2,/,'The configurational',
     *          ' entropy model is probably incorrect.',/)
169   format (/,'**error ver169** cart, imod=',i2,' is an invalid ',
     *          'request')
180   format (/,'**error ver180** too many pseudocompounds ',i8,' gener'
     *       ,'ated subdividing site ',i1,' of solution: ',a,/,
     *        'this error can usually be eliminated by one of the ',
     *        'following actions (best listed first):',//,
     *        2x,'- increase the ',a,' value of the initial_resolution '
     *          ,'keyword in perplex_option.dat')
1803  format (2x,'- restrict the subdivision range for the solution')
1804  format (2x,'- if nonlinear subdivision is specified for the solu'
     *          ,'tion then increase the',/,4x,'resolution parameters '
     *          ,'or change to linear subdivision',/,
     *        2x,'- increase parameter k13 and recompile',/)
1805  format ('NOTE!! ',a,' specifies nonlinear subdivision sch',
     *        'emes, the number of pseudo-',/,'compounds generated by ',
     *        'nonlinear subdivision is independent of initial_resolut',
     *        'ion',/,'unless the non_linear_switch option is ',a,
     *        ' For additional information on',/,'nonlinear subdivi',
     *        'sion refer to the commentary in the header of the ',
     *        'solution model',/,'file. If you do not understand nonli',
     *        'near subdivision use linear subdivision.',/)
181   format (/,'**error ver181** too many reactions,',
     *          ' increase dimension k2 (',i6,')')
182   format (/,'**error ver182** too many invariant points,',
     *           ' increase parameter k2 (',i6,')')
183   format (/,'**error ver183** too many assemblages; increase ',
     *        ' parameter k2 (',i6,'), routine ',a)
197   format (/,'**error ver197** to many components (',i2,'), increase'
     *         ,' parameter k5 (',i2,'), routine:',a,/)
200   format (/,'**error ver200** you are trying to use a fluid ',
     *        'equation of state for an invalid component',/)
204   format ('**error ver204** too many stable assemblages, i',
     *        'ncrease dimension j9 (',i8,')',/)
206   format ('**error ver206** too many univariant assemblages ',
     *        'intersect the edges of the diagram, i',
     *        'ncrease dimension k2 (',i6,')',/)
207   format (/,'**error ver207** the value of the stretching ',
     *          ' parameter (',g13.6,')',/,'for solution ',a,
     *          ' is invalid (<1) for transform subdivision,',/,
     *          'check section 4 of PERPLEX documentation.',/)
208   format (/,'**error ver208** too many phases on one side of a',/
     *        ' reaction.',/,'Do not use the full reaction',
     *        ' equation option (',a,').')
227   format (/,'**error ver227** in solution model ',a,' a DQF ',
     *          'correction is specified for endmember: ',i2,/,
     *          'DQF corrections can only be made on the idependent ',
     *          'endmembers of a solution model',/)
999   format (/,'**error vertex** unspecified error ier=',i3,/,
     *        ' real=',g15.7,/,' i=',i12,/,' char=',a)
      end

      subroutine warn (ier,realv,int,char)
c---------------------------------------------------------------------
c write warning message and continue execution

c generic warnings:

c 49  - future instances of warning int will not be repeated.
c 99  - just dump char
c 100 - use int (i3) as error number and dump char.
c 999 - unspecified real, int, char dump.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      integer ier,int

      double precision realv
 
      character char*(*)

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)
c----------------------------------------------------------------------
      if (ier.eq.1) then 
         write (*,1) 
      else if (ier.eq.2) then 
         write (*,2) realv
      else if (ier.eq.3) then 
         write (*,3)
      else if (ier.eq.4) then 
         write (*,4) char
      else if (ier.eq.5) then
         write (*,5) 
      else if (ier.eq.6) then
         write (*,6) 
      else if (ier.eq.7) then
         write (*,7) 
      else if (ier.eq.8) then
         write (*,8) h8
      else if (ier.eq.9) then
         write (*,9) char
      else if (ier.eq.10) then
         write (*,10) int, realv, char
      else if (ier.eq.11) then
         write (*,11) char
      else if (ier.eq.12) then
         write (*,12) char
      else if (ier.eq.13) then
         write (*,13) char, char
      else if (ier.eq.14) then
         write (*,14) char
      else if (ier.eq.15) then
         write (*,15)
      else if (ier.eq.16) then
         write (*,16) char
      else if (ier.eq.17) then
         write (*,17) char, char
      else if (ier.eq.18) then
         write (*,18) realv
      else if (ier.eq.19) then
         write (*,19) 
      else if (ier.eq.20) then
         write (*,20)
      else if (ier.eq.21) then
         write (*,21) realv, char
      else if (ier.eq.22) then
         write (*,22) realv, char
      else if (ier.eq.23) then
         write (*,23) char     
      else if (ier.eq.24) then
         write (*,24) realv
      else if (ier.eq.25) then 
         write (*,25) int, char
      else if (ier.eq.26) then 
         write (*,26) char
      else if (ier.eq.27) then 
         write (*,27) int
      else if (ier.eq.28) then
         write (*,28)
      else if (ier.eq.29) then
         write (*,29) char
      else if (ier.eq.30) then
         write (*,30) char
      else if (ier.eq.31) then 
         write (*,31)
      else if (ier.eq.32) then
         write (*,32) char
      else if (ier.eq.33) then
         write (*,33) char
      else if (ier.eq.34) then
         write (*,34) char
      else if (ier.eq.35) then
         write (*,35) char,realv
      else if (ier.eq.36) then 
         write (*,36) realv, char 
      else if (ier.eq.37) then
         write (*,37)  
      else if (ier.eq.38) then
         write (*,38) 
      else if (ier.eq.39) then
         write (*,39) 
      else if (ier.eq.40) then
         write (*,40)
      else if (ier.eq.41) then
         write (*,41) char, int
         call prtptx
      else if (ier.eq.42) then
         write (*,42)
      else if (ier.eq.43) then
         write (*,43) char
      else if (ier.eq.44) then
         write (*,44) char
      else if (ier.eq.45) then
         write (*,45) char
      else if (ier.eq.46) then 
         write (*,46) realv, char, char
      else if (ier.eq.47) then
         write (*,47) int, realv
      else if (ier.eq.48) then 
         write (*,48) 
      else if (ier.eq.49) then 
         write (*,49) int, char
      else if (ier.eq.50) then
         write (*,50) char
      else if (ier.eq.51) then 
         write (*,51) char
      else if (ier.eq.52) then 
         write (*,52) char
      else if (ier.eq.53) then 
         write (*,53) realv
      else if (ier.eq.54) then 
         write (*,54)
      else if (ier.eq.55) then 
         write (*,55) char
      else if (ier.eq.56) then 
         write (*,56) char
      else if (ier.eq.57) then
         write (*,57) char
      else if (ier.eq.58) then

         write (*,58)
         write (*,582)
         if (lopt(49)) write (*,584)
         write (*,585)
         if (lopt(32)) write (*,586)
         write (*,413)
         write (*,580) char

58    format (/,'**warning ver058** exhausted memory during adaptive o',
     *       'ptimization, the result will',/,'be reported as a failed',
     *       ' optimization, if such failures are excessive the ',/,
     *       'problem may be eliminated by the ',
     *       'following actions (best listed first):',/)
580   format (2x,'- increase parameter ',a,' and recompile ',
     *           'Perple_X')
582   format (2x,'- set keep_auto to T or default ',
     *           'in perplex_option.dat')
584   format (2x,'- set refinement_switch to F ',
     *           'in perplex_option.dat')
585   format (2x,'- reduce refinement_points (< c+2, > 0) ',
     *           'in perplex_option.dat')
586   format (2x,'- set aq_solvent_solvus to F ',
     *           'in perplex_option.dat')
413   format (2x,'- simplify the calculation, e.g., eliminate ',
     *           'components and/or simplify solution models')

      else if (ier.eq.589) then
         write (*,589)
      else if (ier.eq.59) then
         write (*,59) char
      else if (ier.eq.60) then
         write (*,60) char
         if (int.eq.1) then 
            write (*,601) char
         else 
            write (*,602)
         end if
      else if (ier.eq.61) then
         write (*,61) char
      else if (ier.eq.62) then
         write (*,62) char
      else if (ier.eq.63) then
         write (*,63)
      else if (ier.eq.64) then
         write (*,64) realv
      else if (ier.eq.68) then
         write (*,68)
      else if (ier.eq.72) then
c                                 generic warning, also 99
         write (*,72) char
      else if (ier.eq.73) then
         write (*,73) char, realv, int
      else if (ier.eq.74) then
         write (*,74)
      else if (ier.eq.79) then
         write (*,79) char
      else if (ier.eq.87) then
         write (*,87)
      else if (ier.eq.88) then
         write (*,88)
      else if (ier.eq.89) then
         write (*,89)
      else if (ier.eq.90) then
         write (*,90) 
      else if (ier.eq.91) then
         write (*,91)
      else if (ier.eq.92) then 
         write (*,92) int, l7, char, (l7 - 1)/2**(grid(3,2)-1) + 1
      else if (ier.eq.99) then
         write (*,99) char
      else if (ier.eq.100) then
         write (*,100) int, char
      else if (ier.eq.106) then
         write (*,106) char
      else if (ier.eq.108) then
         write (*,108)
      else if (ier.eq.109) then
         write (*,109)
      else if (ier.eq.113) then
         write (*,113) int
      else if (ier.eq.114) then
         write (*,114)
      else if (ier.eq.172) then
         write (*,172) 
      else if (ier.eq.173) then
         write (*,173) 
      else if (ier.eq.175) then
         write (*,175) char,ier,realv
      else if (ier.eq.176) then
         write (*,176) char, iopt(21)
      else if (ier.eq.177) then
         write (*,177) nopt(5)
      else if (ier.eq.205) then
         write (*,205) int
         write (*,900)
      else if (ier.eq.228) then 
         write (*,228) char, realv, int, char
      else
         write (*,999) ier, char, realv, int
      end if

1     format (/,'**warning ver001** the amount of a saturated phase is'
     *       ,' < 0, this indicates that',/,'the specified amount of a '
     *       ,'saturated component is inadequate to saturate the system'
     *       ,/)
2     format (/,'**warning ver002** the molar amount of a phase is ',
     *        g12.3,
     *        ' (< -zero_mode) this may be',/,'indicative of numeric ',
     *        'instability',/)
3     format (/,'**warning ver003** the solution model file is ',
     *         ' inconsistent with this',/,
     *         'this version of Perple_X. Update the file and/or ',
     *         'Perple_X',/)
4     format (/,'**warning ver004** the data includes ',a,' values, '
     *      ,'probably because bad_number',/,'in perplex_option.dat = '
     *      ,'NaN, these values will be replaced by zeros. To avoid ',/
     *      ,'this behavior set bad_number to a numeric value or use a'
     *      ,' plotting program capable',/,'of handling NaNs, e.g., ',
     *       'MatLab or PYWERAMI',/)
5     format (/,'**warning ver005** fluid components are specified',
     *        ' as thermodynamic AND as either',/,'saturated phase',   
     *      ' or saturated components; almost certainly a BAD idea.',/)
6     format (/,'**warning ver006** fluid components are specified',
     *        ' as both thermodynamic AND',/,'saturated ',    
     *        'components; almost certainly a BAD idea.',/)
7     format (/,'**warning ver007** fluid components are specified as'
     *       ,' a saturated phase component',/,'AND as a thermodynamic', 
     *        'or saturated component; almost certainly a BAD idea.',/)
8     format ('**warning ver08** too exclude more phases ',
     *        'increase paramter h8 (',i3,')')
9     format ('**warning ver009** unable to deconstruct transition,'
     *       ,/,'data for ',a,' will not be output.')
10    format (/,'**warning ver010** not able to traverse ',
     *          'the entire  extent of equilibrium',/,'(',i6,')',
     *          ' at v(3)=',g12.6,/,
     *          'this error can usually be avoided by increasing the ',
     *          'finite difference',/,'increment delt(iv(1)) or delv',
     *          '(iv(2)), as defined on card 6 of',/,'the file on n2.',
     *          ' In routine:',a,/)
11    format (/,'**warning ver011** ',a,' has > 1 transition with dp/d',
     *          'T ne 0 and may not be',/,'treated correctly.')
12    format (/,'**warning ver012** ',a,' has a transition ',
     *          ' with dp/dT < 0 and may not be treated ',/,
     *          ' correctly.')
13    format (/,'**warning ver013** because the total amount of the com'
     *         ,'mponents in ',a,'is <= 0',/,'it will be rejected from '
     *         ,'this calculation although it is a legitimate phase.',/,
     *          'To prevent this rejection transform the data base comp'
     *          'onents (e.g., using CTRANSF)',/,'so that the total amo'
     *         ,'unt of the components in ',a,' is > 0.',/)
14    format (/,'**warning ver014** You can not redefine the ',
     *          'saturated phase component:',a,/,'To circumvent this ',
     *          'restriction use CTRANSF to make a data base with the',/
     *         ,'the desired component transformations',/)
15    format (/,'**warning ver015** if you select > 1 saturated ',
     *          'component, then the order you',/,'enter the ',
     *          'components determines the saturation heirarchy and may'
     *          ,' effect your',/,'results (see Connolly 1990).',/)
16    format (/,'**warning ver016** ',a,' has been rejected because it',
     *       ' has no associated volumetric EoS.',/,'To override this ',
     *        'behavior set auto_exclude to false or add an ',
     *        'association.',/)
17    format (/,'**warning ver017** ',a,' is a relict equipartition so',
     *        'lution model. The use of',/,'such models in Perple_X 6.',
     *        '9.1+ may result in erratic, though formally correct,',/,
     *        'phase field geometries, to avoid this problem replace ',
     *        a,' with an up-to-date model.',/)
18    format (/,'**warning ver018** the value of the default dependen',
     *         't variable (',g14.6,') for the following',/,
     *         'equilibrium was inconsistent with the an earlier ',
     *         'determination of the invariant condition',/,
     *         'and will be reset. This may cause the curve to ',
     *         'kink near the invariant point',/)
19    format ('**warning ver019** you must specify at least ',
     *        'one thermodynamic component, try again',/)
20    format ('**warning ver020** sfol2')
21    format ('**warning ver021** xmax (',g12.6,') > 1 for '
     *         ,' solution model ',a,/,' xmax will be reset to 1',
     *        /,' see documentation, section 4.0.')
22    format ('**warning ver022** xmin (',g12.6,') < 0 for '
     *         ,' solution model ',a,/,' xmin will be reset to 1',
     *        /,' see documentation, section 4.0.')
23    format ('**warning ver023** xmin > xmax for solution ',a,/,
     *        'xmin will be set to xmax NO PSEUDOCOMPOUNDS WILL BE',
     *        ' GENERATED.',/,'see documentation, section 4.0',/)
24    format (/,'**warning ver024** wway, increment refined out of',
     *          ' range (',g8.1,')',/,'before the stable',
     *          ' extension of the equilibria was located')
25    format (/,'**warning ver025** ',i1,' endmembers for ',a,
     *          ' The solution will not be considered.',/)
26    format ('**warning ver026** only one endmember for ',a,
     *          ' The solution will not be considered.')
27    format (/,'**warning ver027** only ',i2,' user defined ',
     *      'compositions permitted.',/,'do multiple runs with WERAMI',
     *      'or redimension common block comps.',/)
28    format (/,'**warning ver028** minimization failed, ill-',
     *        'conditioned?',/)
29    format ('**warning ver029** programming error, routine ',a,/)
30    format (/,'**warning ver030** Because of missing endmembers, ',
     *        'or that the',/,'subdivision',
     *        ' scheme specified for solution model ',a,/,'is too',
     *        ' restrictive, there are no valid compositions for', 
     *        ' this model.',/)
31    format (/,'**warning ver031** this choice is disabled because ',
     *        'the dependent_potentials',/,'keyword is missing or off',
     *        ' in perplex_option.dat, to use this choice set the',/,
     *        'keyword to on and re-run VERTEX.',/)
32    format ('**warning ver032** fixed activity option requested',
     *          ' for ',a,/,'This option is disabled, the',
     *          ' solution will not be considered.')
33    format ('**warning ver033** missing endmembers for ',a,/,
     *        'The model may be recast in > 1 way for',
     *        ' the endmember subset.',/,'To control this choice',
     *        ' eliminate undesired endmembers.')
34    format ('**warning ver034** ',a,' could not be recast as',
     *          ' a simpler model.',/,'The solution will not be',
     *          ' considered. Add the missing endmembers or eliminate'
     *          ,/,'additional endmembers to allow this model.',/)
35    format (/,'**warning ver035** ',a,' is only for pure fluids',
     *        /,' XCO2 will be reset to: ',f4.2,/)
36    format ('**warning ver021** xinc (',g12.6,') < 0 for'
     *         ,' solution model ',a,/,'xinc will be reset to 1.'
     *         ,' see documentation, section 4.0. ',/)
37    format (/,'**warning ver37** you will not be able to plot the ',
     *       'results of this',/,'calculation with PSVDRAW. PSVDRAW ',
     *       'only plots ternary composition diagrams.',/)
38    format (/,'**warning ver38** you will not be able to plot the ',
     *       'results of this',/,'calculation with PSVDRAW. PSVDRAW ',
     *       'only plots mixed-variable diagrams for',/,'a binary ',
     *       'system with one independent potential variable.',/)
39    format (/,'**warning ver39** PSVDRAW will plot the results of ',
     *       'this calculation as a',/,'projected section, such plots ',
     *       'may be difficult to interpret. To plot',/,
     *       'pseudosections as a an explicit function of a systems ',
     *       'composition use',/, 'gridded minimization.',/)
40    format (/,'**warning ver040** you have configured a ',
     *       'problem with only one independent variable.',/)
41    format (/,'**warning ver040** ',a,' occurs ',i1,' times in an as',
     *       'semblage, properties will be',/,'reported for its last i',
     *       'nstance. To avoid this problem do not use multi-property',
     *       /,'choices to extract the properties of this solution.',/)
42    format (/,'**warning ver042** an optimization failed due ',
     *          'to numerical instability',/,
     *          'or because the phases of the system do not span ',
     *          'its bulk composition.',//,
     *          4x,'In the 1st case (best solutions listed first):',/,
     *          8x,'set intermediate_savrpc and intermediate_savrpc to',
     *             ' T and/or',/,
     *          8x,'increase (sic) optimization_precision and/or',/,
     *          8x,'increase (sic) replicate_threshold and/or',/,
     *          8x,'increase (sic) rep_dynamic_threshold.'/,
     *          4x,'see: www.perplex.ch/perplex_options.html for ',
     *          'explanation.',//,
     *          4x,'In the 2nd case: ',
     *          'change the bulk composition or add phases.',/)
43    format (/,'**warning ver043** ',a,' is the base 10 log of the',
     *       ' activity, values > 0 imply',/
     *       'supersaturation with respect to the reference species.',/,
     *       'Specify a different value (Y/N)?')
44    format (/,'**warning ver044** ordinarily ',a,' should be in the ',
     *        'range [0,1].'/,'Specify a different value (Y/N)?')
45    format (/,'**warning ver045** ',a,' involves a nonlinear EoS.',/,
     *          ' Output properties are stoichiometric averages.',/)
46    format (/,'**warning ver046** temperature (',g12.6,' K) is out',
     *        ' of range for endmember/phase: ',a,/,a,
     *        ' will be destabilized at this condition. In some cases',
     *        ' this problem can be corrected by',/,'setting ',
     *        'Anderson_Gruneisen to T in the Perple_X option file.',/)
47    format (/,'**warning ver047** univariant field ',i6,' terminates',
     *        ' at an invariant field',/,'that could not be located ',
     *         'within the tolerance specified in the thermodynamic',/,
     *         'data file (PTOL= ',g12.6,').',/)
48    format (/,'**warning ver048** fluid phase pseudocompound data ',
     *         'does not include',/,' volumetric properties (SWASH).',/)
49    format (/,'**warning ver049** warning ',i3,' will not be repeated'
     *         ,' for future instances of this problem.',/,
     *          'currently in routine: ',a,//,
     *          'To override this limit on the number of warnings set ',
     *          'warn_no_limit to T',/)
50    format (/,'**warning ver050** reformulating prismatic ',
     *          'solution: ',a,' because of missing endmembers. ',
     *        /,'(reformulation can be controlled explicitly ',
     *          'by excluding additional endmembers).',/)
51    format (/,'**warning ver051** cannot make ',a,' because of ',
     *          'missing data or an'
     *       ,/,'invalid definition in the thermodynamic data file.',/)
52    format (/,'**warning ver052** rejecting ',a,'; excluded or '
     *       ,'invalid composition.',/)
53    format (/,'**warning ver053** the failure rate during speciation',
     *          ' (o/d) calculations is ',f5.1,'%.',/,
     *          'A high failure rate may cause failed optimizations. ',
     *          'Usually the failure rate can be',/,'reduced by ',
     *          'increasing speciation_max_it in ',
     *          'perplex_option.dat',/)
54    format (/,'**warning ver054** property choices 25, 36, and 38 are'
     *         ,' not allowed in combination',/,'with other property '
     *         ,'choices',/)
55    format (/,'**warning ver055** possible compositions of ',a,' lie',
     *        ' entirely within the saturated',/,'component composition'
     *         ,' space. the compositions will not be considered.',/,
     *         'If this is problematic, then eliminate the component ',
     *         'saturation constraints',/,'or use convex.',/)
56    format (/,'**warning ver056** the EoS specified for ',a,' by the',
     *        ' hybrid_EoS option will be',/,'overridden by the EoS sp',
     *        'ecified in the problem definition file. To prevent this',
     *      /,'behavior set the GFSM option to True.',/)
57    format (/,'**warning ver044** ordinarily ',a,' should > 0, value',
     *          's <= 0 may cause numerical',/,
     *          'instability. Specify a different value (Y/N)?')
589   format (/,'**warning ver589** wway, the equilibrium of the '
     *         ,'following reaction',/,'is inconsistent with the ',
     *          'invariant equilibrium.',/)
59    format (/,'**warning ver059** endmember ',a,
     *        ' has invalid site populations.',/)
60    format (/,'**warning ver060** Tait EoS conditions out of range ',
     *          'for endmember/phase: ',a,'(at T=',g14.6,' K)')
601   format ('base v1 is negative, endmember/phase ',a,' will be ',
     *        'destabilized.',/)
602   format ('base v2 will be zeroed on the assumption that it is ',
     *        'insignificant.',/)
61    format (/,'**warning ver061** the data includes NaN values, '
     *      ,'probably because bad_number',/,'in perplex_option.dat = '
     *      ,'NaN, these values will be replaced by zeros. To avoid ',/,
     *       'this behavior set bad_number to a numeric value or use a',
     *       ' plotting program capable',/,'of handling NaNs, e.g., ',
     *       'MatLab or PYWERAMI.',//,'program/routine: ',a,/)
62    format (/,'**warning ver062** ',a,' is an electrolytic fluid, th',
     *        'e default value of ',/,'aq_bad_results has been changed',
     *        ' from err to 101 to allow fractionation to completely',/,
     *        'deplete solute components from the condensed phase asse',
     *        'mblage',/)
63    format (/,'**warning ver063** wway, invariant point on an edge?',
     *        /)
64    format (/,'**warning ver064** AQSOLV failed to converge on ionic',
     *       ' stength, currently =',f7.1,', this',/,
     *       'usually occurs at conditions where the HKF assumptions a',
     *       're invalid. Increasing speciation_max_it',/,
     *       'may facilitate convergence.',/)
68    format (/,'**warning ver068** degenerate initial assemblage in ',
     *          'COFACE, this should never occur',/,'if you see this ',
     *          'message please report the problem',/)
72    format (/,'**warning ver072** ',a,/)
73    format (/,'**warning ver073** an invariant point has been ',
     *          'skipped, routine: ',a,/,
     *          'decreasing DTOL (',g9.3,') in the thermodynamic ', 
     *          'data file for variable ',i1,/,
     *          'may eliminate this problem',/)
74    format (/,'**warning ver074** no new equilibria identified,',
     *          ' if degenerate segments have',/,' been skipped',
     *          ' increase the computational reliability level.',/)
79    format (/,'**warning ver079** univeq failed on an edge for ',
     *          'the following equilibrium.',/,' Probable cause is ',
     *          'extreme independent variable limits (e.g., xco2=0)',/
     *          ' or poor convergence criteria ',
     *          'in the thermodynamic data file. In routine:',a,/)
87    format (/,'**warning ver087** wway-univeq did not converge ',
     *          'when div was refined',/)
88    format (/,'**warning ver088** SMPLX converged to a non-unique ',
     *        'solution.',/,3x,'Probable cause: system composition ',
     *        'coincides with that of ',//,3x,'a compound or a ',
     *        'tieline between compounds.',//,3x,'This may lead to ',
     *        'erratic results.',//,3x,'To avoid this problem ',
     *        'perturb the systems composition.',/)
89    format (//,'**warning ver089** BUILD you did not request',
     *        'plot file output.',/,' You will not be able to process',
     *        ' the results of the requested calculation.',//)
90    format (/,'**warning ver090** an optimization failed. This may i',
     *        'ndicate an infeasible bulk',/,'composition or that the ',
     *        'LP_max_iteration option value is too small.',//,
     *        'In the former case optimization will fail for the bulk ',
     *        'composition at all physical',/,'conditions and the prob',
     *        'lem can only be remedied by increasing the range of',/,
     *        'compositions spanned by the possible phases of the ',
     *        'system.',/)
91    format (/,'**warning ver091** optimization failed. Change ',
     *        'minimnization method',/)
92    format (/,'**warning ver092** you have requested ',i4,
     *        ' grid points. Current',/,'dimensioning is for ',
     *        i4,' points. To obtain the requested resolution',/,
     *        'increase parameter L7 and recompile; or reduce the ',
     *        'required resolution via',/,'the ',a,' keyword in ',
     *        'perplex_option.dat',/,'The program will continue ',
     *        'with an effective grid resolution of ',i4,
     *        ' points.')
99    format (/,'**warning ver099** ',a,/)
100   format (/,'**warning ver',i3,'** ',a,/)
106   format ('**warning ver106** programming error in ',a)
108   format (/,'**warning ver108** wway, a phase field with the '
     *         ,'following',/,' reaction is stable on both ',
     *          'sides of an invariant point',/,' this error can ',
     *          'usually be avoided by increasing the finite ',/,
     *          ' difference increment delt(iv(1)) or delv',
     *          '(iv(2)), defined',
     *           /,' on card 6 of the thermodynamic data file',/)
109   format (/,'**warning ver109** you may ',
     *        'have assigned a mobile component as an independent ',/,
     *        ' variable without defining the component',/)
113   format (/,'**warning ver113** maximum variance for equilibrium',
     *        ' tracing [the variance keyword in ',/,
     *        'perplex_option.dat] must be > 0, but is ',i2,
     *        '. Set to 1 for the current calculation',/)
114   format (/,'**warning ver114** the default increment of an ',
     *        'independent variable is <',/,'1 percent of ',
     *        'its range, this is usually inefficient',/)
172   format (/,'**warning ver172** you cannot use this equation of',
     *          ' state with Y(CO2)',/,' as an indepedent variable, ',
     *          ' pick another one:',/)
173   format (/,'**warning ver173** invalid buffer choice ',/)
175   format (/,'**warning ver175** speciation routine ',a,' did',
     *          ' not converge ',/,' possibly due to graphite super-',
     *          'saturation. ier = ',i1,' real = ',g16.5,/)
176   format (/,'**warning ver176** fluid equation of state routine ',
     *        a,' did not converge.',/,'If execution continues this ',
     *        'may lead to incorrect results. To avoid this ',/,
     *        'problem increase speciation_max_it (',i4,') in ',
     *        'perplex_option.dat or choose',/,
     *        'a different equation of state.',//,
     *        'NOTE: at compositional extremes fluid speciation ',
     *        'calculations may require ',/,
     *        'thousands of iterations',/)
177   format (/,'**warning ver177** Invalid fluid speciation. ',
     *          'Reducing speciation tolerance (',g14.6,') in ',
     *          'perplex_option.dat',/,'may resolve this problem',/)
205   format (/,'**warning ver205** too many new phase assemblages, ',
     *        'found by routine newhld',/,'increase dimension j9 (',
     *        i8,')',/)
228   format (/,'**warning ver228** in solution model ',a,' negative ',
     *          'composition (',g12.6,') for component ',i2,/,
     *          'this may indicate an incorrect stoichiometric ',
     *          'dependent endmember definition.',//,
     *          'This warning is issued only for the first negative ',
     *          'composition of ',a,/)
900   format ('the calculation may be incomplete !!!!',/)
999   format (/,'**warning unspecified** ier =',i3,' routine ',a6
     *         ,' r = ',g12.6,' int = ',i9,/)
      end

      subroutine rmakes (jopt)
c----------------------------------------------------------------------
c rmakes is called by topn2 to read make definitions of thermodynamic
c entities, these entities are defined as a linear combination of 
c exisiting entities as defined in the thermodynamic file, with an 
c optional pressure/temperature dependent DQF correction. The format
c assumes data on one line of less than 240 characters, the expected format
c is

c name = num_1 * name_1 + num_2 * name_2 ....+ num_int * name_int
c dnum_1 dnum_2 dnum_3

c where i num_j is a number or fraction (i.e., two numbers separated by a 
c '/') and name_j is the name of the int existing entities. 
c and the dqf correction to the entity 'name' is
c Gdqf(J/mol) = dnum_1 + T[K]*dnum_2 + P[bar]*dnum_3

c end_of_data is either a "|" or the end of the record.

c make definitions are preceded by the keyword:

c begin_makes 

c and truncated by the keyword:

c end_makes

c if jopt > 3, data is echoed to LUN n8 (for ctransf/actcor).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, ier, iscan, i, nreact, jopt

      external iscan

      double precision rnum, nums(m3)

      character tname*8, name*8, rec*(lchar), tag*3

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      call readcd (n2,ier,.true.)
      if (ier.ne.0) goto 90 
c                                 echo data for ctransf/actcor
      if (jopt.gt.3) write (n8,'(400a)') chars(1:com)

      nmak = 0 

      write (rec,'(400a)') chars
      read (rec,'(a3)') tag

      do while (tag.ne.'end')   

         nmak = nmak + 1
         if (nmak.gt.k16) call error (55,mkcoef(1,1),nmak,'RMAKES')
c                                 get first name
         ibeg = 1
         call readnm (ibeg,iend,com,ier,tname)
         if (ier.ne.0) goto 90
c                                 find start of data marker '='
         ibeg = iscan (1,com,'=') + 1
c                                 the rest of the data should
c                                 consist of coefficients followed
c                                 by names
         nreact = 0 

         do while (ibeg.lt.com) 
c                                 find the number
            call readfr (rnum,ibeg,iend,com,ier)
            if (ier.eq.2) then 
c                                 ier = 2 = a read error
               goto 90
            else if (ier.eq.1) then 
c                                 ier = 1, end-of-definition
               exit 
            end if 
c                                 find the name
            call readnm (ibeg,iend,com,ier,name)
            if (ier.ne.0) goto 90

            nreact = nreact + 1
            if (nreact.gt.k17) call error (56,mkcoef(1,1),nmak,'RMAKES')

            mkcoef(nmak,nreact) = rnum 
            mknam(nmak,nreact) = name
           
         end do

         if (nreact+1.gt.k17) call error (56,mkcoef(1,1),nmak,'RMAKES')
         mknam(nmak,nreact+1) = tname
         mknum(nmak) = nreact
c                                 now the dqf
         call readcd (n2,ier,.true.)
         if (ier.ne.0) goto 90
c                                 echo data for ctransf/actcor 
         if (jopt.gt.3) write (n8,'(400a)') chars(1:com)
c                                 read the DQF coefficients
         ibeg = 1
         call redlpt (nums,ibeg,iend,ier) 
         if (ier.ne.0) goto 90

         do i = 1, m3 
            mdqf(nmak,i) = nums(i)
         end do 
c                                 start next make definition
         call readcd (n2,ier,.true.)
         write (rec,'(400a)') chars
         read (rec,'(a3)') tag
c                                 echo data for ctransf/actcor
         if (jopt.gt.3) write (n8,'(400a)') chars(1:com)

c                                 reject excluded makes
         do i = 1, ixct
            if (tname.eq.exname(i)) then 
               nreact = nreact - 1
               exit 
            end if
         end do 

      end do 

      goto 99

90    write (*,1000) chars(1:com)
      stop
      
1000  format (/,'**error ver200** READMK bad make definition in the',
     *        ' thermodynamic data file',/,'currently reading: ',/
     *        ,400a)

99    end 

      subroutine readnm (ibeg,iend,siz,ier,name)
c----------------------------------------------------------------------
c readnm looks for the first word in a record chars, ibeg is the index
c of the 1st letter, iend is the index of the last letter.

c input 
c         ibeg - starting index for search
c         siz  - end index for search
c output
c         ibeg - starting index of word
c         iend - end index of word
c         ier  - error code
c         name - word
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, iscan, iscnlt, ier, siz

      external iscan, iscnlt

      character name*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      ier = 0 
c                                 find start of name
      ibeg = iscnlt (ibeg,siz,' ') 
c                                 find next blank
      iend = iscan (ibeg,siz,' ') - 1
c                                 initialize to be safe:
      name = ' '

      if (iend - ibeg.le.7) then

         write (name,'(20a)') chars(ibeg:iend)

      else 
c                                 can't be a valid name, save it
c                                 anyway in case it's a tag
         write (name,'(20a)') chars(ibeg:ibeg+7)
         ier = 4

      end if 

      ibeg = iend + 1

      end 

      subroutine readcd (nloc,ier,strip)
c----------------------------------------------------------------------
c readcd - read 240 column card image from unit 9, strip out unwanted
c characters if strip. ier = 1 no card found.
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      logical strip

      integer ier, iscan, ict, i, iscnlt, ibeg, nloc

      external iscan, iscnlt

      character card*(lchar)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      ier = 0 

      ibeg = 0
  
      com = 0 

      card = ' '

      do while (ibeg.ge.com) 

         read (nloc,'(a)',end=90) card

         if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1
c                                 '|' in first column
            if (com.eq.0) cycle
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')

         end if 

      end do 
c                                 there is a non-blank data character
      if (strip) then 

         ict = 1

         do i = 2, com
c                                 strip out '+' and '*' chars
            if (chars(i).eq.'+'.or.chars(i).eq.'*') chars(i) = ' '
c                                 eliminate blanks after '/' and '-'
c                                 and double blanks
            if ((chars(ict).eq.'/'.and.chars(i  ).ne.' ') .or. 
     *          (chars(ict).eq.'-'.and.chars(i  ).ne.' ') .or.
     *          (chars(ict).eq.' '.and.chars(i  ).ne.' ') .or.
     *          (chars(ict).ne.'-'.and.chars(ict).ne.'/'.and.
     *           chars(ict).ne.' ') ) then
                ict = ict + 1
                chars(ict) = chars(i)
            end if

         end do 

         com = ict

      else
c                                 scan backwards to the last non-blank
         com = iscnlt (com,1,' ')

      end if

      goto 99

90    ier = 3

99    end

      subroutine readfr (rnum,ibeg,iend,siz,ier)
c----------------------------------------------------------------------
c readfr looks for a number or two numbers separated by a backslash / in
c array elements chars(iend:ibeg), the latter case is interpreted as a ratio. 
c the result is returned as num
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision rnum, rnum1 

      integer ibeg, iend, iback, ier, iscan, iscnlt, siz

      external iscan, iscnlt

      character num*30

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ier = 0 
c                                 now find start of a number
      ibeg = iscnlt (ibeg,siz,' ')  
c                                 find backslash
      iback = iscan (ibeg,siz,'/') - 1
c                                 find next blank
      iend = iscan (ibeg,siz,' ') - 1
c                                 three cases:
      if (iend.ge.com) then

         ier = 1
         goto 99 

      else if (iback.gt.iend) then
c                                 no fraction
         if (iend-ibeg+1.gt.30) goto 90
c                                 first constant
         write (num,'(30a)') chars(ibeg:iend)
         read (num,*,err=90) rnum

      else 
c                                 fraction write numerator
         if (iback+1-ibeg.gt.30) goto 90
c                                 first number
         write (num,'(30a)') chars(ibeg:iback)       
         read (num,*,err=90) rnum
c                                 second number 

         if (iend-iback-1.gt.30) goto 90
         write (num,'(30a)') chars(iback+2:iend)
         read (num,*,err=90) rnum1

         rnum = rnum/rnum1

      end if 

      ibeg = iend + 1

      goto 99

90    ier = 2

99    end 

      subroutine redfr0 (rnum,ibeg,iend,ier)
c----------------------------------------------------------------------
c redfr0 looks for a number or two numbers separated by a backslash / in
c that array chars(iend:ibeg), the latter case is interpreted as a ratio. 
c the result is returned as rnum. differs from readfr in that redfr0
c expects ibeg/iend are known on input.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision rnum, rnum1 

      integer ibeg, iend, iback, ier, iscan

      external iscan

      character num*30

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ier = 0 

c                                 find backslash
      iback = iscan (ibeg,iend,'/') - 1
c                                 two cases:
      if (iback.ge.iend) then
       
         iback = iscan(ibeg,iend,' ') - 1
c                                 no fraction
         if (iback-ibeg+1.gt.30) goto 90
c                                 simple number
         write (num,'(30a)') chars(ibeg:iback)
         read (num,*,err=90) rnum

      else 
c                                 fraction write numerator
         if (iback+1-ibeg.gt.30) goto 90
c                                 first number
         write (num,'(30a)') chars(ibeg:iback)       
         read (num,*,err=90) rnum
c                                 second number 
         if (iend-iback-1.gt.30) goto 90
         write (num,'(30a)') chars(iback+2:iend)      
         read (num,*,err=90) rnum1

         rnum = rnum/rnum1

      end if 

      return

90    ier = 2

      end

      subroutine getnam (name,ids)
c----------------------------------------------------------------------
c subroutine to retrieve phase name corresponding to index ids
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids

      character name*14

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c-----------------------------------------------------------------------
      if (ids.lt.0) then
c                                 simple compound:
         name = names(-ids)

      else  
c                                 solution phases:
         if (iopt(24).eq.0.or.lname(ids).eq.'unclassified') then
c                                 use model name
            name = fname(ids)

         else if (iopt(24).eq.1) then
c                                 use phase abbreviation
            name = aname(ids) 

         else 
c                                 use full name
            name = lname(ids)

         end if 

      end if 
      
      end 

      subroutine sopen 
c-----------------------------------------------------------------------
c simple file open for data echo programs (e.g., actcor, rewrit, ctransf)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character*100 n2name

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
c                                 first the thermo data file
      do 

         call fopen2 (2,n2name)
 
         if (iam.eq.6) then 
            write (*,1070) 'ctransf.dat'
            open (n8,file='ctransf.dat')
         else if (iam.eq.9) then 
            write (*,1070) 'actcor.dat'
            open (n8,file='actcor.dat')
         else if (iam.eq.10) then 
            write (*,1070) 'new_'//n2name
            open (n8,file='new_'//n2name)
         end if 

         exit

      end do 
 
1070  format (/,'Output will be written to file: ',a,/)
 
      end

      block data
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision ptx
      integer ipt2
      common/ cst32 /ptx(l5),ipt2

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)

      character*80 com
      common/delet/com 

      integer hs2p
      double precision hsb
      common/ cst84 /hsb(i8,4),hs2p(6)

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)

      double precision vrt
      integer irt
      logical sroot
      common/ rkroot /vrt,irt,sroot

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer jd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,jd
c-----------------------------------------------------------------------
      data hs2p/4, 5, 18, 19, 20, 21/

      data iff,ipt2,goodc,badc,badinv/3*0,6*0d0,h9*0,h9*0/
c
      data mu, uf/ k8*0d0, 2*0d0/

      data r/8.3144126d0/

      data gflu, sroot/ 2*.false./

      data com/' '/
c                                 tags for thermo data i/o
      data strgs/'G0 ','S0 ','V0 ','c1 ','c2 ','c3 ','c4 ','c5 ','c6 ',
     *           'c7 ','b1 ','b2 ','b3 ','b4 ','b5 ','b6 ','b7 ','b8 ',
     *           'b9 ','b10','c8 ','c9 ','c10','c11',
     *           'Tc ','B  ','p  ','v  ','cs1','cs2','cs3','cs4'/
      data mstrg/'m0','m1','m2','k0','k1','k2'/
      data dstrg/'d1','d2','d3','d4','d5','d6','d7','d8','d9'/
      data tstrg/'t1 ','t2 ','t3 ','t4 ','t5 ','t6 ','t7 ','t8 ','t9 ',
     *           't10','t11','t12','t13','t14','t15'/
      data e16st/'G0 ','S0 ','V0 ','Cp0','w ','q ','a1 ','a2 ','a3 ',
     *           'a4 ','c1 ','c2 ','HOH'/
c     data estrg/'eG0','eS0','eV0','ec1','ec2','ec3','ec4','ec5','ec6',
c    *           'ec7','eb1','eb2','eb3','eb4','eb5','eb6','eb7','eb8'/
c                                 tags for interaction coefficients (Redlich-Kister polynomial)
      data wstrg/'w0 ','wT ','wP ','wP1','wP2','wP0'/
c                                 fluid eos species
      data specie /
     *      'H2O ','CO2 ','CO  ','CH4 ','H2  ','H2S ','O2  ',
     *      'SO2 ','COS ','N2  ','NH3 ','O   ','SiO ','SiO2',
     *      'Si  ','C2H6','HF  '/

      data times,btime,etime/90*0d0/
c                                 na1 must be initalized because the
c                                 array element caq(jd,na1) is used to test
c                                 for aq speciation on output.
      data na1/1/

      end

      subroutine getphi (name,aq,eof)
c----------------------------------------------------------------------
c read phase data from the thermodynamic data file from lun N2, assumes
c topn2 has read the header section of the file.

c on input aq is a flag which determines if solute species data is
c accepted (ieos = 15-16).
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i, it, j, ier

      double precision ct

      logical eof, aq

      character key*22, val*3, name*8,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
      eof = .false.

      do 

         call redcd1 (n2,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (ier.lt.0) then
 
            eof = .true.
            exit

         else if (ier.gt.0) then
                 
            call error (23,ct,i,name)

         end if 
c                                 name          
         read (key,'(a)',iostat=ier) name
         if (ier.ne.0) exit
c                                 on the off chance of a loose end keyword
         if (key.eq.'end') cycle
c                                 EoS
         read (nval2,*,iostat=ier) ieos
         if (ier.ne.0) exit    
c                                 look for comments
c        write (commnt,'(80a)') chars(com:com+79)
c                                 composition
         call formul (n2)
c                                 thermodynamic data
         call indata (n2)
c                                 do component transformation if
c                                 itrans is not zero
         if (itrans.gt.0) then
 
            do i = 1, itrans
               it = ictr(i)
               if (comp(it).ne.0d0.and.ctrans(it,i).ne.0d0) then
c                                 ct is how much of the new
c                                 component is in the phase.
                  ct =  comp(it) / ctrans(it,i)
 
                  do j = 1, icmpn
                     comp(j) = comp(j) - ct * ctrans(j,i)
                  end do 
 
                  comp(it) = ct
               end if 
            end do
         end if

         if (.not.aq.and.(ieos.eq.15.or.ieos.eq.16)) cycle

c                                 a data writing program, don't mess
c                                 with ieos
         if (iam.ne.6.and.iam.ne.9) then
c                                 standard form with no volumetric EoS, 
c                                 reset ieos internally:
             if (ieos.gt.0.and.ieos.lt.5.and.thermo(3,k10).eq.0d0)
     *          ieos = 0

         end if 
 
         exit 

      end do

      end

      subroutine indata (lun)
c----------------------------------------------------------------------
c called by getphi to decompile thermodynamic data cards read from lun.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lun, ier, iscan, iscnlt, i, j, ibeg, iend, ic2p(k4), jkind

      character key*22, values*80, strg*80

      double precision var

      logical ok

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      double precision emodu
      common/ cst318 /emodu(k15)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(k4),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer ic
      common/ cst42 /ic(k0)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer iam
      common/ cst4 /iam

      save ic2p
      data ic2p/31,32,22,1,2,3,4,5,6,7,12,13,14,15,16,17,18,19,20,21,8,
     *          9,10,11,23,24,25,26,27,28,29,30/
c-----------------------------------------------------------------------
c                                 initialize data
c                                 flag for t-dependent disorder
      idiso = 0
c                                 flag for mock-lambda transitions
      ilam = 0 
c                                 counter of mock-lambda transitions
      jlam = 0 
c                                 flag for bulk + shear moduli
      ikind = 0
c                                 flag for just shear
      jkind = 0
c                                 hsc conversion 
      hsc(k10) = .false.
c                                 standard thermo parameters
      do i = 1, k4
         thermo(i,k10) = 0d0
      end do 
c                                 shear modulus
      do i = 1, k15
         emodu(i) = 0d0
      end do
c                                 lamda transitions
      do j = 1, m6
         do i = 1, m7
            tm(i,j) = 0d0
         end do
      end do
c                                 t-dependent disorder
      do i = 1, m8
         td(i) = 0d0
      end do 

      do 
c                                 find a data card
         call redcd0 (lun,ier,key,values,strg)
         if (ier.ne.0) call error (23,tot,ier,strg) 

         ibeg = 1

         if (key.eq.'end') then 

            exit 

         else if (key.eq.'transition') then 

            ibeg = iscan (ibeg,com,'=') + 1
            ibeg = iscnlt (ibeg,com,' ')
            iend = iscan (ibeg+1,com,'=') + 1
c                                 write ilam data to values
            write (values,'(80a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg)
c                                 ilam as read is the counter, code
c                                 currently assumes the data is entered
c                                 sequentially, therefore this isn't necessary.
            read (values,*,iostat=ier) ilam
            if (ier.ne.0) call error (23,tot,ier,strg)
c                                 next get the type flag jlam.
            ibeg = iend
            iend = iscnlt (ibeg,com,'9')

            write (values,'(80a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg)
            read (values,*,iostat=ier) jlam
            if (ier.ne.0) call error (23,tot,ier,strg) 
c                                 position for next keyword
            ibeg = iend

         end if
c                                 read remaining keywords and values
c                                 from card
         do 

            key = ''
c                                 locate end of keyword
            if (ibeg.ge.com) exit 
            iend = iscan (ibeg,com,'=') - 1
            if (iend.ge.com) exit
c                                 write keyword
            write (key,'(22a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg) 
c                                 locate data
            ibeg = iscnlt (iend+2,com,' ')
            iend = iscan (ibeg,com,' ')
c                                 write data 
            write (values,'(80a)',iostat=ier) chars(ibeg:iend)
            if (ier.ne.0) call error (23,tot,ier,strg) 
c                                 shift pointer to next key
            ibeg = iscnlt(iend,com,' ')
c                                 assign data
            ok = .false.
c                                 =====================================
c                                 thermo data 
            if (ieos.eq.12.or.ieos.eq.14.or.ieos.eq.17) then
c                                 calphad format
               do i = 1, k4
                  if (key.eq.strgs(i)) then 
                     read (values,*,iostat=ier) thermo(ic2p(i),k10)
                     if (ier.ne.0) call error (23,tot,ier,key) 
                     ok = .true.
                     exit 
                  end if 
               end do

            else if (ieos.eq.16) then 
c                                 DEW/HKF aqueous data
               do i = 1, 13
                  if (key.eq.e16st(i)) then 
                     read (values,*,iostat=ier) thermo(i,k10)
                     if (ier.ne.0) call error (23,tot,ier,key) 
                     ok = .true.
                     exit 
                  end if 
               end do

            else 
c                                 generic thermo data 
               do i = 1, 21

                  if (key.eq.strgs(i)) then 

                     read (values,*,iostat=ier) thermo(i,k10)
                     if (ier.ne.0) call error (23,tot,ier,strg) 
                     ok = .true.
                     exit

                  else if (key.eq.'GH') then

                     read (values,*,iostat=ier) thermo(1,k10)
                     if (ier.ne.0) call error (23,tot,ier,strg)
                     hsc(k10) = .true.

                     if (hscon.and.iam.ne.5) then 
c                                 convert HSC G0 to SUP G0
                        do j = 1, icomp
                           thermo(1,k10) = thermo(1,k10) 
     *                                   + tr*comp(ic(j))*sel(j)
                        end do

                     else if (hscon) then 

                        do j = 1, icmpn
                           thermo(1,k10) = thermo(1,k10) 
     *                                   + tr*comp(j)*sel(j)
                        end do

                     end if 
 
                     ok = .true.
                     exit

                  end if 

               end do 

            end if 

            if (ok) cycle
c                                 =====================================
c                                 shear mod data 
            do i = 1, 6

               if (key.eq.mstrg(i)) then 
c                                 set shear/bulk mod flag
                  if (i.lt.4) then
c                                 shear
                     ikind = 1

                  else if (i.gt.3) then 
c                                 bulk
                     jkind = 1

                  end if 
                  
                  read (values,*,iostat=ier) emodu(i)
                  if (ier.ne.0) call error (23,tot,ier,strg) 
                  ok = .true.

                  exit

               end if

            end do

            if (ok) cycle
c                                 =====================================
c                                 explicit temperature-dependent disorder data 
            do i = 1, m8
               if (key.eq.dstrg(i)) then 
c                                 set disorder flag
                  idiso = 1
                  read (values,*,iostat=ier) td(i)
                  if (ier.ne.0) call error (23,tot,ier,strg) 
                  ok = .true.
                  exit 
               end if 
            end do 

            if (ok) cycle
c                                 =====================================
c                                 mock-lambda transition data 
            do i = 1, m7
               if (key.eq.tstrg(i)) then 
                  read (values,*,iostat=ier) tm(i,ilam)
                  if (ier.ne.0) call error (23,tot,ier,strg) 
                  ok = .true.
                  exit 
               end if 
            end do 

            if (ok) cycle

            call error (9,var,i,key)

         end do

      end do
c                                 reset ikind flag:
c                                 ikind = 0, no explicit moduli
c                                 ikind = 1 and jkind = 0 set ikind = 1 => just shear
c                                 ikind = 1 and jkind = 1 set ikind = 2 => both
c                                 ikind = 0 and jkind = 1 set ikind = 3 => just bulk
      ikind = ikind + jkind
      if (ikind.eq.1.and.jkind.eq.1) ikind = 3

      end

      subroutine getkey (n,ier,key,values,strg)
c----------------------------------------------------------------------
c getkey calls redcd0 and outputs error message 21 on error.
c----------------------------------------------------------------------
      implicit none
 
      integer n, ier 

      character key*(*), values*(*), strg*(*)
c----------------------------------------------------------------------

      call redcd0 (n,ier,key,values,strg)

      if (ier.ne.0) call error (21,0d0,n,strg)

      end 

      subroutine redcd0 (lun,ier,key,values,strg)
c----------------------------------------------------------------------
c this routine seeks a non-blank card that contains data, i.e., something
c other than a comment (text preceded by a "|") character. 
c the first word of the data is saved as key, the remaining words are
c are saved in values, and the complete data is saved as strg.
c the full record (including comments) is saved in chars.

c the card is also loaded into chars with:

c  length - position of last non-blank character
c  com    - position of the comment character
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer lun, ier, iscan, iscnlt, ibeg, iend, iblank

      character card*(lchar), key*(*), values*(*), strg*(*)

      external iscan, iscnlt

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      ier = 0 
      key = ' '

      do 

         read (lun,'(a)',iostat=ier) card

         if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,lchar,'|') - 1
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')
c                                 find the next blank
            iblank = iscan (ibeg,com,' ')
c                                 len < ibeg => only comments
            if (ibeg.ge.com) cycle
c                                 full record length
            length = iscnlt (lchar,1,' ')

            exit 

         else if (ier.ne.0) then 

            exit 

         end if 

      end do 

      if (ier.eq.0) then 
c                                 find end of keyword 
         iend = ibeg + 1
         iend = iscan (iend,lchar,' ') - 1
         if (iend.gt.22) iend = 22
c                                 load chars into key
         write (key,'(22a)') chars(ibeg:iend)
c                                 now the values
         ibeg = iscnlt (iend+1,lchar,' ') 

         if (ibeg.lt.com) then 

            iend = iscnlt (com,ibeg,' ')
            if (iend-ibeg.gt.79) iend = ibeg + 79
c                                 load chars into value
            write (values,'(80a)') chars(ibeg:iend)
c                                 load chars into strg
            if (iend.gt.80) iend = 80
            write (strg,'(80a)') chars(1:iend)
        
         else
c                                 no values
            strg = key

         end if 

      end if 

      end 

      subroutine formul (lun)
c----------------------------------------------------------------------
c formul reads a text formula and coverts the result to the composition
c array comp. allowed component names are read from xcmpnt array to 
c avoid problems associated with component transformations. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lun, len0, len1, ier, iscan, i, ibeg, iend

      character key*22, values*80, strg*80, ctemp*5

      logical ok

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)
c-----------------------------------------------------------------------
      do i = 1, icmpn
         comp(i) = 0d0
      end do

      call getkey (lun,ier,key,values,strg)

      if (ier.ne.0) call error (23,0d0,i,strg)

      ibeg = 1
      iend = iscan (ibeg,lchar,' ') - 1

      do 
c                                 find the "(" and ")"
         len0 = iscan (ibeg,iend,'(') 
         len1 = iscan (len0,iend,')')
c                                 write the name and number
         write (ctemp,'(5a)')   chars(ibeg:len0-1)
c                                 identify the component
         ok = .false.

         do i = 1, icmpn

            if (xcmpnt(i).eq.ctemp) then
               call redfr0 (comp(i),len0+1,len1-1,ier)
               if (ier.eq.0) ok = .true.
               exit 
            end if 

         end do      

         if (.not.ok) call error (23,0d0,i,strg)

         if (len1.eq.iend) exit

         ibeg = len1 + 1

      end do                  

      end 

      subroutine outdat (lun,id,option)
c----------------------------------------------------------------------
c lun    - output LUN
c id     - phase pointer
c option - source of compositional data 

c called by frendly, vertex, ctransf, actcor, rewrit. for processed data 
c requires unver and unlam to recover original data; unlam puts the 
c transition data into the local array tm; other data is the primary 
c arrays (cp or comp[see option below], thermo, therdi)

c if option = 0  then formula for entity id is created from the 
c    composition array comp and name array cmpnt.
c if option = 1 then formula for entity id created from the 
c    composition array cp and name array cmpnt via pointer ic.
c if option = 2  then formula for entity id created from the 
c    composition array cp0 and name array cmpnt.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lun, i, j, ibeg, iend, id, option, jcomp, siz

      character text(14)*1

      double precision var, dg

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      integer ic
      common/ cst42 /ic(k0)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer eos
      common/ cst303 /eos(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      character*80 commnt
      common/delet/commnt

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)
c-----------------------------------------------------------------------
c                                 =====================================
c                                 name & EoS
      write (lun,*) 
      read (names(id),'(8a)') chars(1:8)
      ibeg = 9
      var = eos(id)
      call outthr (var,' EoS',4,ibeg) 

      if (commnt.ne.' ') then 
         chars(ibeg) = '|'
         read (commnt,'(80a)') chars(ibeg+1:ibeg+80)
         ibeg = ibeg + 80
      end if 

      write (lun,'(400a)') chars(1:ibeg)
c                                 =====================================
c                                 formula
      ibeg = 1
      iend = 0

      if (option.eq.1) then 
         jcomp = icomp
      else 
          jcomp = icmpn
      end if 

      dg = 0d0 

      do i = 1, jcomp

         if (option.eq.0) then
            var = comp(i)
         else if (option.eq.1) then 
            var = cp(i,id)
         else if (option.eq.2) then 
            var = cp0(i,id)
         end if 

         if (var.ne.0) then 
c                                 load text name
            iend = ibeg + cl(ic(i)) - 1

            read (cmpnt(ic(i)),'(5a)') chars(ibeg:iend)
c                                 left parenthesis
            chars(iend + 1) = '('
c                                 get number
            call numtxt (var,text,siz)
c                                 load number into chars
            ibeg = iend + 2
            iend = ibeg + siz - 1

            do j = ibeg, iend
               chars(j) = text(j-ibeg+1)
            end do
c                                get the delta g HSC correction
            dg = dg + var*sel(i)

            chars(j) = ')'
 
            ibeg = iend + 2

         end if 

      end do 
c                                 write the formula
      write (lun,'(400a)') chars(1:iend+1)
c                                 =====================================
c                                 thermo data
      if (eos(id).eq.16) then 
c                                 HKF aqueous electrolyte data (13 values)
         ibeg = 1
 
         do i = 1, 5
            call outthr (thermo(i,id),e16st(i),3,ibeg)
         end do

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

         ibeg = 1
 
         do i = 6, 10
            call outthr (thermo(i,id),e16st(i),3,ibeg)
         end do

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

         ibeg = 1
 
         do i = 11, 13
            call outthr (thermo(i,id),e16st(i),3,ibeg)
         end do

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      else if (eos(id).eq.12.or.eos(id).eq.14.or.eos(id).eq.17) then 

         call error (72,r,id,'routine OUTDAT is not programmed to '//
     *              'output CALPHAD data tags.')

      else 

         ibeg = 1

         if (hscon.and.hsc(id)) then
c                                 convert back to HSC apparent G
            call outthr (thermo(1,id) - tr*dg,'GH',2,ibeg)

         else if (hsc(id)) then 
c                                 direct output of HSC apparent G
            call outthr (thermo(1,id),'GH',2,ibeg)

         else

            call outthr (thermo(1,id),strgs(1),2,ibeg)

         end if
 
         do i = 2, 3
            call outthr (thermo(i,id),strgs(i),2,ibeg)
         end do
c                                 write G,S,V
         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)
c                                 c1->c7 of thermo data
         ibeg = 1
  
         do i = 4, 10
            call outthr (thermo(i,id),strgs(i),2,ibeg)
         end do
c                                 write c1->c7
         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)
c                                 b1->b10 of thermo data
         ibeg = 1

         do i = 11, 18
            call outthr (thermo(i,id),strgs(i),2,ibeg)
         end do
c                                 write b1->b8, c8
         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      end if 
c                                 =====================================
c                                 shear/bulk modulus
      ibeg = 1

      do i = 1, 6
         call outthr (emod(i,id),mstrg(i),2,ibeg)
      end do

      if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)
c                                 =====================================
c                                 disorder parameters
      if (idis(id).ne.0) then

         ibeg = 1

         do i = 1, 8
            call outthr (td(i),dstrg(i),2,ibeg)
         end do 

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      end if 
c                                 =====================================
c                                 transition parameters
      do i = 1, lct(id)

         ibeg = 1
         var = i
         call outthr (var,'transition',10,ibeg)
         var = ltyp(id)
         call outthr (var,'type',4,ibeg) 

         do j = 1, m7
            call outthr (tm(j,i),tstrg(j),3,ibeg)
         end do 

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

      end do 

      write (lun,1000)

1000  format ('end',/)

      end

      subroutine outthr (num,strg,siz,ibeg)
c----------------------------------------------------------------------
c output prettified data.

c    num - the numeric data 
c    strg - a text tag for the data
c    siz  - length of strg
c    ibeg - pointer to the location for the data in the output array (chars)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision num

      character strg*(*), text(14)*1

      integer i, ibeg, iend, siz, len0, jend

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      if (num.ne.0d0.or.strg.eq.'EoS') then 
c                                 pad with 2 blanks, if not at line begining
         if (ibeg.gt.1) then
            chars(ibeg) = ' '
            ibeg = ibeg + 1
         end if

         iend = ibeg + siz - 1

         read (strg,'(14a)') chars(ibeg:iend)
c                                 trim out trailing blanks
         jend = ibeg

         do i = ibeg + 1, iend
            if (chars(i).eq.' ') cycle
            jend = jend + 1
         end do

         iend = jend

         chars(iend+1) = ' '
         chars(iend+2) = '='
         chars(iend+3) = ' '

         call numtxt (num,text,len0)

         do i = 1, len0
            chars(iend+3+i) = text(i)
         end do

         chars(iend+3+i) = ' '
         chars(iend+4+i) = ' '

         ibeg = iend + 4 + i

      end if 

      end 

      subroutine numtxt (num,text,siz)
c----------------------------------------------------------------------
c convert a g14.7e2 number to simplest possible text
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision num

      character text(*)*1, strg*14

      logical dec

      integer i, siz, inum, ier, ibeg, iend, jscnlt, jscan

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

c----------------------------------------------------------------------
      inum = int(num)

      siz = 14

      if (dabs(num-inum).lt.zero) then 
c                                 the number can be represented as 
c                                 an integer
         write (strg,'(i14)',iostat=ier) inum

      else 

         write (strg,'(g14.7E2)',iostat=ier) num

      end if

      read (strg,'(14a)') text(1:siz)

      ibeg = jscnlt (1,siz,' ',text)
      iend = jscan (ibeg,siz,' ',text) - 1
c                                 shift text left
      siz = 0 

      dec = .true.

      do i = ibeg, iend

         siz = siz + 1
         text(siz) = text(i)

         if (text(siz).gt.'A') dec = .false. 

      end do 
c                                 pruning:
      if (text(1).eq.'0') then
c                                 cut leading zero/+
         do i = 1, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1 

      else if (text(1).eq.'-'.and.text(2).eq.'0') then
c                                 cut leading zero
         do i = 2, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1

      end if

      if (dec) then 
c                                decimal number
         iend = jscan (1,siz,'.',text)
c                                reduce len to cut trailing zeroes
         if (iend.lt.siz) siz = jscnlt (siz,iend,'0',text)

      else if (num-inum.ne.0d0) then 
c                                 find the E char
         iend = jscnlt (1,siz,'A',text)
         ibeg = jscnlt (iend-1,1,'0',text) + 1
         inum = iend - ibeg
c             
         do i = ibeg, siz - inum
            text(i) = text(i + inum)
         end do   

         siz = siz - inum
c                                 the E character is now at
         ibeg = iend - inum   

         if (text(ibeg+1).eq.'+') then

            inum = 1 
            if (text(ibeg+2).eq.'0') inum = 2
c                                 delete superfluous + and 0
            do i = ibeg+1, siz - inum
               text(i) = text(i + inum)
            end do
          
            siz = siz - inum

         else if (text(ibeg+1).eq.'-') then
c                                 delete superfluous 0
            if (text(ibeg+2).eq.'0') then
               do i = ibeg+2, siz - 1
                  text(i) = text(i + 1)
               end do
           
               siz = siz - 1

            end if 

         end if 
       
      end if 

      end 

      subroutine znmtxt (num,text,siz)
c----------------------------------------------------------------------
c convert a f7.3 number to simplest possible text, 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision num

      character text(*)*1, strg*7

      integer i, siz, inum, ier, ibeg, iend, jscnlt, jscan

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      inum = int(num)

      siz = 7

      if (dabs(num-inum).lt.zero) then 
c                                 the number can be represented as 
c                                 an integer
         write (strg,'(i7)',iostat=ier) inum

      else 

         write (strg,'(f7.4)',iostat=ier) num

      end if

      read (strg,'(7a)') text(1:siz)

      ibeg = jscnlt (1,siz,' ',text)
      iend = jscan (ibeg,siz,' ',text) - 1
c                                 shift text left
      siz = 0 

      do i = ibeg, iend

         siz = siz + 1
         text(siz) = text(i)

      end do 
c                                 pruning:
      if (text(1).eq.'0') then
c                                 cut leading zero/+
         do i = 1, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1 

      else if (text(1).eq.'-'.and.text(2).eq.'0') then
c                                 cut leading zero
         do i = 2, siz - 1
            text(i) = text(i + 1)
         end do

         siz = siz - 1

      end if

      do i = siz + 1, 7
         text(i) = ' '
      end do

      iend = jscan (1,siz,'.',text)
c                                reduce len to cut trailing zeroes
c     if (iend.lt.siz) siz = jscnlt (siz,iend,'0',text)

      end

      subroutine fopen1 
c-----------------------------------------------------------------------
c fopen1 gets the project name and opens the problem definition file
c    n1name = project_name.dat
c iam is a flag indicating the calling program:
c    4 - build
c    1 - vertex
c    2 - meemum
c   13 - unsplt, global
c   14 - unplst, local
c------------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character y*1,n1name*100

      integer ierr

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iam
      common/ cst4 /iam

      integer jx, jy, lev, xn, yn
      common/ cst58 /jx, jy, lev, xn, yn
c-----------------------------------------------------------------------
      do 
c                                 get the root for all output files
c                                 except if unsplt-local
         if (iam.ne.14) then 

            if (iam.eq.4) then 
c                                 BUILD
               write (*,1040)
c                                 readrt loads the root into prject
               call readrt

            else  
c                                 VERTEX, MEEMUM, and plotting programs
               write (*,1030)
c                                 Amir #1
               call readrt

            end if 

         end if 
c                                 make the problem definition file name
         call mertxt (n1name,prject,'.dat',0)

         if (iam.eq.4) then 

            write (*,1070) n1name
c                                 BUILD
            open (n1,file=n1name,iostat=ierr,status='new')

            if (ierr.ne.0) then
c                                 name exists
               write (*,1050) n1name
               read (*,'(a)') y

               if (y.eq.'Y'.or.y.eq.'y') then 
c                                 overwrite it
                  open (n1,file=n1name)

               else
c                                 try again 
                  cycle 

               end if 

            end if
         
         else
c                                 Amir #2
c           prject = 'amir_mantle_input'
c                                 VERTEX, MEEMUM, UNSPLT
            open (n1,file=n1name,iostat=ierr,status='old')

            if (ierr.ne.0) then
c                                 name does not exist
               write (*,1080) n1name
               read (*,'(a)') y
               if (y.eq.'Y'.or.y.eq.'y') then 
c                                 try again
                  cycle 

               else 
c                                 quit
                  stop 
 
               end if 

            end if

            if (iam.eq.13) then
c                                 unsplt, read my_project.spt
               call mertxt (tfname,prject,'.spt',0)
               open (n8,file=tfname,iostat=ierr,status='old')
               if (ierr.ne.0) then
c                                 file does not exist
                  call error (12,0d0,ierr,tfname)
               end if 

               read (n8,*,iostat=ierr) jx
               if (ierr.ne.0) call error (12,0d0,ierr,tfname)
               read (n8,*,iostat=ierr) jy
               if (ierr.ne.0) call error (12,0d0,ierr,tfname)

            end if
         end if

         exit 

      end do 

1030  format (/,'Enter the project name (the name assigned ',
     *        'in BUILD) [default = my_project]:')
1040  format (/,'Enter a name for this project (the name',
     *        ' will be used as the',/,'root for all output file names)'
     *       ,' [default = my_project]:')
1050  format (/,'The file: ',a,/,'exists, overwrite it (y/n)?')
1070  format (/,'The problem definition file will be named: ',a) 
1080  format (/,'**warning ver191** no problem definition file named: ',
     *       a,/,'Run BUILD to create the file or change project names.'
     *       ,//,'Enter a different project name (y/n)?')

      end 

      subroutine readrt
c----------------------------------------------------------------------
c readrt - read file name root from terminal
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'
 
      integer kscan, iscnlt, ierr, siz

      character*100 prject, tfname
      common/ cst228 /prject,tfname

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      do 

         read (*,'(a)') prject

         if (prject.ne.' ') then 

            read (prject,'(100a)') chars(1:100)
c                                 find end of name ' '
            length = iscnlt (100,1,' ') 
c                                 check length
            if (length.gt.90) then 
               write (*,1010) 
               cycle 
            end if 
c                                 look for path characters / or \
            siz = kscan (100,1,'/')
            if (siz.eq.0) siz = kscan (100,1,'\')

            if (siz.eq.length) then 
               write (*,1030)
               cycle
            end if 
c                                 check if directory is valid
            if (siz.ne.0) then

               write (tfname,'(100a)') chars(1:siz)
               call mertxt (tfname,tfname,'delete_me',0)

               open (n1,file=tfname,iostat = ierr)
               close (n1,status='delete')

               if (ierr.ne.0) then 
                  write (*,1040)
                  cycle 
               end if  
c                                 mertxt uses chars, so re-read chars
               read (prject,'(100a)') chars(1:100)

            end if
c                                 look for illegal "." character
            if (kscan(siz+1,length,'.').lt.length) then 
               write (*,1000)
               cycle 
            end if 
c                                 look for illegal " " character
            if (kscan(siz+1,length,' ').lt.length) then 
               write (*,1020)
               cycle
            end if 

         else

            prject = 'my_project'

         end if 

         exit

      end do 

1000  format (/,'file/project names cannot include . characters, ',
     *          'try again',/)
1010  format (/,'file/project names must be < 91 characters, '
     *         ,'try again',/)
1020  format (/,'file/project names cannot include blanks, ',
     *          'try again',/)
1030  format (/,'file/project names cannot end with a / or \ character',
     *        ', try again',/)
1040  format (/,'the path specified in your project name is invalid,',
     *          ' check that all the ',/,
     *          'directories in the path exist, try again.',/)
      end

      subroutine getrt
c----------------------------------------------------------------------
c getrt - extracts root file name from the full file name in tfname
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'
 
      integer kscan, siz

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      read (tfname,'(100a)') chars(1:100)
c                                 find end of name ' '
      length = kscan (1,100,' ') - 1
c                                 look for dot character
      siz = kscan (length,1,'.') - 1

      if (siz.le.0) siz = length

      write (prject,'(100a)') chars(1:siz)

      end

      subroutine fopen2 (jam,name)
c-----------------------------------------------------------------------
c fopen2 - choose and open a thermodynamic data file on unit n2, jam 
c indicates behavior required by the calling program:
c  0 - name passed as argument, error if not found.
c  1 - BUILD, queries for name and writes it to N1
c  2 - queries for name
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*100 name, y*1, ddata*14, text*140

      integer ierr, jam

      data ddata/'hp02ver.dat   '/
c-----------------------------------------------------------------------

      do 

         if (jam.ne.0) then 
            write (*,1000)
            read (*,'(a)') name
            if (name.eq.' ') name = ddata
         end if 

         open (n2,file=name,iostat=ierr,status='old')

         if (ierr.ne.0) then
c                                 system could not find the file
            if (jam.eq.0) call error (120,0d0,n2,name)
c                                 if not vertex allow another try
            write (*,1010) name
            read (*,'(a)') y

            if (y.ne.'Y'.and.y.ne.'y') then
               write (*,1060)
               stop
            end if             
c                                 try again
            cycle

         end if
 
         if (jam.ne.1) exit 
c                                 BUILD, echo name to n1: 
         call mertxt (text,name,'thermodynamic data file',5)
         write (n1,'(a)') text 

         exit

      end do 
 
1000  format (/,'Enter thermodynamic data file name',
     *          ' [default = hp02ver.dat]:')
1010  format (/,'**warning ver191** FOPEN2 cannot find file:',/,a
     *         ,//,'try again (y/n)?')
1060  format (/,'O.K., I quit too.')

      end

      subroutine mertxt (text,text1,text2,nblank)
c----------------------------------------------------------------------- 
c mertxt - merge text two strings cutting trailing/leading blanks, with 
c          nblank characters between the two strings if text1 eq 
c          blank then pad by an initial 40 blanks.

c     input  - text1, text2 - character strings
c              nblank - number of blanks between the strings in text
c     output - text - character string 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer nchar1, nchar2, nblank

      character text*(*), text1*(*), text2*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      chars(1:lchar) = ' '
c                                 strip leading blanks in text1 and
c                                 get pointer to end of string
      call leblnk (text1,1,nchar1)

      if (nchar1.le.0) then 
c                                 text1 is blank
         nchar1 = 40

      else 
c                                 put nblank blanks between 1st and 
c                                 2nd strings, this is necessary despite
c                                 initialization because leblnk may 
c                                 shift strings left.
         chars(nchar1+1:nchar1 + nblank) = ' '

      end if
c                                 nchar1 points to the first empty char
      nchar1 = nchar1 + nblank + 1 
c                                 strip leading blanks from string in text2 and
c                                 get pointer to end of string in chars
      call leblnk (text2,nchar1,nchar2)

      text = ' '

      if (nchar2.gt.len(text)) call error (10,0d0,len(text),text2)

      write (text,'(400a)') chars(1:nchar2)

      end

      subroutine gettrn (jopt)
c----------------------------------------------------------------------
c jopt = 3 -> build
c jopt = 5 -> ctransf
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,jopt,ict,ier,jscan
 
      character*5 pname, rname, y*1

      double precision sum, ssum

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec
c-----------------------------------------------------------------------
c                                 recombine components:
      do 
         write (*,1030)
         write (*,1040) (cmpnt(i), i = 1, icmpn)
         write (*,1050)
         read (*,'(a)') y
         if (y.ne.'Y'.and.y.ne.'y') exit

         write (*,1060)
         read (*,'(a)') pname
         if (pname.eq.' ') exit
c                                 get the identity of the real comp
c                                 to be replaced.
50       write (*,1070) pname
         read (*,'(a)') rname

         do i = 1, icmpn
            if (rname.eq.cmpnt(i)) then
c                                 matches a name, check if it's a 
c                                 special component
               do j = 1, ispec

                  if (i.ne.idspe(j)) cycle 
c                                 matches special component id(j) 
                  if (jopt.eq.3) then
c                                 don't allow build users 
c                                 to transform saturated
c                                 phase components
                     call warn (14,atwt(1),i,cmpnt(i))
                     goto 60
                  else 
c                                 ctransf, ask the user if the 
c                                 new component will be a special 
c                                 component
                     write (*,1010) cmpnt(i),pname
                     read (*,'(a)') y
                     if (y.eq.'y'.or.y.eq.'Y') cycle
                     idspe(j) = 0 

                  end if 
               end do 

               icout(1) = i
               goto 70

            end if
         end do 
 
60       write (*,1080) 
         write (*,1040) (cmpnt(i), i = 1, icmpn)
         goto 50
c                                 get the identities of the other 
c                                 components in the new component:
70       itrans = itrans + 1
         ict = 1
         if (itrans.gt.k0) call error (999,atwt(1),ict,'GETTRN')
       
         write (*,4050) k5-1,pname
30       read (*,'(a)') rname
         if (rname.eq.'     ') goto 80

         do i = 1, icmpn
            if (rname.eq.cmpnt(i)) then 
               ict = ict + 1
               icout(ict) = i
               goto 30
            end if
         end do 
c                                 no match, try again message
         write (*,2300)
         goto 30
c                                 get the component stoichiometries:
80       write (*,4030) (cmpnt(icout(i)),i=1,ict)
         write (*,4040) pname

         do 
            read (*,*,iostat=ier) (ctrans(icout(i),itrans), i= 1, ict)
            if (ier.eq.0) exit
            call rerr
         end do 
 
         write (*,1100) pname,(ctrans(icout(i),itrans),
     *                      cmpnt(icout(i)), i = 1, ict)
         write (*,1110)
         read (*,'(a)') y
 
         if (y.eq.'y'.or.y.eq.'Y') then
            sum = 0d0
            ssum = 0d0             
            do i = 1, ict
               sum = sum + ctrans(icout(i),itrans) * atwt(icout(i))
               ssum = ssum + ctrans(icout(i),itrans) * sel(icout(i))
            end do 
            atwt(icout(1)) = sum
            sel(icout(1)) = ssum
            cmpnt(icout(1)) = pname
            cl(icout(1)) = jscan(1,5,' ',pname) - 1
            tcname(itrans) = pname
            ictr(itrans) = icout(1)
         else
            itrans = itrans - 1
            write (*,1000) 
         end if

      end do 

1000  format ('Try again.')
1010  format (/,a,' is a possible saturated phase component. Is ',
     *        'the new component ',a,/,'also a possible saturated ',
     *        'phase component (Y/N)?')
1030  format (/,'The current data base components are:')
1040  format (12(1x,a))
1050  format ('Transform them (Y/N)? ')
1060  format ('Enter new component name, < 6 characters,',
     *          ' left justified: ')
1070  format ('Enter old component to be replaced',
     *          ' with ',a,': ')
1080  format ('Select the component from the set: ')
1100  format (1x,a,' = ',6(f6.2,1x,a),/,9x,6(f6.2,1x,a))
1110  format ('Is this correct (Y/N)? ')
2300  format (/,'You made a mistake, try again.',/
     *          'Check spelling and upper/lower case matches.',/)
4030  format ('Enter stoichiometric coefficients of:',/,
     *        2x,12(a,1x))
4040  format ('in ',a,' (in above order): ')
4050  format ('Enter other components (< ',i2,') in ',a,' 1 per',
     *        ' line, <enter> to finish:')

      end

      subroutine topn2 (option)
c----------------------------------------------------------------------
c topn2 reads the header of the thermodynamic data file, if option > 3 
c echoes data to n8
c
c     option  calling program
c       0     vertex, meemum (all calls from clib.f)
c       1     frendly 
c       2     build (2nd, 3rd calls)
c       3     build (first call) 
c       4     actcor
c       5     ctransf
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character tag*4, string*140, key*22, values*80, strg*80

      integer option, i, j, ier, iscan

      double precision sum, ssum

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
  
      character vname*8, xname*8
      common/ csta2 /xname(k5),vname(l2)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct
c-----------------------------------------------------------------------
      rewind n2
c                                 frendly or actcor is reading
c                                 the header, no transformations
      if (option.eq.1.or.option.eq.4) itrans = 0
c                                 test for old (pre 4/2010) data file format:
      read (n2,*,iostat=ier) i
      if (ier.eq.0) call error (8,r,i,dname)

      rewind n2
c                                 database name
      call getkey (n2,ier,key,values,strg)

      dname = strg
c                                 extrinsic variable names & reference values
c                                 read "begin"
      call getkey (n2,ier,key,values,strg)

      do i = 1, l2

         call getkey (n2,ier,key,values,strg)
c                                 over ride default mobile component naming to 
c                                 allow build to use fugacities or activities.
         if (option.gt.3.or.i.lt.4) read (key,'(a8)') vname(i)
         read (values,*) v(i), delt(i)

      end do
c                                 reset the chemical potential tolerances if the
c                                 potential is to be replaced by activity/fugacity
      do i = 1, ipot
         if (jv(i).gt.3) then
            if (imaf(jv(i)-3).ne.1) delt(jv(i)) = delt(jv(i))/1e4
         end if
      end do 
c                                 set log p variable name.
      if (icopt.gt.4.and.lopt(14)) vname(1) = 'log[P,b]'
c                                 set log p variable name.
      if (icopt.gt.4.and.lopt(37)) vname(3) = 'log[X_f]'
c                                 read end key
      call getkey (n2,ier,key,values,strg)
c                                  set reference conditions
      pr = v(1)
      tr = v(2)
c                                  this block of code probably never gets executed?
      if (option.lt.4) then 
         if ((ifug.ge.10.and.ifug.le.12).or.
     *       ifug.eq.15.or.ifug.eq.17.or.ifug.eq.18) then 
            vname(3) = ' X(O) ' 
         else if (ifug.eq.25) then 
            vname(3) = 'Y(CO2)*'
         else if (ifug.eq.13) then 
            vname(3) = 'X(H2)'
         end if 
      end if  
c                                 read tolerance dtol
      call getkey (n2,ier,key,values,strg)

      read (values,*) dtol
      dtol = -dabs(dtol)
c                                 set utol and ptol, the 
c                                 tolerances (in energy units) for determination of
c                                 the stability of divariant, univariant
c                                 and invariant equilibria or reactions.
c                                 utol must be smaller than -utol
c                                 ptol must be > 2*-dtol
      utol = -dtol/1d1
      ptol = -dtol*3d0 

      hscon = .false.
      oxchg = .false.

      do i = 1, k0
         sel(i) = 0d0
         cox(i) = 0d0
      end do

      do
c                                 component names, formula weights, and, optionally, 3rd law s_elements
         call getkey (n2,ier,key,values,strg)
c                                 look for optional HSC_conversion key
         if (key.eq.'HSC_conversion') then 

            hscon = .true.

         else if (key.eq.'reference_oxidation_st') then

            oxchg = .true.

         else if (key.eq.'begin_components') then 

            exit 

         else

            call error (72,utol,i,'invalid thermodynamic data file '//
     *                             'keyword '//key)

         end if

      end do 


      icmpn = 0 

      do

         call getkey (n2,ier,key,values,strg)
         
         if (key.eq.'end_components') exit

         icmpn = icmpn + 1
c                                 get component string length
         cl(icmpn) = iscan(1,length,' ') - 1

         if (hscon.and.oxchg) then

            read (strg,*,iostat=ier) cmpnt(icmpn), atwt(icmpn), 
     *                   sel(icmpn), cox(icmpn), dispro(icmpn)
            if (ier.ne.0) then
               dispro(icmpn) = .false.
               read (strg,*) cmpnt(icmpn), atwt(icmpn), sel(icmpn), 
     *                       cox(icmpn)
            end if

         else if (hscon) then

            read (strg,*,iostat=ier) cmpnt(icmpn), atwt(icmpn), 
     *                               sel(icmpn), dispro(icmpn)
            if (ier.ne.0) then
               dispro(icmpn) = .false.
               read (strg,*) cmpnt(icmpn), atwt(icmpn), sel(icmpn)
            end if

         else

            read (strg,*,iostat=ier) cmpnt(icmpn), atwt(icmpn), 
     *                               dispro(icmpn)

            if (ier.ne.0) then
               dispro(icmpn) = .false.
               read (strg,*) cmpnt(icmpn), atwt(icmpn)
            end if

         end if

      end do
c                                 save old names for component transformations
c                                 in ctransf or build.
      do i = 1, k0
         xcmpnt(i) = cmpnt(i)
      end do
c                                 for programs that read the composition 
c                                 vector, check that the cp-dimension (k5)
c                                 is adequate
c     if (option.ne.0.and.option.ne.2.and.icmpn.gt.k5) 
c     *    call error (197,r,icmpn,'TOPN2')
c                                 read special components, override by lopt(63)
      lopt(7) = .false.
  
      call getkey (n2,ier,key,values,strg)

      if (key.eq.'begin_special_componen') then

         ispec = 0 

         do 
 
            call getkey (n2,ier,key,values,strg)
         
            if (key.eq.'end_special_components') exit

            if (.not.lopt(63)) then 
               do j = 1, icmpn
                  if (key.eq.cmpnt(j)) then 
                     ispec = ispec + 1
                     idspe(ispec) = j
                     lopt(7) = .true.
                     exit
                  end if 
               end do 
            end if
         end do 

      else
c                                 no saturated phase constraint possible
c                                 set lopt(7) to false. 
         backspace(n2)

      end if 

      if (option.eq.3.or.option.eq.5) then 
c                                 get transformations if build or ctransf
         call gettrn (option)
c                                 check special components
         if (lopt(7)) then 
            j = 0
            do i = 1, ispec
               if (idspe(i).ne.0) then
                  j = j + 1
                  idspe(j) = idspe(i)
               end if
            end do 
            ispec = j
            if (ispec.eq.0) lopt(7) = .false.
         end if 

      else if (option.ne.2) then 
c                                 substitute transformed component names
c                                 and compute the new formula wieghts, this
c                                 is done in gettrns for ioption 3/5.
         do i = 1, itrans

            cmpnt(ictr(i)) = tcname(i)

            sum = 0d0
            ssum = 0d0 

            do j = 1, icmpn
               sum =  sum  + ctrans(j,i) * atwt(j)
               ssum = ssum + ctrans(j,i) * sel(j)
            end do
 
            atwt(ictr(i)) = sum
            sel(ictr(i))  = ssum

         end do 

      end if 

      if (option.gt.3) then 
c                                 echo formatted header data for ctransf/actcor:
         write (n8,1000) 
         write (n8,'(a,a,/)') dname,' |<= data base title'

         write (n8,'(a,a)') 'begin_standard_variables |<= name (<9 ',
     *                      'characters), reference value, tolerance'
         do i = 1, l2 
            write (n8,'(a8,1x,f7.2,3x,g6.1E1)') vname(i),v(i),delt(i)
         end do

         write (n8,'(a,/)') 'end_standard_variables'

         write (n8,'(a,g6.1E1,a,/)') 'tolerance  ',dtol,
     *         '  |<= DTOL for unconstrained minimization, energy units'

         if (hscon) then
            write (n8,'(a,//,a)') 'HSC_conversion |<= tag enabling HSC '
     *                          //'to SUP apparent energy conversion, '
     *                          //'requires elemental entropies in the '
     *                          //'component list below',
     *                         'begin_components | < 6 chars, '//
     *                         'molar weight (g), elemental entropy (R)'
            write (n8,'(a5,2x,f9.4,3x,f9.4)') (cmpnt(i),atwt(i),sel(i),
     *                                        i = 1, icmpn)

         else 

            write (n8,'(a)') 'begin_components | < 6 chars, '//
     *                       'molar weight (g)'
            write (n8,'(a5,1x,f9.4)') (cmpnt(i),atwt(i), i = 1, icmpn)

         end if 

         write (n8,'(a,/)') 'end_components'

         if (lopt(7)) then 
            write (n8,'(a)') 'begin_special_components'
            do i = 1, ispec
               write (n8,'(a)') cmpnt(idspe(i))
            end do 
            write (n8,'(a,/)') 'end_special_components'
         end if  

      end if 
c                                 read and echo unformatted comments and make data                            
      do 

         read (n2,'(a)',iostat=ier) string
         if (ier.ne.0) call error (21,r,i,dname)
         read (string,'(a)') tag

         if (option.gt.3) then
            call mytrim (string)
            write (n8,'(400a)') chars(1:length)
         end if 

         if (string.eq.'begin_makes'.and.option.lt.4) then
 
            call rmakes (option)

            cycle 

         else if (tag.ne.'end') then     

            cycle

         else

            exit

         end if

      end do  

1000  format (/,' | comments are indicated by the | character.',/,
     *     ' | check for warnings at the end of the header section.',/) 

      end 

      double precision function dnan()
c----------------------------------------------------------------------
c george helffrich's function to make a safe NaN
c----------------------------------------------------------------------
      implicit none 

      integer i4

      character c8(8)

      double precision val

      equivalence (c8,i4,val)
c----------------------------------------------------------------------
      i4 = 1

      if (ichar(c8(1)).eq.1) then
C                                 Little-endian
         c8(8)=char(127)
         c8(7)=char(248)
         c8(6)=char(0)
         c8(5)=char(0)
         c8(4)=char(0)
         c8(3)=char(0)
         c8(2)=char(0)
         c8(1)=char(0)

      else
C                                 Big-endian
         c8(1)=char(127)
         c8(2)=char(248)
         c8(3)=char(0)
         c8(4)=char(0)
         c8(5)=char(0)
         c8(6)=char(0)
         c8(7)=char(0)
         c8(8)=char(0)

      endif

      dnan = val

      end

      subroutine plblrb (typ)
c----------------------------------------------------------------------
c write a blurb on plotting options for the calculation
c
c type 1 - tab, ctr, table format
c type 2 - plt, psvdraw format
c type 3 - pts, pspts format
c type 4 - phm, phemgp table format
c----------------------------------------------------------------------
      implicit none

      integer typ
c----------------------------------------------------------------------

      if (typ.eq.1) then
c                                 2d - tab format
         write (*,1000) 
         write (*,1010)
      else if (typ.eq.2) then 
c                                 plt format
         write (*,1020)
      else if (typ.eq.3) then
c                                 pts format
         write (*,1030)
      else if (typ.eq.4) then 
c                                 1d - tab format
         write (*,1000) 
         write (*,1040)
      end if 

1000  format (/,'The tabulated data from this calculation can be ',
     *          'plotted with:',/)
1010  format (5x,'PSTABLE - a Perple_X plotting program',
     *      /,5x,'PYWERAMI - github.com/ondrolexa/pywerami',
     *      /,5x,'PERPLE_X_PLOT - a MATLAB plotting script',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *           'for details of the table format refer to:',/,
     *      /,5x,'perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format',
     *           '.txt',/)
1020  format (/,'The output from this calculation can be plotted with ',
     *          'PSVDRAW',/)
1030  format (/,'The output from this calculation can be plotted with ',
     *          'PSPTS or converted to',/,'table/plot format with ',
     *          'PT2CURV',/)
1040  format (5x,'pstable - a Perple_X plotting program',
     *      /,5x,'perple_x_plot - a Matlab plotting script',
     *      /,5x,'spread-sheet programs, e.g., Excel',//,
     *          'for details of the table format refer to:',/,
     *      /,5x,'perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format',
     *           '.txt',/)

      end

c-------------------------------------------------------------------
c text manipulation routines
c
c evntually all such rountines will be in one block or file, this is
c not the case now. 6/20/2011. 
c-------------------------------------------------------------------


      subroutine unblnk (text)
c------------------------------------------------------------------- 
c unblnk - subroutine to remove blanks from text
 
c     text - character string 
c     length - length of unblanked character string, 0 on input 
c             if unknown.
c-------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar

      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c-------------------------------------------------------------------
      nchar = len(text)
 
      read (text,'(400a)') (chars(i), i = 1, nchar)
c                                 scan for blanks:
      length = 0

      do i = 1, nchar
         if (chars(i).eq.' ') cycle 
         length = length + 1
         chars(length) = chars(i)
      end do 

      write (text,'(400a)') (chars(i), i = 1, length)

      end

      subroutine enblnk (text)
c---------------------------------------------------------------------- 
c enblnk - scan text to first blank and cut remaining text.
 
c     text - character string 
c     jchar - length of unblanked character string, 0 on input 
c             if unknown.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ict,nchar
 
      character text*(*), bitsy(lchar)*1 

      nchar = len(text) 
      read (text,'(400a)') (bitsy(i), i = 1, nchar)
c                                 scan for blanks:
      ict = 0

      do i = 1, nchar
         if (bitsy(i).eq.' ') exit
         ict = ict + 1
      end do 

      text = ' '
      write (text,'(400a)') (bitsy(i), i = 1, ict)

      end

      subroutine reblnk (text)
c------------------------------------------------------------------- 
c reblnk - subroutine to replace blanks followed by a character
c          with an underscore.
 
c     text - character string 
c     jchar - length of unblanked character string, 0 on input 
c             if unknown.
c----------------------------------------------------------------------
      implicit none 

      integer i,ict
 
      character*8 text, bitsy(8)*1 
 
      read (text,'(400a)') bitsy
c                                 scan for blanks:
      ict = 0

      do i = 1, 7

         if (i.eq.1.and.bitsy(i).eq.' ') cycle

         if (bitsy(i).eq.' '.and.bitsy(i+1).ne.' ') then
            ict = ict + 1
            bitsy(ict) = '_'
         else if (bitsy(i).eq.' ') then 
            cycle
         else
            ict = ict + 1
            bitsy(ict) = bitsy(i)
         end if

      end do 
 
      bitsy(ict+1) = bitsy(8)

      write (text,'(400a)') (bitsy(i), i = 1, ict + 1)

      end

      subroutine ftext (ist,iend)

c subprogram to filter blanks from text 

      implicit none

      include 'perplex_parameters.h'

      integer ist,iend,i,itic,igot,jend

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      itic = ist - 1
      igot = 0
      jend = iend 

      do i = ist, iend-1

         if (chars(i).eq.' '.and.chars(i + 1).eq.' '.or. 
     *       chars(i).eq.' '.and.chars(i + 1).eq.')'.or. 
     *       chars(i).eq.' '.and.chars(i + 1).eq.'('.or.
     *       igot.eq.0.and.chars(i).eq.' ') cycle
         if (i.gt.ist) then
            if (chars(i-1).eq.'-'.and.chars(i).eq.' ') cycle 
         end if 

         itic = itic + 1
         igot = 1
         chars(itic) = chars(i)

      end do 

      if (chars(iend).ne.' ') then
         itic = itic + 1
         chars(itic) = chars(iend)
      end if 

      iend = itic + 1

      do i = iend, jend 
         chars(i) = ' '
      end do 
 
      end

      subroutine leblnk (text,ibeg,nchar)
c----------------------------------------------------------------------
c leblnk - scan text and strip leading blanks load into chars from
c position ibeg, nchar is the last non blank in char
c non-blank character in stripped string
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ibeg, nchar, ist
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      nchar = len(text) + ibeg -1 
      if (nchar.gt.lchar) nchar = lchar

      read (text,'(400a)') chars(ibeg:nchar)
c                                 find last non-blank   
      
      do i = ibeg, nchar
         if (chars(i).gt.' ') exit
      end do  

      ist = i 

      if (ist.gt.nchar) then 
c                                 all blank
         nchar = 0

      else 

         if (ist.gt.ibeg) then 
c                                 shift chars ist-1 left
            do i = ist, nchar
               chars(i+ibeg-ist) = chars(i)
            end do
         end if  
c                                 find last non-blank
         nchar = nchar + ibeg - ist

         do i = nchar, ibeg, -1
            if (chars(i).gt.' ') exit
         end do 

         nchar = i 

      end if 

      end

      integer function iscan (ibeg,iend,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of char in chars(ibeg..iend), if not
c found iscan is ?; kscan does the forward/inverse scans and sets value
c for bad ranges
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character char*1

      integer ibeg, iend

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      do iscan = ibeg, iend

         if (chars(iscan).eq.char) exit

      end do 

      end 

      integer function kscan (ibeg,iend,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of char in chars(ibeg..iend), if not
c found iscan is iend + inc
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character char*1

      integer ibeg, iend, inc

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      if (ibeg.gt.iend) then 
         inc = -1
      else 
         inc = 1
      end if 

      do kscan = ibeg, iend, inc 

         if (chars(kscan).eq.char) exit

      end do 
c                                 normal loop exit kscan = iend + inc
      end 

      integer function iscnlt (ibeg,iend,char)
c----------------------------------------------------------------------
c iscan finds the first occurence of a character in chars(ibeg..iend) that
c is greater than char. assuming ascii collating sequence +/- < 0 < a
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character char*1

      integer ibeg, iend, inc

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      if (ibeg.le.iend) then 
         inc = 1
      else 
         inc = -1
      end if 

      do iscnlt = ibeg, iend, inc

         if (chars(iscnlt).gt.char) exit

      end do 

      end 

      integer function jscnlt (ibeg,iend,char,chars)
c----------------------------------------------------------------------
c iscan finds the first occurence of a character in chars(ibeg..iend) that
c is greater than char. assuming ascii collating sequence +/- < 0 < a
c----------------------------------------------------------------------
      implicit none

      character char*1, chars(*)*1

      integer ibeg, iend, inc
c----------------------------------------------------------------------

      if (ibeg.le.iend) then 
         inc = 1
      else 
         inc = -1
      end if 

      do jscnlt = ibeg, iend, inc

         if (chars(jscnlt).gt.char) exit

      end do 

      end 

      integer function jscan (ibeg,iend,char,chars)
c----------------------------------------------------------------------
c jscan finds the first occurence of char in chars(ibeg..iend)
c----------------------------------------------------------------------
      implicit none

      character char*1, chars(*)*1

      integer ibeg, iend
c----------------------------------------------------------------------
      do jscan = ibeg, iend

         if (chars(jscan).eq.char) exit

      end do 

      end 

      subroutine mytrim (text)
c----------------------------------------------------------------------
c mytrim - scan text and delete trailing blank characters.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c---------------------------------------------------------------------- 
      nchar = len(text) 

      read (text,'(400a)') chars(1:nchar)
c                                find last non-blank
      length = 1 
      
      do i = 1, nchar
         if (chars(i).gt.' ') length = i
      end do

      end 

      subroutine deblnk (text)
c----------------------------------------------------------------------
c deblnk - scan text and delete multiple blank characters, strip
c out sequential + - or - + operators, trailing -/+ operators, 
c leading blanks, and leading + operators. also deletes blanks
c preceding punctuation (',' and ':' and ';').
 
c     text - character string 

c no test is made for a blank string or a string of "+" signs.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar

      logical strip
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c---------------------------------------------------------------------- 
      nchar = len(text) 

      read (text,1000) chars(1:nchar)
c                                find last non-blank
      length = 1 
      
      do i = 1, nchar
         if (chars(i).gt.' ') length = i
      end do

      nchar = length
c                                 kill any trailing +/- or ','
      if (chars(nchar).eq.'+'.or.chars(nchar).eq.'-'.or.
     *    chars(nchar).eq.',') nchar = nchar - 1
         
c                                 scan for first non blank/+ character:
      length = 0 
      
      do i = 1, nchar
         if (chars(i).eq.' '.or.chars(i).eq.'+') cycle
         length = i
         exit 
      end do 
c                                 shift everything right
      if (length.gt.1) then 

         length = length - 1
         
         do i = length + 1, nchar
            chars(i-length) = chars(i)
         end do 

         nchar = nchar - length

      end if 

      length = 1
      
      do i = 2, nchar
c                                 strip out double blanks
         if ((chars(i).eq.' '.and.chars(i+1).eq.' ').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.':').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.';').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.',').or.
     *       (chars(i).eq.' '.and.chars(i+1).eq.')')) cycle

         length = length + 1
         chars(length) = chars(i)

      end do

      nchar = length

      if (nchar.eq.1) return
c                                 strip put + - and - + strings
      strip = .false.

      do i = 1, nchar - 2

         if (chars(i).eq.'+'.and.chars(i+1).eq.'-'.or.
     *       chars(i).eq.'-'.and.chars(i+1).eq.'+') then

             chars(i) = '-'
             chars(i+1) = ' '
             strip = .true.

         else if (chars(i).eq.'+'.and.chars(i+2).eq.'-'.or.
     *            chars(i).eq.'-'.and.chars(i+2).eq.'+') then
c                                allow +/- or -/+
             if (chars(i+1).eq.'/') cycle 

             chars(i) = '-'
             chars(i+2) = ' '
             strip = .true.

         end if 

      end do 
c                                 special cases:
      if (nchar.gt.2) then 
         if (chars(nchar).eq.'*'.and.chars(nchar-1).eq.' '.and.
     *                               chars(nchar-2).eq.',') then
             chars(nchar-2) = '*'
             chars(nchar) = ' '
         end if
      end if

      if (strip) then 
c                                 strip out new double blanks
         length = 1

         do i = 2, nchar

            if (chars(i).eq.' '.and.chars(i-1).eq.' ') cycle

            length = length + 1
            chars(length) = chars(i)

         end do 

      end if 

      write (text,1000) chars(1:length),(' ',i = length+1, len(text))
 
1000  format (400a)

      end

      subroutine getstg (text)
c---------------------------------------------------------------------- 
c getstg - subroutine to get first non-blank string  
c          
c     text - character string on input, first non-blank strg on output
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar, ist
 
      character text*(*)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      nchar = len(text) 
      if (nchar.gt.lchar) nchar = lchar

      read (text,'(400a)') chars(1:nchar)
c                                 scan for blanks:
      ist = 1

      do i = 1, nchar
         if (chars(i).eq.' ') cycle 
         ist = i
         exit 
      end do 

      do i = ist, nchar
         if (chars(i).ne.' ') cycle 
         nchar = i-1
         exit 
      end do 

      text = ' '

      write (text,'(400a)') chars(ist:nchar)

      end

      subroutine redlpt (coeffs,ibeg,iend,ier)
c----------------------------------------------------------------------
c redlpt - read coefficients of a linear p-t function:

c                f = c0 + c1 T + c2 P

c from chars array. 

c on input

c    ibeg - the first possible location of the data

c assumes one of two formats:

c    c0 c1 c2

c or

c    c0 +/- ci Tag ....

c where the first character of a < 9 character Tag is used to identify P or T coefficients
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ibeg, iend, ier, iscan, iscnlt, itag

      double precision coeffs(3)

      external iscan, iscnlt

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      do i = 2, 3
         coeffs(i) = 0d0
      end do 
c                                 scan for an equals sign from ibeg
      iend = iscan (ibeg,com,'=') + 1
      if (iend.lt.com) ibeg = iend
c                                 get the first number
      ibeg = iscnlt (ibeg,com,' ') 

      call readfr (coeffs(1),ibeg,iend,com,ier)
      if (ier.ne.0.or.iend+1.ge.com) return

      ibeg = iend + 2
      itag = ibeg
c                                 try reading as though no tags 
c                                 are present (pre-6.7.3)
      do i = 2, 3

         call readfr (coeffs(i),ibeg,iend,com,ier)
         if (ier.ne.0) exit

      end do 

      if (ier.eq.0) return

      do i = 2, 3
         coeffs(i) = 0d0
      end do 
c                                 if an error, numbs/tags must be present
c                                 locate the number
      ibeg = itag
      iend = iscan (ibeg,com,' ') 
c                                 locate the first character of the tag
      itag = iend + 1

      if (chars(itag).eq.'T'.or.chars(itag).eq.'t') then
         i = 2
      else if (chars(itag).eq.'P'.or.chars(itag).eq.'p') then
c                                 must be c2, but check for invalid tag
         i = 3
      else 
         ier = 1
         return
      end if 
c                                 read the number
      call readfr (coeffs(i),ibeg,iend,com,ier)
c                                 the next number, if present begins at
      ibeg = iscan (itag,com,' ') + 1
      iend = iscan (ibeg,com,' ')

      if (ier.ne.0.or.iend.ge.com) return 
c                                 swap indexes
      if (i.eq.2) then 
          i = 3
       else 
          i = 2
      end if 
c                                 read the second tag
      call readfr (coeffs(i),ibeg,iend,com,ier)

      end 

      subroutine nanchk (x,y,text)
c----------------------------------------------------------------------
c nanchk - check x-y coordinate pair for bad_number and resets to 0.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical warn1

      character*(*) text

      double precision x, y

      save warn1
      data warn1/.true./
c----------------------------------------------------------------------
      if (warn1.and.(isnan(x).or.isnan(y))) then 
c                                 
         call warn (61,x,1,text)
         warn1 = .false.

      end if 

      if (isnan(x)) x = 0d0
      if (isnan(y)) y = 0d0

      end 

      subroutine begtim (itime)
c----------------------------------------------------------------------
c begin timer itime
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer itime
c----------------------------------------------------------------------
      call CPU_TIME(btime(itime))

      end 

      subroutine endtim (itime,output,chars)
c----------------------------------------------------------------------
c end timer itime with optional output
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical output

      character chars*(*)

      integer itime
c----------------------------------------------------------------------

      call CPU_TIME(etime(itime))
      times(itime) = times(itime) + (etime(itime)-btime(itime))

      if (output) then 

      write (*,'(/,a,3(2x,g14.7))') chars,
     *                              times(itime), 
     *                              etime(itime)-btime(itime)
      write (666,'(/,a,3(2x,g14.7))') chars,
     *                              times(itime), 
     *                              etime(itime)-btime(itime)

      end if 

      end 

      subroutine inqopn (lun,fname)
c-----------------------------------------------------------------------
c subroutine to open temporary (deleteable) files named fname on unit
c lun. this routine is used to avoid timing conflicts that arise to the
c time lag between a runtime close request and its implementation by 
c WINDOWS.
c-----------------------------------------------------------------------
      include 'perplex_parameters.h'

      integer ier, lun

      logical lopen, lname

      character fname*(*)
c-----------------------------------------------------------------------
      open (lun, file = fname, iostat = ier, status = 'new')

      if (ier.ne.0) then

         open (lun, file = fname, iostat = ier)

         if (ier.ne.0) then 
c                                 this case is that WINDOWS hasn't yet
c                                 implemented a previous close request.
            write (*,'(2(/,a))') '**error ver099** unable to open '
     *                            //fname
     *      ,'check that the file is not being used by another program.'
            write (*,'(/,a,i3)') 'IOSTAT = ',ier

            inquire (lun, opened=lopen, named=lname, name=fname)

            if (lopen) then

               write (*,'(a,i3,a)') 'system or programming error: LUN ',
     *                              lun,'is already open'
               if (lname) write (*,'(a)') 'and attached to file: ',
     *                    fname
               call errdbg ('please report this error')

            end if

         else
c                                 this case is just that the file hasn't
c                                 been closed
            close (lun, status = 'delete')
            open (lun, file = fname)

         end if

      end if

      end

      subroutine outsei
c-----------------------------------------------------------------------
c output details of modulus computation for endmembers and solutions
c if lopt(50).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      character stag*9, ptag*8, made*12

      logical stx, notstx, lstx, lmake

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer eos
      common/ cst303 /eos(k10)

      integer ipoint,kphct,imyn
      common/ cst60  /ipoint,kphct,imyn

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      integer make
      common / cst335 /make(k10)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character prject*100,tfname*100
      common/ cst228 /prject,tfname
c-----------------------------------------------------------------------
      call mertxt (tfname,prject,'_seismic_data.txt',0)

      call inqopn (n8,tfname)

      stx = .false.
      notstx = .false.
      lmake = .false.

      write (n8,1233) valu(19),nopt(6),lopt(17),valu(15),nopt(1),
     *                valu(14),lopt(20),lopt(4),.false.,lopt(65),
     *                nopt(65)

      write (n8,1030)

      write (n8,'(/,a)') 'Endmembers and stoichiometric compounds:'
      write (n8,1000)

      do i = istct, ipoint

         if (eos(i).eq.5.or.eos(i).eq.6) then 
            stx = .true.
            lstx = .true.
         else
            lstx = .false.
            notstx = .true.
         end if

         if (iemod(i).eq.3) then

            ptag = 'explicit'
            stag = 'missing'

         else if (iemod(i).eq.0) then

            ptag = 'implicit'
            stag = 'missing'

         else if (iemod(i).eq.1) then

            ptag = 'implicit'
            stag = 'explicit'

         else if (iemod(i).eq.2) then

            ptag = 'explicit'
            stag = 'explicit'

         end if

         if (.not.lopt(17)) ptag = 'implicit'

         if (iopt(16).eq.1) then

            if (stag.eq.'missing') stag = 'Poisson'

         else if (iopt(16).eq.2) then 

            stag = 'Poisson'

         end if

         if (lstx) then
            if (iemod(i).ge.1) stag = 'implicit'
            ptag = 'implicit'
         end if

         made = ' '

         if (make(i).ne.0) then 
            made = 'made entity*'
            lmake = .true.
         end if 

         if (ifp(i).ne.0) stag = 'fluid'

         write (n8,1050) names(i),ptag,stag,made

      end do

      if (lmake) write (n8,1040)

      if (isoct.gt.0) then

         if (stx.and.notstx) write (n8,1010) 

         write (n8,'(/,a)') 'Solutions:'
         write (n8,1000)

         do i = 1, isoct

            if (.not.pmod(i).and..not.smod(i)) then

               ptag = 'implicit'
               stag = 'missing'

            else if (pmod(i).and..not.smod(i)) then

               ptag = 'explicit'
               stag = 'missing'

            else if (.not.pmod(i).and.smod(i)) then

               ptag = 'implicit'
               stag = 'explicit'

            else if (pmod(i).and.smod(i)) then

               ptag = 'explicit'
               stag = 'explicit'

            end if

            if (lopt(17).and.pmod(i)) ptag = 'explicit'

            if (iopt(16).gt.0.and.(.not.smod(i).or.iopt(16).eq.2))
     *                                stag = 'Poisson'

            if (stx) then
               ptag = 'implicit'
               if (stag.eq.'explicit') stag = 'implicit*'
            end if

            if (fp(i)) stag = 'fluid'

            if (lname(i).eq.'liquid') stag = 'liquid'

            write (n8,1050) fname(i),ptag,stag

         end do

         if (stx) write (n8,1020)

      end if 

      close (n8)

1000  format (/,20x,'  Bulk Mod    Shear Mod ',/,
     *          20x,'  ---------   ---------')
1010  format (/,'**warning ver119** this computation mixes inconsistent'
     *         ,' thermodynamic data',/,'the following table may not be'
     *         ,' reliable.',/)
1020  format (/,'*computed as the Reuss average of the implicit endmemb'
     *       ,'er shear moduli.')
1030  format ('In the tables below: implicit moduli are calculated rigo'
     *       ,'rously from the EoS,',/,'explicit moduli are computed '
     *       ,'from empirical functions provided in the',/
     *       ,'thermodynamic data file.',/)
1040  format (/,'*explicit moduli of made endmembers are computed as a '
     *       ,'linear combination of ',/,'the real endmembers specified'
     *       ,' in the corresponding make definition.',/)
1050  format (6x,a10,6x,a8,4x,a9,4x,a)
1233  format (/,'Seismic wavespeed computational options:',//,
     *        4x,'bounds                  ',a3,7x,'[VRH] HS',/,
     *        4x,'vrh/hs_weighting        ',f3.1,7x,'[0.5] 0->1',/,
     *        4x,'explicit_bulk_modulus   ',l1,9x,'[T] F',/,
     *        4x,'poisson_ratio           ',a3,7x,'[on] all off; ',
     *        'Poisson ratio = ',f4.2,/,
     *        4x,'seismic_output          ',a3,7x,'[some] none all',/,
     *        4x,'poisson_test            ',l1,9x,'[F] T',/,
     *        4x,'Anderson-Gruneisen      ',l1,9x,'[F] T',/,
     *        4x,'Tisza_test              ',l1,9x,'[F] T',/,
     *        4x,'fluid_shear_modulus     ',l1,9x,'[T] F',/,
     *        4x,'phi_d                   ',f4.2,6x,'[0.36] 0->1',/)

      end

      subroutine outgrd (loopx,loopy,jinc,lun,ind2)
c----------------------------------------------------------------------
c output grid data to the plot file, called by VERTEX and UNSPLT
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer loopx,loopy,jinc,i,j,jst,kst,kd,ltic,iend,lun,jgrd(l7),
     *        ind1, ind2, ier

      logical ext 

      character string*(lchar), name*170, text*3

      integer igrd
      common/ cst311 /igrd(l7,l7)

      integer iam
      common/ cst4 /iam

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c----------------------------------------------------------------------
      if (lun.ne.n4) then 
c                                 write interim result file list
         call mertxt (name,prject,'.irf',0)
         inquire (lun, opened = ext)

         open (lun, file = name, iostat = ier, position = 'append')

         if (ier.eq.0) then 

            if (refine) then 
               ind1 = 1
            else 
               ind1 = 0
            end if 

            write (lun,*) ind1, ind2

            close (lun)

c                                 writing interim blk file
            write (text,'(a,i1,i1)') '_',ind1, ind2
            call mertxt (name,prject,text,0)
            call mertxt (name,name,'.blk',0)

            open (lun, file = name)

            rewind (n5)
c                                 the length of text should be able to 
c                                 handle format 1010 in outbl1
            do

               read (n5,'(a)',end=99) name
               write (lun,'(a)') name

            end do

99          close (lun)
            backspace (n5)
c                                 and the interim plt file
            call mertxt (name,prject,text,0)
            call mertxt (name,name,'.plt',0)
            open (lun, file = name)

         else 

            write (*,1000)
            return

         end if

      end if

      if (iam.eq.1) then
         write (lun,*) loopx, loopy, jinc
      else 
c                                 unsplt, set jinc to flag for werami
         write (lun,*) loopx, loopy, -1
      end if
c                                 fill in grid
      do i = 1, loopx, jinc

         if (i.ne.1.and.igrd(i,1).eq.0) then 
            jgrd(1) = igrd(i-jinc,1)
         else 
            jgrd(1) = igrd(i,1)
         end if 

         kst = 1

20       jst = kst

         if (i.ne.1.and.igrd(i,jst).eq.0) then
            jgrd(jst) = igrd(i-jinc,jst)
         else 
            jgrd(jst) = igrd(i,jst)
         end if 

         kd = jgrd(jst)

         ltic = -1

         do j = jst, loopy

            if (i.ne.1.and.igrd(i,j).eq.0) then 
               jgrd(j) = igrd(i-jinc,j)
            else 
               jgrd(j) = igrd(i,j)
            end if 

            if (jgrd(j).eq.0.or.jgrd(j).eq.kd) then
               ltic = ltic + 1
               if (j.eq.loopy) write (lun,*) ltic,kd
            else 
               write (lun,*) ltic,kd
               kst = j
               goto 20
            end if 
         end do 
      end do
c                                 write assemblage list
      write (lun,*) iasct

      do i = 1, iasct
         write (lun,*) iavar(1,i),iavar(2,i),iavar(3,i)
         write (lun,*) (idasls(j,i), j = 1, iavar(3,i))
      end do 

      if (lun.ne.n4) then
c                                 close interim plt file
         close (lun)

      else if (io3.eq.0) then 
c                                 write assemblages to print file
         write (n3,'(/,1x,a,a,/)') 'Stable assemblages identified ',
     *                         'by assemblage index:'
         do i = 1, iasct
            call psbtxt (i,string,iend)
            write (n3,'(i4,a,400a)') i,' - ',chars(1:length)
         end do

      end if

1000  format (/,'**warning ver999** OS fileio error, the irf file is'
     *        ,' corrupt and interim results',/,'for this calculation'
     *        ,' will be unreadable unless the irf file is edited',/)

      end

      subroutine getvar  
c--------------------------------------------------------------------
c getvar makes a list of variables to be used for i/o:

c if icopt = 10 -> using nodal coordinates else, 

c if icopt =  9/11 -> using 2d frac coordinates else:

c one-dimensional diagram (oned = .true.) then 

c the vertical (real) axis is variable iv(2), the horizontal axis
c is dummy.

c two-dimensional diagram (oned = .false.) and not icopt = 9, then 

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      integer iam
      common/ cst4 /iam

      logical pzfunc
      integer ilay,irep,npoly,ord
      double precision abc0,vz,iblk
      common/ cst66 /abc0(0:mord,mpol),vz(6),iblk(lay,k5),ilay,
     *               irep(lay),npoly,ord,pzfunc

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)   

      character vname*8, xname*8
      common / csta2 /xname(k5),vname(l2)  

      logical oned
      common/ cst82 /oned

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------

      if (icopt.eq.7.and.fileio) then 
c                                 1d-fractionation with file input,
c                                 use nodal coordinates:
         vnm(1) = 'node #'
         vmn(1) = 1
         vmx(2) = 1d0
         vmn(2) = 0d0 
         vmx(1) = loopy 
         oned = .true.
         
         jvar = ipot + 1

         do i = 2, jvar
            vnm(i) = vname(jv(i-1))
         end do  

      else if (icopt.lt.9) then 

         jvar = ipot

         if (idep.ne.0) jvar = ipot + 1

         if (icont.eq.1) then 

            do i = 1, jvar
               vnm(i) = vname(jv(i))
               vmx(i) = vmax(jv(i))
               vmn(i) = vmin(jv(i))
               var(i) = vmin(jv(i))
            end do   

         else 

            if (icont.eq.2) then 

               jvar = jvar + 1

               vnm(1) = ' X(C1)  '
               vmx(1) = 1d0
               vmn(1) = 0d0

               do i = 2, jvar
                  vnm(i) = vname(jv(i-1))
                  vmx(i) = vmax(jv(i-1))
                  vmn(i) = vmin(jv(i-1))
                  var(i) = vmin(jv(i-1))
               end do   
 
            else 

               jvar = jvar + 2

               vnm(1) = ' X(C1)  '
               vmx(1) = 1d0
               vmn(1) = 0d0

               vnm(2) = ' X(C2)  '
               vmx(2) = 1d0
               vmn(2) = 0d0

               do i = 3, jvar
                  vnm(i) = vname(jv(i-2))
                  vmx(i) = vmax(jv(i-2))
                  vmn(i) = vmin(jv(i-2))
                  var(i) = vmin(jv(i-2))
               end do   

            end if 

         end if 

         if (oned) then 
c                                 make a fake y-axis for 1-d plots
            vmx(2) = 1d0
            vmn(2) = 0d0

         end if

      else if (icopt.eq.9) then 
c                                using non-thermodynamic coordinate frame
         vmn(1) = vz(4)
         vmx(1) = vz(5)

         if (iam.ne.1) then
c                                 vertex sets the number of y-nodes in frac2d
            ncol = loopy
         else 
c                                 pssect/werami get the number of nodes from the plt file
            loopy = ncol

         end if 

         if (flsh) then
c                                  flush calculations: 
            vnm(1) = 'Q,kg/m^2'
            vnm(2) = 'dz,m   '
c                                  set the base to 
            vmn(2) = vz(1)/2d0
            vmx(2) = vmn(2) + dfloat(ncol-1)*vz(1)

         else
c                                  frac2d calculations.
            vnm(1) = 'z0,m'
            vnm(2) = 'dz,m'
c                                  set y = 0 ti be the top
            vmx(2) = -vz(1)/2d0
            vmn(2) = vmx(2) - dfloat(ncol-1)*vz(1)

         end if

         jvar = 4

         do i = 3, 4
            vnm(i) = vname(jv(i-2))
         end do

      else if (icopt.eq.12) then 

         vnm(1) = 'n,alqt. '
         vnm(2) = 'node#      '

         vmn(2) = 1d0
         vmx(2) = dfloat(iopt(36)) + 1d0
         var(2) = 1d0
 
         vmn(1) = 0d0
         vmx(1) = nopt(36)*dfloat(iopt(36))
         var(1) = 0d0

         v(1) = vmin(1)
         v(2) = vmin(2)

         jvar = ipot + 2

         do i = 3, jvar
            vnm(i) = vname(jv(i-2))
            vmx(i) = vmax(jv(i-2))
            vmn(i) = vmin(jv(i-2))
            var(i) = vmin(jv(i-2))
         end do

      end if 

      end

      subroutine  mkcomp (jcomp,ids)
c----------------------------------------------------------------
c mkcomp makes the jcomp'th user defined compositional variable
c the first k5 compositions are reserved for chsprp, the remaining
c k5 compositions are for solvus testing

c the solution ids is associated with the composition.

c   jcx  - the number of components to define the numerator of
c          the composition.
c   jcx1 - the number of components to define the denominator of
c          the composition.
c   icps - the indices of the components (1..jcx,jcx+1...jcx1).
c   rcps - the cofficients on the compenents as indexed by icps.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*5 y*1, units*13, text*195, what*9, sym*1

      integer jcomp, ier, i, ids, kount

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character cname*5
      common/ csta4  /cname(k5)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer icps, jcx, jcx1, kds
      logical stol, savg, spec
      double precision rcps, a0
      common/ comps /rcps(k7,2*k5),a0(k7,2),icps(k7,2*k5),jcx(2*k5),
     *               jcx1(2*k5),kds(2*k5),stol(i11),savg(i11),spec(2*k5)

      integer spct
      double precision ysp
      character*8 spnams
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)
c----------------------------------------------------------------------
c                                choose components vs species
      write (*,1000) fname(ids)
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') then 
         spec(jcomp) = .true.
         what = ' species'
      else
         spec(jcomp) = .false.
         what = 'component'
      end if 
c                                set units for composition
      if (spec(jcomp)) then
         units = 'mole fraction'
         sym = 'y'
      else if (iopt(2).eq.0) then
         units = 'molar amount '
         sym = 'n'
      else 
         units = ' mass amount '
         sym = 'm'
      end if 
c                                get the composition to be contoured
10    if (lopt(22)) then
c                                with moronic constant 
         write (*,1100) sym, sym, sym, units, what, what
      else 
         write (*,1110) sym, sym, sym, units, what, what
c                                zero the constant
         a0(jcomp,1) = 0d0
         a0(jcomp,2) = a0(jcomp,1)
      end if 
  
      do 

         if (spec(jcomp)) then 
            write (*,1030) what,'numerator',k5+1
         else 
            write (*,1030) what//'s','numerator',k5+1
         end if 

         read (*,*,iostat=ier) jcx(jcomp)

         if (ier.ne.0.or.jcx(jcomp).lt.1) then
            write (*,1020)
            cycle 
         end if 

         exit

      end do 
c                                define the numerator
      do 

         write (*,1040) what,'numerator'

         if (spec(jcomp)) then
            write (*,1010) (i,spnams(i,ids), i = 1, spct(ids))
            kount = spct(ids)
         else 
            write (*,1010) (i,cname(i), i = 1, icomp)
            kount = icomp
         end if 

         read (*,*,iostat=ier) (icps(i,jcomp),rcps(i,jcomp), 
     *                                     i = 1, jcx(jcomp))
         do i = 1, jcx(jcomp)
            if (icps(i,jcomp).lt.1.or.icps(i,jcomp).gt.kount) then
               ier = 1
               exit 
            end if 
         end do 

         if (ier.ne.0) then
            write (*,1020)
            cycle 
         end if 

         exit 

      end do  

      if (lopt(22)) then 
         write (*,1050) 'a1'
         call rdnumb (a0(jcomp,1),0d0,i,0,.true.)
      end if 
c                                define the denominator
      do 

         if (spec(jcomp)) then 
            write (*,1030) what,'denominator',k5+1-jcx(jcomp)
         else 
            write (*,1030) what//'s','denominator',k5+1-jcx(jcomp)
         end if 

         write (*,1140)
         read (*,*,iostat=ier) jcx1(jcomp)

         if (ier.ne.0.or.jcx1(jcomp).lt.0) then
            write (*,1020)
            cycle 
         end if 
 
         jcx1(jcomp) = jcx(jcomp) + jcx1(jcomp)
        
         exit 

      end do 

      if (jcx1(jcomp).gt.jcx(jcomp)) then 

         do 

            write (*,1040) what,'denominator'

            if (spec(jcomp)) then
               write (*,1010) (i,spnams(i,ids), i = 1, spct(ids))
            else 
               write (*,1010) (i,cname(i), i = 1, icomp)
            end if 

            read (*,*,iostat=ier) (icps(i,jcomp),rcps(i,jcomp), 
     *                                 i = jcx(jcomp)+1, jcx1(jcomp))

            do i = jcx(jcomp)+1, jcx1(jcomp)
               if (icps(i,jcomp).lt.1.or.icps(i,jcomp).gt.icomp) then
                  ier = 1
                  exit 
               end if 
            end do 

            if (ier.ne.0) then
               write (*,1020)
               cycle 
            end if 

            if (lopt(22)) then 

               write (*,1050) 'a2'
               call rdnumb (a0(jcomp,2),0d0,i,0,.true.)
c                                show the user the composition: 
               write (*,1070)   

               if (spec(jcomp)) then
                  write (text,1120) a0(jcomp,1),(rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = 1, jcx(jcomp))
               else           
                  write (text,1130) a0(jcomp,1),
     *                         (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                                          i = 1, jcx(jcomp))
               end if 

            else 

               write (*,1070)  

               if (spec(jcomp)) then
                  write (text,1120) (rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = 1, jcx(jcomp))
               else          
                  write (text,1120) (rcps(i,jcomp),cname(icps(i,jcomp)),
     *                                            i = 1, jcx(jcomp))
               end if 

            end if  

            call deblnk (text)
            write (*,1150) text    
            write (*,*) '   divided by '

            if (lopt(22)) then 

               if (spec(jcomp)) then
                  write (text,1130) a0(jcomp,1),(rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               else  
                  write (text,1130) a0(jcomp,2),
     *                        (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               end if 

            else

               if (spec(jcomp)) then
                  write (text,1120) (rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               else  
                  write (text,1120) (rcps(i,jcomp),cname(icps(i,jcomp)),
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               end if 

            end if 

            call deblnk (text)
            write (*,1150) text 

            exit 

         end do 

      else 

         if (lopt(22)) then 

            if (spec(jcomp)) then
               write (text,1130) a0(jcomp,1),(rcps(i,jcomp),
     *                           spnams(icps(i,jcomp),ids),
     *                           i = 1, jcx(jcomp))
            else  
               write (text,1130) a0(jcomp,1),
     *                           (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                           i = 1, jcx(jcomp))
            end if 

         else

            if (spec(jcomp)) then
               write (text,1120) (rcps(i,jcomp),
     *                           spnams(icps(i,jcomp),ids),
     *                           i = 1, jcx(jcomp))
            else  
               write (text,1120) (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                                       i = 1, jcx(jcomp))
            end if 

         end if  

         call deblnk (text)
         write (*,1080) text 

      end if 
 
      write (*,1090)
      read (*,'(a)') y
      if (y.eq.'y'.or.y.eq.'Y') goto 10

      kds(jcomp) = ids

1000  format (/,'Define the composition in terms of the species/endmem',
     *          'bers of ',a,' (y/n)?',//,'Answer no to define a ',
     *          'composition in terms of the systems components.',/,
     *          'Units (mass or molar) are controlled by the ',
     *          'composition keyword in',/,'perplex_option.dat.')
1010  format (2x,i2,' - ',a)
1020  format (/,'Invalid input, try again:',/)
1030  format (/,'How many ',a,' in the ',a,' of the',
     *          ' composition (<',i2,')?')
1040  format (/,'Enter ',a,' indices and weighting factors for the '
     *        ,a,':')
1050  format (/,'Enter the optional constant ',a,' [defaults to 0]:')
1070  format (/,'The compositional variable is:')
1080  format (/,'The compositional variable is: ',a,/)
1090  format ('Change it (y/n)?')
1100  format (/,'Compositions are defined as a ratio of the form:',/,
     *        4x,'[a1 + Sum {w(i)*',a,'(i), i = 1, c1}] / ',
     *           '[a2 + Sum {w(i)*',
     *        a,'(i), i = c2, c3}]',/,15x,
     *        a,'(j)   = ',a,' of ',a,' j',/,15x,
     *        'w(j)   = weighting factor of ',a,' j (usually 1)',/,
     *    15x,'a1, a2 = optional constants (usually 0)')
1110  format (/,'Compositions are defined as a ratio of the form:',/,
     *        4x,' Sum {w(i)*',a,'(i), i = 1, c1} / Sum {w(i)*',
     *        a,'(i), i = c2, c3}',/,15x,
     *        a,'(j)   = ',a,' of ',a,' j',/,15x,
     *        'w(j)   = weighting factor of ',a,' j (usually 1)')
1120  format (15('+',1x,f4.1,1x,a5,1x))
1130  format (f4.1,1x,15('+',1x,f4.1,1x,a5,1x))
1140  format ('Enter zero to use the numerator as a composition.')
1150  format (/,a,/)  

      end

      subroutine rnam1 (iex,xnam,what)
c----------------------------------------------------------------------
c read a solution name (what = 0) compound name (what = 1) or either
c (what = 2) from console, return
c iex = -id if a compound
c iex = ikp if a solution
c iex = 0 if invalid choice
c----------------------------------------------------------------------
      implicit none

      integer iex, what

      character*10 xnam
c----------------------------------------------------------------------
      iex = 0

      do 

         if (what.eq.0) then 
            write (*,1040) 'solution' 
         else if (what.eq.1) then 
            write (*,1040) 'compound' 
         else 
            write (*,1040) 'solution or compound' 
         end if  

         read (*,'(a)') xnam

         call matchj (xnam,iex)

         if (iex.ne.0) exit 
         write (*,1100) xnam

      end do 

1040  format (/,'Enter ',a,' (left justified): ')
1100  format (/,'No such entity as ',a,', try again: ')

      end


      subroutine lpwarn (idead,char)
c----------------------------------------------------------------------
c write warning messages from lpnag as called by routine 'char',
c set flags ier and idead, the optimization is a total fail if
c idead set to 1.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer idead, iwarn91, iwarn42, iwarn90, iwarn01, iwarn02, 
     *        iwarn03, iwarn58

      character char*(*)

      double precision c

      save iwarn91, iwarn42, iwarn90, iwarn01, iwarn02, iwarn03, 
     *     iwarn58

      data iwarn91, iwarn42, iwarn90, iwarn01, iwarn02, iwarn03, 
     *     iwarn58/7*0/
c----------------------------------------------------------------------
c                                             look for errors
      if (idead.eq.2.or.idead.gt.4.and.idead.lt.8.and.
     *                      (lopt(64).or.iwarn91.lt.6)) then 
c                                             unbounded solution, or
c                                             other programming error.
         call warn (91,c,idead,char)

         call prtptx

         iwarn91 = iwarn91 + 1
         if (iwarn91.eq.5.and..not.lopt(64)) 
     *                        call warn (49,c,91,'LPWARN')

      else if (idead.eq.3.and.
     *                      (lopt(64).or.iwarn42.lt.6)) then 
c                                             no feasible solution
         call warn (42,c,idead,char)

         call prtptx

         iwarn42 = iwarn42 + 1
         if (iwarn42.eq.6.and..not.lopt(64)) 
     *                        call warn (49,c,42,'LPWARN')

      else if (idead.eq.4.and.
     *                      (lopt(64).or.iwarn90.lt.6)) then 
c                                             iteration count exceeded,
c                                             probable cause no feasible
c                                             solution.
         call warn (90,c,idead,char) 
         iwarn90 = iwarn90 + 1
         if (iwarn90.eq.5) call warn (49,c,90,'LPWARN')

      else if ((lopt(64).or.iwarn58.lt.11)
     *         .and.(idead.eq.58.or.idead.eq.59)) then 

         if (idead.eq.58) then 
            call warn (58,c,k21,char)
         else 
            call warn (58,c,k25,char)
         end if

         call prtptx

         iwarn58 = iwarn58 + 1

         if (iwarn58.eq.10.and..not.lopt(64)) 
     *                        call warn (49,c,58,'LPWARN')

      else if (idead.eq.101.and.
     *        (lopt(64).or.iwarn01.lt.10).and.lopt(32)) then

          iwarn01 = iwarn01 + 1
          call warn (100,c,101,'under-saturated solute-component.'
     *              //' To output result set aq_bad_result to 102')

          if (iwarn01.eq.10.and..not.lopt(64)) 
     *                        call warn (49,c,101,'LPWARN')

      else if (idead.eq.102.and.
     *        (lopt(64).or.iwarn02.lt.10).and.lopt(32)) then

         iwarn02 = iwarn02 + 1
         call warn (100,c,102,'pure and impure solvent phases '//
     *             'coexist within solvus_tolerance. '//
     *             'To output result set aq_bad_result to 101')

         call prtptx

          if (iwarn02.eq.10.and..not.lopt(64)) 
     *                        call warn (49,c,102,'LPWARN')

      else if (idead.eq.103.and.
     *        (lopt(64).or.iwarn03.lt.10).and.lopt(32)) then

         iwarn03 = iwarn03 + 1
         call warn (100,c,103,'pure and impure solvent phases '//
     *              'coexist. To output result set aq_bad_result.')

         call prtptx

          if (iwarn03.eq.10.and..not.lopt(64)) 
     *                        call warn (49,c,103,'LPWARN')

      end if

      end 

      subroutine muwarn (quit,iter)
c----------------------------------------------------------------------
c if optimization succeeds but no chemical potentials write warning that
c iteration will be aborted and a low quality result ouput
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer iwarn, iter

      logical quit

      save iwarn

      data iwarn/0/
c----------------------------------------------------------------------
      quit = .true.

      if (iwarn.lt.9.or.lopt(64)) then

         iwarn = iwarn + 1

         write (*,106) iter

         call prtptx

         if (iwarn.eq.10) call warn (49,0d0,106,'MUWARN')

      end if

106   format (/,'**warning ver106** chemical potentials could not be ',
     *       'determined after ',i2,' iterations.',/,
     *       'Iteration has been aborted and the ',
     *       'low quality result output.',/)

      end

      subroutine prtptx 
c----------------------------------------------------------------------
c print current independent variable values to console, intended for 
c warning messages.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      character tag*8

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont
c----------------------------------------------------------------------
      write (*,'(a,/)') 'Current conditions:'

      do i = 2, icont

         if (i.eq.2) then
            tag = 'X(C1)   '
         else
            tag = 'X(C2)   '
         end if

         write (*,1000) tag,cx(i-1)

      end do

      do i = 1, ipot
         write (*,1000) vname(iv(i)),v(iv(i))
      end do

      write (*,'(/)')

1000  format (5x, a,' = ',g14.7)

      end 

      subroutine psbtxt (id,string,iend)
c----------------------------------------------------------------------
c subprogram to write a text labels for bulk composition output 
c id identifies the assemblage

      implicit none

      include 'perplex_parameters.h'

      character string*(*), pname*14

      integer i, ist, iend, id, ids

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      iend = 0

      string = ' '

      ist = 1

      do i = 1, lchar
         chars(i) = ' '
      end do

      do i = 1, iavar(3,id)
             
         ids = idasls(i,id)

         call getnam (pname,ids) 

         ist = iend + 1
         iend = ist + 14
         read (pname,'(400a)') chars(ist:iend)

         call ftext (ist,iend)

      end do 

      write (string,'(400a)') chars(1:iend) 

      length = iend

      end 

      subroutine assort (jlist,ksol,np)
c----------------------------------------------------------------------
c sort solution phase indices in an assemblage so that the phases are 
c in the same order as entered in the input file.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical quit

      integer i, j, jd, np, ifnd, jlist(k5), ksol(k19,k19)
c----------------------------------------------------------------------
         ifnd = 0
         quit = .false.

         do i = 1, isoct
c                                 original solution position
            jd = solptr(i)

            do j = 1, np
c                                 look through the list to see if the solution 
c                                 occurs, if it does, count it and set pointer
c                                 in jlist
               if (jd.eq.ksol(j,1)) then

                  ifnd = ifnd + 1
                  jlist(ifnd) = j
                  if (ifnd.eq.np) then
                     quit = .true.
                     exit
                  end if

               end if

            end do

            if (quit) exit

         end do

      end

      subroutine fopen (n2name,prt,n9name,err)
c-----------------------------------------------------------------------
c open files for subroutine input1.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, err, tic 

      integer ier

      character n2name*100, prt*3, name*100, n9name*100

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iam
      common/ cst4 /iam

      save first

      data first/.true./
c----------------------------------------------------------------------
c                                 open thermodynamic data file
      call fopen2 (0,n2name)

      tic = .false.
      err = .false.

      if (iam.eq.3.or.iam.eq.7.or.iam.eq.14) then
c                                 use existing plt/blk files
c                                 iam - 14 - unsplt (local)

c                                 plt/blk files for werami/pssect opened 
c                                 later by redplt to allow interim results
         if (iam.eq.14) then 
c                                 open the plot file
            call mertxt (name,prject,'.plt',0)

            open (n4, file = name, iostat = ier, status = 'old')

            if (ier.ne.0) err = .true.
c                                 open assemblage file
            call mertxt (name,prject,'.blk',0)

            open (n5, file = name, iostat = ier, status = 'old')

            if (ier.ne.0) err = .true.

         end if

      else if (iam.eq.1.or.iam.eq.2.or.iam.eq.13.or.iam.eq.15) then 
c                                 iam -  1 - vertex
c                                 iam -  2 - meemum
c                                 iam - 13 - unsplt (global)
c                                 iam - 15 - convex

         if (first) then
            tic = .true.
            call mertxt (name,prject,'.dat',0)
            write (*,1160) name
            write (*,1170) n2name
         end if
c                                 open print/plot files if requested
         if (prt.ne.' '.and.prt.ne.'no_'.and.iam.ne.13) then 

            io3 = 0 
            call mertxt (name,prject,'.prn',0)
            open (n3, file = name)

         else

            io3 = 1
            name = 'none requested'

         end if

         if (first.and.iam.ne.2) then
c                                 plt output file
            io4 = 0
            call mertxt (name,prject,'.plt',0)
            if (iam.ne.13) write (*,1180) name

            open (n4, file = name, iostat = ier, status = 'new')
            if (ier.ne.0) then 
               open (n4, file = name)
               close (n4, status = 'delete')
               open (n4, file = name)
            end if

            write (*,1190) name

            if (iam.ne.15) then 
c                                 blk output file
               call mertxt (name,prject,'.blk',0)
               open (n5, file = name, iostat = ier, status = 'new')
               if (ier.ne.0) then 
                  open (n5, file = name)
                  close (n5, status = 'delete')
                  open (n5, file = name)
               end if

               write (*,1220) name

            end if

         else if (iam.ne.15) then 

            rewind (n5)

         end if

      else

         call error (999,0d0,n9,'oops fopen')

      end if 

      if (n9name.ne.' ') then

         io9 = 0 
c                                 open solution model file
         open (n9,file = n9name,iostat = ier,status = 'old')
         if (ier.ne.0) call error (120,0d0,n9,n9name)

         if (tic) write (*,1210) n9name

      else

         io9 = 1
         if (tic) write (*,1210) 'not requested'

      end if

      first = .false.

1160  format (/,'Reading problem definition from file: ',a)
1170  format ('Reading thermodynamic data from file: ',a)
1180  format ('Writing print output to file: ',a)
1190  format ('Writing plot output to file: ',a)
1210  format ('Reading solution models from file: ',a)
1220  format ('Writing phase assemblage data to file: ',a)

      end

      subroutine setau1
c----------------------------------------------------------------------
c setau1 sets autorefine dependent parameters. called by vertex, werami,
c pssect, convex, and meemum.

c output is set to false if autorefine mode is not auto (i.e., iopt(6) = 2) 
c or it is auto and in the second cycle.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character y*1, badnam(h9)*10

      integer ibad2, ibad1, igood, i, j, ier

      character n10nam*100,n11nam*100,n8nam*100

      character prject*100,tfname*100
      common/ cst228 /prject,tfname
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
      refine = .false.
c                                 only use autorefine if solutions
c                                 are present and it is requested.
      if (isoct.ne.0) then 

         call mertxt (n10nam,prject,'.arf',0)
         open (n10, file = n10nam, iostat = ier, status = 'old')

         call mertxt (n8nam,prject,'.tof',0)

         if (iam.eq.1.or.iam.eq.2.or.iam.eq.15) then
c                                 VERTEX, MEEMUM, or CONVEX:
            if (iam.eq.1.or.iam.eq.15) call inqopn (n8,n8nam)

            ibad1 = 0 

            if (ier.ne.0.and.(iam.eq.1.or.iam.eq.15)) then 
c                                 no auto_refine data
               open (n10, file = n10nam, status = 'unknown')

            else if (ier.eq.0.and.(iam.eq.1.or.iam.eq.15)) then 

               if (iam.eq.15) then 
                  read (n10,*,iostat=ier) ibad1, ibad2, igood
                  if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)
               end if 
c                                 changed to .and. 9/10/19
               if (iopt(6).ne.2.and.outprt) write (*,1030) n10nam

               if (iopt(6).eq.1) then 
c                                 manual mode, allow reinitialization
c                                 or suppression.
                  write (*,1060) 
                  read (*,'(a)') y

                  if (y.eq.'y'.or.y.eq.'Y') then

                     iopt(6) = 0

                  else 

                     refine = .true.

                  end if

                  outprt = .true.
 
               else if (outprt) then  
c                                 second cycle of automated mode
                  refine = .true.

               end if  

               write (n8,*) refine

            else if (ier.eq.0.and.iam.eq.2.and.iopt(6).ne.0) then 
c                                 MEEMUM, ask the user if he wants
c                                 to use the data 
               write (*,'(/,a,a,/,a)') 'Auto-refine data exists from a',
     *                  ' previous calculation with VERTEX.',
     *                   'Do you want MEEMUM to use this data (y/n)?'
               read (*,'(a)') y

               if (y.ne.'y'.and.y.ne.'Y') then

                  iopt(6) = 0

               else 

                  refine = .true.
                  iopt(6) = 1

                  write (*,1030) n10nam

               end if

            else if (ier.ne.0.and.iam.eq.2) then 

               iopt(6) = 0

            end if 
c                                 set cycle dependent parameters
            if (refine) then 

               i = 2

            else 

               i = 1

            end if
c                                 auto solvus_tolerance, only relevant for CONVEX
            if (lopt(9).and.iam.eq.15) nopt(8) = 1.5d0*rid(3,i)

         else if (iam.eq.13) then
c                                 the global level of unsplt, which should generate 
c                                 neither arf or tof files; at this point if ier is
c                                 zero, the arf file has been successfully opened:
c                                 kill and close it:
            if (ier.ne.0) close (n10, status = 'delete')
c                                 open and kill the tof file
            open (n8, file = n8nam, status = 'unknown')
            close (n8, status = 'delete')
c                                 open and kill the irf file
            call mertxt (n8nam,prject,'.irf',0)
            open (n8, file = n11nam, iostat=ier, status = 'unknown')
            close (n8,status = 'delete')

         else
c                                 werami/pssect open and read tof file
            open (n8, file = n8nam, iostat=ier, status = 'old')

            if (ier.eq.0) then 
c                                 read flag that indicates if auto-refine
c                                 or exploratory stage parameters.
               read (n8,*,iostat=ier) refine

            else 

               call errdbg ('missing *.tof file')

            end if 

         end if
c                                 only want *_auto_refine.txt for the exploratory
c                                 stage. VERTEX or CONVEX:
         if (refine) then

            lopt(11) = .false.
            close (n11)

         else if (iam.eq.1.or.iam.eq.15.and.lopt(11)) then
c                                 user friendly text version of the exploratory stage
c                                 auto_refine file:
            call mertxt (n11nam,prject,'_auto_refine.txt',0)
            open (n11, file = n11nam, status = 'unknown')
c                                 write blurb
            write (n11,1000) 'www.perplex.ethz.ch/perplex/faq/warning_'
     *                       //'ver991_relax_solution_model_limits.txt'

         end if 

      end if 

      close (n8)
c                                 just to be sure
      if (iopt(6).eq.0) refine = .false.

      if (refine.and.iam.eq.15) then 
c                                 CONVEX: reject solution models that were 
c                                 not found to be stable and set parameters 
c                                 that depend on refinement
         ibad2 = 0 

         do 50 i = 1, isoct

            do j = 1, ibad1
               if (fname(i).eq.badnam(j)) then
                  if (iam.eq.1.or.iam.eq.15) write (*,1070) fname(i)
                  goto 50
               end if 
            end do 

            ibad2 = ibad2 + 1
            fname(ibad2) = fname(i)

50       continue 

         isoct = ibad2 

         write (*,'(/)')

      end if





      if (iopt(6).eq.2.and..not.refine) then
c                                 this means it must be in the exploratory
c                                 stage
         outprt = .false.

      else

         outprt = .true.

      end if

      if (.not.(iopt(6).eq.2.and.refine).and.iopt(34).ne.0.and.
     *    iam.eq.1) then
c                                  initialize (i.e., delete prior) intermediate 
c                                  results file 
         call mertxt (n11nam,prject,'.irf',0)
         open (1000, file = n11nam, iostat=ier, status = 'unknown')
         close (1000,status = 'delete')

      end if

1000  format (//,'NOTE: this file echoes the auto-refine data after ',
     *       'the exploratory stage. If',/,'the composition of a phase',
     *       ' has been relaxed (**warning ver991**) during this stage,'
     *    /,'best practice is to modify the appropriate subdivision sch'
     *      ,'eme* and repeat the',/,'exploratory stage calculation un'
     *      ,'til the warnings are eliminated. This process can be',
     *     /,'expedited by setting the auto_refine option = man or off',
     *    //,'For a summary of the compositional ranges at the end of',
     *       ' the auto-refine stage refer',/,'to the console output.'
     *   ,//,
     *       '*refer to the header section of the solution model file',
     *       'for explanation of subdivision schemes',//,  
     *       'and:',//,a,//'for additional information.',//)
1030  format (/,'Reading data for auto-refinement from file: ',a,/)
1060  format ('Suppress or reinitialize auto-refinement (y/n)?')
1070  format ('Eliminating solution model: ',a,' in auto-refinement.')

      end 

      subroutine setau2
c----------------------------------------------------------------------
c setau2 sets/resets autorefine parameters after the solution models have
c been read. setau1 must be called first.

c output is set to true if autorefine mode is auto (i.e., iopt(6) = 2) 
c but no solutions are present (isoct = 0). 
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,index

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      logical oned
      common/ cst82 /oned

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c-----------------------------------------------------------------------
      if (isoct.eq.0) then 
     
         index = 2
         outprt = .true.

      else if (.not.outprt) then

         index = 1

      else 

          if (refine) then

             index = 2

          else 

             index = 1

          end if 

      end if 
c                                 set auto-refine dependent parameters
      if (icopt.eq.5) then 
c                                 gridded minimization
         if (oned) then 
            jlow = grid(4,index)
            loopx = 1
         else 
            jlow = grid(2,index)
            loopx = grid(1,index) 
         end if

         jlev = grid(3,index) 
          
      else if (icopt.gt.5) then 
c                                 1d/2d phase fractionation
         jlow = grid(4,index)

      else if (icopt.eq.1) then 
c                                 schreinemakers diagrams

c                                 max variance of curves to be traced
          isudo = grid(5,index)
c                                 default variable tracing increment
          do i = 1, 2
             dv(iv(i)) = (vmax(iv(i)) - vmin(iv(i)))*rid(1,index)
          end do 

      else if (icopt.eq.3) then 
c                                 mixed variable diagrams 

c                                 no variance restriction
          isudo = 99
c                                 default search increment
          dv(iv(1)) = (vmax(iv(1)) - vmin(iv(1)))*rid(1,index)

      end if 

      end 

      subroutine input1 (first,err)
c-----------------------------------------------------------------------
c input1 reads data from a file on unit n1, this data controls the
c computational options and is modified frequently.

c iam - indicates calling program 1 - vertex
c                                 2 - meemum
c                                 3 - werami
c                                13 - unsplt, global call
c                                14 - unsplt, local call
c                                 any other values no output
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      logical eof, first, err

      character*100 blank*1,string(3)*8,rname*5,name*8,strg*80,n2name,
     *              n9name,y*1,sname*10,prt*3,plt*3

      integer idum, nstrg, i, j, k, ierr, icmpn, jcont, kct

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      character*100 cfname
      common/ cst227 /cfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      character*162 title
      common/ csta8 /title(4)

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character*8 xname,vname
      common/ csta2 /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5) 

      integer icp2
      common/ cst81 /icp2

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      double precision buf
      common/ cst112 /buf(5)

      integer ivfl
      common/ cst102 /ivfl

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer ibuf,hu,hv,hw,hx 
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      logical oned
      common/ cst82 /oned

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      character eoscmp*8
      common/ cst98 /eoscmp(2)

      integer iam
      common/ cst4 /iam

      save blank
      data blank/' '/
c-----------------------------------------------------------------------
c                             output = .false. then in 1st cycle of
c                             autorefine.
      if (.not.outprt) then 
c                                 read computational option file
         call fopen1
      
      else 
c                                 create the file name
         call mertxt (tfname,prject,'.dat',0)
         open (n1, file = tfname, iostat = ierr, status = 'old')
         if (ierr.ne.0) call error (120,r,n1,tfname)

      end if 
c                                 begin reading input:

c                                 read name of thermodynamic data file
      read (n1,'(a)') n2name
      call enblnk (n2name)
c                                 read print and graphic file names
      read (n1,'(a)') prt

      read (n1,'(a)') plt

      read (n1,'(a)') n9name
      call enblnk (n9name)
c
      do i = 1, 4
         title(i) = ' '
      end do 
c                                 read title for the calculation:
      read (n1,'(a)') title(1)
c                                 read computational option or option file name
c                                 use error condition to determine which:
      read (n1,'(a)') tfname
c                                 get first non-blank string 
      call getstg (tfname)

      read (tfname,'(i2)',iostat=ierr) icopt 

      if (ierr.eq.0) then 
c                                 if no error, old version
         tfname = 'perplex_option.dat'

      else
c                                 new version, read icopt
         read (n1,*,err=998) icopt

      end if 
c                                 if fractionation path from data 
c                                 file, get name:
      fileio = .false.

      if (icopt.eq.10.or.icopt.eq.11) then 

         fileio = .true.

         read (n1,'(a)') cfname
         call enblnk (cfname)

         if (icopt.eq.10) then 
            icopt = 7
         else 
            icopt = 9
         end if

      end if 
c                                 if meemum, override whatever computational option
c                                 is set in the input file. 
      if (iam.eq.2) icopt = 5
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum

      read (n1,*,err=998) itrans
      read (n1,*,err=998) icmpn
c                                 read new component definitions:
      do i = 1, itrans
         read (n1,'(a,1x,i2)') tcname(i), ictr(i)
         read (n1,*) (ctrans(j,i), j = 1, icmpn)
      end do

      read (n1,*,err=998) iwt
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum  
      read (n1,*,err=998) idum
c                                 read code for choice of fluid equation
c                                 of state from terminal. 
      read (n1,*,err=998) ifug
      
      if (ifug.eq.8 .or.ifug.eq.10.or.ifug.eq.12.or.ifug.eq.16.or.
     *    ifug.eq.17.or.ifug.eq.19.or.ifug.eq.20.or.ifug.eq.24.or.
     *    ifug.eq.25) then

         read (n1,*,err=998) ibuf,hu,dlnfo2,elag

      else if (ifug.eq.6 .or.ifug.eq.7 .or.ifug.eq.11.or.ifug.eq.18.or.
     *         ifug.eq.21.or.ifug.eq.22.or.ifug.eq.23) then

        call error (72,0d0,0,' the input file specifies a disabled '//
     *                       'or ivalid internal fluid EoS')

      end if 

      if (ibuf.eq.5) read (n1,*,err=998) buf

      if (hu.eq.1) then 
c                                 hardwired fluid EoS endmember names
         eoscmp(1) = 'H2      '
         eoscmp(2) = 'O2      '

      else 

         eoscmp(1) = 'H2O     '
         eoscmp(2) = 'CO2     '

      end if 
c                                 no dependent variable
      iind = 0 
c                                 dummy variable
      read (n1,*,err=998) idum
c                                 idum is just a 1d/2d flag for 
c                                 gridded minimization, for backwards 
c                                 compatibility set the to 2d if > 2 or < 1.
      if (idum.eq.1) then 
         oned = .true.
      else
         oned = .false.
      end if 

      read (n1,*,err=998) idep
      read (n1,*,err=998) c0,c1,c2,c3,c4

      if (idep.eq.1) then 
         iind = 2
      else if (idep.eq.2) then 
         iind = 1
      end if 
c                                 decode thermodynamic components
c                                 read to the beginning of the component list
      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 count (icp) and save names (cname)
      icp = 0
      jbulk = 0

      do 

         read (n1,'(a,a)') rname,strg

         if (rname.eq.'end t') then 
c                                 finished, check for no components
            if (icp.eq.0) then
               write (*,*) 'No thermodynamic components'
               goto 998
            else if (icopt.eq.5.and.jbulk.lt.icp) then 
               write (*,*) 'All thermodynamic components must be ',
     *                     'constrained.'
               goto 998
            end if 
        
            exit 

         else if (rname.eq.blank) then 
 
            cycle

         else

            icp = icp + 1
            cname(icp) = rname
c                                 encode a graphics names for the
c                                 compositional variables, this is kind of
c                                 pointless, but it looks good.
            write (xname(icp),'(a,a,a)') 'x(',rname,')'
c                                 unblank the name
            call unblnk (xname(icp))
            if (icp.gt.k5) call error (197,r,icp,'INPUT1')

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) icont

         if (icopt.eq.12) then 
            k = 2
         else 
            k = icont
         end if 

         if (k.ne.0) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, k)
         end if 

      end do

      icp1 = icp + 1
      icp2 = icp + 2

      hcp = icp
c                                 decode saturated components
c                                 isat is the saturated component counter
      isat = 0
      io2  = 0 

      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 

      do 

         read (n1,'(a,a)') rname,strg
         if (rname.eq.blank) cycle 

         if (rname.eq.'end s') then 

            icomp = icp + isat
            exit 

         else if (rname.eq.blank) then 

            cycle 

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) jcont

         if (jcont.ne.0) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, jcont)
         end if

         isat = isat + 1
         if (isat.gt.h5) call error (15,r,i,'BUILD')
         cname(icp+isat) = rname
         if (rname.eq.'O2') io2 = isat

      end do 
c                                 decode saturated phase components
      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 ifct is the saturated phase component counter
      ifct = 0

      do 

         read (n1,'(a)') rname

         if (rname.eq.'end s') then 
            icomp = icomp + ifct
            exit 
         else if (rname.eq.blank) then 
            cycle 
         end if 
      
         ifct = ifct + 1
         if (ifct.gt.2) call error (44,r,i,'BUILD')
c                                 save the component if only one
c                                 for use in input2.
         if (ifct.eq.1) zname = rname
         cname(icomp+ifct) = rname
      end do 
c                                  decode mobile components
c                                  jmct - mobile component counter
      jmct = 0 
      ifact = 0
      jmuct = 0 

      do 

         call rdstrg (n1,nstrg,string,eof)

         if (eof) then 

            goto 998

         else if (string(1).eq.'begin') then

            cycle 

         else if (string(1).eq.'end') then

            icomp = icomp + jmct
            exit 

         else 

            read (string(1),'(a5)') rname
            jmct = jmct + 1
            if (jmct.gt.2) call error (45,r,i,'BUILD')
            cname(icomp+jmct) = rname

            if (nstrg.eq.1) then 
c                                 old format, create variable name
               write (vname(3+jmct),'(a,a)') 'mu_',rname
               imaf(jmct) = 1
               jmuct = jmuct + 1

            else 
c                                 new format
               read (string(2),'(a1)') y
               vname(3+jmct) = string(2)
               afname(jmct) = string(3)

               if (y.eq.'m') then 
c                                 chemical potential
                  imaf(jmct) = 1
                  jmuct = jmuct + 1

               else if (y.eq.'f') then 

                  imaf(jmct) = 2

               else if (y.eq.'a') then 

                  imaf(jmct) = 3

               end if 

               if (imaf(jmct).gt.1) ifact = ifact + 1 

            end if 
               
         end if 

      end do 
c                             the ifct flag can probably be set later if fluid
c                             is in the thermodynamic composition space.   
      jfct = icp + isat 
c                             jprct+1..icomp -> (jmct.ne.0) mobile components 
      jprct = icomp - jmct
c                             kbulk counter used for aq speciation which allows
c                             saturated + mobile components
      kbulk = jbulk + jmct
c                             excluded phases
      ixct = 0
c                             decode excluded phases
      do 
         read (n1,'(a)',end=998) name
         if (name.eq.'begin ex') exit
      end do

      do 

        read (n1,'(a)') name

         if (name.eq.'end excl') then
            exit
         else if (name.eq.blank) then
            cycle 
         end if 

         ixct = ixct + 1
         if (ixct.gt.h8) call error (13,r,i,'BUILD')
         exname(ixct) = name

      end do  
c                             solution phases:
      do 
         read (n1,'(a)',end=998) sname
         if (sname.eq.'begin solu') exit
      end do
c                             isoct - solution phase counter,
c                             io9 is a flag = 0 no solution file
      isoct = 0

      do 

         read (n1,'(a)') sname
 
         if (sname.eq.'end soluti') then 
            if (io9.eq.1) isoct = 0 
            exit 
         else if (sname.eq.blank) then 
            cycle  
         end if 

         isoct = isoct + 1
         if (isoct.gt.h9) call error (25,r,i,'BUILD')
         fname(isoct) = sname

      end do  
c                             read the maximum pressure, temper-
c                             ature, xco2, u1, and u2; the minimum
c                             pressure temperature, xco2, u1, and u2;
c                             and the default pressure, temperature,
c                             xco2, and chemical
c                             potential increments use kelvins, bars and
c                             joules as units (if no mobile components
c                             enter two zeroes for each read).
      read (n1,*,err=998) vmax
      read (n1,*,err=998) vmin
      read (n1,*,err=998) dv
c                             read the default indices of the
c                             dependent, independent, and secondary
c                             independent intensive variables, p = 1,
c                             t = 2, and xco2 = 3, respectively.
      read (n1,*,err=998) (iv(i), i = 1, 5)
c                             check variable ranges are consistent,
c                             variable iv(1), UNLESS:
c                             normal calculations w/o limits
      if (icopt.ne.0.and.icopt.ne.4.and.icopt.ne.12.and.
c                             fractionation calculations
     *    icopt.ne.7.and.icopt.le.9.
c                              MEEMUM
     *    and.iam.ne.2) then

         if (iv(1).eq.3.and.ifct.eq.0) call error (110,r,i,'I')

         if (iv(1).eq.3.and.ifct.eq.1) then 

            if (icopt.ne.7.and.iv(2).ne.3) call error (111,r,i,'I')

         end if 

         if (vmin(iv(1)).ge.vmax(iv(1)).and.icopt.lt.5) then 

            call error (112,r,i,'less than or equal')

         else if (vmin(iv(1)).eq.vmax(iv(1)).and.
     *            (icont.lt.3.or.oned.and.icont.eq.1)) then

            call error (112,r,i,'equal')

         end if 

         if (vname(iv(1)).eq.blank) call error (116,vmin(1),i,'I')

      end if
c                             variable iv(2):
      if (iam.ne.2.and.(icopt.eq.1.or.
     *                  (icopt.eq.5.and.icont.eq.1.and..not.oned))) then

         if (iv(2).eq.3.and.ifct.eq.0) call error (110,r,i,'INPUT1')

         if (iv(2).eq.3.and.ifct.eq.1) call error (111,r,i,'INPUT1')

         if (icopt.eq.1) then 

            if (vmin(iv(2)).ge.vmax(iv(2))) call error (112,r,i,
     *                                            'less than or equal')

         else 

            if (vmin(iv(2)).eq.vmax(iv(2))) call error (112,r,i,'equal')

         end if 

         if (vname(iv(2)).eq.blank) call error (116,r,i,'INPUT1')

      end if
c                             if a chemical potential is specified as an
c                             independent variable (iv(1-3)), check if
c                             the variable is defined:
      kct = 0

      do i = 1, 3
         if (iv(i).gt.3) kct = kct + 1
      end do 
c                             identify the variable used to determine
c                             which phases lie on the left hand side
c                             of a reaction equation.
      if (icopt.eq.3) then
         ivfl = iv(1)
      else if (iv(1).eq.2.or.iv(2).eq.2) then
c                             choose T
         ivfl = 2
      else if (iv(1).eq.1.or.iv(2).eq.1) then
c                             no T, so choose P
         ivfl = 1
      else
c                             no P or T, choose independent V
         if (iv(2).ne.3) then
            ivfl = iv(2)
         else
            ivfl = iv(1)
         end if
      end if
c                             ok, now find out which variables are
c                             dummies and store the indexes of the
c                             non-dummy variables in jv.
      ipot = 0

      do i = 1, 5
c                             variables v(1) (p) and v(2) (t) are
c                             only dummies if idep is set.
         if ((iv(i).ne.idep.or.icopt.eq.7.or.icopt.eq.9).and.
     *       (iv(i).eq.1.or.iv(i).eq.2)) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variable v(3) is a dummy if ifct = 0:
         else if ((iv(i).eq.3).and.ifct.gt.0) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variables v(4) and v(4) are dummies if
c                             imyn = 1:
         else if (jmct.ne.0) then
            if (iv(i).eq.4) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            else if (iv(i).eq.5.and.jmct.eq.2) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            end if
         end if

      end do 
c                                 if dependent variable add to jv list, could
c                                 increment ipot, but maybe it's better not to.
      if (idep.ne.0) jv(ipot+1) = idep
c                                 set convergence criteria for routine univeq
      if (icopt.le.3) then 

         call concrt

      else if (icopt.eq.12) then 
c                                 0-d infiltration
         read (n1,*,err=998) iopt(36), nopt(36)

      end if 

      if (icopt.ne.0) close (n1)
c                                 open files requested in input
      call fopen (n2name,prt,n9name,err)
c                                 err only set for unsplt (iam.eq.14)
      if (err) return
c                                 read auxilliary input for 2d fractionation
      if (icopt.eq.9) call rdain
c                                 get runtime parameters
      if (first.or.(.not.first).and.(.not.outprt).or.iam.eq.13) 
     *   call redop1 (first,tfname)

      goto 999
c                                 archaic error trap
998   call mertxt (n2name,prject,'.dat',0)
      call error (27,r,i,n2name)

999   end

      subroutine setblk
c-----------------------------------------------------------------------
c for gridded minimization setblk computes the bulk composition
c and initializes the arrays for lpsol.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x0

      integer i,j

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)
c-----------------------------------------------------------------------
      x0 = 1d0

      if (lopt(1)) then 
c                                 closed composition
         do j = 1, icont-1
            x0 = x0 - cx(j)
         end do 

      end if 

      do j = 1, jbulk
         cblk(j) = x0*dblk(1,j)
      end do 
         
      do j = 1, jbulk
         do i = 2, icont 
            cblk(j) = cblk(j) + cx(i-1)*dblk(i,j)
         end do 
      end do
c                                 modify cblk here to change the 
c                                 composition before minimization.
      ctotal = 0d0 
c                                 get total moles to compute mole fractions             
      do i = 1, hcp
         ctotal = ctotal + cblk(i)
      end do

      do i = 1, hcp 
         b(i) = cblk(i)/ctotal
      end do

      end 

      subroutine inblnk (text,char)
c----------------------------------------------------------------------
c inblnk - scan text to last '/' or '\' and insert char after.
 
c     text - character string 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nchar
 
      character text*(*), bitsy(lchar)*1, char*1 
c----------------------------------------------------------------------
      nchar = len(text) 
      read (text,1000) (bitsy(i), i = 1, nchar)
c                                 scan for blanks:

      do i = nchar,1,-1
c                                 this line may cause problems
c                                 on some operating systems that 
c                                 recognize the backslash as an escape
c                                 character.
         if (bitsy(i).eq.'/') goto 10
         bitsy(i+1) = bitsy(i)
      end do 

      i = 0

10    bitsy(i+1) = char

      write (text,1000) (bitsy(i), i = 1, nchar)
 
1000  format (400a)

      end

      subroutine matchj (unnown,itis)
c----------------------------------------------------------------------
 
c matchj - subroutine to determine if the string unnown is a valid
c          solution or compound name.
 
c   itis = -id if compound
c   itis = ikp if solution 
c   itis = 0 if invalid
c----------------------------------------------------------------------
      implicit none

      integer i, itis
 
      character*10 unnown
 
      include 'perplex_parameters.h'

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c---------------------------------------------------------------------- 
 
      itis = 0

      do i = 1, isoct
         if (unnown.eq.fname(i)) then
             itis = i
             goto 99
         end if
      end do

      do i = 1, iphct
         if (unnown.eq.names(i)) then
            itis = -i
            goto 99
         end if
      end do 

99    end

      subroutine maktit 
c-----------------------------------------------------------------------
c create a title for graphics output, the title consists of the 
c calculation title + saturation hierarchy (provided one is 
c specified) and is the first two elements of title (csta8).
c if icopt = 1 or 3, also adds a blurb about reaction convention.

c title is max 3 lines, but four lines are written to be consistent
c with old plot file formats written by frendly, pt2curv etc.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      character*162 title
      common/ csta8 /title(4)

      character*8 vname,xname     
      common/ csta2  /xname(k5),vname(l2)

      integer ivfl
      common/ cst102 /ivfl

      character*5 cname
      common/ csta4 /cname(k5)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
      do i = 2, 4
         title(i) = ' '
      end do                              
c                               saturated and buffered component names:
      if (isat.gt.0) then 
         write (title(2),1070) (cname(i+icp), i= 1, isat)
      else 
         write (title(2),1000) ' '
      end if 
c                                 reaction convention
      if (icopt.eq.1.or.icopt.eq.3) write (title(3),1080) vname(ivfl)

      do i = 1, 3
         call deblnk (title(i))
      end do 

1000  format (a)
1070  format ('Component saturation hierarchy: ',7(a,1x))
1080  format ('Reaction equations are written with the high ',
     *         a,'assemblage to the right of the = sign')

      end

      subroutine rdain
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for 2d fractionation 
c calculations, called by VERTEX and WERAMI
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical dynam, titrat, qfile, err

      integer i,j,k,ier

      double precision zlayer

      character*100 name

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc 

      logical pzfunc
      integer ilay,irep,npoly,ord
      double precision abc0,vz,iblk
      common/ cst66 /abc0(0:mord,mpol),vz(6),iblk(lay,k5),ilay,
     *               irep(lay),npoly,ord,pzfunc

      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),b(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      character*100 cfname
      common/ cst227 /cfname
c-----------------------------------------------------------------------
c                                 look for input data from a file 
c                                 of type aux
      call mertxt (name,prject,'.aux',0)
      open (n8,file=name,status='old',iostat=ier)

      if (ier.ne.0) call error (51,vz(1),ilay,name)

      call mertxt (name,prject,'.fld',0)
      open (n12,file=name)
c                                 set the number of independent variables
c                                 to 1 (the independent path variable must
c                                 be variable jv(1), and the dependent path
c                                 variable must be jv(2), the path variables
c                                 can only be pressure and temperature
      ipot = 1

      jbulk = icp
c                                 true => flush model, ~true => subducting column
      read (n8,*) dynam
c                                 this is just a trick to avoid changing the variable 
c                                 names, the i/o needs to be rewired for dynam + titrat
      flsh = .not.dynam
c                                 true => basal mass flux
      read (n8,*) titrat
c                                 true => anneal the column
      read (n8,*) anneal
c                                 true => don't output console info on the fractionated phase
      read (n8,*) short
c                                 true => p-t field from file/internal function
      read (n8,*) pzfunc
c                                 Perple_X assumes upward directed depth, but to 
c                                 make the input intuitive, the input is specified
c                                 in downward coordinates, hence the sign changes 
c                                 below:
c                                 thickness of a box in column
      read (n8,*) vz(1)
c                                 gradient in variable jv(1) with z, jv(1)
c                                 is the independent variable, for subduction
c                                 this is logically pressure, i.e., dp(bar)/dz(m)
      read (n8,*) vz(2)
c                                 z if flush or or dzmax if not flush
      if (flsh) read (n8,*) vz(3)
c                                 value of the x-coordinate at the origin
      vz(4) = 0d0
      if (.not.flsh) read (n8,*) vz(4)
c                                 max value of the x-coordinate
      read (n8,*) vz(5)
c                                 Perple_X assumes an upward directed column coordinate
c                                 but in frac2d 
c                                 make the input intuitive, the input is specified
c                                 in downward coordinates, hence the sign changes 
c                                 below:

      if (flsh) then
c                                 specification n t-z points to fit n-1^th order
c                                 polynomial 
         read (n8,*) npoly

         if (npoly.gt.mpol) call error (72,b(1),i,'too many t-z '/
     *                     /'coordinates increase mpol in common cst66')

         do i = 1, npoly

            read (n8,*) b(i), a(i,1)

            do j = 2, npoly - 1
               a(i,j) = a(i,1)**j
            end do

            a(i,j) = 1d0

         end do

         call factor (a,k8,npoly,ipvt,err)

         if (.not.err) call subst (a,k8,ipvt,npoly,b,err)

         if (err) call error (72,b(1),i,'degenerate t-z'//
     *                                     ' coordinates, FRAC2D')
         do i = 1, npoly
            abc0(1,i) = b(i)
         end do

      else
c                                 now we need a path function for the dependent
c                                 variable, here we take a function defined in
c                                 terms of the absolute depth of the top of the
c                                 column (z0) and the relative depth (dz) within
c                                 the column
         if (.not.pzfunc) then
c                                 slab dip (degree)
            read (n8,*) vz(6)
c                                 number of geothermal polynomials
            read (n8,*) npoly
            if (npoly.gt.mpol) call error (72,b(1),i,'too many '/
     *                       /'geotherms increase mpol in common cst66')
c                                 order of geothermal polynomials
            read (n8,*) ord
            if (ord.gt.mord) call error (72,b(1),i,'geothermal '/
     *      /'polynomial order too high, increase mord in common cst66')

            do i = 1, npoly
c                                 depth in column for the i'th geotherm
              read (n8,*) abc0(ord+1,i) 
c                                 convert orthogonal depth to vertical depth
              abc0(ord+1,i) = abc0(ord+1,i) / 
     *                        dcos(vz(6)*.1745329252d-1)
c                                 polynomial coefficients for the geotherm
              read (n8,*) (abc0(j,i), j = 0, ord)

            end do

         end if

      end if 
c                                 get the initial global composition array
c                                 consisting of ibox compositions defined 
c                                 in terms of icp components. this read
c                                 statement assumes that H2O an CO2 (if 
c                                 thermodynamic components) are the 1st and
c                                 2nd components (if present). 
      ilay = 0
      ncol = 0
c                                 number of nodes with appended composition
c                                 end of data indicated by zero 
      do 

         read (n8,*) zlayer

         if (zlayer.eq.0) exit 

         ilay = ilay + 1

         if (ilay.eq.lay) call error (72,b(1),i, 
     *                               'increase lay in common cst66')

         read (n8,*) (iblk(ilay,i),i=1,icp)

         irep(ilay) = idint(zlayer/vz(1)) 

         ncol = ncol + irep(ilay)

         if (ncol.gt.maxbox) call error (72,b(1),i, 
     *                            'increase maxbox in common cst66')

      end do
c                                 read aliquot composition
      if (titrat) then 

         if (ilay+1.eq.lay) call error (72,b(1),i, 
     *                               'increase lay in common cst66')

         read (n8,*) qfile

         if (qfile) then
            write (*,*) 'oink'
            call errpau
         else 
            read (n8,*) (iblk(ilay+1,i),i=1,icp)
         end if 
      end if 

      close (n8)
c                                 two cases, file input or analytical
      if (fileio) then 
c                                 file input of nodal p-t coordinates
         open (n8,file=cfname,status='old',iostat=ier)
c                                 read header info
         read (n8,*) i, nrow

         if (ncol*nrow.gt.k2) call error (72,b(1),i,'too many'/
     *      /' coordinates, increase k2 to ncol*nrow in routine FRAC2D')

         if (i.ne.ncol) call error (72,b(1),i,'the number of '//
     *     'nodes in a column specified in: '//cfname//' must equal the'
     *   //' number of nodes specified in the aux file.')

         do i = 1, nrow

            k = (i-1) * ncol

            do j = 1, ncol
               read (n8,*) vn(k+j,1),vn(k+j,2)
            end do 

         end do

         close (n8)

      end if

      end

      subroutine fr2dpt (p0,dz)
c----------------------------------------------------------------------
c subroutine to set p-t variables from i-j coordinates in 2d-fractionation
c calculations, called by VERTEX and WERAMI
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical err

      integer i,j

      double precision p0, z0, dz, z2, z3, z4, z5, z6, t0, t1, t2,aa,bb

      logical pzfunc
      integer ilay,irep,npoly,ord
      double precision abc0,vz,iblk
      common/ cst66 /abc0(0:mord,mpol),vz(6),iblk(lay,k5),ilay,
     *               irep(lay),npoly,ord,pzfunc

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23  /a(k8,k8),b(k8),ipvt(k8),idv(k8),iophi,idphi,
     *                iiphi,iflg1

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
      if (fileio) then 
c                                convert p0-dz coordinate to nodal 
c                                values
         i = idint((p0 - vmn(1))/dvr(1)) + 1
         j = ncol + idint(dz/vz(1))

         v(1) = vn((i-1)*ncol + j, 1)
         v(2) = vn((i-1)*ncol + j, 2)

      else if (pzfunc) then
c                                 this could be made a lot more efficient by
c                                 making the quadratic coeffs once for each column
         z0 = p0/1d3
         z2 = z0*z0
         z3 = z2*z0
         z4 = z3*z0
         z5 = z4*z0
         z6 = z5*z0

         t2 = -0.1099312D-6*z4 +0.5065153D-4*z3 -0.3902580D-2*z2 
     *        +0.3024415D0 *z0 +0.8107985D3

          if (z0.lt.75d0) then
c                                t0, t1 shallow
             t0 = 0.1255734D-5*z5 -0.2000554D-3*z4 +0.1180485D-1*z3 
     *           -0.3163565D0 *z2 +0.6026698D1 *z0 + 0.276185544D3
             t1 = 0.1409099D-4*z4 -0.1603057D-2*z3 + 0.5553760D-1*z2 
     *           +0.2762566D0 *z0 +0.4401928241D3
          else if (z0.lt.78.99d0) then
c                                t0 deep
             t0 = -0.2059655D-9*z6 +0.2323113D-6*z5 - 0.1076535D-3*z4 
     *            +0.2625959D-1*z3 -0.3566382D1 *z2 + 0.2582593D3 *z0 
     *            -0.6916326D4
c                                t1 shallow
             t1 = 0.1409099D-4*z4 -0.1603057D-2*z3 + 0.5553760D-1*z2 
     *           +0.2762566D0 *z0 +0.4401928241D3
          else
c                                t0, t1 deep
             t0 = -0.2059655D-9*z6 +0.2323113D-6*z5 - 0.1076535D-3*z4 
     *            +0.2625959D-1*z3 -0.3566382D1 *z2 + 0.2582593D3 *z0 
     *            -0.6916326D4

             t1 = -0.3998088D-6*z4 +0.3672092D-3*z3 - 0.1290587D0*z2 
     *            +0.2181334D2 *z0 -0.5161647D3
          end if

         aa = -t1 / 272d0 + t2 / 850d0 + t0 / 400d0
         bb = -dsqrt(2d0) * (64d0*t2 - 625d0*t1 + 561d0*t0)/6800d0

         v(1) = (p0 - dz) * vz(2)

         v(2) = aa*dz**2/1d6 - bb*dz/1d3 + t0

      else if (flsh) then 

         z0 = (vz(3) -dz)
         v(1) = z0 * vz(2)
         v(2) = abc0(1,npoly)

         do i = 1, npoly-1
            v(2) = v(2) + abc0(1,i) * z0 ** i
         end do

      else
c                                 compute the npoly t-corrdinates
         do i = 1, npoly 
c                                 b - geotherm t
            b(i) = abc0(0,i)
c                                 depth for geotherm
            z0 = p0 + abc0(ord+1,i)

            do j = 1, ord
               b(i) = b(i) + abc0(j,i) * z0**j
            end do

            do j = 1, npoly-1
               a(i,j) = z0**j
            end do 

            a(i,j) = 1d0

         end do

         call factor (a,k8,npoly,ipvt,err)

         if (.not.err) call subst (a,k8,ipvt,npoly,b,err)

         if (err) call error (72,b(1),i,'degenerate t-z'//
     *                                     ' coordinates, FRAC2D')
c                                  true depth is 
         z0 = p0 - dz
c                                  pressure is
         v(1) = z0 * vz(2)
c                                  temperature is
         v(2) = b(npoly)

         do i = 1, npoly-1
            v(2) = v(2) + b(i) * z0** i
         end do

      end if

      end

      subroutine tabhed (n,vmn,dv,nv,nvar,n5name,n6name)
c----------------------------------------------------------------------
c  write header for nvar-dimensional table output
c     vmn - minimum values of the indendpendent variables
c      dv - the increments
c      nv - number of nodes
c     lop - werami option
c    nvar - number of independent variables
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nv(2), nvar, ivar, n

      double precision vmn(2), dv(2)

      character*100 n6name, n5name, vname(l3)*14

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      character vnm*8
      common/ cxt18a /vnm(l3)

      integer iam
      common/ cst4 /iam

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short
c------------------------------------------------------------------------
c                                 generate a file name and
c                                 open the file on n5
      if (iam.eq.1) then
c                                 frac2d file names: on input n5name
c                                 contains the string to be appended
c                                 to the project name
         call fopenv (n,n5name)
      else
         call fopenn (n,nvar,n5name,n6name)
      end if
c                                 initialize max/min record
      do i = 1, iprop
         prmx(i) = -1d99
         prmn(i) = 1d99
      end do
c                                 set flag for cumulative properties
c                                 generated by allprp/allmod
      first = .true.
c                                 ctr file tab format header
      write (n,'(a)') '|6.6.6'
c                                 a title
      write (n,'(a)') n5name
c                                 number of indepedent variables
      write (n,*) nvar
c                                 independent variable names,
c                                 value, increment & nodes
      do i = 1, nvar
         write (n,'(a)') vnm(i)
         write (n,*) vmn(i)
         write (n,*) dv(i)
         write (n,*) nv(i)
      end do
c                                 number of pseudo-dependent variables,
      ivar = 2
      if (icopt.eq.7.and.fileio) then
         ivar = 3
      else if (icopt.eq.9.and.iam.eq.1) then
         ivar = 1
      else
         ivar = 2
      end if
c                                 convert a8 names to a14
      do i = 1, ivar
         vname(i) = vnm(i)
         call unblnk (vname(i))
      end do
c                                 output the dependent variable counter and list
      if (kcx(1).eq.999) then
c                                  phemgp file
         write (n,*) ivar + iprop + 2
         write (n,'(200(a20,1x))') 'Name','Counter',
     *                             (vname(i), i = 1, ivar),
     *                             (dname(i),i = 1, iprop)

      else
c                                  tab file
         if (lopt(15).or.nvar.eq.1) then
c                                  with pseudo-dependent variables
            write (n,*) ivar + iprop
            write (n,'(200(a14,1x))') (vname(i), i = 1,ivar),
     *                                (dname(i), i = 1,iprop)


         else
c                                  terse format
            write (n,*) iprop
            write (n,'(200(a14,1x))') (dname(i), i = 1,iprop)

         end if

      end if

      end

      subroutine fopenn (n,dim,n5name,n6name)
c----------------------------------------------------------------------
c decide on file name and type for werami output files, open n5name on
c LUN n, n6name is opened if necessary by prphed or modhed
c dim - integer 1 = 1d, 2 = 2d, 0 = text echo
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ier, dim, n

      character*100 n5name, n6name, num*4

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      character*100 prject,tfname
      common/ cst228 /prject,tfname
c----------------------------------------------------------------------
c                                 make plot file
      do i = 1, 1000
c                                 loop to find an unused name made of
c                                 project + "i"
         write (num,'(a1,i3)') '_',i

         call unblnk (num)

         call mertxt (tfname,prject,num,0)

         if ((kop(1).eq.36.or.kop(1).eq.38).and.kcx(1).eq.999) then
c                                 phemgp format
            call mertxt (n5name,tfname,'.phm',0)

         else

            if (dim.eq.0) then
               call mertxt (n5name,tfname,'.txt',0)
            else
               call mertxt (n5name,tfname,'.tab',0)
            end if

            if (kop(1).eq.25) call mertxt (n6name,tfname,'.plt',0)

         end if

         open (n, file=n5name, status='new', iostat=ier)
c                                 presume error means a file with name
c                                 n5name already exists
         if (ier.eq.0) exit

         if (i.gt.999) call error (999,0d0,i,tfname)

      end do

      if (dim.eq.0) write (*,1000) n5name

1000  format (/,'Console output will be echoed in file: ',a,/)

      end

      subroutine fopenv (n,string)
c----------------------------------------------------------------------
c decide on file name and type for werami output files, open n5name on
c LUN n, n6name is opened if necessary by prphed or modhed
c dim - integer 1 = 1d, 2 = 2d, 0 = text echo
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ier, n

      character string*(*)

      character*100 prject,tfname
      common/ cst228 /prject,tfname
c----------------------------------------------------------------------

      call mertxt (tfname,prject,string,0)

      call mertxt (string,tfname,'.tab',0)

      open (n, file=string, status='replace', iostat=ier)

      if (ier.ne.0) call error (999,0d0,i,
     *             'file '//tfname//' is in use by another application')

      end

      subroutine factor (a,lda,n,ipvt,error)
c-----------------------------------------------------------------------
c factor is a subroutine which calculates the triangular
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
c       error- false if a is of rank = n, and true if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error

      integer ipvt(*), lda, i,j,k,ip1,n,istr,nm1

      double precision a(lda,*),d(lda),rmax,tmax,temp,ratio
c-----------------------------------------------------------------------
      error = .false.
c                            initialize ipvt,d
      do i = 1, n

         ipvt(i) = i
         rmax = 0d0

         do j = 1, n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.nopt(50)) then
            error = .true.
            return
         end if

         d(i) = rmax

      end do
c                            begin decomposition:
      nm1 = n-1
c
      do i = 1, nm1
         ip1 = i+1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/d(i)
         istr = i
         do j = ip1, n
            tmax = dabs(a(j,i))/d(j)
            if (tmax.gt.rmax) then
               rmax = tmax
               istr = j
            end if
         end do

         if (dabs(rmax).lt.nopt(50)) then
            error = .true.
            return
         end if
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then
            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = d(istr)
            d(istr) = d(i)
            d(i) = temp

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
            do k = ip1,n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do
         end do
      end do

      if (dabs(a(n,n)).lt.nopt(50)) error = .true.

      end


      subroutine subst (a,lda,ipvt,n,b,error)
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
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error

      integer lda, ipvt(*), ip, i, j, n, ii

      double precision a(lda,*), b(*), x(lda), sum
c----------------------------------------------------------------------
c                                 solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)
      do i = 2, n

         sum = 0d0

         do j = 1, i - 1
            sum = a(i,j)*x(j)+sum
         end do

         ip = ipvt(i)
         x(i) = b(ip)-sum

      end do
c                                 solve ux = y for x:
      if (a(n,n).eq.0d0) then
c                                 this check should be superfluous,
c                                 but reopt requires it. should check
c                                 what's with factor.
         error = .true.
         return
      end if

      x(n) = x(n)/a(n,n)

      do ii = 1, n - 1

         i = n-ii

         sum = 0d0

         do j = i + 1, n
            sum = a(i,j)*x(j)+sum
         end do

         if (a(i,i).eq.0d0) then
c                                 as above.
            error = .true.
            return
         end if

         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)

      end do
      b(n) = x(n)

99    end

      subroutine concrt
c---------------------------------------------------------------------
c concrt determines convergence criteria and limits for reasonable
c solutions for the routine univeq. the array delt is also used in
c the routine slope.

c input:  dv - array of default search increments.
c         vmax,vmin - arrays containing the maximum and minimum values
c           of the independent and dependent intensive variables.
c         iv - array indexes variables in vmax, vmin, and output arrays
c output: ulim, blim - arrays with upper and lower limits for reasonable
c           solutions (vmax,vmin+/-4*dv)
c         delt - array containing the finite difference increments
c           (vmax-vmin/10**5), delt also serves as the test value
c           for convergence.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision ddv

      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)
c---------------------------------------------------------------------
      do i = 1, l2
         if (dv(i).lt.0d0) call error (34,dv(i),i,'CONCRT')
         if (i.eq.3) then
            ulim(i) = vmax(i)
            blim(i) = vmin(i)
         else if (i.eq.1.or.i.eq.2) then
            ulim(i) = vmax(i)+dv(i)
            blim(i) = vmin(i)-dv(i)
            if (blim(i).lt.0d0) blim(i) = 1d0
         else
            ulim(i) = vmax(i)+dv(i)
            blim(i) = vmin(i)-dv(i)
         end if
         ddv = vmax(i)-vmin(i)
         if (ddv.lt.0d0) call error (35,ddv,i,'CONCRT')
      end do

      end

      subroutine conver (g,s,v,a,b,c,d,e,f,gg,c8,
     *                   b1,b2,b3,b4,b5,b6,b7,b8,
     *                   b9,b10,b11,b12,b13,tr,pr,r,ieos)
c---------------------------------------------------------------------
c convert thermodynamic equation of state from a pr tr reference state
c to minimize references to to reference conditions and constants.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ieos

      double precision g,s,v,a,b,c,d,e,f,gg,c8,b1,b2,b3,b4,b5,b6,b7,b8,
     *                 b9,b10,b11,b12,b13,tr,pr,n,v0,k00,k0p, dadt0,
     *                 gamma0,q0,etas0,g0,g0p,r,c1,c2, alpha0, beta0,
     *                 yr,theta,psi,eta

      double precision emodu
      common/ cst318 /emodu(k15)
c                                constants for anderson electrolyte extrapolation (ieos = 15)
      save alpha0, beta0, dadt0
      data alpha0, beta0, dadt0 /25.93d-5,45.23d-6,9.5714d-6/
c                                constants for hkf electrolyte formulation (ieos = 16)
      save psi, theta, yr, eta
      data psi, theta, yr, eta/2600d0, 228d0, -5.79865d-5, 694656.968d0/
c----------------------------------------------------------------------
c                                first conditional reformulates and returns for eos:
c                                      1, 5, 6, 11, 12, 15, 16, 101, 102
c                                reformulates and continues to second conditional for
c                                      eos < 100 and special cases (see final conditional)
      if (ieos.eq.1) then
c                                G(P,T) polynomial forms, e.g., Helgeson et al 1978 (AJS)
c                                Berman 1988 (J Pet).
         g  = g
c                                add the integral of sdt at tr,pr:
     *       + s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *       - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *       - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *       + c8 / 4d0 * tr**4
c                                add the constant components of the
c                                integral -vdp at (t,pr):
     *       - v * pr + b2 * tr * pr + b4 * pr * pr / 2d0
c                                sign on b7 corrected June 16, 2004. PJ Gorman.
     *       - b6 * pr**3/3d0 - b7 * tr * tr * pr
c                                S* (T)
         s  = a - b2*pr - s + a*dlog(tr) + b*tr - c/tr/tr/2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f/tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 /3d0
c                                b7 term added June 16, 2004. PJ Gorman.
     *        + b7 * 2d0 * pr * tr
c                                V* (P)
         v  = v - b2 * tr - b4 * pr
c                                b6 term added June 16, 2004. PJ Gorman.
     *          + b6 * pr * pr
c                                sign on b7 corrected June 16, 2004. PJ Gorman.
     *          + b7 * tr * tr
c                                a*  (-T log T)
c        a  = a
c                                b*  (-T*T)
         b  = b7 * pr + b / 2d0
c                                c*  (-1/T)
         c  = c / 2d0
c                                d*  (-T**3)
         d  = d / 6d0
c                                e*  (sqrt T)
         e  = 4d0 * e
c                                f*  (log T)
c        f  = f
c                                gg* (-1/T**2)
         gg = gg/6d0
c                                c8* T**4
         c8 = c8/12d0
c                                b2* (PT)
         b2 = b2
c                                b7 term added June 16, 2004. PJ Gorman.
     *          - 2d0 * b7 * tr
c                                b3* ((P-Pr) exp (T / c1))
c        b3 = b3
c                                b4* (P**2)
         b4 = b4/2d0
c                                b6 term added June 16, 2004. PJ Gorman.
     *          - b6 * pr
c                                b5* (exp (P/ c2))
c        b5 = c2 * b5
c                                b6* (P**3)
         b6 = b6/3d0
c                                b7* (P*T**2)
c        b7 = b7

         return

      else if (ieos.eq.5.or.ieos.eq.6) then
c                              Mie-Gruneisen Solid Models:
         if (ieos.eq.5) then
c                              stixrude & bukowinski JGR '93 +
c                              stixrude & lithgow-bertelloni 2005a (JGR)
            n = s

         else
c                              stixrude & lithgow-bertelloni GJI '05
            n = -s

         end if

         v0     = -v
         k00    = a
         k0p    = b
         gamma0 = d
         q0     = e
         etas0  = f
         g0     =  emodu(1)
         g0p    =  emodu(2)
c                                 nr9
         b1 = 9d0*n*r
c                                 c1
         b2 = 9d0*k00*v0
c                                 c2
         b3 = k0p/2d0-2d0
c                                 c3
         b4 = 3d0*b2*b3
c                                 aii
         b5 = 6d0*gamma0
c                                 aiikk
         b6 = -12d0*gamma0 + 36d0*gamma0**2 - 18d0*q0*gamma0
c                                 as
         b7 = -(gamma0 + etas0)
c                                 aiikk2
         b8 = b6/2d0
c                                 aii2
         b9 = b5/2d0
c                                 nr9t0
         b10 = b1*tr
         b11 = (3d0*k00*g0p-5d0*g0)
         b12 = ((6d0*g0p-24d0+4.5d0*k0p)*k00-14d0*g0)

         return

      else if (ieos.eq.11) then
c                                Mie-Gruneisen Stixrude liquid Model:
c                                G0 = F0
c                                S0 = S0 => S0 - Cv
c                                V0 = V0
c                                a  = Cv
c                                b  = K0 => 4.5d0*K0*V0
c                                c  = K0' => 4.5d0*K0*V0*(K'-4)
c                                d  = y0 => y0 - y'
c                                e  = y'
c                                f  = T0
c                                --- dependent ---
c                                gg = (S0-Cv-Cv*y0)*T0
c                                b1 = Cv*(y0+ln(T0))-S0+Cv
c                                b2 = ln(v0)

         gg = (s - a - a*d)*f
         b1 = a*(d+dlog(f)) - s + a
         b2 = dlog(v)

         s = s - a
         b = 4.5d0*b*v
         c = b*(c-4d0)
         d = d - e

         return

      else if (ieos.eq.12.or.ieos.eq.14.or.ieos.eq.17) then
c                                 calphad format, don't do anything.
         return

      else if (ieos.eq.15) then
c                                 H&P aqueous species, rearrange constants
c                                 and load into thermo(10-14) => (gg,b1,b2,b3,b4)
c                                 the HOH flag is in b3, move to b11
         b11 = b3

         gg = -s + tr*b + (a - b*tr)/tr/dadt0*alpha0
         b1 =  (a - b*tr)/tr/dadt0
         b2 = -b/2d0
         b3 = tr*(s - b/2d0*tr) + g - pr*v
     *        + (a-b*tr)/tr/dadt0*(beta0*pr-alpha0*tr)
         b4 = v - (a-b*tr)/tr/dadt0*beta0

         return

      else if (ieos.eq.16) then
c                                 coming in b3 is a flag (HOH) for H+, OH-, returned
c                                 as b11 (thermo(21))
        b11 = b3
c                                 DEW/HKF aqueous species
c                                 psi, theta, epsr, yr are generic parameters.
c                                 coming in HKF species parameters loaded as:
c                                 g, s, v,   a, b, c,  d,   e,  f, gg, b1, b2
c                                 (actually v and cp0 are not loaded)
c                                 and correspond to (HKF notation):
c                                 g, s, v, cp0, w,  q, a1, a2, a3, a4, c1, c2
c                                 the compound constants on output will be
c        b2 = -s + c1*dlog(tr) + c1 + w*yr + dlog(tr/(tr-theta))*c2/theta**2 => b8 in HKF_G.mws
         b3 = -s + b1*dlog(tr) + b1 + b*yr
     *            + dlog(tr/(tr-theta))*b2/theta**2
c        b3 = (-w*yr-c1+s)*tr + (-1/epsilonr+1)*w - a1*pr - a2*ln(psi+pr) + g + c2/theta => b9
c        may 2, 2017 modified to:
c        b3 = (-w*yr-c1+s)*tr + w - a1*pr - a2*ln(psi+pr) + g + c2/theta => b9
         b4 = (-b*yr-b1+s)*tr + b - d*pr
     *                        - e*dlog(psi+pr) + g + b2/theta
c        b4 = -a3*pr-a4*ln(psi+pr) => b10
         b5 = -f*pr - gg*dlog(psi+pr)
c        b5 = -c2/(tr-theta)/theta => b11
         b6 = -b2/(tr-theta)/theta
c        b6 = c2/theta^2 => b12
         b7 = b2/theta**2
c        b7 = -c1-c2/theta^2
         b8 = -(b1+b7)
c                                 the reference condition born radius (thermo 19)
         if (b.ne.0d0.or.c.ne.0d0) then 
            b9 = 5d9 * eta * c**2 / (1.622323167d9 * eta * c + 5d9 * b)
         else 
            b9 = 0d0
         end if

         return
c                                 remaining standard forms have caloric polynomial
      else if (ieos.le.202.or.ieos.eq.604.or.ieos.eq.605.or.
     *         ieos.eq.606.or.ieos.eq.700.or.ieos.eq.701.or.
     *         ieos.eq.702) then
c                                G(Pr,T) polynomial
         g  = g
     *       + s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *       - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *       - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *       + c8 / 4d0 * tr**4
         s  = a - s + a * dlog(tr) + b * tr - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 /3d0

         b  = b / 2d0
         c  = c / 2d0
         d  = d / 6d0
         e  = 4d0 * e
         gg = gg/6d0
         c8 = c8/12d0
c                                 fluid special case, this is sloppy
         if (ieos.gt.100.and.ieos.lt.120.or.ieos.eq.201.or.ieos.eq.202)
     *      return

      end if
c                                 -------------------------------------
c                                 mechanical eos:
      if (ieos.eq.3) then
c                                 for Ghiorso et al.'s PMELTS formulation
c                                 b6 (K0) is 0 and alpha is a constant to an
c                                 arbitrary reference, therefore parameters are left as is.
      else if (ieos.eq.7) then
c                                 exponential polynomial on volume, e.g., Haas et al 1986,
c                                 Gottschalk 1997. b3 = alpha, -b8 = beta:
         b1 = -v / b8 / dexp (b3*tr)

      else if (ieos.eq.8.or.ieos.eq.9) then
c                                 The HP Tait EoS:
         if (ieos.eq.8) then
c                                 HP 2010 Full Tait EoS
c                                 on input b1=alpha0, b5=theta, b6=k0, b7=k", b8=k'
c                                 thermal pressure coefficients:
c                                 alpha0/k0/xi0*theta -> b1
            b1 = 1d0/b5*b1*b6*tr**2/dexp(b5/tr)*(dexp(b5/tr)-1d0)**2
            b9 = 1d0/(dexp(b5/tr)-1d0)
c                                 Tait a parameter -> b6
            c1 = (1d0+b8)/(1d0+b8+b6*b7)
c                                 Tait b parameter -> b7
            c2 = b8/b6 - b7/(1d0+b8)
c                                 Tait c parameter -> b8 = 1-c
            b8 = 1d0 - (1d0+b8+b6*b7)/(b8**2+b8-b6*b7)
            b7 = c2
            b6 = c1
            b10 = b7*b8

         else
c                                 True tait used for melt endmembers
c                                 kt = b6 + b5*(T-Tr)
c                                 vt = v0*exp(b1*(T-Tr))
            b9 = 1d0 + b8
            b10 = b8*b9
            b11 = b7/b9
c                                 tait parameters computed as f(T)
         end if

      else if (ieos.eq.10) then
c                                 ideal gas, could make a reference pressure correction here.
      else if (ieos.eq.13) then
c                                 1) alpha = b1 + b2*T + b3/T + b4/T^2
c                                    which is reformulated here to
c                                    int(alpha,T=Tr..Tf) = b13 + b1*T + b2*T^2 + b3*ln(T) + b4/T
         b2 = b2/2d0
         b4 = -b4
         b13 = -(b1*tr + b2*tr*tr + b3*dlog(tr) + b4/tr)
c                                 2)  K' is a f(T,Tr)
      else if (b8.ne.0d0) then
c                                 All remaining forms (ieos = 2, 4, >100) assume:
c                                 1) alpha = b1 + b2*T + b3/T + b4/T^2 + b5/sqrt(T)
c                                    which is reformulated here to
c                                    int(alpha,T=Tr..Tf) = b13 + b1*T + b2*T^2 + b3*ln(T) + b4/T + b5*sqrt(T)
         b2 = b2/2d0
         b4 = -b4
         b5 = 2d0*b5
         b13 = -(b1*tr + b2*tr*tr + b3*dlog(tr) + b4/tr + b5*dsqrt(tr))
c                                 old parameter test forms:

c                                 2) if lopt(4) not true, then the isothermal bulk modulus is
c                                    K = b6 + b7*(T-Tr)
c                                    which is reformulated to
c                                    K = b6 + b7*T
         if (.not.lopt(4)) b6 = b6 - b7*tr
c                                 operation solving constants for the Murnaghan (ieos=2), these are not
c                                 used for the BM3.
         b9 = (1d0-1d0/b8)
         b10 = b8*pr
         b12 = b8 - 1d0
c                                 anderson-gruneisen parameter is assumed = K' (abs(b8)) except for
c                                 special EoS forms
         if (ieos.gt.300) then
            b11 = -s
         else
c                                 special EoS, anderson-gruneisen stored in
c                                 s-position.
            b11 = dabs(b8)
         end if

      end if

      end

      double precision function depvar (var)
c--------------------------------------------------------------------
c depvar computes the dependent variable from the independent variable
c var
c--------------------------------------------------------------------
      implicit none

      double precision var

      integer iind,idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      depvar = c0 + var*(c1 + var*(c2 + var*(c3 + var*c4)))

c      depvar = c0 + var*(c1 + var*(c2 + var*(c3 + var*(c4 + c5*var))))

      end


      subroutine makecp (inames,mnames,first)
c----------------------------------------------------------------------
c makecp reads the thermodynamic to compute the composition vector of
c made phases, called by vertex. programs without composition checking
c use smakcp.

c output to console if first = .true.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer inames, i, j, k,ict, id, incomp(k0), jct, mmeos(k16*k17)

      logical inph(k16*k17), inmk(k16), eof, good, first

      double precision mcp(k16*k17,k0)
      character name*8, mnames(k16*k17)*8

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 make a list of all definition names:
      inames = 0

      do i = 1, nmak

         inmk(i) = .true.

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do

            if (k.le.inames) cycle

            inames = inames + 1
            mnames(inames) = mknam(i,j)

         end do

      end do
c                                 array to find valid makes:
      do i = 1, inames
         inph(i) = .false.
      end do
c                                 now get the composition vectors for
c                                 the mnames phases:
      do

         call getphi (name,.false.,eof)

         if (eof) exit

         do i = 1, inames
            if (name.eq.mnames(i)) then

               do j = 1, icmpn
                  mcp(i,j) = comp(j)
               end do

               inph(i) = .true.
               mmeos(i) = ieos

               exit
            end if
         end do

      end do
c                                 find valid makes:
      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames

               if (mnames(k).eq.mknam(i,j).and.(.not.inph(k))) then

                  inmk(i) = .false.

c                 if (first) then
c                    if (iam.ne.3.and.iam.lt.6.or.iam.eq.15)
c    *                  call warn (51,tot,icmpn,mknam(i,mknum(i)+1))
c                 end if

                  exit

               else if (mnames(k).eq.mknam(i,j)) then

                  mkind(i,j) = k
                  meos(i) = mmeos(k)

               end if

            end do

            if (.not.inmk(i)) exit

         end do

      end do
c                                 compute the composition for each
c                                 made entitity and check if it's
c                                 valid
      do i = 1, nmak

         if (inmk(i)) then

            name = mknam(i,mknum(i)+1)

            do j = 1, icmpn
               mcomp(i,j) = 0d0
            end do

            do j = 1, mknum(i)
               id = mkind(i,j)
               do k = 1, icmpn
                  mcomp(i,k) = mcomp(i,k) + mkcoef(i,j)*mcp(id,k)
               end do
            end do
c                                 test the composition vector
c                                 is it a normal phase (i.e.,
c                                 non-zero thermodynamic components)
            do k = 1, icmpn
               comp(k) = mcomp(i,k)
            end do

            call chkphi (1,name,good)

            if (good) then

               mksat(i) = .false.

            else
c                                 no thermo componnents, sat comps?
               call chkphi (0,name,good)

               mksat(i) = .true.

               if (.not.good) then

                  inmk(i) = .false.

c                 call warn (52,tot,icmpn,mknam(i,mknum(i)+1))

                  cycle

               end if

            end if

         end if
      end do
c                                 clean up arrays:
      ict = 0

      do i = 1, icmpn
         incomp(i) = 0
      end do


      do i = 1, nmak

         if (inmk(i)) then

            ict = ict + 1

            mksat(ict) = mksat(i)

            mknum(ict) = mknum(i)

            meos(ict) = meos(i)

            do j = 1, mknum(ict)+1
               mknam(ict,j) = mknam(i,j)
               mkind(ict,j) = mkind(i,j)
            end do

            do j = 1, mknum(ict)
               mkcoef(ict,j) = mkcoef(i,j)
            end do

            do j = 1, k17
               mdqf(ict,j) = mdqf(i,j)
            end do

            do j = 1, icmpn
               mcomp(ict,j) = mcomp(i,j)
c                                get list of used components
               if (mcomp(ict,j).ne.0d0.and.incomp(j).eq.0) incomp(j) = j

            end do

         end if

      end do

      jct = 0

      do j = 1, icmpn

         if (incomp(j).ne.0) then
            jct = jct + 1
            incomp(jct) = incomp(j)
         end if

      end do

      nmak = ict
c                                remake list of phases required for
c                                makes:
      inames = 0

      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do

            if (k.le.inames) cycle

            inames = inames + 1
            mnames(inames) = mknam(i,j)

         end do

      end do

      if (nmak.gt.0.and.first.and.(iam.lt.3.or.iam.eq.15)) then
c                                write list of valid makes:
         write (*,1010)
         do j = 1, nmak, 6
            k = j + 5
            if (k.gt.nmak) k = nmak
            write (*,1000) (mknam(i,mknum(i)+1),i=j,k)
         end do 
         write (*,'(80(''-''))')

      end if

1000  format (4x,6(a,2x))
1010  format (/,80('-'),/,'Summary of make-definition entities:',/)

      end

      subroutine smakcp (inames,mnames)
c----------------------------------------------------------------------
c smakcp reads the thermodynamic data file to compute the composition
c vector of made phases without composition checking. vertex uses
c makecp
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer inames, i, j, k, ict, id, incomp(k0), jct

      logical inph(k16*k17), inmk(k16), eof

      double precision mcp(k16*k17,k0)

      character name*8, mnames(k16*k17)*8

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer iam
      common/ cst4 /iam

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
c----------------------------------------------------------------------
c                                 make a list of all definition names:
      inames = 0

      do i = 1, nmak

         inmk(i) = .true.

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do

            if (k.le.inames) cycle

            inames = inames + 1
            mnames(inames) = mknam(i,j)

         end do

      end do
c                                 array to find valid makes:
      do i = 1, inames
         inph(i) = .false.
      end do
c                                 now get the composition vectors for
c                                 the mnames phases:
      do

         call getphi (name,.false.,eof)

         if (eof) exit

         do i = 1, inames
            if (name.eq.mnames(i)) then

               do j = 1, icmpn
                  mcp(i,j) = comp(j)
               end do

               inph(i) = .true.

               exit
            end if
         end do

      end do
c                                 find valid makes:
      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames

               if (mnames(k).eq.mknam(i,j).and.(.not.inph(k))) then
                  inmk(i) = .false.
c                 if (iam.ne.3.and.iam.lt.6.or.iam.eq.15)
c    *               call warn (51,tot,icmpn,mknam(i,mknum(i)+1))
                  exit
               else if (mnames(k).eq.mknam(i,j)) then
                  mkind(i,j) = k
               end if

            end do

            if (.not.inmk(i)) exit

         end do

      end do
c                                 compute the composition for each
c                                 made entitity and check if it's
c                                 valid
      do i = 1, nmak

         if (inmk(i)) then

            do j = 1, icmpn
               mcomp(i,j) = 0d0
            end do

            do j = 1, mknum(i)
               id = mkind(i,j)
               do k = 1, icmpn
                  mcomp(i,k) = mcomp(i,k) + mkcoef(i,j)*mcp(id,k)
               end do
            end do
c                                 test the composition vector
c                                 is it a normal phase (i.e.,
c                                 non-zero thermodynamic components)
            do k = 1, icmpn
               comp(k) = mcomp(i,k)
            end do

         end if
      end do
c                                 clean up arrays:
      ict = 0

      do i = 1, icmpn
         incomp(i) = 0
      end do


      do i = 1, nmak

         if (inmk(i)) then

            ict = ict + 1

            mknum(ict) = mknum(i)

            do j = 1, mknum(ict)+1
               mknam(ict,j) = mknam(i,j)
            end do

            do j = 1, mknum(ict)
               mkcoef(ict,j) = mkcoef(i,j)
            end do

            do j = 1, k17
               mdqf(ict,j) = mdqf(i,j)
            end do

            do j = 1, icmpn
               mcomp(ict,j) = mcomp(i,j)
c                                get list of used components
               if (mcomp(ict,j).ne.0d0.and.incomp(j).eq.0) incomp(j) = j

            end do

         end if

      end do

      jct = 0
      do j = 1, icmpn
         if (incomp(j).ne.0) then
            jct = jct + 1
            incomp(jct) = incomp(j)
         end if
      end do

      nmak = ict

c     if (nmak.gt.0) write (*,1010) (cmpnt(incomp(j)),j=1,jct)
c                                remake list of phases required for
c                                makes:
      inames = 0

      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do

            if (k.le.inames) cycle

            inames = inames + 1
            mnames(inames) = mknam(i,j)

         end do
c                                write list of valid makes:
c        write (*,1000) mknam(i,mknum(i)+1),(mcomp(i,incomp(j)),j=1,jct)

      end do

c     write (*,1020)

c1000  format (a,1x,15(f5.2,1x))
c1010  format (/,'Summary of valid make definitions:',//,10x,15(a,1x))
c1020  format (/)

      end

      subroutine chkphi (ichk,name,good)
c-----------------------------------------------------------------------
c ichk = 0 and 2 -> test for saturated entities
c ichk = 1 and 3 -> test for non-saturated entities
c ichk = 4 -> look for phases that consist entirely of constrained components
c ichk > 1  and not 4 -> do not compare against excluded list (for make definitions).
c ichk = 5 allow phases that include only mobile/saturated components (aq species).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 name

      integer ichk,i,j

      logical good

      integer iam
      common/ cst4 /iam

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ic
      common/ cst42 /ic(k0)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ixct,ifact
      common/ cst37 /ixct,ifact

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname
c-----------------------------------------------------------------------
      good = .true.
c                               reject the data if excluded in input:
      if (ichk.lt.2.or.ichk.ge.4) then
         do i = 1, ixct
            if (name.eq.exname(i)) goto 90
         end do
      end if
c                               reject phases which have a component
c                               not included in the input:
      do i= 1, icmpn
         if (icout(i).eq.0.and.comp(i).ne.0d0) goto 90
      end do

      if (ichk.eq.5) return
c                               reject phases with negative/zero compositions
      tot = 0d0

      do j = 1, icmpn

         if (comp(j).lt.0d0.and.comp(j).gt.-nopt(50)) then
            comp(j) = 0d0
         end if

         tot = tot + comp(j)

      end do

      if (tot.eq.0d0) then

         if (iam.eq.1.or.iam.eq.15.or.iam.eq.2) then

            call warn (13,tot,j,name)

            if (lopt(56)) call wrnstp

         end if

         goto 90

      end if
c                               check for GFSM fluid species when saturated phase
c                               and saturated component constraints are in use.




c                               do a check to make sure that the phase does
c                               not consist of just mobile components
      tot = 0d0

      do j = 1, icp + isat + ifct
         tot = tot + comp(ic(j))
      end do

      if (jmct.gt.0.and.ichk.ne.4.and.tot.eq.0d0) goto 90

      if (lopt(5).and.iam.ne.5.and.iam.ne.6) then
c                               auto_exclude any phase without a pressure EoS from the
c                               thermodynamic-saturated component data space, the
c                               more general exclusion would preclude the use of such data
c                               unless it is explicitly identified as the basis for a fugacity
c                               variable.
          if (tot.ne.0d0.and.ieos.eq.0) then
c                               got a bad operator, check that it doesn't match a
c                               special component
             good = .false.

             if (lopt(7)) then

                do i = 1, ispec
                   if (name.eq.cmpnt(idspe(i))) then
                      good = .true.
                      exit
                   end if
                end do

             end if

             if (.not.good) then

                if (ichk.eq.1.and.(iam.eq.1.or.iam.eq.15)
     *                        .or.iam.eq.2.or.iam.eq.4)
     *                                 call warn (16,tot,j,name)

                goto 90

             end if

          end if

       end if
c                               the following is not executed by build:
c                               if ichk ne 0 reject phases that consist entirely
c                               of saturated components:
      if (ichk.eq.0.or.ichk.eq.2) return
c                               phases with null composition, saved if ichk = 4
c                               otherwise rejected.
      tot = 0d0

      do j = 1, icp
         tot = tot + comp(ic(j))
      end do

      if (tot.lt.0d0.and.ichk.eq.1.and.lopt(60)) then
c                               reject phases such as H2 = O2/2 - H2O
         if (iam.eq.1.or.iam.eq.2.or.iam.eq.15) then

            call warn (13,tot,j,name)

            if (lopt(56)) call wrnstp

         end if

         goto 90

      else if (tot.eq.0d0.and.ichk.ne.4) then

         goto 90

      else if (tot.ne.0d0.and.ichk.eq.4) then

         goto 90

      else if (ichk.eq.4.and.tot.eq.0d0) then
c                               reject a null phase if it contains only
c                               saturated components, since these phases
c                               are already saved in the sat list.
         tot = 0d0

         do j = icp1, icp + isat + ifct + jmct
            tot = tot + comp(ic(j))
         end do

         if (tot.eq.0d0) goto 90

      end if

      return

90    good = .false.

      end

      logical function findph (igood)
c-----------------------------------------------------------------------
c findph checks if a phase loaded by getphi consists entirely of
c component igood, if it does it returns true.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer igood, i

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
c-----------------------------------------------------------------------

      if (comp(igood).ne.0d0) then
c                                 the phase has the component,
c                                 does it have others?
         do i= 1, icmpn
            if (i.ne.igood.and.comp(i).ne.0d0) then
               findph = .false.
               return
            end if
         end do

         findph = .true.

      else

         findph = .false.

      end if

      end



      double precision function gtrans (id,j)
c-----------------------------------------------------------------------
c gtrans computes the reference pressure free energy of a compound
c aboves its jth transition.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,j

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      gtrans = therlm(12,j,id)
c                                 -sdt
     *       + t * (therlm(3,j,id) - therlm(5,j,id) * dlog(t)
     *       - t * (therlm(6,j,id) + therlm(8,j,id) * t))
     *       - (therlm(7,j,id) + therlm(11,j,id) / t) / t
     *       + therlm(9,j,id) * dsqrt(t) + therlm(10,j,id)*dlog(t)

      end

      double precision function gclpht (id,j)
c-----------------------------------------------------------------------
c gclpht computes the reference pressure free energy of a compound
c aboves its jth transition (SGTE data type)
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,j

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      gclpht = therlm(5,j,id) + therlm(6,j,id)*t +
     *         therlm(7,j,id)*t*dlog(t) + therlm(8,j,id)/t +
     *         therlm(9,j,id)/t**2 + therlm(10,j,id)/t**3 +
     *         therlm(11,j,id)/t**9 + therlm(12,j,id)*t**2
c                                 added terms for eleanor, for SGTE 1st therlm indices > 3
c                                 are incremented by 1 relative to input array tm index
c                                 parameters are otherwise unmodified. thus therlm(13...)
c                                 is tm(12...) = t12. JADC, 12/3/2017.
     *        + therlm(13,j,id)*t**3 + therlm(14,j,id)*dsqrt(t)
     *        + therlm(15,j,id)*dlog(t)

      end

      subroutine lamhel (p,t,g,vdp,ld,lct)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of helgeson et al 1978 (AJS).

c     there is something seriously wrong in this routine!!!!

c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        lct = number of transitions
c                        g = initial free energy

c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ld,lct,i,jtran
      double precision t,g,gtrans,vdp,trtp,p,dt,pstar

      external gtrans

      double precision therlm,therdi
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
      if (t.lt.therlm(1,1,ld)) return

      do i = 1, lct
         if (t.lt.therlm(1,i,ld)) then

            if (i.eq.1) then
c                                 T<T lowest transition, ignore
c                                 possibility of < 0 clapeyron slope
c                                 and exit:
                return
             else
c                                 at the i-1 th transition:
                jtran = i - 1
                exit
             end if
          else
             jtran = i
          end if
      end do

      g = gtrans(ld,jtran) + vdp
c                                 add vtrans*dp terms, this is
c                                 only set up for one transition
c                                 with a non-zero clapeyron.
c                                 should write warning.
      if (therlm(2,1,ld).eq.0d0) return

      trtp = therlm(1,1,ld) + (p-pr)/therlm(2,1,ld)
      dt = t - therlm(1,1,ld)

      if (t .gt. trtp) then
c                                 1 bar polymorph is stable at p,
c                                 p*thermo(3,id) is 0 bar standard
c                                 state pdv polymorph contribution &
c                                 therlm(4,j,ld) is Delta V of trans
c                                 contribution.
         pstar = dt*therlm(2,1,ld) + pr

         g = g + (p-pstar)*therlm(4,1,ld)

      else
c                                 the 1 bar polymorph isn't stable.
         g = g + dt * therlm(2,1,ld) * therlm(4,1,ld)

      end if

      end

      subroutine calpht (t,g,ld,lct)
c---------------------------------------------------------------------
c     calculate the extra energy of a standard cp transition.

c     input variables

c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        lct = number of transitions
c                        g = initial free energy

c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ld, lct, i, jtran
      double precision t, g, gclpht

      external gclpht

      double precision therlm,therdi
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
      if (t.lt.therlm(1,1,ld)) return

      do i = 1, lct
         if (t.lt.therlm(1,i,ld)) then

            if (i.eq.1) then
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
                return
             else
c                                 at the i-1 th transition:
                jtran = i - 1
                exit
             end if
          else
             jtran = i
          end if
      end do

c                                    SGTE data format
      g = gclpht (ld,jtran)

      end

      subroutine lamqtz (p,t,g,ld,id)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of helgeson et al 1978 (AJS). eq. (114) (corrected for the
c     misplaced parenthesis and the sign on c sub alpha), this is
C     probably bogus for coesite note that only one transition
c     is allowed.

c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        g = initial free energy

c     constants
c                Tr = 298.15 K
c                Pr = 1. bar
c                S(lope) = 38.5 bar/K
c                ba = 0.066 j/bar
c                aa = 549.82 K
c                ca = -4.973d-6 j/bar/bar
c                trt = 848. K

c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id,ld

      double precision p,t,ps,g,pdv

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision trt, tr, pr, s, aa, ba, ca , vb

      save trt, tr, pr, s, aa, ba, ca , vb

      data trt, tr, pr, s, aa, ba, ca, vb / 848., 298.15, 1d0,
     *                  38.5, 549.82, 0.066, -4.973d-6, 2.372 /

c      trtp = trt + (p-pr) / s

      if (t.gt.trt) then
         ps = (t-trt) * therlm(2,1,ld) + pr
      else
         ps = pr
      end if
c                                 if above the ref P transition T
c                                 do the cp integration:
      if (t.gt.trt) g = therlm(8,1,ld) + (p-ps) * thermo(3,id)
     *                - therlm(3,1,ld) * (t-trt)
     *     + therlm(5,1,ld) * (t - trt - t * dlog(t/trt))
     *     - (therlm(7,1,ld) + therlm(6,1,ld) * t * trt * trt)
     *     * (t - trt)**2 / 2d0 / t / trt / trt

c                                 now add in pdv terms to
c                                 the free energy, note that
c                                 pdv term for the ref polymorph
c                                 is already in:

      pdv  =  vb * (ps-pr)
     *      - ca * ( (2d0 * pr * (p-ps) - (p * p - ps * ps) ) / 2d0
     *              + s * (t-tr) * (p-ps) )
     *      + s * (ba + aa*ca*s) * (t-tr)
     *                           * dlog ((aa+p/s) / (aa + ps/s))

      g = g + pdv

      end

      subroutine lamubc (p,t,gspk,k,lct)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of berman and brown (1985, contribs. miner. petro.)

c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        k = pointer to phase in therlm

c     returned - gspk - energy because of lambda spike

c     modified from perkins et al 1986 (computers and geosciences).
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j,lct,k

      double precision gspk,ctrans,aspk2,bspk2,ct2,ct3,a1,b1,c1,t92,t93,
     *                 tr92,tr93,dhspk,dsspk,t9,tr,teq,tq1bar,p,tr9,
     *                 dvdtr,dvdp,dstr,abspk,t

      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)
c---------------------------------------------------------------------
      gspk=0d0

      do j = 1, lct

         tq1bar = therlm(3,j,k)

         if (tq1bar.eq.0d0) cycle

            tr    = therlm(7,j,k)
            teq   = therlm(4,j,k) * (p-1d0) + tq1bar
            ctrans = tq1bar - teq
            tr9    = tr - ctrans

            if (t.lt.tr9) cycle

            aspk2 = therlm(1,j,k)
            bspk2 = therlm(2,j,k)
            dvdtr = therlm(5,j,k)
            dvdp  = therlm(6,j,k)
            dstr  = therlm(8,j,k) / tq1bar
            abspk = therlm(9,j,k)

            if (t .gt. teq) then
               t9 = teq
            else
               t9 = t
            end if

            ct2 = ctrans * ctrans
            ct3 = ct2 * ctrans

            a1 = aspk2 * ctrans + 2d0 * abspk * ct2 + bspk2 * ct3
            b1 = aspk2 + 4d0 * abspk * ctrans + 3d0 * bspk2 * ct2
            c1 = 2d0 * abspk + 3d0 * ctrans * bspk2

            t92 = t9 * t9
            t93 = t9 *t92
            tr92 = tr9 * tr9
            tr93 = tr92 * tr9

            dhspk = a1 * (t9 - tr9)
     *            + b1 * (t92 - tr92) / 2d0
     *            + c1 * (t93 - tr93) / 3d0
     *            + bspk2 * (t9*t93 - tr93*tr9) / 4d0

            dsspk = a1 * (dlog(t9) - dlog(tr9))
     *            + b1 * (t9 - tr9)
     *            + c1 * (t92 - tr92) / 2d0
     *            + bspk2 * (t93 - tr93) / 3d0

            gspk = gspk - (t9 * dsspk) + dhspk

            if (t.gt.teq) gspk = gspk - (dstr + dsspk ) * (t-teq)

            gspk = gspk + dvdtr * (p-1d0) * (t9-tr)
     *                  + dvdp * ((p*p-1d0)/2d0 - (p-1d0))

      end do

      end

      subroutine disord (gval,id)
c----------------------------------------------------------------------
c compute t-dependent disorder contribution, g is the
c gibbs energy of disordering, therdi(8,id) is the t of onset of
c disordering, therdi(9,id) is the t of complete disorder.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision gval,dh,tt,ds,trr

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      trr = therdi(8,id)
      if (t .lt. trr) return

      tt  = t
      if (t.gt.therdi(9,id)) tt = therdi(9,id)

      dh = therdi(1,id) * (tt - trr)
     *     + 2d0 * therdi(2,id) * (dsqrt(tt) - dsqrt(trr))
     *     - therdi(3,id) * (1d0 / tt - 1d0 / trr)
     *     + therdi(5,id) * dlog(tt/trr)
     *     + therdi(6,id) * (tt*tt - trr*trr) / 2d0
     *     + therdi(7,id) * (tt**3 - trr**3) / 3d0

      ds = therdi(1,id) * dlog(tt/trr)
     *     - 2d0 * therdi(2,id) * (tt**(-0.5d0) - trr**(-0.5d0))
     *     - therdi(3,id) * (1d0/tt/tt - 1d0/trr/trr) / 2d0
     *     - therdi(5,id) * (1d0/tt - 1d0 / trr)
     *     + therdi(6,id) * (tt - trr)
     *     + therdi(7,id) * (tt*tt - trr*trr) / 2d0

      gval = gval + dh - (t * ds)

      if (therdi(4,id).ne.0d0) gval = gval + dh/therdi(4,id) * (p - pr)

      end

      subroutine lamla0 (dg,intvdp,ld)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using  the
c     Landau model as implemented INCORRECTLY in thermocalc pre-ds6.

c     The correct implementation of the HP98 Landau model is given
c     by function lamla1.

c     input variables

c                        t = temperature in k
c                        p = pressure in bar
c                        ld = pointer to phase in therlm
c                        intvdp = the vdp integral of the phase

c     returned - dg - difference in free energy due to the transition

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ld

      double precision dg,tc,tc0,q2,intvdp

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      tc0 = therlm(1,1,ld)
      tc = tc0 + therlm(3,1,ld)*(p-pr)

      if (t.lt.tc) then
c                                 the hp98 form is
         q2 = dsqrt(1d0-t/tc)

      else

         q2 = 0d0

      end if
c                                 See landau_d55.mws
c                                                JADC Jan 26, 2012.
      dg = therlm(2,1,ld) *
c                                 This is the ds5 PX version
     *   ( (t-tc)*q2*0.6666667d0 - therlm(8,1,ld)*t + therlm(4,1,ld) )
c                                 This is the ds5 TC version
c    * ((tc0/3d0*q2**2 + (t-tc))*q2 + therlm(7,1,ld) - t*therlm(8,1,ld))
c                                 + int(vt,p)
     *     + therlm(6,1,ld)*intvdp

      end

      subroutine lamla1 (dg,intvdp,ld)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using  the
c     Landau model (Holland and Powell '98) but as apparently
c     implemented for DS6 (i.e., HP 2010)

c     input variables

c                        t = temperature in k
c                        p = pressure in bar
c                        ld = pointer to phase in therlm
c                        intvdp = the vdp integral of the phase

c     returned - dg - difference in free energy due to the transition
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ld

      double precision dg,tc,tc0,q2,intvdp

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      tc0 =  therlm(1,1,ld)
      tc = tc0 + therlm(3,1,ld)*(p-pr)

      if (t.lt.tc) then
c                                 partially disordered:
         q2 = dsqrt((tc-t)/tc0)

      else

         q2 = 0d0

      end if
c                                 This is the hp98 version, differs
c                                 from what's in the TC code such that
c                                 dGPX - dGTC = (tc0-tc)*q2^3/3.
c                                 See landau_d60.mws
c                                                JADC Jan 26, 2012.
c      dg = therlm(2,1,ld) *
c     *  (therlm(7,1,ld) + t*(q2 - therlm(8,1,ld)) + tc*(q2**3/3d0 - q2))
c        + vdp...
c                                 TC version, according to hans vrijmoed
c                                 e-mail 9/10/12 this is correct, should
c                                 check again....
      dg = therlm(2,1,ld) *
     *   (therlm(7,1,ld) + t*(q2 - therlm(8,1,ld))
     *                   - tc*q2 + tc0*q2**3/3d0)
c                                 + int(vt,p)
     *     + therlm(6,1,ld)*intvdp

      end

      double precision function lamla2 (ld)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using the
c     Putnis Landau O/D model as implemented by Stixrude 2021.

c     in contrast to lamla0 and lamla1 the reference state is the
c     low temperature phase.

c                        ld = pointer to phase in therlm (lmda(id))
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ld

      double precision tc,tc0,q2

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      tc0 = therlm(1,1,ld)
      tc = tc0 + therlm(3,1,ld)*(p-pr)

      if (t.lt.tc) then
c                                 partially disordered:
         q2 = dsqrt((tc-t)/tc0)

      else

         q2 = 0d0

      end if

      lamla2 = therlm(2,1,ld) * (
     *         (t-tc)*(q2 - 1d0) + tc0*(q2**3 - 1d0)/3d0 )

      end

      subroutine lamla4 (dg,ld)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using
c     the Landau model as implemented by Stxrude & Lithgow-Bertelloni
c     GJI 2005, THIS IS WRONG, likely why Stx 21 flipped the reference
c     state, need to check this.

c     input variables

c                        t = temperature in k
c                        p = pressure in bar
c                        ld = pointer to phase in therlm

c     returned - dg - difference in free energy due to the transition
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ld

      double precision dg,tc,tc0,q2,vlan

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      tc0 =  therlm(1,1,ld)
      tc = tc0 + therlm(3,1,ld)*(p-pr)

      if (t.lt.tc) then
c                                 partially disordered:
         q2 = dsqrt((tc-t)/tc0)
         vlan = therlm(2,1,ld) * therlm(3,1,ld) * 
     *          ((t - tc0 - p * therlm(3,1,ld))/(tc0*q2) - q2)/2

      else

         q2 = 0d0
         vlan = 0d0

      end if

      dg = therlm(2,1,ld) *
     *     (therlm(7,1,ld) + t*(q2 - therlm(8,1,ld))
     *                     - tc*q2 + tc0*q2**3/3d0) - p * vlan

      end

      double precision function gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)
c-----------------------------------------------------------------------
c gfesi0 computes the G for BCC FeSi alloy once the speciation has been
c computed if function gfesi. See FeSiBCC.mws.

c    y  - the bulk Fe mole fraction
c    g1 - free energy of Bcc Fe, without Gmag
c    g2 - free energy of Bcc Si
c-----------------------------------------------------------------------
      implicit none

      double precision g2, g12, y, x, w0, w1, w2, xy, yx, x1, rt,
     *                 gord
c-----------------------------------------------------------------------
      yx  = 2d0*y - x
      x1  = 1d0 - x
      xy  = 1d0 - 2d0*y + x

      gfesi0  = ( dlog(x/x1*xy/yx)*x/2d0
     *          + dlog(yx/xy)*y
     *          + dlog(xy*x1)/2d0)*rt
     *      - g12*yx*x
     *      - 64d0*w2*y**4
     *      + 16d0*(8d0*w2 - w1)*y**3
     *      + 4d0 *(6d0*w1 - 20d0*w2 - w0)*y**2
     *      + 2d0 *(8d0*w2 + gord + w0 - 4d0*w1 - g2)*y + g2

      end

      double precision function gfesi1 (y,x,w0,w1,w2,rt)
c-----------------------------------------------------------------------
c gfesi0 computes the G - Gmech for BCC FeSi alloy once the speciation has been
c computed if function gfesi. See FeSiBCC.mws.

c    y  - the bulk Fe mole fraction
c    g1 - free energy of Bcc Fe, without Gmag
c    g2 - free energy of Bcc Si
c-----------------------------------------------------------------------
      implicit none

      double precision y, x, w0, w1, w2, xy, yx, x1, rt, gcon, gex
c-----------------------------------------------------------------------
      yx  = 2d0*y - x
      x1  = 1d0 - x
      xy  = 1d0 - 2d0*y + x

      gfesi1  = ( dlog(x/x1*xy/yx)*x/2d0
     *          + dlog(yx/xy)*y
     *          + dlog(xy*x1)/2d0)*rt
     *          + (((-64d0*w2*y + 128d0*w2 - 16d0*w1)*y
     *              + 24d0*w1 - 80d0*w2 - 4d0*w0)*Y + 4d0*x*w0 + 2d0*w0
     *              + 16d0*w2 - 8d0*w1)*y - 2d0*x**2*w0


      gcon = ( dlog(x/x1*xy/yx)*x/2d0
     *          + dlog(yx/xy)*y
     *          + dlog(xy*x1)/2d0)*rt


      gex  =
     *          + (((-64d0*w2*y + 128d0*w2 - 16d0*w1)*y
     *              + 24d0*w1 - 80d0*w2 - 4d0*w0)*Y + 4d0*x*w0 + 2d0*w0
     *              + 16d0*w2 - 8d0*w1)*y - 2d0*x**2*w0


      end

      subroutine dgfesi (dg,d2g,y,x,g12,rt)
c-----------------------------------------------------------------------
c dgfesi first and second derivatives of gfesi with respect to the ordered
c species concentration (x). After Lacaze & Sundman 1990, see FeSiBCC.mws.

c    y  - the bulk Fe mole fraction
c-----------------------------------------------------------------------
      implicit none

      double precision y, x, xy, yx, x1, rt, dg, d2g, g12
c-----------------------------------------------------------------------
      yx  = 2d0*y - x
      x1  = 1d0 - x
      xy  = 1d0 - 2d0*y + x

      dg  = -2d0*(y - x)*g12 + dlog(x*xy/x1/yx)*rt/2d0

      d2g =  2d0*g12 + (xy/x1/yx + x/x1/yx + x*xy/x1**2/yx
     *               + x*xy/x1/yx**2)/x/xy*x1*yx*rt/2d0

      end

      double precision function gfesic (y1,y3,y4,g1,g2,g3,g4,id)
c---------------------------------------------------------------------
c fesic4 returns the free energy change for Fe-Si-C alloy after
c Lacaze & Sundman 1990.

c    id = 30 -> BCC
c    id = 31 -> FCC

c    y1..y4 - mole fractions of Fe, Si, FeC and SiC, respectively
c    g1..g4 - free energies of Fe, Si, FeC and SiC, respectively
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision g1, g2, g3, g4, y1, y3, y4, gmech,
     *                 logu, logx, gconf, gex, x, y, u, v, gmag

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
c     x is the site fraction of Fe on the first site, y is the site
c     fraction of Si on the first site, u is the site fraction of C
c     on the second site, v is the site fraction of vacancies on the
c     second site

      x = y1 + y3
      u = y3 + y4
      y = 1d0 - x
      v = 1d0 - u

      gmech = x * v * g1 + y * v * g2 + x * u * g3 + y * u * g4

      if (x.gt.0d0.and.x.lt.1d0) then

         logx = (x * dlog(x) + y * dlog(y))

      else

         logx = 0d0

      end if

      if (u.gt.0.and.u.lt.1d0) then

         logu = (u * dlog(u) + v * dlog(v))

      else

         logu = 0d0

      end if

      if (id.eq.30) then
c                                 BCC
         gconf = r*t*(logx + 3d0*logu)

         gex = x*y*v*(-0.153138560d6 + 0.4648d2 * t - 0.92352d5 * x +
     *         0.92352d5*y + 0.62240d5*(x-y)**2) + 0.78866d5 *x*y*u
     *         - 0.190d3*x*u*v*t + gmag(1d0)

      else if (id.eq.31) then
c                                 FCC
         gconf = r * t * (logx + logu)

         gex = x*y*v* (-0.1252477d6 + 0.41116d2 * t - 0.1427076d6 * x
     *         + 0.1427076D6 * y + 0.899073D5 * (x - y) ** 2) +
     *         x * y * u * (0.1432199d6 + 0.3931d2 * t - 0.2163205d6 * x
     *         + 0.2163205D6 * y) - 0.34671d5 * x * u * v

      end if

      gfesic = gmech + gconf + gex

      end

      double precision function gfecr1 (y,g1,g2)
c-----------------------------------------------------------------------
c     gfecr1 returns the free energy change for BCC Fe-Cr alloy after
c     Andersson and Sundman, 1987, updated by Xiong et al. 2011
c     (solution model id = 32)
c     The only reason this model needs to be built in is due to the
c     continuous change in magnetic Tc and B through the FM-AFM transition.
c     Otherwise it could be a normal solution model

c    y      - mole fractions of Fe
c    g1, g2 - free energies of Fe-bcc and Cr-bcc

c                                      G. Helffrich, ELSI, 8 Apr. 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision g1, g2, y, gmech,
     *                 gconf, gex, gmag2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      gmech = y * g1 + (1-y) * g2

      if (y.lt.1d0.and.y.gt.0d0) then

          gconf = r * t * (y*dlog(y) + (1-y)*dlog(1-y))
      else
          gconf = 0d0

      end if

c     gex below is from Xiong et al. (2011) and uses an alternative magnetic
c     model than is implemented in gmag2/gmags.
      gex = (1 - y) * y * (0.2421206D5 - 0.15507D2 * t +
     *      (1 - 2 * y) * (0.166469D4 + 0.286D0 * t) +
     *      ((1 - 2 * y) ** 2) * (-0.1325088D5 + 0.8252D1 * t))
c     gex below is from Andersson & Sundman (1987) whose magnetic model is
c     consistent with gmag2/gmags.  This doesn't work due to unknown errors
c     in the Andersson & Sundman (1987) data listed in the article.
c     gex = y*(1d0-y)*(20 500d0 - 9.68d0*t)

      gfecr1 = gmech + gconf + gex + gmag2(y)

      end

      double precision function gmag2 (x)
c-----------------------------------------------------------------------
c gmag2 returns the magnetic contribution to G for BCC Fe in FeCr alloy
c after Andersson and Sundman (1987).
c     x - bulk mole fraction of Fe

c                                      G. Helffrich, ELSI, 8 Apr. 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision b,tc,x,gmags,xfe,xcr

      if (x.eq.0d0) then

         gmag2 = 0d0

      else

c        tc = 0.13545D4 * x - 0.3115D3 + (1d0 - x) * x *
c    *        (0.2200D4 - 0.11D4 * x)
c        b = 0.2228D1 * x - 0.8D-2 - 0.85D0 * (1d0 - x) * x

         xfe = x
         xcr = 1d0-x
         tc = 1043d0 * xfe + (-311.5d0) * xcr +
     *        xfe * xcr * (1650d0 + 550d0*(xcr-xfe))
         b = 2.22d0 * xfe + (-0.008d0) * xcr + xfe*xcr*(-0.008d0)
         gmag2 = gmags (tc,b,0.40d0)

      end if

      end

      double precision function gmet (id)
c----------------------------------------------------------------
c function reads SGTE data format for reference Gibbs free energy
c and evaluates thermal and pressure EoS
c
c polynomial for Gibbs free energy function
c g = a + b*T + c*T*lnT + d/T + e/T**2 + f/T**3 + g/T**9 +
c         h*T**2 + i*T**3 + j*T**4 + k*T**7
c
c EOS after Brosh et al., 2007, 2008:
c v0 volume at pr,tr; nn number of atoms; gam0 Gr�neisen parameter;
c tet0 Enstein temperature; b1,dd1,b0,dd0 fitting coefficients;
c Bo bulk modulus; Bpo pressure derivative of bulk modulus
c                                 -------------------------------
c in terms of coefficient tags nastia's polynomial is

c g = c1 + c2*T + c3*T*lnT + c4/T + c5/T**2 + c6/T**3 + c7/T**9 +
c         c8*T**2 + c9*T**3 + c10*T**4 + c11*T**7

c to accomodate saxena & eriksson (2015, Fe-S) with minimal effort
c for eleanor i use the G0 and S0 tags for the sqrt(T) and ln(T)
c SGTE terms. creating a set of tags specifically for SGTE would
c be preferable, but require rewriting nastia's data (which might,
c or might not, be worthless), so the SGTE polynomial is now
c
c g = c1 + c2*T + c3*T*lnT + c4/T + c5/T**2 + c6/T**3 + c7/T**9
c        + c8*T**2 + c9*T**3 + c10*T**4 + c11*T**7 + G0*sqrt(T)
c        + S0*ln(T)
c
c G0 and S0 are loaded into thermo(31...) and thermo(32...) via
c the ic2p pointer array.

c i didn't follow the cp/s/h stuff through the contorted 'cst1'
c 'cst2' mess to see if it's ever relevant, if it is then the
c code must be modified to include the appropriate G0/S0 terms.
c if that stuff really is necessary it's probably better to clean
c this up rather than to apply more plasters.

c                                           JADC, 12/3/2017
c-----------------------------------------------------------------
      implicit none
      include 'perplex_parameters.h'

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer id, nn

      double precision a,b,c,d,e,f,g,h,i,j,k
      double precision v0,gam0,tet0,tet02,gam02,b1,dd1,b0,dd0,Bo,Bpo
      double precision gsgte,gsgte0,sr,hr,cpr
      double precision gqh,gqhr,sqhr,hqhr,cpqhr
      double precision intp,difc,gbrosh
      double precision colcom,harter
      double precision tc,tc1,t0,beta,pp,ff,gmagn,vv,cst1,cst2
c----------------------------------------------------------------------
c                             allocate polynomial coeffiecients for
c                             reference Gibbs free energy function
      a = thermo(1,id)
      b = thermo(2,id)
      c = thermo(3,id)
      d = thermo(4,id)
      e = thermo(5,id)
      f = thermo(6,id)
      g = thermo(7,id)
      h = thermo(8,id)
      i = thermo(9,id)
      j = thermo(10,id)
      k = thermo(11,id)
c                             allocate coefficients for EoS
c                             and magnetic term
      gam0 = thermo(12,id)
      nn = idint(thermo(13,id))
      tet0 = thermo(14,id)
      b0 = thermo(15,id)
      dd0 = thermo(16,id)
      b1 = thermo(17,id)
      dd1 = thermo(18,id)
      Bo = thermo(19,id)
      Bpo = thermo(20,id)
      v0 = thermo(22,id)
      tc = thermo(23,id)
      beta = thermo(24,id)
      pp = thermo(25,id)
      vv = thermo(26,id)
c                           thermodynamic properties calculated at tr=0 K
c                           cst1: hsgte(at reference T) - hqh(at ref T) -
c                                 - tr/2*(Cpsgte(at ref T) - Cpqh(at ref T)) = hsgte(at reference T)
c                           cst2: (sqh(at ref T) - ssgte(at ref T)) +
c                                 + (Cpsgte(at ref T)-Cpqh(at ref T)) = - ssgte(at ref T)
      cst1 = thermo(27,id)
      cst2 = thermo(28,id)
c                           additional Grueneisen parameter and Einstein temperature for graphite
      gam02 = thermo(29,id)
      tet02 = thermo(30,id)

c                          read SGTE data
      gsgte = a + b*t + c*t*dlog(t) + d/t + e/t**2 + f/t**3
     *                + g/t**9 + h*t**2 + i*t**3 + j*t**4 + k*t**7
     *                + thermo(31,id)*dsqrt(t) + thermo(32,id)*dlog(t)
c                          check for transitions:
      if (ltyp(id).ne.0) call calpht (t,gsgte,lmda(id),lct(id))

c                          quasi-harmonic term
      if (nn.eq.0) then
c                          quasi-harmonic term for graphite
          gqh = 1d0*r*t*dlog(1d0 - dexp(-tet0/t)) +
     *          2d0*r*t*dlog(1d0 - dexp(-tet02/t))

      else

          gqh = 3d0*nn*r*t*dlog(1d0 - dexp(-tet0/t))

      end if

c                             interpolation function
      intp = 1d0/(1D0 + b1)*(b1 + dsqrt(1D0 +
     *              2D0*b1*(1D0 + dd1)*p/Bo))*dexp(0.10D1/b1 -
     *              1d0/b1*dsqrt(1D0 + 2D0*b1*(1D0 + dd1)*p/Bo))

c                    evaluate S, Cp, and H at reference T and P
c                    at tr and pr defined in the thermodynamic data file (tr=298.15 K)
      if (cst1.eq.0d0.or.cst2.eq.0d0) then

          gsgte0 = a + b*tr + c*tr*dlog(tr) + d/tr + e/tr**2 + f/tr**3 +
     *                 g/tr**9 + h*tr**2 + i*tr**3 + j*tr**4 + k*tr**7
c                                            sr = -dg/dt at reference t
          sr = -b - c*dlog(tr) - c + d/tr**2 + 2d0*e/tr**3 +
     *              3d0*f/tr**4 + 0.9D1*g/tr**10 - 2d0*h*tr -
     *              3d0*i*tr**2 - 0.4D1*j*tr**3 - 0.7D1*k*tr**6
c                                            hr = gr+tr*sr at reference t
          hr = gsgte0 + tr*sr
c                                            cpr = t*ds/dt at reference t
          cpr = -c - 2d0*d/tr**2 - 6d0*e/tr**3 - 12d0*f/tr**4
     *          - 9d1*g/tr**10
     *          - 2d0*h*tr - 6d0*i*tr**2 - 12d0*j*tr**3 - 42d0*k*tr**6

          gqhr = 3d0*nn*r*tr*dlog(1d0 - dexp(-tet0/tr))

          sqhr = 3d0*nn*r*tet0/tr/(dexp(tet0/tr) - 1d0) -
     *                     3d0*nn*r*dlog(1d0 - dexp(-tet0/tr))

          hqhr = 3d0*nn*r*tet0/(dexp(tet0/tr) - 1d0)

          cpqhr = 3d0*nn*r*tet0**2/tr**2*dexp(-tet0/tr)/
     *                    (1d0 - dexp(-tet0/tr))**2

c                         Cp(SGTE) - Cp(QH) with low temperature correction
c                         for the case tr=298.15
          if (t.lt.tr) then

             difc = t**2/(2D0*tr)*(cpr - cpqhr)

          else

             difc = -(gsgte - hr + t*sr) + (gqh - hqhr + t*sqhr)
     *                + (t-tr/2D0)*(cpr-cpqhr)

          end if


      else

         difc = -gsgte + gqh + cst1 + t*cst2

      end if


      gbrosh = colcom(Bo,v0,Bpo,p) +
     *         harter(nn,r,t,p,tet0,tet02,Bo,b0,dd0,gam0,gam02) -
     *         gqh + difc*(1D0 - intp)


c                     magnetic contribution using Inden-Hillert-Jarl model

      if (tc.eq.0D0.or.pp.eq.0D0) then
c                     no magnetic contribution
         gmagn = 0D0

      else
c                     pressure dependence of Tc (as modelled for cementite)
         if (vv.eq.0D0) then
            tc1 = tc
         else
            tc1 = tc*dexp(vv*p)
         end if

         t0 = t/tc1

         ff = 0d0

         if (pp.eq.0.28D0) then
c                           fcc,hcp metals and cementite
c                           for fcc-Fe beta is a negative number,
c                           therefore, gmagn can't be calculated
c                           dlog(<0)

            if (t0.lt.1D0) then
                  ff = 1d0 - 0.8603387544D0/t0 -
     #                0.1744912404D0*t0**3 -
     #                0.7755166236D-2*t0**9 - 0.1744912404D-2*t0**15
            else
                  ff = -0.4269022681D-1/t0**5 -
     #                  0.1355245296D-2/t0**15 -
     #                  0.2846015121D-3/t0**25
            end if

          else if (pp.eq.0.4D0) then
c                           bcc metals
             if (t0.lt.1D0) then
                 ff = 1d0 - 0.9052993829D0/t0 -
     #                0.1530083464D0*t0**3 -
     *                0.6800370949D-2*t0**9 -
     *                0.1530083464D-2*t0**15
             else
                 ff = -0.6417312080D-1/t0**5 -
     *                 0.2037241930D-2 /t0**15 -
     *                 0.4278208053D-3/t0**25
             end if

          end if

          gmagn = r*t*dlog(beta+1D0)*ff

      end if

      gmet = gsgte + gbrosh + gmagn

      end function gmet

      double precision function colcom (Bo,v0,Bpo,p)
c------------------------------------------------------------------------
c integral for cold compression path used in EOS Brosh et al., 2007, 2008
c------------------------------------------------------------------------
      implicit none

      double precision          Bo,V0,Bpo,p
      double precision          a,I4XC,G42,G41,G4L,G4M1,G4PC,n

         n = 4d0
         a = (n-1D0)/(3.*Bpo-1d0)
         I4XC = 1d0 - a + a*(1d0 + n/a*p/Bo/3d0)**(1D0/n)

         G42 = 0.15D1*Bpo**3 - 0.6D1*Bpo**2 +
     *                         0.8D1*Bpo - 0.3555555555D1

         G41 = -9D0*Bpo**3 + 27D0*Bpo**2 -
     *                       24D0*Bpo + 0.5333333333D1

         G4L = 9D0*Bpo**3 - 18D0*Bpo**2 +
     *                       9D0*Bpo - 0.1333333333D1

         G4M1 = 3D0*Bpo**3 - 3D0*Bpo**2 +
     *                           Bpo - 0.111111111D0

         G4PC = G42*I4XC**(-2) + G41*I4XC**(-1) - G4L*dlog(I4XC) +
     *          G4M1*I4XC - G42 - G41 - G4M1

        colcom = Bo*v0*G4PC

      end function colcom

      double precision function harter (nn,r,t,p,tet0,tet02,Bo,b0,dd0,
     *                                  gam0,gam02)
c---------------------------------------------------------------------
c integrated "thermal" volume used in EOS Brosh et al., 2007, 2008
c---------------------------------------------------------------------
      implicit none

      double precision          r,t,p,tet0,tet02,Bo,b0,dd0,gam0,gam02
      double precision          b,tet,tet1,tet2,I2XT,G2XT
      integer                   m, nn

         m = 2
         b = (m-1D0)/(3d0*b0-1d0)

         I2XT = 1d0 - b + b*(1d0 + m/b*(1D0 + dd0)*
     *          P/Bo/3d0)**(1D0/m)

         G2XT = (4.5D0*b0 - 3D0)*I2XT**(-2) + (-9D0*b0+3D0)*I2XT**(-1)
     *                                                      + 4.5D0*b0

         if (nn.eq.0) then

             tet1 = tet0*dexp(gam0/(1D0+dd0)*G2XT)
             tet2 = tet02*dexp(gam02/(1D0+dd0)*G2XT)

             harter = 1D0*r*t*dlog(1D0-dexp(-tet1/t)) +
     *               2D0*r*t*dlog(1D0-dexp(-tet2/t))

         else

             tet = tet0*dexp(gam0/(1D0+dd0)*G2XT)

             harter = 3D0*nn*r*t*dlog(1D0-dexp(-tet/t))

         end if

      end function harter

      double precision function gterm2(id)
c----------------------------------------------------------------
c function evaluates 2nd part of modified EoS for liquid C after Brosh
c
c EOS includes 2 terms: G1 (gsgte + gbrosh) and G2; G2 converges to
c a constant at high p; the term is added to the eos of liquid C
c to model graphite-like behavior of liquid at lower p and diamond-
c like at high p
c
c parameters for the second part of eos: v02, b20 (B2), b2(B'2),
c tet2 (theta2), al2 (alpha2); Murnaghan EoS was used for the
c additional term
c in thermodynamic data file: c1 - v02; c2 - b20; c3 - b2; c4 - theta2;
c c5 - alpha2
c-----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer id

      double precision tet2, b2,t1, t2, b20
c----------------------------------------------------------------------
c                             allocate coefficients for 2nd term in EoS

      b2 = thermo(3,id)
      b20 = thermo(2,id)
      tet2 = thermo(4,id)

      t1  = dexp(- b2*thermo(5,id) * (t - dlog(1d0+t/tet2)*tet2) )
      t2  = 1d0-1d0/b2

      gterm2 = b20*thermo(1,id)/(b2-1d0)*((t1+b2*p/b20)**t2 - t1**t2)

      end function gterm2

      double precision function getstr (x,y,bad)
c----------------------------------------------------------------------
c compute the conformal stretching paramter value that will give resolution y
c for conformal coordinate x
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer it

      double precision x, y, stx, st, f, df, dst, st2
c----------------------------------------------------------------------
      st = y

      it = 0

      bad = .false.

      do

         st2 = st + 2d0
         stx =  (st2 / st) ** x

         f = (st * (-st2 + y) * stx + st2 * (st +
     *        y)) / (st * stx  + st2)

         df = (-(stx ** 2) * st ** 2 + 4d0 * (1d0 + st) * (x
     *        - 1d0) * stx + st2 ** 2) / (st * stx + st2) ** 2

         dst = -f/df

         if (st + dst.lt.0d0) dst = -st / 2d0

         st = st + dst

         it = it + 1

         if (it.gt.iopt(21)) then
            bad = .true.
            exit 
         end if

         if (dabs(dst).lt.1d-3*y) exit

      end do

      getstr = st

      end


      subroutine dgfes (dg,d2g,y,x,rt,g00,g01,g02,g04,g10,g20,g30)
c-----------------------------------------------------------------------
c dgfes first and second derivatives of gfes with respect to the ordered
c species concentration (x).

c The expressions for dg and d2g are unwieldy, mainly due to the
c composition-dependent coordination in the modified QC model.
c To (greatly) reduce the length of the expressions I've directly
c incorporated model values of coordination numbers, viz:
c   ZFeFe = ZSS = 6;  ZFeS = ZSFe = 2

c    y  - the bulk Fe mole fraction
c-----------------------------------------------------------------------
      implicit none

      double precision y, x, rt, pre, pre2, lgt, bod, bod2, frcs, dg,
     *                 d2g, g00, g01, g02, g04, g10, g20, g30
c-----------------------------------------------------------------------

      pre = 3d0/(32d0*(1d0 + 2d0*x)**2d0)


      lgt = -48d0*rt*dlog( (2d0 + x - 2d0*y - 4d0*x*y)/
     *                              (2d0*(1d0 + x - y - 2d0*x*y)**2d0) )
     *        + 32d0*rt*dlog( -1d0*x/(2d0*(-1d0 - x + y + 2d0*x*y)*
     *                                             (y - x + 2d0*x*y)) )
     *        - 48d0*rt*dlog( (-3d0*x + 2d0*y + 4d0*x*y)/
     *                                      (2d0*(y - x + 2*x*y)**2d0) )


      bod = 16d0*(g00 + g10 + g20 + g30 + g01*y - g10*y - 2d0*g20*y
     *             - 3d0*g30*y + g02*y**2 + g20*y**2 + 3d0*g30*y**2
     *             - g30*y**3 + g04*y**4)
     *       + 16d0*x*(g10 + 2d0*g20 + 3d0*g30 - 6d0*g02*y - 4d0*g10*y
     *           - 10d0*g20*y - 18d0*g30*y + 8d0*g02*y**2 + 8d0*g20*y**2
     *           + 27d0*g30*y**2 - 12d0*g04*y**3 - 12d0*g30*y**3
     *           + 16d0*g04*y**4 + g01*(4d0*y - 3d0))
     *       + (4d0*x**2)*(4d0*g10 + 11d0*g20 + 21d0*g30 - 16d0*g10*y
     *           - 64d0*g20*y - 153d0*g30*y  + 162d0*g04*y**2
     *           + 8d1*g20*y**2 + 324d0*g30*y**2 - 48d1*g04*y**3
     *           - 192d0*g30*y**3 + 352d0*g04*y**4
     *           + 4d0*g01*(-3d0+4d0*y) + g02*(27d0-96d0*y+8d1*y**2))
     *       + (8d0*x**3)*(7d0*g30 + 2d0*g20*(1d0 - 4d0*y)**2
     *           + 2d0*g02*(3d0 - 4d0*y)**2 - 108d0*g04*y - 66d0*g30*y
     *           + 54d1*g04*y**2 + 192d0*g30*y**2 - 864d0*g04*y**3
     *           - 16d1*g30*y**3 + 448d0*g04*y**4)
     *       + (x**4)*(-12d0*g30*(-1d0 + 4d0*y)**3
     *                     + (g04*(-3d0 + 4d0*y)**3d0)*(-15d0 + 68d0*y))
     *       + (8d0*g04*x**5)*(3d0 - 4d0*y)**4

      dg = pre*(bod + lgt)


      pre2 = 3d0/(32d0*(1d0 + 2d0*x)**3d0)

      frcs = 3d0*(3d0 + x - (7d0 + 6d0*x)*y + (4d0 + 8d0*x)*y**2)/
     *                 ((-1d0 - x + y + 2d0*x*y)*(-2 - x + (2 + 4*x)*y))
     *      + 3d0*(3d0*x - (1d0 + 1d1*x)*y + (4d0 + 8d0*x)*y**2)/
     *                   ((-x + y + 2d0*x*y)*(-3d0*x + (2d0 + 4d0*x)*y))
     *      - 2d0*(x**2 + y - 4d0*y*x**2 + (4d0*x**2 - 1d0)*y**2)/
     *                    (x*(-1d0 - x + y + 2d0*x*y)*(y - x + 2d0*x*y))

      bod2 = -16d0*(4d0*g00 + 3d0*g01 + 3d0*g10 + 2d0*g20 + g30
     *                                          + 6d0*g02*y + 2d0*g20*y
     *              + 6d0*g30*y - 4d0*g02*y**2 - 4d0*g20*y**2
     *                                                  - 15d0*g30*y**2
     *              + 12d0*g04*y**3 + 8d0*g30*y**3 - 12d0*g04*y**4)
     *       + 24d0*x*(g20*(1d0 - 4d0*y)**2 + g02*(3d0 - 4d0*y)**2
     *              + 6d0*g04*((3d0 - 4d0*y)**2)*y**2
     *                           - 3d0*g30*((1d0 - 4d0*y)**2)*(y - 1d0))
     *       + 24d0*(x**2)*(7d0*g30 + 2d0*g20*(1d0 - 4d0*y)**2
     *                                      + 2d0*g02*(3d0 - 4d0*y)**2
     *              - 108d0*g04*y - 66d0*g30*y + 54d1*g04*y**2
     *                                                + 192d0*g30*y**2
     *              - 864d0*g04*y**3 - 16d1*g30*y**3 + 448d0*g04*y**4)
     *       + 4d0*(x**3)*(8d0*g02*(3d0 - 4d0*y)**2
     *            + 8d0*((1d0 - 4d0*y)**2)*(g20 + 5d0*g30 - 11d0*g30*y)
     *              + 3d0*g04*((3d0 - 4d0*y)**2)*
     *                                   (15d0 - 104d0*y + 128d0*y**2))
     *       + 12d0*(x**4)*(3d0*g04*((4d0*y - 3d0)**3)*(12d0*y - 5d0)
     *                                      - 4d0*g30*(4d0*y - 1d0)**3)
     *       + 48d0*g04*(x**5)*(3d0 - 4d0*y)**4

      d2g = pre2*(bod2 + 16d0*(1d0 + 2d0*x)*rt*frcs - 4d0*lgt)

      end

      double precision function
     *              gfes0 (y,x,g1,g2,rt,g00,g01,g02,g04,g10,g20,g30)
c-----------------------------------------------------------------------
c Called by gfes

      implicit none

      double precision y, x, g1, g2, rt, g00, g01, g02, g04, g10,
     *                 g20, g30, gmech, tdscnf, gex
c-----------------------------------------------------------------------

      gmech = y*g2 + (1d0 - y)*g1

      tdscnf =  (rt/(2d0 + 4d0*x)) *
     *            (2d0*(1d0 + 2d0*x)*(y - 1d0)*dlog(1d0 - y) - 2d0*
     *                                          (1d0 + 2d0*x)*y*dlog(y)
     *             - 6d0*x*dlog(-(x/(2d0*(-1d0-x+y+2d0*x*y)*
     *                                                (-x+y+2d0*x*y))))
     *             + 3d0*(-2d0-x+(2d0+4d0*x)*y)*
     *                    dlog((2d0+x-2d0*(1d0+2d0*x)*y)/
     *                                      (2d0*(1d0+x-y-2d0*x*y)**2))
     *             - 3d0*(-3d0*x+(2d0+4d0*x)*y)*
     *                    dlog((-3d0*x+(2d0+4d0*x)*y)/
     *                                       (2d0*(-x+y+2d0*x*y)**2)))

      gex = ( 3d0*x/(4d0*(8d0 + 16d0*x)) ) *
     *            (16d0*g00 + 8d0*g10*(2d0 + x - 2d0*(1d0 + 2d0*x)*y)
     *             + 4d0*g20*(2d0 + x - 2d0*(1d0 + 2d0*x)*y)**2
     *                              + 8d0*g01*((2d0 + 4d0*x)*y - 3d0*x)
     *             + 4d0*g02*((2d0 + 4d0*x)*y - 3d0*x)**2
     *                               + g04*((2d0 + 4d0*x)*y - 3d0*x)**4
     *             - 2d0*g30*((2d0 + 4d0*x)*y - x - 2d0)**3)

      gfes0 = gmech - tdscnf + gex

      end


      double precision function gmet2 (id)
c----------------------------------------------------------------
c EOS of Saxena & Eriksson (2015)a on Fe-S.
c
c This is a variant of Brosh et al., 2007, 2008;
c      see Saxena & Eriksson (2015)b on the pure Fe system.
c BUT in S&E (2015)a, the quasi-harmonic terms omit N (a.p.f.u.).
c This makes no theoretical sense, but lets us reproduce the
c results of the Fe-S models.
c
c For the cold compression part, I've used only coefficient
c c4, not the full c2-c5 formulation. This should make no
c significant difference for Earth pressures (Fei & Brosh 2014,
c confirmed by ecrg in Mathematica). To implement the full version
c (Brosh et al 2007 App A), would need to provide atomic numbers.
c
c Coded by Eleanor CR Green (ecrg) April 2018.
c                                 -------------------------------
c
c function reads SGTE data format for reference Gibbs free energy
c and evaluates thermal and pressure EoS
c
c to accomodate saxena & eriksson (2015, Fe-S) with minimal effort
c for eleanor i use the G0 and S0 tags for the sqrt(T) and ln(T)
c SGTE terms. creating a set of tags specifically for SGTE would
c be preferable, but require rewriting nastia's data (which might,
c or might not, be worthless), so the SGTE polynomial is now
c
c g = c1 + c2*T + c3*T*lnT + c4/T + c5/T**2 + c6/T**3 + c7/T**9
c        + c8*T**2 + c9*T**3 + c10*T**4 + c11*T**7 + G0*sqrt(T)
c        + S0*ln(T)
c
c G0 and S0 are loaded into thermo(31...) and thermo(32...) via
c the ic2p pointer array.

c v0 volume at pr,tr; nn number of atoms; gam0 Gr�neisen parameter;
c tet0 Enstein temperature; b1,dd1,b0,dd0 fitting coefficients;
c Bo bulk modulus; Bpo pressure derivative of bulk modulus

c                                           JADC, 12/3/2017
c-----------------------------------------------------------------
      implicit none
      include 'perplex_parameters.h'

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer id, nn, cc, ndum, cqh

      double precision a,b,c,d,e,f,g,h,i,j,k,l,m
      double precision v0,gam0,tet0,b1,dd1,b0,dd0,Bo,Bpo,beta,pp,tc
      double precision x4t,gamNx,gamN1,ccc,gc,x2t,gqh,tet
      double precision ib,ifunc
      double precision gmagn,tau,ff,dterm
      double precision gsgte,gsgte0,sr,hr,cpr
      double precision gqhp0,sqhr,hqhr,cpqhr,difc
      double precision xn,gamN


c----------------------------------------------------------------------

c                             allocate coefficients for EoS
c                             and magnetic term
      gam0 = thermo(12,id)
      nn = idint(thermo(13,id))
      tet0 = thermo(14,id)
      b0 = thermo(15,id)
      dd0 = thermo(16,id)
      b1 = thermo(17,id)
      dd1 = thermo(18,id)
      Bo = thermo(19,id)
      Bpo = thermo(20,id)
      v0 = thermo(22,id)
      tc = thermo(23,id)
      beta = thermo(24,id)
      pp = thermo(25,id)

c --------------------------------

c                             cold compression term (c4 only)
      cc = 4
      x4t = xn (cc,Bo,Bpo,p)
      gamNx = gamN (cc,x4t,Bpo)
      gamN1 = gamN (cc,1d0,Bpo)
      ccc = gamNx - gamN1
      gc = Bo*v0*ccc

c -----------
c                             quasi-harmonic term
      ndum = 1    ! dummy a.p.f.u.: S&E15 use n0=1 for all end-members
      cqh = 2   ! coefft in Brosh et al (2007) Xn function
      x2t = xn (cqh,Bo/(1d0+dd0),b0,p)
      gamNx = gamN (cqh,x2t,b0)
      gamN1 = gamN (cqh,1d0,b0)
      tet = tet0 * dexp( gam0/(1d0 + dd0) * (gamNx - gamN1) )
      gqh = 3d0*ndum*r*t*dlog(1d0 - dexp(-tet/t))


c -----------
c                             interpolating function:
c                             interpolate between SGTE cp at 1atm
c                             and QH model at high P.
      ib = dsqrt(1d0 + 2d0*b1*(1d0 + dd1)*p/Bo)
      ifunc = 1d0/(1d0+b1) * (b1+ib) * dexp( (1d0-ib) / b1 )

c -----------
c                             Inden-Hillert-Jarl magnetic contribution
c                                          (I assume; S&E15 don't say)
      if (tc.eq.0D0.or.pp.eq.0D0) then
c                             no magnetic contribution
          gmagn = 0D0
      else if (tc.lt.0D0) then
          gmagn = 0D0
c               neglect magnetic contribution: small at T of interest,
c               avoids -ve Log
      else
          tau = t/tc
          dterm = 0.4604444444d0 + 0.7318935837d0*(1d0/pp - 1d0)
          if (tau.lt.1D0) then
            ff = 1d0 - (79d0/(140d0*tau*pp) + 0.9537223340*(1d0/pp-1d0)*
     *               (tau**3/6d0 + tau**9/135d0 + tau**15/600d0))/dterm
          else
            ff = -(1d-1/tau**5 + 3.1746031746d-3/tau**(15) +
     *                  6.6666666666d-4/tau**(25))/dterm
          end if
          gmagn = r*t*dlog(beta+1D0)*ff
      end if

c -----------
c                             (Gqh(t,p0) - Gsgte(t,p0)) component

c                             allocate polynomial coeffiecients for
c                             reference Gibbs free energy function
      a = thermo(1,id)
      b = thermo(2,id)
      c = thermo(3,id)
      d = thermo(4,id)
      e = thermo(5,id)
      f = thermo(6,id)
      g = thermo(7,id)
      h = thermo(8,id)
      i = thermo(9,id)
      j = thermo(10,id)
      k = thermo(11,id)
      l = thermo(31,id)
      m = thermo(32,id)

c                          read SGTE data
      gsgte = a + b*t + c*t*dlog(t) + d/t + e/t**2 + f/t**3
     *                + g/t**9 + h*t**2 + i*t**3 + j*t**4 + k*t**7
     *                + l*dsqrt(t) + m*dlog(t)
c                          check for transitions:
      if (ltyp(id).ne.0) call calpht (t,gsgte,lmda(id),lct(id))

c                    sgte g,h,s,cp values at reference T
      gsgte0 = a + b*tr + c*tr*dlog(tr) + d/tr + e/tr**2 + f/tr**3 +
     *            g/tr**9 + h*tr**2 + i*tr**3 + j*tr**4 + k*tr**7
     *           + l*dsqrt(tr) + m*dlog(tr)
      sr = -b - c*dlog(tr) - c + d/tr**2 + 2d0*e/tr**3 +
     *            3d0*f/tr**4 + 0.9D1*g/tr**10 - 2d0*h*tr -
     *            3d0*i*tr**2 - 0.4D1*j*tr**3 - 0.7D1*k*tr**6 -
     *            m/tr - 0.5D0/dsqrt(tr)
      hr = gsgte0 + tr*sr
      cpr = -c - 2d0*d/tr**2 - 6d0*e/tr**3 - 12d0*f/tr**4 - 9d1*g/tr**10
     *        - 2d0*h*tr - 6d0*i*tr**2 - 12d0*j*tr**3 - 42d0*k*tr**6
     *        + m/tr + 0.25d0/dsqrt(tr)

c                    quasi-harmonic model values
      gqhp0 = 3d0*ndum*r*t*dlog(1d0 - dexp(-tet0/t))  ! at t,p0
      sqhr = 3d0*ndum*r*tet0/tr/(dexp(tet0/tr) - 1d0) -
     *                     3d0*ndum*r*dlog(1d0 - dexp(-tet0/tr))
      hqhr = 3d0*ndum*r*tet0/(dexp(tet0/tr) - 1d0)
      cpqhr = 3d0*ndum*r*tet0**2/tr**2*dexp(-tet0/tr)/
     *                    (1d0 - dexp(-tet0/tr))**2

      if (t.lt.tr) then
             difc = t**2/(2D0*tr)*(cpr - cpqhr)
      else
             difc = -(gsgte - hr + t*sr) + (gqhp0 - hqhr + t*sqhr)
     *                + (t-tr/2D0)*(cpr-cpqhr)
      end if

c                             assemble G
      gmet2 = gc + gqh + gsgte - gqhp0 + difc*(1d0 - ifunc) + gmagn

      end function gmet2

      double precision function gamN (n,xt,Bpo)
c----------------------------------------------------------------
c     Used in the Brosh et al (2007) equation of state.

      implicit none
      integer n,k,bin(n+1)
      double precision xt,Bpo,an,kk(n+1),kr,bnk,dk
c -----------
c                             binomial coeffts for n=2->n=5
      if (n.eq.2) then
            bin = (/ 1,2,1 /)
      else if (n.eq.3) then
            bin = (/ 1,3,3,1 /)
      else if (n.eq.4) then
            bin = (/ 1,4,6,4,1 /)
      else if (n.eq.5) then
            bin = (/ 1,5,10,10,5,1 /)
      else
            write (*,*) 'rlib:gamN: illegal n'
            stop
      end if

      an = (real(n)-1d0)/(3d0*Bpo-1d0)

      do k = 0, n
          kr = real(k)
          bnk = real(bin(k+1))
          if (k.eq.3) then
                dk = -3d0*dlog(xt)
          else
                dk = xt**(3d0-kr) * kr/(kr-3d0)
          end if
          kk(k+1) = bnk * (an - 1d0)**(n-k) * dk
      end do

      gamN = 3d0 / (an**(n-1) * real(n)) * sum(kk)

      end function gamN

      double precision function xn (n,Bo,Bpo,p)
c----------------------------------------------------------------
c     Used in the Brosh et al (2007) equation of state.

      implicit none
      double precision an,Bo,Bpo,p
      integer n
c -----------
      an = (n-1d0)/(3d0*Bpo-1d0)
      xn = 1d0/(1d0 - an + an*(1d0 + n/(3*an) * p/Bo)**(1d0/real(n)))
      end function xn

      double precision function gmags (tc,b,pee)
c-----------------------------------------------------------------------
c gmags returns the magnetic contribution to G parameterized as in
c Sundman 1991 J. Ph. Equil. v12, 127-140.
c     tc - transition temperature
c     b - Bohr magneton value
c     pee - structural magnetic parameter
c
c Not described or referenced in the original ref, but in an obscure reference,
c viz. Hertzman and Sundman (1982), CALPHAD v6. 67-80, negative tc means the
c material is antiferromagnetic, and the magnetic strength will also be a
c negative number by convention.  As explained in Xiong et al. (2012)
c CALPHAD v39, 11-20, to get actual values for evaluation, one examines the
c structural parameter pee.  For bcc structure p = 0.4, while for fcc & hcp
c structures p = 0.28.  The convention is for bcc (p=0.4) antiferromagnetic
c structures, divide tc and b by -1; for fcc & hcp (p=0.28) antiferromagnetic
c structures, divide tc and b by -3.

c                                      G. Helffrich, ELSI, 8 Apr. 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision a0,a1,t5,t15,t25,f1,f0,f3,f9,f15

      parameter (a0=518d0/1125d0, a1=11692d0/15975d0)
      parameter (t5=1d0/10d0, t15=1d0/315d0, t25=1d0/1500d0)
      parameter (f1=79d0/140d0, f0=474d0/497d0, f3=1d0/6d0,
     *           f9=1d0/135d0, f15=1d0/600d0)

      double precision a,b,t0,tc,f,pee,bc

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      if (tc.lt.0d0) then
         if (pee.lt.0.4d0) then
            bc = -b/3
            t0 = -3*t/tc
         else
            bc = -b
            t0 = -t/tc
         endif
      else
         t0 = t/tc
         bc = b
      endif

      a = a0 + a1*(1d0/pee - 1d0)

      if (t0.lt.1d0) then

         f = t - (f1*tc/pee + t * f0 * (1d0/pee - 1d0) *
     *       (f3 + t0**6 * (f9 + t0**6 * f15)) * t0**3) / a

      else

         f = -t*(t5 + (t15 + t25 / t0**10) / t0**10) / t0**5 / a

      end if

      gmags = r*f*dlog(bc+1)

      end

      double precision function gmag (x)
c-----------------------------------------------------------------------
c gmag returns the magnetic contribution to G for BCC Fe in FeSi alloy.
c after Lacaze & Sundman 1990.
c     x - bulk mole fraction of Fe
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision b,t0,tc,f,x

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      if (x.eq.0d0) then
         gmag = 0d0
         return
      end if

      tc = ((-0.1008D4 * x + 0.1512D4) * x + 0.539D3) * x

      t0 = t/tc

      if (t0.lt.1d0) then

         f = 1d0 - 0.905299383D0 / t0 - (0.153008346D0
     *       + (0.680037095D-2 + 0.153008346D-2 * t0**6) * t0**6)
     *       * t0**3

      else

         f = -(0.641731208D-1 + (0.203724193D-2 + 0.42782080051D-3 /
     *        t0**10) / t0**10) / t0**5

      end if

      b = 2.22d0 * x

      gmag = r*t*f*dlog(b+1)

      end

      double precision function hserfe (t)
c-----------------------------------------------------------------------
c hserfe returns the hser(fe) function of Lacaze & Sundman 1990.
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.1811d0) then
         hserfe  = 1224.83d0 + (124.134d0 -23.514d0*dlog(t)
     *                  + (-.439752d-2-.5892691d-7*t)*t)*t + 77358.5d0/t
      else
         hserfe = -25384.451d0 + (299.31255d0 - 46d0*dlog(t))*t
     *          + 2.2960305e31/t**9
      end if

      end

      double precision function hsersi (t)
c-----------------------------------------------------------------------
c hserfe returns the hser(si) function of Lacaze & Sundman 1990.
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.1687d0) then
         hsersi = -8162.61d0 + ((137.227d0-22.8318d0*dlog(t)) +
     *                     (-.191129d-2-.355178d-8*t)*t)*t+176667d0/t
      else
         hsersi = -9457.64d0 + t*(167.272d0 - 27.196d0*dlog(t))
     *        - .420369e31/t**9
      end if

      end


      double precision function fefcc (t)
c-----------------------------------------------------------------------
c fefcc returns the gfefcc function after Andersson and Sundman, 1987
c for calculation of Fe-Cr phase diagram
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.1811d0) then
         fefcc  = -0.23757d3 + 0.132416d3 * t - 0.246643d2 * t * dlog(t)
     *            - 0.375752d-2 * t ** 2 - 0.589269d-7 * t ** 3 +
     *              0.773585d5 / t

      else
         fefcc = -0.27098266d5 + 0.30025256d3 * t - 0.46d2 * t * dlog(t)
     *           + 0.278854d32 / t ** 9

      end if

      end

      double precision function crbcc (t)
c-----------------------------------------------------------------------
c crbcc returns the gcrbcc function after Andersson and Sundman, 1987
c for calculation of Fe-Cr phase diagram
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.2180d0) then
         crbcc = -0.885193d4 + 0.15748d3 * t - 0.26908d2 * t * dlog(t) +
     *            0.189435d-2 * t ** 2 - 0.147721d-5 * t ** 3 +
     *            0.139250d6 / t

      else
         crbcc = -0.34864d5 + 0.34418d3 * t - 0.50d2 * t * dlog(t) -
     *            0.288526d33 / t ** 9

      end if


      end

      double precision function hserc (t)
c-----------------------------------------------------------------------
c hserc returns the reference Gibbs energy of C
c-----------------------------------------------------------------------
      implicit none
      double precision t


      if (t.ge.1d-2.and.t.lt.103d0) then
         hserc = -0.104914084D4 - 0.9009204D-1 * t - 0.275D-4 * t ** 3

      else if (t.ge.103d0.and.t.le.350d0) then
         hserc = -0.98825091D3 - 0.739898691D1 * t + 0.176583D1 * t *
     #           dlog(t) - 0.1706952D-1 * t ** 2

      else
         hserc = -0.17368441D5 + 0.17073D3 * t - 0.243D2 * t * dlog(t)
     #           - 0.4723D-3 * t ** 2 + 0.2562600D7 / t -
     #             0.2643D9 / t ** 2 + 0.12D11 / t ** 3

      end if
      end

      double precision function strtch (y)
c----------------------------------------------------------------------
c get the x-value from the unstretched coordinate y
c----------------------------------------------------------------------
      implicit none

      double precision y, t

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      t = bpm**(1d0-y)

      strtch = (bp1 - bm1*t)/(1d0 + t)

      end

      double precision function unstch (x)
c----------------------------------------------------------------------
c get the y-value from the stretched coordinate x
c----------------------------------------------------------------------
      implicit none

      double precision x

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      unstch = 1d0 - dlog((bp1-x)/(bm1+x))/lbpm

      end

      subroutine assptx
c---------------------------------------------------------------------
c assptx assigns the current value of v1 and v2 to the array ptx and
c increments the ordinate counter ipt2.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      ipt2 = ipt2 + 2

      if (ipt2.gt.l5) ipt2 = l5

      ptx(ipt2-1) = v(iv1)
      ptx(ipt2)   = v(iv2)

      end

      subroutine wrnstp
c---------------------------------------------------------------------
c wrnstp terminate execution on iner
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character y*1

      write (*,'(a)') 'Continue execution despite this warning (Y/N)?'

      if (lopt(56)) then
c                                 read the choice
         read (*,'(a)') y

         if (y.ne.'y'.and.y.ne.'Y') then
c                                 quit on warning
            stop

         else
c                                 blurb on override
            write (*,1000)

         end if

      else
c                                 automatically continue
         write (*,1010)

      end if

1000  format (/,'To automatically answer interactive warnings affirmat',
     *        'ively, set warn_interactive',/,'to false.',/)
1010  format (/,'**warning ver536** the preceding interactive warning ',
     *        'was automatically answered Y',/,'because warn_interacti',
     *        've has been set to F, this is often bad practice',/)

      end

