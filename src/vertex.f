
c Please do not distribute any part of this source.
 
c Copyright (c) 1987-2020 by James A. D. Connolly, Institute for Geochemistry
c & Petrology, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.

      program vertx        
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *  vertex.sep.10.1998  *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------
 
c references:
 
c       connolly jad, kerrick dm (1987) bounds:an algorithm and program
c         for calculating composition diagrams, calphad 11, 1-57.
c       connolly jad (1990) multivariable phase diagrams: an algorithm
c         based on generalized thermodynamics, am j sci 290, p 666-718.
 
c-----------------------------------------------------------------------
c file structure:
 
c i/o logical unit numbers are specified in perplex_parameters.h 
c-----------------------------------------------------------------------
c the necessity for double precision real
c variables has not been established, users interested in reducing 
c memory requirements may wish to experiment with single precision 
c variables.
c-----------------------------------------------------------------------
c functional(1) tolerance(2) constant(3) or zero (4) dependent
c     subprograms:                         (this list is not current)

c                 testit (4)      sreset (2)      asschk (4)
c itestc   (4)    loadit (1)      flipit (2)      wway   (2)
c                 solod  (1,3)    schk   (4)      concrt (2,3)
c miscib   (4)    newass (4)      pchk   (2,4)    univeq (2,3)
c balanc   (4)    uproj  (1,2,3)  grxn   (1)      gphase (1)
c mrk      (1,3)  hsmrk  (1,3)    newrap (1,3)    fug    (1,3)

c units, tolerances, and constants:

c function specific
c routines for fluid properties (mrk, trkmrk, hsmrk) require units of
c kelvins (v(2)) and bars (v(3)) and joules.
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------

c parameters are assigned in "perplex_parameter.h"

c-----------------------------------------------------------------------
      include 'perplex_parameters.h'

      logical first, err

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer iam
      common/ cst4 /iam

      save err,first
      data err,first/.false.,.true./
c----------------------------------------------------------------------- 
c                                 iam indicates the Perple_X program
c                                    iam = 1  - vertex
c                                    iam = 2  - meemum, MC_fit
c                                    iam = 3  - werami 
c                                    iam = 4  - build 
c                                    iam = 5  - frendly
c                                    iam = 6  - ctransf
c                                    iam = 7  - pssect 
c                                    iam = 8  - psvdraw
c                                    iam = 9  - actcor
c                                    iam = 10 - rewrite
c                                    iam = 11 - fluids
c                                    iam = 13 - unsplt (global)
c                                    iam = 14 - unsplt (local)
c                                    iam = 15 - convex
      iam = 1
c                                 perplexwrap.f flags
      getInput = .true.
      sWarn = .false.
c                                 read input files, etc
      call iniprp
c                                 -------------------------------------
c                                 start the total timer (30)
      if (lopt(61)) call begtim (30)

      if (.not.refine) then
c                                 two-stage calculation,
c                                 inform user of 1st stage
         write (*,1000) 'exploratory'

      else

         write (*,1000) 'auto-refine'
c                                 header info for print and graphics files
c                                 title page for print file:
         if (io3.ne.1) call outtit

      end if
c                                 do the calculation
      call docalc
c                                 output ranges etc compositions if p2yx inversion
      if (lopt(11)) call outlim
c                                 output autorefine arf file and load
c                                 rpcs into static array (routine reload)
      call outarf

      if (iopt(6).ne.2) then
c                                 quitting after exploratory stage:
c                                 close n4/n5, delete interim results,
c                                 first is a dummy.
         call interm (.true.,first)

      else
c                                 start auto-refine stage
         outprt = .true.
         first = .false.
c                                 set refine to indicate the stage
         call setau1
c                                 set grid parameters
         call setau2
c                                 suppress output to print file, why?
c        io3 = 1
c                                 close/open prt/plt/blk
         if (io3.eq.0) then
c                                 prt output file
            call mertxt (tfname,prject,'.prn',0)
            call inqopn (n3,tfname)
            call outtit

         end if
c                                 plt output file
         call mertxt (tfname,prject,'.plt',0)
         call inqopn (n4,tfname)
c                                 blk output file
         call mertxt (tfname,prject,'.blk',0)
         call inqopn (n5,tfname)

         write (*,'(80(''-''))')
         write (*,1000) 'auto-refine'
c                                 load the former dynamic compositions
c                                 into the static arrays if manual
         if (iopt(6).eq.1) call reload (refine)
c                                 repeat the calculation
         call docalc
c                                 output ranges etc compositions if p2yx inversion
         if (lopt(11)) call outlim
c                                 output arf file if auto-re-refine
         if (lopt(55)) call outarf
c                                 clean up intermediate results
         call interm (outprt,err)

      end if

      if (lopt(61)) call cumtim
c                                 end of job msg
      write (*,1020) prject
c                                 dynamic objective function calls
      write (*,*) count

1000  format ('** Starting ',a,' computational stage **',/)
1020  format (80('-'),//,'End of job: ',a,//,80('-'),/)

      end

      subroutine cumtim
c----------------------------------------------------------------------
c output cumulative time for:

c           static LP optimization (timer 13)
c           dynamic LP optimization (timer 14)
c           successive QP optimization (timer 15)
c           total time (timer 30)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer n

      double precision tt
c----------------------------------------------------------------------
c                                 the total time is in etime(30)
      call CPU_TIME(etime(30))

      call mertxt (tfname,prject,'.tim',0)

      open (993,file=tfname)

      n = 6

      tt = (times(1) + times(13) + times(14) + times(15))

      do

         write (n,1000)

         write (n,1010) 'Static G calculation ',
     *                  times(1)/60,times(1)/etime(30)*1d2
         write (n,1010) 'Dynamic G calculation',
     *                  times(2)/60.,times(2)/etime(30)*1d2
         write (n,1010) 'Static LP            ',
     *                  times(13)/60.,times(13)/etime(30)*1d2
         write (n,1010) 'Dynamic LP           ',
     *                  times(14)/60.,times(14)/etime(30)*1d2
         write (n,1010) 'Successive QP        ',
     *                  (times(15)-times(2))/60.,
     *                  (times(15)-times(2))/etime(30)*1d2
         write (n,1010) 'Total of above       ',
     *                  tt/60.,tt/etime(30)*1d2
         write (n,1010) 'Total elapsed time   ',
     *                  etime(30)/60.,1d2
         if (n.ne.6) then 

            write (n,1020)
            write (n,'(2x,a,/)') 'Other statistics:'
            write (n,1030) 'Good SLP minimizations ',rcount(4)
            write (n,1030) 'Bad SLP minimizations  ',rcount(5)
            write (n,1030) 'Good SQP minimizations ',rcount(2)
            write (n,1030) 'Bad SQP minimizations  ',rcount(3)
            write (n,1030) 'SQP G evaluations      ',count
            write (n,1020)

            exit

         end if

         n = 993

      end do

1000  format (80('-')/,5x,'Timing',20x,'min.',9x,'% of total',/)
1010  format (2x,a21,3x,g14.5,7x,f5.1)
1020  format (80('-'))
1030  format (5x,a,1x,i16)

      end 

      subroutine docalc
c----------------------------------------------------------------------
c do the exploratory or autorefine stage calculation requested by 
c vertex
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c-----------------------------------------------------------------------
c                                 initialize potentials
      call inipot
c                                 initialize the bulk
      call iniblk

      if (icopt.eq.2) then
c                                 liquidus calculation.
         call liqdus

      else if (icopt.ge.0.and.icopt.le.4.or.icopt.eq.8) then

         call error (72,0d0,0,'you must run CONVEX for this type '//
     *                        'of calculation')

      else if (icopt.eq.5) then 
c                                 optimization on a 2-d grid.
         call wav2d1

      else if (icopt.eq.7) then 
c                                 fractionation on a 1-d path.
         call frac1d

      else if (icopt.eq.12) then 
c                                 0-d fractionation/titration
         call titrat

      else if (icopt.eq.9) then 
c                                 fractionation on a 2-d path (space-time) 
         call frac2d

      else
c                                 disabled stability field calculation
         call error (32,0d0,k2,'MAIN')

      end if

      end

      subroutine frac1d
c-----------------------------------------------------------------------
c a main program template to illustrate how to call the minimization 
c procedure. the example here is 1d fractionation.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, idead, ier, it

      logical pttrue(l2), readyn

      double precision errr(k5)

      external readyn

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical fileio, flsh, anneal, verbos, siphon, colcmp, usecmp
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,verbos,siphon,
     *                usecmp, colcmp

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer inc,jpot
      common/ cst101 /inc(l2),jpot

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer fmode,ifrct,ifr
      logical gone
      common/ frct1 /fmode,ifrct,ifr(k23),gone(k5)

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1
c-----------------------------------------------------------------------

      iasct = 0
      ibulk = 0
      it = 0  

      loopy = jlow
c                                 get phases to be fractionated
      call frname 

      open (n0-1,status='scratch')
c                                 patch to initialize unused potentials
      pttrue(1:l2) = .false.

      do i = 1, ipot
         do j = 1, jpot
            if (jv(i).eq.j) pttrue(j) = .true.
         end do
      end do

      do j = 1, jpot
         if (pttrue(j)) then 
            write (*,1080) vname(j),v(j)
            if (v(j).eq.0d0) call error (72,v(j),j,
     *     'the sectioning value of '//vname(j)//' cannot be zero.')
         end if
      end do 

1080  format ('The sectioning value for ',a,' is read from the list ',
     *        'list of minimum potential',/,'values at the end of your',
     *        ' problem definition file, its value is currently: ',
     *        g14.7)
c                                 two cases fileio: input from
c                                 file, else analytical path function
      if (fileio) then 

         open (n8,file=cfname,status='old',iostat=ier)

         if (ier.ne.0) call error (6,v(1),i,cfname)

         j = 0

         do 

            read (n8,*,iostat=ier) (v(jv(i)), i = 1, ipot)
c                                 echo to scratch
            write (n0-1,*) v

            if (jmct.gt.0) call subinc

            if (ier.ne.0) exit

            j = j + 1
            it = it + 1 

            ctotal = 0d0
c                                 get total moles to compute mole fractions             
            do i = 1, icp
               ctotal = ctotal + cblk(i)
            end do

            do i = 1, icp 
               b(i) = cblk(i)/ctotal
            end do
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
            call lpopt (1,j,idead)
c                                 fractionate the composition:
            call fractr (idead,outprt,.true.,errr)

            if (it.gt.99) then 
               write (*,'(i5,a)') j,' optimizations completed...'
               it = 0
            end if 

         end do 

         close (n8)

         loopy = j 

      else
c                                 initialize sectioning variables
         call setvar 
c                                 loopx = 99 indicates interactive
c                                 mode, otherwise traverse on  
         if (loopx.eq.99) then 

            do 

               write (*,1070) (vname(jv(i)), i = 1, ipot)
               read (*,*) (v(jv(i)), i = 1, ipot)
c                                 echo to scratch
               write (n0-1,*) v

               if (v(jv(1)).eq.0d0) exit
c                                 the systems molar composition is in 
c                                 the array cblk.
               ctotal = 0d0
c                                 get total moles to compute mole fractions             
               do i = 1, icp
                  ctotal = ctotal + cblk(i)
               end do

               do i = 1, icp 
                  b(i) = cblk(i)/ctotal
               end do
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
               call lpopt (1,j,idead)

               write (*,1000)
               do i = 1, jbulk
                  write (*,1010) cname(i), cblk(i)
               end do 
c                                 the bulk composition can be modified here, 
c                                 e.g., to add back previously fractionated
c                                 phases.
               write (*,1060) 

               if (readyn()) then 
                  write (*,1020) 
                  read (*,*) (dcomp(i), i = 1, jbulk)
                  do i = 1, jbulk
                     cblk(i) = cblk(i) + dcomp(i)
                  end do 
               end if 

            end do

         else
c                                 automated mode:
            do j = 1, loopy

               ctotal = 0d0
c                                 get total moles to compute mole fractions             
               do i = 1, icp
                  ctotal = ctotal + cblk(i)
               end do

               do i = 1, icp 
                  b(i) = cblk(i)/ctotal
               end do

               call setvr0 (j,j)
c                                 echo to scratch
               write (n0-1,*) v
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
               call lpopt (1,j,idead)
c                                 fractionate the composition:
               call fractr (idead,outprt,.true.,errr)

               it = it + 1

               if (it.gt.99) then 
                  write (*,'(i5,a)') j,' optimizations completed...'
                  it = 0
               end if 
 
            end do 

         end if

      end if 
c                                 output 
      if (outprt.and.io4.eq.0) call outgrd (1,loopy,1,n4,0)
c                                 close fractionation data files
      do i = -1, ifrct
         close (n0+i)
      end do

      close (n0-1,status='delete')

1000  format (/,'Composition is now:',/)
1010  format (1x,a,1x,g14.6)
1020  format (/,'Enter molar amounts of the components to be added ',
     *        '(ordered as above:')
1050  format (a)
1060  format (/,'Modify composition (y/n)?')
1070  format (/,'Enter (zeroes to quit) ',7(a,1x))

      end

      subroutine titrat
c-----------------------------------------------------------------------
c a main program template to illustrate how to call the minimization 
c procedure. the example here is 1d fractionation.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, idead, it

      double precision iblk(k5), errr(k5)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical fileio, flsh, anneal, verbos, siphon, colcmp, usecmp
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,verbos,siphon,
     *                usecmp, colcmp

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1

      integer fmode,ifrct,ifr
      logical gone
      common/ frct1 /fmode,ifrct,ifr(k23),gone(k5)
c-----------------------------------------------------------------------

      iasct = 0
      ibulk = 0
      it = 0
      verbos = .true.

      if (iopt(36)+1.gt.l7) then
         call warn (92,nopt(1),iopt(36)+1,'TITRAT')
         iopt(36) = l7 - 1
      end if

      do i = 1, icp
c                                 the molar composition of the infiltrant
         iblk(i) = dblk(2,i)

      end do
c                                 get phases to be fractionated
      call frname 

      do j = 1, iopt(36) + 1

         it = it + 1 

         write (*,'(/,a,i3,/)') 'aliquot ',j - 1

         ctotal = 0d0
c                                 get total moles to compute mole fractions             
         do i = 1, icp
            ctotal = ctotal + cblk(i)
         end do

         do i = 1, icp 
            b(i) = cblk(i)/ctotal
         end do
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
         call lpopt (1,j,idead)
c                                 fractionate the composition:
         call fractr (idead,outprt,verbos,errr)
c                                  add the infiltrant
         do i = 1, icp 
            cblk(i) = cblk(i) + nopt(36) * iblk(i)
         end do 

         if (it.gt.99) then 
            write (*,'(i5,a)') j,' optimizations completed...'
            it = 0
         end if 

      end do 

      close (n8)
c                                 output 
      if (outprt) then 

        call outgrd (1, iopt(36) + 1, 1,n4,0)
c                                 close fractionation data files
         do i = 1, ifrct
            close (n0+i)
         end do

         close (n0)

      end if

1000  format (/,'Composition is now:',/)
1010  format (1x,a,1x,g14.6)
1020  format (/,'Enter molar amounts of the components to be added ',
     *        '(ordered as above:')
1050  format (a)
1060  format (/,'Modify composition (y/n)?')
1070  format (/,'Enter (zeroes to quit) ',7(a,1x))

      end 

      subroutine frac2d
c-----------------------------------------------------------------------
c a subprogram template to illustrate how to do 2-d (space-time) 
c fractionation. 

c modified oct 15, 2005 for more rational input and variable gradients.
c modified oct 19, 2011 for file input
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical readyn

      character n6name*100, n5name*100, tmp*8

      integer i,j,k,l,m,idead,two(2),lun,iox,itop(lay),
     *        layer(maxbox),ibot,minus, ier, nblen

      double precision gblk(maxbox,k5),cdcomp(k5,lay),vox(k5),rho,zbox,
     *                 tot,lcomp(k5,lay),cmass(k5),cfmass(k5),area,
     *                 imass(k5),errr(k5),icerr(k5),ccerr(k5),cccomp(k5)

      external readyn, nblen

      double precision x, y
      common/ cxt46 /x, y

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)  

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1

      logical pzfunc
      integer ilay,irep,npoly,ord
      double precision abc0,vz,iblk
      common/ cst66 /abc0(0:mord,mpol),vz(6),iblk(lay,k5),ilay,
     *               irep(lay),npoly,ord,pzfunc

      logical fileio, flsh, anneal, verbos, siphon, colcmp, usecmp
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,verbos,siphon,
     *                usecmp, colcmp

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, tirst
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               tirst,kfl(i11),tname

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      logical first

      save first, vox, iox

      data first/.true./
c-----------------------------------------------------------------------
c                                 initialization
      iasct = 0
      ibulk = 0

      if (jmct.gt.0) call errdbg 
     *               ('Frac2d not set up for mobile components')

      if (flsh) then 
c                                 the y coordinate increases downward in frac2d
c                                 and increases upward in flush.
         minus = -1
      else
         minus = 1
      end if 
c                                 jlow set by 1dpath keyword in perplex_option.dat
c                                 is the number of increments in the x direction.
      loopx = jlow

      if (first) then 
c                                 set first to prevent re-reading of the 
c                                 input in auto_refine
         first = .false.
c                                 get the phase to be fractionated
         call frname 
c                                 check for consistent input if fileio
         if (fileio) then 

c           if (loopx.ne.nrow) then 
c               write (*,'(2(/,a,i4,a,a))') 
c     *         '** error ** the number of columns (',nrow,
c     *         ') specified in the coordinate file must equal the',
c     *         'number of z increments (',loopx,
c     *        ')specified in the aux file.'
c              
              loopx = nrow

c          end if

         end if
c                                 work out the oxide stoichiometries for excess O 
c                                 computation
         iox = 0
         vox = 0d0

         do i = 1, icp
            if (cname(i).eq.'C') then 
               vox(i) = 1d0
            else if (cname(i).eq.'H2') then 
               vox(i) = 0.5d0
            else if (cname(i).eq.'Si') then 
               vox(i) = 1d0
            else if (cname(i).eq.'Al') then 
               vox(i) = 0.75d0
            else if (cname(i).eq.'Fe') then 
               vox(i) = 1d0
            else if (cname(i).eq.'Mg') then 
               vox(i) = 1d0
            else if (cname(i).eq.'Ca') then 
               vox(i) = 1d0
            else if (cname(i).eq.'Na') then 
               vox(i) = 0.25d0
            else if (cname(i).eq.'K') then 
               vox(i) = 0.25d0
            else if (cname(i).eq.'S2') then 
               vox(i) = 2d0
            else if (cname(i).eq.'O2') then 
               vox(i) = 0d0
               iox = i
            else 
               write (*,*) 'Unidentified component for O deficit ',
     *                     'computation: ',cname(i)
            end if 
         end do 
      else 
c                                 get the phase to be fractionated
         call frname 
c                                 NOTE if not fileio, then jlow must not change
         if (.not.fileio)  nrow = jlow

      end if
c                                 check resolution dependent dimensions
      if (loopx*ncol.gt.k2) then
         write (*,*) ' parameter k2 must be >= ncol*loopx'
         write (*,*) ' increase parameter k2 for routine DUMMY1'
         write (*,*) ' or increase box size (vz(1)) or decrease'
         write (*,*) ' number of path increments (loopx) or try'
         write (*,*) ' the large parameter version of VERTEX'
         write (*,*) ' k2 = ',k2 
         write (*,*) ' ncol * loopx = ',ncol*loopx
         stop
      end if 

      ncol = 0 

      do i = 1, ilay

         do j = 1, irep(i)

            ncol = ncol + 1

            layer(ncol) = i

            do k = 1, icp1
               gblk(ncol,k) = iblk(i,k)
            end do

            itop(i) = ncol

         end do

      end do

      if (colcmp) then

         write (*,'(a)') 'Read column compositions from file (y/n)?'

         if (readyn()) then
            
            do

               write (*,'(a)') 
     *               'Enter file name (e.g., my_project_90.cmp):'
               read (*,*) tfname

               open (n8,file=tfname,iostat=ier,status='old')

               if (ier.eq.0) then 

                  exit

               else 

                  write (*,'(a,/,a)') 
     *           'File does not exist or is locked by another process.',
     *           'Try again (y/n)?'

                 if (readyn()) then 
                    cycle
                 else
                    write (*,'(a)') 
     *              'Coninuing with initialization from *.aux file'
                    exit
                 end if

               end if

            end do

            if (ier.eq.0) then 
c                                 good to go
               if (anneal) then
                  write (*,'(/,a)')
     *            'anneal is inconsistent with colcmp, turn off'//
     *            ' anneal (y/n)?'
                  if (readyn()) anneal = .false.
               end if

               read (n8,'(g12.6,a)') x

               if (x.gt.vmn(1)) then 
                  write (*,'(/,a,f6.0,/,a,f6.0,/,a)') 
     *            'column compositions were computed at z0 = ',x,
     *            '> current (*.aux file) coordinate is ',vz(4),
     *            'continue without correcting the aux file (y/n)?'
               
                  if (.not.readyn()) call errdbg ('quitting')
               end if

               do k = 1, ncol
                  read (n8,*) i, (gblk(k,j),j=1,icp)
               end do

               close (n8)

            end if

         end if

      end if
c                                 organize the coordinate frame and thermodynamic variables
      call getvar
c                                 initialize coordinate frame and sectioning variables
      call setvar 
c                                 set up stuff for tab file output, this is the only routine other
c                                 than WERAMI that writes tab files.
      two(1) = loopx
c                                 number of variables in table
      iprop = 2*icp1
      
      do j = 1, icp
         write (dname(j),'(a14)') cname(j)
         write (dname(j+icp1),'(a14)') cname(j)(1:nblen(cname(j)))
     *                                                      //'_{cum}'
      end do

      do j = 1, icp1
         cmass(j) = 0d0
         cfmass(j) = 0d0 
         imass(j) = 0d0
         ccerr(j) = 0d0 
      end do 

      dname(icp1) = 'O2_{def}'
      dname(2*icp1) = 'O2_{cum-def}'

      lun = n0 + 100
c                                 initialization for the top of
c                                 each layer
      write (*,'(a)') 'Layer  bot_node_dz,m   top_node_dz,m  '//
     *                '  bot_dz,m  top_dz,m   thickness,m'

      ibot = 0

      do j = 1, ilay

         do i = 1, icp1
            lcomp(i,j) = 0d0 
            cdcomp(i,j) = 0d0
         end do 

         write (*,'(2x,i1,2(7x,f8.0),5x,2(2x,f8.0),5x,f7.1)') j,
     *            vmn(2) + ibot * vz(1),
     *            vmn(2) + (itop(j)-1) * vz(1),
     *            vmn(2) + (ibot-0.5d0) * vz(1),
     *            vmn(2) + (itop(j)-0.5d0) * vz(1),
     *            (itop(j) - ibot) * vz(1)

         ibot = itop(j)

         write (n5name,'(a,i1)') '_cumulative_change_layer_',j
         call tabhed (lun + j,vmn(1),dvr(1),two,1,n5name,n6name)

      end do

c                                 number of variables in table
      iprop = icp

      do j = 1, ilay
         write (n5name,'(a,i1)') '_average_comp_layer_',j
         call tabhed (lun + ilay + j,vmn(1),dvr(1),two,1,n5name,n6name)
      end do

      n5name = '_cumulative_change_column'
      call tabhed (lun + 2*ilay + 1,vmn(1),dvr(1),two,1,n5name,n6name)

      write (*,'(/)')

      if (anneal) then

         x = vmn(1)
c                                 annealing removes the fractionated phase at
c                                 the initial condition and renormalizes to 1000 g
c                                 of material.
         do k = 1, ncol
c                                 k-loop varies depth within the column (dz)
            y = vmn(2) + dfloat(k-1) * vz(1)  
c                                 set the p-t conditions
            call fr2dpt (x,y)
c                                 now put the bulk composition from the global
c                                 array into the local array and get the total
c                                 number of moles (ctotal)
            ctotal = 0d0
c                                 get total moles to compute mole fractions             
            do i = 1, icp1
               dcomp(i) = 0d0
               cblk(i) = gblk(k,i)
               if (cblk(i).lt.nopt(11)) cblk(i) = 0d0
               ctotal = ctotal + cblk(i)
            end do

            do i = 1, icp 
               b(i) = cblk(i)/ctotal
            end do

            call lpopt (1,k,idead)

            if (idead.eq.101) then 

               write (*,*) 'orca pa duddle'

            else if (idead.ne.0.and.verbos) then 
               write (*,2000) x,y,layer(k),k,j,v(1),v(2),
     *                        (cblk(i),i=1,icp)
               write (*,2010) (itop(i),i=1,ilay)
            end if 

            call fractr (idead,.false.,verbos,errr)

            do i = 1, icp 
c                                 subtract the fluid from the current composition
               gblk(k,i) = gblk(k,i) - dcomp(i)
c                                 by not applying the zero_bulk threshold to the 
c                                 global array near zero components may accumulate
c                                 to become significant
               if (gblk(k,i).lt.0d0) gblk(k,i) = 0d0

            end do
c                                 reset carbon
c           gblk(k,2) = iblk(layer(k),2)
c                                 get the mass of the intial composition and
c                                 renormalize
            tot = 0d0 

            do i = 1, icp 
c                                 local mass
               tot = tot + gblk(k,i)*atwt(i)
            end do
c
            do i = 1, icp 
c                                 reset to 1000 g.
               gblk(k,i) = gblk(k,i)*1d3/tot
               lcomp(i,layer(k)) = lcomp(i,layer(k))
     *                           + gblk(k,i)/dfloat(irep(layer(k)))
            end do

            if (verbos) write (*,'(f8.0,12(f6.3,1x))') y,
     *                                           (gblk(k,i),i=1,icp)

         end do

         if (verbos) then 

            write (*,'(/)')

            do i = 1, ilay
               write (*,'(i1,1x,12(f6.3,1x))') i,(lcomp(j,i),j=1,icp)
            end do

          end if 

         if (flsh) then 
c                                 renormalize the aliquot to 1kg
            tot = 0d0 

            do i = 1, icp 
c                                 local mass
               tot = tot + iblk(ilay+1,i)*atwt(i)
            end do

            tot  = 1d3/tot
            rho  = 1d4*vz(2)
            zbox = vz(1)
            area = 1d0/rho/zbox
            tot  = tot*area

            do i = 1, icp 
               dblk(1,i) = iblk(ilay+1,i)*tot
            end do

            if (verbos) then 
               write (*,'(/,a,12(f10.3,1x))') 'kg/m2 aliquot ',
     *                                      (dblk(1,i),i=1,icp)
            end if

         end if
c                                 end of annealling = T section.
      else if (flsh) then
c                                 not annealing and doing titrate, set the 
c                                 aliquot to whatever the user specified
         do i = 1, icp 
            dblk(1,i) = iblk(ilay+1,i)
         end do

      end if
c                                 get total mass for conservation test
      do k = 1, ncol
         do i = 1, icp1
            imass(i) = imass(i) + gblk(k,i)
         end do 
      end do
c                                 loopx is the number of steps along the subduction
c                                 path:
      do j = 1, loopx
c                                 initialize column mass for conservation test
         do i = 1, icp1
            cmass(i) = 0d0
            icerr(i) = 0d0
            cccomp(i) = 0d0
         end do 
c                                 initialize avg layer comp
         do l = 1, ilay
            do m = 1, icp1
               lcomp(m,l) = 0d0
            end do
         end do

         if (verbos) then 
            write (*,*) '##########################################'
            write (*,'(/,a,i4,a,i4/)') 'Column ',j,' of ',loopx
            write (*,*) '##########################################'
         end if 
c                                 j loop varies the pressure at the top of the 
c                                 column (p0), set p0:
         x = vmn(1) + dfloat(j-1)*dvr(1)
c                                 for each box in our column compute phase 
c                                 relations
         do k = 1, ncol
c                                 k-loop varies depth within the column (dz)
            y = vmn(2) + dfloat(k-1) * vz(1)
c                                 set the p-t conditions
            call fr2dpt (x,y)
c                                 now put the bulk composition from the global
c                                 array into the local array and get the total
c                                 number of moles (ctotal)
            ctotal = 0d0
c                                 get total moles to compute mole fractions
            do i = 1, icp1
               dcomp(i) = 0d0
               cblk(i) = gblk(k,i)
c                                 apply the zero_bulk filter only to the working 
c                                 composition, this allows near zero components to 
c                                 accumulate without destabilizing the optimization
               if (cblk(i).lt.nopt(11)) cblk(i) = 0d0

               ctotal = ctotal + cblk(i)

            end do

            do i = 1, icp 
               b(i) = cblk(i)/ctotal
            end do
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
            if (io3.eq.0) write (n3,*) 'step ',j,' on path, box ',k

            call lpopt (j,k,idead)

            if (idead.eq.101) then 

               write (*,*) 'orca pa duddle ',j,k

            else if (idead.ne.0) then 

               write (*,*) idead
c                                 write failure info to fld file:
               write (n12,2000) x,y,layer(k),k,j,v(1),v(2),
     *                          (cblk(i),i=1,icp)
               write (n12,2010) (itop(i),i=1,ilay)

               write (*,2000) x,y,layer(k),k,j,v(1),v(2),
     *                        (cblk(i),i=1,icp)
               write (*,2010) (itop(i),i=1,ilay)

            end if 

            call fractr (idead,outprt,verbos,errr)
c                                 at this point we've computed the stable
c                                 assemblage at each point in our column
c                                 and could do mass transfer, etc etc 
c                                 here we'll simply fractionate the fluid 
c                                 phase
            if (iox.ne.0) dcomp(icp1) = dcomp(iox)

            do i = 1, icp 
c                                 subtract the fluid from the current composition
               gblk(k,i) = gblk(k,i) - dcomp(i)
c                                 by not applying the zero_bulk threshold to the 
c                                 global array near zero components may accumulate
c                                 to become significant
               if (gblk(k,i).lt.0d0) gblk(k,i) = 0d0

               if (.not.siphon) then 
c                                 add the fluid to the overlying node
                  if (k.lt.ncol) gblk(k+1,i) = gblk(k+1,i) + dcomp(i)
               end if

c                                 oxygen deficit and cumulative change
               if (iox.ne.0d0) dcomp(icp1) = dcomp(icp1) 
     *                                     - vox(i)*dcomp(i)
c                                 average layer comp
               lcomp(i,layer(k)) = lcomp(i,layer(k)) + gblk(k,i)/
     *                             dfloat(irep(layer(k)))
               if (.not.siphon) then
c                                 save layer specific results
                  do l = 1, ilay
c                                 cumulative change
                     if (k.eq.itop(l)) cdcomp(i,l) =
     *                                 cdcomp(i,l) + dcomp(i)
                  end do

               else
c                                 siphoning, add all fluid to
c                                 cccomp, cdcomp
                  cccomp(i) = cccomp(i) + dcomp(i)
                  cdcomp(i,ilay) = cdcomp(i,ilay) + dcomp(i)
               end if

            end do

            if (.not.siphon) then

               do l = 1, ilay
c                                 cumulative change
                  if (k.eq.itop(l)) then 
                     cdcomp(icp1,l) = cdcomp(icp1,l) + dcomp(icp1)
                     write (lun + l,'(200(g13.6,1x))') 
     *                                         x,(dcomp(i),i=1,icp1),
     *                                         (cdcomp(i,l),i=1,icp1)
                  end if

               end do

            else

               cccomp(icp1) = cccomp(icp1) + dcomp(icp1)
               cdcomp(icp1,ilay) = cdcomp(icp1,ilay) + dcomp(icp1)
               write (lun + ilay,'(200(g13.6,1x))') 
     *                                         x,(cccomp(i),i=1,icp1),
     *                                         (cdcomp(i,ilay),i=1,icp1)

            end if

            do i = 1, icp
c                                 mass being lost from the column
               if (k.eq.ncol) cfmass(i) = cfmass(i) + dcomp(i)

               if (errr(i).gt.nopt(11)) then 
c                                 instantaneous fractionated column mass error
                  icerr(i) = icerr(i) + errr(i)
               end if 
c                                 instantaneous column mass
               cmass(i) = cmass(i) + gblk(k,i)
            end do 
c                                 end of the k index loop
         end do

         do i = 1, icp
            ccerr(i) = ccerr(i) + icerr(i)
         end do 

         if (flsh) then 
            write (*,'(/,a,f9.0)') 'Average Layer Compositions at '
     *                             //'QINT,kg/m^2 = ',x
         else
            write (*,'(/,a,f9.0)') 'Average Layer Compositions at '
     *                             //'z0,m = ',x
         end if 

         write (*,'(5x,12(3x,a,3x))') (cname(i),i=1,icp)

         do l = ilay, 1, -1
c                                 average layer compositions
            write (lun + ilay + l,'(20(g13.6,1x))') 
     *                                       x,(lcomp(i,l),i=1,icp)
            write (*,'(i1,1x,12(f10.5,1x))') l,(lcomp(i,l),i=1,icp)

         end do

c                                 conservation tests:
         write (*,'(/,a)') 'Instantaneous and cumulative erro'
     *   //'r in column molar mass'
         write (*,'(2x,12(f10.5,1x))') (icerr(i),i=1,icp)
         write (*,'(2x,12(f10.5,1x))') (ccerr(i),i=1,icp)

         write (*,'(/,a)') 'Cumulative molar mass-loss by fractionation'
         if (.not.siphon) then
            write (*,'(2x,12(f10.5,1x))') (cfmass(i),i=1,icp)
         else
            write (*,'(2x,12(f10.5,1x))') (cdcomp(i,ilay),i=1,icp)
         end if

         if (flsh) then 

            write (*,'(/,a)') 'Cumulative molar mass-gain by titration'
            write (*,'(2x,12(f10.5,1x))') (x*dblk(1,i),i=1,icp)

            write (*,'(/,a)') 'Within-column net gain in molar mass'
            write (*,'(2x,12(f10.5,1x))') (x*dblk(1,i)-cfmass(i),
     *                i=1,icp)

            write (lun + 2*ilay + 1,'(20(g13.6,1x))') 
     *                                       x,(x*dblk(1,i)-cfmass(i),
     *                i=1,icp)
         else

            write (lun + 2*ilay + 1,'(20(g13.6,1x))') 
     *                                       x,(-cfmass(i),i=1,icp)

         end if 

         write (*,'(/,a)') 'Instantaneous column molar '
     *   //'mass '
         write (*,'(2x,12(f10.5,1x))') (cmass(i),i=1,icp)

         do i = 1, icp 
            errr(i) = cmass(i)+cfmass(i)
            if (flsh) errr(i) = errr(i) - x*dblk(1,i)
         end do 

         write (*,'(/,a)') 'Instantaneous column molar '
     *                   //'mass + mass-loss - mass-gain'
         write (*,'(2x,12(f10.5,1x))') (errr(i),i=1,icp)

         write (*,'(/,a)') 'Initial column molar mass'
         write (*,'(2x,12(f10.5,1x))') (imass(i),i=1,icp)

         write (*,'(/,a)') 'Molar mass imbalance'
         write (*,'(2x,12(f10.0,1x))') (imass(i)-errr(i),i=1,icp)

         write (*,'(/)')
c                                 add the aliquot to the base of the column
         if (flsh) then 
            do i = 1, icp 
               gblk(1,i) = gblk(1,i) + dvr(1) * dblk(1,i)
            end do
         end if 
c                                 end of j index loop
      end do 

      do j = 1, 2*ilay + 1
         close (lun + j)
      end do

      if (outprt) call outgrd (loopx,ncol,1,n4,0)

      if (colcmp) then 
c                                 dump column compositions
         write (tmp,'(''_'',i3,''.cmp'')') idint(x/1000)
         call unblnk (tmp)
         call mertxt (tfname,prject,tmp,0)

         open (n8,file=tfname)

         write (n8,'(g12.6,a)') x,' <- final Z0 coordinate'

         do k = 1, ncol
            write (n8,'(i4,1x,15(g12.6,1x))') k, (gblk(k,j),j=1,icp)
         end do

         close (n8)

      end if

      write (*,'(/,a)') 'NOTE: use nodal coordinates for layer'//
     *                    ' boundaries in PSSECT and WERAMI.'
      write (*,'(/,a)') 'Layer  bot_node_dz,m   top_node_dz,m  '//
     *                  '  bot_dz,m   top_dz,m  thickness,m'

      ibot = 0

      do j = 1, ilay

         write (*,'(2x,i1,2(7x,f8.0),5x,2(2x,f8.0),5x,f7.1)') j,
     *            vmn(2) + ibot * vz(1),
     *            vmn(2) + (itop(j)-1) * vz(1),
     *            vmn(2) + (ibot-0.5d0) * vz(1),
     *            vmn(2) + (itop(j)-0.5d0) * vz(1),
     *            (itop(j) - ibot) * vz(1)
         ibot = itop(j)

      end do

2000  format (/,' failed at z0-dz = ',2(g14.7,1x),' layer ',i1,' node '
     *       ,i3,' column ',i3,/,' p-t-c ',2(g14.7,1x),/,12(g14.7,1x))
2010  format (' layer boundaries are ',5(i3,1x))
2020  format (/,' did not fail at z0-dz = ',2(g14.7,1x),' l ',i1,' n '
     *       ,i3,' column ',i3,/,' p-t-c ',2(g14.7,1x),/,12(g14.7,1x))

      end

      subroutine lpopt (i,j,idead)
c-----------------------------------------------------------------------
c lpopt - is a shell to maintain backwards compatibility used for gridded
c minimization calculations. it calls lpopt0 to do a minimization, then 
c does some minor bookkeeping for node i,j. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, idead
c-----------------------------------------------------------------------
c                                 check for positive bulk, this
c                                 is for icont = 2 with closed
c                                 compositions space (lopt(1) = T), george's
c                                 indexing should eliminate the 
c                                 possibility.
      call chkblk (idead)

      if (idead.eq.0) call lpopt0 (idead)
c                                 if idead = 0 optimization was ok
      if (idead.ne.0) idead = k2

      call isgood (i,j,idead)

      end 

      subroutine isgood (i,j,idead)
c-----------------------------------------------------------------------
c isgood - for sucessful gridded minimization sort and index the 
c          assemblage. for bad minimizations flag the grid/assemblage
c          pointer. in either case increment count stats.
c          if idead > 0, should be either k2 or k2-1 which get mapped to
c          iap(idead) = k3 or k3-1
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,idead
c-----------------------------------------------------------------------
      if (idead.eq.0) then 
c                                 all systems go
         rcount(4) = rcount(4) + 1
c                                 at this point the compositions of
c                                 the np solutions are in cp3, ctot3, x3 indexed
c                                 by np and the original indices of the 
c                                 ncpd compounds are -kkp(np+1..ntot), and
c                                 the molar amounts of the phases are in amt.
         call sorter (igrd(i,j),i,j)

      else

         rcount(5) = rcount(5) + 1

         igrd(i,j) = idead
         iap(idead) = idead + k3-k2

      end if  

      end 

      subroutine stblk1 (i,j,loopx,loopy,idead)
c-----------------------------------------------------------------------
c stblk1:

c 1) generates 2d bulk composition at (i,j)
c 2) if closed space (lopt(1)), tests that the composition is bounded 
c 3) checks for degeneracy and negative compositions.

c for out-of-bound and negative compositions the node is assigned 
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, loopx, loopy, idead

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c------------------------------------------------------------------------
      idead = 0
c                                 compute bulk coordinate at node i,j
      cx(1) = (i-1)/dfloat(loopx-1)
      cx(2) = (j-1)/dfloat(loopy-1)

      if (lopt(1)) then 
c                                 do a quick check for closed composition
c                                 if bad, assign bad result and return
         if (cx(1)+cx(2).gt.r1) then
            igrd(i,j) = k2
            idead = 2
            return
         end if

      end if 
c                                 compute the real compositions
      call setblk
c                                 checks for degeneracy and out-of-bounds
c                                 compositions; this check obviates the 
c                                 need for fancy loop indexing.
      call chkblk (idead)

      if (idead.ne.0) igrd(i,j) = k2

      end

      subroutine lpopt1 (idead,statik)
c-----------------------------------------------------------------------
c lpopt1 - does optimization for george's liquidus search without saving
c or finalizing the results. results must be finalized by calling rebulk
c and sorter, the call to rebulk needs to know is static or not.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, idead, inc, lphct, jter, lpprob

      double precision ax(k5),x(k1),oldt,oldp,gtot,
     *                 tol,oldx,clamda(k1+k5)

      logical quit, statik

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer jphct,istart
      common/ cst111 /jphct,istart

      double precision g
      common/ cst2 /g(k1)

      integer hcp,idv
      common/ cst52 /hcp,idv(k7)

      integer tphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision wmach
      common/ cstmch /wmach(10)

      double precision bl,bu
      common/ cstbup /bl(k1+k5),bu(k1+k5)

      save ax, x, clamda
c-----------------------------------------------------------------------
c                                 check for positive bulk, this
c                                 is for icont = 2 with closed
c                                 compositions space (lopt(1) = T), george's
c                                 indexing should eliminate the 
c                                 possibility.
      idead = 0

      statik = .true.

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
      do k = 1, jpoint
         g2(k) = c(k)
      end do 
c                                 load the bulk into the constraint array
      bl(jphct+1:jphct+icp) = b(1:icp)
      bu(jphct+1:jphct+icp) = b(1:icp)

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

      else
c                                 save lphct to recover static solution if
c                                 no refinement 
         lphct = jphct 
c                                 find discretization points for refinement
         call yclos1 (x,clamda,jphct,quit)
c                                 returns quit if nothing to refine
         if (quit) then 
c                                 final processing, .true. indicates static
            statik = .true.
c                                 for recovery of a previous optimization 
c                                 result by savlst, call savpa, currently
c                                 this is only done by liqdus
            call savpa (statik)

         else
c                                 initialize refinement point pointers
            hkp(1:ipoint) = 0 
c                                 reoptimize with refinement
            call reopt (idead,gtot)

            if (idead.eq.0) then

               statik = .false.

            else if (idead.eq.-1) then
c                                 hail mary
               jphct = lphct
               idead = 0

               call yclos0 (x,is,jphct) 

            end if 

         end if 

      end if

      t = oldt
      p = oldp
      xco2 = oldx

c                                 if idead = 0 optimization was ok
      if (idead.eq.0) then 

         rcount(4) = rcount(4) + 1

      else

         rcount(5) = rcount(5) + 1

      end if

      end 

      subroutine amihot (i,j,jhot,jinc)
c----------------------------------------------------------------------
c check if cell with lower left corner at node ij is homogeneous
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jinc, i, j, jhot
c----------------------------------------------------------------------
      jhot = 1

      if (iap(igrd(i,j)).eq.iap(igrd(i,j+jinc)).and.
     *    iap(igrd(i,j)).eq.iap(igrd(i+jinc,j+jinc)).and.
     *    iap(igrd(i,j)).eq.iap(igrd(i+jinc,j))) jhot = 0

c     if (lopt(18).and.iap(igrd(i,j)).eq.k3) jhot = 1
     
      end  

      subroutine filler (i,j,kinc)
c--------------------------------------------------------------------
c fill nodes between identical edges or diagonals, this could
c be modified to a less biased (symmetrical) strategy as in 
c routine aminot. 
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer kinc, i, j, ii, jj, k
c----------------------------------------------------------------------
      if (kinc.eq.1) return 

      if (igrd(i,j).eq.igrd(i+kinc,j+kinc)) then
c                                right diagonal
c        if (igrd(i,j).ne.k2.or.(.not.lopt(18))) then 
         do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
            if (igrd(i+k,j+k).eq.0) igrd(i+k,j+k) = igrd(i,j)
         end do
c        end if  

      else if (igrd(i+kinc,j).eq.igrd(i,j+kinc)) then
c                                left diagonal
c        if (igrd(i+kinc,j).ne.k2.or.(.not.lopt(18))) then 
         do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
            if (igrd(i+k,j+kinc-k).eq.0) 
     *          igrd(i+k,j+kinc-k) = igrd(i,j+kinc)
         end do 
c        end if 

      end if
c                                bottom and top edge
      do jj = j, j+kinc, kinc
c        if (lopt(18).and.igrd(i,jj).eq.k2) then 
c           cycle 
         if (igrd(i,jj).eq.igrd(i+kinc,jj)) then 
            do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
               if (igrd(i+k,jj).eq.0) igrd(i+k,jj) = igrd(i,jj)
            end do 
         end if
      end do 
c                                left and right edges
      do ii = i, i+kinc, kinc
         if (igrd(ii,j).eq.igrd(ii,j+kinc)) then 
            do k = 1, kinc-1
c                                the conditional prevents 
c                                overwriting of identities if 
c                                compression is off and searching
c                                for true boundaries.
               if (igrd(ii,j+k).eq.0)  igrd(ii,j+k) = igrd(ii,j)
            end do 
         end if
      end do 

      end

      subroutine liqdus
c--------------------------------------------------------------------
c liqdus does liquidus finding on a 1- or 2-d compositional grid.
c at each compositional step across the grid, it does a hunt for the 
c temperature at which the last solid disappears.

c the phase assemblage at node(i,j) of the grid is identified by the 
c pointer to a phase assemblage.

c George Helffrich, 5/23.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical init

      integer kinc, kinc2, kinc21, icent, jcent, ie, je, ii, jj,
     *        iic, iil, jjc, jjl, klow,
     *        jtic, ktic, icell, ihot, jhot, hhot, khot,
     *        i, j, k, h, hh, kk, ll, idead,
     *        nblen

      integer iind(4), jind(4), iiind(4,2), jjind(4,2),
     *        icind(4), jcind(4), ieind(5), jeind(5),
     *        jinc(l8), lhot(4),
     *        hotij(l7,2), kotij(l7,2)

      external nblen

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      save init
      data init/.true./

      save iind, jind, iiind, jjind, icind, jcind
      data iind, jind   /0,0,1,1, 0,1,1,0/
      data iiind, jjind /0,0,1,1,1,1,2,2, 1,1,2,0,0,2,1,1/
      data icind, jcind /0,1,2,1, 1,2,1,0/
      data ieind, jeind /0,0,2,2,0, 0,2,2,0,0/
c-----------------------------------------------------------------------
      if (init) call initlq
c                               initialize assemblage counter
      iasct = 0 
      ibulk = 0 
c                               load arrays for lp solution
      call initlp 
c                               jlow is the number of nodes
c                               at the lowest level, the number of
c                               nodes is
c                               first level:
      loopy = (jlow-1) * 2**(jlev-1) + 1 

      loopx = (loopx-1) * 2**(jlev-1) + 1 

      if (loopx.gt.l7) then
         call warn (92,v(iv1),loopx,'x_node')
         klow = (l7 - 1)/2**(jlev-1)
         loopx = klow * 2**(jlev-1) + 1 
      end if  

      if (loopy.gt.l7) then
         call warn (92,v(iv1),loopy,'y_node')
         klow = (l7 - 1)/2**(jlev-1)
         loopy = klow * 2**(jlev-1) + 1 
      end if
c                               initialize igrd (this is critical 
c                               for auto_refine).
      do j = 1, loopy
         do i = 1, loopx
            igrd(i,j) = 0
         end do 
      end do 
c                               could check here if loopx*loopy, the
c                               theoretical max number of assemblages
c                               is > k2, but in practice the number of
c                               assemblages << k2, so only test when 
c                               actually set.
      do j = 1, k2
         iap(j) = 0 
      end do 
c                               increments at each level
      do j = 1, jlev
         jinc(j) = 2**(jlev-j)
      end do 

      kinc = jinc(1)
      jinc1 = kinc

      call setvar

      ktic = 0

      if (init) then
         write (*,1050) 'Beginning',
     *                whatlq(1:nblen(whatlq)),vname(iv1)(1:1),
     *                nopt(2),unitlq(1:nblen(unitlq))
      else
         write (*,1050) 'Continuing',
     *                whatlq(1:nblen(whatlq)),vname(iv1)(1:1),
     *                nopt(2),unitlq(1:nblen(unitlq))
      end if
c                              now traverse compositional grid:
c                              lower triangle; upper is symmetric across diag.
      do i = 1, loopx, kinc

         do j = 1, loopy, kinc
c                                 set bulk, check limits and degeneracy
            call stblk1 (i,j,loopx,loopy,idead) 

            if (idead.ne.0) cycle
c                                 look for the liquidus:
            call fndliq (i,j,ktic,idead)

         end do

      end do
c                                 get hot points
      ihot = 0 
      kinc2 = kinc/2
      kinc21 = kinc2 + 1

      do i = 1, loopx - kinc, kinc
         do j = 1, loopy - kinc, kinc
            if (igrd(i,j+kinc).le.0 .or. igrd(i+kinc,j).le.0) then
               jhot = 1
            else
               call amihot (i,j,jhot,kinc)
            end if
            if (jhot.ne.0) then 
               ihot = min(l7,ihot + 1)
               hotij(ihot,1) = i
               hotij(ihot,2) = j 
c                                 cell is heterogeneous
c                                 fill in homogeneous diagonals
c                                 and edges
               if (iopt(18).ne.0.and.kinc.gt.1) call filler (i,j,kinc)
            else 
c                                 cell is homogeneous
c                                 fill in entire cell
               if (kinc.gt.1) call aminot (i,j,kinc,kinc2,kinc21)
            end if 
         end do 
      end do

      if (ihot.eq.0) goto 10

      ktic = (loopx/kinc+1)*(loopy/kinc+1)
c                              now refine on all higher levels:
      do k = 2, jlev
c                              set new hot cell counter
         khot = 0 
         jtic = 0 
c                              now working on new level
         kinc = jinc(k)
         kinc2 = kinc/2
c
         write (*,1065) ihot,k
c                              flush stdout for paralyzer
         flush (6)
c                              compute assemblages at refinement
c                              points
         do h = 1, ihot
c                              cell corner
            iic = hotij(h,1) 
            jjc = hotij(h,2)
c                              first compute the central node of
c                              the hot cell
            icent = iic + kinc
            jcent = jjc + kinc
c                              forget cells already on grid diagonal
            if (jjc.ge.loopy-(iic-1)-kinc) cycle

            if (igrd(icent,jcent).eq.0) then

               call stblk1 (icent,jcent,loopx,loopy,idead) 

               if (idead.eq.0) call fndliq (icent,jcent,
     *                                      jtic,idead)

            end if 
c                              now determine which of the diagonals
c                              has a change
            hhot = 0 

            do hh = 1, 4

               i = iic + iind(hh)*2*kinc
               j = jjc + jind(hh)*2*kinc
               lhot(hh) = 0

               if (iap(igrd(i,j)).ne.iap(igrd(icent,jcent))) then
c                              cell is hot
                  khot = khot + 1
                  hhot = hhot + 1
                  kotij(khot,1) = iic + iind(hh)*kinc
                  kotij(khot,2) = jjc + jind(hh)*kinc
                  lhot(hh) = 1
c                              compute assemblages at new nodes
                  do kk = 1, 2

                     ii = iic + iiind(hh,kk)*kinc
                     jj = jjc + jjind(hh,kk)*kinc

                     if (igrd(ii,jj).eq.0) then

                        call stblk1 (ii,jj,loopx,loopy,idead) 

                        if (idead.eq.0) call fndliq (ii,jj,jtic,idead)

                     end if 
                  end do 
               end if 
            end do 
c                              if less than 3 hot sub-cells check
c                              edges
            if (hhot.lt.4.and.hhot.gt.1) then 

               do hh = 1, 4
c                              index the edge node
                  ii = iic + icind(hh)*kinc
                  jj = jjc + jcind(hh)*kinc

                  if (igrd(ii,jj).ne.0) then 
c                              could have a second hot cell, check
c                              both corners
                     icell = hh

                     do kk = 1, 2
                        ie = iic + ieind(hh+kk-1)*kinc
                        je = jjc + jeind(hh+kk-1)*kinc
                        icell = icell + kk - 1
                        if (icell.gt.4) icell = 1

                        if (iap(igrd(ii,jj)).ne.iap(igrd(ie,je)).and.
     *                     lhot(icell).eq.0) then 
c                               new cell
                           khot = min(l7,khot + 1)
                           hhot = hhot + 1
                           lhot(icell) = 1
c                                cell index is 
                           ii = iic + iind(icell)*kinc
                           jj = jjc + jind(icell)*kinc
                           kotij(khot,1) = ii
                           kotij(khot,2) = jj
c                                compute assemblage at cell nodes
                           do ll = 1, 4

                              iil = ii + iind(ll)*kinc
                              jjl = jj + jind(ll)*kinc

                              if (igrd(iil,jjl).eq.0) then

                                 call stblk1 (iil,jjl,loopx,loopy,idead)

                                 if (idead.eq.0) call fndliq (iil,jjl,
     *                                   jtic,idead)

                              end if 
                           end do  
                        end if
                     end do 
                  end if 
               end do 
            end if      

            do hh = 1, 4

               i = iic + iind(hh)*kinc
               j = jjc + jind(hh)*kinc
               if (i.lt.loopx.and.j.lt.loopy) then 
                  if (lhot(hh).eq.0) then 
c                                fill cold cells
                     call aminot1 (icent,jcent,i,j,kinc)
                  else 
c                                fill hot cells
                     if (iopt(18).ne.0) call filler (i,j,kinc)
                  end if 
               end if 
            end do 
         
         end do

         write (*,1070) k,jtic
         ktic = ktic + jtic

         write (*,1080) ktic,(loopx/kinc)*(loopy/kinc+1)/2*18

         if (khot.eq.0.or.k.eq.jlev) exit 
c                             now switch new and old hot list
         ihot = khot
         do i = 1, khot
            hotij(i,1) = kotij(i,1)
            hotij(i,2) = kotij(i,2)
         end do 

      end do 

      write (*,1060) rcount(5),loopx*(loopy+1)/2
      write (*,1080) rcount(4),loopx*(loopx+1)/2*18

10    if (outprt) call outgrd (loopx,loopy,1,n4,0)

      init = .false.

1040  format (2(i4,1x),a,a)
1050  format (/,3(a,1x),'refinement ',
     *        'to +/-',f6.2,1x,a,1x,'tolerance.',/)
1060  format (/,i6,' grid cells of ',i6,
     *             ' failed liquidus/solidus search.',/)
1065  format (/,i6,' grid cells to be refined at grid level ',i1)
1070  format (/,7x,'refinement at level ',i1,' involved ',i6,
     *             ' minimizations')
1080  format (i6,' minimizations required of the ',
     *        'theoretical limit of ',i6)
1090  format (a,7x,'...working (',i6,' minimizations done)',a,$)

      end 

      subroutine fndliq (i,j,ktic,idead)
c--------------------------------------------------------------- 
c fndliq iterates on element (i,j) in the grid to locate the liquidus or
c solidus assemblage to within tolerance tol.
c opts are search options, encoded as two bits:
c    x 0 = finding liquidus
c    x 1 = finding solidus
c    0 x = T search (high -> liquid, low -> solid)
c    1 x = P search (low -> liquid, high -> solid)
c ktic is an iteration counter.
c liq(1:nliq) is a list of liquid phases.
c on return, idead ne 0 if there is a failure to find a liquidus/solidus
c assemblage.
c 
c clsliq returns: type = 0 if no liquid
c                 type = 1 if liquid + solid
c                 type = 2 if liquid only
c
c type of search: iv1 = 1 -> P
c                 iv1 = 2 -> T

c George Helffrich, 5/23

c the trick here is that for liquidus calculations the last subliquidus
c result needs to be saved, whereas for solidus calculations the last
c supersolidus result should be saved.
c--------------------------------------------------------------- 
      implicit none

      include 'perplex_parameters.h'

      logical sol, pl, statik, abort

      integer i, j, k, l, ktic, idead

      double precision tlo,thi

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c--------------------------------------------------------------- 
      sol = mod(opts,2) .eq. 1
      pl = opts/2 .eq. 1

      tlo = vmin(iv1)
      thi = vmax(iv1)
c                              check the assemblage at the minimum
      v(iv1) = vmin(iv1)
c                              update dependent variables, if any
      call incdp0
c                              chkblk returns idead ~0 only if composition
c                              is out of bounds, george's loops should never
c                              generate this case.
      call chkblk (idead)

      if (idead.ne.0) then 

         write (*,*) 'outta bounds?',i,j

         return

      end if

      ktic = ktic + 1
      if (0.eq.mod(ktic,500)) write (*,1090) cr,ktic

      call lpopt1 (idead,statik)

      if (idead .ne. 0) then

         write (*,1020) 'low',vname(iv1),i,j

         call isgood (i,j,k2)

         return

      end if

      call clsliq (l)

      if (iv1.eq.2 .and.
     *    ((sol .and. l.ne.0) .or. (.not.sol .and. l.eq.2))) then
c                              only warn if past exploratory phase
         call liqwrn (i,j,'no solids','lowest')

         call isgood (i,j,k2-1)

         return

      end if

      if (iv1.eq.1 .and. l.eq.0) then
c                              only warn if past exploratory phase
         call liqwrn (i,j,'no liquid','lowest')

         call isgood (i,j,k2-1)

         return

      end if
c                              check the assemblage at the maximum
      v(iv1) = vmax(iv1)
c                              update dependent variables, if any
      call incdp0
c                              do the optimization
      ktic = ktic + 1

      if (0.eq.mod(ktic,500)) write (*,1090) cr,ktic

      call lpopt1 (idead,statik)

      if (idead.ne.0) then

         write (*,1020) 'high',vname(iv1),i,j

         call isgood (i,j,k2)

         return

      end if

      call clsliq (l)

      if (iv1.eq.2 .and. ((.not.sol .and. l.ne.2) .or.
     *                     (sol .and. l.eq.0))) then
c                              only warn if past exploratory phase
         call liqwrn (i,j,'solids','highest')

         call isgood (i,j,k2-1)

         return

      else if (iv1.eq.1 .and. l.ne.0) then
c                              only warn if past exploratory phase
         call liqwrn (i,j,'liquid','highest')

         call isgood (i,j,k2-1)

         return

      end if
c                                 iterate by narrowing interval to 1/2**16
c                                 or to uncertainty < nopt(2)
      do k = 1, 16

        ktic = ktic + 1

        v(iv1) = (tlo+thi)/2
c                                 update dependent variables, if any
        call incdp0
c                                 do the optimization
        call lpopt1 (idead,statik)

        if (idead .ne. 0) exit

        call clsliq (l)

        if (pl) then

           if (l .eq. 0) then
c                                 s: upper bound, save solid
              thi = v(iv1)

           else if (l .eq. 1) then
c                                 s+l: save solid
              if (sol) then
c                                 solidus: lower bound
                 tlo = v(iv1)

              else
c                                 liquidus: upper bound
                 thi = v(iv1)

              end if

           else
c                                 l: lower bound
              tlo = v(iv1)

           end if

        else

           if (l .eq. 2) then
c                                 if all liquid, upper bound
              thi = v(iv1)

           else if (l .ge. 1) then
c                                 s+l: save solid
              if (sol) then
c                                 solidus: upper bound
                 thi = v(iv1)

              else
c                                 liquidus: lower bound
                 tlo = v(iv1)

              end if

           else
c                                 lower bound, keep solid assemblage
              tlo = v(iv1)

           end if

        end if

        if (0.eq.mod(ktic,500)) write (*,1090) cr,ktic

        if (thi - tlo .lt. nopt(2)) exit
c                                 save the last wrong-side result
        if (sol.and.l.eq.0 .or. .not.sol.and.l.ne.2) then

           call savlst (.false.,statik,l)

        end if

      end do
c                                 final processing:
      if (idead.ne.0) then
c                                 fndliq failed
         if (refine) write (*,1020) 'liquidus grid',i,j
c                                 here's an opportunity to set
c                                 a bad value for the temperature.
         idead = k2
      else 

         if (sol.and.l.ne.0 .or. .not.sol.and.l.eq.2) then
c                                 on the solid side of the solidus or
c                                 the liquid side of the liquidus, back
c                                 off to the last L+S result
             call savlst (.true.,statik,l)

         end if
c                                 finalize and save liquidus assemblage
         call rebulk (abort,statik)
c                                 rebulk can set abort, but this would
c                                 be for electrolytic fluids
         if (abort) idead = 99

      end if

      call isgood (i,j,idead)

1020  format (/,'**Unable to define',2(1x,a),' assemblage: ind.',
     *        2(1x,i5))
1090  format (a,7x,'...working (',i6,' minimizations done)')

      end

      subroutine liqwrn (i,j,a,b)
c----------------------------------------------------------------------
c warn of failed initial condition tests in fndliq
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character a*(*),b*(*), assmb*128, text*240

      integer i, j, l, nblen

      external nblen

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
      call smptxt (assmb,l)

      write (text,1010) i, j, a, b, vname(iv1)(1:1), assmb(1:l)

      call deblnk (text)

      write (*,'(/,a)') text(1:nblen(text))
      write (*,1020) vname(iv1)

1010  format ('**warning ver327**',2(1x,i5),' has ',a,
     *          ' at',2(1x,a),': ',a)
1020  format (/,2x,'Possible causes for this problem include:',/,
     *          4x,'1 - an unduly restricted search range for ',a,/,
     *          4x,'2 - stability of melt endmembers not ',
     *             'included in the melt model.',/)

      end 

      subroutine savlst (recov,statik,l)
c----------------------------------------------------------------------
c save (~recov) or recover (recov) the information necessary to 
c reconstruct the previous optimization result during iteration.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical recov, statik

      integer l, i, ids

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      logical xstic
      integer xnpt, xjdv, xl, xlcoor, xlkp
      double precision xamt, twrong, xycoor
      common/ cstlst /xamt(k19), twrong, xycoor(k22), xlkp(k19),
     *                xjdv(k19), xlcoor(k19), xnpt, xstic, xl

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c----------------------------------------------------------------------
      if (.not.recov) then

         xl = l
         twrong = v(iv1)
         xstic = statik
         xnpt = npt

         do i = 1, npt

            xjdv(i) = jdv(i)
            xamt(i) = amt(i)
            xlkp(i) = lkp(i)
c                                 some how compounds can get into
c                                 the dynamic list, i.e., jdv(i) may point to 
c                                 a compound at i > jpoint. this should be prevented
            if (jdv(i).gt.jpoint.and.lkp(i).lt.0) then 
               write (*,*) 'oinkers ',jdv(i), jkp(jdv(i)), lkp(i)
            end if

            if (lkp(i).lt.0) cycle

            ids = lkp(i)
            xlcoor(i) = lcoor(i)
            xycoor(lcoor(i)+1:lcoor(i)+nstot(ids)) = 
     *                       ycoor(lcoor(i)+1:lcoor(i)+nstot(ids))

         end do

      else

         statik = xstic
         npt = xnpt

         do i = 1, npt
c                                 need to reset ctot2 etc?
            jdv(i) = xjdv(i)
            amt(i) = xamt(i)
            ids = xlkp(i)
            lkp(i) = ids

            if (ids.lt.0) cycle

            jkp(jdv(i)) = ids

            lcoor(i) = xlcoor(i)
            ycoor(lcoor(i)+1:lcoor(i)+nstot(ids)) = 
     *                       xycoor(lcoor(i)+1:lcoor(i)+nstot(ids))

         end do

      end if

      end

      subroutine smptxt (string,iend)
c----------------------------------------------------------------------
c subprogram to write a text labels for bulk composition output 
c id identifies the assemblage
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character string*(*), char14*14

      integer i, ist, iend

      integer length,com
      character chars*1, card*400
      common/ cst51 /length,com,chars(400),card

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
c----------------------------------------------------------------------
      iend = 0

      string = ' '

      ist = 1

      do i = 1, lchar
         chars(i) = ' '
      end do

      do i = 1, npt

         call getnam (char14,jkp(jdv(i)))

         ist = iend + 1
         iend = ist + 14
         read (char14,'(400a)') chars(ist:iend)

         call ftext (ist,iend)

      end do 

      write (string,'(400a)') chars(1:iend) 

      length = iend

      end 

      subroutine clsliq (type)
c----------------------------------------------------------------------
c classify grid item one of three ways:
c type = 0 if no liquid
c type = 1 if liquid + solid
c type = 2 if liquid only
c--------------------------------------------------------------- 
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, type

      logical is, liq, sol

      integer ikp
      common/ cst61 /ikp(k1)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c--------------------------------------------------------------- 
      type = 0
      liq = .false.
      sol = .false.

      do i = 1, npt
c                                 the problem with this test is liqlst may
c                                 contain negative index of a compound or the
c                                 positive index of a solution. test for both
c                                 cases.
         if (lkp(i).lt.0) then

            ids = ikp(-lkp(i))
            id = lkp(i)

         else

            ids = lkp(i)
            id = 0

         end if

         do j = 1, nliq

            if (liqlst(j).gt.0) then
c                                 compare solution model indices
               is = ids.eq.liqlst(j)

            else
c                                 compare compound indices
               is = id.eq.liqlst(j)

            end if

            if (is) exit

         end do

         liq = liq .or. is
         sol = sol .or. .not. is

      end do

      if (liq) then

         if (sol) then
            type = 1
         else
            type = 2
         end if

      end if

      end

      subroutine wav2d1
c--------------------------------------------------------------------
c wav2d does constrained minimization on a 2 dimensional multilevel
c grid (ith column of a 2-d grid), lowest resolution is the default, 
c and phase boundaries are located at the highest level

c the phase assemblage at node(i,j) of the grid is identified by the 
c pointer to a phase assemblage.

c the compound grid data igrd is output to the "plot" file.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer kinc, jinc(l8), iind(4), jind(4), iiind(4,2),htic,
     *        jjind(4,2), hotij(l7*l7,2), kotij(l7*l7,2), icind(4),
     *        jcind(4), lhot(4), ieind(5), jeind(5), i, ihot,
     *        iil, jjl, ll, je, ie, icell,
     *        ii,jj,kk,hh,hhot,jcent,icent,jjc,iic,h,jtic,khot,k,
     *        ktic,j,jhot,klow,kinc2,kinc21,idead

      double precision dinc, tot

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1

      save iind, jind, iiind, jjind, icind, jcind

      data iind, jind   /0,0,1,1, 0,1,1,0/
      data iiind, jjind /0,0,1,1,1,1,2,2, 1,1,2,0,0,2,1,1/
      data icind, jcind /0,1,2,1, 1,2,1,0/
      data ieind, jeind /0,0,2,2,0, 0,2,2,0,0/
c-----------------------------------------------------------------------
c                               initialize assemblage counter
      iasct = 0 
      ibulk = 0 
c                               call old routine for 1d grid
      if (loopx.eq.1) then
         call wavgrd
         return
      end if 
c                               jlow is the number of nodes
c                               at the lowest level, the number of
c                               nodes is
      klow = jlow - 1
c                               first level:
      loopy = klow * 2**(jlev-1) + 1 

      loopx = (loopx-1) * 2**(jlev-1) + 1 

      if (loopy.gt.l7) then
         call warn (92,v(iv1),loopy,'y_node')
         klow = (l7 - 1)/2**(jlev-1)
         loopy = klow * 2**(jlev-1) + 1 
      end if

      if (loopx.gt.l7) then
         call warn (92,v(iv1),loopx,'x_node')
         klow = (l7 - 1)/2**(jlev-1)
         loopx = klow * 2**(jlev-1) + 1 
      end if
c                               initialize igrd (this is critical 
c                               for auto_refine).
      igrd = 0
c                               could check here if loopx*loopy, the
c                               theoretical max number of assemblages
c                               is > k2, but in practice the number of
c                               assemblages << k2, so only test when 
c                               actually set.
      iap = 0
c                               increments at each level
      do j = 1, jlev
         jinc(j) = 2**(jlev-j)
      end do 

      kinc = jinc(1)
      jinc1 = kinc

      call setvar
c                               init progress info
      dinc = 1d2/real((loopx-1)/kinc + 1)
      tot = 0d0

c     if (lopt(28)) call begtim (11)
c                               do all points on lowest level
      do i = 1, loopx, kinc
         do j = 1, loopy, kinc
            call setvr0 (i,j)
            call lpopt (i,j,idead)
         end do
c                               progress info
         tot = tot + dinc
         write (*,1030) tot,char(13)
c                               flush stdout for paralyzer
         flush (6)

      end do

      write (*,1030)

c     if (lopt(28)) call endtim (11,.true.,'low level grid')
c                               output interim plt file
      if (iopt(34).ne.0) call outgrd (loopx,loopy,kinc,1000,1)
c                               get hot points
      ihot = 0 
      kinc2 = kinc/2
      kinc21 = kinc2 + 1

      do i = 1, loopx-1, kinc
         do j = 1, loopy-1, kinc
            if (igrd(i,j+kinc).le.0 .or. igrd(i+kinc,j).le.0) then
               jhot = 1
            else
               call amihot (i,j,jhot,kinc)
            end if
            if (jhot.ne.0) then 
               ihot = ihot + 1
               hotij(ihot,1) = i
               hotij(ihot,2) = j 
c                               cell is heterogeneous
c                               fill in homogeneous diagonals
c                               and edges
               if (iopt(18).ne.0.and.kinc.gt.1) call filler (i,j,kinc)
            else 
c                               cell is homogeneous
c                               fill in entire cell
               if (kinc.gt.1) call aminot (i,j,kinc,kinc2,kinc21)
            end if 
         end do 
      end do

      if (ihot.eq.0) goto 10 

      ktic = (loopx/kinc+1)*(loopy/kinc+1)

      if (jlev.gt.1) write (*,1050) 

      htic = 0 
c                              now refine on all higher levels:
      do k = 2, jlev
c                              set new hot cell counter
         khot = 0 
         jtic = 0 
c                              now working on new level
         kinc = jinc(k)
         kinc2 = kinc/2
c
         write (*,1060) ihot,k
c                               flush stdout for paralyzer
         flush (6)

c        if (lopt(28)) call begtim (12)
c                              compute assemblages at refinement
c                              points
         do h = 1, ihot
c                              cell corner
            iic = hotij(h,1) 
            jjc = hotij(h,2)
c                              first compute the central node of
c                              the hot cell
            icent = iic + kinc
            jcent = jjc + kinc
            if (igrd(icent,jcent).eq.0) then 
               call setvr0 (icent,jcent)
               call lpopt (icent,jcent,idead)
               jtic = jtic + 1
               htic = htic + 1
            end if 
c                              now determine which of the diagonals
c                              has a change
            hhot = 0 

            do hh = 1, 4
            
               i = iic + iind(hh)*2*kinc
               j = jjc + jind(hh)*2*kinc
               lhot(hh) = 0 

               if (iap(igrd(i,j)).ne.iap(igrd(icent,jcent))) then 
c                              cell is hot
                  khot = khot + 1
                  hhot = hhot + 1
                  kotij(khot,1) = iic + iind(hh)*kinc
                  kotij(khot,2) = jjc + jind(hh)*kinc
                  lhot(hh) = 1
c                              compute assemblages at new nodes
                  do kk = 1, 2
                     ii = iic + iiind(hh,kk)*kinc
                     jj = jjc + jjind(hh,kk)*kinc
                     if (igrd(ii,jj).eq.0) then 
                        call setvr0 (ii,jj)
                        call lpopt (ii,jj,idead)
                        jtic = jtic + 1
                        htic = htic + 1
                     end if 
                  end do 
               end if 
            end do 
c                              if less than 3 hot sub-cells check
c                              edges
            if (hhot.lt.4.and.hhot.gt.1) then 

               do hh = 1, 4
c                              index the edge node
                  ii = iic + icind(hh)*kinc
                  jj = jjc + jcind(hh)*kinc
                  if (igrd(ii,jj).ne.0) then 
c                              could have a second hot cell, check
c                              both corners
                     icell = hh

                     do kk = 1, 2
                        ie = iic + ieind(hh+kk-1)*kinc
                        je = jjc + jeind(hh+kk-1)*kinc
                        icell = icell + kk - 1
                        if (icell.gt.4) icell = 1
                        
                        if (iap(igrd(ii,jj)).ne.iap(igrd(ie,je)).and.
     *                     lhot(icell).eq.0) then 
c                               new cell
                           khot = khot + 1
                           hhot = hhot + 1
                           lhot(icell) = 1
c                                cell index is 
                           ii = iic + iind(icell)*kinc
                           jj = jjc + jind(icell)*kinc
                           kotij(khot,1) = ii
                           kotij(khot,2) = jj
c                                compute assemblage at cell nodes
                           do ll = 1, 4
                              iil = ii + iind(ll)*kinc
                              jjl = jj + jind(ll)*kinc
                              if (igrd(iil,jjl).eq.0) then
                                 call setvr0 (iil,jjl)
                                 call lpopt (iil,jjl,idead)
                                 jtic = jtic + 1
                                 htic = htic + 1
                              end if 
                           end do  
                        end if
                     end do 
                  end if 
               end do 
            end if      

            do hh = 1, 4

               i = iic + iind(hh)*kinc
               j = jjc + jind(hh)*kinc
               if (i.lt.loopx.and.j.lt.loopy) then 
                  if (lhot(hh).eq.0) then 
c                                fill cold cells
                     call aminot1 (icent,jcent,i,j,kinc)
                  else 
c                                fill hot cells
                     if (iopt(18).ne.0) call filler (i,j,kinc)
                  end if 
               end if 
            end do 
         
            if (htic.gt.500) then 
               write (*,1090) jtic, char(13)
               htic = 0 
            end if
 
         end do

         write (*,1070) k,jtic
         ktic = ktic + jtic

         write (*,1080) ktic,(loopx/kinc+1)*(loopy/kinc+1)

c        if (lopt(28)) call endtim (12,.true.,'nth level grid')

         if (khot.eq.0.or.k.eq.jlev) exit 
c                             now switch new and old hot list
         ihot = khot
         do i = 1, khot
            hotij(i,1) = kotij(i,1)
            hotij(i,2) = kotij(i,2)
         end do
c                               output interim plt file
         if (iopt(34).ne.0) call outgrd (loopx,loopy,kinc,1000,k)

      end do
c                                 ouput grid data
10    if (outprt) call outgrd (loopx,loopy,1,n4,0)

1030  format (f5.1,'% done with low level grid.',a,$)
1050  format (//,'Beginning grid refinement stage.',/)
1060  format (i6,' grid cells to be refined at grid level ',i1)
1070  format (7x,'refinement at level ',i1,' involved ',i6,
     *        ' minimizations')
1080  format (i6,' minimizations required of the ',
     *        'theoretical limit of ',i7)
1090  format (7x,'...working (',i6,' minimizations done)',a,$)

      end 

      subroutine wavgrd
c----------------------------------------------------------------------
c wavgrd does constrained minimization on a 1 dimensional multilevel
c grid (ith column of a 2-d grid), lowest resolution is the default, 
c and phase boundaries are located at the highest level

c the phase assemblage at node(i,j) of the grid is identified by the 
c pointer igrd(i,j), igrd is a pointer to a phase assemblage.
c the compound grid data igrd is output to the "plot" file.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer kinc,jinc(l8),i,j,k,klast,kmax,kd,idead,jtic,jmin,
     *        jmax,klev,klow,ltic 

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1
c-----------------------------------------------------------------------
c                               jlow is the number of nodes
c                               at the lowest level, the number of
c                               nodes is
      klow = jlow - 1
c                               first level:
      loopy = klow * 2**(jlev-1) + 1 

      if (loopy.gt.l7) then
         call warn (92,v(iv1),loopy,'WAVGRD')
         klow = (l7 - 1)/2**(jlev-1)
         loopy = klow * 2**(jlev-1) + 1 
      end if        
c                               initialize igrd/iap (this is critical 
c                               for auto_refine).
      do j = 1, loopy
         igrd(1,j) = 0
         iap(j) = 0 
      end do        
c                                increments
      do j = 1, jlev
         jinc(j) = 2**(jlev-j)
      end do 

      jinc1 = jinc(1)
c                                initialize variables
      call setvar 
c                                initialize grid vars
      kinc = jinc(1)
      klev = 1

      i = 1

      j = 1
      jmax = 1 
      jmin = 0
      jtic = 0
      ltic = 0 

      do while (j.le.loopy) 

         if (igrd(i,j).eq.0) then
c                                call to incvr2 moved here (formerly
c                                at end of loop), april 5, 2006. JADC
            call setvr0 (j,1)
            call lpopt (i,j,idead)
            jtic = jtic + 1
c                                progress information
            ltic = ltic + 1
            if (ltic.eq.20) then 
               write (*,'(i5,a)') jtic,' optimizations done...'
               ltic = 0
            end if

         end if

         kd = iap(igrd(i,j))
c                                jmax should always be the last known node
c                                at lowest level of resolution
         if (j.gt.jmax) then 
            jmax = j
            jmin = j - jinc(1) 
            kmax = kd
         end if 

         if (j.eq.1) then
 
            klast = kd
               
         else if (kd.ne.klast.and.klev.ne.jlev) then
c                                crossed a boundary and not
c                                at max resolution, increment resolution
            klev = klev + 1
            kinc = -jinc(klev)

         else if (kd.eq.klast.and.klev.ne.jlev) then 
c                                backed into field?
            if (kd.eq.kmax) then
c                                go to next low res node
               kinc = jinc(1)
               j = jmax
               klev = 1
               goto 10
            end if 
            
            klev = klev + 1
            kinc = jinc(klev)

         else if (klev.eq.jlev) then
c                                boundary has been located at max
c                                resolution
              
c                                jmin is the minimum value at which 
c                                subsequent searchs can go
            do k = j, jmax - 1
               if (igrd(i,k).eq.0) then 
                  jmin = k - 1
                  goto 20 
               end if 
            end do 
          
20          klast = iap(igrd(i,jmin))
c                                check if there is a grid point with
c                                this assemblage at a higher level
            do k = jmin, jmax - 1

               if (igrd(i,k).eq.0) cycle
               if (iap(igrd(i,k)).eq.klast) jmin = k

            end do 

            if (klast.eq.kmax.or.j.eq.jmax-1) then
c                                the grid can be filled in between j and jmax
c                                go on to the next low level grid point:
               klast = iap(igrd(i,jmax))
               kinc = jinc(1)
               klev = 1
               j = jmax
                  
            else 
c                                find the level to search for next boundary
c                                must be > level 1
               do k = 2, jlev
                  if (jmax-jmin.gt.jinc(k)) then
                     j = jmax
                     klev = k
                     kinc = -jinc(k)
                     goto 10 
                  end if
               end do 
            end if 
         end if 

10       j = j + kinc

         if (j.lt.jmin+1) then 
            j = j - kinc
            klev = klev + 1
            kinc = -jinc(klev)
            goto 10 
         end if

      end do
c                                 output graphics data
      if (outprt) call outgrd (i,loopy,jinc(jlev),n4,0)

      end 

      subroutine aminot (i,j,kinc,kinc2,kinc21)
c--------------------------------------------------------------------
c fill in nodes of a homogeneous cell, called only for lowest level
c grid, aminot1 is called for higher level cells. kinc is always
c even and > 1. modfied to fill from each corner if kinc > 2 in case
c compression is off. JADC sept '04. 

c kinc   - the current grid increment
c kinc2  - kinc/2 (integer division)
c kinc21 - kinc2 + 1
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, kinc, ii, jj, kinc2, kinc21
c                                 the commented version fills
c                                 with just assemblage i,j this
c                                 would be fine if compression is
c                                 on, but is biased otherwise
c      do ii = i, i + kinc - 1
c         do jj = j, j + kinc - 1
c                                 the conditional prevents 
c                                 overwriting of identities if 
c                                 compression is off and searching
c                                 for true boundaries.
c            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i,j)
c         end do
c      end do 
c                                 first fill from lower left (i,j)
      do ii = i, i + kinc2 
         do jj = j, j + kinc2
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i,j)
         end do
      end do
c                                 then from lower right (i+kinc,j)
      do ii = i + kinc21, i + kinc
         do jj = j, j + kinc2
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i+kinc,j)
         end do
      end do
c                                 then from upper left (i,j+kinc)
      do ii = i, i + kinc2
         do jj = j + kinc21, j + kinc
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i,j+kinc)
         end do
      end do
c                                 then from upper right (i+kinc,j+kinc)
      do ii = i + kinc21, i + kinc
         do jj = j + kinc2 + 1, j + kinc
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(i+kinc,j+kinc)
         end do
      end do
     
      end  

      subroutine aminot1 (icent,jcent,i,j,kinc)
c--------------------------------------------------------------------
c fill in nodes of a homogeneous cell, called only for high level
c grids, aminot is called for lowest level cells. kinc is always
c even and > 1.  JADC sept '04. 

c kinc   - the current grid increment
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, ii, jcent, icent, jj, kinc

      do ii = i, i + kinc
         do jj = j, j + kinc
c                                 the conditional prevents 
c                                 overwriting of identities if 
c                                 compression is off and searching
c                                 for true boundaries.
            if (igrd(ii,jj).eq.0) igrd(ii,jj) = igrd(icent,jcent)
         end do
      end do 
     
      end 

      subroutine frname 
c----------------------------------------------------------------------
c frname - prompt for names of phases to be fractionated

c          ifrct - number of phases to be fractionated
c          ifr(ifract) - 0 if cpd, else ikp (solution pointer)
c          jfr(ifract) - cpd pointer
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical first 

      integer jam, jfrct, i

      double precision numb

      character phase(k23)*10, name*100

      integer fmode,ifrct,ifr
      logical gone
      common/ frct1 /fmode,ifrct,ifr(k23),gone(k5)
 
      save first, phase

      data first/.true./
c----------------------------------------------------------------------
      if (first) then 

         first = .false.
c                                 choose mode: 0 - don't fractionate
c                                              1 - specified phases
c                                              2 - all solids
         write (*,1030)
         call rdnumb (numb,numb,fmode,0,.false.)

         if (fmode.eq.1) then 
c                                 get phases to be fractionated 
            ifrct = 1

            do 

               write (*,1040)
               read (*,'(a)') phase(ifrct)

               if (phase(ifrct).eq.' ') exit

               call matchj (phase(ifrct),ifr(ifrct))

               if (ifr(ifrct).eq.0) then

                  write (*,1100) phase(ifrct)
                  cycle

               else if (ifr(ifrct).gt.0) then

                  if (ksmod(ifr(ifrct)).eq.39.and.lopt(32).and.
     *                iopt(22).eq.0) then
c                                 fractionating an electrolytic fluid:
c                                 override solid component depletion
c                                 error trap in yclos2 to allow output,
c                                 no examples where this does anything.
                     lopt(71) = .false.

                     call warn (62,numb,ifrct,phase(ifrct))

                  end if

               end if

               ifrct = ifrct + 1

               if (ifrct.gt.k23) call error (1,0d0,ifrct,'k23')

            end do 

            ifrct = ifrct - 1

         else

            ifrct = 0 

         end if 

      else if (fmode.eq.1) then 
c                                 must be in autorefine cycle, make
c                                 new phase list from old list:
         jfrct = ifrct
         ifrct = 0  
         
         do i = 1, jfrct 

            call matchj (phase(i),jam)

            if (jam.eq.0) cycle 
       
            ifrct = ifrct + 1
            ifr(ifrct) = jam 

         end do

      else

         ifrct = 0 

      end if

      if (fmode.ne.0) then 
c                                 initialize "gone" flags
         do i = 1, k5
            gone(i) = .false.
         end do 

c                                 open fractionation files
         call mertxt (name,prject,'_fractionated_bulk.dat',0)
         open (n0,file=name,status='unknown')
         write (*,1050)

         do i = 1, ifrct 

            call fropen (i,phase(i))

         end do

      end if 

1030  format (/,'Choose computational mode:',/,
     *       5x,'0 - no fractionation [default]',/,
     *       5x,'1 - fractionate specified phases',/,
     *       5x,'2 - fractionate all phases other than liquid'/)
1040  format (/,'Enter the name of a phase to be fractionated',
     *        /,'(left justified, <cr> to finish): ')
1050  format (/,'The fractionated bulk composition will be ',
     *          'written to file: fractionated_bulk.dat',/)
1100  format (/,'No such entity as ',a,', try again: ')

99    end 

      subroutine fropen (i,phase)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nblen

      character phase*10

      external nblen
c-----------------------------------------------------------------------
 
      tfname = '_'//phase//'.dat'

      call unblnk (tfname)

      call mertxt (tfname,prject,tfname,0)

      write (*,1010) phase(1:nblen(phase)), tfname(1:nblen(tfname))

      open (n0+i,file=tfname,status='unknown')

1010  format (/,'The fractionated amount and composition of ',a,/,
     *          'will be written to file: ',a,/)

      end

      subroutine fractr (idead,output,verbos,errr)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character gname*10, lgname*22, phase*14

      external gname, lgname

      integer i,j,k,idead,ier

      logical there(k23), warn, output, quit, liquid, verbos

      double precision mass(k23), tmass, x, errr(k5)
 
      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer fmode,ifrct,ifr
      logical gone
      common/ frct1 /fmode,ifrct,ifr(k23),gone(k5)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c----------------------------------------------------------------------- 
c                                 fractionation effects:
      do i = 1, jbulk
         dcomp(i) = 0d0
         errr(i) = 0d0
      end do 

      do i = 1, ifrct
         there(i) = .false.
      end do 

      if (idead.eq.0) then
c                                 error can only be non-zero for lagged speciation
         do i = 1, jbulk
            errr(i) = -cblk(i)
         end do 
c                                 optimization suceeded
c                                 compute mass fractions in case
c                                 threshold based fractionation:
         tmass = 0d0 
c                                 get mass fractions:
         do j = 1, ntot

            mass(j) = 0d0

            do k = 1, jbulk
c                                 error is the back-calculated bulk - input bulk
               errr(k) = errr(k) + amt(j)*cp3(k,j)
               mass(j) = mass(j) + amt(j)*cp3(k,j)*atwt(k)
            end do

            tmass = tmass + mass(j)

         end do

         do j = 1, ntot
            mass(j) = mass(j)/tmass
         end do

         if (fmode.eq.2) then 
c                               check if liquid is stable
            liquid = .false.

            do j = 1, ntot

               if (lgname(phase,kkp(j)).eq.'liquid') then 
                  liquid = .true. 
                  exit
               end if

            end do

            if (liquid) then 
c                               check if the solids are already
c                               in the fractionation list
               do j = 1, ntot

                  if (lgname(phase,kkp(j)).eq.'liquid') cycle

                  quit = .false.

                  do i = 1, ifrct

                     if (ifr(i).eq.kkp(j)) then
                        quit = .true.
                        exit 
                     end if

                  end do
c                                phase j is in the list so its 
c                                file is alread open, cycle 
                  if (quit) cycle
c                                else open a new file
                  ifrct = ifrct + 1

                  ifr(ifrct) = kkp(j)

                  call fropen (ifrct,phase)
c                                paste in previous coordinates
                  rewind (n0-1)

                  do 
c                                the danger in doing this is that 
c                                v may not be written to scratch with 
c                                full precision
                     read (n0-1,*,iostat=ier) v

                     if (ier.ne.0) then 
c                               presumably eof, reset the file pointers
                        backspace (n0+ifrct)

                        backspace (n0-1)

                        exit 

                     else 

                        write (n0+ifrct,1200) 
     *                                    (v(jv(i)),i=1,ipot),nopt(7),
     *                                    (nopt(7),k=1,jbulk)

                     end if

                  end do

               end do

            end if

         end if

         if (fmode.eq.1.or.fmode.eq.2.and.liquid) then 
c                                 do fractionation
            do i = 1, ifrct

               do j = 1, ntot

                  if (kkp(j).eq.ifr(i)) then 
c                                 the phase to be fractionated
c                                 is present, remove from bulk
                     there(i) = .true.
c                                 simple back calculated speciation
c                                 is being used for fractionation:
                     if (lopt(67)) call aqrxdo (j,-1,.false.)

                     if (amt(j).lt.0d0) amt(j) = 0d0

                     if (mass(j).gt.nopt(33)) then 
c                                 mass fraction exceeds upper threshold
                        x = (mass(j) - nopt(32))/mass(j)

                        if (verbos) then 
c                                 write to console
                           write (*,1185) (vname(iv(k)),v(iv(k)),
     *                                                       k=1,ipot)
                           write (*,1175) (cblk(k),k=1,jbulk)
                           write (*,1190) amt(j),gname(ifr(i)),
     *                                   (amt(j)*cp3(k,j),k=1,jbulk)
                        end if

                        do k = 1, jbulk 
                           dcomp(k) = dcomp(k) + x*amt(j)*cp3(k,j)
                        end do

c                                 write to file
                        if (output) write (n0+i,1200) 
     *                                       (v(jv(k)),k=1,ipot),amt(j),
     *                                       (amt(j)*cp3(k,j),k=1,jbulk)
                     else 
c                                 write to console
                        write (*,1185) (vname(iv(k)),v(iv(k)),k=1,ipot)

                        write (*,1180) gname(ifr(i)),mass(j),nopt(33)

                     end if 

                  end if

               end do

            end do 

         end if 
c                                 write output for fractionated phases 
c                                 that are NOT stable
         do i = 1, ifrct

            if (.not.there(i)) then 
c                                 console output:
c                                 write to console
               if (fmode.eq.1.and.verbos) then 

                  write (*,1185) (vname(iv(k)),v(iv(k)),k=1,ipot)
                  write (*,1210) gname(ifr(i))

               end if 
c                                 write to file
               if (output) write (n0+i,1200) (v(iv(k)),k=1,ipot),
     *                                nopt(7),(nopt(7),k=1,jbulk) 

            end if 

         end do

         if (fmode.eq.2.and..not.liquid) then 

            write (*,1185) (vname(iv(k)),v(iv(k)),k=1,ipot)

            write (*,'(a)') 'liquid not stable, no phases will be '//
     *                      'fractionated'
         end if 

      else 
c                                 optimization failed
         if (idead.lt.100) then
            write (*,1030) (vname(jv(i)),v(jv(i)), i = 1, ipot)
         else 
            write (*,1035) (vname(jv(i)),v(jv(i)), i = 1, ipot)
         end if
c                                 
         do i = 1, ifrct
c                                 write bad_number to fractionation file
            if (output) write (n0+i,1200) (v(iv(k)),k=1,ipot),nopt(7),
     *                                    (nopt(7),k=1,jbulk)
         end do

      end if 

      warn = .false.
c                                 remove fractionated mass from bulk,
c                                 warn on complete depletion of a component
      do i = 1, jbulk

         cblk(i) = cblk(i) - dcomp(i)

         if (cblk(i).le.nopt(11)) then 

            if (.not.gone(i)) then
               warn = .true.
               gone(i) = .true. 
               write (*,1000) cname(i)
            end if 

            if (cblk(i).lt.nopt(11)) cblk(i) = 0d0

         end if 

      end do

      if (output.and.ifrct.gt.0) 
     *   write (n0,1200) (v(iv(i)), i = 1, ipot), (cblk(k),k=1,jbulk)
 
      if (warn) then 
         do i = 1, jbulk 
            if (.not.gone(i)) write (*,1010) cname(i),cblk(i)
         end do
      end if  

1000  format (3(/,5('*** WARNING ***')),
     *     //,'Fractionation has eliminated ',a,' from the bulk ',
     *        'composition.',//,
     *        'This may destabilize the minimization algorithm. ',
     *        'If excessive',/,'minimization errors occur, restart the',
     *        ' fractionation calculation at the',/,'current point on ',
     *        'the fractionation path with the current molar bulk ', 
     *        'composition:',/)
1010  format (4x,a,2x,g12.6)  
1030  format (/,'optimization failed at:',//,5(1x,a8,'=',g12.6))
1035  format (/,'optimization reported as failed due to EoS problem at:'
     *       ,//,5(1x,a8,'=',g12.6))
1175  format ('current molar bulk composition is:',/,15(1x,g12.6))
1180  format (a,' is stable, but its mass fraction (',f5.3,') is below '
     *      ,'or at the upper fractionation threshold (',f5.3,').')
1185  format (/,'At ',5(a,'=',g12.6,' '))
1190  format ('fractionating ',g12.6,
     *        ' moles of ',a,'; will change bulk by:',/,15(1x,g12.6))
1200  format (21(1x,g12.6))
1210  format (a,' is not stable.')

      end 

      character*22 function lgname (char14,ids)
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character char14*(*)

      integer ids

      integer ikp
      common/ cst61 /ikp(k1)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c----------------------------------------------------------------------
      call getnam (char14,ids)

      if (ids.lt.0) then 

         if (ikp(-ids).gt.0) then 
            lgname = lname(ikp(-ids))
         else if (ids.gt.0) then 
            lgname = char14
         end if 

      else 

         lgname = lname(ids)

      end if 

      end 
