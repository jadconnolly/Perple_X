c Please do not distribute any part of this source.
 
c Copyright (c) 1987-2023 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved. 

      program MC_fit
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------- 
c                                 MC_fit uses the MEEMUM iam flag value
      iam = 2
c                                 perplexwrap.f flags
      getInput = .true.
      sWarn = .false.
c                                 initialization, read files etc.
      call iniprp
c                                 open inversion problem file
      call opnimc
c                                 do the inversion
      call invers

      end

      subroutine opnimc 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      integer ier
c----------------------------------------------------------------------- 
c                                 open inversion problem file
      call mertxt (tfname,prject,'.imc',0)
      open (n8,file=tfname,status='old',iostat=ier)

      if (ier.ne.0) call errdbg 
     *   ('can''t open assemblage composition file: '//tfname)
c----------------------------------------------------------------------- 
c                                 IMC file input Section 1
c                                 read compositional input units
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'molar'.or.key.eq.'mass') then

         if (key.eq.'molar') then
            iopt(2) = 0
         else
            iopt(2) = 1
         end if

      else

         call errdbg ('expecting mass or molar tag, found: '
     *               //key)

      end if
c                                 read problem type flag
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'invptx'.or.key.eq.'invxpt') then

         if (key.eq.'invxpt') then
            invxpt = .true.
         else
            invxpt = .false.
         end if

      else

         call errdbg ('expecting invxpt or invptx tag, found: '
     *               //key)

      end if
c                                 read algorithm flag
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'minim'.or.key.eq.'grid') then

         if (key.eq.'minim') then
            mcgrid =  .false.
         else
            mcgrid = .true.
         end if

         if (mcgrid.and..not.invxpt) call errdbg ('grid search not '
     *                                          //'allowed with invptx')

      else

         call errdbg ('expecting minim or grid tag, found: '
     *               //key)

      end if
c                                 best model criterion (both output, but
c                                 choice selects the scores output to 
c                                 n6). 
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'likelihood'.or.key.eq.'bayes') then

         if (key.eq.'bayes') then
            bayes =  .true.
         else
            bayes = .false.
         end if

      else

         call errdbg ('expecting likelihood or bayes tag, found: '
     *               //key)

      end if
c                                 george's normalization, etc
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'JADC'.or.key.eq.'GRH') then

         if (key.eq.'GRH') then
            GRH =  .true.
         else
            GRH = .false.
         end if

      else

         call errdbg ('expecting GRH or JADC tag, found: '
     *               //key)

      end if
c                                 output only improved results
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'better'.or.key.eq.'all') then

         if (key.eq.'better') then
            better =  .true.
         else
            better = .false.
         end if

      else

         call errdbg ('expecting better or all tag, found: '
     *               //key)

      end if
c                                 output optimization counter
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'vital_sign'.or.key.eq.'quiet') then

         if (key.eq.'vital_sign') then
            vital =  .true.
         else
            vital = .false.
         end if

      else

         call errdbg ('expecting vital_sign or quiet tag, found: '
     *               //key)

      end if
c                                 make new seed for random number generator
c                                 george's normalization, etc
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'seed'.or.key.eq.'no_seed') then

         if (key.eq.'seed') then
            seed =  .true.
         else
            seed = .false.
         end if

      else

         call errdbg ('expecting seed or no_seed tag, found: '
     *               //key)

      end if
c----------------------------------------------------------------------- 
c                                 IMC file input Section 2
c                                 quadratic/linear objective function
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'quadratic'.or.key.eq.'linear') then

         if (key.eq.'linear') then
            grhobj =  .true.
         else
            grhobj = .false.
         end if

      else

         call errdbg ('expecting quadratic or linear tag, found: '
     *               //key)

      end if

      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) oktol

      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) invtol

      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) wcomp

      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) wextra

      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) wmiss
c----------------------------------------------------------------------- 
c                                 IMC file input Section 3
c                                 Monte-Carlo logistical variables:


c     call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

c     ier = index('TtFfGg',key(1:1))
c     if (ier.eq.0) stop '***Bad problem type, check random value'
c     random(1) = (ier-1)/2

c     if (val.ne.' ') then
c        ier = index('TtFf',val(1:1))
c        if (ier.eq.0) stop '***Bad problem type, check random value'
c        ier = (ier-1)/2
c        random(1) = random(1) + 10*ier
c     end if
c                                 number of starting guesses used for each
c                                 Nelder-Meade inverse problem, George also
c                                 somehow converts his input to mtry
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) mtry
c                                 error evaluation loop counter
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) nunc(1)

      if (val.ne.' ') then
         read (val,*) nunc(2)
         if (nunc(2).gt.5) then
            write(*,'(2a,/,a)')
     *     '**warning** 2nd nunc value in MC inversion parameter',
     *     ' is probably too big -','  will slow search; set to 5'
           nunc(2) = 5
         end if
      else
         nunc(2) = 0
      end if
c----------------------------------------------------------------------- 
c                                 IMC file input Section 3
c                                 Nelder-Meade parameters
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) simplx
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) frac
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) conchk
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) jprint
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) iquad
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) kcount

      if (wmiss.gt.oktol .or.wextra.gt.oktol .or.wcomp.gt.oktol) then
         write(*,'(2a,/,a)')
     *     '**warning** oktol value in MC inversion parameters',
     *     ' is probably too small -','  check wcomp, wextra, wmiss'
      end if

      end

      subroutine invers 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k

      logical pdat

      double precision bstx(l2+k5), x(200,l2+k5), sx(l2+k5), ss(l2+k5)
c----------------------------------------------------------------------- 
c                                 best model statistics for error evaluation
      call mertxt (tfname,prject,'.bst',0)
      open (n7,file=tfname)
      call mertxt (tfname,prject,'.bay',0)
      open (n9,file=tfname)
c                                 initialize drand
      if (seed) call random_seed
c                                 get best model, 1st argument sets 
c                                 random perturbation off, 2nd sets 
c                                 output of all sucessful optimizations
      call bstmod (.false.,.true., bstx)

      x(1,1:nparm) = bstx(1:nparm)

      pdat = .true.

      do i = 1, nunc(1)

         call opnimc 
c                                 suppress any grid search at this point
         if (random(1).ge.2) then
            random(1) = random(1)/10
            if (random(1).eq.0 .and. i.eq.1) call random_seed
            pdat = .false.
         end if

         call bstmod (pdat,.false.,bstx)

         x(1+i,1:nparm) = bstx(1:nparm)

         sx(1:nparm) = 0d0
         ss(1:nparm) = 0d0

         do j = 1, 1 + i

            do k = 1, nparm
c                                 sum x
               sx(k) = sx(k) + x(j,k)
            end do

         end do

         do j = 1, 1 + i

            do k = 1, nparm
c                                 sum of squares
               ss(k) = ss(k) + (x(j,k) - sx(k)/(i+1d0))**2

            end do

         end do

         do k = 1, nparm
c                                 standard deviation
            ss(k) = dsqrt(ss(k)/(i))

         end do

         write (n7,'(20(1pg12.6,1x))') ss(1:nparm)

      end do

      close (n7)

      stop

      end 

      subroutine bstmod (randm, n6out, x)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, n, icount, ifault, jcount, nfree, lu, 
     *        j, k, igood, ibest, id, jbest, pnum(l2+k5)

      logical readyn, bad, randm, n6out, improv

      double precision var(l2+k5), objf, bstx(l2+k5), sx(l2+k5), 
     *                 x0(l2+k5), step(l2+k5), bstobj, x(*), 
     *                 bstvar(l2+k5), ssp, bay, bstbx(l2+k5), bstbay

      external readyn, mcobj2

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)
c----------------------------------------------------------------------- 
c                                 flag to make mcobj print extended
c                                 ouput for invptx
      fprint = .false.

      if (n6out) then 
c                                 output file
         call mertxt (tfname,prject,'.out',0)
c                                 in case the user has set print, detach the file
         close (n6)
         open (n6,file=tfname)

         write (n6,*) 'tol/frac/simplx',invtol, frac, simplx
         write (n6,'(80(''-''))')

      end if 

      if (invxpt) then
c                                 parameter inversion:
c                                 read experimental data and inversion candidates
c                                 if randm, then experimental data is perturbed within
c                                 its uncertainty.
         call mcxpt (n6out,randm)

         if (random(1).ge.2) mtry = random(3)-1
c                                 parameter max - min
         n = 0
c                                 compounds
         do i = 1, mccpd
c                                 for each coefficient
            do j = 1, mcpct(i)
               n = n + 1
               plow(n) = cprng(i,j,1)
               pdelta(n) = cprng(i,j,2) - cprng(i,j,1)
               pnum(n) = nint(cprng(i,j,3))
            end do
         end do
c                                 solutions
         do i = 1, mcsol
c                                 for each term
            do j = 1, mctrm(i)
c                                 for each coefficient
               do k = 1, mccoef(i,j) 
                  n = n + 1
                  plow(n) = sprng(i,j,k,1)
                  pdelta(n) = sprng(i,j,k,2) - sprng(i,j,k,1)
                  pnum(n) = nint(sprng(i,j,k,3))
               end do
            end do
         end do

         nfree = n

      else
c                                 ptx inversion:
c                                 read phase compositions
         call mccomp
c                                 set george's "weight"
         xptpt(1,5) = 1d0

         n = ipot + cextra + xptnph(1) - 1

         nfree = ipot + cextra

         do k = 1, n
            if (k.le.ipot) then
               plow(k) = vmin(k)
               pdelta(k) = vmax(k) - vmin(k)
            else
               plow(k) = 0d0
               pdelta(k) = 1d0
            end if

         end do

      end if
c                                 initialize scaled coordinate
      if (mcgrid) then
         sx = 0d0
      else
         sx = 0.5d0
      end if
c                                 unscale step size for search
      step(1:n) = frac*pdelta(1:n)

      nparm = n
      ibest = 0
      igood = 0
      ifault = 0
      bstobj = 1d99
      bstbay = 1d99
      jbest = 0

      do i = 1, mtry
c                                 counter in mcobj2
         optct = 0
c                                 unscale sx
         do j = 1, n
            x(j) = plow(j) + sx(j)*pdelta(j)
         end do
c                                 save starting coordinate
         x0(1:n) = x(1:n)

         if (mcgrid) then
c                                 grid seach:
            call mcobj2 (x,objf,bad)

         else
c                                 Nelder-Meade search:
c                                 initialize icount in case of failure
            icount = 0

            call minim (x, step, n, objf, kcount, jprint, invtol, 
     *                  conchk, iquad, simplx, var, mcobj2, icount, 
     *                  jcount, ifault, oktol, nfree)

         end if
c                                 count and check if model improved
         if (ifault.le.2) call savbst (x,var,objf,bstobj,bstx,bay,
     *                    bstbay,bstbx,bstvar,n,i,ibest,jbest,igood)
         
         improv = i.eq.ibest.or.i.eq.jbest.or..not.better
c                                 print loop
         lu = 6
         consol = .false.

         do
c                                 output some stats
            if (ifault.gt.2.or.(ifault.gt.0.and.objf.gt.oktol)) then
c                                 minim has failed
               write (lu,1020) ifault, icount, igood, i

            else if (mcgrid.and.lu.eq.6) then 
c                                 george's minimal console output
               write (*,'(a,i7,1h/,i7,1x,a,2(1x,1pg12.6),a,$)')
     *                'Try ',i,random(3),'best:',bstobj,bstbay,char(13)

            else if (.not.invxpt.and.improv) then
c                                thermobarometry
               write (lu,'(/80(''-''))')
               write (lu,2000) i
               write (lu,2070) (vname(j),x(j), j = 1, ipot)
c                                write likelihood results
               write (lu,2010) objf
               
               if (i.eq.ibest) then
                  write (lu,2020)
               else
                  write (lu,2030) bstobj, ibest
               end if
c                                 write Bayeseian results
               write (lu,2040) bay
               if (i.eq.jbest) then
                  write (lu,2050)
               else
                  write (lu,2060) bstbay, jbest
               end if
c                                  unnecessary obj call for print output
               fprint = .true.
               call mcobj2 (x,objf,bad)
               fprint = .false.

           else if (invxpt.and.
     *                  (.not.mcgrid.or.(mcgrid.and.improv))) then
c                                   george only outputs improved results         
               write (lu,'(/80(''-''))')

               if (iquad.gt.0.and.jcount.gt.0.and..not.mcgrid) then
                  write (lu,1125) i, igood, icount, jcount, objf, 
     *                           bstobj, ibest
               else
                  write (lu,1120) i, igood, icount, objf, bstobj, ibest
               end if

               write (lu,1130) ssp, bay, bstbay, jbest
               write (lu,1080) sx(1:n)
               write (lu,1085) x0(1:n)
               write (lu,1030) x(1:n)

               if (lu.eq.n6.and.invxpt) then 

                  write (lu,'(/,a,i7,a,/)') 'Scores for try = ',i,
     *                                   ' follow:'

                  do id = 1, mxpt

                     if (xptpt(id,5).eq.1d0) then
                        write (lu,1010)
     *                   id, xptnam(id),' score = ',scores(id)
                     else
                        write (lu,1010)
     *                   id, xptnam(id),' score = ',scores(id),
     *                   xptpt(id,5)
                     end if

                  end do

               end if

               write (lu,'(/80(''-''))')

            end if

            if (lu.eq.n6.or..not.n6out) exit

            lu = n6
            consol = .true.

         end do

         if (.not.mcgrid) then
c                               new random starting point
            do j = 1, n
               call random_number (sx(j))
            end do

         else
c                               new grid search point
            icount = random(2)
            random(2) = icount + 1
            
            do j = 1, n
               sx(j) = dble(mod(icount,pnum(j)))/max(1,pnum(j)-1)
               icount = icount / pnum(j)
            end do

         end if

      end do

      if (bayes) then
         x(1:n) = bstx(1:n)
         objf = bstobj
      else
         x(1:n) = bstbx(1:n)
         objf = bstbay
      end if
c                               write best model to *.bst and *.bay
      write (n7,1000) bstx(1:n), bstobj
      write (n9,1000) bstbx(1:n), bstbay

      if (n6out) close (n6)

1000  format (20(1pg13.6,1x),a)
1010  format (i3,1x,2a,g12.6,:,3h * ,g8.3)
1020  format (80('-'),/,'Search FAILED, ifault = ',i3,', icount = ',
     *       i4,', igood = ',i4,', mtry = ',i7,/,80('-'))
1030  format ('Final coordinates: ',20(1pg13.6,1x))
1050  format (/,'Number of function evaluations: ',i5,', igood = ',i3,/)
1080  format ('Initial normalized coordinates: ',20(1pg12.6,1x))
1085  format ('Initial coordinates: ',20(1pg13.6,1x))
1120  format (/,'Try ',i7,', successes so far = ',i7,/,
     *          'Objective function evaluations this try = ',i5,/,
     *          'Last objective function value this try OBJ = ',g12.6,/,
     *          'Best OBJ so far = ',g12.6,
     *          ' obtained on try ',i7,/)
1125  format (/,'Try ',i7,', successes so far = ',i7,/,
     *          'Objective function evaluations this try = ',i5,/,
     *          ' + ',i4,' evaluations for quadratic surface fitting',/,
     *          'Last objective function value this try OBJ = ',g12.6,/,
     *          'Best OBJ so far = ',g12.6,
     *          ' obtained on try ',i7,/)
1130  format (/,'Scaled parameter SSP = ',g12.6,/,
     *          'Bayes score SSP * OBJF = ',g12.6,/,
     *          'Best Bayes score so far = ',g12.6,
     *          ' obtained on try ',i7,/)

2000  format (/,'Try ',i5,' has converged at:',/)
2010  format (/,'Likelihood score for this Try: ',g12.6)
2020  format ('This is the best likelihood score obtained so far.')
2030  format ('Best likelihood score (',g12.6,
     *        ') so far was obtained on Try',i5)
2040  format (/,'Bayesian score for this Try: ',g12.6)
2050  format ('This is the best Bayesian score obtained so far.')
2060  format ('Best Bayesian score (',g12.6,
     *        ') so far was obtained on Try',i5)
2070  format (29x,a8,' = ',g12.6)

      end


      subroutine mccomp
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for MC thermobarometry
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical randm, bad

      integer ier, nblen, i

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      double precision comp(k5), ecomp(k5)

      external nblen

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c-----------------------------------------------------------------------
c                                to allow invptx to use the same objective 
c                                function and i/o structure as invxpt the 
c                                multi-experiment arrays of invxpt are used
c                                here even though invptx treats only one 
c                                observation (mxpt = 1). 
      mxpt = 1
c                                 number of phase compositions
      cxpt = 0
c                                 unmeasured component counter
      cextra = 0
c                                 look for thermobarometry problem in
c                                 my_project.imc
      call mertxt (tfname,prject,'.imc',0)
      open (n8,file=tfname,status='old',iostat=ier)

      if (ier.ne.0) call errdbg 
     *   ('can''t open assemblage composition file: '//
     *   tfname(1:nblen(tfname)))

      if (random(1).ge.2) call errdbg
     *   ('can''t use grid search option (yet)')
c                                 -------------------------------------
c                                 IMC file input Section 6
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.ne.'begin_assemblage') call errdbg (
     *   'expecting begin_assemblage keyword, found '//key)
c                                 IMC file input Section 6
c                                 read sample name
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.ne.'sample_name') call errdbg (
     *   'expecting sample_name keyword, found '//key)

      read (strg,*) xptnam(mxpt)
c                                 pressure range
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'pressure_range') then

         read (strg1,*) vmin(1), vmax(1)

      else

         call errdbg ('expecting pressure_range tag, found: '
     *               //key)

      end if
c                                 temperature range
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'temperature_range') then

         read (strg1,*) vmin(2), vmax(2)

      else

         call errdbg ('expecting temperature_range tag, found: '
     *               //key)

      end if
c                                 read optional unmeasured component
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'unmeasured_component') then

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.ne.'begin_comp') call errdbg ('expecting begin_comp'//
     *                            ' tag, found :'//key)

         call gtcomp (comp,ecomp,randm,'comp',bad,.false.)

         if (bad) call errdbg ('error reading unmeasured_component '//
     *                         'composition')

         do i = 1, kbulk
            xptc(cxpt+i) = comp(i)
            xpte(cxpt+i) = ecomp(i)
         end do
c                                 pointer to the composition of phase nph in expt mexpt 
         xptptr(mxpt,k5) = cxpt
c                                 variable counter increment
         cextra = 1
c                                 increment composition pointer
         cxpt = cxpt + kbulk

      else if (key.eq.'phase_name') then

         backspace (n8)

      else

         call errdbg ('expecting unmeasured_component or phas_name tag,'
     *               //'found: '//key)

      end if
c                                 read the assemblage data
      call gtassmb (comp,ecomp,randm,'assemblage',bad)

      close (n8)

      end

      double precision function pertrb (num,err)
c-----------------------------------------------------------------------
c a function to add random error to num
c-----------------------------------------------------------------------
      implicit none

      double precision num, err, x
c
      call random_number (x)

      if (err.le.0d0) then

         pertrb = num

      else 

         pertrb = num + 2d0*err*(x - 0.5d0)

         if (pertrb.lt.0d0) then
            pertrb = 0d0
            write (*,*) 'oinky poinky'
         end if

      end if

      end 

      subroutine mcxpt (n6out, randm)
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for MC inversion of exptal
c data. parameters:
c     n6out - echo parameter output to n6
c     randm - randomly perturb bulk compositions before inversion
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, bad, randm, n6out, done, eof

      integer i, j, k, ier, id, nblen

      integer*8 ngrid

      double precision err, pertrb, comp(k5), ecomp(k5)

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40, str(3)*8

      external pertrb, nblen

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)
c-----------------------------------------------------------------------
c                                 read solutions and compounds to be perturbed
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.ne.'begin_phase_free_param') 
     *                     call errdbg ('invalid data, last read '//key)

      mccpd = 0
      mcsol = 0

      ok = .false.
      bad = .false.
c                                 loop to read inversion phases
      do

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'end_phase_free_paramet') exit
c                                 phase name
         read (key,'(a)') tname
c                                 check name
         call matchj (tname,id)

         if (id.eq.0) then

            write(*,1040) tname(1:nblen(tname))

            do
               call redcd1(n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
               if (ier.ne.0 .or. key.eq.'end_list') exit
            end do

            cycle

         end if

         if (id.lt.0) then

            ok = .true.
            mccpd = mccpd + 1
            mcid(mccpd) = -id
            mcpct(mccpd) = 0

            do
c                                 read free parameters for compounds
               call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

               if (key.eq.'end_list') exit

               if (key.ne.'parameter') call errdbg ('expecting '//
     *                    'parameter tag for '//tname//' found '//key)

               mcpct(mccpd) = mcpct(mccpd) + 1

               if (val.eq.'a') then 
                  mcpid(mccpd,mcpct(mccpd)) = 1
               else if (val.eq.'b') then 
                  mcpid(mccpd,mcpct(mccpd)) = 2
               else if (val.eq.'c') then 
                  mcpid(mccpd,mcpct(mccpd)) = 3
               else if (val.eq.'K') then 
                  mcpid(mccpd,mcpct(mccpd)) = 4
               else if (val.eq.'K''') then 
                  mcpid(mccpd,mcpct(mccpd)) = 5
               else if (val.eq.'V0') then 
                  mcpid(mccpd,mcpct(mccpd)) = 6
               else if (val.eq.'HJT' .or. val.eq.'HJB') then 
                  if (lct(-id).gt.1) then
                     write(*,1020) tname,'multiple'
                     call errpau
                  else if (lct(-id).ne.1) then
                     write(*,1020) tname,'no'
                     call errpau
                  else if (ltyp(-id).ne.8) then
                     write(*,1010) tname,'has no H&J transition'
                     call errpau
                  end if
                  if (val.eq.'HJT') then
                     mcpid(mccpd,mcpct(mccpd)) = 7
                  else
                     mcpid(mccpd,mcpct(mccpd)) = 8
                  end if
               else 
                  call errdbg ('invalid parameter name for ' //
     *                         tname(1:nblen(tname)) // ': "' //
     *                         val(1:nblen(val)) // '"')
               end if
c                                 read parameter range and grid search count
               call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

               if (key.ne.'range') call errdbg ('expecting '//
     *                    'range tag for '//tname//' found '//key)
               read (strg1,*) cprng(mccpd,mcpct(mccpd),1),
     *                        cprng(mccpd,mcpct(mccpd),2)
               read (strg1,*,iostat=ier) err,err,
     *                        cprng(mccpd,mcpct(mccpd),3)
               if (ier.ne.0)  cprng(mccpd,mcpct(mccpd),3) = 2d0
               if (cprng(mccpd,mcpct(mccpd),3).le.0d0) 
     *                        cprng(mccpd,mcpct(mccpd),3) = 2d0

            end do

         else if (id.gt.0) then
c                                 read free parameters for solutions
c                                 for Redlich-Kister models, the syntax is
c                                 parameter 1 = w0 value,
c                                 parameter 2 = wT value,
c                                 parameter 3 = wP (even though this would
c                                    be the index of the Brosh wP0 parameter,
c                                    it is fudged later on).

            if (.not.bad) mcsol = mcsol + 1
            mcids(mcsol) = id
c                                initialize term counter
            mctrm(mcsol) = 0
c                                initialize coefficient counter
            mccoef(mcsol,1:m1) = 0

            done = .false.

            do
c                                 read free excess terms for solutions
               call readcd (n8,ier,.true.)

               if (chars(1).eq.'e') exit

               if (chars(1).ne.'W') call errdbg ('invalid excess term'
     *                                           //' for '//tname)
               call redtrm (bad)

               if (bad) then
c                                 invalid endmember, read to end_list
c                                 and cycle
                  do
                     call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,
     *                            strg1)
                     if (key.eq.'end_list') exit
                  end do

                  cycle

               end if
c                                 at this point we've read a legitimate term
c                                 now check that it exists in the model:
c                                 ------------------------------------------
c                                 indices are now in isub(1,1:iord) look for 
c                                 the term in the real model
c                                 checks vary depending on type: extyp(id)
c                                 = 0 margules; = 1 redlich-kister; = 2 van laar
               do j = 1, jterm(id)

                  if (extyp(id).ne.1 .and. rko(j,id).ne.iord) cycle

                  ok = .true.

                  do i = 1, iord
c                                 now check the endmembers match, assume same
c                                 ordering as in the solution model
                     if (extyp(id).eq.1) then
c                                 redlich-kister
                        ok = ok .and. jend(id,2+i).eq.isub(1,i)
                     else
c                                 margules, try match on jth term
                        ok = ok .and.
     *                     isub(1,i).eq.jend(id,2+jsub(i,j,id))
                     end if

                  end do

                  if (ok) exit

               end do
              
               if (.not.ok) then
c                                term does not exist in solution model
                  print*,(chars(i),i=1,com),'term does not exist'
                  cycle 
               end if
c                                the term exists:
c                                increment term counter
               mctrm(mcsol) = mctrm(mcsol) + 1
c                                save the term pointer
               mcj(mcsol,mctrm(mcsol)) = j
c                                now find which coefficients of the
c                                term are free
               do 
c                                 loop to read free coefficients in term
                  call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,
     *                         strg1)
c                                 key may be either "parameter"
c                                 or the name of a new term "W(..." or
c                                 "Wk(..." or "end_list"

                  if (key.eq.'end_list') then
                     if (mccoef(mcsol,mctrm(mcsol)).eq.0) call errdbg(
     *                                 'a term with no coefficients?')
                     done = .true.

                     exit

                  else if (key(1:1).eq.'W') then
c                                 backspace the input to use redtrm
                     backspace(n8)

                     exit

                  else if (key.ne.'parameter') then 

                      call errdbg ('invalid coefficient or parameter '//
     *                             'specification for '//tname)

                  end if
c                                 count the coefficient
                  mccoef(mcsol,mctrm(mcsol)) = 
     *                   mccoef(mcsol,mctrm(mcsol)) + 1
c                                 read the coefficient range and grid numbers
                  read (strg1,*) 
     *           mccoid(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol))),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),1),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),2)
                  read (strg1,*,iostat=ier) i,err,err,
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
                  if (ier.ne.0) 
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
     *               = 2d0
                  if (
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
     *               .le.0d0)
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
     *               = 2d0

               end do
c                                 flag to make gall recalculate
c                                 staticredt compounds
               mcflag(id) = .true.

               if (done) exit

            end do

         end if

      end do
c                                 holy schmoly! if you followed that mess
c                                 you deserve a medal. count the parameters
c                                 and echo them if requested
      if (n6out) write(n6,*) 'Inverted parameters:'
      nparm = 0 
c                                 first for endmembers
      do i = 1, mccpd
         if (n6out) then
c                                 echo parameters to output
            id = mcid(i)
            do j = 1, mcpct(i)
               k = mcpid(i,j)
               if (k .eq. 1) then
                  key = 'delta-G'
               else if (k .eq. 2) then
                  key = 'delta-S'
               else if (k .eq. 3) then
                  key = 'delta-V'
               else if (k .eq. 4) then
                  key = 'K'
               else if (k .eq. 5) then
                  key = 'K'''
               else if (k .eq. 6) then
                  key = 'V0'
               else if (k .eq. 7) then
                  key = 'H&J-T'
               else if (k .eq. 8) then
                  key = 'H&J-B'
               else
                  call errdbg('**internal error** MCXPT encoding')
               end if
               if (random(1).lt.2) then
                  write(n6,1050) nparm+j,names(id),key(1:nblen(key)),
     *               cprng(i,j,1), cprng(i,j,2)
               else
                  write(n6,1050) nparm+j,names(id),key(1:nblen(key)),
     *               cprng(i,j,1), cprng(i,j,2), int(cprng(i,j,3))
               end if
            end do
            
         end if
c                                 not so hard to count
         nparm = nparm + mcpct(i)

      end do
c                                 now solutions
      do i = 1, mcsol
c                                 for each term, classify type
         id = mcids(i)
         if (extyp(id).eq. 0) then
            val = 'W'
         else if (extyp(id).eq. 1) then
            val = 'Wk'
         else if (extyp(id).eq. 2) then
            call errdbg('**internal error** no van Laar solns yet')
            val = 'WV'
         else
            call errdbg('**internal error** MCXPT extyp(.) invalid')
         end if

         do j = 1, mctrm(i) 
c                                 grunt work to count and name the term
            if (n6out) then
c                                 echo parameters to output
               if (extyp(id).eq.1) then
c                                 Redlich-Kister names always pairwise
                  k = jend(id,2+1)
                  strg = names(k)(1:nblen(names(k))) //
     *               ' ' // names(jend(id,2+2))
               else if (extyp(id).eq.0) then
c                                 Margules names depend on order of solution
                  strg = ' '
                  do k = 1, rko(j,id)
                     strg(1+nblen(strg):) = ' ' //
     *               names(jend(id,2+jsub(k,j,id)))
                  end do
                  strg = strg(2:)
               else if (extyp(id).eq. 2) then
c                 nothing yet for van Laar - who knows how to input them?
               end if

               do k = 1, mccoef(i,j)
c                                 which parameters are getting varied?
                  if (random(1).lt.2) then
                     write(n6,1060) nparm+k, fname(id),
     *                  val(1:nblen(val)), strg(1:nblen(strg)), j,
     *                  mccoid(i,j,k),
     *                  sprng(i,j,k,1), sprng(i,j,k,2)
                  else
                     write(n6,1060) nparm+k, fname(id),
     *                  val(1:nblen(val)), strg(1:nblen(strg)), j,
     *                  mccoid(i,j,k),
     *                  sprng(i,j,k,1), sprng(i,j,k,2),
     *                  int(sprng(i,j,k,3))
                  end if
               end do

            end if

            nparm = nparm +  mccoef(i,j)

         end do
      end do

      if (n6out) write (n6,'(80(''-''))')

      if (nparm.eq.0) call errdbg ('no free parameters! no free lunch!')

c                                 if grid search, set up counts
      if (random(1).ge.2) then
         ngrid = 1
         do i = 1, mccpd
            do j = 1, mcpct(i)
               ngrid = ngrid*nint(cprng(i,j,3))
            end do
         end do
         do i = 1, mcsol
            do j = 1, mctrm(i)
               do k = 1, mccoef(i,j) 
                  ngrid = ngrid*nint(sprng(i,j,k,3))
               end do
            end do
         end do
         if (ngrid.ge.2d0**32) 
     *      call errdbg ('Excessive grid search combinations, coarsen')
         random(2) = 0
         random(3) = ngrid
         write(*,*) 'Grid search will explore ',ngrid,' combinations.'
      end if
c                                 number of xpts
      mxpt = 0
c                                 number of phase compositions
      cxpt = 0
c                                 flag for observations with 
c                                 extraneous components
      bad = .false.
c                                 loop to read observation data
      do

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 assume EOF on error
         if (ier.ne.0) exit

         if (key.ne.'begin_exp') call errdbg (
     *                           'expecting begin_exp, found '//key)
c                                an experiment, only increment
c                                counter if previous result was ok
         if (.not.bad) mxpt = mxpt + 1

         if (mxpt.gt.l11) call errdbg ('too many expts, increase l11')
c                                 reset bad, bad will be set to true
c                                 if a bulk or phase composition includes
c                                 a missing component. 
         bad = .false.
c                                initialize solution counters
         msolct(mxpt,1:isoct) = 0
c                                read name
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         xptnam(mxpt) = key
c                               read expt p,t, save any error
         call rdstrg (n8,i,str,eof)

         read (str(1),*) xptpt(mxpt,1)
         read (str(2),*) err
         xptpt(mxpt,3) = err

         if (randm) xptpt(mxpt,1) = pertrb (xptpt(mxpt,1),err)

         call rdstrg (n8,i,str,eof)

         read (str(1),*) xptpt(mxpt,2)
         read (str(2),*) err
         xptpt(mxpt,4) = err
c                               weight
         xptpt(mxpt,5) = 1d0

         if (randm) xptpt(mxpt,2) = pertrb (xptpt(mxpt,2),err)
c                               read bulk composition
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.ne.'begin_bulk') call errdbg (
     *                           'expecting begin_bulk, found '//key)
c                                 read bulk composition, normalize, and
c                                 purturb
         call gtcomp (comp,ecomp,randm,'bulk',bad,.true.)

         do i = 1, icomp
            xptblk(mxpt,i) = comp(i)
         end do
c                                 read optional weight
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (ier.eq.0 .and. key.eq.'weight') then
            read (val, *, iostat=ier) xptpt(mxpt,5)
            if (ier.ne.0 .or. xptpt(mxpt,5).lt.0d0)
     *         call errdbg ('bad weight: '//val)
         else
            backspace(n8)
         end if
c                                 read phase names, modes
c                                 and compositions
         call gtassmb (comp,ecomp,randm,'exp',bad)

         if (bad) then

            call mertxt (tfname,prject,'.dat',0)
            write (*,1000) xptnam(mxpt), tfname(1:nblen(tfname))
            if (grh) write (*,1030)

         end if
c                                 next experiment
      end do

      close (n8)

1000  format (/,'warning ver502** observation: ',a,/,'has been rejecte',
     *          'd because it includes a component not specified in: ',
     *          a,//,80('-'))
1010  format (/,a,1x,a)
1020  format (/,a,1x,'has ',a,' magnetic transitions')
1030  format ('or because it has > 1 composition with no uncertainty')
1040  format (/,'**warning ver502** parameter ',a,' skipped - ',
     *          'not a compound or solution',//,80('-'))
1050  format (i2,2(1x,a),2(1x,1pg12.5),1x,i3)
1060  format (i2,1x,a,1x,a,1h(,a,1h),i2,' parameter',i2,2(1x,1pg12.5),
     *        1x,i3)

      end

      subroutine gtassmb (comp,ecomp,randm,tag,bad)
c-----------------------------------------------------------------------
c read phase data between between begin_\\tag/end_\\tag keywords. 

c      data includes:

c           phase names
c           phase modes and error (optional)
c           phase compositions and error for solutions (via gtcomp).

c if invxpt returns bad if expt includes a component not specified 
c    in the problem definition file or if grh and missing compositional
c    uncertainties

c if invptx no return on bad data.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, bad, randm, absent(k5)

      integer i, ids, nph, ier

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40, tag*(*)

      double precision comp(*), ecomp(*), pertrb

      external pertrb

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c----------------------------------------------------------------------
      nph = 0
      ok = .true.
      absent = .true.
      msolct(mxpt,1:isoct) = 0

      do
c                                 now read phase data
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'end_'//tag) then

            ok = .false.
            exit

         end if
c                                 phase name
         read (strg,'(a)') tname

         nph = nph + 1
c                                 check name
         call matchj (tname,ids)
c                                 instead of making this an error 
c                                 could just set bad to reject the xpt
         if (ids.eq.0) call errdbg ('no such entity as: '//tname)
c                                 all clear, save id
         xptids(mxpt,nph) = ids
c                                 look for optional modal data, only 
c                                 invptx, so 1d arrays
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'phase_mode') then
c                                 modal data
            read (nval1,*) pmode(nph)
            read (nval2,*) emode(nph)

         else
c                                 no modal data, initialize
            backspace (n8)
            pmode(nph) = -1d0
            emode(nph) = 1d0

         end if
c                                 if compound don't read composition
         if (ids.lt.0) then
            xptnph(mxpt) = nph
            cycle
         end if
c                                 counters in case same solution > 1 time
         if (msolct(mxpt,ids).eq.0) then
c                                 first case
            msolct(mxpt,ids) = 1
            msloc(mxpt,msolct(mxpt,ids)) = nph

         else

            msolct(mxpt,ids) = msolct(mxpt,ids) + 1
            msloc(mxpt,msolct(mxpt,ids)) = nph

         end if
c                                 get solution composition
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.ne.'begin_comp') call errdbg (
     *                           'expecting begin_comp, found: '//key)

         call gtcomp (comp,ecomp,randm,'comp',bad,.true.)

         do i = 1, kbulk
            
            if (comp(i).ne.0d0) absent(i) = .false.
c                                 should loop if cextra can be > 1
            if (absent(i).and.cextra.gt.0) then 
               if (xptc(xptptr(1,k5)+i).ne.0) absent(i) = .false. 
            end if

            xptc(cxpt+i) = comp(i)
            xpte(cxpt+i) = ecomp(i)
            
         end do
c                                 pointer to the composition of phase nph in expt mexpt 
         xptptr(mxpt,nph) = cxpt

         xptnph(mxpt) = nph
c                                 increment composition pointer
         if (.not.bad) cxpt = cxpt + kbulk

      end do

      if (invxpt) then



      else 

         if (nph.lt.2) call errdbg ('input must specify > 1 phase')

         ok = .true.

         do i = 1, kbulk
            if (absent(i)) then
               write (*,1000) cname(i)
               ok = .false.
            end if
         end do

         write (*,'(/)')

         if (lopt(56).and..not.ok) call wrnstp

      end if

1000  format (/,'**warning ver777** the problem definition file specif',
     *        'ies a component (',a,') that',/,'is absent from the obs',
     *        'erved phase assemblage specified in the *.imc file. To',/
     *       ,'avoid bad practice delete the absent component.')

      end 

      subroutine gtcomp (comp,ecomp,randm,tag,bad,norm)
c-----------------------------------------------------------------------
c read compositional data and uncertainities between begin_\\tag/end_\\tag 
c keywords. converts mass units to molar units if iopt(2) = 1.
c if randm perturn composition within uncertainty.
c returns bad if the composition includes a component not specified in
c the problem definition file. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, bad, randm, norm

      integer i, j, ier

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40, tag*(*)

      double precision tot, comp(*), ecomp(*), pertrb

      external pertrb
c----------------------------------------------------------------------
      comp(1:kbulk) = 0d0
      ecomp(1:kbulk) = 1d0

      do

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'end_'//tag) exit

         ok = .false.
c                                 should test for component name replication?
         do i = 1, kbulk
            if (key.eq.cname(i)) then
               ok = .true.
               exit
            end if
         end do

         if (.not.ok) then

            if (invxpt) then 
c                                 set bad to reject the experiment but continue
c                                 reading data in case there are valid experiments
               bad = .true.

            else
c                                  thermobarometry, no point in continuing.
                call errdbg ('phase contains a component:'//key//' not '
     *          //'specified in the problem defintiion file')

            end if

         end if

         read (strg1,*) comp(i), ecomp(i)
c                                 convert mass input to molar units
         if (iopt(2).eq.1) then

            comp(i) = comp(i)/atwt(i)
            ecomp(i) = ecomp(i)/atwt(i)

         end if

         if (randm) comp(i) = pertrb (comp(i),ecomp(i))

      end do

      tot = 0d0
c                                 george's 0 counter
      j = 0

      do i = 1, kbulk
         if (ecomp(i).eq.0d0) j = j + 1
         tot = tot + comp(i)
      end do
      
      if (norm) then 

         if (grh) then 
c                                 george assumes normalized input, but has 
c                                 this suspect test. a phase with 
c                                 fewer components than the system will
c                                 fail, and what's with the > 1?
            if (j.gt.1 .and. j.ne.kbulk) bad = .true.
            xpterr(mxpt) = j

         else

            do i = 1, kbulk
c                                 normalize composition and error for
c                                 scoring:
               comp(i) = comp(i)/tot
               ecomp(i) = ecomp(i)/tot

            end do

         end if

      end if

      end

      subroutine redtrm (bad)
c----------------------------------------------------------------------
c redtrm - read excess function term of the form 

c        W(name-name-...)

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer ibeg, jend, ier, iscan, imax, iscnlt, id

      character name*8

      external iscnlt, iscan

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------
      ier = 0
      bad = .false.
c                                 find expansion type
      if (chars(2).eq.'k'.or.chars(2).eq.'K') then
         xtyp = 1
      else
         xtyp = 0
      end if
c                                 find brackets
      ibeg = iscan (1,com,'(') + 1
      imax = iscan (1,com,')') - 1

      if (ibeg.gt.com.or.imax.gt.com) then
         bad = .true.
         return 
       end if 
c                                 data found
      iord = 0

      do while (ibeg.lt.imax)

         call readnm (ibeg,jend,imax,ier,name)
         if (ier.ne.0) call errdbg ('wroink! excess term')

         iord = iord + 1
         if (iord.gt.m2) call error (49,wg(1,1),m2,name)

         call matchj (name,id)

         if (id.eq.0) then
            write (*,1000) chars(1:com)
            write (*,1010) name
            bad = .true.
            return
         end if

         isub(1,iord) = -id

      end do

1000  format (/,'Warning: term ',60a)
1010  format ('includes an endmember (',a,') not currently in the ',
     *        'model, the term will be rejected',/)

      end

      subroutine mcsetv (x)
c-----------------------------------------------------------------------
c a subprogram to set potential variables from scaled coordinates during
c MC inversion.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision x(*)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)
c-----------------------------------------------------------------------
c                                 convert the scaled potential variables back
c                                 to the normal variables:
      do i = 1, ipot
         v(jv(i)) = vmin(jv(i)) + x(i) * (vmax(jv(i)) - vmin(jv(i)))
      end do

      end 

      subroutine mcsetb (x)
c-----------------------------------------------------------------------
c set bulk composition for MC thermobarometry, assumes normalized 
c compositions
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x(*)

      integer i, j, k

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
c-----------------------------------------------------------------------
      cblk = 0d0
      ctotal = 0d0

      do i = 1, xptnph(1)
c                                 composition pointer
         cxpt = xptptr(1,i)
c                                 compositional var pointer
         k = ipot + cextra + i

         if (i.lt.xptnph(1)) then 
            ctotal = ctotal + x(k)
         else
            x(k) = 1d0 - ctotal
         end if

         do j = 1, kbulk
            cblk(j) = cblk(j) + x(k) * xptc(cxpt+j)
         end do

      end do
c                                 modify cblk here to change the 
c                                 composition before minimization.
      if (cextra.gt.0) then
         
         k = ipot + cextra
         
         cxpt = xptptr(1,k5)

         do j = 1, kbulk
            cblk(j) = cblk(j) + x(k) * xptc(cxpt+j)
         end do

      end if

      ctotal = 1d0

      b(1:kbulk) = cblk(1:kbulk) 

      end

      subroutine mcstb2 (id,ip,it,ner)
c-----------------------------------------------------------------------
c set bulk composition for MC experiment inversion
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id, ip, it, ner

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c-----------------------------------------------------------------------
      ctotal = 0d0

      do i = 1, kbulk
         cblk(i) = xptblk(id,i)
      end do
c                                 get total moles to compute mole fractions 
      do i = 1, kbulk
         ctotal = ctotal + cblk(i)
      end do

      do i = 1, kbulk
         b(i) = cblk(i) / ctotal
      end do
c                                 set p-t
      v(1) = xptpt(id,1)
      v(2) = xptpt(id,2)
c                                 add error contribution for p & t
c                                 ***WARNING***
c                                 assumes v(3) & v(4) are unused potentials;
c                                 only p (1) and t (2) valid.
c                                 xptpt(.,3:4) are p & t uncertainties, not
c                                 potentials.
c                                 ***WARNING***
      if (ip.ne.0) v(1) = v(1) + ip * xptpt(id,3)/ner
      if (it.ne.0) v(2) = v(2) + it * xptpt(id,4)/ner

      end

      subroutine minim (p,step,nop,func,max,jprint,stopcr,nloop,iquad,
     *                  simp,var,functn,neval,meval,ifault,oktol,nfree)
c----------------------------------------------------------------------
c     a program for function minimization using the simplex method.
c     the minimum found will often be a local, not a global, minimum.
c
c     for details, see nelder & mead, the computer journal, january 1965
c
c     programmed by d.e.shaw,
c     csiro, division of mathematics & statistics
c     p.o. box 218, lindfield, n.s.w. 2070
c
c     with amendments by r.w.m.wedderburn
c     rothamsted experimental station
c     harpenden, hertfordshire, england
c
c     further amended by alan miller,
c     csiro, division of mathematics & statistics
c     private bag 10, clayton, vic. 3168
c
c     arguments:

c     p()     = input, starting values of parameters
c               output, final values of parameters
c     step()  = input, initial step sizes
c     nop     = input, no. of parameters, incl. any to be held fixed
c     func    = output, the function value corresponding to the final
c               parameter values
c     max     = input, the maximum no. of function evaluations allowed
c     iprint  = input, print control parameter
c                     < 0 no printing
c                     = 0 printing of parameter values and the function
c                         value after initial evidence of convergence.
c                     > 0 as for iprint = 0 plus progress reports after
c                         every iprint evaluations, plus printing for the
c                         initial simplex.
c     stopcr  = input, stopping criterion
c     nloop   = input, the stopping rule is applied after every nloop
c               function evaluations.
c     iquad   = input, = 1 if the fitting of a quadratic surface is required
c                      = 0 if not
c     simp    = input, criterion for expanding the simplex to overcome
c               rounding errors before fitting the quadratic surface.
c     var()   = output, contains the diagonal elements of the inverse of
c               the information matrix.
c     functn  = input, name of the user's subroutine - arguments (p,func)
c               which returns the function value for a given set of
c               parameter values in array p. must be declared external in 
c               the calling program.
c     ifault  = output, = 0 for successful termination
c                       = 1 if maximum no. of function evaluations exceeded
c                       = 2 if information matrix is not +ve semi-definite
c                       = 3 if nop < 1
c                       = 4 if nloop < 1

c       advice on usage:

c       if the function minimized can be expected to be smooth in the vicinity
c       of the minimum, users are strongly urged to use the quadratic-surface
c       fitting option.   this is the only satisfactory way of testing that the
c       minimum has been found.   the value of simp should be set to at least
c       1000 times the rounding error in calculating the fitted function.
c       e.g. in double precision on a micro- or mini-computer with about 16
c       decimal digit representation of floating-point numbers, the rounding
c       errors in calculating the objective function may be of the order of
c       1.e-12 say in a particular case.   a suitable value for simp would then
c       be 1.e-08.   however, if numerical integration is required in the
c       calculation of the objective function, it may only be accurate to say
c       1.e-05 and an appropriate value for simp would be about 0.1.
c       if the fitted quadratic surface is not +ve definite (and the function
c       should be smooth in the vicinity of the minimum), it probably means
c       that the search terminated prematurely and you have not found the
c       minimum.

c       n.b. p, step and var (if iquad = 1) must have dimension at least nop
c            in the calling program.
c       the dimensions below are for a maximum of 20 parameters.
c       the dimension of bmat should be at least nop*(nop+1)/2.
c
c       latest revision - 11 august 1991
c----------------------------------------------------------------------
c                                 original code used implicit typing
      implicit none

      external functn

      logical bad

      integer nop, lout, jprint, ifault, nloop, nap, loop, iflag, i, 
     *        irow, j, np1, imax, imin, iquad, neval, max, i1, i2, j1, 
     *        k, l, ii, ij, nullty, irank, jj, ijk, meval, nfree

      double precision zero, one, two, three, half, a, b, c, func, 
     *                 fnap, fnp1, savemn, test, simp, a0, stopcr,
     *                 rmax, ymin, hmax, hmin, hstar, hstst, hstd, 
     *                 hmean, p(nop),step(nop),var(nop), oktol

      double precision g(21,20),h(21),pbar(20),pstar(20),pstst(20),
     *                 aval(20),bmat(210),pmin(20),vc(210),temp(20)

      data zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/, half/0.5d0/
c
c     a = reflection coefficient, b = contraction coefficient, and
c     c = expansion coefficient.
c
      data a,b,c/1.d0, 0.5d0, 2.d0/
c
c     set lout = logical unit no. for output
c
      data lout/6/
c
c     if progress reports have been requested, print heading
c
      if(jprint.gt.0) write(lout,1000) jprint
 1000 format(' progress report every',i4,' function evaluations'/,
     1  ' eval.  func.',15x,'parameter values')

c     check input arguments

      ifault=0
      if(nop.le.0) ifault=3
      if(nloop.le.0) ifault=4
      if(ifault.ne.0) return

c     set nap = no. of parameters to be varied, i.e. with step.ne.0

      nap=0
      loop=0
      iflag=0

      do i=1,nop
        if(step(i).ne.zero) nap=nap+1
      end do

c     if nap = 0 evaluate function at the starting point and return

      if(nap.gt.0) go to 30

      call functn(p,func,bad)
      if (bad) ifault = 99
      return

c     set up the initial simplex

   30 do i=1,nop
         g(1,i)=p(i)
      end do

      irow=2

      do i=1,nop

        if(step(i).eq.zero) cycle

        do j=1,nop
           g(irow,j)=p(j)
        end do

        g(irow,i)=p(i)+step(i)
        irow=irow+1

      end do

      np1 = nap + 1
      neval = 0
      meval = 0

      do i=1,np1

        do j=1,nop
           p(j)=g(i,j)
        end do

        call functn(p,h(i),bad)

        if (bad) then 
           ifault = 99
           return
        end if

        neval=neval+1

        if(jprint.le.0) cycle

        write(lout,1010) neval,h(i),(p(j),j=1,nop)
 1010   format(/i4, 2x, g12.5, 2x, 5g12.5, 3(/20x, 5g12.5))

      end do

c     start of main cycle.

c     find max. & min. values for current simplex (hmax & hmin).

  100 loop=loop+1
      imax=1
      imin=1
      hmax=h(1)
      hmin=h(1)

      do i=2,np1

        if(h(i).le.hmax) go to 110
        imax=i
        hmax=h(i)

        cycle

  110   if(h(i).ge.hmin) cycle

        imin=i
        hmin=h(i)

      end do

c     find the centroid of the vertices other than p(imax)

      do i=1,nop
         pbar(i)=zero
      end do

      do i=1,np1
         if(i.eq.imax) cycle
         do j=1,nop
            pbar(j)=pbar(j)+g(i,j)
         end do 
      end do

      do j=1,nop
         fnap = nap
         pbar(j)=pbar(j)/fnap
      end do

c     reflect maximum through pbar to pstar,
c     hstar = function value at pstar.

      do i = 1, nop

         pstar(i)=a*(pbar(i)-g(imax,i))+pbar(i)

         if ((i.gt.nfree.and.pstar(i).gt.1d0).or.
     *       (i.gt.nfree.and.pstar(i).lt.0d0)) then 
            
            if (pstar(i).gt.1d0.and.pbar(i).lt.1d0) then
               pstar(i) = 1d0
            else if (pstar(i).lt.0d0.and.pbar(i).gt.0d0) then
               pstar(i) = 0d0               
            else 
              ifault = 99
              return
            end if 
            
         end if

      end do

      call functn(pstar,hstar,bad)

      neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

      if(jprint.le.0) go to 180
      if(mod(neval,jprint).eq.0) write(lout,1010) neval,hstar,
     1  (pstar(j),j=1,nop)

c     if hstar < hmin, reflect pbar through pstar,
c     hstst = function value at pstst.

  180 if (hstar.ge.hmin) go to 220

      do i = 1, nop

         pstst(i)=c*(pstar(i)-pbar(i))+pbar(i)

         if ((i.gt.nfree.and.pstst(i).gt.1d0).or.
     *       (i.gt.nfree.and.pstst(i).lt.0d0)) then 

            if (pstst(i).gt.1d0.and.pbar(i).lt.1d0) then
               pstst(i) = 1d0
            else if (pstst(i).lt.0d0.and.pbar(i).gt.0d0) then
               pstst(i) = 0d0               
            else 
              ifault = 99
              return
            end if 

         end if

      end do

      call functn(pstst,hstst,bad)

      neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

      if(jprint.le.0) go to 200
      if(mod(neval,jprint).eq.0) write(lout,1010) neval,hstst,
     1  (pstst(j),j=1,nop)
c
c     if hstst < hmin replace current maximum point by pstst and
c     hmax by hstst, then test for convergence.
c
  200 if(hstst.ge.hmin) go to 320

      do i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstst(i)
      end do

      h(imax)=hstst
      go to 340

c     hstar is not < hmin.
c     test whether it is < function value at some point other than
c     p(imax).   if it is replace p(imax) by pstar & hmax by hstar.

  220 do i=1,np1
         if(i.eq.imax) cycle
         if(hstar.lt.h(i)) go to 320
      end do

c     hstar > all function values except possibly hmax.
c     if hstar <= hmax, replace p(imax) by pstar & hmax by hstar.

      if(hstar.gt.hmax) go to 260

      do i=1,nop
         if(step(i).ne.zero) g(imax,i)=pstar(i)
      end do

      hmax=hstar
      h(imax)=hstar

c     contracted step to the point pstst,
c     hstst = function value at pstst.

  260 do i=1,nop
         pstst(i)=b*g(imax,i) + (1.d0-b)*pbar(i)
      end do

      call functn(pstst,hstst,bad)

      neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

      if(jprint.le.0) go to 280
      if(mod(neval,jprint).eq.0) write(lout,1010) neval,hstst,
     1  (pstst(j),j=1,nop)

c     if hstst < hmax replace p(imax) by pstst & hmax by hstst.

  280 if(hstst.gt.hmax) go to 300

      do i=1,nop
         if(step(i).ne.zero) g(imax,i)=pstst(i)
      end do

      h(imax)=hstst

      go to 340

c     hstst > hmax.
c     shrink the simplex by replacing each point, other than the current
c     minimum, by a point mid-way between its current position and the
c     minimum.

  300 do i=1,np1

         if(i.eq.imin) cycle

         do j=1,nop
            if(step(j).ne.zero) g(i,j)=(g(i,j)+g(imin,j))*half
            p(j)=g(i,j)
         end do

         call functn(p,h(i),bad)

         neval = neval+1

         if (bad) then 
            ifault = 99
            return
         end if

         if(jprint.le.0) cycle
         if(mod(neval,jprint).eq.0) write(lout,1010) neval,h(i),
     *              (p(j),j=1,nop)

      end do

      go to 340

c     replace maximum point by pstar & h(imax) by hstar.

  320 do i=1,nop
         if(step(i).ne.zero) g(imax,i)=pstar(i)
      end do

      h(imax)=hstar

c     if loop = nloop test for convergence, otherwise repeat main cycle.

  340 if(loop.lt.nloop) go to 100

c     calculate mean & standard deviation of function values for the
c     current simplex.

      hstd=zero
      hmean=zero

      do i=1,np1
         hmean=hmean+h(i)
      end do

      fnp1 = np1
      hmean=hmean/fnp1

      do i=1,np1
         hstd=hstd+(h(i)-hmean)**2
      end do

      hstd=sqrt(hstd/float(np1))

c     if the rms > stopcr, set iflag & loop to zero and go to the
c     start of the main cycle again.

      if(hstd.le.stopcr.or.neval.gt.max) go to 410

      iflag=0
      loop=0
      go to 100

c     find the centroid of the current simplex and the function value there.

  410 do i=1,nop

         if(step(i).eq.zero) cycle
         p(i)=zero

         do j=1,np1
            p(i)=p(i)+g(j,i)
         end do

         fnp1 = np1
         p(i)=p(i)/fnp1

      end do

      call functn(p,func,bad)

      neval=neval+1

      if (bad) then 
         ifault = 99
         return
      end if

      if (func.gt.oktol) then
         write (*,'(a,1x,g12.6)') 'Aborting, bad objf: ',h(i)
         ifault = 2
         return
      end if

      if(jprint.le.0) go to 390
      if(mod(neval,jprint).eq.0) write(lout,1010) neval,func,
     1  (p(j),j=1,nop)

c     test whether the no. of function values allowed, max, has been
c     overrun; if so, exit with ifault = 1.

  390 if(neval.le.max) go to 420
      ifault=1
      if(jprint.lt.0) return
      write(lout,1020) max
 1020 format(' no. of function evaluations exceeds ',i5)
      write(lout,1030) hstd
 1030 format(' rms of function values of last simplex =',g14.6)
      write(lout,1040)(p(i),i=1,nop)
 1040 format(' centroid of last simplex =',4(/1x,1p6g13.5))
      write(lout,1050) func
 1050 format(' function value at centroid =',g14.6)
      return

c     convergence criterion satisfied.
c     if iflag = 0, set iflag & save hmean.
c     if iflag = 1 & change in hmean <= stopcr then search is complete.

  420 if(jprint.lt.0) go to 430
      write(lout,1060)
 1060 format(/' evidence of convergence')
      write(lout,1040)(p(i),i=1,nop)
      write(lout,1050) func
  430 if(iflag.gt.0) go to 450
      iflag=1
  440 savemn=hmean
      loop=0
      go to 100
  450 if(abs(savemn-hmean).ge.stopcr) go to 440
      if(jprint.lt.0) go to 460
      write(lout,1070) neval
 1070 format(//' minimum found after ',i5,' function evaluations')
      write(lout,1080)(p(i),i=1,nop)
 1080 format(' minimum at',4(/1x,1p6g13.6))
      write(lout,1090) func
 1090 format(' function value at minimum =',g14.6)
  460 if(iquad.le.0) return
c-------------------------------------------------------------------

c     quadratic surface fitting

      if(jprint.ge.0) write(lout,1110)
 1110 format(/' quadratic surface fitting about supposed minimum'/)

c     expand the final simplex, if necessary, to overcome rounding
c     errors.

      do i=1,np1

  470   test=abs(h(i)-func)

        if(test.ge.simp) cycle

        if (func.gt.oktol) then
           write (*,'(a,1x,g12.6)') 'Aborting, bad objf: ',func
           ifault = 2
           return
        end if

        do j=1,nop

           if(step(j).ne.zero) g(i,j)=(g(i,j)-p(j))+g(i,j)

           pstst(j)=g(i,j)

           if (isnan(g(i,j))) then

              write (*,'(a,1x,g12.6)') 'Aborting, bad coord: ',g(i,j)
              ifault = 2
              return

           end if

        end do

        call functn(pstst,h(i),bad)

        if (bad) then 
           ifault = 99
           return
        end if

        if (h(i).gt.oktol) then
           write (*,'(a,1x,g12.6)') 'Aborting, bad objf: ',h(i)
           ifault = 2
           return
        end if

        meval = meval+1
c                                 quick fix: this segment goes into an 
c                                 infinite loop if g(i,j) doesn't change
        if (meval.gt.max) then
           write (*,*) 'Aborting, infinite loop during surface fitting'
           ifault = 2
           return
        end if 

        go to 470

      end do

c     function values are calculated at an additional nap points.

      do i = 1, nap

         i1=i+1

         do j=1,nop
            pstar(j)=(g(1,j)+g(i1,j))*half
         end do

         call functn(pstar,aval(i),bad)

         if (bad) then 
            ifault = 99
            return
         end if

         if (aval(i).gt.oktol) then
            write (*,'(a,1x,g12.6)') 'Aborting, bad objf: ',h(i)
            ifault = 2
            return
         end if

         meval= meval + 1

      end do
c
c     the matrix of estimated second derivatives is calculated and its
c     lower triangle stored in bmat.
c
      a0=h(1)
      do 540 i=1,nap
        i1=i-1
        i2=i+1
        if(i1.lt.1) go to 540
        do 530 j=1,i1
          j1=j+1
          do 520 k=1,nop
  520     pstst(k)=(g(i2,k)+g(j1,k))*half
          call functn(pstst,hstst,bad)

        if (bad) then 
           ifault = 99
           return
        end if

        if (hstst.gt.oktol) then
           write (*,'(a,1x,g12.6)') 'Aborting, bad objf: ',h(i)
           ifault = 2
           return
        end if

          meval = meval+1
          l=i*(i-1)/2+j
          bmat(l)=two*(hstst+a0-aval(i)-aval(j))
  530   continue
  540 continue
      l=0
      do 550 i=1,nap
        i1=i+1
        l=l+i
        bmat(l)=two*(h(i1)+a0-two*aval(i))
  550 continue

c     the vector of estimated first derivatives is calculated and
c     stored in aval.

      do 560 i=1,nap
        i1=i+1
        aval(i)=two*aval(i)-(h(i1)+three*a0)*half
  560 continue
c
c     the matrix q of nelder & mead is calculated and stored in g.
c
      do 570 i=1,nop
  570 pmin(i)=g(1,i)
      do 580 i=1,nap
        i1=i+1
        do 580 j=1,nop
        g(i1,j)=g(i1,j)-g(1,j)
  580 continue
      do 590 i=1,nap
        i1=i+1
        do 590 j=1,nop
          g(i,j)=g(i1,j)
  590 continue
c
c     invert bmat
c
      call syminv(bmat,nap,bmat,temp,nullty,ifault,rmax)
      if(ifault.ne.0) go to 600
      irank=nap-nullty
      go to 610
  600 if(jprint.ge.0) write(lout,1120)
 1120 format(/' matrix of estimated second derivatives not +ve defn.'/
     1  ' minimum probably not found'/)
      ifault=2
      return

c     bmat*a/2 is calculated and stored in h.

  610 do 650 i=1,nap
        h(i)=zero
        do 640 j=1,nap
          if(j.gt.i) go to 620
          l=i*(i-1)/2+j
          go to 630
  620     l=j*(j-1)/2+i
  630     h(i)=h(i)+bmat(l)*aval(j)
  640   continue
  650 continue

c     find the position, pmin, & value, ymin, of the minimum of the
c     quadratic.

      ymin=zero
      do 660 i=1,nap
  660 ymin=ymin+h(i)*aval(i)
      ymin=a0-ymin
      do 670 i=1,nop
        pstst(i)=zero
        do 670 j=1,nap
  670 pstst(i)=pstst(i)+h(j)*g(j,i)
      do 680 i=1,nop
  680 pmin(i)=pmin(i)-pstst(i)
      if(jprint.lt.0) go to 682
      write(lout,1130) ymin,(pmin(i),i=1,nop)
 1130 format(' minimum of quadratic surface =',g14.6,' at',
     1  4(/1x,6g13.5))
      write(lout,1150)
 1150 format(' if this differs by much from the minimum estimated',
     1  1x,'from the minimization,'/
     2  ' the minimum may be false &/or the information matrix may be',
     3  1x,'inaccurate'/)

c     calculate true function value at the minimum of the quadratic.

  682 meval = meval + 1
      call functn(pmin, hstar,bad)

        if (bad) then 
           ifault = 99
           return
        end if

c     if hstar < func, replace search minimum with quadratic minimum.

      if (hstar .ge. func) go to 690
      func = hstar
      do 684 i = 1, nop
  684 p(i) = pmin(i)
      write(lout, 1140) func
 1140 format(' true func. value at minimum of quadratic = ', g14.6/)
c
c     q*bmat*q'/2 is calculated & its lower triangle stored in vc
c
  690 do 760 i=1,nop
        do 730 j=1,nap
          h(j)=zero
          do 720 k=1,nap
            if(k.gt.j) go to 700
            l=j*(j-1)/2+k
            go to 710
  700       l=k*(k-1)/2+j
  710       h(j)=h(j)+bmat(l)*g(k,i)*half
  720     continue
  730   continue
        do 750 j=i,nop
          l=j*(j-1)/2+i
          vc(l)=zero
          do 740 k=1,nap
  740     vc(l)=vc(l)+h(k)*g(k,j)
  750   continue
  760 continue

c     the diagonal elements of vc are copied into var.

      j=0
      do 770 i=1,nop
        j=j+i
        var(i)=vc(j)
  770    continue
      if(jprint.lt.0) return
      write(lout,1160) irank
 1160 format(' rank of information matrix =',i3/
     1  ' generalized inverse of information matrix:-')
      ijk=1
      go to 880
  790 continue
      write(lout,1170)
 1170 format(/' if the function minimized was -log(likelihood),'/
     1  ' this is the covariance matrix of the parameters'/
     2  ' if the function was a sum of squares of residuals'/
     3  ' this matrix must be multiplied by twice the estimated',
     4  1x'residual variance'/' to obtain the covariance matrix.'/)
      call syminv(vc,nap,bmat,temp,nullty,ifault,rmax)
c
c     bmat now contains the information matrix
c
      write(lout,1190)
 1190 format(' information matrix:-'/)
      ijk=3
      go to 880
c
c     calculate correlations of parameter estimates, put into vc.
c
  800 ijk=2
      ii=0
      ij=0
      do 840 i=1,nop
        ii=ii+i

        if(vc(ii).gt.zero) then
          vc(ii)=one/sqrt(vc(ii))
        else 
          vc(ii)=zero
        end if

        jj=0
        do 830 j=1,i-1
          jj=jj+j
          ij=ij+1
          vc(ij)=vc(ij)*vc(ii)*vc(jj)
  830   continue
        ij=ij+1
  840 continue
      write(lout,1200)
 1200 format(/' correlation matrix:-')
      ii=0
      do 850 i=1,nop
        ii=ii+i
        if(vc(ii).ne.zero) vc(ii)=one
  850 continue
      go to 880
  860 write(lout,1210) meval
 1210 format(/' a further',i4,' function evaluations have been used'/)
      return

c     pseudo-subroutine to print vc if ijk = 1 or 2, or
c     bmat if ijk = 3.

  880 l=1
  890 if(l.gt.nop) go to (790,860,800),ijk
      ii=l*(l-1)/2
      do 910 i=l,nop
        i1=ii+l
        ii=ii+i
        i2=min(ii,i1+5)
        if(ijk.eq.3) go to 900
        write(lout,1230)(vc(j),j=i1,i2)
        go to 910
  900   write(lout,1230)(bmat(j),j=i1,i2)
  910 continue
 1230 format(1x,6g13.5)
      write(lout,1240)
 1240 format(/)
      l=l+6
      go to 890
      end

      subroutine syminv(a,n,c,w,nullty,ifault,rmax)
c
c     algorithm as7, applied statistics, vol.17, 1968.
c
c     arguments:-
c     a()     = input, the symmetric matrix to be inverted, stored in
c               lower triangular form
c     n       = input, order of the matrix
c     c()     = output, the inverse of a (a generalized inverse if c is
c               singular), also stored in lower triangular.
c               c and a may occupy the same locations.
c     w()     = workspace, dimension at least n.
c     nullty  = output, the rank deficiency of a.
c     ifault  = output, error indicator
c                     = 1 if n < 1
c                     = 2 if a is not +ve semi-definite
c                     = 0 otherwise
c     rmax    = output, approximate bound on the accuracy of the diagonal
c               elements of c.  e.g. if rmax = 1.e-04 then the diagonal
c               elements of c will be accurate to about 4 dec. digits.
c
c     latest revision - 18 october 1985
c
c*************************************************************************
c
      implicit double precision (a-h, o-z)

      integer n, nrow, nullty, ifault, nn, irow, ndiag, l, i, icol, 
     *        jcol, mdiag, j, k

      double precision zero, one, rmax, x
      double precision a(*),c(*),w(n)

      data zero/0.d0/, one/1.d0/
c
      nrow=n
      ifault=1
      if(nrow.le.0) go to 100
      ifault=0
c
c     cholesky factorization of a, result in c
c
      call chola(a,nrow,c,nullty,ifault,rmax,w)
      if(ifault.ne.0) go to 100
c
c     invert c & form the product (cinv)'*cinv, where cinv is the inverse
c     of c, row by row starting with the last row.
c     irow = the row number, ndiag = location of last element in the row.
c
      nn=nrow*(nrow+1)/2
      irow=nrow
      ndiag=nn
   10 if(c(ndiag).eq.zero) go to 60
      l=ndiag
      do 20 i=irow,nrow
        w(i)=c(l)
        l=l+i
   20 continue
      icol=nrow
      jcol=nn
      mdiag=nn
   30 l=jcol
      x=zero
      if(icol.eq.irow) x=one/w(irow)
      k=nrow
   40 if(k.eq.irow) go to 50
      x=x-w(k)*c(l)
      k=k-1
      l=l-1
      if(l.gt.mdiag) l=l-k+1
      go to 40
   50 c(l)=x/w(irow)
      if(icol.eq.irow) go to 80
      mdiag=mdiag-icol
      icol=icol-1
      jcol=jcol-1
      go to 30
c
c     special case, zero diagonal element.
c
   60 l=ndiag
      do 70 j=irow,nrow
        c(l)=zero
        l=l+j
   70 continue
c
c      end of row.
c
   80 ndiag=ndiag-irow
      irow=irow-1
      if(irow.ne.0) go to 10
  100 return
      end

      subroutine chola(a, n, u, nullty, ifault, rmax, r)
c
c     algorithm as6, applied statistics, vol.17, 1968, with
c     modifications by a.j.miller
c
c     arguments:-
c     a()     = input, a +ve definite matrix stored in lower-triangular
c               form.
c     n       = input, the order of a
c     u()     = output, a lower triangular matrix such that u*u' = a.
c               a & u may occupy the same locations.
c     nullty  = output, the rank deficiency of a.
c     ifault  = output, error indicator
c                     = 1 if n < 1
c                     = 2 if a is not +ve semi-definite
c                     = 0 otherwise
c     rmax    = output, an estimate of the relative accuracy of the
c               diagonal elements of u.
c     r()     = output, array containing bounds on the relative accuracy
c               of each diagonal element of u.
c
c     latest revision - 18 october 1985
c
c*************************************************************************
c
      implicit double precision (a-h, o-z)

      integer n, irow, l, icol, m, ifault, nullty, i, j, k

      double precision a(*),u(*),r(n), w, eta, zero, five, rmax, rsq
c
c     eta should be set equal to the smallest +ve value such that
c     1.0 + eta is calculated as being greater than 1.0 in the accuracy
c     being used.
c
      data eta/1.d-16/, zero/0.d0/, five/5.d0/
c
      ifault=1
      if(n.le.0) go to 100
      ifault=2
      nullty=0
      rmax=eta
      r(1)=eta
      j=1
      k=0
c
c     factorize column by column, icol = column no.
c
      do 80 icol=1,n
        l=0
c
c     irow = row number within column icol
c
        do 40 irow=1,icol
          k=k+1
          w=a(k)
          if(irow.eq.icol) rsq=(w*eta)**2
          m=j
          do 10 i=1,irow
            l=l+1
            if(i.eq.irow) go to 20
            w=w-u(l)*u(m)
            if(irow.eq.icol) rsq=rsq+(u(l)**2*r(i))**2
            m=m+1
   10     continue
   20     if(irow.eq.icol) go to 50
          if(u(l).eq.zero) go to 30
          u(k)=w/u(l)
          go to 40
   30     u(k)=zero
          if(abs(w).gt.abs(rmax*a(k))) go to 100
   40   continue
c
c     end of row, estimate relative accuracy of diagonal element.
c
   50   rsq=sqrt(rsq)
        if(abs(w).le.five*rsq) go to 60
        if(w.lt.zero) go to 100
        u(k)=sqrt(w)
        r(i)=rsq/w
        if(r(i).gt.rmax) rmax=r(i)
        go to 70
   60   u(k)=zero
        nullty=nullty+1
   70   j=j+icol
   80 continue
      ifault=0
c
  100 return
      end

      double precision function score (kd,id,j)
c-----------------------------------------------------------------------
c a function to evaluate the distance between oberved and candidate 
c compositions for MC, parameters are:
c     kd - phase indicator
c     id - experiment number
c     j  - phase number in experiment
c-----------------------------------------------------------------------
      include 'perplex_parameters.h'

      integer id, kd, j, k, l, m

      double precision total, term
c-----------------------------------------------------------------------
      total = 0d0
      score = 0d0
      m = 0
c                                 normalization depends on whether any
c                                 compositional uncertainty given
      do l = 1, kbulk

         total = total + pcomp(l,kd)
         if (xpte(xptptr(id,j)+l).eq.0d0) m = l

      end do

      if (m.eq.0) then
c                                 calculate residual
         do l = 1, kbulk

            k = xptptr(id,j)+l

            term = dabs(pcomp(l,kd)/total - xptc(k)) / xpte(k)

            if (.not.grhobj.or.term.lt.1d0) term = term**2

            score = score + term

         end do

      else if (m.eq.1) then
c                                 only one has no uncertainty - normalize
c                                 to it, calculate weighted residual
         total = xptc(m) / pcomp(m,kd)

         do l = 1, kbulk

            if (l.eq.m) cycle

            k = xptptr(id,j)+l

            term = dabs((pcomp(l,kd)*total - xptc(k)) / xpte(k))

            if (.not.grhobj.or.term.lt.1d0) term = term**2

            score = score + term

         end do

      end if

      end

      subroutine x2ther (x)
c-----------------------------------------------------------------------
c a to map the MC search coordinates into thermodynamic arrays
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, n, id, ids, jt

      double precision x(*)

      character name*14

      integer eos
      common/ cst303 /eos(k10)
c-----------------------------------------------------------------------
      n = 0
c                                 compounds
      do i = 1, mccpd
c                                 for each coefficient
         do j = 1, mcpct(i)
            n = n + 1
            id = mcid(i)
            ids = mcpid(i,j)
            if (ids.le.3) then
               mdqf(make(mcid(i)),ids) = x(n)
            else if (ids.eq.4) then
c                                 K - modified compound data directly
               thermo(16,id) = x(n)
            else if (ids.eq.5) then
c                                 K' - modified compound data directly
               thermo(18,id) = x(n)
               if (lopt(4)) thermo(21,id) = abs(x(n))
            else if (ids.eq.6) then
c                                 V0 - modified compound data directly
               thermo(3,id) = x(n)
            else if (ids.eq.7) then
c                                 Tc - modify compound data directly
               therlm(1,1,id) = x(n)
            else if (ids.eq.8) then
c                                 B - modify compound data directly
               therlm(2,1,id) = x(n)
            end if
         end do
      end do
c                                 solutions, this is only for margules
c                                 for rk (extyp(mcids(i)).eq.1), the 
c                                 coefficients are in wkl(m16,m17,m18,h9)
      do i = 1, mcsol
c                                 solution id:
         ids = mcids(i)
c                                 for each term
         do j = 1, mctrm(i)
c                                 term id
            jt = mcj(i,j)
c                                 for each coefficient
            do k = 1, mccoef(i,j)
               n = n + 1
               if (extyp(ids).eq.1) then
c                                 for rk coefficients, funny handling of wP
                  if (k.lt.3) then
                     wkl(k,j,jt,ids) = x(n)
                  else
                     wkl(6,j,jt,ids) = x(n)
                  end if
               else if (extyp(ids).eq.0) then
c                                 for margules coefficients
                  wgl(mccoid(ids,j,k),jt,ids) = x(n)
               else
                  call getnam (name,ids)
                  call errdbg ('unimplemented solution type '//
     *                             'for '//name)
                  
               end if
            end do

         end do
      end do

      end

      subroutine mcobj2 (x,obj,bad)
c-----------------------------------------------------------------------
c a subprogram to evaluate objective function for MC data inversion
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer id, i, neg

      double precision x(*), obj, xptscr, value

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      external xptscr
c-----------------------------------------------------------------------
      do i = 1, nparm
         if (abs(x(i)).gt.1d99) then
            write (*,*) 'woo',x(1:nparm)
            bad = .true.
            return
         end if
      end do

      neg = nunc(2)

      if (invxpt) then
c                                 map the search coordinates to thermodynamic
c                                 parameters
         call x2ther (x)

      else
c                                 set the potentials
         do i = 1, ipot
            v(i) = x(i)
         end do
c                                 set the bulk composition
         call mcsetb (x)

      end if

      obj = 0d0
c                                 loop through observations
      do id = 1, mxpt
c                                 set p, t, bulk
         if (invxpt) call mcstb2 (id,0,0,neg)
c                                 do the optimization and via getloc compute
c                                 system and derivative properties; these
c                                 computations are costly and can be streamlined
c                                 for specific applications.
         call meemum (bad)

         if (bad) then 
            write (*,*) 'oink'
            scores(id) = 1d99
            obj = 1d99
            cycle
         end if
c                                 get score for this p, t, x
         value = xptscr(id)
c                                 if a missing or extra phase, widen search
c                                 within p, t uncertainty range
         if (
     *       neg.gt.0 .and.
     *       (xptpt(id,3).ne.0d0 .or. xptpt(id,4).ne.0d0) .and.
     *       value.gt.wmiss .or. value.gt.wextra/xptnph(id)
     *   ) then
            do i = 1, neg
               if (xptpt(id,3).ne.0d0) then
                  call mcstb2 (id, -i, 0, neg)
                  call meemum (bad)
                  if (.not.bad) value = min(value,xptscr(id))
                  call mcstb2 (id, +i, 0, neg)
                  call meemum (bad)
                  if (.not.bad) value = min(value,xptscr(id))
               end if
               if (xptpt(id,4).ne.0d0) then
                  call mcstb2 (id, 0, -i, neg)
                  call meemum (bad)
                  if (.not.bad) value = min(value,xptscr(id))
                  call mcstb2 (id, 0, +i, neg)
                  call meemum (bad)
                  if (.not.bad) value = min(value,xptscr(id))
               end if
            end do
         end if

         scores(id) = value

         obj = obj + value*xptpt(id,5)

      end do

      end

      double precision function xptscr (id)
c-----------------------------------------------------------------------
c function to score a calculation's agreement with an experiment for MC
c data inversion
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, imout(k5), imin(k5), used(k5), jmiss(k5)

      character bdname(k5)*14

      integer id, jd, ids, i, j, k, kct(k5), ksol(k5,k5), ibest, mpred,
     *        idpred(k5), idextr(k5), jmin(k5), mextra, jdbest, kd, lu

      double precision score, best, res, mode, tot, extra, comp, miss

      external score

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus
c---------------------------------------------------------------------
c                                 compute the observation objective function
      kct = 0
      imout = .true.
      imin = .false.
      used = .false.
      jmiss = .true.
      ksol = 0
      jmin = 0
      optct = optct + 1

      do i = 1, ntot

         do j = 1, xptnph(id)

            if (kkp(i).eq.xptids(id,j)) then

               imout(i) = .false.

               kct(j) = kct(j) + 1
               ksol(j,kct(j)) = i

            end if

         end do

      end do

      extra = 0d0
      comp = 0d0
      miss = 0d0
      mpred = 0
c                                 first extraneous phases:
      do i = 1, ntot
c                                 penalty is wextra * mass_fraction^2
         if (imout(i)) extra = extra + wextra * 
     *                 (props(17,i)*props(16,i)/psys(17))**2

      end do
c                                 potential target phases:
      do j = 1, xptnph(id)

         if (kct(j).eq.0) cycle

         jd = ksol(j,1)
         ids = kkp(jd)

         if (kct(j).eq.1.and.ids.lt.0) then
c                                 a compound, just count
            mpred = mpred + 1
            imin(jd) = .true.
            jmiss(j) = .false.
            jmin(mpred) = j

         else if (kct(j).eq.1.and.msolct(id,ids).eq.1) then
c                                 found solution and no ambiguity
            mpred = mpred + 1
            imin(jd) = .true.
            jmiss(j) = .false.
            jmin(mpred) = j
c                                 compute and add residual
            comp = comp + wcomp * score (jd,id,j)

         else if (kct(j).gt.1) then
c                                 if there are multiple appearances of a
c                                 solution in an assemblage, use the
c                                 best fit
            best = 1d99
            ok = .false.

            do i = 1, kct(j)

               jd = ksol(j,i)

               if (used(jd)) cycle

               res = score (jd,id,j)

               if (res.lt.best) then
                  best = res
                  ibest = i
                  jdbest = jd
                  ok = .true.
               end if

            end do

            if (ok) then

               mpred = mpred + 1

               used(ibest) = .true.
               imin(jdbest) = .true.
               jmiss(j) = .false.
               jmin(mpred) = j

               comp = comp + wcomp * best

            end if

         end if

      end do
c                                 score remaining extraneous phaes
      do j = 1, xptnph(id)

         if (kct(j).lt.2) cycle

         do i = 1, kct(j)

            jd = ksol(j,i)

            if (used(jd)) cycle
c                                 set used to avoid double counting
            used(jd) = .true.
c                                 residual is wextra * mass_fraction^2
            extra = extra + wextra * 
     *                   (props(17,jd)*props(16,jd)/psys(17))**2

         end do

      end do
c                                 missing phase residual
      miss = wmiss * (1d0 - dble(mpred)/xptnph(id))**2
c                                 accumulate scores
      xptscr = miss + extra + comp

      if (.not.invxpt.and.fprint) then
c                                 output optimal P-T and compositions
c                                 for inverse thermo-barometry
c                                 ------------------------------------
c                                 locate predicted and extra phases
         mpred = 0
         mextra = 0

         do i = 1, ntot
            if (imin(i)) then
               mpred = mpred + 1
               idpred(mpred) = i
            else if (imout(i)) then
               mextra = mextra + 1
               idextr(mextra) = i
            end if
         end do
         
         if (consol) then
            lu = 6
         else
            lu = n6
         end if
         
         write (lu,1100) comp, extra, miss

         if (mpred.gt.0) then 
c                                 phases observed and predicted
            write (lu,1000)
            write (lu,1020) (cname(k), k = 1, kbulk)

            do j = 1, mpred

               jd = idpred(j)
               kd = jmin(j)

               mode = props(1,jd)*props(16,jd)/psys(1)*1d2

               tot = 0d0

               do k = 1, kbulk
                  tot = tot + pcomp(k,jd)
               end do

               write (lu,'(a)') pname(jd)
               write (lu,1030) 'predicted*', mode, 
     *                            (pcomp(k,jd)/tot, k = 1, kbulk)
               write (lu,1035) 'observed* ', (xptc(xptptr(1,kd)+k), 
     *                                                     k = 1, kbulk)
               write (lu,1045) 'std resid ', (
     *               (xptc(xptptr(1,kd)+k) - pcomp(k,jd)/tot) 
     *              / xpte(xptptr(1,kd)+k), k = 1, kbulk)

            end do

         end if

         if (mextra.gt.0) then 
c                                 phases predicted but not observed
            write (lu,1040)
            write (lu,1020) (cname(k), k = 1, kbulk)

            do j = 1, mextra

               jd = idextr(j)

               mode = props(1,jd)*props(16,jd)/psys(1)*1d2

               tot = 0d0

               do k = 1, kbulk
                  tot = tot + pcomp(k,jd)
               end do

               write (lu,1050) pname(jd), mode, 
     *                            (pcomp(k,jd)/tot, k = 1, kbulk)

            end do

         end if
c                                 phases observed but not predicted
         mextra = 0

         do j = 1, xptnph(id)

            if (jmiss(j)) then
               mextra = mextra + 1
               call getnam(bdname(mextra),xptids(1,j))
            end if

         end do

         if (mextra.gt.0) write (lu,1060) (bdname(j), j = 1, mextra)

         write (lu,1080)

         do j = 1, kbulk
            write (lu,1090) cname(j), cblk(j), mu(j)
         end do

         write (lu,1010)

      end if
         
      if (.not.fprint)  write (*,'(a,g12.6,1x,i6,60x,a,$)')
     *                  'Score: ', miss + extra + comp, optct, char(13)

1000  format (/,'The following observed phases are predicted:',/)
1010  format (/,'*normalized molar units',/,80('-'))
1020  format (17x,'vol %    ',20(1x,a,2x))
1030  format (2x,a,3x,f7.3,3x,20(f7.4,1x))
1035  format (2x,a,13x,20(f7.4,1x))
1040  format (/,'The following predicted phases are not observed:',/)
1045  format (2x,a,13x,20(f7.3,1x))
1050  format (a,1x,f7.3,3x,20(f7.4,1x))
1060  format (/,'The following observed phases are not predicted:',//,
     *           10(2x,a))
1080  format (/,5x,'       Effective Bulk*      ',
c                12345678901234567890123456789
     *            '  Chemical Potentials (J/mol)')
1090  format (a5,10x,f8.3,20x,g13.6)
1100  format (/,'Components of the Likelihood Score',/,5x,
     *          'Predicted compositions:      ',g14.8,/,5x,
     *          'Extraneous predicted phases: ',g14.8,/,5x,
     *          'Missed observed phases: ',g14.8,/)
1120  format (29x,a8,' = ',g12.6)
      end

      subroutine savbst (x,var,objf,bstobj,bstx,bay,bstbay,bstbx,bstvar,
     *                   n,i,ibest,jbest,igood)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, n, j, igood, ibest, jbest

      double precision var(*), objf, x(*), bstobj, bstvar(*), ssp,
     *                 bay, bstbx(*), bstbay, bstx(*)
c----------------------------------------------------------------------
      igood = igood + 1
c                                 compute ss of parameter deviations
      ssp = 1d0

            do j = 1, n
c                                 center "bayesian" score in interval
               if (pdelta(j).ne.0d0) then
                  ssp = ssp +
     *               ((2*(x(j) - plow(j) - pdelta(j)/2))/pdelta(j))**2
               end if

            end do 
c                                 best "bayesian" score
            bay = ssp * objf

            if (bay.lt.bstbay) then

               jbest = i
               bstbay = bay
               bstbx(1:n) = x(1:n)

            end if 
c                                 save max likelihood result
            if (objf.lt.bstobj) then

               ibest = i
               bstobj = objf
               bstx(1:n) = x(1:n)
               bstvar(1:n) = var(1:n)

            end if

      end