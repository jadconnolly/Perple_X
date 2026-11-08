c Please do not distribute any part of this source.
 
c Copyright (c) 1987-2023 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved. 

      program mc_fit
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j 

      integer iam
      common/ cst4 /iam

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------- 
c                                 mc_fit uses the MEEMUM iam flag value
      iam = 2
c                                 until it get's it's own iam flag use
c                                 mcfit to control i/o etc.
      mcfit = .true.
c                                 perplexwrap.f flags
      getInput = .true.
      sWarn = .false.
c                                 read input normal thermo files, etc
      call iniprp
c                                 -------------------------------------
c                                 at this point it's possible to override 
c                                 options/variable configuration etc

c                                 a saturated phase constraint
      if (ifct.eq.1) then 
c                                 it'll have to be pure, remove the variable 
c                                 from the ipot/iv list
         j = 0

         do i = 1, ipot

            if (iv(i).eq.3) then 
               v(iv(i)) = vmin(iv(i))
            else 
               j = j + 1
               iv(j) = iv(i)
            end if

         end do

         ipot = j

      end if
c                                 set to molar output regardless of 
c                                 option file, mass input units in the imc
c                                 file allowed via lmass flag
      iopt(2) = 0
c                                 -------------------------------------
c                                 open inversion problem file
      call opnimc

      if (invxpt.and.grh) then

         call errdbg ('grh option not implemented')

      else if (invxpt) then

         call invrs1

      else 
c                                inverse thermobarometry
         call invers

      end if 

      end

      subroutine opnimc
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, ier, noerr(k10), nblen
      
      external nblen
c----------------------------------------------------------------------- 
c                                 open inversion problem file
      call mertxt (tfname,prject,'.imc',0)
      open (n8,file=tfname,status='old',iostat=ier)

      if (ier.ne.0) call errdbg 
     *   ('can''t open the inversion input file: '//tfname)
c                                 rewind in case, n8 hasn't been closed
c                                 by mccomp (e.g., thermo-only error analysis)
      rewind (n8)
c                                 first execution, read options
      call mcopts
c                                 consitency checks on options:
c                                 -------------------------------------
c                                 uncertainty analysis error sources, sets jnvrnd :
c                                 invunc = 1 = > perturb all data
c                                 invunc = 2 = > perturb analytical data only
c                                 invunc = 3 = > perturb thermodynamic data only
c
c                                 works in concert with invprt set in subroutine invers
c                                 invprt = F = > not uncertainty analysis(no pertrubations)
c                                 invprt = T = > uncertainty analysis
      if (.not.uncrty.and.(invunc.eq.1.or.invunc.eq.3)) then

         write (*,1020) n2name(1:nblen(n2name))

         call wrnstp

         invunc = 2

      else if (mcfrst.and.(invunc.eq.1.or.invunc.eq.3)) then

         j = 0

         do i = 1, imkend

            if (make(i).ne.0) then
               cycle
            else if (mkptr(i).eq.0) then
               if (deltah(i).eq.0d0) then
                  j = j + 1
                  noerr(j) = i
               end if
            else if (deltah(mkptr(i)).eq.0d0) then
               j = j + 1
               noerr(j) = i
            end if

         end do

         if (j.gt.0) then
            write (*,1000)
            write (*,1010) (names(noerr(i)), i = 1,j)
            write (*,'(/)')
            call wrnstp
         end if

      end if
c                                 override ptry if using best central 
c                                 model for perturbation analysis initial 
c                                 condition.
      if (.not.newstt) ptry = 1
c                                 george's conditions:
      if (wmiss.gt.oktol .or.wextra.gt.oktol .or.wcomp.gt.oktol) then
         write(*,'(2a,/,a)')
     *     '**warning** oktol value in MC inversion parameters',
     *     ' is probably too small -','  check wcomp, wextra, wmiss'
      end if

1000  format (/,'**warning ver503** the following entities have no ',
     *          'associated thermodynamic uncertainty:',/)
1010  format (8(a,1x))
1020  format (/,'**warning ver503** the mc_fit uncertainties option is',
     *       ' set to evaluate thermodynamic',/,'uncertainties, but ', 
     *       a,' does not include thermodynamic uncertainties.',/)

      end

      subroutine invrs1
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k

      double precision bstx(l2+k5), x(200,l2+k5), sx(l2+k5), ss(l2+k5), 
     *                 bstobj, oldobj
c----------------------------------------------------------------------- 
c                                 best model statistics for error evaluation
      call mertxt (tfname,prject,'.bst',0)
      open (n7,file=tfname)
c                                 initialize drand
      if (seed) call random_seed
c                                 get best model, 1st argument sets 
c                                 random perturbation off, 2nd sets 
c                                 output of all sucessful 
      call bstmd1 (.false.,.true., bstx)

      x(1,1:nparm) = bstx(1:nparm)

      oprt = .false.

      do i = 1, nunc(1)

         call opnimc

         call bstmd1 (.true.,.false., bstx)

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

         write (n7,'(20(g12.6,1x))') ss(1:nparm)

      end do

      close (n7)

      stop

1010  format (/,'Try: ',i4,/)
1030  format ('Final coordinates: ',5(g12.6,1x))
1050  format (/,'Number of function evaluations: ',i5,', igood = ',i3,/)
1080  format ('Initial normalized coordinates: ',5(g12.6,1x))
1085  format ('Initial coordinates: ',5(g12.6,1x))
1100  format (i3,' Successes in ',i4,' tries.')
1110  format (/,'Best result so far is try ',i3,', bstobj = ',g12.6,/)
1120  format ('ntry = ',i3,', igood = ',i3,', icount = ',i3,
     *        ', obj = ',g12.6,', ibest = ',i3,', bestobj = ',g12.6)

      end

      subroutine bstmd1 (randm,n6out, bstx)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, n, icount, jcount, ifault, j, k, igood, ibest

      logical readyn, bad, randm, n6out

      double precision var(l2+k5), objf, x(l2+k5), sx(l2+k5), x0(l2+k5),
     *                 tol, step(l2+k5), bstobj, bstx(*), rvar,
     *                 bstvar(l2+k5)

      external readyn, mcobj2

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)
c----------------------------------------------------------------------- 
c                                 read experimental data and inversion candidates
c                                 if randm, then experimental data is perturbed within
c                                 its uncertainty.
      call mcxpt (n6out,randm)
c                                 parameter max - min
      n = 0
c                                 compounds
      do i = 1, mccpd
c                                 for each coefficient
         do j = 1, mcpct(i)
            n = n + 1
            plow(n) = cprng(i,j,1)
            pdelta(n) = cprng(i,j,2) - cprng(i,j,1)
c           pcent(n) = (cprng(i,j,2) + cprng(i,j,1))/2d0
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
c              pcent(n) = (sprng(i,j,k,2) + sprng(i,j,k,1))/2d0
            end do
         end do
      end do
c                                 unscaled step size for search
      step(1:n) = frac*pdelta(1:n)
c                                 initialize scaled coordinate
      sx(1:n) = 0.5d0

      ibest = 0
      igood = 0
      bstobj = 1d99
      oprt = .false.
      fprint = .false.

      if (n6out) then 
c                                 output file
         call mertxt (tfname,prject,'.out',0)
         open (n6,file=tfname)

         write (n6,*) 'tol/frac/simplx',tol, frac, simplx
         write (n6,'(80(''-''))')
      end if 

      do i = 1, mtry

         write (*,1010) i
         write (*,1080) sx(1:n)
c                                 unscale sx
         do j = 1, n
            x(j) = plow(j) + sx(j)*pdelta(j)
         end do 

         x0 = x
c                                 initialize icount in case of failure
         icount = 0

         write (*,1085) x(1:n)

         call minim (x, step, n, objf, kcount, jprint, invtol, 
     *                  conchk, iquad, simplx, covar, mcobj2, icount, 
     *                  jcount, ifault, oktol)

         if (ifault.gt.2.or.(ifault.gt.0.and.objf.gt.1d2)) then 

            write (*,1020) ifault, icount

         else

            igood = igood + 1

            if (objf.lt.bstobj) then 
               ibest = i
               bstobj = objf
               bstx(1:n) = x(1:n)
               bstvar = var
            end if

            if (n6out) then
               write (n6,1120) i, igood, icount, objf, ibest, bstobj
               write (n6,1085) x0(1:n)
               write (n6,1030) x(1:n)
               write (n6,'(80(''-''))')
            end if

            write (*,1030) x(1:n)
            write (*,1050) icount, igood

         end if

         write (*,'(80(''-''))')
c                               new starting point
         do j = 1, n
            call random_number (sx(j))
         end do

      end do
c                               print residuals for best model
      x(1:n) = bstx(1:n)
      oprt = .true.
c                               write best model to *.bst
      write (n7,'(20(g12.6,1x))') bstx(1:n), bstobj

      fprint = .true.
      write (n6,'(80(''-''))')
      write (n6,1100) igood, mtry
      write (n6,1110) ibest, bstobj
      write (n6,'(80(''-''))')

      call mcobj2 (x,objf,rvar,bad)
      fprint = .false.

      write (*,1100) igood, mtry
      write (*,1110) ibest, bstobj

      if (n6out) close (n6)

1010  format (/,'Try: ',i4,/)
1020  format (/,'Minimization FAILED, ifault = ',i3,', icount = ',i4,/)
1030  format ('Final coordinates: ',5(g12.6,1x))
1050  format (/,'Number of function evaluations: ',i5,', igood = ',i3,/)
1080  format ('Initial normalized coordinates: ',5(g12.6,1x))
1085  format ('Initial coordinates: ',5(g12.6,1x))
1100  format (i3,' Successes in ',i4,' tries.')
1110  format (/,'Best result so far is try ',i3,', bstobj = ',g12.6,/)
1120  format ('ntry = ',i3,', igood = ',i3,', icount = ',i4,
     *        ', obj = ',g12.6,', ibest = ',i3,', bestobj = ',g12.6)

      end

      subroutine invers 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical readyn

      integer i, j, k, lu, np1, bfit

      double precision bstx(l12), x(200,l12), sx(l12), ss(l12), 
     *                 bstobj, bstcov(l12**2), oldobj

      external readyn
c----------------------------------------------------------------------- 
c                                 output summary files:
c                                 ---------------------
c                                 central model results
      call mertxt (tfname,prject,'_central.pts',0)
      open (n0,file=tfname)
c                                 perturbed model results
      call mertxt (tfname,prject,'_perturbed.pts',0)
      open (n11,file=tfname)
c                                 legible output, all results
      call mertxt (tfname,prject,'.out',0)
c                                 in case the user has set print, detach the file
      close (n6)
      open (n6,file=tfname)

      write (n6,'(80(''-''))')
c                                 initialize drand
      if (seed) call random_seed

      uncert = .false.

      if (.not.mchot) then 
c                                 get best model, 1st argument sets 
c                                 random perturbation off
c                                 arguments: init, x
         call bstmod (.false., bstx, bfit, bstobj, bstcov, 1)
c                                 set cblk to bstx
         write (n0,1020) 1, bfit, (xinv(j),j=1,nparm), bstobj, 999d0, 
     *                   0d0, 0d0, (cblk(k), k = 1, lbulk)

         x(1,1:nparm) = xinv(1:nparm)
         x(1,nparm+1) = bstobj

      else
c                                 this call is just doing initialization.
c                                 arguments: init, x
         call bstmod (.true., bstx, bfit, bstobj, bstcov, 1)

         if (nunc(1).eq.0) call errdbg ('the value of NUNC in the '//
     *                         'must be > 0 for uncertainty analysis')

         write (*,1000) 'Enter the ',nparm,' coordinates and OBJ of the'
     *                //' model you want uncertainties for:'

         read (*,*) (x(1,j),j=1,nparm+1)
c                                 set the fit-all-data flag
         write (*,1100)
        
         if (readyn()) then 
            bfit = 1
         else 
            bfit = 0
         end if

      end if
c                                 write table header to perturbed pts
      write (n11,1090) vnm(1:nparm),'OBJF','Marker','Misfit',
     *                              'Bayes',cname(1:lbulk)
c                                 write best model result to perturbed file
      write (n11,1020) 1, bfit, (x(1,j),j=1,nparm+1), 999d0, 0d0, 0d0, 
     *                 (0d0, k = 1, lbulk)
c                                 trick to save parm ranges in pts file format
      lu = n0
      do k = 1, 2

         do i = 1, 2
            write (lu,1020) 8+i, 0, 
     *                      (plow(j) + (i-1)*pdelta(j), j=1,nparm),
     *                      bstobj, 999d0, 0d0, 0d0, (0d0, j = 1, lbulk)
         end do

         lu = n11

      end do

      oldobj = bstobj
c                                 uncert turns uncertainty analysis on in general
      uncert = .true.
c                                 invprt turns thermodynamic uncertainty on specifically
      if (invunc.eq.1.or.invunc.eq.3) then 
         invprt = .true.
      else 
         invprt = .false.
      end if

      mpert = 1

      np1 = nparm + 1

      do i = 1, nunc(1)

         write (*,1040) i, nunc(1)
c                                 opnimc opens n8 and reads to the end of section 4
c                                 of the imc file, mccomp then reads section 5 and
c                                 closes (n8), for invptx this is only necessary 
c                                 for inunc < 3
         if (invunc.lt.3) call opnimc
c                                 set number of attempts
         mtry = ptry

         bstx(1:nparm) = x(1,1:nparm)
c                                 arguments: init, x
         call bstmod (.false.,bstx, bfit, bstobj, bstcov, i)

         if (nogood) then 

            lu = 6

            do j = 1, 1
               if (mtry.gt.1) then
                  write (lu,1060) i
                  if (missng) write (lu,1070)
               end if
               lu = n6
            end do

            cycle

         else 

            mpert = mpert + 1

         end if
c                                 bstmod will set xinv to the best model
c                                 results
         x(mpert,1:nparm) = xinv(1:nparm)
         x(mpert,np1) = bstobj
c                                 write perturbed model results
         if (bayes) then
            write (n11,1020) 5, jfit, (x(mpert,j),j=1,np1),
     *                            999d0, lik, bay,
     *                            (cblk(k), k = 1, lbulk)
         else
            write (n11,1020) 2, jfit, (x(mpert,j),j=1,np1),
     *                            999d0, lik, bay,
     *                            (cblk(k), k = 1, lbulk)
         end if

         sx(1:np1) = 0d0
         ss(1:np1) = 0d0

         do j = 1, mpert

            do k = 1, np1
c                                 sum x
               sx(k) = sx(k) + x(j,k)
            end do

         end do

         do j = 2, mpert

            do k = 1, np1
c                                 sum of squares
               ss(k) = ss(k) + (x(j,k) - x(1,k))**2

            end do

         end do

         if (mpert.gt.1) then

            do k = 1, np1
c                                 standard deviation
               ss(k) = dsqrt(ss(k)/(mpert-1d0))

            end do

         end if 

         lu = 6

         do j = 1, 2

            write (lu,'(80(''-''))')
            write (lu,1080) mpert
            write (lu,1025) (vnm(k), k = 1, nparm),'  OBJF  '

            write (lu,1030) 'Central model parameters:',
     *                       (x(1,k), k = 1, np1)
            write (lu,1030) 'Perturbed central model parameter: ',
     *                       (x(mpert,k), k = 1, np1)
            write (lu,1030) 'Mean perturbed parameters:',
     *                       (sx(k)/mpert, k = 1, np1)
            if (mpert.gt.1) then

               write (lu,1030) 'Parameter standard deviation:',
     *                          (ss(k), k = 1, np1)

c              write (lu,1030) 'Standard errors on the mean parameter v'
c    *                         //'alues: ',(ss(k)/(m-1d0), k = 1, nparm)

               if (nmcov.and.bstcov(1).ne.0d0) then

                  write (lu,1030) 'Nelder-Mead standard deviation:',
     *                          (dsqrt(covar(jkdiag(k))), k = 1, nparm)

               end if

            end if

            lu = n6

         end do

      end do


      if (mpert.lt.2.and.ptry.gt.0) then

            lu = 6

            do j = 1, 2
               if (ptry.ge.1.and.nunc(1).gt.0) write (lu,1065)
               lu = n6
            end do

      end if

      write (*,1050)

      close (n0)
      close (n6)
      close (n11)

      stop

1000  format (/,a,i4,a,/)
1010  format (*(1pg13.6,1x))
1020  format (1x,i2,2x,i1,2x,*(1pg13.6,1x))
1025  format (t40,*(3x,a8,3x))
1030  format (a,t40,*(1pg13.6,1x))
1040  format (2(80('-'),/),'Starting uncertainty evaluation loop ',
     *       'iteration ',i3,' of ',i4,' requested.',/)
1050  format (/,80('-'),/,
     *       'The results have been written to *.out and *.pts files. ',
     *       'The *.pts files can be',/,'processed with mc_fit_plot, ',
     *       'Perple_X_plot or PSPTS, the *.out file summarizes the',/,
     *       'central model and perturbation analyses in legible ',
     *       'format.',/)
1060  format (/,'No good solutions found for perturbation try ',i3,
     *          ' of the central model.',/,'When a lack of successful ',
     *          'perturbed models hinders uncertainty assessment ',/,
     *          'possible remedies include:',//,
     *       2x,'1 - doubling NUNC (*.imc file)',/,
     *       2x,'2 - doubling MTRY (*.imc file)')
1065  format (/,'No good solutions found for the perturbation analysis',
     *          '. When a lack of successful ',/,
     *          'perturbed models hinders uncertainty assessment ',/,
     *          'possible remedies include:',//,
     *       2x,'1 - doubling NUNC (*.imc file)',/,
     *       2x,'2 - doubling MTRY (*.imc file)')
1070  format (2x,'3 - changing no_miss to miss_ok (*.imc file)')
1080  format (/,'After ',i4,' successful perturbations:',/)
1090  format ('sym fit',*(3x,a8,3x))
1100  format (/,'Does this model fit all measured analytical data ',
     *          '(Y/N)?',
     *        /,'This is only important if you intend to use the fit-',
     *        'all-data filter in mc_fit_plot')

      end 

      subroutine bstmod (init, xcent, bfit, objf, bstcov, ipert)
c----------------------------------------------------------------------
c init   -> just initialize and return
c uncert -> use x as starting guess
c randm  -> perturb observational data within uncertainty
c xcent  -> on input, if uncert & best_central_model, best central 
c           model coordinates 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, n, icount, ifault, jcount, lu, nblen, mph,
     *        j, k, igood, ibest, jbest, ipert, ii, bfit

      logical readyn, bad, improv, init, centrl

      double precision objf, point5, bstcov(*), xcent(*), x(l12),
     *                 x0inv(l12), step(l12), ssp, rvar,
     *                 bstlx(l12), bstlik, bstlx0(l12), oxinv(l12),
     *                 bstbx0(l12), bstbx(l12), bstbay, pertrb

      character strg*35, strg1*200

      external readyn, mcobj2, nblen, pertrb

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c----------------------------------------------------------------------- 
c                                 flag to make mcobj print extended
c                                 ouput for invptx
      fprint = .false.
c                                 flag to signal whether any acceptable
c                                 models were found 
      nogood = .true.

      centrl = .not.uncert

      if (uncert) then
c                                 doing uncertainty analysis:
         write (strg,'(a,i3)') 'Central model perturbation ',ipert
c                                 all sources or just analytic
         if (invunc.eq.1.or.invunc.eq.2) call mccomp (uncert)

         if (invunc.eq.1.or.invunc.eq.3) then 
c                                 make the thermodynamic data 
c                                 perturbations
            do i = 1, imkend
c                                 i = 1:ipoint -> includes make entities,
c                                 but the increments specified here will
c                                 not be used because gcpd returns directly
c                                 for made endities; however, the makes 
c                                 invoke real data in range i = 1:mkend
c                                 if these are in the list 1:ipoint they are
c                                 made equivalent to the earlier increment.
c                                 this has the consequence that made entities
c                                 are not perturbed indendently.
               if (mkptr(i).eq.0) then
c                                 real data or a make entity that
c                                 has been excluded:
                  hinc(i) = pertrb (1d0,deltah(i))
               else 
c                                 data has already been perturbed
c                                 in the real data list
                  hinc(i) = hinc(mkptr(i))
               end if

            end do

         end if

      else
c                                 doing central model
c                                 read phase compositions
         call mccomp (uncert)
         write (strg,'(a,i3)') 'Central model '

      end if

      mcfrst = .false.
c                                 set george's "weight"
c KA hack
      do i = 1, mxpt
         xptpt(i,5) = 1d0
      end do
c                                 for phase simplex and bulk cases, the first 
c                                 ipot variables are the same
      n = ipot
c                                 phase simplex case
c!!!!!!!!!!!!!!!!!!!! need a local counter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.not.mcbulk) n = ipot + xptnph(1)

      if (unmeas.gt.0) n = n + unmeas
c KA hack
c        if (unmeas.gt.0) n = n + 2 * unmeas + 1
c                                  here's the mess:
c                                  x(1:ipot) are the independent potentials
c                                  x(ipot+1,ipot+unmeas) are the unmeasured components
c                                  and if .not.mcbulk then
c                                  x(ipot+unmeas+1,ipot+unmeas+nph) are the variables
c                                  that will be used to map the measured composition 
c                                  to a nph - 1 simplex using a scale * simplex 
c                                  transformation which requires that x(ipot+unmeas+1:
c                                  nparm) are unconstrained 
c                                  nelder-mead variable indices:
      nmi1 = ipot + 1
      nmi2 = ipot + unmeas
      nmi3 = nmi2 + 1
      nmi5 = nmi2 + xptnph(1)
      nmi4 = nmi2 + xptnph(1) - 1
      mph = xptnph(1)

      do k = 1, nmi5
         if (k.le.ipot) then
            plow(k) = vmin(k)
            pdelta(k) = vmax(k) - vmin(k)
         else
            plow(k) = 0d0
            pdelta(k) = 1d0
         end if
      end do
c                                  initialize ranges for "unconstrained" variables
c                                  guessing scales for nelder-mead variables to
c                                  linear values
      xlow(1:n) = 0d0
      xdelta(1:n) = 1d0
c                                  adjust fornelder-mead box transform variables
      if (sftpot) then
         xlow(1:ipot) = -sftdel/2d0
         xdelta(1:ipot) = sftdel
      end if 

      if (sftc) then 
         xlow(nmi1:nmi2) = -sftdel/2d0
         xdelta(nmi1:nmi2) = sftdel
      end if

      if (sftmax) then 
c                                  nelder-mead scale * softmax transform
c                                  variables
c                                  scale * softmax tranform variables:
         xlow(nmi3:nmi4) = -dlog(dfloat(mph)) - sftdel/2d0
         xlow(nmi5) = -dlog(dfloat(mph-1)) - sftdel/2d0
         xdelta(1:nmi5) = sftdel

      end if
c                                 initialize to linear step
       step(1:n) = frac
c                                 unscale step size for search
      if (sftpot) step(1:ipot) = frac*xdelta(1:ipot)

      if (sftc) step(nmi1:nmi2) = frac*xdelta(nmi1:nmi2)

      if (sftmax) step(nmi3:nmi5) = frac*xdelta(nmi3:nmi5)

      nparm = n
      ibest = 0
      igood = 0
      ifault = 0
      bstlik = 1d99
      bstbay = 0d0
      jbest = 0
      point5 = 0.5d0
c                                 save indices of diagonal elements
c                                 of the nelder-mead covariance matrix
      k = 0

      do j = 1, nparm
         jkdiag(j) = k + j
         k = k + nparm
      end do

      if (init) return

      if (centrl) then
c                                 write header to central pts
         write (n0,1030) vnm(1:nparm),'OBJF','Marker','Misfit',
     *                                'Bayes', cname(1:lbulk)

      end if
c                                 if doing uncertainty on input xcent is the 
c                                 central model x
      do i = 1, mtry
c                                 counter in mcobj2
         optct = 0
c                                 initialize covar(1) to prevent re-use of old 
c                                 covariance results
         covar(1) = 0d0
c                                 invxpt and invptx have different
c                                 initialization, and georges grid search 
c                                 yet another
         if (uncert.and.mtry.eq.1) then
c                                 using the best central model coordinates as
c                                 starting guess and perturbing the 
c                                 observational data
            xinv(1:n) = xcent(1:n)
            call trnsyx (x)

         else if (i.eq.1) then
c                                 doing the 1st try of either central or
c                                 uncertainty analysis. initialize to linear
            x(1:n) = 0.5d0
c                                 initialize x values to range centers
c                                 box transform variables:
            if (sftpot)  x(1:ipot) = 0d0
            if (sftc)  x(nmi1:nmi2) = 0d0
c                                 scale * softmax tranform variables:
            if (sftmax) then 
               x(nmi3:nmi4) = -dlog(dfloat(mph))
               x(nmi5) = -dlog(dfloat(mph-1))
            end if

         else if (uncert.or.(centrl.and.i.gt.1)) then
c                                  doing different starting guess for the 
c                                  central model or doing uncertainty with
c                                  random starting guesses.
            do j = 1, n

               if (sftpot.and.j.le.ipot .or. sftc.and.j.ge.nmi1) then
c                                 sigmoid box transformation x => [xlow,xlow+xdelta]
                  x(j) = xlow(j) + pertrb (point5,1d0) * xdelta(j)
               else
c                                 linear transformation x => [0,1]
                  x(j) = pertrb (point5,1d0)

               end if

            end do

         end if 
c                                 convert x to xinv
c                                 save starting coordinate
         call trnsxy (x)
         x0inv(1:n) = xinv(1:n)
c                                 Nelder-Meade search:
c                                 initialize icount in case of failure
         icount = 0

         call minim (x, step, n, objf, kcount, jprint, invtol, 
     *                  conchk, iquad, simplx, covar, mcobj2, icount, 
     *                  jcount, ifault, oktol)
c                                 check that this is the best x
         oxinv = xinv
         call trnsxy (x)

         do j = 1, n
            if (oxinv(j).ne.xinv(j)) then
               write (*,*) 'whoa ',j,xinv(j) - oxinv(j)
            end if
         end do

         if (ifault.gt.2.or.(ifault.gt.0.and.objf.gt.oktol)) then
c                                 write failure msgs 
            if (uncert) then

               write (strg1,'(a,a,i5,a,i5,a)') strg,
     *            ' of try ',i,
     *            ' did not converge after ',icount,
     *            ' misfit function evaluations.'

            else 

               write (strg1,'(a,a,i5,a,i5,a)') strg,', Try ',i,
     *            ' did not converge after ',icount,
     *            ' misfit function evaluations.'

            end if

            call deblnk (strg1)

            lu = 6

            do 
               write (lu,'(80(''-''))')
               write (lu,'(a)') strg1(1:nblen(strg1))
               write (lu,1080) igood
               if (uncert) write (lu,1085) mpert - 1
               
               write (lu,'(80(''-''))')
               if (lu.eq.n6) exit
               lu = n6
            end do 
c                                  cycle won't work for george
            cycle

         else

            write (strg1,'(a,a,i4,a,i4,a)') strg,' Try ',i,
     *                                        ' converged'
            call deblnk (strg1)
c                                  got a legitimate result
            if (nomiss.and.missng.or.objf.gt.mxobjf) then
c                                  but reject because either no_miss or 
c                                  objf > filter
c                                  -----------------------------------------
c                                  mustft is not properly implemented here
c                                  because it's only set if mcobj2 is called
c                                  with fprint = .true.
c              .or. mustft.and..not.fitdat) then

               lu = 6

               do

                  write (lu,'(80(''-''))')
                  write (lu,'(a,a,a)') 'Although ',strg1(1:nblen(strg1))
     *                           ,' the result will be rejected because'
                  if (objf.gt.mxobjf) then 
                     write (lu,1091)
                  else if (missng) then
                     write (strg1,'(*(a,1x))') (pname(j)
     *                                          (1:nblen(pname(j))),
     *                                          j = 1, ntot)
                     call deblnk (strg1)
                     write (lu,1090) strg1(1:nblen(strg1))
                  else 
                     write (lu,1092)
                  end if

                  write (lu,1093)

                  if (lu.eq.n6) exit
                  lu = n6

               end do 

               cycle

            else

               nogood = .false.

            end if

            if (nmcov.and.ifault.eq.2.and.centrl) then
c                                 warn about failed covariance matrix
c                                 if doing central models
               lu = 6

               do
                  write (lu,1100)
                  if (lu.eq.n6) exit
                  lu = n6
               end do

            end if
c                                 count and check if model improved
c                                 bstlik - minimum misfit
c                                 bstlx  - minimum misfit model x
c                                 bstlx0 - initial x of min misfit model
c                                 bstlco - min misfit covariance/correlation matrix
c                                 bstbay - maximum bayesian probability
c                                 bstbx  - maximum bayesian x
c                                 bstbx0 - initial x of max bayesian
c                                 bstbco - max bayesian covariance/correlation matrix
c                                 bay    - current bayesian probability
            call savbst (objf,bstlik,bstlx,ssp,x0inv,bstlx0,
     *                   bstbx0,bstbay,bstbx,n,i,ibest,jbest,igood)
c                                 write all successful results to 
c                                 *_central.pts file
            if (centrl) then
               if (bayes) then 
                  write (n0,1025) 6, jfit, (xinv(j),j=1,nparm), bay, 
     *                            999d0, objf, bay,
     *                            (cblk(k), k = 1, lbulk)
               else 
                  write (n0,1025) 3, jfit, (xinv(j),j=1,nparm), objf,
     *                            999d0, objf, bay,
     *                            (cblk(k), k = 1, lbulk)
               end if
            end if

            improv = i.eq.ibest.or.i.eq.jbest

         end if
c                                 three print cases:
c                                 bstout = 0 print all converged models
c                                 bstout = 1 print only improved models
c                                 bstout = 2 print only best model, this loop
c                                 handles only the first two cases writes
c                                 to n6 (*.out), n11 (*.pts)
         if (((improv.and.bstout.eq.1).or.bstout.eq.0) .and.centrl) then
c                                 call to prtsum may change cblk (bulk composition)
c                                 but keeps original x
            call prtsum (objf,bstlik,ssp,x0inv,
     *                   bstbay,n,ibest,jbest,igood,jcount,
     *                   icount,strg1)
c                                 write details of scoring to *.out console and n6
            fprint = .true.
            call mcobj2 (x,objf,rvar,bad)
            fprint = .false.

         end if

      end do
c                               write warnings if no good results
      if (centrl) then
c                               use rndm to flag central model
         if (nogood) then
c                               terminate if no good result on central model
            lu = 6

            do j = 1, 2
               write (lu,1060)
               if (missng) write (lu,1070)
               write (lu,'(/)')
               lu = n6
            end do

            stop

         end if

      else
c                               return if no good result on a pertubed model
         if (nogood) return

      end if

      if (bstout.eq.3.and.centrl) then
c                               only writing best model output
c                               write best model to *.out
         do j = 1, 2
c KA hack, this isn't doing anything useful
            do ii = 1, mxpt
               call mcsetb (ii)
            end do
c                                 loop to write to console and *.out
c                                 write notice and stats to console
            call prtsum (objf,bstlik,ssp,x0inv,
     *                   bstbay,n,ibest,jbest,igood,jcount,
     *                   icount,strg1)

         end do

      end if
c                                 set choice for statistics computation
c                                 on return to invers
      if (bayes) then
         xinv(1:n) = bstbx(1:n)
         bstcov(1:n**2) = bstbco(1:n**2)
         bfit = bayfit
         objf = bstbay
      else
         xinv(1:n) = bstlx(1:n)
         bstcov(1:n**2) = bstlco(1:n**2)
         bfit = likfit
         objf = bstlik
      end if

      if (bstout.eq.3.and.centrl) then
c                                 mcobj2 uses the Nelder-Mead variables
            call trnsyx (x)
c                                 write details of scoring to *.out console and n6
            fprint = .true.
            call mcobj2 (x,objf,rvar,bad)
            fprint = .false.

      end if

1010  format (i3,1x,2a,g12.6,:,3h * ,g8.3)
1020  format (80('-'),/,'Search FAILED, ifault = ',i3,', icount = ',
     *       i4,', igood = ',i4,', mtry = ',i7,/,80('-'))
1025  format (1x,i2,2x,i1,2x,*(1pg13.6,1x))
1030  format ('sym fit',*(3x,a8,3x))
1060  format (/,'No good solutions found for your unperturbed observat',
     *          'ional data execution will',/,'be terminated. Possible',
     *          ' remedies:',//,
     *       2x,'1 - double MTRY (*.imc file)')
1070  format (2x,'2 - change no_miss to miss_ok (*.imc file)')
1080  format (i4,' tries have converged so far.')
1085  format (i4,' successful perturbations so far.')
1090  format ('one or more observed phases are not predicted. To use ',
     *        'such results set the ',/,'mc_fit no_missing_phases opti',
     *        'on to F. The predicted assemblage was:',/,1x,a)
1091  format ('the misfit function exceeds the mc_fit misfit_filter_va',
     *        'lue option value.')
1092  format ('the model does not predict all analytical data within ',
     *        'its uncertainty')
1093  format (80('-'))
1100  format (/,'**warning ver505** quadratic surface fitting failed. ',
     *        'The Nelder Mead Covariance',/,'and correlation matrices',
     *        ' will not be available for this result. Increasing',/,
     *        'expansion_criterion may remedy this problem.',/)
2000  format (/,'Try ',i5,' has converged at:',/)
2010  format (/,'Misfit for this Try: ',g12.6)
2020  format ('This is the minimum Misfit obtained so far.')
2030  format ('The minimum Misfit (',g12.6,
     *        ') so far was obtained on Try',i5)
2040  format (/,'Bayesian score for this Try: ',g12.6)
2050  format ('This is the best Bayesian score obtained so far.')
2060  format ('Best Bayesian score (',g12.6,
     *        ') so far was obtained on Try',i5)
2070  format (29x,a8,' = ',g12.6)

      end

      subroutine mccomp (randm)
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for MC thermobarometry
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical randm, bad

      integer ier, nblen

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
c                                 to allow invptx to use the same objective 
c                                 function and i/o structure as invxpt the 
c                                 multi-experiment arrays of invxpt are used
c                                 here even though invptx treats only one 
c                                 observation (mxpt = 1). 
      mxpt = 1
c                                 initialize phase weighting
      invwt = 1d0 
c                                 number of phase compositions
      cxpt = 0
c                                 unmeasured component counter
c KA hack
c10    unmeas = 0
      unmeas = 0
c                                 effect bulk specification
      mcbulk = .false.

      vnm(1) = 'P_bar'
      vnm(2) = 'T_K'
c                                 -------------------------------------
c                                 IMC file input Section 5
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.ne.'begin_assemblage') call errdbg (
     *   'expecting begin_assemblage keyword, found '//key)
c                                 read sample name
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.ne.'sample_name') call errdbg (
     *   'expecting sample_name keyword, found '//key)

      read (strg,*) xptnam(mxpt)
c                                 pressure range
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'pressure_range') then

         read (strg1,*) vmin(1), vmax(1)
         if (mcfrst) write (*,'(a)') 'read pressure range'

      else

         call errdbg ('expecting pressure_range tag, found: '
     *               //key)

      end if
c                                 temperature range
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'temperature_range') then

         read (strg1,*) vmin(2), vmax(2)
         if (mcfrst) write (*,'(a)') 'read temperature range'

      else

         call errdbg ('expecting temperature_range tag, found: '
     *               //key)

      end if

      call redunm
c                                 read optional effective bulk
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'begin_bulk'.or.key.eq.'phase_name') then

         backspace (n8)

      else

         call errdbg ('expecting unmeasured_component, phase_name, or b'
     *              //'egin_bulk tag, found: '//key)

      end if
c                                 read the assemblage data
      call gtassmb (comp,ecomp,randm,'assemblage',bad)

c KA hack
c     if (mxpt.eq.1) then
c         mxpt = mxpt + 1
c        goto 10
c     end if

      close (n8)

      end

      subroutine redunm
c-----------------------------------------------------------------------
c read unmeasured components and corresponding limit expressions
c set meas/unmeas component counters
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer ier, i, j

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      double precision comp(k5), ecomp(k5)
c-----------------------------------------------------------------------
c                                 initialize simple limit flag
      simlim = .false.
      needms = .false.

      do 
c                                 read optional unmeasured components
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'redox_component'.or.key.eq.'unmeasured_component') 
     *                                                             then

            if (mcfrst) write (*,'(a)') 'reading '//key

            unmeas = unmeas + 1

            uncomp(unmeas) = 0

            do i = 1, lbulk
               if (val.eq.cname(i)) then
                  uncomp(unmeas) = i
                  exit
               end if
            end do

            if (uncomp(unmeas).eq.0) call errdbg (val//' is an invalid'/
     *                               /' component name')

            vnm(unmeas+2) = 'C_'//val

            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 look for limiting components
            if (key.eq.'begin_limits') then

               if (mcfrst) write (*,'(a)') 'reading limiting components'

               call gtcomp (comp,ecomp,.false.,'limits',bad,.false.)

               if (bad) call errdbg (
     *                              'error reading limiting components')

               do i = 1, lbulk
                  cmpmin(unmeas,i) = comp(i)
                  cmpmax(unmeas,i) = ecomp(i)
               end do

            else

               call errdbg ('expecting begin_limits tag found: '//key)

            end if

         else
c                                 no limiting component specifications
            backspace (n8)
            exit

         end if
c                                 end of unmeasured component loop
      end do
c                                 make a list of the measured thermodynamic components:
      meas = 0

      do i = 1, jbulk

         bad = .false.

         do j = 1, unmeas
            if (uncomp(j).eq.i) then 
               bad = .true.
               exit
            end if
         end do

         if (bad) cycle

         meas = meas + 1
         cmeas(meas) = i

      end do

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

         pertrb = num *( 1d0 + 2d0*err*(x - 0.5d0) )

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

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer length,com
      character chars*1, card*400
      common/ cst51 /length,com,chars(400),card

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
c                                 reading a parameters for a
c                                 compound
         if (id.lt.0) then

            call redcpd (id)

         else if (id.gt.0) then
c                                 mcsol counts parameter groups,
c                                 rather than solutions
            mcsol = mcsol + 1

            mcids(mcsol) = id
c                                initialize term counter
            mctrm(mcsol) = 0
c                                initialize coefficient counter
            mccoef(mcsol,1:m1) = 0
c
            done = .false.
c                                 read parameter group for a solution
c                                 first read the initial tag to determine
c                                 the group type
c                                 W -> excess function group
c                                 D -> internal DQF group
c                                 e -> end_list -> end of group
            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

            if (key(1:1).eq.'D'.or.key(1:1).eq.'d') then
c                                 DQF group
               call reddqf (id)
c                                 set the group type (excess term)
               mcspid(mcsol) = 2

            else if (key(1:1).eq.'W'.or.key(1:1).eq.'w') then
c                                 Excess group
               call redtrm (id)
c                                 set the group type (excess term)
               mcspid(mcsol) = 1

            else 

               call errdbg ('invalid parameter'//' for '//tname)

            end if 
c                                 flag to make gall recalculate
c                                 static compounds
            mcflag(id) = .true.

         end if

      end do
c                                 holy schmoly! if you followed that mess
c                                 you deserve a medal. count the parameters
c                                 and echo them if requested
      do i = 1, 2
         if (i.eq.1.and.n6out) then
            call parmct (n6out,n6)
         else if (i.eq.2) then
            call parmct (.true.,6)
            exit
         end if
      end do

      if (grh) then
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

            write(*,*) 'Grid search will explore ',ngrid,
     *                 ' combinations.'

         end if

      end if
c                                 read unmeasured component stuff, this
c                                 needs to be called to set component
c                                 counters.
      call redunm
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

         if (mcfrst) write (*,'(a)') 'reading data for experiment '//key
c                               read expt p,t, save any error
         if (grh) then

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

         else

            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

            read (key,*) xptpt(mxpt,1)
            read (val,*) err

            if (randm) xptpt(mxpt,1) = pertrb (xptpt(mxpt,1),err)

            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

            read (key,*) xptpt(mxpt,2)
            read (val,*) err

            if (randm) xptpt(mxpt,2) = pertrb (xptpt(mxpt,2),err)

            xptpt(mxpt,5) = 1d0

         end if
c                               read bulk composition
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.ne.'begin_bulk') call errdbg (
     *                           'expecting begin_bulk, found '//key)

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

1000  format (/,'**warning ver502** observation: ',a,/,'has been rejec',
     *          'ted because it includes a component not specified in: '
     *         ,a,//,80('-'))
1030  format ('or because it has > 1 composition with no uncertainty')
1040  format (/,'**warning ver502** parameter ',a,' skipped - ',
     *          'not a compound or solution',//,80('-'))

      end

      subroutine gtassmb (comp,ecomp,randm,tag,bad)
c-----------------------------------------------------------------------
c read phase data between between begin_//tag/end_//tag keywords. 

c      data includes:

c           optional bulk composition
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

      integer i, j, ids, nph, ier, nblen

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40, tag*(*)

      double precision comp(*), ecomp(*), repurb, tot

      external repurb, nblen

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision cp
      common/ cst12 /cp(k5,k10)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c----------------------------------------------------------------------
      nph = 0
      ok = .true.
      absent = .true.
      msolct(mxpt,1:isoct) = 0
      mcmode(mxpt) = .false.

      if (mcfrst) then
         write (*,'(a)') 'reading assemblage data'
      end if 
c                                 optional bulk composition
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

      if (key.eq.'begin_bulk') then

         mcbulk = .true.

         if (mcfrst) 
     *         write (*,'(a)') 'reading [effective] bulk composition'

         if (kiso) write (*,1030) 
c                                 get bulk composition
         call gtcomp (comp,ecomp,randm,'bulk',bad,.true.)

         if (kiso) write (*,'(/)')

         if (bad) call errdbg ('error reading bulk composition')

         do i = 1, lbulk

            xptc(cxpt+i) = comp(i)
            xpte(cxpt+i) = ecomp(i)

         end do
c                                 pointer to composition
c                                 increment composition pointer
c                                 pointer to the composition of phase nph in expt mexpt 
         blkptr(mxpt) = cxpt

         cxpt = cxpt + lbulk

      else

         if (.not.invxpt) mcbulk = .false.
         backspace (n8)

      end if

      do
c                                 now read phase data
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.ne.'phase_name') then

            if (key.eq.'end_'//tag) exit

            call errdbg ('the end_'//tag//' tag is missing or '
     *                 //'there is invalid data around '//strg1)

         end if
c                                 phase name
         read (strg,'(a)') tname

         if (mcfrst) write (*,'(a,a)') 'reading data for ',tname

         nph = nph + 1
c                                 check name
         call matchj (tname,ids)
c                                 instead of making this an error 
c                                 could just set bad to reject the xpt
         if (ids.eq.0) then 

            call errdbg ('no such entity as: '//tname)

         else if (ids.lt.0) then
c                                 check that the phase doesn't match a
c                                 composant
            do i = 1, isat
               do j = 1, isct(i)
                  if (sids(i,j).eq.-ids) then

                     write (*,1020) tname(1:nblen(tname)),
     *                              cname(icp+i)(1:nblen(cname(icp+i))),
     *                              tname(1:nblen(tname))
                     call errpau

                  end if
               end do
            end do
         end if
c                                 all clear, save id
         xptids(mxpt,nph) = ids
c                                 initialize phase weight
         invwt(mxpt,nph) = 1d0
c                                 and name:
         if (ids.lt.0) then 
            vnm(unmeas+2+nph) = tname
         else 
            vnm(unmeas+2+nph) = aname(ids)
         end if
c                                 get optional inversion weighting 
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'end_'//tag) then

            exit

         else if (key.eq.'inversion_weight') then     

            read (strg,*,iostat=ier) invwt(mxpt,nph)

            if (ier.ne.0) call errdbg ('the '//key(1:nblen(key))//
     *                        ' keyword requires a numeric value')

         else 

            backspace (n8)

         end if 
c                                 look for optional modal data, only 
c                                 invptx, so 1d arrays
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'phase_mode') then
c                                 modal data
            read (strg,*,iostat=ier) pmode(mxpt,nph)
            read (nval1,*,iostat=ier) emode(mxpt,nph)

            if (ier.ne.0) call errdbg ('the '//key(1:nblen(key))//
     *                     ' keyword requires two numeric values')

            mcmode(mxpt) = .true.

            if (pmode(mxpt,nph).le.0d0) call errdbg ('modes must be '//
     *               '> 0, reading: '//tname(1:nblen(tname)))

            if (emode(mxpt,nph).le.0d0) call errdbg ('modal uncertain'//
     *             'ties must be > 0, reading: '//tname(1:nblen(tname)))

            if (randm.and.(invunc.eq.1.or.invunc.eq.2)) then
c                                 perturb modes if doing uncertainty
c                                 analysis

            end if

         else
c                                 no modal data, initialize
            backspace (n8)
            if (mcmode(mxpt)) call errdbg ('no modal data for '//
     *                                    tname(1:nblen(tname))//
     *                        ' if modal data is provided, it '//
     *                        'must be provided for all phases.')

         end if
c                                 if compound don't read composition
         if (ids.lt.0) then

            do i = 1, lbulk
               if (cp(i,-ids).ne.0d0) absent(i) = .false.
            end do

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
c                                write phase name for molar composition output
         if (kiso) write (*,1010) tname(1:nblen(tname))

         call gtcomp (comp,ecomp,randm,'comp',bad,.true.)

         if (kiso) write (*,'(/)')

         do i = 1, lbulk
c                               check for non-zero components
            if (comp(i).ne.0d0) absent(i) = .false.
            xptc(cxpt+i) = comp(i)
            xpte(cxpt+i) = ecomp(i)

         end do
c                                 pointer to the composition of phase nph in expt mexpt 
         xptptr(mxpt,nph) = cxpt
c                                 increment composition pointer
         if (.not.bad) cxpt = cxpt + lbulk

      end do
c                                 phase counter
      xptnph(mxpt) = nph

      if (kiso) stop

      if (mcmode(mxpt)) then
c                                 modes have been provided
         tot = 0d0

         do i = 1, nph
            tot = tot + pmode(mxpt,i)
         end do

         if (tot.gt.1d0+nopt(9).or.tot.lt.1d0-nopt(9)) then

            if (mcfrst) 
     *         write (*,1040) xptnam(mxpt)(1:nblen(xptnam(mxpt))),tot

            do i = 1, nph
               pmode(mxpt,i) = pmode(mxpt,i)/tot
               if (.not.relerr) emode(mxpt,i) = emode(mxpt,i)/tot
            end do

         end if

         if (isat.gt.0) then
c                                 saturated components being used with modal
c                                 data warning:
            write (*,1050)
            call wrnstp

         end if

         if (.not.relerr) then
c                                 convert absolute modal error to relative error
            do i = 1, nph
               emode(mxpt,i) = emode(mxpt,i)/pmode(mxpt,i)
            end do

         end if

         if (randm) then

            do i = 1, nph
c                                 doing uncertainty analysis, perturb the modes
               pmode(mxpt,i) = repurb (pmode(mxpt,i),emode(mxpt,i))

            end do

         end if

      end if

      if (invxpt) then

      else 

         if (nph.lt.2) call errdbg ('input must specify > 1 phase')

         ok = .true.

         do i = 1, meas
c                                 don't bother with saturated/mobile
c                                 components.
            bad = .false.

            do j = 1, icp
c                                 don't bother with unmeasured components
               if (uncomp(j).eq.i) then
                  bad = .true.
                  exit
               end if

            end do

            if (.not.absent(i).or.bad) cycle

            write (*,1000) cname(i)
            ok = .false.

         end do

         write (*,'(/)')

         if (.not.ok) call wrnstp

      end if

1000  format (/,'**warning ver777** the problem definition file specif',
     *        'ies a component (',a,') that',/,'is absent from the obs',
     *        'erved phase assemblage specified in the *.imc file. To',/
     *       ,'avoid bad practice delete the absent component.')
1010  format (/,'Converted molar composition of ',a,
     *          ' follows, with RELATIVE error:',/)
1020  format (/,'**error ver778** the stability of ',a,' is determined',
     *          ' by the saturated component',/,'constraint on ',a,' ei'
     *         ,'ther eliminate the constraint from the *.dat file or',
     *        /,'remove ',a,' from the list of observed phases in the ',
     *          '*.imc file.',/)
1030  format (/,'Converted bulk molar composition of follows, with ',
     *          'RELATIVE error:',/)
1040  format (/,'**warning ver779** the phase modes for ',a,' sum to ',
     *        f10.6,' they will be renormalized.',/)
1050  format (/,'**warning ver780** using modal data when one or more ',
     *          'components is saturated is bad practice.')

      end 

      double precision function repurb (x,xerr)
c-----------------------------------------------------------------------
c hack to prevent negative perturbed values generated by pertrb, a better
c alternative would be to assume a log-normal distribution.
c-----------------------------------------------------------------------
      implicit none

      double precision min, chk, x, xerr, pertrb

       external pertrb
c-----------------------------------------------------------------------
      chk = pertrb (x,xerr)
c                                 the largest negative perturbation
      min = x * (1d0 - xerr)

      if (min.lt.0d0.and.chk.lt.x) then
c                                 the perturbed value chk is on the low
c                                 side of comp(i), rescale so that it's
c                                 between 0 and comp(i)
         repurb = x * (x - chk)/(x - min)

      else

         repurb = chk

      end if

      end 

      subroutine gtcomp (comp,ecomp,randm,tag,bad,norm)
c-----------------------------------------------------------------------
c read compositional data and uncertainities between begin_\\tag/end_\\tag 
c keywords. converts mass units to molar units if lmass.
c if randm perturbs composition within uncertainty.
c returns bad if the composition includes a component not specified in
c the problem definition file. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, bad, randm, norm

      integer i, j, ier

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40, tag*(*)

      double precision tot, comp(*), ecomp(*), repurb

      external repurb
c----------------------------------------------------------------------
      comp(1:lbulk) = 0d0
      ecomp(1:lbulk) = 1d0
      bad = .false.
c                                 limits tag signals that the component 
c                                 names are followed by two molar stoichiometric
c                                 coefficients, tag is also used to shut off
c                                 mass/mol conversions
      if (tag.eq.'limits') ecomp(1:lbulk) = 0d0

      do

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.'end_'//tag) exit

         ok = .false.

         if (tag.eq.'limits'.and.(key.eq.'simple_mole'.or.
     *                            key.eq.'Simple'.or.
     *                            key.eq.'simple_mass')) then
c                                 simple limit expressions in terms of mole
c                                 or mass fractions
             simlim(unmeas) = .true.

             if (key.eq.'simple_mass') then
c                                 simple limit expressions in terms of mass
c                                 fractions
                simmol(unmeas) = .false.   
c                                 flag indicating mass units needed.
                needms = .true.   
     
             else 

                simmol(unmeas) = .true.

             end if

             ok = .true.
             i = uncomp(unmeas)

         else 
c                                 should test for component name replication?
            do i = 1, lbulk
               if (key.eq.cname(i)) then
                  ok = .true.
                  exit
               end if
            end do

         end if

         if (.not.ok) then

            if (invxpt) then 
c                                 set bad to reject the experiment but continue
c                                 reading data in case there are valid experiments
               bad = .true.
c                                 that was a bad idea, issue an error:
               call errdbg ('invalid or misspelled component: '//key)

            else if (tag.eq.'limits') then

                call errdbg ('limit expression includes a component: '//
     *          key//' not specified in the problem definition file')

            else if (tag.eq.'bulk') then

                call errdbg ('bulk composition includes a component: '//
     *          key//' not specified in the problem definition file')

            else

                call errdbg ('phase contains a component: '//key//' not'
     *          //' specified in the problem definition file')

            end if

         end if

         read (strg1,*) comp(i), ecomp(i)

         if (tag.ne.'limits') then
c                                 ecomp > 0 check to not overwrite nerone estimation flag
            if (.not.relerr.and.comp(i).ne.0d0.and.ecomp(i).gt.0d0) then
c                                 convert absolute error to relative
               ecomp(i) = ecomp(i)/comp(i)

            else if (.not.relerr.and.ecomp(i).ne.-1d0) then

               call errdbg ('absolute errors cannot be used on zeroed'//
     *                      ' compositions, switch to relative error')
            end if
c                                 convert mass input to molar units
            if (lmass) comp(i) = comp(i)/atwt(i)
c                                output molar composition and errors
            if (lmass.and.kiso) 
     *                 write (*,1000) cname(i), comp(i), ecomp(i)

         end if

      end do
c                                 norm is set by the calling program and is only false
c                                 for limit expressions.

c                                 the composition is now molar, but not normalized
      if (norm.and..not.grh) then
c                                 this tot is just for the error estimates
         tot = 0d0

         do i = 1, lbulk
            tot = tot + comp(i)
         end do

         do i = 1, lbulk
c                                 cycle on absent or unmeasured components
            if (comp(i).eq.0d0) cycle 
c                                 use error estimate functions if ecomp < 0
c                                 numbers chosen to give 2% for X = 0.5,
c                                 20% for X = 0.05 (or thereabouts). 
            if (ecomp(i).lt.0d0) ecomp(i) = 0.025d0/(comp(i)/tot)**0.7d0
c                                 perturb data
            if (randm) comp(i) = repurb (comp(i), ecomp(i))

         end do

         if (randm) then 
c                                 get new total for normalization
            tot = 0d0

            do i = 1, lbulk
               tot = tot + comp(i)
            end do

         end if

         do i = 1, lbulk
c                                 normalize composition for scoring:
            comp(i) = comp(i)/tot

         end do

      else if (norm.and.grh) then 
c                                 george's 0 counter
         j = 0

         do i = 1, lbulk
            if (ecomp(i).eq.0d0) j = j + 1
            if (randm) comp(i) = repurb (comp(i), ecomp(i))
         end do
c                                 george assumes normalized input, but has 
c                                 this suspect test. a phase with 
c                                 fewer components than the system will
c                                 fail, and what's with the > 1?
         if (j.gt.1 .and. j.ne.lbulk) bad = .true.
         xpterr(mxpt) = j

      end if

1000  format (a5,2x,g12.6,2x,g12.6)

      end

      subroutine redtrm (id)
c----------------------------------------------------------------------
c redtrm - read excess function term of the form 

c        W(name name ...) or Wk(name name ...)

c        jd - endmember pointer
c        id - solution model pointer
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, bad, done

      integer i, j, id

      integer ibeg, ier, iscan, imax, iscnlt, jd, kend

      character name*8

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      external iscnlt, iscan

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer length,com
      character chars*1, card*400
      common/ cst51 /length,com,chars(400),card

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      double precision err
c----------------------------------------------------------------------
      ier = 0
      bad = .false.
      done = .false.
c                                 find expansion type
      if (chars(2).eq.'k'.or.chars(2).eq.'K') then
         xtyp = 1
      else
         xtyp = 0
      end if

      if (xtyp.ne.extyp(id)) call errdbg ('excess term type specified'//
     *                       'for '//tname//'in the *.imc file is inc'//
     *                       'onsisetent with the expansion for '//tname
     *                        //'in the solution model file')

c                                 data found, loop to read free terms:
      do

c                                 find brackets
         ibeg = iscan (1,com,'(') + 1
         imax = iscan (1,com,')') - 1

         if (ibeg.gt.com.or.imax.gt.com) call errdbg ('invalid excess'//
     *                                       ' term format for '//tname)

         iord = 0

         do while (ibeg.lt.imax)

            call readnm (ibeg,kend,imax,ier,name)
            if (ier.ne.0) call errdbg ('invalid excess'//
     *                                 ' term format for '//tname)

            iord = iord + 1

            if (iord.gt.m2) call error (49,wg(1,1),m2,name)

            call matchj (name,jd)

            if (jd.eq.0) call errdbg ('invalid endmember '//name//
     *                ' in excess term specified for '//tname//
     *                ' in the *.imc file')

            isub(1,iord) = -jd

         end do
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
               if (xtyp.eq.1) then
c                                 redlich-kister
                  ok = ok .and. jend(id,2+i).eq.isub(1,i)
               else
c                                 margules, try match on jth term
                  ok = ok .and.isub(1,i).eq.jend(id,2+jsub(i,j,id))
               end if

            end do

            if (ok) exit

         end do
              
         if (.not.ok) then
c                                 term does not exist in solution model
            print*,(chars(i),i=1,com),'term does not exist'
            call errdbg ('quitting in redtrm')

         end if
c                                 the term exists:
c                                 increment term counter
         mctrm(mcsol) = mctrm(mcsol) + 1
c                                 save the term pointer
         mcj(mcsol,mctrm(mcsol)) = j
c                                 find the free coefficients
         do 
c                                 loop to read free coefficients in term
            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 key may be either "parameter"
c                                 or the name of a new term "W(..." or
c                                 "Wk(..." or "end_list"
            if (key.eq.'end_list') then

               done = .true.
               exit

            else if (chars(1).eq.'W'.or.chars(1).eq.'w') then
c                                 starting a new term
               exit

            else if (key.ne.'parameter') then 

               call errdbg ('invalid coefficient or parameter '//
     *                                     'specification for '//tname)

            end if
c                                 count the coefficient
            mccoef(mcsol,mctrm(mcsol)) = 
     *                   mccoef(mcsol,mctrm(mcsol)) + 1

            if (.not.grh) then
c                                 read the parameter index
               read (strg1,*) 
     *           mccoid(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)))

               call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 read the parameter range
               read (strg1,*) 
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),1),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),2)

            else
c                                 read the coefficient range and grid numbers
               read (strg1,*) 
     *           mccoid(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol))),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),1),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),2)

               call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

               read (strg1,*,iostat=ier) i,err,err,
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
               if (ier.ne.0) 
     *            sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
     *                                                             = 2d0
               if (
     *            sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
     *               .le.0d0)
     *            sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),3)
     *               = 2d0

             end if

          end do
c                                 start a new term or exit
          if (done) exit

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
c KA hack
      subroutine mcsetb (ixpt)
c-----------------------------------------------------------------------
c set bulk composition for MC thermobarometry, assumes normalized 
c compositions
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, ixpt

      double precision cmax, cmin, mtotal

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
      if (mcbulk) then
c                                 use effective bulk
c                                 composition pointer
c KA hack
         cxpt = blkptr(ixpt)

         do j = 1, lbulk
            cblk(j) = xptc(cxpt+j)
         end do

      else

         cblk = 0d0
c                                 compute bulk from phase simplex
c KA hack
         do i = 1, xptnph(ixpt)
c                                 composition pointer
c KA hack
            cxpt = xptptr(ixpt,i)

            do j = 1, lbulk

               cblk(j) = cblk(j) + xinv(nmi2+i) * xptc(cxpt+j)

               if (cblk(j).lt.0d0) then 
                  write (*,*) 'lt 0 ',j,xinv(i),cblk(j)
               end if

            end do

         end do

      end if
c                                 get the total number of moles
c                                 and total mass
      ctotal = 0d0

      do j = 1, lbulk
         ctotal = ctotal + cblk(j)
      end do

      if (unmeas.gt.0) then
c                                unmeasured components
         if (needms) then 
c                                need mass fractions
            mtotal = 0d0

            do j = 1, lbulk
               mtotal = mtotal + cblk(j) * atwt(j)
            end do

         end if

         do l = 1, unmeas
c                                 
            cmax = 0d0
            cmin = 0d0
            k = uncomp(l)

            if (simlim(l).and.simmol(l)) then 
c                                simple molar limit expression
               cblk(k) = cmpmin(l,k) + xinv(ipot+l) * 
     *                  (cmpmax(l,k) - cmpmin(l,k)) / ctotal

            else if (simlim(l).and..not.simmol(l)) then 
c                                simple mass limit expression
               cblk(k) = ( cmpmin(l,k) + xinv(ipot+l) * 
     *                    (cmpmax(l,k) - cmpmin(l,k)) * mtotal) /
     *                    atwt(k)

            else 
c                                coupled limit expression
               do j = 1, lbulk
                  if (j.eq.k) cycle
                  cmin = cmin + cblk(j)*cmpmin(l,j)
                  cmax = cmax + cblk(j)*cmpmax(l,j)
               end do
c KA hack
c              cblk(k) = cmin + (cmax - cmin) * x(ipot+l+(ixpt-1))

               cblk(k) = cmin + (cmax - cmin) * xinv(ipot+l)

            end if

         end do
c                                 update ctotal 
         do l = 1, unmeas
            ctotal = ctotal + cblk(uncomp(l))
         end do

      end if
c                                 normalize, should be to icp?
      b(1:lbulk) = cblk(1:lbulk) / ctotal

      end



      subroutine trnsyx (x)
c-----------------------------------------------------------------------
c convert the constrained xinv variables to the unconstrained nelder-mead 
c x variables.xinv are in cst68.

c x(1:nmi2) - box transformation to simple bounds
c x(nmi3:nmi5) - scale x softmax transformation to simplicial 

c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision x(*)
c-----------------------------------------------------------------------

      do i = 1, nmi2

         if (sftpot.and.i.le.ipot .or. sftc.and.i.ge.nmi1) then
c                                 undo sigmoid box transformation

            x(i) = -dlog( - (xinv(i) - plow(i) - pdelta(i)) /
     *                      (xinv(i) - plow(i)) )
         else
c                                 undo linear transformation
            x(i) = (xinv(i) - plow(i))/pdelta(i)
         end if

      end do
c                                 done if using bulk composition
      if (mcbulk) return
c                                 undo softmax:
      x(nmi5) = -dlog((1-xinv(nmi5))/xinv(nmi5))

      do i = nmi3, nmi4

         if (xinv(i).gt.nopt(50)) then
            x(i) = dlog(xinv(i))
         else 
            x(i) = dlog(nopt(50))
         end if

      end do
c                                 end trnsyx
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

      do i = 1, lbulk
         cblk(i) = xptblk(id,i)
      end do
c                                 get total moles to compute mole fractions 
      do i = 1, lbulk
         ctotal = ctotal + cblk(i)
      end do

      do i = 1, lbulk
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


      double precision function score (kd,id,j,fit)
c-----------------------------------------------------------------------
c a function to evaluate the distance between oberved and candidate 
c compositions for MC, parameters are:
c     kd - phase indicator
c     id - experiment number
c     j  - phase number in experiment
c     fit - true if model fits within error
c-----------------------------------------------------------------------
      include 'perplex_parameters.h'

      integer id, kd, j, k, l, m, ll

      logical fit

      double precision total, term
c-----------------------------------------------------------------------
      score = 0d0
      total = 0d0
      m = 0
c                                 normalization depends on whether any
c                                 compositional uncertainty given
      do ll = 1, meas

         l = cmeas(ll)

         total = total + pcomp(l,kd)
         if (xpte(xptptr(id,j)+l).eq.0d0) m = l

      end do

      if (m.eq.0) then
c                                 calculate residual
         do ll = 1, meas

            l = cmeas(ll)

            k = xptptr(id,j)+l

            if (pcomp(l,kd).eq.0d0.and.xptc(k).eq.0d0) cycle

            term = dabs(pcomp(l,kd)/total - xptc(k))

            mresid = mresid + 1
            rss = rss + term**2

            if (term.gt.un2ft*xpte(k)*xptc(k)) fit = .false.

            if (.not.GRH) then

               if (lsqchi.eq.1) then 
c                                 LSQ compositional obj
                  term = term**2
 
               else if (lsqchi.eq.2) then 
c                                 Chi-square compositional obj
                  if (pcomp(l,kd).ne.0d0) then
                     term = term**2 / (pcomp(l,kd)/total)
                  else
c                                 pseudo Chi-square
                     term = term**2 / xptc(k)
                  end if 

               else if (lsqchi.eq.3) then 
c                                 weighted Chi-square compositional obj
c                                 this case arises if someone doesn't 
c                                 put a value for a component that exists
c                                 in the solution model
                  if (xpte(k)*xptc(k).eq.0d0) cycle
c
                  term = (term / (xpte(k)*xptc(k)))**2

               end if

            else if (term.lt.1d0) then 
c                                 grhobj, i'm not sure george ever had
c                                 anything but this              
               term = term**2

            end if

            score = score + invwt(id,j) * term

         end do

      else if (m.eq.1) then
c                                 GRH stuff:
c                                 only one has no uncertainty - normalize
c                                 to it, calculate weighted residual
         total = xptc(m) / pcomp(m,kd)

         do ll = 1, meas

            l = cmeas(ll)

            if (l.eq.m) cycle

            k = xptptr(id,j)+l

            term = dabs((pcomp(l,kd)*total - xptc(k))) 

            if (.not.grhobj.or.term.lt.1d0) term = term**2

            score = score + term

         end do

      end if

      end

      double precision function mscore (kd,id,j,fit)
c-----------------------------------------------------------------------
c a function to evaluate the distance between oberved and candidate 
c modes for MC, parameters are:
c     kd - predicted phase pointer
c     id - experiment number
c     j  - phase number in experiment
c     fit - true if model fits within error
c-----------------------------------------------------------------------
      include 'perplex_parameters.h'

      integer id, kd, j

      logical fit

      double precision term, mode
c-----------------------------------------------------------------------
c                                 predicted volumetric mode:
      mode = props(1,kd)*props(16,kd)/psys(1)
c                                 calculate residual
      term = dabs(mode - pmode(id,j))

      mresid = mresid + 1
      rss = rss + term**2

      if (term.gt.un2ft*emode(id,j)*mode) fit = .false.

      if (lsqchm.eq.1) then 
c                                 LSQ modal obj
         term = term**2
 
      else if (lsqchm.eq.2) then 
c                                 Chi-square modal obj
         term = term**2 / mode

      else if (lsqchm.eq.3) then 
c                                 weighted Chi-square modal obj
c                                 this case arises if someone doesn't 
c                                 put a value for a component that exists
c                                 in the solution model
         term = (term / (emode(id,j)*pmode(id,j)))**2

      end if

      mscore =  invwt(id,j) * term

      end

      subroutine x2ther (x)
c-----------------------------------------------------------------------
c a to map the MC search coordinates into thermodynamic arrays
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, n, id, jd, ids, jt

      double precision x(*)

      character name*14
c-----------------------------------------------------------------------
      n = 0
c                                 compounds
      do i = 1, mccpd
c                                 for each coefficient
         do j = 1, mcpct(i)

            n = n + 1
            id = mcid(i)
            jd = make(id)
            ids = mcpid(i,j)

            if (ids.le.3) then
               mdqf(jd,ids) = x(n)
               if (ids.eq.1.and.qdef(jd)) mdqf(jd,1) = x(1) + 1e6
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
c                                 an excess function
         if (mcspid(i).eq.1) then
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
                     wgl(mccoid(i,j,k),jt,ids) = x(n)

                  else

                     call getnam (name,ids)
                     call errdbg ('unimplemented solution type '//
     *                             'for '//name)

                  end if
c                                end coeffiecient loop
               end do
c                                end excess term loop
            end do

         else if (mcspid(i).eq.2) then
c                                 for each internal DQF
            do j = 1, mctrm(i)
c                                 DQF id
               jt = mcj(i,j)
c                                 for each coefficient
               do k = 1, mccoef(i,j)

                  n = n + 1

                  dqfg(mccoid(i,j,k),jt,ids) = x(n)
c                                 if _q make, decrement the enthalpic
c                                 term, on input this should checked,
c                                 currently it's just assumed internal
c                                 dqf's are _q makes
                  if (mccoid(i,j,k).eq.1) then

                     dqfg(mccoid(i,j,k),jt,ids) =
     *               dqfg(mccoid(i,j,k),jt,ids) - 1d6

                  end if

               end do

            end do

         else

            call errdbg ('something agly')

         end if
c                                end solution parameter group loop (mcsol)
      end do

      end

      subroutine mcobj2 (x, obj, rvar, bad)
c-----------------------------------------------------------------------
c a subprogram to evaluate objective function for MC data inversion
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer id, i, neg

      double precision x(*), obj, xptscr, value, rvar

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
c                                 untransform the nelder-mead variables
         call trnsxy (x)
c                                 set the potentials
         do i = 1, ipot
            v(i) = xinv(i)
         end do

c KA hack, moves mcsetb below
c        pdqf = x(nparm) * 8d4

      end if

      obj = 0d0
c                                 loop through observations
      do id = 1, mxpt
c                                 set p, t, bulk
         if (invxpt) then

            call mcstb2 (id,0,0,neg)

         else 
c                                 set the bulk composition
            call mcsetb (id)

         end if
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
         value = xptscr (id)

         if (isnan(value)) then
            write (*,*) 'oink 33'
         end if

         if (grh) then
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

         end if

         scores(id) = value

         obj = obj + value*xptpt(id,5)

      end do
c                                 normalize misfit to per xpt
      if (.not.grh) obj = obj / float(mxpt)

      write (*,'(a,1pg13.6,1x,i6)') 'obj, optct ',obj,optct

c                                 convert rss to half the residual variance
      if (lsqchi.eq.3) then 
         rvar = 1d0
      else 
         rvar = rss / float(mresid-nparm) / 2d0
      end if

      end

      double precision function xptscr (id)
c-----------------------------------------------------------------------
c function to score a calculation's agreement with an experiment for MC
c data inversion
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, imout(k5), imin(k5), used(k5), jmiss(k5), 
     *        ftbst, outnph(k5), jused(k5), skip, fit, tfit

      character bdname(k5)*14

      integer id, jd, ids, i, j, k, kct(k5), ksol(k5,k5), ibest, mpred,
     *        idextr(k5), jmin(k5), mextra, jdbest, kd, lu, kmin(k5), 
     *        l, nmiss, jj, jbest, jds(h5), istsat, nblen, ld, kdk, kk

      double precision score, best, res, mode, tot, extra, comp, miss, 
     *                 cscor(k5), sres(k5), sres1(k5), scmode, 
     *                 mscor(k5), mscore, sres2(k5), den

      external score, nblen, mscore

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus
c---------------------------------------------------------------------
c                                 compute the observation objective function
c                                 and the residual sum-of-squares
      kct = 0
c                                 rss
      rss = 0d0
c                                 number of residuals
      mresid = 0

      imout = .true.
      imin = .false.
      used = .false.
      jmiss = .true.
      outnph = .false.
      fit = .true.
      ksol = 0
      jmin = 0
      optct = optct + 1

      do i = 1, xptnph(id)

         do j = 1, ntot

            if (kkp(j).eq.xptids(id,i)) then
c                                 kct(j) is the number of times a predicted
c                                 solution model occurs in the observed assemblage
               kct(i) = kct(i) + 1
c                                 ksol(j,...) points to the positions of the
c                                 the predicted solution model in the observed assemblage
               ksol(i,kct(i)) = j

            end if

         end do

      end do
c                                 classify cases
      do j = 1, xptnph(id)

         if (kct(j).eq.0) cycle

         ids = kkp(ksol(j,1))

         if (ids.lt.0) then 

            outnph(j) = .true.

         else if (kct(j).ge.msolct(id,ids)) then
c                                 the predicted sollution occurs .le. times than
c                                 in the observed assemblage, a doudble loop with 
c                                 outer index xptnph will function for matching
            outnph(j) = .true.

         end if

      end do

      extra = 0d0
      comp = 0d0
      miss = 0d0
      scmode = 0d0
      mpred = 0
      mextra = 0
c                                 potential target phases:
      do j = 1, xptnph(id)

         if (kct(j).eq.0.or..not.outnph(j)) cycle

         jd = ksol(j,1)
         ids = kkp(jd)

         if (ids.lt.0) then
c                                 a compound, cycle and count in the next 
c                                 loop
            cycle 

         else if (kct(j).eq.1.and.msolct(id,ids).eq.1) then
c                                 found solution and no ambiguity
            mpred = mpred + 1
c                                 flag that the phase is not missing in the observed list
            jmiss(j) = .false.
c                                 observed list pointer
            jmin(mpred) = j
c                                 predicted list pointer
            kmin(mpred) = jd
c                                 compute composition score
            cscor(mpred) = score (jd,id,j,fit)
c                                 compute and add residual
            comp = comp + wcomp * cscor(mpred)

            if (mcmode(id)) then 
c                                 scoring on modes as well
               mscor(mpred) = mscore (jd,id,j,fit)
c
               scmode = scmode + wmode * mscor(mpred)

            end if
c                                 flag that the phase has been used in the predicted list
            used(jd) = .true.

         else 

            best = 1d99
            ok = .false.

            do i = 1, kct(j)
c                                 jd locates the candidates in the predicted
c                                 assemblage. this loop finds the best match
c                                 between observed phase j and predicted phase
c                                 jd
               jd = ksol(j,i)

               if (used(jd)) cycle

               tfit = .true.

               res = score (jd,id,j,tfit)

               if (res.lt.best) then
                  best = res
                  jdbest = jd
                  ftbst = tfit
                  ok = .true.
               end if
c                                 end phase loop over kct(j)
            end do

            if (ok) then
c                                 the best available match for observed phase j, 
c                                 i.e., predicted phase ibest, has been found.
               mpred = mpred + 1
c                                 used flags that a phase has been matced in the 
c                                 
               used(jdbest) = .true.
c                                 jmiss flags that the observed phase has been matched
               jmiss(j) = .false.
c                                 jmin points to the matched phase in the observed
c                                 assemblage
               jmin(mpred) = j
c                                 kmin points to the matching phase in the predicted
c                                 assemblage
               kmin(mpred) = jdbest

               if (fit) fit = ftbst

               cscor(mpred) = best

               comp = comp + wcomp * best

               if (mcmode(id)) then 
c                                 scoring on modes as well
                  mscor(mpred) = mscore (jdbest,id,j,fit)
c
                  scmode = scmode + wmode * mscor(mpred)

               end if

            else
c                                  the number of observed versions of the solution is
c                                  less than the number of predicted versions
               call errdbg ('wroink wroink')

            end if

         end if

      end do

      jused = .false.

      do k = 1, xptnph(id)
c                                 count compounds in this loop to maintain cpd last ordering for output
         if (kct(k).eq.0) cycle

         ids = kkp(ksol(k,1))

         if (kkp(ksol(k,1)).lt.0) then

            jd = ksol(k,1)
c                                 a compound, just count
            mpred = mpred + 1
            jmiss(k) = .false.
c                                 jmin(i) points to the phase entered 
c                                 in the observed phase list
            jmin(mpred) = k
c                                 kmin(i) points to the phase in the
c                                 predicted assemblage list
            kmin(mpred) = jd

            if (mcmode(id)) then 
c                                 scoring on modes as well
               mscor(mpred) = mscore (jd,id,k,fit)
c
               scmode = scmode + wmode * mscor(mpred)

            end if

            used(jd) = .true.

            cycle

         end if

         if (outnph(k)) cycle
c                                 if here then there are more occurences of the 
c                                 observed solution than the prediced, so need
c                                 an outer over the predicted phases
         do i = 1, kct(k)

            jd = ksol(k,i)
            if (used(jd)) cycle
            ids = kkp(jd)

            best = 1d99
            ok = .false.

            do jj = 1, msolct(id,ids)
c                                 j points to the location of the observed solution
               j = msloc(1,jj)

               if (jused(j)) cycle

               tfit = .true.

               res = score (jd,id,j,tfit)

               if (res.lt.best) then
                  best = res
                  ibest = i
                  jbest = j
                  jdbest = jd
                  ftbst = tfit
                  ok = .true.
               end if
c                                 end of inner solution loop
            end do

            if (ok) then
c                                 the best available match for observed phase j, 
c                                 i.e., predicted phase ibest, has been found.
               mpred = mpred + 1
c                                 used prevents subsequent use of the predicted
c                                 phase for matching observed phases
               used(jdbest) = .true.

               jused(jbest) = .true.
c                                 imin flags that the predicted phase is in, redundant
c                                 with used?
               imin(jdbest) = .true.
c                                 jmiss flags that the observed phase has been matched
               jmiss(j) = .false.
c                                 jmin points to the matched phase in the observed
c                                 assemblage
               jmin(mpred) = j
c                                 kmin points to the matching phase in the predicted
c                                 assemblage
               kmin(mpred) = jdbest

               if (fit) fit = ftbst

               cscor(mpred) = best

               comp = comp + wcomp * best

            end if
c                                  end of outer solution loop.
         end do

      end do

      istsat = 0
c                                  score extraneous predicted phases:
      do i = 1, ntot

         skip = .false.
c                                  skip saturated component phases:
         do j = 1, isat
            if (skp(j).eq.kkp(i)) then
               istsat = istsat + 1
               jds(istsat) = i
               skip = .true.
               exit
            end if
         end do

         if (skip) cycle
c                                  penalty is wextra * mass_fraction^2
         if (.not.used(i)) then

            mextra = mextra + 1
            idextr(mextra) = i
            extra = extra + wextra * 
     *                             (props(17,i)*props(16,i)/psys(17))**2

         end if

      end do
c                                 missing phase residual
      miss = wmiss * (1d0 - dble(mpred)/xptnph(id))
c                                 accumulate scores
      if (mpred.gt.0) then 
         xptscr = miss + extra + (comp + scmode)/mpred
      else
         xptscr = miss + extra
      end if

      if (nomiss.and.miss.gt.0d0) then

         missng = .true.
         fit = .false.

      else

         missng = .false.

      end if
c                                setting jfit here risks mistating whether the final
c                                result out of minim is correctly identified as 
c                                fitting or not fitting the data.... needs to be fixed.
      if (fit) then 
         jfit = 1
       else
         jfit = 0
      end if

      if ((.not.invxpt.and.fprint).or.invxpt) then

         lu = 6

         do l = 1, 2

            if (l.eq.2) lu = n6
c                                 if no_miss is set reject incomplete predictions
            if (missng)  then 

               write (lu,1130)

               cycle

            end if
c                                 set the common version for bstmod
            fitdat = fit
c                                 output optimal P-T and compositions
c                                 for inverse thermo-barometry
c                                 ------------------------------------
            if (.not.mustft) then

               if (fit) then

                  write (lu,1150)

               else

                  write (lu,1151)

               end if

            end if

            write (lu,1100) comp, extra, miss
            if (mcmode(id)) write (lu,1105) scmode

            if (mpred.gt.0) then 
c                                 phases observed and predicted
               write (lu,1000)

               do j = 1, mpred
c                                 jd/kmin points to the phase in the predicted assemblage
                  jd = kmin(j)
c                                 ld/jmin points to the phase in the observed assemblage
                  ld = jmin(j)

                  kd = xptptr(id,ld)

                  mode = props(1,jd)*props(16,jd)/psys(1)*1d2

                  write (lu,1021) pname(jd), 
     *                             (cname(cmeas(k)), k = 1, meas)

                  if (jd.le.np) then 
c                                 a solution
                     tot = 0d0

                     do k = 1, meas
                        tot = tot + pcomp(cmeas(k),jd)
                     end do

                    sres = 0d0
                    sres1 = 0d0
                    sres2 = 0d0

                    do kk = 1, meas

                       k = cmeas(kk)
c                                 skip unused components
                       if (xptc(kd+k).eq.0d0) cycle

                       kdk = kd + k
                       sres(kk) = (xptc(kdk) - pcomp(k,jd)/tot) 
c                                 sqrt chi
                       if (pcomp(k,jd).ne.0d0) then
                          den = pcomp(k,jd)/tot
                       else
c                                 pseudo Chi-square
                          den = xptc(kdk)
                       end if
 
                       if (den.ne.0d0) sres1(kk) = sres(kk)/dsqrt(den)
c                                standardized residual
                       sres2(kk) = sres(kk)/(xpte(kdk)*xptc(kdk))

                     end do

                     write (lu,1030) 'predicted', mode, 
     *                           (pcomp(cmeas(k),jd)/tot, k = 1, meas)

                     if (mcmode(id)) then

                        den = pmode(id,ld)*1d2 - mode

                        write (lu,1030) 'observed', pmode(id,ld)*1d2,
     *                               (xptc(kd+cmeas(k)), k = 1, meas)
                        if (runrm) write (lu,1030) 'residual', den,
     *                                          (sres(k), k = 1, meas)
                        if (rnorm) write (lu,1037) 'std_resid', 
     *                                          den/(emode(id,ld)*1d2),
     *                                          (sres2(k), k = 1, meas)
                        if (rlsq.or.(lsqchm.eq.1.and.rcobj)) 
     *                     write (lu,1037) 'LSQ', den**2,
     *                                         (sres(k)**2, k = 1, meas)
                        if (rchi.or.(lsqchm.eq.2.and.rcobj))
     *                     write (lu,1037) 'Chi', den**2/mode, 
     *                                        (sres1(k)**2, k = 1, meas)
                        if (rwch.or.(lsqchm.eq.3.and.rcobj))  
     *                     write (lu,1037) 'wChi', 
     *                                      (den/(emode(id,ld)*1d2))**2,
     *                                      (sres2(k)**2, k = 1, meas)

                        write (lu,1025) 'Mode Misfit:', mscor(j)

                     else

                        write (lu,1035) 'observed',
     *                                 (xptc(kd+cmeas(k)), k = 1, meas)
                        if (runrm) write (lu,1035) 'residual', 
     *                                           (sres(k), k = 1, meas)
                        if (rnorm) write (lu,1036) 'std_resid.', 
     *                                          (sres2(k), k = 1, meas)
                        if (rlsq.or.(lsqchi.eq.1.and.rcobj)) 
     *                     write (lu,1036) 'LSQ', 
     *                                        (sres(k)**2, k = 1, meas)
                        if (rchi.or.(lsqchi.eq.2.and.rcobj))
     *                     write (lu,1036) 'Chi', 
     *                                        (sres1(k)**2, k = 1, meas)
                        if (rwch.or.(lsqchi.eq.3.and.rcobj))  
     *                     write (lu,1036) 'wChi', 
     *                                        (sres2(k)**2, k = 1, meas)

                     end if

                     write (lu,1025) 'Composition Misfit:', cscor(j)

                  else
c                                 a compound
                     write (lu,1030) 'predicted', mode

                     if (mcmode(id)) then

                        den = pmode(id,ld)*1d2 - mode

                        write (lu,1030) 'observed', pmode(id,ld)*1d2
                        if (runrm) write (lu,1030) 'residual', den
                        if (rnorm) write (lu,1031) 'std_resid.', 
     *                                          den/(emode(id,ld)*1d2)
                        if (rlsq.or.(lsqchm.eq.1.and.rcobj)) 
     *                     write (lu,1031) 'LSQ', den**2
                        if (rchi.or.(lsqchm.eq.2.and.rcobj))
     *                     write (lu,1031) 'Chi', den**2/mode 
                        if (rwch.or.(lsqchm.eq.3.and.rcobj))  
     *                     write (lu,1031) 'wChi', 
     *                                      (den/(emode(id,ld)*1d2))**2

                        write (lu,1025) 'Mode Misfit:', mscor(j)

                     end if

                     write (lu,1025) 
     *                        'Composition Misfit zero by definition.'

                  end if
c                                 end of mpred loop
               end do

               write (lu,1011) 

            end if

            if (isat.gt.0) then

               write (lu,1160)
               write (lu,1020) (cname(k), k = 1, lbulk)
c                                 take care of composants with amt > 0
               do j = 1, istsat

                  jd = jds(j)
                  mode = props(1,jd)*props(16,jd)/psys(1)*1d2

                  tot = 0d0

                  do k = 1, lbulk
                     tot = tot + pcomp(k,jd)
                  end do

                  write (lu,1050) pname(jd), mode, 
     *                            (pcomp(k,jd)/tot, k = 1, lbulk)

               end do
c                                 next composants with amt < 0
               do j = 1, isat

                  if (xskp(j).le.0d0) then

                     jd = -skp(j)
                     tot = 0d0

                     do k = 1, lbulk
                        tot = tot + cp(k,jd)
                     end do

                     write (lu,1050) names(jd)(1:nblen(names(jd)))//'*',
     *                            nopt(7), (cp(k,jd)/tot, k = 1, lbulk)
                  end if

               end do

               if (isat.gt.istsat) write (lu,1170)

            end if

            if (mextra.gt.0) then 
c                                 phases predicted but not observed
               write (lu,1040)
               write (lu,1020) (cname(k), k = 1, lbulk)

               do j = 1, mextra

                  jd = idextr(j)

                  mode = props(1,jd)*props(16,jd)/psys(1)*1d2

                  tot = 0d0

                  do k = 1, lbulk
                     tot = tot + pcomp(k,jd)
                  end do

                  write (lu,1050) pname(jd), mode, 
     *                            (pcomp(k,jd)/tot, k = 1, lbulk)
c                                 end of miss loop
               end do

             end if
c                                 phases observed but not predicted
             nmiss = 0

             do j = 1, xptnph(id)

               if (jmiss(j)) then
                  nmiss = nmiss + 1
                  call getnam(bdname(nmiss),xptids(1,j))
               end if

            end do

            if (nmiss.gt.0) write (lu,1060) (bdname(j), j = 1, nmiss)

c                                 complete bulks and chemical potentials
            write (lu,1080) (cname(j), j = 1, lbulk)

            tot = 0d0

            do j = 1, lbulk
               tot = tot + fbulk(j)
            end do
c                                 bulk:
            write (lu,1081) 'Bulk', (fbulk(j)/tot, j = 1, lbulk)
c                                 phases:
            do j = 1, mpred
c                                 jd/kmin points to the phase in the predicted assemblage
               jd = kmin(j)

               tot = 0d0
               do k = 1, lbulk
                  tot = tot + pcomp(k,jd)
               end do

               write (lu,1081) pname(jd), 
     *                         (pcomp(k,jd)/tot, k = 1, lbulk)

            end do

            write (lu,1082) ('-----', j = 1, lbulk)
            write (lu,1083) 'Potentials',(mu(j), j = 1, lbulk)
            write (lu,1010)

            if (nmcov.and.covar(1).ne.0d0) then
c                                 write covariance/correlation matrix
               write (lu,1180) (vnm(j), j = 1, nparm)

               k = 0

               do j = 1, nparm
                  write (lu,1190) vnm(j),covar(k+1:k+nparm)
                  k = k + nparm
               end do

               write (lu,1200) (dsqrt(covar(jkdiag(j))), j = 1, nparm)

               write (lu,1210)

            else if (nmcov) then
c                                 covariance matrix calculation failed
               write (lu,1220)

            end if

            write (lu,'(/,80(''-''))')
c                                 end double unit print loop
         end do

      end if

c                       write (*,'(a,1pg13.6,1x,i6),a,$)')
c     *         'Misfit: ', miss + extra + comp + scmode, optct, char(13)
c     if (.not.fprint)  write (*,'(a,1pg13.6,1x,i6)')
c    *         'Misfit: ', xptscr, optct


1000  format (/,'The following observed phases are predicted:*')
1010  format (/,'**mole fractions')
1011  format (/,'*mole fractions relative to measured components, compl'
     *         ,'ete compositions below')  
1020  format (t16,'vol %',t26,*(a,5x))
1021  format (/,a,t16,'vol %',t26,*(a,5x))
1025  format (2x,a,t24,1pg13.6e2,2x,a)
1030  format (2x,a,t14,f7.3,3x,*(f8.6,2x))
1031  format (2x,a,t14,1pg8.1e1)
1035  format (2x,a,t24,*(f8.6,2x))
1036  format (2x,a,t23,*(1pg9.1e2,1x))
1037  format (2x,a,t14,1pg8.1e1,1x,*(1pg9.1e2,1x))
1040  format (/,'The following predicted phases** are not observed:',/)
1045  format (2x,a,13x,*(f7.3,1x))
1050  format (a,t14,f7.3,1x,*(2x,f8.6))
1060  format (/,'The following observed phases are not predicted:',//,
     *           10(2x,a))
1080  format (/,'Complete compositions** and chemical potent'
     *         ,'ials (J/mol)',//,t13,*(3x,a,2x))
1081  format (a,t12,*(2x,f8.6))
1082  format (t13,*(3x,a,2x))
1083  format (a,t13,*(1x,f9.0))


*1080  format (/,5x,'       Effective Bulk*      ',
*     *            '  Chemical Potentials (J/mol)')
1090  format (a5,10x,f8.3,20x,g13.6)
1100  format (/,'Components of the Misfit Score',//,5x,
     *          'Predicted compositions',t34,1pg13.6e1,/,5x,
     *          'Extraneous predicted phases:',t34,1pg13.6e1,/,5x,
     *          'Missed observed phases:',t34,1pg13.6e1)
1105  format (5x,'Predicted modes:',t34,1pg13.6e1,/)
1130  format (/,'Converged but because the predicted assemblage was in',
     *          'complete the model will not',/,'be counted. To includ',
     *          'e incomplete predictions set no_miss to miss_ok',//,
     *          80('-'))
1150  format (/,'This model fits all data within observational',
     *          ' uncertainty')
1151  format (/,'This model DOES NOT fit all data within observational',
     *          ' uncertainty')
1160  format (/,'The following phases are in equilibrium with the obse',
     *          'rved assemblage as determined',/,'by the component sa',
     *          'turation constraints:',/)
1170  format (/,'*The indicated phase is stable but undersaturated wit',
     *          'h respect to the effective bulk composition.')
1180  format (/,'Misfit covariance and correlation matrices* and ',
     *          'standard deviations:',/,
     *        10x,*(3x,a8,3x))
1190  format (2x,a8,*(1x,1pg13.6))
1200  format (/,2x,'std_dev ',*(1x,1pg13.6))
1210  format (/,
     *        '*Covariance elements are on and below the left diagonal',
     *        ', corellation elements are',/,'above the diagonal, on-d',
     *        'iagonal correlations are 1 by definition.')
1220  format (/,'**warning ver504** covariance matrix calculation fai',
     *       'led, increase expansion_criterion',/,'or set print_leve',
     *       'l > 0 for additional diagnostic information')

      pause

      end

      subroutine savbst (objf,bstlik,bstlx,ppb,x0inv,bstlx0,
     *                   bstbx0,bstbay,bstbx, n,i,ibest,jbest,igood)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, n, m, j, igood, ibest, jbest

      double precision objf, bstlik, ppb, x0inv(*), bstlx0(*), pre,
     *                 bstbx0(*), bstbx(*), bstbay, bstlx(*), fac
c----------------------------------------------------------------------
      igood = igood + 1
c                                 compute ss of parameter deviations
      ppb = 1d0
c                                 pre-exponential term factor
      pre = 0.398942 * ra2zs
c                                 multiplier on (res/range)
      fac = -0.5d0 * ra2zs**2

      if (ptonly) then 
         m = 2
      else
         m = n
      end if
c                                 center "bayesian" score in interval
      do j = 1, m
c                                 compute the product of the priors
         if (pdelta(j).ne.0d0) then

            ppb = ppb * pre/pdelta(j) * dexp ( fac *
c                                 (resiual / range) ^ 2
     *                            ( (xinv(j) - plow(j) - pdelta(j)/2d0) 
     *                              / pdelta(j) )**2 )

         end if

      end do
c                                 best "bayesian" score
c                                 Likelihood = exp(-objf) strictly 
c                                 for wChi only
      bay = exp(-objf) * ppb
      lik = objf

      if (bay.gt.bstbay) then

         jbest = i
         bayfit = jfit
         bstbay = bay
         bstbx(1:n) = xinv(1:n)
         bstbx0(1:n) = x0inv(1:n)
         bstbco = covar

      end if 
c                                 save min misfit result
      if (objf.lt.bstlik) then

         ibest = i
         likfit = jfit
         bstlik = objf
         bstlx(1:n) = xinv(1:n)
         bstlx0(1:n) = x0inv(1:n)
         bstlco = covar

      end if

      end

      subroutine prtsum (objf,bstlik,ppb,x0inv,
     *                   bstbay,n,ibest,jbest,igood,jcount,
     *                   icount,strg1)
c----------------------------------------------------------------------
c print simple summary of converged model to: console and *.out
c
c objf - should be current model final coordinates
c trnsxy must be called to do the inverse transformation x -> xinv
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character strg1*(*), text1*100

      integer n, j, igood, ibest, jbest, lu, nblen, jcount, icount

      double precision objf, ppb, x0inv(*), ox(l12), bstbay, bstlik

      external nblen
c----------------------------------------------------------------------
c                                 write notice and stats to console
            lu = 6

            ox(1:n) = xinv(1:n)

            do j = 1, 2
c                                 loop to write to console, and *.out
               write (lu,'(80(''-''))')
               write (text1,'(a,a,i4,a)') strg1(1:nblen(strg1)),',',
     *                                    igood,' successes so far.'
               call deblnk (text1)
               write (lu,'(a)') text1(1:nblen(text1))

               if (iquad.gt.0.and.jcount.gt.0) then
                  write (text1,1110) icount, jcount
               else
                  write (text1,1100) icount
               end if 

               call deblnk (text1)
               write (lu,'(/,a)') text1(1:nblen(text1))

               write (lu,1120) objf

               write (text1,1125) bstlik, ibest
               call deblnk (text1)
               write (lu,'(a)') text1(1:nblen(text1))

               write (lu,1130) ppb, bay
               write (text1,1140) bstbay, jbest
               call deblnk (text1)
               write (lu,'(a,/)') text1(1:nblen(text1))

               write (lu,1010) vnm(1:n)
               write (lu,1020) x0inv(1:n)
               write (lu,1030) ox(1:n)

               lu = n6

             end do
c                                 reset x, just in case
      xinv(1:n) = ox(1:n)

1000  format (a,1x,*(1pg13.6,1x),a)
1010  format (20x,*(4x,a8,2x))
1020  format ('Initial coordinates: ',*(1pg13.6,1x))
1030  format ('Final coordinates:   ',*(1pg13.6,1x))
1050  format (/,'Number of function evaluations: ',i5,', igood = ',i3,/)
1090  format ('Try ',i5,', successes so far = ',i5)
1100  format ('Misfit function evaluations this try => ',i5)
1110  format ('Misfit function evaluations this try => ',i5,
     *          ', + ',i5,' evaluations for quadratic surface fitting')
1120  format (/,'Misfit this try = ',1pg13.6)
1125  format ('Best Misfit so far = ',1pg13.6,' obtained on Try ',i5)
1130  format (/,'Prior probability, PP = ',1pg13.6,/,
     *          'Bayes score, PP * exp(-Misfit) = ',g12.6)
1140  format ('Best Bayes score so far = ',1pg13.6,' obtained on Try ',
     *        i5)

      end

      subroutine mcopts
c----------------------------------------------------------------------
c set default parameters for mc_fit and read option file
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier, i, nblen, oquad

      logical inum, rnum

      character key*28, val*3, nval1*12, nval2*12,
     *          nval3*12, strg*40, strg1*40, fnm*80

      double precision x

      external nblen

      save oquad
c----------------------------------------------------------------------
c                                 set all parameters to their defaults
c                                 -------------------------------------
c                                 formerly imc file section 1 variables:
c                                 -------------------------------------
c                                 molar_composition:
c                                 T - assume molar units for input compositions
c                                 F - assume molar units for input compositions
      lmass = .false.
c                                 kiso:
c                                 T - convert mass input to molar and output, terminate before inversion.
c                                 F - do nothing
      kiso  = .false.
c                                 thermobarometry = .not.invxpt, invxpt:
c                                 T - invert experimental data
c                                 F - thermobarometric inversion
      invxpt = .false.
c                                 MINIM_algorithm = .not.mcgrid, mcgrid:
c                                 T - grid search algoirthm
c                                 F - minim algorithm
      mcgrid = .false.
c                                 Bayes:
c                                 T - pseudo-bayesian best model criterion
c                                 F - misfit best model criterion
      bayes = .false.
c                                 prior_probability_PT_only
c                                 T - assume prior probabilities are known only for P and T
c                                 F - assume all ranges are based on prior probabilities
      ptonly = .true.
c                                 ra2zs:
c                                 width of range in terms of the parameter standard deviations
      ra2zs = 4d0
c                                 un2ft:
c                                 multiplier on uncertainty for assessing whether a model fits within error
      un2ft = 2d0
c                                 hot_start:
c                                 F - do uncertainty analysis around best central model
c                                 T - do uncertainty analysis around a user specified model coordinate
      mchot = .false.
c                                 GRH:
c                                 T - use georges normalizations, etc
c                                 F - use my normalization
      GRH = .false.
c                                 model_output:
c                                 0 - output all converged models
c                                 1 - output only improved models
c                                 2 - output only the best model
      bstout = 0
c                                 vital_sign:
c                                 T - output optimization counter
c                                 F - quiet
      vital = .true.
c                                 seed:
c                                 T - seed random number generator => succesive runs of mc_fit will generate different results
c                                 F - don't seed random number generator => succesive runs of mc_fit will yield the same result
      seed = .false.
c                                 must_fit_all_data:
c                                 T - reject models that do not fit all analytical data within uncertainty*uncertainty_multiplier
c                                 F - accept models regardless of whether they fit analytical data
      mustft = .false.
c                                 no_missing_phases:
c                                 T - reject models that fail to predict an observed phase
c                                 F - accept models regardless of whether they predict the observed phases
      nomiss = .true.
c                                 relative_error:
c                                 T - assume input uncertainties are relative error
c                                 F - assume input uncertainties are absolute error
      relerr = .true.
c                                 T - use central model for uncertainty analysis
c                                 F - use random initial guesses for uncertainty analysis
      newstt = .false.
c                                 1 - total uncertainty
c                                 2 - analytical uncertainty only
c                                 3 - thermodynamic uncertainty only
      invunc = 1
c                                 Nelder-Meade covariance matrix
      nmcov = .false.
c                                 add central model Nelder-Meade variance to 
c                                 perturbation analysis variance
      unplus = .true.
c                                 -------------------------------------
c                                 output normalized residuals
      rnorm = .true.
c                                 output unormalized residuals
      runrm = .false.
c                                 output OBJ contribution per component
      rcobj = .true. 
c                                 output LSQ contribution per component
      rlsq = .false.
c                                 output Chi2
      rchi = .false.
c                                 output wChi
      rwch = .false.
c                                 -------------------------------------
c                                 formerly imc file section 2:
c                                 Objective function mode, tolerances, and weights
c                                 -------------------------------------
c                                 1 - least square compositional scoring
c                                 2 - Chi square scoring
c                                 3 - weighted Chi square
      lsqchi = 2
c
      objnam(1) = 'LSQ'
      objnam(2) = 'Chi'
      objnam(3) = 'wChi'
c                                 and for modes (same scheme):
      lsqchm = 2
c                                 maximum value of the objective function, 
c                                 larger values trigger an error during search (or should)
      oktol = 1d5
c                                 maximum misfit for an acceptable result
      mxobjf = oktol
c                                 the variance of the objective function for an acceptable
c                                 solution
      invtol = 1d-4
c                                 sftmx = softmax * scale for simplicial variables
      sftmax = .true.
c                                 sftpot = soft interior for potanetials
      sftpot = .true.
c                                 sftc = soft interior for unmeasured components
      sftc = .true.
c                                 sftdel = width of the approximately linear region of the
c                                 softmax/soft-interion transformations
      sftdel = 6.0
c                                 weighting factor for the compositional term of the 
c                                 objective function. for wChi this is reset to 1d-4
c                                 if the user hasn't overridden wcomp
      wcomp = 1d0
c                                 weighting factor for the modal term of the 
c                                 objective function
      wmode = 1d0
c                                 weighting factor for the extraneous phase term of the 
c                                 objective function
      wextra = 1d4
c                                 weighting factor for the missing phase term of the 
c                                 objective function
      wmiss = 1d3
c                                 -------------------------------------
c                                 formerly imc file section 3:
c                                 Monte-Carlo logistical variables
c                                 -------------------------------------
c                                 number of starting guesses used to find the best fit central model
      mtry = 100
c                                 number of perturbations used to evaluate uncertainty of central model
      nunc(1) = 100
c                                 george's extra grid search parameter
      nunc(2) = 5
c                                 number of tries for each perturbation
      ptry = 1
c                                 -------------------------------------
c                                 formerly imc file section 4:
c                                 Nelder-Mead MINIM algorithm variables
c                                 -------------------------------------
c                                 criterion used to decide whether final simplex is expanded 
c                                 before quadratic fitting, should be at least 1000 x rounding error, 
c                                 generally a lot more. 
      simplx = 1d-3
c                                 fractional step size
      frac = 0.1
c                                 how often to check for convergence
      conchk = 10
c                                 print control: 0 => normal, < 0 => silent, > 0 => verbose
      jprint = -1
c                                 quadratic surface fitting to improve location of the minimum, 
c                                 0 => don't fit, 1 => fit
      iquad = 0
      oquad = 0
c                                 max # of objective evaluations per try, recommended nparm*1024
      kcount = 1000
c----------------------------------------------------------------------
c                                 read mc_fit option file name from *.imc file (on n8)
      call redfnm (n8,ier,fnm)

      if (ier.ne.0) call errdbg ('expecting the mc_fit option file na'//
     *              'me or ''none'', found: '//fnm)

      if (fnm.eq.'none') then

         if (mcfrst) write (*,1010) tfname(1:nblen(tfname))
         return

      else 

         open (n9, file = fnm, iostat = ier, status = 'old')

         if (ier.ne.0) then 
            call errdbg ('could not find the specified mc_'//
     *                   'fit option file: '//fnm)
         else 
            if (mcfrst) write (*,1050) fnm(1:nblen(fnm))
         end if

      end if
c----------------------------------------------------------------------
c                                 sift through the option file
      do

         call redcd1 (n9,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 ier ne 0 = eof
         if (ier.ne.0) exit
c                                 cycle on default specification
         if (strg.eq.'default') cycle
c                                 set print flags
         inum = .false.
         rnum = .false.
c                                 -------------------------------------
c                                 if here we have a keyword and value
         if (key.eq.'molar_composition_input') then

            if (val.eq.'F') lmass = .true.

         else if (key.eq.'kiso') then

            if (val.eq.'T') kiso = .true.

         else if (key.eq.'softmax') then

            if (val.eq.'F') sftmax = .false.

         else if (key.eq.'soft_interior_potentials') then

            if (val.eq.'F') sftpot = .false.

         else if (key.eq.'soft_interior_C') then

            if (val.eq.'F') sftc = .false.

         else if (key.eq.'soft_delta') then

            read (strg,*) sftdel
            rnum = .true.

         else if (key.eq.'thermobarometry') then

            if (val.eq.'F') invxpt = .true.

         else if (key.eq.'MINIM_algorithm') then

            if (val.eq.'F') mcgrid  = .true.

         else if (key.eq.'Bayes') then

            if (val.eq.'T') bayes = .true.

         else if (key.eq.'prior_probability_PT_only') then

            if (val.eq.'F') ptonly = .false.

         else if (key.eq.'range_to_sigma') then

            read (strg,*) ra2zs
            rnum = .true.

         else if (key.eq.'uncertainty_multiplier'
     *        .or.key.eq.'sigma_multiplier') then
c                                 uncertainty_multiplier deprecated in 7.2.6
            read (strg,*) un2ft
            rnum = .true.

         else if (key.eq.'must_fit_all_data') then

            if (val.eq.'T') mustft = .true.

         else if (key.eq.'perturbation_hot_start') then
c                                 ~mchot => do perturbation analysis automatically after finding best central model
c                                 mchot => do perturbation analysis around a user specified inversion coordinate
            if (val.eq.'T') mchot = .true.

         else if (key.eq.'GRH') then

            if (val.eq.'T') GRH = .true.

         else if (key.eq.'model_output') then

            if (val.eq.'imp') then 

               bstout = 1

            else if (val.eq.'bes') then 

               bstout = 2

            else if (val.eq.'all') then 

               bstout = 0
           
            else 

               call badopt (key,strg)

            end if

         else if (key.eq.'residuals') then

            if (val.eq.'T') runrm = .true.

         else if (key.eq.'residuals_standardized') then

            if (val.eq.'F') rnorm = .false.

         else if (key.eq.'score_misfit') then 

            if (val.eq.'F') rcobj = .false.

         else if (key.eq.'score_LSQ') then 

            if (val.eq.'T') rlsq = .true.

         else if (key.eq.'score_Chi') then 

            if (val.eq.'T') rchi = .true.

         else if (key.eq.'score_wChi') then 

            if (val.eq.'T') rwch = .true.

         else if (key.eq.'vital_sign') then

            if (val.eq.'F') vital = .false.

         else if (key.eq.'seed') then

            if (val.eq.'T') seed = .true.

         else if (key.eq.'no_missing_phases') then

            if (val.eq.'F') nomiss = .false.

         else if (key.eq.'relative_error') then

            if (val.eq.'F') relerr = .false.

         else if (key.eq.'best_central_model') then
c                                 ~newstt => do perturbation analysis around best central model
c                                 newstt => do perturbation analysis using random starting coordinates
            if (val.eq.'F') newstt = .true.

         else if (key.eq.'uncertainties') then

            if (val.eq.'all') then 

               invunc = 1

            else if (val.eq.'ana') then

               invunc = 2

            else if (val.eq.'the') then

               invunc = 3

            else 

               call badopt (key,strg)

            end if

         else if (key.eq.'Nelder-Mead_covariance') then

            if (strg.eq.'T') nmcov = .true.

         else if (key.eq.'uncertainty_plus') then

            if (strg.eq.'F') unplus = .false.

         else if (key.eq.'objective_composition'.or.
     *            key.eq.'misfit_composition') then

            if (strg.eq.'LSQ') then

               lsqchi = 1

            else if (strg.eq.'Chi') then

               lsqchi = 2

            else if (strg.eq.'wChi') then

               lsqchi = 3

            else 

               call badopt (key,strg)

            end if

         else if (key.eq.'objective_mode'.or.key.eq.'misfit_mode') then

            if (strg.eq.'LSQ') then

               lsqchm = 1

            else if (strg.eq.'Chi') then

               lsqchm = 2

            else if (strg.eq.'wChi') then

               lsqchm = 3

            else 

               call badopt (key,strg)

            end if

         else if (key.eq.'max_misfit_value'.or.
     *            key.eq.'max_OBJF_value') then

            read (strg,*) oktol
            rnum = .true.

         else if (key.eq.'misfit_filter_value') then

            read (strg,*) mxobjf
            rnum = .true.

         else if (key.eq.'convergence_threshold') then

            read (strg,*) invtol
            rnum = .true.

         else if (key.eq.'w_composition') then

            read (strg,*) wcomp
            rnum = .true.
c           lwcomp = .true.

         else if (key.eq.'w_extra') then

            read (strg,*) wextra
            rnum = .true.

         else if (key.eq.'w_missing') then

            read (strg,*) wmiss
            rnum = .true.

         else if (key.eq.'number_of_tries') then

            read (strg,*) mtry
            inum = .true.

         else if (key.eq.'number_of_perturbations') then

            read (strg,*) nunc(1)
            inum = .true.

         else if (key.eq.'number_of_perturbation_tries') then

            read (strg,*) ptry
            inum = .true.

         else if (key.eq.'expansion_criterion') then

            read (strg,*) simplx
            rnum = .true.

         else if (key.eq.'step_size') then

            read (strg,*) frac
            rnum = .true.

         else if (key.eq.'convergence_testing') then

            read (strg,*) conchk
            inum = .true.

         else if (key.eq.'print_level') then

            read (strg,*) jprint
            inum = .true.

         else if (key.eq.'quadratic_fitting') then

            if (val.eq.'T') iquad = 1
            oquad = 1

         else if (key.eq.'max_misfit_evaluations'.or.
     *            key.eq.'max_OBJF_evaluations') then

            read (strg,*) kcount
            inum = .true.

         else if (key.ne.'|') then

            call error (72,nopt(1),iopt(1),key//' is not a valid mc_fit'
     *                 //' option file keyword and must be deleted '
     *                 //'or corrected.')

         end if
c                                 reset wcomp "default" is wChi
c        if (.not.lwcomp.and.lsqchi.eq.3) wcomp = 1d-4
c                                 print modified options:
         if (mcfrst) then
            if (inum) then
               read (strg,*) i
               write (*,1030) key, i
            else if (rnum) then 
               read (strg,*) x
               write (*,1040) key, x
            else
               write (*,1020) key, strg
            end if
         end if

      end do

      close (n9)
c                                 set dependent options
      if (nmcov.and.mcfrst) then
c                                 nelder-mead covariance currently 
c                                 requires quadratic surface fitting:
         if (iquad.eq.0) then 
            iquad = 1
            key = 'quadratic_fitting'
            write (*,1070) key,'T*'
            write (*,1060)
         end if

      else if (.not.mcfrst) then
c                                 doing uncertainty analysis, shut
c                                 off quadratic fitting unless it was
c                                 explicitly set.
         iquad = oquad

      end if

1000  format ('Context specific options are echoed in: ',a)
1010  format (/,'**warning** ',a,' does not specify an mc_fit option',/,
     *       ' all options will be set'/,'to default, refer to mc_fit_',
     *       'option.dat for default option values',/)
1020  format ('mc_fit option: ',a,' has been set to ',a)
1030  format ('mc_fit option: ',a,' has been set to ',i6)
1040  format ('mc_fit option: ',a,' has been set to ',g12.6)
1050  format (/,'Reading mc_fit options from: ',a,/)
1060  format (/,'*Quadratic fitting is costly. Do not use it ',
     *      '(i.e., set Nelder-Mead_covariance',/,' to F) unless you',
     *      ' want the Nelder-Mead parameter correlation matrix and/or',
     *      /,' Nelder-Mead uncertainty estimates.',/,80('-'))
1070  format ('mc_fit option: ',a,' has been automatically reset to ',a)

      end 

      subroutine badopt (key,strg)
c----------------------------------------------------------------------
c error diagnostic for bad values for an option keyword
c----------------------------------------------------------------------
      implicit none

      character key*28, strg*40
c----------------------------------------------------------------------
      write (*,1000) strg, key

1000  format (/,'**error ver324** Invalid value: ',a,
     *        /,'for mc_fit option: ',a,/)

      call errpau

      end 

      subroutine reddqf (id)
c----------------------------------------------------------------------
c read solution specific non-excess function parameters for mc_fit data
c inversionts. 

c        DQF(name) or dqf(name)

c        jd - endmember pointer
c        id - solution model pointer
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, done

      integer j, id

      integer ibeg, ier, iscan, imax, iscnlt, jd, kend

      character name*8

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      external iscnlt, iscan

      integer length,com
      character chars*1, card*400
      common/ cst51 /length,com,chars(400),card

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c----------------------------------------------------------------------
      ier = 0
      done = .false.

      do
c                                 find brackets
         ibeg = iscan (1,com,'(') + 1
         imax = iscan (1,com,')') - 1

         if (ibeg.gt.com.or.imax.gt.com) call errdbg ('invalid DQF '//
     *                                            'format for '//tname)

         call readnm (ibeg,kend,imax,ier,name)

         if (ier.ne.0) call errdbg ('invalid DQF format for '//tname)

         call matchj (name,jd)

         if (jd.eq.0) call errdbg ('invalid endmember '//name//
     *                             ' in DQF specified for '//tname//
     *                             ' in the *.imc file')
         ok = .false.

         do j = 1, jdqf(id)
c                                 check for a match with internal dqf
            if (-jd.eq.jend(id,2+jndq(j,id))) then
               ok = .true.
               exit
            end if

         end do

         if (.not.ok) call errdbg ('endmember '//name//' is not in '//
     *      'or has been ejected from '//tname)
c                                 the dqf exists
         mctrm(mcsol) = mctrm(mcsol) + 1
c                                 save the term pointer
         mcj(mcsol,mctrm(mcsol)) = j
c                                 find the free coefficients
         do 
c                                 loop to read free coefficients in the DQF
            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 key may be either "parameter"
c                                 or a new DQF "DQF(name)" or "end_list"
            if (key.eq.'end_list') then

               done = .true.
               exit

            else if (chars(1).eq.'D'.or.chars(1).eq.'d') then
c                                 starting a new DQF
               exit

            else if (key.ne.'parameter') then 

               call errdbg ('invalid DQF coefficient '//
     *                      'specification for '//tname)

            end if
c                                 count the coefficient
            mccoef(mcsol,mctrm(mcsol)) = 
     *                   mccoef(mcsol,mctrm(mcsol)) + 1

c                                 read the parameter index
            read (strg1,*) 
     *           mccoid(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)))

            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 read the parameter range
            read (strg1,*)
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),1),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),2)

          end do
c                                 start a new term or exit
          if (done) exit

      end do

      end

      subroutine redcpd (id)
c----------------------------------------------------------------------
c read inversion parameters for compunds, called by mcexpt.

c        id - endmember pointer
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, nblen, ier

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      double precision err

      external nblen

      integer length,com
      character chars*1, card*400
      common/ cst51 /length,com,chars(400),card

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)
c----------------------------------------------------------------------

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

1010  format (/,a,1x,a)
1020  format (/,a,1x,'has ',a,' magnetic transitions')

      end 

      subroutine parmct (n6out,lun)
c----------------------------------------------------------------------
c count the parameters and echo them if requested
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical n6out

      integer i, j, k, id, jd, nblen, lun

      character key*22, val*3, strg*40

      external nblen

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer length,com
      character chars*1, card*400
      common/ cst51 /length,com,chars(400),card
c----------------------------------------------------------------------
      if (n6out) write(lun,*) 'Inverted parameters:'

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
                  write(lun,1050) nparm+j,names(id),key(1:nblen(key)),
     *               cprng(i,j,1), cprng(i,j,2)
               else
                  write(lun,1050) nparm+j,names(id),key(1:nblen(key)),
     *               cprng(i,j,1), cprng(i,j,2), int(cprng(i,j,3))
               end if
            end do
            
         end if
c                                 not so hard to count
         nparm = nparm + mcpct(i)

      end do
c                                 solution parameter groups:
      do i = 1, mcsol

         id = mcids(i)

         if (mcspid(i).eq.1) then
c                                 excess term group
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
     *                     ' '// names(jend(id,2+2))
                     else if (extyp(id).eq.0) then
c                                 Margules names depend on order of solution
                     strg = ' '

                     do k = 1, rko(j,id)
                        strg(1+nblen(strg):) = ' ' //
     *                      names(jend(id,2+jsub(k,j,id)))
                     end do

                     strg = strg(2:)

                  else if (extyp(id).eq. 2) then

                     call errdbg ('no counting for van laar')

                  end if

                  do k = 1, mccoef(i,j)
c                                 which parameters are getting varied?
                     if (random(1).lt.2) then
                        write(lun,1060) nparm+k, fname(id),
     *                     val(1:nblen(val)), strg(1:nblen(strg)), j,
     *                     mccoid(i,j,k),
     *                     sprng(i,j,k,1), sprng(i,j,k,2)
                     else
                        write(lun,1060) nparm+k, fname(id),
     *                     val(1:nblen(val)), strg(1:nblen(strg)), j,
     *                     mccoid(i,j,k),
     *                     sprng(i,j,k,1), sprng(i,j,k,2),
     *                     int(sprng(i,j,k,3))
                     end if

                  end do

               end if

               nparm = nparm +  mccoef(i,j)

            end do

         else if (mcspid(i).eq.2) then

            do j = 1, mctrm(i)
c                                 grunt work to count and name the term
               if (n6out) then

                  jd = jend(id,2+mcj(i,j))

                  strg = 'DQF('//names(jd)(1:nblen(names(jd)))//')'

                  do k = 1, mccoef(i,j)
c                                 which parameters are getting varied?
                     if (random(1).lt.2) then
                        write(lun,1070) nparm+k, fname(id),
     *                     strg(1:nblen(strg)), j,
     *                     mccoid(i,j,k),
     *                     sprng(i,j,k,1), sprng(i,j,k,2)
                     else
                        write(lun,1070) nparm+k, fname(id),
     *                     strg(1:nblen(strg)), j,
     *                     mccoid(i,j,k),
     *                     sprng(i,j,k,1), sprng(i,j,k,2),
     *                     int(sprng(i,j,k,3))
                     end if
c                                 end optional print loop
                  end do

               end if

               nparm = nparm +  mccoef(i,j)

            end do

         end if

      end do

      if (n6out) write (lun,'(80(''-''))')

      if (nparm.eq.0) call errdbg ('no free parameters! no free lunch!')

1050  format (i2,2(1x,a),2(1x,1pg12.5),1x,i3)
1060  format (i2,1x,a,1x,a,1h(,a,1h),i2,' parameter',i2,2(1x,1pg12.5),
     *        1x,i3)
1070  format (i2,1x,a,1x,a,i2,' parameter',i2,2(1x,1pg12.5),
     *        1x,i3)

      end
