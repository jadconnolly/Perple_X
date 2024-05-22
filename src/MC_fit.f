c Please do not distribute any part of this source.
 
c Copyright (c) 1987-2023 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved. 

      program MC_fit
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical invprb

      integer conchk, kcount, iprint, iquad, ntry

      double precision tol, simplx, frac

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------- 
c                                 MC_inv uses the MEEMUM iam flag value
      iam = 2
c                                 initialization, read files etc.
      call iniprp
c                                 open inversion problem file
      call opnimc (invprb,ntry,tol,simplx,frac,conchk,iprint,iquad,
     *             kcount)

      if (invprb) then

         call invxpt (ntry,tol,simplx,frac,conchk,iprint,iquad,kcount)

      else

         call invptx (ntry,tol,simplx,conchk,iprint,iquad,kcount)

      end if

      end

      subroutine opnimc (invprb,ntry,tol,simplx,frac,conchk,iprint,
     *                   iquad, kcount)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical invprb

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40

      integer conchk, kcount, iprint, iquad, ntry, ier

      double precision tol, simplx, frac
c----------------------------------------------------------------------- 
c                                 open inversion problem file
      call mertxt (tfname,prject,'.imc',0)
      open (n8,file=tfname,status='old',iostat=ier)

      if (ier.ne.0) call errdbg 
     *   ('can''t open assemblage composition file: '//tfname)
c                                 read problem type flag
      invprb = .true.

      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      if (key.ne.'T'.and.key.ne.'t') invprb = .false.
c                                 make new seed for random number generator
      random = .true.

      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      if (key.ne.'T'.and.key.ne.'t') random = .false.
c                                 number of starting guesses used for each
c                                 Nelder-Meade inverse problem
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) ntry
c                                 error evaluation loop counter
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) nunc
c                                 read Nelder-Meade parameters
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) oktol
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) tol
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) simplx
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) frac
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) conchk
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) iprint
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) iquad
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) kcount
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) wcomp
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) wextra
      call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
      read (key,*) wmiss

      end

      subroutine invptx (ntry,tol,simplx,conchk,iprint,iquad,kcount)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, n, conchk, kcount, icount, ifault, iprint,
     *        iquad, j, igood, ntry, ibest

      logical readyn

      character amount*6

      double precision var(l2+k5), objf, mcobjf, x(l2+k5), 
     *                 tol, step(l2+k5), simplx, bstobj

      external readyn, mcobjf, mcobj1

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      double precision atwt
      common/ cst45 /atwt(k0) 

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      character*5 cname
      common/ csta4 /cname(k5)

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------- 
c                                 iwt is set by input, it is only used below to determine
c                                 whether to convert weight compositions to molar. the 
c                                 computations are done solely in molar units. 
      amount = 'molar '

      if (iwt.eq.1) amount = 'weight'
c                                 read phase compositions
         call mccomp

         n = ipot+mphase-1

         x(1:n) = 0.5d0

         tol = 1d-8
         step(1:n) = 5d-1
         conchk = 10
         iprint = -1
         iquad = 1
         ibest = 0
         igood = 0
         bstobj = 1d99
         simplx = 1d-8
         ntry = 100
c                                 max number of objf evaluations
         kcount = 10000
c                                 initialize drand
         call random_seed

         do i = 1, ntry

            write (*,1080) x(1:n)

            call mcsetv (x)

            write (*,1090) (v(jv(j)),j= 1, ipot)
c                                 MINIM R O'neil
c                                 www.scilab.org/sites/default/files/neldermead.pdf
c                                 code: lib.stat.cmu.edu/apstat/47
            call minim (x, step, n, objf, kcount, iprint, tol, 
     *                  conchk, iquad, simplx, var, mcobj1, icount, 
     *                  ifault, oktol)

c https://people.sc.fsu.edu/~jburkardt/f77_src/asa047/asa047.html
c           call nelmin (mcobjf, n, xini, xopt, objf, tol, step, conchk,
c    *                kcount, icount, numres, ifault)

            if (ifault.ne.0) then 

               write (*,1020) ifault,igood,icount

            else

               igood = igood + 1

               if (objf.lt.bstobj) then 
                  ibest = i
                  bstobj = objf
               end if

               write (*,1030) x(1:n)
               call mcsetv (x)
               write (*,1040) (v(jv(j)),j= 1, ipot)
               write (*,1050) icount, igood
            end if
c                               new starting point
            do j = 1, n
               call random_number (x(j))
            end do

         end do

      stop


1080  format ('Initial normalized coordinates: ',
     *        5(g12.6,1x))
1090  format (10x,'initial potentials: ',5(g12.6,1x))
1020  format ('Minimization failed IFAULT = ',i3,' igood = ',i3)
1030  format ('Final coordinates: ',5(g12.6,1x))
1040  format ('Final potentials: ',5(g12.6,1x))
1050  format ('Number of function evaluations: ',i5,
     *        ' igood = ',i3,' icount = ',i5)
1000  format (/,'Interactively enter bulk compositions (y/n)?',/,
     *          'If you answer no, MEEMUM uses the bulk composition',
     *         ' specified in the input file.',/)
1010  format (/,'Enter value of bulk compositional variable X(C',i1,'):'
     *       )
1060  format (/,'Enter ',a,' amounts of the components:')
1070  format (/,'Enter (zeroes to quit) ',7(a,1x))

      end 

      subroutine invxpt (ntry,tol,simplx,frac,conchk,iprint,iquad,
     *                   kcount)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, conchk, kcount, iprint,
     *        iquad, j, k, ntry

      logical invprb

      double precision tol, simplx, frac, bstx(l2+k5), 
     *                 x(200,l2+k5), sx(l2+k5), ss(l2+k5)
c----------------------------------------------------------------------- 
c                                 best model statistics for error evaluation
      call mertxt (tfname,prject,'.bst',0)
      open (n7,file=tfname)
c                                 initialize drand
      if (random) call random_seed
c                                 get best model, 1st argument sets 
c                                 random perturbation off, 2nd sets 
c                                 output of all sucessful optimizations
      call bstmod (.false.,.true.,ntry,tol,simplx,frac,conchk,iprint,
     *             iquad, kcount, bstx)

      x(1,1:nparm) = bstx(1:nparm)

      oprt = .false.

      do i = 1, nunc

         call opnimc (invprb,ntry,tol,simplx,frac,conchk,iprint,iquad,
     *                kcount)

         call bstmod (.true.,.false.,ntry,tol,simplx,frac,conchk,iprint,
     *                iquad, kcount, bstx)

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

1020  format (/,'Minimization FAILED, ifault = ',i3,', icount = ',i3,/)
1030  format ('Final coordinates: ',5(g12.6,1x))
1050  format (/,'Number of function evaluations: ',i5,', igood = ',i3,/)
1080  format ('Initial normalized coordinates: ',5(g12.6,1x))
1085  format ('Initial coordinates: ',5(g12.6,1x))
1100  format (i3,' Successes in ',i4,' tries.')
1120  format ('ntry = ',i3,', igood = ',i3,', icount = ',i3,
     *        ', obj = ',g12.6,', ibest = ',i3,', bestobj = ',g12.6)

      end 

      subroutine bstmod (randm,n6out,ntry,tol,simplx,frac,conchk,iprint,
     *                   iquad, kcount, bstx)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, n, conchk, kcount, icount, ifault, iprint,
     *        iquad, j, k, igood, ntry, ibest, id, jbest

      logical readyn, bad, randm, n6out

      double precision var(l2+k5), objf, x(l2+k5), sx(l2+k5), x0(l2+k5),
     *                 tol, step(l2+k5), simplx, bstobj, frac, bstx(*), 
     *                 bstvar(l2+k5), plow(l2+k5), pdelta(l2+k5), ssp,
     *                 bay, bstbx(l2+k5), bstbay

      external readyn, mcobj2

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------- 
c                                 read experimental data and inversion candidates
c                                 if randm, then experimental data is perturbed within
c                                 its uncertainty.
      call mcxpt (randm)
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
      bstbay = 1d99
      jbest = 0

      if (n6out) then 
c                                 output file
         call mertxt (tfname,prject,'.out',0)
         open (n6,file=tfname)

         write (n6,*) 'tol/frac/simplx',tol, frac, simplx
         write (n6,'(80(''-''))')

      end if 

      do i = 1, ntry
c                                 unscale sx
         do j = 1, n
            x(j) = plow(j) + sx(j)*pdelta(j)
         end do 

         x0 = x
c                                 initialize icount in case of failure
         icount = 0

         call minim (x, step, n, objf, kcount, iprint, tol, 
     *               conchk, iquad, simplx, var, mcobj2, icount, 
     *               ifault, oktol)

         call mcobj2 (x,objf,bad)

         if (ifault.gt.2.or.(ifault.gt.0.and.objf.gt.oktol)) then

            write (*,1020) ifault, icount, igood, ntry

         else

            igood = igood + 1
c                                 compute ss of parameter deviations
            ssp = 0d0

            do j = 1, n
               ssp = ssp + ((x(j) - plow(j))/pdelta(j))**2
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
               bstvar = var
            end if

            write (*,1120) i, igood, icount, objf, bstobj, ibest
            write (*,1130) ssp, bay, bstbay, jbest
            write (*,1080) sx(1:n)
            write (*,1085) x0(1:n)
            write (*,1030) x(1:n)

            if (n6out) then

               write (n6,1120) i, igood, icount, objf, bstobj, ibest
               write (n6,1130) ssp, bay, bstbay, jbest
               write (n6,1085) x0(1:n)
               write (n6,1030) x(1:n)
               write (n6,'(/,a,i3,a,/)') 'Scores for try = ',i,
     *                                   ' follow:'
               do id = 1, mxpt
                  write (n6,'(i3,1x,a,g12.6,1x)') 
     *                   id, xptnam(id)//' score =',scores(id)
               end do

               write (n6,'(/80(''-''))')

            end if

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

      call mcobj2 (x,objf,bad)

      if (n6out) close (n6)

1020  format (/,'Minimization FAILED, ifault = ',i3,', icount = ',i4,
     *          ', ntry = ',i4,', igood = ',i4,/)
1030  format ('Final coordinates: ',20(g12.6,1x))
1050  format (/,'Number of function evaluations: ',i5,', igood = ',i3,/)
1080  format ('Initial normalized coordinates: ',20(g12.6,1x))
1085  format ('Initial coordinates: ',20(g12.6,1x))
1120  format (/,'Try ',i3,', successes so far = ',i3,/,
     *          'Objective function evaluations this try = ',i4,/,
     *          'Last objective function value this try OBJ = ',g12.6,/,
     *          'Best OBJ so far = ',g12.6,
     *          ' obtained on try ',i3,/)
1130  format (/,'Scaled parameter SSP = ',g12.6,/,
     *          'Bayes score SSP * OBJF = ',g12.6,/,
     *          'Best Bayes score so far = ',g12.6,
     *          ' obtained on try ',i3,/)

      end 

      subroutine mccomp
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for MC thermobarometry
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, ier

      character c

      double precision total

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision atwt
      common/ cst45 /atwt(k0)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname
c-----------------------------------------------------------------------
c                                 look for phase compositions in
c                                 my_project.imc
      call mertxt (tfname,prject,'.imc',0)
      open (n8,file=tfname,status='old',iostat=ier)

      if (ier.ne.0) call errdbg 
     *   ('can''t open assemblage composition file: '//tfname)

      if (iopt(2).eq.1) then
         write (*,*) 'resetting composition_phase option to mol'
         iopt(2) = 0d0
      end if

      mphase = 0

      do

         read (n8,'(a)',iostat=ier) tname
c                                 presumed EOF
         if (ier.ne.0) exit
c                                 filter out comments
         read (tname,'(a)') c
         if (c.eq.'|'.or.c.eq.' ') cycle
c                                 got a live one
         mphase = mphase + 1
c                                 check name
         call matchj (tname,pids(mphase))

         if (pids(mphase).eq.0) 
     *      call errdbg ('no such entity as: '//tname)
c                                 all clear
         read (n8,*,iostat=ier) 
     *                        pmode(mphase), (pblk(mphase,j),j=1,kbulk),
     *                        emode(mphase), (eblk(mphase,j),j=1,kbulk)

         if (ier.ne.0) call errdbg ('invalid data format for: '//tname)
c                                 convert to molar if mass units
         if (iwt.eq.1) then 

            do j = 1, icomp
               pblk(mphase,j) = pblk(mphase,j)/atwt(j)
               eblk(mphase,j) = eblk(mphase,j)/atwt(j)
            end do 

         end if

      end do

      if (mphase.lt.2) call errdbg ('input must specify > 1 phase')
c                                 normalize to the icomp (>= icp) components
      do i = 1, mphase

         total = 0d0

         do j = 1, icomp

            total = total + pblk(i,j)

         end do

         do j = 1, icomp

            pblk(i,j) = pblk(i,j) / total
            eblk(i,j) = eblk(i,j) / total

         end do

      end do

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

      subroutine mcxpt (randm)
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for MC inversion of exptal
c data.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok, bad, randm, done

      integer i, j, k, nph, ier, id, ids

      double precision err, pertrb, tot

      character key*22, val*3, nval1*12, nval2*12, nval3*12,
     *          strg*40, strg1*40, char*1

      external pertrb

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      integer jend
      common/ cxt23 /jend(h9,m14+2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character cname*5
      common/ csta4  /cname(k5)

      double precision atwt
      common/ cst45 /atwt(k0)

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m1,h9),
     *               jsub(m2,m1,h9)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c-----------------------------------------------------------------------
      if (iopt(2).eq.1) then
         write (*,*) 'resetting composition_phase option to mol'
         iopt(2) = 0d0
      end if
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
               else 
                  call errdbg ('invalid parameter name for '//tname//
     *                         ': '//val)
               end if
c                                 read parameter range
               call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

               if (key.ne.'range') call errdbg ('expecting '//
     *                    'range tag for '//tname//' found '//key)
               read (strg1,*) cprng(mccpd,mcpct(mccpd),1),
     *                        cprng(mccpd,mcpct(mccpd),2)

            end do

         else if (id.gt.0) then

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
               do j = 1, jterm(id)

                  if (rko(j,id).ne.iord) cycle

                  ok = .true.

                  do i = 1, iord
c                                 now check the endmembers match, assume same
c                                 ordering as in the solution model
                     if (isub(1,i).ne.jend(id,2+jsub(i,j,id))) then
c                                 if .not.ok, no match on jth term
                        ok = .false.
                        exit

                     end if

                  end do

                  if (ok) exit

               end do
              
               if (.not.ok) then
c                                term does not exist in solution model
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
c
                  read (strg1,*) 
     *           mccoid(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol))),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),1),
     *           sprng(mcsol,mctrm(mcsol),mccoef(mcsol,mctrm(mcsol)),2)

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
      nparm = 0 
c                                 first for endmembers
      do i = 1, mccpd
c                                 not so bad
         nparm = nparm + mcpct(i)
      end do
c                                 now solutions
      do i = 1, mcsol
c                                 for each term
         do j = 1, mctrm(i) 
            nparm = nparm +  mccoef(i,j)
         end do
      end do

      if (nparm.eq.0) call errdbg ('no free parameters! no free lunch!')
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

         if (key.ne.'begin_exp') then
            call errdbg ('invalid data, last read '//key)
         end if
c                                an experiment, only increment
c                                counter if previous result 
c                                was ok
         if (.not.bad) mxpt = mxpt + 1

         if (mxpt.gt.l11) call errdbg ('too many expts, increase l11')
c
         bad = .false.
c                                initialize bulk
         xptblk(mxpt,1:icp) = 0
c                                initialize solution counters
         msolct(mxpt,1:isoct) = 0
c                                read name
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         xptnam(mxpt) = key
c                               read expt p,t
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         read (key,*) xptpt(mxpt,1)
         read (val,*) err

         if (randm) xptpt(mxpt,1) = pertrb (xptpt(mxpt,1),err)

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         read (key,*) xptpt(mxpt,2)
         read (val,*) err

         if (randm) xptpt(mxpt,2) = pertrb (xptpt(mxpt,2),err)
c                               read bulk composition
         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.ne.'begin_bulk') then
            call errdbg ('invalid data, last read '//key)
         end if

         do

            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

            if (key.eq.'end_bulk') exit

            ok = .false.

            do i = 1, icp
               if (key.eq.cname(i)) then
                  ok = .true.
                  exit
               end if
            end do

            if (.not.ok) bad = .true.

            read (strg1,*) xptblk(mxpt,i), err
            if (randm) xptblk(mxpt,i) = pertrb (xptblk(mxpt,i),err)

         end do

         nph = 0
         ok = .true.

         do
c                                 now read phase data
            call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

            if (key.eq.'end_exp') then
               ok = .false.
               exit
            end if
c                                 phase name
            read (key,'(a)') tname

            nph = nph + 1
c                                 check name
            call matchj (tname,ids)

            xptids(mxpt,nph) = ids

            if (ids.eq.0) call errdbg ('no such entity as: '//tname)
c                                 all clear,
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

            if (key.ne.'begin_comp') 
     *         call errdbg ('invalid data, last read '//key)

            do

               call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)

               if (key.eq.'end_comp') exit

               ok = .false.

               do i = 1, icp
                  if (key.eq.cname(i)) then
                     ok = .true.
                     exit
                  end if
               end do

               if (.not.ok) bad = .true.

               read (strg1,*) xptc(cxpt+i), err

               if (randm) xptc(cxpt+i) = pertrb (xptc(cxpt+i),err)

            end do
c                                 normalize bulk
            tot = 0d0

            do i = 1, icp
               tot = tot + xptc(cxpt+i)
            end do

            do i = 1, icp
               xptc(cxpt+i) = xptc(cxpt+i)/tot
            end do
c                                 pointer to the composition of phase nph in expt mexpt 
            xptptr(mxpt,nph) = cxpt

            xptnph(mxpt) = nph
c                                 increment composition pointer
            if (.not.bad) cxpt = cxpt + icp

         end do

         if (bad) then 
            call mertxt (tfname,prject,'.dat',0)
            write (*,1000) xptnam(mxpt), tfname
         end if
c                                 next experiment
      end do

      close (n8)

1000  format (/,'warning ver502** observation: ',a,/,'has been rejecte',
     *          'd because it includes a component not specified in: ',
     *          a,//,80('-'))

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
     *        'term will be rejected',/)

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
c set bulk composition for MC thermobarometry
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
      cblk(1:kbulk) = 0d0
      ctotal = 0d0 

      do i = 1, mphase

         k = ipot + i

         if (i.lt.mphase) then 
            ctotal = ctotal + x(k)
         else
            x(k) = 1d0 - ctotal
         end if

         do j = 1, kbulk
            cblk(j) = cblk(j) + x(k) * pblk(i,j)
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
         b(i) = cblk(i) / ctotal
      end do

      end

      subroutine mcstb2 (id)
c-----------------------------------------------------------------------
c set bulk composition for MC experiment inversion
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer is
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1),is(k1+k5)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c-----------------------------------------------------------------------
      cblk(1:kbulk) = 0d0
      ctotal = 0d0

      do i = 1, hcp
         cblk(i) = xptblk(id,i)
      end do
c                                 get total moles to compute mole fractions 
      do i = 1, hcp
         ctotal = ctotal + cblk(i)
      end do

      do i = 1, hcp 
         b(i) = cblk(i) / ctotal
      end do
c                                 set p-t
      v(1) = xptpt(id,1)
      v(2) = xptpt(id,2)

      end

      subroutine mcobj1 (x,obj,bad)
c-----------------------------------------------------------------------
c a subprogram to evaluate objective function for MC thermobarometry
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad, ok, imout(k5), imin(k5)

      integer i, j, l, kct(k5), ksol(k5,k5)

      double precision x(*), obj, total, mpred

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

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)
c-----------------------------------------------------------------------
c                                 convert the scaled potential variables back
c                                 to the normal variables:
      do i = 1, ipot
         v(jv(i)) = vmin(jv(i)) + x(i) * (vmax(jv(i)) - vmin(jv(i)))
      end do
c                                 set the bulk composition
      call mcsetb (x)
c                                 do the optimization and via getloc compute system
c                                 and derivative properties, these computations are costly 
c                                 and can be streamlined for specific applications.
      call meemum (bad)
c                                 compute the objective function
      kct(1:mphase) = 0

      ok = .false. 

      imout(1:ntot) = .true.
      imin(1:mphase) = .false.

      do i = 1, ntot

         do j = 1, mphase

            if (kkp(i).eq.pids(j)) then

               ok = .true.
               imout(i) = .false.

               kct(j) = kct(j) + 1
               ksol(j,kct(j)) = i

            end if

         end do

      end do

      if (.not.ok) then
c                                 none of the target phases found
c                                 set obj to ridiculous value
         obj = 1d99

      else

         obj = 0d0
         mpred = 0
c                                 first extraneous phases:
         do i = 1, ntot

            if (imout(i)) then 
c                                 residual is wextra * mass_fraction^2
               obj = obj + 
     *               wextra * (props(17,i)*props(16,i)/psys(17))**2

            end if

         end do
c                                 potential target phases:
         do j = 1, mphase

            if (kct(j).eq.1) then

               mpred = mpred + 1d0
c                                 found phase and no ambiguity
               if (kkp(ksol(j,1)).gt.0) then
c                                 it's a solution compute and add residual
                  total = 0d0
c                                 normalization
                  do l = 1, icomp
                     total = total + pcomp(l,ksol(j,1))
                  end do
c                                 residual
                  do l = 1, icomp
                     obj = obj +
     *               wcomp * (pcomp(l,ksol(j,1))/total - 
     *                       pblk(j,l))**2
                  end do

               end if

            else if (kct(j).gt.1) then 

               mpred = mpred + 1d0
c                                 found phase and ambiguity
               call errdbg ('need more')

            end if

         end do
c                                 missing phase residual
         obj = obj + wmiss * (1d0 - mpred/mphase)**2

      end if

      end

      subroutine minim (p,step,nop,func,max,iprint,stopcr,nloop,iquad,
     *                  simp,var,functn,neval,ifault,oktol)
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
c      the dimension of bmat should be at least nop*(nop+1)/2.
c
c       latest revision - 11 august 1991
c----------------------------------------------------------------------
c                                 original code used implicit typing
      implicit none

      external functn

      logical bad

      integer nop, lout, iprint, ifault, nloop, nap, loop, iflag, i, 
     *        irow, j, np1, imax, imin, iquad, neval, max, i1, i2, j1, 
     *        k, l, ii, ij, nullty, irank, jj, ijk

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
      if(iprint.gt.0) write(lout,1000) iprint
 1000 format(' progress report every',i4,' function evaluations'/,
     1  ' eval.  func.',15x,'parameter values')
c
c     check input arguments
c
      ifault=0
      if(nop.le.0) ifault=3
      if(nloop.le.0) ifault=4
      if(ifault.ne.0) return
c
c     set nap = no. of parameters to be varied, i.e. with step.ne.0
c
      nap=0
      loop=0
      iflag=0
      do 10 i=1,nop
        if(step(i).ne.zero) nap=nap+1
   10 continue
c
c     if nap = 0 evaluate function at the starting point and return
c
      if(nap.gt.0) go to 30
      call functn(p,func,bad)
      if (bad) ifault = 99
      return
c
c     set up the initial simplex
c
   30 do 40 i=1,nop
   40 g(1,i)=p(i)
      irow=2
      do 60 i=1,nop
        if(step(i).eq.zero) go to 60
        do 50 j=1,nop
   50   g(irow,j)=p(j)
        g(irow,i)=p(i)+step(i)
        irow=irow+1
   60 continue
      np1=nap+1
      neval=0
      do 90 i=1,np1
        do 70 j=1,nop
   70   p(j)=g(i,j)
        call functn(p,h(i),bad)

        if (bad) then 
           ifault = 99
           return
        end if

        neval=neval+1
        if(iprint.le.0) go to 90
        write(lout,1010) neval,h(i),(p(j),j=1,nop)
 1010   format(/i4, 2x, g12.5, 2x, 5g12.5, 3(/20x, 5g12.5))
   90 continue
c
c     start of main cycle.
c
c     find max. & min. values for current simplex (hmax & hmin).
c
  100 loop=loop+1
      imax=1
      imin=1
      hmax=h(1)
      hmin=h(1)
      do 120 i=2,np1
        if(h(i).le.hmax) go to 110
        imax=i
        hmax=h(i)
        go to 120
  110   if(h(i).ge.hmin) go to 120
        imin=i
        hmin=h(i)
  120 continue
c
c     find the centroid of the vertices other than p(imax)
c
      do 130 i=1,nop
  130 pbar(i)=zero
      do 150 i=1,np1
        if(i.eq.imax) go to 150
        do 140 j=1,nop
  140   pbar(j)=pbar(j)+g(i,j)
  150 continue
      do 160 j=1,nop
      fnap = nap
  160 pbar(j)=pbar(j)/fnap
c
c     reflect maximum through pbar to pstar,
c     hstar = function value at pstar.
c
      do 170 i=1,nop
  170 pstar(i)=a*(pbar(i)-g(imax,i))+pbar(i)
      call functn(pstar,hstar,bad)
      neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

      if(iprint.le.0) go to 180
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,hstar,
     1  (pstar(j),j=1,nop)
c
c     if hstar < hmin, reflect pbar through pstar,
c     hstst = function value at pstst.
c
  180 if(hstar.ge.hmin) go to 220
      do 190 i=1,nop
  190 pstst(i)=c*(pstar(i)-pbar(i))+pbar(i)
      call functn(pstst,hstst,bad)
      neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

      if(iprint.le.0) go to 200
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,hstst,
     1  (pstst(j),j=1,nop)
c
c     if hstst < hmin replace current maximum point by pstst and
c     hmax by hstst, then test for convergence.
c
  200 if(hstst.ge.hmin) go to 320
      do 210 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstst(i)
  210 continue
      h(imax)=hstst
      go to 340
c
c     hstar is not < hmin.
c     test whether it is < function value at some point other than
c     p(imax).   if it is replace p(imax) by pstar & hmax by hstar.
c
  220 do 230 i=1,np1
        if(i.eq.imax) go to 230
        if(hstar.lt.h(i)) go to 320
  230 continue
c
c     hstar > all function values except possibly hmax.
c     if hstar <= hmax, replace p(imax) by pstar & hmax by hstar.
c
      if(hstar.gt.hmax) go to 260
      do 250 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstar(i)
  250 continue
      hmax=hstar
      h(imax)=hstar
c
c     contracted step to the point pstst,
c     hstst = function value at pstst.
c
  260 do 270 i=1,nop
  270 pstst(i)=b*g(imax,i) + (1.d0-b)*pbar(i)
      call functn(pstst,hstst,bad)
      neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

      if(iprint.le.0) go to 280
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,hstst,
     1  (pstst(j),j=1,nop)
c
c     if hstst < hmax replace p(imax) by pstst & hmax by hstst.
c
  280 if(hstst.gt.hmax) go to 300
      do 290 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstst(i)
  290 continue
      h(imax)=hstst
      go to 340
c
c     hstst > hmax.
c     shrink the simplex by replacing each point, other than the current
c     minimum, by a point mid-way between its current position and the
c     minimum.
c
  300 do 315 i=1,np1
        if(i.eq.imin) go to 315
        do 310 j=1,nop
          if(step(j).ne.zero) g(i,j)=(g(i,j)+g(imin,j))*half
          p(j)=g(i,j)
  310   continue
        call functn(p,h(i),bad)
        neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

        if(iprint.le.0) go to 315
        if(mod(neval,iprint).eq.0) write(lout,1010) neval,h(i),
     1              (p(j),j=1,nop)
  315 continue
      go to 340
c
c     replace maximum point by pstar & h(imax) by hstar.
c
  320 do 330 i=1,nop
        if(step(i).ne.zero) g(imax,i)=pstar(i)
  330 continue
      h(imax)=hstar
c
c     if loop = nloop test for convergence, otherwise repeat main cycle.
c
  340 if(loop.lt.nloop) go to 100
c
c     calculate mean & standard deviation of function values for the
c     current simplex.
c
      hstd=zero
      hmean=zero
      do 350 i=1,np1
  350 hmean=hmean+h(i)
      fnp1 = np1
      hmean=hmean/fnp1
      do 360 i=1,np1
  360 hstd=hstd+(h(i)-hmean)**2
      hstd=sqrt(hstd/float(np1))
c
c     if the rms > stopcr, set iflag & loop to zero and go to the
c     start of the main cycle again.
c
      if(hstd.le.stopcr.or.neval.gt.max) go to 410
      iflag=0
      loop=0
      go to 100
c
c     find the centroid of the current simplex and the function value there.
c
  410 do 380 i=1,nop
        if(step(i).eq.zero) go to 380
        p(i)=zero
        do 370 j=1,np1
  370   p(i)=p(i)+g(j,i)
        fnp1 = np1
        p(i)=p(i)/fnp1
  380 continue

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

      if(iprint.le.0) go to 390
      if(mod(neval,iprint).eq.0) write(lout,1010) neval,func,
     1  (p(j),j=1,nop)
c
c     test whether the no. of function values allowed, max, has been
c     overrun; if so, exit with ifault = 1.
c
  390 if(neval.le.max) go to 420
      ifault=1
      if(iprint.lt.0) return
      write(lout,1020) max
 1020 format(' no. of function evaluations exceeds',i5)
      write(lout,1030) hstd
 1030 format(' rms of function values of last simplex =',g14.6)
      write(lout,1040)(p(i),i=1,nop)
 1040 format(' centroid of last simplex =',4(/1x,6g13.5))
      write(lout,1050) func
 1050 format(' function value at centroid =',g14.6)
      return
c
c     convergence criterion satisfied.
c     if iflag = 0, set iflag & save hmean.
c     if iflag = 1 & change in hmean <= stopcr then search is complete.
c
  420 if(iprint.lt.0) go to 430
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
      if(iprint.lt.0) go to 460
      write(lout,1070) neval
 1070 format(//' minimum found after',i5,' function evaluations')
      write(lout,1080)(p(i),i=1,nop)
 1080 format(' minimum at',4(/1x,6g13.6))
      write(lout,1090) func
 1090 format(' function value at minimum =',g14.6)
  460 if(iquad.le.0) return
c-------------------------------------------------------------------
c
c     quadratic surface fitting
c
      if(iprint.ge.0) write(lout,1110)
 1110 format(/' quadratic surface fitting about supposed minimum'/)
c
c     expand the final simplex, if necessary, to overcome rounding
c     errors.
c
      neval=0

      do 490 i=1,np1
  470   test=abs(h(i)-func)

        if(test.ge.simp) go to 490

        if (func.gt.oktol) then
           write (*,'(a,1x,g12.6)') 'Aborting, bad objf: ',func
           ifault = 2
           return
        end if

        do 480 j=1,nop

           if(step(j).ne.zero) g(i,j)=(g(i,j)-p(j))+g(i,j)

           pstst(j)=g(i,j)

           if (isnan(g(i,j))) then

              write (*,'(a,1x,g12.6)') 'Aborting, bad coord: ',g(i,j)
              ifault = 2
              return

           end if

  480   continue

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

        neval=neval+1
c                                 quick fix: this segment goes into an 
c                                 infinite loop if g(i,j) doesn't change
        if (neval.gt.max) then
           write (*,*) 'Aborting, infinite loop during surface fitting'
           ifault = 2
           return
        end if 

        go to 470

  490 continue
c
c     function values are calculated at an additional nap points.
c
      do 510 i=1,nap
        i1=i+1
        do 500 j=1,nop
  500   pstar(j)=(g(1,j)+g(i1,j))*half
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

        neval=neval+1

  510 continue
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

          neval=neval+1
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
c
c     the vector of estimated first derivatives is calculated and
c     stored in aval.
c
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
  600 if(iprint.ge.0) write(lout,1120)
 1120 format(/' matrix of estimated second derivatives not +ve defn.'/
     1  ' minimum probably not found'/)
      ifault=2
      return
c
c     bmat*a/2 is calculated and stored in h.
c
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
c
c     find the position, pmin, & value, ymin, of the minimum of the
c     quadratic.
c
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
      if(iprint.lt.0) go to 682
      write(lout,1130) ymin,(pmin(i),i=1,nop)
 1130 format(' minimum of quadratic surface =',g14.6,' at',
     1  4(/1x,6g13.5))
      write(lout,1150)
 1150 format(' if this differs by much from the minimum estimated',
     1  1x,'from the minimization,'/
     2  ' the minimum may be false &/or the information matrix may be',
     3  1x,'inaccurate'/)
c
c     calculate true function value at the minimum of the quadratic.
c
  682 neval = neval + 1
      call functn(pmin, hstar,bad)

        if (bad) then 
           ifault = 99
           return
        end if
c
c     if hstar < func, replace search minimum with quadratic minimum.
c
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
c
c     the diagonal elements of vc are copied into var.
c
      j=0
      do 770 i=1,nop
        j=j+i
        var(i)=vc(j)
  770    continue
      if(iprint.lt.0) return
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
  860 write(lout,1210) neval
 1210 format(/' a further',i4,' function evaluations have been used'/)
      return
c
c     pseudo-subroutine to print vc if ijk = 1 or 2, or
c     bmat if ijk = 3.
c
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
c compositions for MC
c-----------------------------------------------------------------------
      include 'perplex_parameters.h'

      integer id, kd, j, l

      double precision total

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c-----------------------------------------------------------------------
      total = 0d0
      score = 0d0
c                                 normalization
      do l = 1, icomp
         total = total + pcomp(l,kd)
      end do
c                                 residual
      do l = 1, icomp
         score = score + (pcomp(l,kd)/total - 
     *                    xptc(xptptr(id,j)+l))**2
      end do

      end

      subroutine x2ther (x)
c-----------------------------------------------------------------------
c a to map the MC search coordinates into thermodynamic arrays
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, n, ids, jt

      double precision x(*)
c-----------------------------------------------------------------------
      n = 0
c                                 compounds
      do i = 1, mccpd
c                                 for each coefficient
         do j = 1, mcpct(i)
            n = n + 1
            mdqf(make(mcid(i)),mcpid(i,j)) = x(n)
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
               wgl(mccoid(ids,j,k),jt,ids) = x(n)
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

      logical bad, ok, imout(k5), imin(k5), used(k5)

      integer id, jd, ids, i, j, kct(k5), ksol(k5,k5), ibest

      double precision x(*), lobj, mpred, obj, score, best, res

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

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)
c-----------------------------------------------------------------------
      do i = 1, nparm
         if (x(i).lt.-1d99.or.x(i).gt.1d99) then
            write (*,*) 'woo',x(1:nparm)
            bad = .true.
            return
         end if
      end do
c                                map the search coordinates to thermodynamic
c                                parameters
      call x2ther (x)
c                                 compositional residual weight
      wcomp  = 1d0
c                                 extra phase amount residual weight
      wextra = 1d1
c                                 missing phase residual weight
      wmiss  = 1d1

      obj = 0d0
c                                 loop through observations
      do id = 1, mxpt
c
         scores(id) = 1d99
c                                 set p, t, bulk
         call mcstb2 (id)
c                                 do the optimization and via getloc compute system
c                                 and derivative properties, these computations are costly 
c                                 and can be streamlined for specific applications.
         call meemum (bad)

         if (bad) then 
            write (*,*) 'oink'
         end if

c        call calpr0 (6)
c                                 compute the observation objective function
         kct(1:xptnph(id)) = 0

         ok = .false. 

         imout(1:ntot) = .true.
         used = .false.
         imin(1:mphase) = .false.
         kct = 0
         ksol = 0

         do i = 1, ntot

            do j = 1, xptnph(id)

               if (kkp(i).eq.xptids(id,j)) then

                  ok = .true.
                  imout(i) = .false.

                  kct(j) = kct(j) + 1
                  ksol(j,kct(j)) = i

               end if

            end do

         end do

         lobj = 0d0
         mpred = 0
c                                 first extraneous phases:
         do i = 1, ntot

            if (imout(i)) then 
c                                 residual is wextra * mass_fraction^2
               lobj = lobj + 
     *                wextra * (props(17,i)*props(16,i)/psys(17))**2

            end if

         end do
c                                 potential target phases:
         do j = 1, xptnph(id)

            if (kct(j).eq.0) cycle

            jd = ksol(j,1)
            ids = kkp(jd)

            if (kct(j).eq.1.and.ids.lt.0) then
c                                 a compound, just count
               mpred = mpred + 1d0

            else if (kct(j).eq.1.and.
     *               msolct(id,ids).eq.1) then
c                                 found solution and no ambiguity
               mpred = mpred + 1d0
c                                 compute and add residual
               lobj = lobj + wcomp * score (jd,id,j)

            else if (kct(j).gt.1) then

               best = 1d99
               ok = .false.

               do i = 1, kct(j)

                  jd = ksol(j,i)

                  if (used(jd)) cycle

                  res = score (jd,id,j)

                  if (res.lt.best) then
                     best = res
                     ibest = i
                     ok = .true.
                  end if

               end do

               if (ok) then

                  used(ibest) = .true.
                  mpred = mpred + 1d0
                  lobj = lobj + wcomp * best

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
               lobj = lobj + wextra * 
     *                      (props(17,jd)*props(16,jd)/psys(17))**2

            end do

         end do
c                                 missing phase residual
         lobj = lobj + wmiss * (1d0 - mpred/xptnph(id))
c                                 accumulate scores
         obj = obj + lobj

         scores(id) = lobj

c        write (*,'(i3,1x,2(g12.6,1x,a))') id, lobj, xptnam(id)

c        if (oprt) write (n6,'(i3,1x,a,g12.6,1x)') 
c    *                   id, xptnam(id)\\' score =', lobj

      end do

c     write (*,'(8(g12.6,1x))') obj, x(1:nparm)

      end