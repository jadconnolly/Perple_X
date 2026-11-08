 
c MLIB - a library of subroutines for mc_fit, predominantly Nelder-Mead
c simplex-based function minimization.

      subroutine minim (p,step,nop,func,max,jprint,stopcr,nloop,
     *              iquad,simp,covar,functn,neval,meval,ifault,oktol)
c----------------------------------------------------------------------
c     a program for function minimization using the simplex method.
c     the minimum found will often be a local, not a global, minimum.
c
c     the version of minim here has been modified for mc_fit and may
c     not function identically to the original code.

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
c     nbdd    = input no of bounded variables for stink hack
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
c     stopcr  = input, stopping criterion, rms of obj function on the vertices
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
     *        k, l, ii, ij, nullty, irank, jj, meval

      double precision zero, one, two, three, half, a, b, c, func, 
     *                 fnap, fnp1, savemn, test, simp, a0, stopcr,
     *                 rmax, ymin, hmax, hmin, hstar, hstst, hstd, 
     *                 hmean, p(nop),step(nop),oktol, stink, rvar, fvar,
     *                 g(21,20),h(21),pbar(20),pstar(20),pstst(20),
     *                 aval(20),bmat(210),pmin(20),vc(210),temp(20),
     *                 covar(*)

      double precision xinv
      common/ cst67 /xinv(20)

      external stink

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
c-----------------------------------------------------------------------
c     if progress reports have been requested, print heading
c
      if(jprint.gt.0) write(lout,1000) jprint
 1000 format('progress report every',i4,' function evaluations'/,
     1       'eval.  func.',15x,'parameter values')

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

      call functn (p,func,rvar,bad)

      if (bad) then
         write (*,*) 'hard abort e'
         ifault = 99
      end if

      return

c     set up the initial simplex

30    do i=1,nop
         g(1,i)=p(i)
      end do

      irow=2

      do i=1,nop

        if(step(i).eq.zero) cycle

        do j=1,nop
           g(irow,j)=p(j)
        end do

        g(irow,i) = p(i)+step(i)

        call stink1 (g(irow,i),p(i),i,ifault)

        if (ifault.eq.99) then 
           write (*,*) 'hard abort e stinknew '
           return
        end if

        irow = irow + 1

      end do

      np1 = nap + 1
      neval = 0
      meval = 0

      do i = 1, np1

         do j = 1, nop
            p(j) = g(i,j)
         end do

         call functn (p,h(i),rvar,bad)

         if (bad) then
            write (*,*) 'hard abort d'
            ifault = 99
            return
         end if

        neval=neval+1

        if(jprint.le.0) cycle

        write(lout,1010) neval,h(i),(xinv(j),j=1,nop)

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

         pstar(i) = a*(pbar(i)-g(imax,i)) + pbar(i)

         call stink1 (pstar(i),pbar(i),i,ifault)

         if (ifault.eq.99) then
            write (*,*) 'hard abort c oldstink '
            return
          end if

      end do

      call functn (pstar,hstar,rvar,bad)

      neval=neval+1

      if (bad) then
         write (*,*) 'hard abort b oldstink '
        ifault = 99
        return
      end if

      if (jprint.gt.0.and.mod(neval,jprint).eq.0) then 
         write (lout,1010) neval,hstar, (xinv(j),j=1,nop)
      end if

c     if hstar < hmin, reflect pbar through pstar,
c     hstst = function value at pstst.

      if (hstar.ge.hmin) go to 220

      do i = 1, nop

         pstst(i)=c*(pstar(i)-pbar(i))+pbar(i)

         call stink1 (pstst(i),pbar(i),i,ifault)

         if (ifault.eq.99) then
            write (*,*) 'hard abort a'
            return
         end if

      end do

      call functn(pstst,hstst,rvar,bad)

      neval=neval+1

        if (bad) then 
           write (*,*) 'hard abort a'
           ifault = 99
           return
        end if

      if (jprint.gt.0 .and. mod(neval,jprint).eq.0) then
         write(lout,1010) neval,hstst,(xinv(j),j=1,nop)
      end if
c
c     if hstst < hmin replace current maximum point by pstst and
c     hmax by hstst, then test for convergence.
c
      if(hstst.ge.hmin) go to 320

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

      call functn(pstst,hstst,rvar,bad)

      neval=neval+1

        if (bad) then 
           ifault = 99
           return
        end if

      if (jprint.gt.0 .and.mod(neval,jprint).eq.0) then 
         write(lout,1010) neval,hstst, (xinv(j),j=1,nop)
      end if

c     if hstst < hmax replace p(imax) by pstst & hmax by hstst.

      if (hstst.gt.hmax) go to 300

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

         call functn(p,h(i),rvar,bad)

         neval = neval+1

         if (bad) then 
            ifault = 99
            return
         end if

         if (jprint.gt.0 .and. mod(neval,jprint).eq.0) then 
            write(lout,1010) neval,h(i), (xinv(j),j=1,nop)
         end if

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

      hstd=sqrt(dmax1(0d0,hstd)/np1)

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

      call functn(p,func,fvar,bad)

      neval=neval+1

      if (bad) then 
         ifault = 99
         return
      end if

      if (func.gt.oktol) then
         write (*,'(a,1x,g12.6)') 'Hard abort, bad objf: ',func
         ifault = 99
         return
      end if

      if (jprint.gt.0 .and. mod(neval,jprint).eq.0) then
         write(lout,1010) neval,func, (xinv(j),j=1,nop)
      end if

c     test whether the no. of function values allowed, max, has been
c     overrun; if so, exit with ifault = 1.

      if (neval.le.max) go to 420

      ifault = 1

c     return minimum found so far

      if (func.gt.h(imin)) then
         p(1:nop) = g(imin,1:nop)
         call trnsxy (p)
         func = h(imin)
      end if

      if (jprint.lt.0) return

      write(lout,1020) max
      write(lout,1030) hstd
      write(lout,1040) (xinv(i),i=1,nop)
      write(lout,1050) 'search limit', func

      return

c     convergence criterion satisfied.
c     if iflag = 0, set iflag & save hmean.
c     if iflag = 1 & change in hmean <= stopcr then search is complete.

420   if (jprint.gt.0) then
 
         write(lout,1060)
         write(lout,1040) (xinv(i),i=1,nop)
         write(lout,1050) 'centroid',func

      end if

      if(iflag.gt.0) go to 450

      iflag=1

440   savemn = hmean
      loop = 0

      go to 100

450   if (abs(savemn-hmean).ge.stopcr) go to 440

c     if func at last centroid was higher than minimum, use minimum

      if (func.gt.h(imin) .and. iquad.le.0) then

         if (jprint.gt.0) write (lout,1090)
         p(1:nop) = g(imin,1:nop)
         call trnsxy (p)
         func = h(imin)

      end if

      if (jprint.gt.0) then

         write(lout,1070) neval, (p(i), i = 1, nop)
         write(lout,1050) 'minimum',func

      else if (jprint.eq.0) then

         write(lout,1080) neval, func

      end if

      if (iquad.le.0) return
c-------------------------------------------------------------------

c     quadratic surface fitting

      if(jprint.ge.0) write(lout,1110)

c     expand the final simplex, if necessary, to overcome rounding
c     errors.

      do i=1,np1

  470   test=abs(h(i)-func)

        if(test.ge.simp) cycle

        if (func.gt.oktol) then
           write (*,'(a,1x,g12.6)') 'Abort during quad fit, bad objf: ',
     *                              func
           ifault = 99

           return

        end if

        do j=1,nop

           if (step(j).ne.zero) then 
              test = (g(i,j)-p(j))+g(i,j)

              call stink1 (pstst(i),pbar(i),i,ifault)

              if (ifault.eq.99) then
                 write (*,*) 'qabort 1'
                 ifault = 2
                 return
              end if

              g(i,j) = test

           end if

           pstst(j) = g(i,j)

           if (isnan(g(i,j))) then

              write (*,'(a,1x,g12.6)') 
     *              'Aborting during quad fit, bad coord: ',g(i,j)
              ifault = 2

              return

           end if

        end do

        call functn (pstst,h(i),rvar,bad)

        if (bad) then
           write (*,'(a,1x,g12.6)') 'Abort during quad, bad opt'
           ifault = 99
           return
        end if

        if (h(i).gt.oktol) then
           write (*,'(a,1x,g12.6)') 'Aborting quad, bad objf: ',h(i)
           ifault = 99
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

         call functn(pstar,aval(i),rvar,bad)

         if (bad) then 
            ifault = 99
            return
         end if

         if (aval(i).gt.oktol) then
            write (*,'(a,1x,g12.6)') 'Abort quad, bad objf: ',h(i)
            ifault = 99
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

        do j=1,i1
 
           j1=j+1

           do k=1,nop
              pstst(k)=(g(i2,k)+g(j1,k))*half
           end do

           call functn(pstst,hstst,rvar,bad)

           if (bad) then 
              ifault = 99
              return
           else if (hstst.gt.oktol) then
              write (*,'(a,1x,g12.6)') 'Abort quad, bad objf: ',h(i)
              ifault = 2
              return
           end if

           meval = meval+1
           l=i*(i-1)/2+j
           bmat(l)=two*(hstst+a0-aval(i)-aval(j))

        end do

  540 continue

      l=0

      do i=1,nap
         i1=i+1
         l=l+i
         bmat(l)=two*(h(i1)+a0-two*aval(i))
      end do

c     the vector of estimated first derivatives is calculated and
c     stored in aval.

      do i=1,nap
        i1=i+1
        aval(i)=two*aval(i)-(h(i1)+three*a0)*half
      end do
c
c     the matrix q of nelder & mead is calculated and stored in g.
c
      do i=1,nop
         pmin(i)=g(1,i)
      end do

      do i=1,nap

        i1=i+1

        do j=1,nop
           g(i1,j)=g(i1,j)-g(1,j)
         end do

      end do 

      do i=1,nap
         i1=i+1

         do j=1,nop
           g(i,j)=g(i1,j)
         end do

      end do
c
c     invert bmat
c
      call syminv(bmat,nap,bmat,temp,nullty,ifault,rmax)

      if (ifault.ne.0) then
c                                 keeping the old min, reset xinv
         call trnsxy (p)

         if(jprint.ge.0) write(lout,1120)
         ifault=2
         return

      end if

      irank=nap-nullty

c     bmat*a/2 is calculated and stored in h.

      do 650 i=1,nap
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

      do i=1,nap
         ymin = ymin+h(i)*aval(i)
      end do

      ymin=a0-ymin

      do i = 1,nop

         pstst(i) = zero

         do j = 1,nap

           pstst(i) = pstst(i)+h(j)*g(j,i)

         end do
      end do

      do i = 1,nop
         pmin(i) = pmin(i) - pstst(i)
      end do

      if (jprint.ge.0) then

         call trnsxy (pmin)

         write(lout,1130) ymin, (xinv(i),i=1,nop)

         write(lout,1150)

      end if

c     calculate true function value at the minimum of the quadratic.

      meval = meval + 1

      call functn (pmin, hstar, rvar, bad)

        if (bad) then 
           write (*,*) 'qabort 2'
           ifault = 99
           return
        end if

c     if hstar < func, replace search minimum with quadratic minimum.

      if (hstar .lt. func) then

         func = hstar
         fvar = rvar

         do i = 1, nop
            p(i) = pmin(i)
         end do

         write(lout, 1140) func

      end if
c
c     q*bmat*q'/2 is calculated & its lower triangle stored in vc
c
c                                 keeping the old min, reset xinv
      call trnsxy (p)

      do i = 1, nop

         do j = 1, nap

          h(j)=zero

          do k = 1, nap

            if(k.gt.j) go to 700
            l=j*(j-1)/2+k
            go to 710
  700       l=k*(k-1)/2+j
  710       h(j)=h(j)+bmat(l)*g(k,i)*half

          end do

       end do

        do j=i,nop

           l=j*(j-1)/2+i
           vc(l)=zero

           do k=1,nap
              vc(l)=vc(l)+h(k)*g(k,j)
           end do

        end do

      end do
c                                 covariance matrix, use fvar is 1/2 the residual variance
c                                 for wChi (=1), LSQ, and Chi (dubious if this is the right
c                                 factor)
      k = 0
      l = 1

      do i = 1, nop

         do j = 1+(i-1)*nop, l+(i-1)*nop
            k = k + 1
            covar(j) = fvar * vc(k)
         end do

         l = l + 1

      end do

c                                 compute the information matrix, inv(covariance)
      call syminv(vc,nap,bmat,temp,nullty,ifault,rmax)
c                                 compute the correlation matrix
      ii=0
      ij=0

      do i=1,nop

        ii=ii+i

        if(vc(ii).gt.zero) then
          vc(ii)=one/sqrt(vc(ii))
        else 
          vc(ii)=zero
        end if

        jj=0

        do j=1,i-1

          jj=jj+j
          ij=ij+1
          vc(ij)=vc(ij)*vc(ii)*vc(jj)
        end do

        ij=ij+1

      end do
c                                 correlation matrix:
c                                 vc stores the elements in the lower diagonal
c                                 put them in the upper diagonal of covmat
      k = 0

      do i = 1, nop
         do j = 1, i

            k = k + 1
            if (i.eq.j) cycle
            covar((j-1)*nop+i) = vc(k)

         end do
      end do

1020  format('No. of function evaluations exceeds ',i5)
1030  format('RMS of function values of last simplex =',g14.6)
1040  format('Centroid of last simplex =',4(/1x,1p6g13.5))
1050  format('Function value at ',a,' =',g14.6)
1060  format(/'Evidence of convergence')
1070  format(//'Minimum found after ',i5,' function evaluations at',
     *                                                4(/1x,1p6g13.6))
1080  format(//'Minimum after ',i5,' function evaluations =',1pg13.6)
1090  format(/'Centroid not the minimum, using lowest vertex')
1110  format(/'Quadratic surface fitting about supposed minimum.'/)
1120  format(/'Matrix of estimated second derivatives not +ve defn.'/
     *        'Minimum probably not found.'/)
1130  format('Minimum of quadratic surface =',g14.6,' at',
     *       4(/1x,6g13.5))
1140  format('True func. value at minimum of quadratic = ', g14.6/)
1150  format('If this differs by much from the minimum estimated',
     *  1x,'from the minimization,'/
     *  'the minimum may be false &/or the information matrix may be',
     *  1x,'inaccurate'/)

      end

      subroutine syminv(a,n,c,w,nullty,ifault,rmax)
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
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
