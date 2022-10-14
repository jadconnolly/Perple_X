
c a library of BLAS level 2 subroutines. Modified subroutines are named
c by appending a digit to the original name or in the case of 6 char
c names by replacing the final character by a digit (1).

c BLAS is a freely-available software package at www.netlib.org/blas/

c JADC, 5/2020.

      subroutine dscal (n, alpha, x, incx)
c----------------------------------------------------------------------
c blas routine
c     x := alpha*x
c----------------------------------------------------------------------
      implicit none

      double precision alpha, x(*)

      integer incx, n, ix
c----------------------------------------------------------------------
      if (n.gt.0) then

         if (alpha.eq.0d0) then

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = 0d0
            end do

         else if (alpha.ne.1d0) then

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = alpha*x(ix)
            end do 

         end if

      end if

      end

      double precision function ddot1 (n, x, incx, y)
c----------------------------------------------------------------------
c modified blas routine ddot:
c 1) incy dropped as argument and fixed to 1

      implicit none

      integer i, ix, iy, incx, n

      double precision x(*), y(*), sum
c----------------------------------------------------------------------
      sum = 0d0

      if (n.gt.0) then

         if ((incx.eq.1).and.(incx.gt.0)) then

            do ix = 1, 1 + (n - 1)*incx, incx
               sum = sum + x(ix)*y(ix)
            end do

         else

            iy = 1    
            if (incx.gt.0) then

               do ix = 1, 1 + (n - 1)*incx, incx
                  sum = sum + x(ix)*y(iy)
                  iy  = iy  + 1
               end do

            else

               ix = 1 - (n - 1)*incx

               do i = 1, n
                  sum = sum + x(ix)*y(iy)
                  ix  = ix  + incx
                  iy  = iy  + 1
               end do 

            end if
         end if
      end if

      ddot1 = sum

      end

      subroutine sssq (n, x, incx, scale, sumsq)
c----------------------------------------------------------------------
c blas routine sssq:

      implicit none

      integer incx, n, ix

      double precision scale, sumsq, x(*), absxi
c----------------------------------------------------------------------
      if (n.gt.0) then

         do ix = 1, 1 + (n-1)*incx, incx

            if (x(ix).ne.0d0) then

               absxi = dabs(x(ix))

               if (scale.lt.absxi) then
                  sumsq = 1d0 + sumsq*(scale/absxi)**2
                  scale = absxi
               else
                  sumsq = sumsq + (absxi/scale)**2
               end if

            end if

         end do

      end if

      end

      double precision function dnrm2 (n, x, incx)
c----------------------------------------------------------------------
      implicit none

      integer incx, n

      double precision x(*), norm, scale, ssq, snorm
c----------------------------------------------------------------------
      if (n.lt.1) then

         norm  = 0d0

      else if (n.eq.1) then

         norm  = dabs(x(1))

      else

         scale = 0d0 
         ssq   = 1d0

         call sssq (n, x, incx, scale, ssq)

         norm  = snorm (scale, ssq)

      end if

      dnrm2 = norm

      end

      integer function idamx1 (n, x)
c----------------------------------------------------------------------
c modified blas routine idamax
c drops incx argument (set to 1)
c----------------------------------------------------------------------
      implicit none

      double precision x(*), xmax

      integer i, imax, ix, n
c----------------------------------------------------------------------
      if (n.gt.0) then

         imax = 1

         if (n.gt.1) then

            xmax = dabs(x(1))
            ix   = 1

            do i = 2, n

               ix = ix + 1

               if (xmax.lt.dabs(x(ix))) then
                  xmax = dabs(x(ix))
                  imax = i
               end if

            end do

         end if

      else

         imax = 0

      end if

      idamx1 = imax

      end


      subroutine dswap (n, x, incx, y, incy)
c----------------------------------------------------------------------
c blas routine 
c----------------------------------------------------------------------
      implicit none

      integer incx, incy, n, i, ix, iy

      double precision x(*), y(*), temp 
c----------------------------------------------------------------------

      if (n.gt.0) then

         if ((incx.eq.incy).and.(incy.gt.0)) then

            do iy = 1, 1 + (n - 1)*incy, incy
               temp    = x(iy)
               x(iy) = y(iy)
               y(iy) = temp
            end do

         else

            if (incx.ge.0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if

            if (incy.gt.0) then

               do iy = 1, 1 + (n - 1)*incy, incy
                  temp    = x(ix)
                  x(ix) = y(iy)
                  y(iy) = temp
                  ix      = ix      + incx
               end do 

            else

               iy = 1 - (n - 1)*incy

               do i = 1, n
                  temp    = x(ix)
                  x(ix) = y(iy)
                  y(iy) = temp
                  iy      = iy      + incy
                  ix      = ix      + incx
               end do 

            end if
         end if
      end if

      end

      subroutine dgemv (trans, m, n, alpha, a, lda, x, beta, y)
c----------------------------------------------------------------------
c modified blas routine dgemv
c assumes incy = incx = 1
c eliminates trans = 'c'

c  dgemv  does one of the matrix-vector operations

c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,

c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.

c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*a'*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.

c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.

c  alpha  - double precision.
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.

c  a      - double precision array of dimension (lda, n).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients.
c           unchanged on exit.

c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max(1, m).
c           unchanged on exit.

c  x      - double precision array of dimension at least
c           (1 + (n - 1)*abs(incx)) when trans = 'n' or 'n'
c           and at least
c           (1 + (m - 1)*abs(incx)) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.

c  beta   - double precision.
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.

c  y      - double precision array of dimension at least
c           (1 + (m - 1)*abs(incy)) when trans = 'n' or 'n'
c           and at least
c           (1 + (n - 1)*abs(incy)) otherwise.
c           before entry with beta non-zero, the incremented array y
c           must contain the vector y. on exit, y is overwritten by the
c           updated vector y.

c----------------------------------------------------------------------
      implicit none

      double precision alpha, beta
      integer lda, m, n
      character        trans*1

      double precision a(lda, *), x(*), y(*)



      double precision temp, temp1, temp2, temp3, temp4
      integer i, info, j, jx, kx, ky, lenx, leny, m4, n4

c     test the input parameters.

      info = 0
      if     (.not.(trans.eq.'n'.or.trans.eq.'n').and.
     *         .not.(trans.eq.'t'.or.trans.eq.'t').and.
     *         .not.(trans.eq.'c'.or.trans.eq.'c')    ) then
         info = 1
      else if (m.lt.0) then
         info = 2
      else if (n.lt.0) then
         info = 3
      else if (lda.lt.max(1, m)) then
         info = 6
      end if

      if (info.ne.0) return

c     quick return if possible.

      if ((m.eq.0).or.(n.eq.0).or.
     *    ((alpha.eq.0d0).and.(beta.eq.1d0))) return

c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.

      if ((trans.eq.'n'.or.trans.eq.'n')) then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if

         kx = 1
         ky = 1

c     start the operations. in this version the inner loops are all
c     equivalent to axpy operations.

c     first form  y := beta*y.

      if (beta.ne.1d0)then

            if (beta.eq.0d0)then
               do 10, i = 1, leny
                  y(i) = 0d0
   10          continue
            else
               do 20, i = 1, leny
                  y(i) = beta*y(i)
   20          continue
            end if

      end if

      if (alpha.eq.0d0) return
      jx = kx

      if ((trans.eq.'n'.or.trans.eq.'n')) then

c        form  y := alpha*a*x + y.

c**** u n r o l l   t o   d e p t h   4 ********************************
            n4 = 4*(n/4)
            do 60, j = 1, n4, 4
               temp1 = alpha*x(jx)
               temp2 = alpha*x(jx + 1)
               temp3 = alpha*x(jx + 2)
               temp4 = alpha*x(jx + 3)
               if (temp1.ne.0d0.or.temp2.ne.0d0.or.temp3.ne.0d0.or.
     *             temp4.ne.0d0) then
                  do 50, i = 1, m
                     y(i) = ((((y(i) + temp1*a(i, j))
     *                        + temp2*a(i, j + 1))
     *                        + temp3*a(i, j + 2))
     *                        + temp4*a(i, j + 3))
   50             continue
               end if
               jx = jx + 4
   60       continue
c**** clean-up loop ****************************************************
            do 80, j = n4 + 1, n, 1
               temp = alpha*x(jx)
               if (temp.ne.0d0) then
                  do 70, i = 1, m
                     y(i) = y(i) + temp*a(i, j)
   70             continue
               end if
               jx = jx + 1
   80       continue

      else

c        form  y := alpha*a'*x + y.

c**** u n r o l l   t o   d e p t h   4 ********************************
            m4 = 4*(m/4)
            do 140, j = 1, m4, 4
               temp1 = alpha*x(jx)
               temp2 = alpha*x(jx + 1)
               temp3 = alpha*x(jx + 2)
               temp4 = alpha*x(jx + 3)
               if (temp1.ne.0d0.or.temp2.ne.0d0.or.temp3.ne.0d0.or.
     *             temp4.ne.0d0)then
                  do 130, i = 1, n
                     y(i) = ((((y(i) + temp1*a(j, i))
     *                        + temp2*a(j + 1, i))
     *                        + temp3*a(j + 2, i))
     *                        + temp4*a(j + 3, i))
  130             continue
               end if
               jx = jx + 4
  140       continue
c**** clean-up loop ****************************************************
            do 160, j = m4 + 1, m, 1
               temp = alpha*x(jx)
               if (temp.ne.0d0)then
                  do 150, i = 1, n
                     y(i) = y(i) + temp*a(j, i)
  150             continue
               end if
               jx = jx + 1
  160       continue

      end if

      end

      subroutine daxpy1 (n, alpha, x, incx, y)
c----------------------------------------------------------------------
c modified blas routine daxpy

      implicit none

      integer incx, n, i, ix, iy

      double precision x(*), y(*), alpha
c----------------------------------------------------------------------
      if (n.gt.0) then

         if (alpha.ne.0d0) then

            if (incx.eq.1.and.incx.gt.0) then

               do ix = 1, 1 + (n - 1)*incx, incx
                  y(ix) = alpha*x(ix) + y(ix)
               end do

            else
              
               iy = 1

               if (incx.gt.0) then

                  do ix = 1, 1 + (n - 1)*incx, incx
                     y(iy) = alpha*x(ix) + y(iy)
                     iy      = iy            + 1
                  end do

               else

                  ix = 1 - (n - 1)*incx

                  do i = 1, n
                     y(iy) = alpha*x(ix) + y(iy)
                     ix      = ix            + incx
                     iy      = iy            + 1
                  end do 

               end if
            end if
         end if
      end if

      end

      subroutine daxpy (n, alpha, x, incx, y,incy)
c----------------------------------------------------------------------
c blas routine daxpy

      implicit none

      integer incx, incy, n, i, ix, iy

      double precision x(*), y(*), alpha
c----------------------------------------------------------------------
      if (n.gt.0) then

         if (alpha.ne.0d0) then

            if (incx.eq.incy.and.incx.gt.0) then

               do ix = 1, 1 + (n - 1)*incx, incx
                  y(ix) = alpha*x(ix) + y(ix)
               end do

            else

               if (incy.ge.0) then
                  iy = 1
               else
                  iy = 1 - (n - 1)*incy
               end if

               if (incx.gt.0) then

                  do ix = 1, 1 + (n - 1)*incx, incx
                     y(iy) = alpha*x(ix) + y(iy)
                     iy      = iy            + incy
                  end do

               else

                  ix = 1 - (n - 1)*incx

                  do i = 1, n
                     y(iy) = alpha*x(ix) + y(iy)
                     ix      = ix            + incx
                     iy      = iy            + incy
                  end do 

               end if
            end if
         end if
      end if

      end

      double precision function adivb (a,b,fail)
c----------------------------------------------------------------------
c blas routine sdiv, note coding error! 

c orginal code was using 

c flmin = 1/flmax ~= wmach(3) = epsmch
c flmax = wmach(7) = huge()
c----------------------------------------------------------------------
      implicit none

      logical fail

      double precision a, b, absb, absa, div

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      absa = dabs(a)
      absb = dabs(b)

      if (absa.le.wmach(3)) then
c                                 a is < machine eps
         div = 0d0

         if (absb.le.wmach(3)) then
c                                 b is also < machine eps
            fail = .true.
         else
            fail = .false.
         end if

      else
c                                 a is > eps
         if (absb.le.wmach(3)) then
c                                 b is < eps, div = sign(a)*huge
            div  =  dsign (wmach(8),a)
c                                 this doesn't make sense, but it's the
c                                 original code. jadc, jul 4, 2019.
c                                 changed to sqrt(huge) [wmach(8)] because
c                                 -huge causes a fpe if a < 0. jadc, 6/24/2020.
            fail = .true.

         else
c                                 b is > eps
            if (absb.ge.1d0) then
c                                 |b| > 1
               fail = .false.

               if (absa.ge.absb*wmach(3)) then
                  div = a/b
               else
                  div = 0d0
               end if

            else
c                                 |b| < 1
               if (absa.le.absb*wmach(7)) then

                  fail = .false.
                  div  =  a/b

               else

                  fail = .true.
                  div  = wmach(7)

                  if (((a.lt.0d0).and.(b.gt.0d0)).or.
     *                ((a.gt.0d0).and.(b.lt.0d0))) div = -div

               end if

            end if

         end if

      end if

      adivb = div

      end

      subroutine dcopy (n, x, incx, y, incy)
c----------------------------------------------------------------------
c blas routine dcopy
c----------------------------------------------------------------------
      implicit none

      integer i, ix, iy, incx, incy, n

      double precision x(*), y(*)
c----------------------------------------------------------------------
      if (n.gt.0) then
         if (incx.eq.incy.and.incy.gt.0) then
            do iy = 1, 1 + (n - 1)*incy, incy
               y(iy) = x(iy)
            end do
         else
            if (incx.ge.0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incy.gt.0) then
               do iy = 1, 1 + (n - 1)*incy, incy
                  y(iy) = x(ix)
                  ix      = ix      + incx
               end do 
            else
               iy = 1 - (n - 1)*incy
               do i = 1, n
                  y(iy) = x(ix)
                  iy      = iy      + incy
                  ix      = ix      + incx
               end do 
            end if
         end if
      end if

      end

      subroutine dtrsv (uplo, trans, diag, n, a, lda, x, incx)
c----------------------------------------------------------------------
c blas routine dtrsv

c  dtrsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.

c  parameters
c  ==========

c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:

c              uplo = 'u' or 'u'   a is an upper triangular matrix.

c              uplo = 'l' or 'l'   a is a lower triangular matrix.

c           unchanged on exit.

c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   a'*x = b.

c           unchanged on exit.

c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:

c              diag = 'u' or 'u'   a is assumed to be unit triangular.

c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.

c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.

c  a      - double precision array of dimension (lda, n).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.

c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max(1, n).
c           unchanged on exit.

c  x      - double precision array of dimension at least
c           (1 + (n - 1)*abs(incx)).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.

c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c----------------------------------------------------------------------
      implicit none 
      integer incx, lda, n
      character*1 diag, trans, uplo
      double precision a(lda, *), x(*)
      double precision temp
      integer i, info, ix, j, jx, kx
      logical            nounit
c----------------------------------------------------------------------

      info = 0
      if     (.not.(uplo .eq.'u'.or.uplo .eq.'u').and.
     *         .not.(uplo .eq.'l'.or.uplo .eq.'l')    ) then
         info = 1
      else if (.not.(trans.eq.'n'.or.trans.eq.'n').and.
     *         .not.(trans.eq.'t'.or.trans.eq.'t').and.
     *         .not.(trans.eq.'c'.or.trans.eq.'c')    ) then
         info = 2
      else if (.not.(diag .eq.'u'.or.diag .eq.'u').and.
     *         .not.(diag .eq.'n'.or.diag .eq.'n')    ) then
         info = 3
      else if (n.lt.0) then
         info = 4
      else if (lda.lt.max(1, n)) then
         info = 6
      else if (incx.eq.0) then
         info = 8
      end if
      if (info.ne.0) return

c     quick return if possible.

      if (n.eq.0) return
c
      nounit = (diag.eq.'n'.or.diag.eq.'n')

c     set up the start point in x if the increment is not unity. this
c     will be  (n - 1)*incx  too small for descending loops.

      if (incx.le.0) then
         kx = 1 - (n - 1)*incx
      else if (incx.ne.1) then
         kx = 1
      end if

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.

      if ((trans.eq.'n'.or.trans.eq.'n')) then

c        form  x := inv(a)*x.

         if ((uplo.eq.'u'.or.uplo.eq.'u')) then
            if (incx.eq.1) then
               do 20, j = n, 1, -1
                  if (x(j).ne.0d0) then
                     if (nounit) x(j) = x(j)/a(j, j)
                     temp = x(j)
                     do 10, i = j - 1, 1, -1
                        x(i) = x(i) - temp*a(i, j)
   10                continue
                  end if
   20          continue
            else
               jx = kx + (n - 1)*incx
               do 40, j = n, 1, -1
                  if (x(jx).ne.0d0) then
                     if (nounit) x(jx) = x(jx)/a(j, j)
                     temp = x(jx)
                     ix   = jx
                     do 30, i = j - 1, 1, -1
                        ix      = ix      - incx
                        x(ix) = x(ix) - temp*a(i, j)
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if (incx.eq.1) then
               do 60, j = 1, n
                  if (x(j).ne.0d0) then
                     if (nounit)
     *                  x(j) = x(j)/a(j, j)
                     temp = x(j)
                     do 50, i = j + 1, n
                        x(i) = x(i) - temp*a(i, j)
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if (x(jx).ne.0d0) then
                     if (nounit) x(jx) = x(jx)/a(j, j)
                     temp = x(jx)
                     ix   = jx
                     do 70, i = j + 1, n
                        ix      = ix      + incx
                        x(ix) = x(ix) - temp*a(i, j)
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else

c        form  x := inv(a')*x.

         if ((uplo.eq.'u'.or.uplo.eq.'u')) then
            if (incx.eq.1) then
               do 100, j = 1, n
                  temp = x(j)
                  do 90, i = 1, j - 1
                     temp = temp - a(i, j)*x(i)
   90             continue
                  if (nounit)
     *               temp = temp/a(j, j)
                  x(j) = temp
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x(jx)
                  ix   = kx
                  do 110, i = 1, j - 1
                     temp = temp - a(i, j)*x(ix)
                     ix   = ix   + incx
  110             continue
                  if (nounit)
     *               temp = temp/a(j, j)
                  x(jx) = temp
                  jx      = jx   + incx
  120          continue
            end if
         else
            if (incx.eq.1) then
               do 140, j = n, 1, -1
                  temp = x(j)
                  do 130, i = n, j + 1, -1
                     temp = temp - a(i, j)*x(i)
  130             continue
                  if (nounit)
     *               temp = temp/a(j, j)
                  x(j) = temp
  140          continue
            else
               kx = kx + (n - 1)*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x(jx)
                  ix   = kx
                  do 150, i = n, j + 1, -1
                     temp = temp - a(i, j)*x(ix)
                     ix   = ix   - incx
  150             continue
                  if (nounit) temp = temp/a(j, j)
                  x(jx) = temp
                  jx      = jx   - incx
  160          continue
            end if
         end if
      end if

      end

      subroutine dtrsv1 (trans, n, a, lda, x)
c----------------------------------------------------------------------
c modified dtrsv blas routine
c----------------------------------------------------------------------
      implicit none

      integer lda, n,  i, j

      character*1 trans

      double precision a(lda,*), x(*), temp
c----------------------------------------------------------------------
      if (n.eq.0) return

      if (trans.eq.'n') then

               do j = n, 1, -1
                  if (x(j).ne.0d0) then
                     x(j) = x(j)/a(j,j)
                     temp = x(j)
                     do i = j - 1, 1, -1
                        x(i) = x(i) - temp*a(i,j)
                     end do 
                  end if
               end do 

      else
               do j = 1, n
                  temp = x(j)
                  do i = 1, j - 1
                     temp = temp - a(i,j)*x(i)
                  end do 
                  temp = temp/a(j,j)
                  x(j) = temp
               end do 

      end if

      end

      subroutine dger1 (m, n, x, y, a, lda)
c----------------------------------------------------------------------
c modified blas routine dger:
c 1) alpha dropped as an argument and fixed to -1
c 2) incx, incy dropped arguments and fixed to 1
c----------------------------------------------------------------------
      implicit none

      integer lda, m, n, i, j, jy

      double precision a(lda,*), x(*), y(*), alpha, temp
c----------------------------------------------------------------------
      if (m.eq.0.or.n.eq.0) return

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.

      jy = 1
      alpha = -1d0

      do j = 1, n
            if (y(jy).ne.0d0) then
               temp = alpha*y(jy)
               do i = 1, m
                  a(i,j) = a(i,j) + x(i)*temp
               end do
            end if
            jy = jy + 1
      end do 

      end

      subroutine dtrmv (uplo, trans, diag, n, a, lda, x, incx)
c----------------------------------------------------------------------
c  dtrmv  does one of the matrix-vector operations

c     x := a*x,   or   x := a'*x,

c  where x is n element vector and a is an n by n unit, or non-unit,
c  upper or lower triangular matrix.

c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:

c              uplo = 'u' or 'u'   a is an upper triangular matrix.

c              uplo = 'l' or 'l'   a is a lower triangular matrix.

c           unchanged on exit.

c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:

c              trans = 'n' or 'n'   x := a*x.

c              trans = 't' or 't'   x := a'*x.

c              trans = 'c' or 'c'   x := a'*x.

c           unchanged on exit.

c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:

c              diag = 'u' or 'u'   a is assumed to be unit triangular.

c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.

c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least 0d0.
c           unchanged on exit.

c  a      - double precision array of dimension (lda, n).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.

c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max(1, n).
c           unchanged on exit.

c  x      - double precision array of dimension at least
c           (1 + (n - 1)*abs(incx)).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.

c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be 0d0.
c           unchanged on exit.
c----------------------------------------------------------------------
      implicit none

      integer incx, lda, n
      character*1 diag, trans, uplo
      double precision a(lda, *), x(*)
      double precision temp
      integer i, info, ix, j, jx, kx
      logical            nounit
c----------------------------------------------------------------------

      info = 0
      if     (.not.(uplo .eq.'u').and.
     *         .not.(uplo .eq.'l')    ) then
         info = 1
      else if (.not.(trans.eq.'n').and.
     *         .not.(trans.eq.'t').and.
     *         .not.(trans.eq.'c')    ) then
         info = 2
      else if (.not.(diag .eq.'u').and.
     *         .not.(diag .eq.'n')    ) then
         info = 3
      else if (n.lt.0) then
         info = 4
      else if (lda.lt.max(1, n)) then
         info = 6
      else if (incx.eq.0) then
         info = 8
      end if

      if (info.ne.0) return

      if (n.eq.0) return

      nounit = (diag.eq.'n'.or.diag.eq.'n')

c     set up the start point in x if the increment is not unity. this
c     will be  (n - 1)*incx  too small for descending loops.

      if (incx.le.0) then
         kx = 1 - (n - 1)*incx
      else if (incx.ne.1) then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if ((trans.eq.'n'.or.trans.eq.'n')) then
c
c        form  x := a*x.
c
         if ((uplo.eq.'u'.or.uplo.eq.'u')) then
            if (incx.eq.1) then
               do 20, j = 1, n
                  if (x(j).ne.0d0) then
                     temp = x(j)
                     do 10, i = 1, j - 1
                        x(i) = x(i) + temp*a(i, j)
   10                continue
                     if (nounit)
     *                  x(j) = x(j)*a(j, j)
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if (x(jx).ne.0d0) then
                     temp = x(jx)
                     ix   = kx
                     do 30, i = 1, j - 1
                        x(ix) = x(ix) + temp*a(i, j)
                        ix      = ix      + incx
   30                continue
                     if (nounit)
     *                  x(jx) = x(jx)*a(j, j)
                  end if
                  jx = jx + incx
   40          continue
            end if
         else
            if (incx.eq.1) then
               do 60, j = n, 1, -1
                  if (x(j).ne.0d0) then
                     temp = x(j)
                     do 50, i = n, j + 1, -1
                        x(i) = x(i) + temp*a(i, j)
   50                continue
                     if (nounit)
     *                  x(j) = x(j)*a(j, j)
                  end if
   60          continue
            else
               kx = kx + (n - 1)*incx
               jx = kx
               do 80, j = n, 1, -1
                  if (x(jx).ne.0d0) then
                     temp = x(jx)
                     ix   = kx
                     do 70, i = n, j + 1, -1
                        x(ix) = x(ix) + temp*a(i, j)
                        ix      = ix      - incx
   70                continue
                     if (nounit)
     *                  x(jx) = x(jx)*a(j, j)
                  end if
                  jx = jx - incx
   80          continue
            end if
         end if
      else

c        form  x := a'*x.

         if ((uplo.eq.'u'.or.uplo.eq.'u')) then
            if (incx.eq.1) then
               do 100, j = n, 1, -1
                  temp = x(j)
                  if (nounit)
     *               temp = temp*a(j, j)
                  do 90, i = j - 1, 1, -1
                     temp = temp + a(i, j)*x(i)
   90             continue
                  x(j) = temp
  100          continue
            else
               jx = kx + (n - 1)*incx
               do 120, j = n, 1, -1
                  temp = x(jx)
                  ix   = jx
                  if (nounit)
     *               temp = temp*a(j, j)
                  do 110, i = j - 1, 1, -1
                     ix   = ix   - incx
                     temp = temp + a(i, j)*x(ix)
  110             continue
                  x(jx) = temp
                  jx      = jx   - incx
  120          continue
            end if
         else
            if (incx.eq.1) then
               do 140, j = 1, n
                  temp = x(j)
                  if (nounit)
     *               temp = temp*a(j, j)
                  do 130, i = j + 1, n
                     temp = temp + a(i, j)*x(i)
  130             continue
                  x(j) = temp
  140          continue
            else
               jx = kx
               do 160, j = 1, n
                  temp = x(jx)
                  ix   = jx
                  if (nounit)
     *               temp = temp*a(j, j)
                  do 150, i = j + 1, n
                     ix   = ix   + incx
                     temp = temp + a(i, j)*x(ix)
  150             continue
                  x(jx) = temp
                  jx      = jx   + incx
  160          continue
            end if
         end if
      end if

      end

      subroutine ssrot1 (n, alpha, x, c, s)
c----------------------------------------------------------------------
c modified blas routine ssrotg (ssrotg)
c----------------------------------------------------------------------
      implicit none
      double precision c(*), s(*), x(*), alpha
      integer n, i, ix
c----------------------------------------------------------------------
      if (n.gt.0) then

            ix = 1

         do i = 1, n - 1
                  call srotg1 (x(ix + 1), x(ix), c(i), s(i))
                  s(i)  = -s(i)
                  x(ix) = -x(ix)
                  ix      =  ix + 1
         end do 
               call srotg1 (alpha, x(ix), c(n), s(n))
               s(n)  = -s(n)
               x(ix) = -x(ix)

      end if

      end

      subroutine ssrotg (pivot, direct, n, alpha, x, incx, c, s)
c----------------------------------------------------------------------
c  ssrotg generates the parameters of an orthogonal matrix p such that
c
c     when   pivot = 'f' or 'f'   and   direct = 'f' or 'f'
c     or     pivot = 'v' or 'v'   and   direct = 'b' or 'b'
c
c        p*(alpha) = (beta),
c          (  x  )   (  0 )
c
c     when   pivot = 'f' or 'f'   and   direct = 'b' or 'b'
c     or     pivot = 'v' or 'v'   and   direct = 'f' or 'f'
c
c        p*(  x  ) = (  0 ),
c          (alpha) = (beta)
c
c  where alpha is a scalar and x is an n element vector.
c
c  when  pivot = 'f' or 'f'  (fixed pivot)
c  and  direct = 'f' or 'f'  (forward sequence) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p(n)*p(n - 1)*...*p(1)
c
c     where p(k) is a plane rotation matrix for the (1, k + 1) plane
c     designed to annihilate the kth element of x.
c
c  when  pivot = 'v' or 'v'  (variable pivot)
c  and  direct = 'b' or 'b'  (backward sequence) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p(1)*p(2)*...*p(n)
c
c     where p(k) is a plane rotation matrix for the (k, k + 1) plane
c     designed to annihilate the kth element of x.
c
c  when  pivot = 'f' or 'f'  (fixed pivot)
c  and  direct = 'b' or 'b'  (backward sequence) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p(1)*p(2)*...*p(n)
c
c     where p(k) is a plane rotation matrix for the (k, n + 1) plane
c     designed to annihilate the kth element of x.
c
c  when  pivot = 'v' or 'v'  (variable pivot)
c  and  direct = 'f' or 'f'  (forward sequence) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p(n)*p(n - 1)*...*p(1)
c
c     where p(k) is a plane rotation matrix for the (k, k + 1) plane
c     designed to annihilate the kth element of x.
c
c  the routine returns the cosine, c(k), and sine, s(k) that define
c  the matrix p(k), such that the two by two rotation part of p(k),
c  r(k), has the form
c
c     r(k) = ( c(k)  s(k)).
c              (-s(k)  c(k))
c
c  on entry, alpha must contain  the scalar alpha and on exit, alpha is
c  overwritten by beta. the cosines and sines are returned in the arrays
c  c and s and the vector x is overwritten by the tangents of the plane
c  rotations (t(k) = s(k)/c(k)).
c----------------------------------------------------------------------
      implicit none 

      double precision alpha
      integer incx, n
      character direct*1, pivot*1
      double precision c(*), s(*), x(*)

      integer i, ix
c----------------------------------------------------------------------
      if (n.gt.0) then
         if ((direct.eq.'b').or.(direct.eq.'B')) then
            ix = 1 + (n - 1)*incx
            if ((pivot.eq.'v').or.(pivot.eq.'V')) then
               do 10, i = n, 2, -1
                  call srotg1 (x(ix-incx), x(ix), c(i), s(i))
                  ix = ix - incx
   10          continue
               call srotg1 (alpha, x(ix), c(1), s(1))
            else if ((pivot.eq.'f').or.(pivot.eq.'F')) then

c              here choose c and s so that
c
c                 (alpha) := ( c  s)*(alpha )
c                 (  0  )    (-s  c) (x(i))
c
c              which is equivalent to
c
c                 (  0  ) := (c  -s)*(x(i))
c                 (alpha)    (s   c) (alpha )
c
c              and so need to return  s(i) = -s  in order to make
c              r(i) look like
c
c                 r(i) = ( c(i)  s(i)).
c                          (-s(i)  c(i))
c
               do 20, i = n, 1, -1
                  call srotg1 (alpha, x(ix), c(i), s(i))
                  s(i)  = -s(i)
                  x(ix) = -x(ix)
                  ix      =  ix      - incx
   20          continue
            end if
         else if ((direct.eq.'f').or.(direct.eq.'F')) then
            ix = 1
            if ((pivot.eq.'v').or.(pivot.eq.'V')) then

c              here choose c and s so that

c                 (x(i + 1)) := ( c  s)*(x(i + 1))
c                 (   0      )    (-s  c) (x(i)    )

c              which is equivalent to

c                 (   0      ) := (c  -s)*(x(i)    )
c                 (x(i + 1))    (s   c) (x(i + 1))
c
c              and so need to return  s(i) = -s  in order to make
c              r(i) look like

c                 r(i) = ( c(i)  s(i)).
c                          (-s(i)  c(i))

               do 30, i = 1, n - 1
                  call srotg1 (x(ix+incx), x(ix), c(i), s(i))
                  s(i)  = -s(i)
                  x(ix) = -x(ix)
                  ix      =  ix      + incx
   30          continue
               call srotg1  (alpha, x(ix), c(n), s(n))
               s(n)  = -s(n)
               x(ix) = -x(ix)
            else if ((pivot.eq.'f').or.(pivot.eq.'F')) then
               do 40, i = 1, n
                  call srotg1 (alpha, x(ix), c(i), s(i))
                  ix = ix + incx
   40          continue
            end if
         end if
      end if

      end

      subroutine sgesrc (side, pivot, direct, m, n, k1, k2, c, s, a,
     *                   lda)
c----------------------------------------------------------------------
c  sgesrc  does the transformation
c
c     a := p*a,   when   side = 'l' or 'l'  ( left-hand side)
c
c     a := a*p',  when   side = 'r' or 'r'  (right-hand side)
c
c  where a is an m by n matrix and p is an orthogonal matrix, consisting
c  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
c  determined by the parameters pivot and direct as follows:
c
c     when  pivot  = 'v' or 'v'  (variable pivot)
c     and   direct = 'f' or 'f'  (forward sequence) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p(k2 - 1)*...*p(k1 + 1)*p(k1),
c
c        where  p(k)  is a plane rotation matrix for the  (k, k + 1)
c        plane.
c
c     when  pivot  = 'v' or 'v'  (variable pivot)
c     and   direct = 'b' or 'b'  (backward sequence) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p(k1)*p(k1 + 1)*...*p(k2 - 1),
c
c        where  p(k)  is a plane rotation matrix for the  (k, k + 1)
c        plane.
c
c     when  pivot  = 't' or 't'  (top pivot)
c     and   direct = 'f' or 'f'  (forward sequence) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p(k2 - 1)*p(k2 - 2)*...*p(k1),
c
c        where  p(k)  is a plane rotation matrix for the (k1, k + 1)
c        plane.
c
c     when  pivot  = 't' or 't'  (top pivot)
c     and   direct = 'b' or 'b'  (backward sequence) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p(k1)*p(k1 + 1)*...*p(k2 - 1),
c
c        where  p(k)  is a plane rotation matrix for the (k1, k + 1)
c        plane.
c
c     when  pivot  = 'b' or 'b'  (bottom pivot)
c     and   direct = 'f' or 'f'  (forward sequence) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p(k2 - 1)*p(k2 - 2)*...*p(k1),
c
c        where  p(k)  is a  plane rotation  matrix  for the  (k, k2)
c        plane.
c
c     when  pivot  = 'b' or 'b'  (bottom pivot)
c     and   direct = 'b' or 'b'  (backward sequence) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p(k1)*p(k1 + 1)*...*p(k2 - 1),
c
c        where  p(k)  is a  plane rotation  matrix  for the  (k, k2)
c        plane.
c
c  c(k) and s(k)  must contain the  cosine and sine  that define the
c  matrix  p(k).  the  two by two  plane rotation  part of the  matrix
c  p(k), r(k), is assumed to be of the form
c
c     r(k) = ( c(k)  s(k)).
c              (-s(k)  c(k))
c
c  if m, n or k1 are less than unity,  or k2 is not greater than k1,  or
c  side = 'l' or 'l'  and  k2  is greater than  m, or  side = 'r' or 'r'
c  and  k2  is greater than  n,  then an  immediate return  is effected.
c----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, m, n
      character*1 direct, pivot, side
      double precision a(lda, *), c(*), s(*)
      double precision aij, ctemp, stemp, temp
      integer i, j
      logical            left, right
c----------------------------------------------------------------------

      left = (side.eq.'l').or.(side.eq.'L')
      right = (side.eq.'r').or.(side.eq.'R')
      if ((min(m, n, k1).lt.1).or.(k2.le.k1).or.
     *    ((left).and.(k2.gt.m)).or.
     *    ((right).and.(k2.gt.n)))return
      if (left) then
         if ((pivot.eq.'v').or.(pivot.eq.'V')) then
            if ((direct.eq.'f').or.(direct.eq.'F')) then
               do 20 j = 1, n
                  aij = a(k1, j)
                  do 10 i = k1, k2 - 1
                     temp = a(i + 1, j)
                     a(i, j) = s(i)*temp + c(i)*aij
                     aij = c(i)*temp - s(i)*aij
   10             continue
                  a(k2, j) = aij
   20          continue
            else if ((direct.eq.'b').or.(direct.eq.'B')) then
               do 40 j = 1, n
                  aij = a(k2, j)
                  do 30 i = k2 - 1, k1, -1
                     temp = a(i, j)
                     a(i + 1, j) = c(i)*aij - s(i)*temp
                     aij = s(i)*aij + c(i)*temp
   30             continue
                  a(k1, j) = aij
   40          continue
            end if
         else if ((pivot.eq.'t').or.(pivot.eq.'T')) then
            if ((direct.eq.'f').or.(direct.eq.'F')) then
               do 60 j = 1, n
                  temp = a(k1, j)
                  do 50 i = k1, k2 - 1
                     aij = a(i + 1, j)
                     a(i + 1, j) = c(i)*aij - s(i)*temp
                     temp = s(i)*aij + c(i)*temp
   50             continue
                  a(k1, j) = temp
   60          continue
            else if ((direct.eq.'b').or.(direct.eq.'B')) then
               do 80 j = 1, n
                  temp = a(k1, j)
                  do 70 i = k2 - 1, k1, -1
                     aij = a(i + 1, j)
                     a(i + 1, j) = c(i)*aij - s(i)*temp
                     temp = s(i)*aij + c(i)*temp
   70             continue
                  a(k1, j) = temp
   80          continue
            end if
         else if ((pivot.eq.'b').or.(pivot.eq.'B')) then
            if ((direct.eq.'f').or.(direct.eq.'F')) then
               do 100 j = 1, n
                  temp = a(k2, j)
                  do 90 i = k1, k2 - 1
                     aij = a(i, j)
                     a(i, j) = s(i)*temp + c(i)*aij
                     temp = c(i)*temp - s(i)*aij
   90             continue
                  a(k2, j) = temp
  100          continue
            else if ((direct.eq.'b').or.(direct.eq.'B')) then
               do 120 j = 1, n
                  temp = a(k2, j)
                  do 110 i = k2 - 1, k1, -1
                     aij = a(i, j)
                     a(i, j) = s(i)*temp + c(i)*aij
                     temp = c(i)*temp - s(i)*aij
  110             continue
                  a(k2, j) = temp
  120          continue
            end if
         end if
      else if (right) then
         if ((pivot.eq.'v').or.(pivot.eq.'V')) then
            if ((direct.eq.'f').or.(direct.eq.'F')) then
               do 140 j = k1, k2 - 1
                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 130 i = 1, m
                        temp = a(i, j + 1)
                        a(i, j + 1) = ctemp*temp - stemp*a(i, j)
                        a(i, j) = stemp*temp + ctemp*a(i, j)
  130                continue
                  end if
  140          continue
            else if ((direct.eq.'b').or.(direct.eq.'B')) then
               do 160 j = k2 - 1, k1, -1
                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 150 i = m, 1, -1
                        temp = a(i, j + 1)
                        a(i, j + 1) = ctemp*temp - stemp*a(i, j)
                        a(i, j) = stemp*temp + ctemp*a(i, j)
  150                continue
                  end if
  160          continue
            end if
         else if ((pivot.eq.'t').or.(pivot.eq.'T')) then
            if ((direct.eq.'f').or.(direct.eq.'F')) then
               do 180 j = k1 + 1, k2
                  ctemp = c(j - 1)
                  stemp = s(j - 1)
                  if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
                     do 170 i = 1, m
                        temp = a(i, j)
                        a(i, j) = ctemp*temp - stemp*a(i, k1)
                        a(i, k1) = stemp*temp + ctemp*a(i, k1)
  170                continue
                  end if
  180          continue
            else if ((direct.eq.'b').or.(direct.eq.'B')) then
               do 200 j = k2, k1 + 1, -1
                  ctemp = c(j - 1)
                  stemp = s(j - 1)
                  if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
                     do 190 i = m, 1, -1
                        temp = a(i, j)
                        a(i, j) = ctemp*temp - stemp*a(i, k1)
                        a(i, k1) = stemp*temp + ctemp*a(i, k1)
  190                continue
                  end if
  200          continue
            end if
         else if ((pivot.eq.'b').or.(pivot.eq.'B')) then
            if ((direct.eq.'f').or.(direct.eq.'F')) then
               do 220 j = k1, k2 - 1
                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 210 i = 1, m
                        temp = a(i, j)
                        a(i, j) = stemp*a(i, k2) + ctemp*temp
                        a(i, k2) = ctemp*a(i, k2) - stemp*temp
  210                continue
                  end if
  220          continue
            else if ((direct.eq.'b').or.(direct.eq.'B')) then
               do 240 j = k2 - 1, k1, -1
                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)
                     stemp = s(j)
                     do 230 i = m, 1, -1
                        temp = a(i, j)
                        a(i, j) = stemp*a(i, k2) + ctemp*temp
                        a(i, k2) = ctemp*a(i, k2) - stemp*temp
  230                continue
                  end if
  240          continue
            end if
         end if
      end if

      end


      subroutine sgesr1 (side, direct, m, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
c modified blas routine sgesrc.
c----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, m, n, i, j

      character*1 direct, side

      double precision a(lda, *), c(*), s(*), aij, ctemp, stemp, temp

      logical            left, right

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      left = (side.eq.'l')
      right = (side.eq.'r')

      if ((min(m, n, k1).lt.1).or.(k2.le.k1).or.
     *    (left .and.k2.gt.m).or.(right.and.k2.gt.n)) return

      if (left) then

            if (direct.eq.'f') then

               do j = 1, n
                  aij = a(k1,j)
                  do i = k1, k2 - 1
                     temp = a(i + 1,j)

                     if (c(i).lt.wmach(3)) c(i) = 0d0

                     a(i,j) = s(i)*temp + c(i)*aij
                     aij = c(i)*temp - s(i)*aij
                  end do 
                  a(k2,j) = aij
               end do

            else if (direct.eq.'b') then

               do j = 1, n
                  aij = a(k2,j)

                  do i = k2 - 1, k1, -1
                     temp = a(i,j)
                     if (c(i).lt.wmach(3)) c(i) = 0d0
                     a(i + 1,j) = c(i)*aij - s(i)*temp
                     aij = s(i)*aij + c(i)*temp
                  end do

                  a(k1,j) = aij
               end do

            end if

      else if (right) then

            if (direct.eq.'f') then

               do 140 j = k1, k2 - 1

                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)

                     if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

                     stemp = s(j)

                     do 130 i = 1, m
                        temp = a(i, j + 1)
                        a(i, j + 1) = ctemp*temp - stemp*a(i,j)
                        a(i,j) = stemp*temp + ctemp*a(i,j)
  130                continue

                  end if

  140          continue

            else if (direct.eq.'b') then

               do j = k2 - 1, k1, -1

                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then

                     ctemp = c(j)

                     if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

                     stemp = s(j)

                     do i = m, 1, -1
                        temp = a(i,j + 1)
                        a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                        a(i,j) = stemp*temp + ctemp*a(i,j)
                     end do

                  end if

               end do

            end if

      end if

      end

      subroutine sgeapr (side, trans, n, perm, k, b, ldb)
c----------------------------------------------------------------------
c  sgeapr does one of the transformations
c
c     b := p'*b   or   b := p*b,   where b is an m by k matrix,
c
c  or
c
c     b := b*p'   or   b := b*p,   where b is a k by m matrix,
c
c  p being an m by m permutation matrix of the form
c
c     p = p(1, index(1))*p(2, index(2))*...*p(n, index(n)),
c
c  where  p(i, index(i)) is the permutation matrix that interchanges
c  items i and index(i). that is p(i, index(i)) is the unit matrix
c  with rows and columns  i and  index(i)  interchanged. of course, if
c  index(i) = i  then  p(i, index(i)) = i.
c
c  this  routine is intended  for use in conjunction with auxiliary
c  routines  that  do  interchange  operations,  such  as  sorting.

c  side   - character*1.
c  trans
c           on entry,  side  (left-hand side, or right-hand side)  and
c           trans  (transpose, or no transpose)  specify the operation
c           to be performed as follows.
c
c           side = 'l' or 'l'   and   trans = 't' or 't'
c
c              do the operation   b := p'*b.
c
c           side = 'l' or 'l'   and   trans = 'n' or 'n'
c
c              do the operation   b := p*b.
c
c           side = 'r' or 'r'   and   trans = 't' or 't'
c
c              do the operation   b := b*p'.
c
c           side = 'r' or 'r'   and   trans = 'n' or 'n'
c
c              do the operation   b := b*p.
c
c           unchanged on exit.

c  n      - integer.

c           on entry, n must specify the value of n.  n must be at least
c           zero.  when  n = 0  then an  immediate  return  is effected.

c           unchanged on exit.

c  perm   - real             array of dimension at least (n).

c           before  entry,  perm  must  contain  the  n indices  for the
c           permutation matrices. index(i) must satisfy

c              1 .le. index(i) .le. m.

c           it is usual for index(i) to be at least i, but this is not
c           necessary for this routine. it is assumed that the statement
c           index = perm(i)  returns the correct integer in  index, so
c           that,  if necessary,  perm(i)  should contain a real value
c           slightly larger than  index.

c           unchanged on exit.

c  k      - integer.

c           on entry with  side = 'l' or 'l',  k must specify the number
c           of columns of b and on entry with  side = 'r' or 'r', k must
c           specify the number of rows of  b.  k must be at least  zero.
c           when  k = 0  then an immediate return is effected.

c           unchanged on exit.
c
c  b      - real  array  of  dimension (ldb, ncolb),  where  ncolb = k
c           when  side = 'l' or 'l'  and  ncolb = m  when  side = 'r' or
c           'r'.

c           before entry  with  side = 'l' or 'l',  the  leading  m by k
c           part  of  the  array   b  must  contain  the  matrix  to  be
c           transformed  and before  entry with  side = 'r' or 'r',  the
c           leading  k by m part of the array  b must contain the matrix
c           to  be  transformed.  on exit,   b  is  overwritten  by  the
c           transformed matrix.

c  ldb    - integer.

c           on entry,  ldb  must specify  the  leading dimension  of the
c           array  b  as declared  in the  calling  (sub) program.  when
c           side = 'l' or 'l'   then  ldb  must  be  at  least  m,  when
c           side = 'r' or 'r'   then  ldb  must  be  at  least  k.
c           unchanged on exit.
c----------------------------------------------------------------------
      implicit none

      integer k, ldb, n
      character side*1, trans*1
      double precision perm(*), b(ldb, *)
      logical            left, null, right, trnsp
      integer i, j, l
      double precision temp
c----------------------------------------------------------------------

      if (min(n, k).eq.0) return

      left = (side.eq.'l').or.(side.eq.'l')
      right = (side.eq.'r').or.(side.eq.'r')
      null = (trans.eq.'n').or.(trans.eq.'n')
      trnsp = (trans.eq.'t').or.(trans.eq.'t')
      if (left) then
         if (trnsp) then
            do 20 i = 1, n
               l = perm(i)
               if (l.ne.i) then
                  do 10 j = 1, k
                     temp = b(i, j)
                     b(i, j) = b(l, j)
                     b(l, j) = temp
   10             continue
               end if
   20       continue
         else if (null) then
            do 40 i = n, 1, -1
               l = perm(i)
               if (l.ne.i) then
                  do 30 j = 1, k
                     temp = b(l, j)
                     b(l, j) = b(i, j)
                     b(i, j) = temp
   30             continue
               end if
   40       continue
         end if
      else if (right) then
         if (trnsp) then
            do 60 j = n, 1, -1
               l = perm(j)
               if (l.ne.j) then
                  do 50 i = 1, k
                     temp = b(i, j)
                     b(i, j) = b(i, l)
                     b(i, l) = temp
   50             continue
               end if
   60       continue
         else if (null) then
            do 80 j = 1, n
               l = perm(j)
               if (l.ne.j) then
                  do 70 i = 1, k
                     temp = b(i, l)
                     b(i, l) = b(i, j)
                     b(i, j) = temp
   70             continue
               end if
   80       continue
         end if
      end if

      end

      subroutine sgeap1 (n, perm, k, b, ldb)
c----------------------------------------------------------------------
c modified blas routine sgeapr
c----------------------------------------------------------------------
      implicit none

      integer k, ldb, n, i, j, l

      double precision perm(*), b(ldb, *), temp
c----------------------------------------------------------------------
      if (min(n, k).eq.0) return

      do i = 1, n

         l = idint(perm(i))

         if (l.ne.i) then

            do j = 1, k
               temp = b(i,j)
               b(i,j) = b(l,j)
               b(l,j) = temp
            end do

         end if

      end do

      end

      subroutine sgrfg (n, alpha, x, incx, tol, zeta)
c----------------------------------------------------------------------
c  sgrfg generates details of a generalized householder reflection such
c  that

c     p*(alpha) = (beta),   p'*p = i.
c       (  x  )   (  0 )

c  p is given in the form

c     p = i - (zeta)*(zeta  z'),
c             (  z )

c  where z is an n element vector and zeta is a scalar that satisfies

c     1.0 .le. zeta .le. sqrt(2.0).

c  zeta is returned in zeta unless x is such that

c     max(abs(x(i))) .le. max(eps*abs(alpha), tol)

c  where eps is the relative machine precision and tol is the user
c  supplied value tol, in which case zeta is returned as 0.0 and p can
c  be taken to be the unit matrix.

c  beta is overwritten on alpha and z is overwritten on x.
c  the routine may be called with  n = 0  and advantage is taken of the
c  case where  n = 1.
c----------------------------------------------------------------------
      implicit none 

      double precision alpha, tol, zeta
      integer incx, n
      double precision x(*)
      double precision beta, scale, ssq

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if (n.lt.1) then
         zeta = 0d0
      else if ((n.eq.1).and.(x(1).eq.0d0)) then
         zeta = 0d0
      else

c        treat case where p is a 2 by 2 matrix specially.
c
         if (n.eq.1) then
c
c           deal with cases where  alpha = zero  and
c           abs(x(1)) .le. max(eps*abs(alpha), tol)  first.
c
            if (alpha.eq.0d0) then
               zeta   =  1d0
               alpha  = dabs (x(1))
               x(1) = -dsign (1d0, x(1))
            else if (dabs(x(1)).le.dmax1(wmach(3)*dabs(alpha),tol)) then
               zeta   =  0d0
            else
               if (dabs(alpha).ge.dabs(x(1))) then
                  beta = dabs(alpha) *dsqrt(1 + (x(1)/alpha)**2)
               else
                  beta = dabs(x(1))*dsqrt(1 + (alpha/x(1))**2)
               end if
               zeta = dsqrt((dabs(alpha) + beta)/beta)
               if (alpha.ge.0d0)
     *            beta = -beta
               x(1) = -x(1)/(zeta*beta)
               alpha  = beta
            end if
         else

c           now p is larger than 2 by 2.

            ssq   = 1d0
            scale = 0d0
            call sssq (n, x, incx, scale, ssq)

c           treat cases where  scale = 0d0,
c           scale .le. max(eps*abs(alpha), tol)  and
c           alpha = 0d0  specially.
c           note that  scale = max(abs(x(i))).

            if ((scale.eq.0d0).or.
     *          (scale.le.dmax1(wmach(3)*dabs(alpha), tol))) then
               zeta  = 0d0
            else if (alpha.eq.0d0) then
               zeta  = 1d0
               alpha = scale*dsqrt(ssq)
               call dscal (n, -1/alpha, x, incx)
            else
               if (scale.lt.dabs(alpha)) then
                  beta = dabs(alpha)*dsqrt(1 + ssq*(scale/alpha)**2)
               else
                  beta = scale *dsqrt(ssq + (alpha/scale)**2)
               end if
               zeta = dsqrt((beta + dabs(alpha))/beta)
               if (alpha.gt.0d0)
     *            beta = -beta
               call dscal (n, -1/(zeta*beta), x, incx)
               alpha = beta
            end if
         end if
      end if

      end

      subroutine sgrfg1 (n, alpha, x, zeta)
c----------------------------------------------------------------------
c modified blas sgrfg routine
c----------------------------------------------------------------------
      implicit none

      double precision alpha, zeta, beta, scale, ssq, x(*)

      integer n

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if (n.lt.1) then

         zeta = 0d0

      else if (n.eq.1.and.x(1).eq.0d0) then

         zeta = 0d0

      else
c        treat case where p is a 2 by 2 matrix specially.

         if (n.eq.1) then

c           deal with cases where  alpha = 0d0  and
c           dabs(x(1)) .le. max(eps*abs(alpha), tol)  first.

            if (alpha.eq.0d0) then
               zeta   =  1d0
               alpha  =  dabs (x(1))
               x(1) = -dsign(1d0, x(1))
            else if (dabs(x(1)).le.max(wmach(3)*dabs(alpha),0d0)) then
               zeta   =  0d0
            else

               if (dabs(alpha).ge.dabs(x(1))) then
                  beta = dabs(alpha)*dsqrt(1d0 + (x(1)/alpha)**2)
               else
                  beta = dabs(x(1))*dsqrt(1d0 + (alpha/x(1))**2)
               end if

               zeta = dsqrt((dabs(alpha) + beta)/beta)
               if (alpha.ge.0d0) beta = -beta
               x(1) = -x(1)/(zeta*beta)
               alpha  = beta

            end if
         else

c           now p is larger than 2 by 2.

            ssq   = 1d1
            scale = 0d0

            call sssq (n, x, 1, scale, ssq)

c           treat cases where  scale = 0d0,
c           scale .le. max(eps*abs(alpha), tol)  and
c           alpha = 0d0  specially.
c           note that  scale = max(abs(x(i))).

            if (scale.eq.0d0.or.
     *          scale.le.max(wmach(3)*dabs(alpha),0d0)) then

               zeta  = 0d0

            else if (alpha.eq.0d0) then

               zeta  = 1d0
               alpha = scale*dsqrt(ssq)
               call dscal (n, -1d0/alpha, x, 1)

            else

               if (scale.lt.dabs(alpha)) then
                  beta = dabs(alpha)*dsqrt(1d0 + ssq*(scale/alpha)**2)
               else
                  beta = scale       *dsqrt(ssq +   (alpha/scale)**2)
               end if

               zeta = dsqrt((beta + dabs(alpha))/beta)
               if (alpha.gt.0d0) beta = -beta
               call dscal (n,-1d0/(zeta*beta), x, 1)
               alpha = beta

            end if

         end if

      end if

      end

      subroutine srotg1 (a, b, c, s)
c----------------------------------------------------------------------
c  modifed blas routine srotg (srotgc/srotg1): c is always returned as non-negative and b 
c  is overwritten by the tangent of the angle that defines the plane rotation.
c----------------------------------------------------------------------
      double precision a, b, c, s, t, abst

      logical            fail
      double precision adivb
      external adivb

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if (b.eq.0d0) then

         c  = 1d0
         s  = 0d0

      else

         t  = adivb (b,a,fail)

         abst = dabs(t)

         if (abst.lt.wmach(4)) then
c                                 wmach(4) = sqrt(epsmch)
            c = 1d0
            s = t

         else if (abst*wmach(4).gt.1d0) then

            c = 1d0/abst
            s = sign (1d0,t)

         else

            c = 1d0 / dsqrt(1d0 + abst**2)
            s = c*t

         end if

         a  = c*a + s*b
         b  = t

      end if

      end

      subroutine sutsrs (side, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
c  sutsrs applies a  given sequence  of  plane rotations  to either  the
c  left,  or the right,  of the  n by n  upper triangular  matrix  u  to
c  transform  u  to an upper spiked matrix. the rotations are applied in
c  planes k1 up to k2.
c
c  the upper spiked matrix, h, is formed as
c
c     h = p*u,   when   side = 'l' or 'L',  (left-hand side)
c
c  where p is an orthogonal matrix of the form
c
c     p = p(k1)*p(k1 + 1)*...*p(k2 - 1),
c
c  p(k) being a plane rotation matrix for the (k, k2) plane, and is
c  formed as
c
c     h = u*p',   when   side = 'r' or 'R',  (right-hand side)
c
c  where p is an orthogonal matrix of the form
c
c     p = p(k2 - 1)*...*p(k1 + 1)*p(k1),
c
c  p(k)  being a  plane rotation matrix for the  (k1, k + 1)  plane.
c
c  the cosine and sine that define  p(k), k = k1, k1 + 1, ..., k2 - 1,
c  must be  supplied  in  c(k) and s(k) respectively. the two by two
c  rotation part of p(k), r(k), is assumed to have the form
c
c     r(k) = ( c(k)  s(k)).
c              (-s(k)  c(k))
c
c  the matrix  u must be supplied in the n by n leading upper triangular
c  part of the array  a, and this is overwritten by the upper triangular
c  part of h.
c
c  when  side = 'l' or 'l'  then a  row spike  is  generated  in  h  and
c  when  side = 'r' or 'r'  then a  column spike is generated. for a row
c  spike the elements  h(k2, k)  and for a  column spike  the elements
c  h(k + 1, k1), k = k1, k1 + 1, ..., k2 - 1, are returned in  s(k).
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c----------------------------------------------------------------------
      integer k1, k2, lda, n
      character*1 side
      double precision a(lda, *), c(*), s(*)
      double precision aij, ctemp, spike, stemp, temp
      integer i, j
c----------------------------------------------------------------------
      if ((min(n, k1).lt.1).or.(k2.le.k1).or.
     *   (k2.gt.n))return
      if ((side.eq.'l').or.(side.eq.'L')) then

c        apply the plane rotations to columns n back to k2.

         do 20 j = n, k2, -1
            temp = a(k2, j)
            do 10 i = k2 - 1, k1, -1
               aij = a(i, j)
               a(i, j) = s(i)*temp + c(i)*aij
               temp = c(i)*temp - s(i)*aij
   10       continue
            a(k2, j) = temp
   20    continue

c        form  the spike  and apply the rotations in columns  (k2 - 1)
c        back to k1.

         do 40 j = k2 - 1, k1, -1
            spike = -s(j)*a(j, j)
            a(j, j) = c(j)*a(j, j)
            do 30 i = j - 1, k1, -1
               aij = a(i, j)
               a(i, j) = s(i)*spike + c(i)*aij
               spike = c(i)*spike - s(i)*aij
   30       continue
            s(j) = spike
   40    continue
      else if ((side.eq.'r').or.(side.eq.'R')) then

c        apply the  plane rotations to columns  (k1 + 1) up to k2  and
c        form the spike.

         do 70 j = k1 + 1, k2
            ctemp = c(j - 1)
            stemp = s(j - 1)
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               do 50 i = 1, k1
                  temp = a(i, k1)
                  a(i, k1) = stemp*a(i, j) + ctemp*temp
                  a(i, j) = ctemp*a(i, j) - stemp*temp
   50          continue
               do 60 i = k1 + 1, j - 1
                  spike = s(i - 1)
                  s(i - 1) = stemp*a(i, j) + ctemp*spike
                  a(i, j) = ctemp*a(i, j) - stemp*spike
   60          continue
               s(j - 1) = stemp*a(j, j)
               a(j, j) = ctemp*a(j, j)
            end if
   70    continue
      end if

      end


      subroutine susqr(side, n, k1, k2, c, s, a, lda)
c  susqr restores an upper spiked matrix  h to upper triangular form by
c  applying a sequence of plane rotations, in planes  k1 up to k2,  from
c  either the left, or the right.

c  the matrix  h is assumed to have non-zero elements only in the spiked
c  positions, h(k2, k) for a row spike and h(k + 1, k1) for a column
c  spike, k = k1, k1 + 1, ..., k2 - 1, and these must be supplied in the
c  elements s(k).

c  when  side = 'l' or 'l'  ( left-hand side)

c     h  is  assumed  to have a  row spike  and is restored to the upper
c     triangular matrix  r as

c        r = p*h,

c     where p is an orthogonal matrix of the form

c        p = p(k2 - 1)*...*p(k1 + 1)*p(k1),

c     p(k)  being a  plane rotation  matrix for the  (k, k2)  plane.

c  when  side = 'r' or 'r'  (right-hand side)

c     h  is assumed to have a  column spike and is restored to the upper
c     triangular matrix r as

c        r = h*p',

c     where p is an orthogonal matrix of the form
c
c        p = p(k1)*p(k1 + 1)*...*p(k2 - 1),
c
c     p(k) being a plane rotation matrix for the  (k1, k + 1) plane.
c
c  the  two by two  rotation  part of  p(k),  q(k),  is of  the form
c
c     q(k) = ( c(k)  s(k))
c              (-s(k)  c(k))
c
c  and  c(k) and s(k) are returned in the kth elements of the arrays
c  c and s respectively.
c
c  the upper triangular part of the matrix  h must be supplied in the  n
c  by n  leading upper triangular part of  a, and this is overwritten by
c  the upper triangular matrix r.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      integer k1, k2, lda, n
      character*1 side
      double precision a(lda, *), c(*), s(*)
      double precision aij, ctemp, spike, stemp, temp
      integer i, j

      if ((min(n, k1).lt.1).or.(k2.le.k1).or.
     *   (k2.gt.n))return
      if ((side.eq.'l').or.(side.eq.'L')) then

c        restore h to upper triangular form by annihilating the elements
c        in  the  spike  of  h.  the  jth rotation  is  chosen  so  that

c        (h(j, j)) := ( c  s)*(h(j , j)).
c        (    0    )    (-s  c) (h(k2, j))

c        apply the rotations in columns k1 up to (k2 - 1).

         do 20 j = k1, k2 - 1
            spike = s(j)
            do 10 i = k1, j - 1
               aij = a(i, j)
               a(i, j) = s(i)*spike + c(i)*aij
               spike = c(i)*spike - s(i)*aij
   10       continue

c           set up the rotation.

            call srotg1(a(j, j), spike, c(j), s(j))
   20    continue

c        apply the rotations to columns k2 up to n.

         do 40 j = k2, n
            temp = a(k2, j)
            do 30 i = k1, k2 - 1
               aij = a(i, j)
               a(i, j) = s(i)*temp + c(i)*aij
               temp = c(i)*temp - s(i)*aij
   30       continue
            a(k2, j) = temp
   40    continue
      else if ((side.eq.'r').or.(side.eq.'R')) then

c        restore h to upper triangular form by annihilating the spike of
c        h. the jth rotation is chosen so that

c           (h(j, j)) := ( c  s)*(h(j, j) ),
c           (    0    )    (-s  c) (h(j, k1))

c        which can be expressed as

c           (0  h(j, j)) := (h(j, k1)  h(j, j))*( c  s).
c                                                         (-s  c)

c        thus return  c(j) = c  and  s(j) = -s  to make the plane
c        rotation matrix look like

c           q(j) = ( c(j)  s(j)).
c                    (-s(j)  c(j))

         do 70 j = k2, k1 + 1, -1
            call srotg1(a(j, j), s(j - 1), ctemp, stemp)
            stemp = -stemp
            s(j - 1) = stemp
            c(j - 1) = ctemp
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               do 50 i = j - 1, k1 + 1, -1
                  spike = s(i - 1)
                  s(i - 1) = stemp*a(i, j) + ctemp*spike
                  a(i, j) = ctemp*a(i, j) - stemp*spike
   50          continue
               do 60 i = k1, 1, -1
                  temp = a(i, k1)
                  a(i, k1) = stemp*a(i, j) + ctemp*temp
                  a(i, j) = ctemp*a(i, j) - stemp*temp
   60          continue
            end if
   70    continue
      end if

      end

      subroutine sload (n, const, x, incx)
c----------------------------------------------------------------------
c blas routine sload (f06fbf)
c----------------------------------------------------------------------
      implicit none

      double precision const, x(*)

      integer incx, n, ix
c----------------------------------------------------------------------
      if (n.gt.0) then

         if (const.ne.0d0) then

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = const
            end do

         else

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = 0d0
            end do

         end if

      end if

      end

      subroutine sscmv (n, alpha, x, incx, y, incy)
c----------------------------------------------------------------------
c  sscmv does the operation
c
c     y := alpha*x
c----------------------------------------------------------------------
      implicit none

      double precision alpha
      integer incx, incy, n
      double precision x(*), y(*)
      integer i, ix, iy
c----------------------------------------------------------------------
      if (n.gt.0) then
         if (alpha.eq.0d0.and.incy.ne.0) then
            call sload(n, 0d0, y, abs(incy))
         else
            if (incx.eq.incy.and.incx.gt.0) then
               do 10, ix = 1, 1 + (n - 1)*incx, incx
                  y(ix) = alpha*x(ix)
   10          continue
            else
               if (incy.ge.0) then
                  iy = 1
               else
                  iy = 1 - (n - 1)*incy
               end if
               if (incx.gt.0) then
                  do 20, ix = 1, 1 + (n - 1)*incx, incx
                     y(iy) = alpha*x(ix)
                     iy      = iy            + incy
   20             continue
               else
                  ix = 1 - (n - 1)*incx
                  do 30, i = 1, n
                     y(iy) = alpha*x(ix)
                     ix      = ix            + incx
                     iy      = iy            + incy
   30             continue
               end if
            end if
         end if
      end if

      end

      integer function isrank (n, x, incx, tol)
c----------------------------------------------------------------------
c  isrank finds the first element of the n element vector x for which
c
c     abs(x(k)).le.tol*max(abs(x(1)), ..., abs(x(k - 1)))
c
c  and returns the value (k - 1) in the function name isrank. if no
c  such k exists then isrank is returned as n.
c
c  if tol is supplied as less than zero then the value epsmch, where
c  epsmch is the relative machine precision, is used in place of tol.
c----------------------------------------------------------------------
      implicit none

      double precision tol
      integer     incx, n
      double precision x(*)
      double precision tl, xmax
      integer     ix, k

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------

      k = 0
      if (n.ge.1) then
         ix = 1
         if (tol.lt.0d0) then
c                                this should be (1/2)*b**(1-p)
c                                or maybe b**(1-p)
            tl = wmach(8)
         else
            tl = tol
         end if
         xmax = dabs(x(ix))
c
   10    if   (k.lt.n) then
            if (dabs(x(ix)).le.tl*xmax) go to 20
            xmax = max(xmax, dabs(x(ix)))
            k    = k  + 1
            ix   = ix + incx
            go to 10
         end if

      end if

   20 isrank = k

      end

      subroutine f06qzz (hess,n,k1,k2,c,s,a,lda)

c  f06qzz  either applies a  given sequence  of  plane rotations  to the
c  right of the n by n reverse lower triangular matrix t, to transform t
c  to a  reverse lower hessenberg matrix  h, or restores a reverse lower
c  hessenberg matrix h to reverse lower triangular form t, by applying a
c  sequence of plane rotations from the right.
c
c  the rotations are applied  in planes k1 up to k2.
c
c  when   hess = 'c' or 'c',   (create),  then   the   reverse   lower
c  hessenberg matrix, h, is formed as
c
c     h = t*p',
c
c  where p is an orthogonal matrix of the form
c
c     p = p(k2 - 1)*...*p(k1 + 1)*p(k1),
c
c  p(k) being a plane rotation matrix for the  (k, k + 1) plane. the
c  cosine and sine that define p(k), k = k1, k1 + 1, ..., k2 - 1, must
c  be  supplied  in  c(k)  and  s(k)  respectively.  the  two by two
c  rotation part of p(k), r(k), is assumed to have the form
c
c     r(k) = ( c(k)  s(k)).
c              (-s(k)  c(k))
c
c  the matrix  t must be supplied in the n by n reverse lower triangular
c  part  of the array  a,  and this is overwritten by the  reverse lower
c  triangular part of  h.
c
c  the super-diagonal elements of  h, h(n - k, k), are returned in the
c  elements s(k),  k = k1, k1 + 1, ..., k2 - 1.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c
c  when   hess = 'r' or 'r',   (remove),  then   the   reverse   lower
c  hessenberg matrix  h  is  assumed  to  have  non-zero  super-diagonal
c  elements  in  positions  h(n - k, k),  k = k1, k1 + 1, ..., k2 - 1,
c  only and  h(n - k, k) must be supplied in  s(k). h is restored to
c  the reverse lower triangular matrix t as
c
c     t = h*p',
c
c  where p is an orthogonal matrix of the form
c
c     p = p(k1)*p(k1 + 1)*...*p(k2 - 1),
c
c  p(k) being a plane rotation for the  (k, k + 1) plane. the cosine
c  and  sine  that  define  p(k)  are  returned  in  c(k) and s(k)
c  respectively.  the  two by two  rotation part of  p(k),  r(k), is
c  of the form
c
c     r(k) = ( c(k)  s(k)).
c              (-s(k)  c(k))
c
c  the reverse lower triangular part of the matrix h must be supplied in
c  the  n by n  reverse  lower  triangular  part  of  a,   and  this  is
c  overwritten by the reverse triangular matrix t.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c
c  when   n = 7, k1 = 2 and k2 = 5   then  t  and  h  are  of  the  form
c
c     t = (0  0  0  0  0  0  x),   h = (0  0  0  0  0  0  x).
c         (0  0  0  0  0  x  x)        (0  0  0  0  x  x  x)
c         (0  0  0  0  x  x  x)        (0  0  0  x  x  x  x)
c         (0  0  0  x  x  x  x)        (0  0  x  x  x  x  x)
c         (0  0  x  x  x  x  x)        (0  x  x  x  x  x  x)
c         (0  x  x  x  x  x  x)        (0  x  x  x  x  x  x)
c         (x  x  x  x  x  x  x)        (x  x  x  x  x  x  x)
c----------------------------------------------------------------------
      integer k1, k2, lda, n
      character*1 hess
      double precision a(lda,*), c(*), s(*)
      double precision ctemp, stemp, suph, temp
      integer i, j
c----------------------------------------------------------------------

      if ((min(n,k1).lt.1).or.(k2.le.k1).or.(k2.gt.n)) return
      if ((hess.eq.'c').or.(hess.eq.'C')) then

c        apply  the  plane rotations  to  columns  k1  up to  (k2 - 1)
c        and  form   the  additional  super-diagonal  elements,  storing
c        h(n - j, j) in s(j).

         do 40 j = k1, k2 - 1
            if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
               stemp = s(j)
               ctemp = c(j)
               s(j) = stemp*a(n-j,j+1)
               a(n-j,j+1) = ctemp*a(n-j,j+1)
               do 20 i = n - j + 1, n
                  temp = a(i,j+1)
                  a(i,j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   20          continue
            end if
   40    continue
      else if ((hess.eq.'r').or.(hess.eq.'R')) then
c
c        restore  h to reverse lower triangular form by annihilating the
c        super-diagonal elements of  h.  the  jth rotation  is chosen so
c        that
c
c          (h(n - j, n - j)) := ( c  s)*(h(n - j, n - j    )),
c          (        0        )    (-s  c) (h(n - j, n - j - 1))

c        which can be expressed as

c           (0  h(n - j, n - j)) :=

c               (h(n - j, n - j - 1)  h(n - j, n - j))*( c  s).
c                                                            (-s  c)

c        thus return  c(j) = c  and  s(j) = -s  to make the plane
c        rotation matrix look like

c           r(j) = ( c(j)  s(j)).(-s(j)  c(j))

         do 80 j = k2 - 1, k1, -1
            suph = s(j)
            call srotg1(a(n-j,j+1),suph,ctemp,stemp)
            stemp = -stemp
            s(j) = stemp
            c(j) = ctemp
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               do 60 i = n - j + 1, n
                  temp = a(i,j+1)
                  a(i,j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   60          continue
            end if
   80    continue
      end if

      end

      subroutine scond (n, x, incx, xmax, xmin)
c----------------------------------------------------------------------
c blas routine scond
c----------------------------------------------------------------------
      implicit none

      double precision xmax, xmin, x(*)
      integer incx, n, ix
c----------------------------------------------------------------------
      if (n.lt.1) then

         xmax = 0d0
         xmin = 0d0

      else

         xmax = dabs(x(1))
         xmin = xmax

         do ix = 1 + incx, 1 + (n - 1)*incx, incx
            xmax = max(xmax,dabs(x(ix)))
            xmin = min(xmin,dabs(x(ix)))
         end do

      end if

      end

      subroutine suhqr (side, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
c blas routine suhqr
c----------------------------------------------------------------------

      implicit none

      integer k1, k2, lda, n, i, j

      character*1 side

      double precision a(lda,*), c(*), s(*), aij, ctemp, stemp, subh, 
     *                 temp

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if ((min(n, k1).lt.1).or.(k2.le.k1).or.(k2.gt.n))return

      if (side.eq.'l') then

         do 20 j = k1, n
            aij = a(k1,j)
            do 10 i = k1, min(j, k2) - 1
               temp = a(i + 1,j)
               a(i,j) = s(i)*temp + c(i)*aij
               aij = c(i)*temp - s(i)*aij
   10       continue
            if (j.lt.k2) then

               subh = s(j)
               call srotg1(aij, subh, c(j), s(j))
               a(j,j) = aij
            else
               a(k2,j) = aij
            end if
   20    continue

      else if (side.eq.'r') then

         do 40 j = k2 - 1, k1, -1

            subh = s(j)

            call srotg1(a(j + 1, j + 1), subh, ctemp, stemp)

            stemp = -stemp

            if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

            s(j) = stemp
            c(j) = ctemp

            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then

               do 30 i = j, 1, -1
                  temp = a(i, j + 1)
                  a(i, j + 1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   30          continue

            end if

   40    continue

      end if

      end

      subroutine sutsrh (side, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
c blas routine
c----------------------------------------------------------------------
      implicit none

      integer k1, k2, lda, n, i, j

      character side*1

      double precision a(lda, *), c(*), s(*), aij, ctemp, stemp, temp

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if (min(n,k1).lt.1.or.k2.le.k1.or.k2.gt.n) return

      if (side.eq.'l') then

         do 20 j = n, k1, -1

            if (j.ge.k2) then
               aij = a(k2,j)
            else
               aij = c(j)*a(j,j)
               s(j) = -s(j)*a(j,j)
            end if

            do 10 i = min(k2,j) - 1, k1, -1
               temp = a(i,j)
               a(i + 1,j) = c(i)*aij - s(i)*temp
               aij = s(i)*aij + c(i)*temp
   10       continue
            a(k1,j) = aij
   20    continue

      else if (side.eq.'r') then

         do j = k1, k2 - 1

            if (c(j).ne.1d0.or.s(j).ne.0d0) then
               stemp = s(j)
               ctemp = c(j)

               if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

               do i = 1, j

                  temp = a(i, j + 1)
c                                        could check for underflow.
c                                 added following line 11/06
                  if (dabs(a(i,j)).lt.wmach(3)) a(i,j) = 0d0                            
                  a(i, j + 1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)

               end do

               s(j) = stemp*a(j + 1, j + 1)
               a(j + 1, j + 1) = ctemp*a(j + 1, j + 1)

            end if
         end do 
      end if

      end

      double precision function snorm (scale, ssq)
c----------------------------------------------------------------------
c original code was using

c flmax = wmach(7) = huge
c flmin = 1/flmax ~= wmach(3) = epsmch
c----------------------------------------------------------------------
      implicit none

      double precision scale, ssq, sqt, norm

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      sqt = dsqrt(ssq)

      if (scale.lt.wmach(7)/sqt) then
         norm = scale*sqt
      else
         norm = wmach(7)
      end if

      snorm = norm

      end

      subroutine sutsqr (side, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
c  do the transformation
c
c     r := p*u*q'  when  side = 'l' or 'l'  ( left-hand side)
c
c     r := q*u*p'  when  side = 'r' or 'r'  (right-hand side),
c
c  where  u and r  are  n by n  upper  triangular  matrices,   p  is  an
c  orthogonal matrix,  consisting of a given sequence of plane rotations
c  to be  applied  in  planes  k1 to k2,  and  q  is  a  unitary  matrix
c  consisting of a sequence of plane rotations, applied in planes  k1 to
c  k2,  chosen to make  r  upper triangular.
c
c  when  side = 'l' or 'l'  then  p  is  given  as a  sequence of  plane
c  rotation matrices
c
c     p = p(k2 - 1)*...*p(k1 + 1)*p(k1),
c
c  where  p(k) is a plane rotation matrix for the  (k, k + 1) plane.
c  in this case the matrix q is given as
c
c     q = q(k2 - 1)*...*q(k1 + 1)*q(k1),
c
c  where  q(k) is a plane rotation matrix for the  (k, k + 1) plane.
c
c  when  side = 'r' or 'r'  then  p  is  given  as a  sequence of  plane
c  rotation matrices
c
c     p = p(k1)*p(k1 + 1)*...*p(k2 - 1),
c
c  where  p(k) is a plane rotation matrix for the  (k, k + 1) plane.
c  in this case the matrix q is given as
c
c     q = q(k1)*q(k1 + 1)*...*q(k2 - 1),
c
c  where  q(k) is a plane rotation matrix for the  (k, k + 1) plane.
c
c  the  upper  triangular  matrix  u  must  be  supplied  in the  n by n
c  leading upper triangular part of  a,  and this  is overwritten by the
c  upper triangular matrix  r.  the cosine  and  sine  that  define  the
c  plane rotation matrix  p(k)  must be supplied in  c(k) and s(k)
c  respectively,  and  the two by two rotation part of  p(k),  t(k),
c  is assumed to be of the form
c
c     t(k) = ( c(k)  s(k)).
c              (-s(k)  c(k))
c
c  the cosine  and  sine that define  q(k)  are overwritten on  c(k)
c  and  s(k)  respectively and the two by two rotation part of  q(k)
c  will have the form of  t(k)  above.
c
c  if  n or k1  are less  than  unity, or  k1  is not  less than  k2, or
c  k2  is greater than  n  then an immediate return is effected.
c----------------------------------------------------------------------
      integer k1, k2, lda, n
      character        side*1
      double precision a(lda, *), c(*), s(*)
      double precision aij, ctemp, fill, stemp, temp
      integer i, i1, j
c----------------------------------------------------------------------
      if ((min(n, k1).lt.1).or.(k2.le.k1).or.
     *   (k2.gt.n))return
      if ((side.eq.'l').or.(side.eq.'L')) then
c
c        apply the left-hand transformations,  column by column,  to the
c        triangular part of  u,  but not to  anywhere  that would  cause
c        fill.
c
         do 20 j = k1 + 1, n
c
c           apply  p(k1) ... p(j - 1)  to column j.
c
            aij = a(k1, j)
            do 10 i = k1, min(j - 1, k2 - 1)
               a(i, j) = s(i)*a(i + 1, j) + c(i)*aij
               aij = c(i)*a(i + 1, j) - s(i)*aij
   10       continue
            a(i, j) = aij
   20    continue
c
c           now apply each  left-hand tranformation  to form the fill-in
c           elements and apply a  right-hand transformation to eliminate
c           the fill-in element.
c
         do 40 j = k1, k2 - 1
c
c           apply  p(j)  to the jth diagonal element  and the  fill-in
c           position.
c
            fill = -s(j)*a(j, j)
            a(j, j) = c(j)*a(j, j)
c
c           now  set up  the rotation  q(j)  to eliminate the  fill-in
c           element,  and  apply  q(j)  to  the  jth  and  (j + 1)th
c           columns.
c
            call srotg1(a(j + 1, j + 1), fill, ctemp, stemp)
            c(j) = ctemp
            s(j) = -stemp
            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then
               stemp = -stemp
               do 30 i = 1, j
                  temp = a(i, j + 1)
                  a(i, j + 1) = ctemp*temp - stemp*a(i, j)
                  a(i, j) = stemp*temp + ctemp*a(i, j)
   30          continue
            end if
   40    continue
      else if ((side.eq.'r').or.(side.eq.'R')) then
c
c        intermingle the  left and right hand transformations so that
c        at the kth step form
c
c           a := q(k)*a*p(k)'.
c
c        first  apply  the  transformations  in  columns  k2 back to k1.
c
         do 60 j = k2 - 1, k1, -1
c
c           first apply  p(j).
c
            if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
               ctemp = c(j)
               stemp = s(j)
               do 50 i = 1, j
                  temp = a(i, j + 1)
                  a(i, j + 1) = ctemp*temp - stemp*a(i, j)
                  a(i, j) = stemp*temp + ctemp*a(i, j)
   50          continue
c
c              next form the fill-in element  a(j + 1, j)  by applying
c              p(j).
c
               fill = s(j)*a(j + 1, j + 1)
               a(j + 1, j + 1) = c(j)*a(j + 1, j + 1)
c
c              now set up the rotation  q(j)  to eliminate the fill-in
c              element.
c
               call srotg1(a(j, j), fill, c(j), s(j))
            end if
   60    continue
c
c        finally  apply  q(k2 - 1) ... q(k1)  to columns  n  back to
c        (k1 + 1).
c
         do 80 j = n, k1 + 1, -1
            i1 = min(k2, j)
            aij = a(i1, j)
            do 70 i = i1 - 1, k1, -1
               temp = a(i, j)
               a(i + 1, j) = c(i)*aij - s(i)*temp
               aij = s(i)*aij + c(i)*temp
   70       continue
            a(k1, j) = aij
   80    continue
      end if

      end

      subroutine sdscl (n, d, incd, x, incx)
c-----------------------------------------------------------------------
c  sdscl does the operation
c
c     x := diag(d)*x
c-----------------------------------------------------------------------
      implicit none

      integer incd, incx, n
      double precision d(*), x(*)
      integer i, id, ix

      double precision wmach
      common/ cstmch /wmach(10)
c-----------------------------------------------------------------------
      if (n.gt.0) then
         if (incd.eq.0.and.incx.ne.0) then
            call dscal(n, d(1), x, abs(incx))
         else if (incd.eq.incx.and.incd.gt.0) then
            do 10, id = 1, 1 + (n - 1)*incd, incd
               x(id) = d(id)*x(id)
   10       continue
         else
            if (incx.ge.0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incd.gt.0) then
               do 20, id = 1, 1 + (n - 1)*incd, incd
                  x(ix) = d(id)*x(ix)
                  ix      = ix              + incx
   20          continue
            else
               id = 1 - (n - 1)*incd
               do 30, i = 1, n
                  x(ix) = d(id)*x(ix)
                  id      = id              + incd
                  ix      = ix              + incx
   30          continue
            end if
         end if
      end if

      end

      subroutine sgeqrp (pivot,m,n,a,lda,zeta,perm,work,ifail)

c  sgeqrp  finds  a  qr factorization  of  the  real  m by n  matrix  a,
c  incorporating  column interchanges,  so that  a  is reduced to  upper
c  triangular form  by means of  orthogonal transformations  and  column
c  permutations.
c
c  2. description
c     ===========
c
c  the m by n matrix a is factorized as
c
c     a = q*(r)*p'      when   m.gt.n,
c           (0)
c
c     a = q*r*p'          when   m = n,
c
c     a = q*(r  x)*p'   when   m.lt.n,
c
c  where  q  is an  m by m  orthogonal matrix, r  is  a  min(m, n)  by
c  min(m, n)  upper triangular matrix and  p is an  n by n permutation
c  matrix.
c
c  the  factorization  is  obtained  by  householder's  method. the  kth
c  transformation matrix, q(k), which is used  to introduce zeros into
c  the kth column of a is given in the form
c
c     q(k) = (i     0  ),
c              (0  t(k))
c
c  where
c
c     t(k) = i - u(k)*u(k)',
c
c     u(k) = (zeta(k)),
c              (   z(k))
c
c  zeta(k)  is a scalar and  z(k)  is an  (m - k)  element vector.
c  zeta(k) and z(k)  are chosen to annhilate the elements  below the
c  triangular part of  a.
c
c  the vector  u(k) is returned in the kth element of  zeta and in the
c  kth column of a, such that zeta(k) is in zeta(k) and the elements
c  of  z(k) are in  a(k + 1, k), ..., a(m, k).  the elements of  r
c  are returned in the upper triangular part of  a.
c
c  q is given by
c
c     q = (q(p)*q(p - 1)*...*q(1))',   p = min(m, n).
c
c  two options are available for the column permutations. in either case
c  the column for which the  sub-diagonal elements are to be annihilated
c  at the  kth step is chosen from the remaining (n - k + 1)  columns.
c  the  particular column chosen as the pivot column is either that  for
c  which  the  unreduced  part  (elements k onwards)  has the  largest
c  euclidean  length, or  is that for  which the ratio of the  euclidean
c  length  of the  unreduced part  to the  euclidean length of the whole
c  column is a maximum.

c  pivot  - character*1.
c
c           on  entry, pivot  specifies  the  pivoting  strategy  to  be
c           performed as follows.
c
c           pivot = 'c' or 'c'   (column interchanges)
c
c              column  interchanges  are  to be  incorporated  into  the
c              factorization, such that the  column whose unreduced part
c              has  maximum  euclidean  length  is chosen  as the  pivot
c              column at each step.
c
c           pivot = 's' or 's'   (scaled column interchanges)
c
c              scaled  column interchanges  are to be  incorporated into
c              the  factorization, such  that the  column for which  the
c              ratio  of the  euclidean  length of the unreduced part of
c              the column to the original euclidean length of the column
c              is a maximum is chosen as the  pivot column at each step.
c
c           unchanged on exit.
c
c  m      - integer.
c
c           on entry, m  must specify the number of rows of a. m must be
c           at  least  zero. when  m = 0  then  an  immediate return  is
c           effected.
c
c           unchanged on exit.
c
c  n      - integer.
c
c           on entry, n  must specify the number of columns of a. n must
c           be  at least zero. when  n = 0  then an immediate return  is
c           effected.
c
c           unchanged on exit.
c
c  a      - real             array of dimension (lda, n).
c
c           before entry, the leading  m by n  part of the array  a must
c           contain the matrix to be factorized.
c
c           on  exit, the  min(m, n) by min(m, n)  upper  triangular
c           part of a will contain the upper triangular matrix r and the
c           m by min(m, n)  strictly lower triangular part of  a  will
c           contain details  of the  factorization  as  described above.
c           when m.lt.n then the remaining m by (n - m) part of a will
c           contain the matrix x.
c
c  lda    - integer.
c
c           on entry, lda  must  specify  the  leading dimension of  the
c           array  a  as declared in the calling (sub) program. lda must
c           be at least  m.
c
c           unchanged on exit.
c
c  zeta   - real             array of dimension at least (n).
c
c           on exit,  zeta(k)  contains the scalar  zeta(k)  for the
c           kth  transformation.  if  t(k) = i  then  zeta(k) = 0.0,
c           otherwise  zeta(k)  contains  zeta(k) as described above
c           and  zeta(k) is always in the range  (1.0, sqrt(2.0)).
c           when n.gt.m the elements  zeta(m + 1), zeta(m + 2), ...,
c           zeta(n)  are used as internal workspace.
c
c  perm   - integer array of dimension at least  min(m, n).
c
c           on exit, perm  contains details of the permutation matrix p,
c           such  that  perm(k) = k  if no  column interchange occured
c           at  the  kth  step  and  perm(k) = j, (k .lt. j .le. n),
c           if columns  k  and  j  were  interchanged  at the  kth step.
c           note that there are  min(m, n) permutations.
c
c  work   - real array of dimension at least (2*n).
c
c           used as internal workspace.
c
c           on exit, work(j), j = 1, 2, ..., n, contains the euclidean
c           length of the jth column of the permuted matrix a*p'.
c
c  ifail  - integer.
c
c           before entry,  ifail  must contain one of the values -1 or 0
c           or 1 to specify noisy soft failure or noisy hard failure  or
c           silent soft failure. (see chapter p01 for further details.)
c
c           on  successful exit, ifail  will be  zero,  otherwise  ifail
c           will  be set to   -1  indicating that an input parameter has
c           been  incorrectly supplied. see the next section for further
c           details.
c
c  4. diagnostic information
c     ======================
c
c  ifail = -1
c
c     one or more of the following conditions holds:
c
c        pivot .ne. 'c' or 'c' or 's' or 's'
c        m     .lt. 0
c        n     .lt. 0
c        lda   .lt. m
c
c  if  on  entry,  ifail  was either  -1 or 0  then  further  diagnostic
c  information  will  be  output  on  the  error message  channel.
      implicit none

      double precision lamda
      parameter (lamda=1.0d-2)
      character*6       srname
      parameter (srname='sgeqrp')
      integer ifail, lda, m, n
      character*1 pivot
      double precision a(lda,*), work(*), zeta(*)
      integer perm(*)
      double precision eps, maxnrm, norm, temp, tol
      integer j, jmax, k, la
      double precision dnrm2
      external          dnrm2

      double precision wmach
      common/ cstmch /wmach(10)

c     compute eps and the initial column norms.

      if (min(m,n).eq.0) then
         ifail = 0
         return
      end if
c                                this should be (1/2)*b**(1-p)
c                                or maybe b**(1-p)
      eps = wmach(8)

      do 20 j = 1, n
         work(j) = dnrm2(m,a(1,j),1)
         work(j+n) = work(j)
   20 continue
c
c     do the factorization. tol is the tolerance for f06frf.
c
      la = lda
      do 120 k = 1, min(m,n)
c
c        find the pivot column.
c
         maxnrm = 0d0
         jmax = k
         if (pivot.eq.'c') then
            do 40 j = k, n
               if (work(j+n).gt.maxnrm) then
                  maxnrm = work(j+n)
                  jmax = j
               end if
   40       continue
         else
            do 60 j = k, n
               if (work(j).gt.0d0) then
                  if (k.le.1) then
                     jmax = j
                     go to 80
                  else if ((work(j+n)/work(j)).gt.maxnrm) then
                     maxnrm = work(j+n)/work(j)
                     jmax = j
                  end if
               end if
   60       continue
   80       continue
         end if
         perm(k) = jmax
         if (jmax.gt.k) then
            call dswap(m,a(1,k),1,a(1,jmax),1)
            temp = work(k)
            work(k) = work(jmax)
            work(jmax) = temp
            work(jmax+n) = work(k+n)
         end if
         tol = eps*work(k)
         if (k.lt.m) then
c
c           use a householder reflection to zero the kth column of a.
c           first set up the reflection.
c
            call sgrfg (m-k,a(k,k),a(k+1,k),1,tol,zeta(k))

            if (k.lt.n) then
               if (zeta(k).gt.0d0) then
                  if ((k+1).eq.n) la = m - k + 1
c
c                 temporarily store beta and put zeta(k) in a(k, k).
c
                  temp = a(k,k)
                  a(k,k) = zeta(k)
c
c                 now do the operation  a := q(k)*a.
c
c                 let  b  denote  the bottom  (m - k + 1) by (n - k)
c                 part of a.
c
c                 first  form   work = b'*u.  (work  is  stored  in the
c                 elements zeta(k + 1), ..., zeta(n).)
c
                  call dgemv ('t',m-k+1,n-k,1d0,a(k,k+1),la,
     *                       a(k,k),0d0,zeta(k+1))
c
c                 now form  b := b - u*work'.
c
                  call dger1 (m-k+1,n-k,a(k,k),zeta(k+1),a(k,k+1),la)

c                 restore beta.
c
                  a(k,k) = temp
               end if
c
c              update  the  unreduced  column  norms.  use  the  linpack
c              criterion for when to recompute the norms, except that we
c              retain  the original column lengths throughout  and use a
c              smaller lamda.
c
               do 100 j = k + 1, n
                  if (work(j+n).gt.0d0) then
                     temp = dabs(a(k,j))/work(j+n)
                     temp = max((1d0+temp)*(1d0-temp),0d0)
                     norm = temp
                     temp = 1d0 + lamda*temp*(work(j+n)/work(j))**2
                     if (temp.gt.1d0) then
                        work(j+n) = work(j+n)*dsqrt(norm)
                     else
                        work(j+n) = dnrm2(m-k,a(k+1,j),1)
                     end if
                  end if
  100          continue
            end if
         end if
  120 continue
c
c     set the final  zeta  when  m.le.n.
c
      if (m.le.n) zeta(m) = 0d0
c
      ifail = 0

99999 format ('    the input parameters contained ',i2,' error(s)')
      end

      subroutine sutsr1 (side,n,k1,k2,s,a,lda)

c  sutsrh applies a  sequence  of  pairwise interchanges to either  the
c  left,  or the right,  of the  n by n  upper triangular matrix  u,  to
c  transform u to an  upper hessenberg matrix. the interchanges are
c  applied in planes k1 up to k2.
c
c  the upper hessenberg matrix, h, is formed as
c
c     h = p*u,    when   side = 'l' or 'l',  ( left-hand side)
c
c  where p is a permutation matrix of the form
c
c     p = p(k1)*p(k1 + 1)*...*p(k2 - 1)
c
c  and is formed as
c
c     h = u*p',   when   side = 'r' or 'r',  (right-hand side)
c
c  where p is a permutation matrix of the form
c
c     p = p(k2 - 1)*...*p(k1 + 1)*p(k1),
c
c  p(k) being a pairwise interchange for the  (k, k + 1) plane.
c  the  two by two
c  interchange part of p(k), r(k), is assumed to have the form
c
c     r(k) = (0  1).
c              (1  0)
c
c  the matrix  u must be supplied in the n by n leading upper triangular
c  part of the array  a, and this is overwritten by the upper triangular
c  part of  h.
c
c  the  sub-diagonal elements of  h, h(k + 1, k),  are returned in the
c  elements s(k),  k = k1, k1 + 1, ..., k2 - 1.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      implicit none

      integer k1, k2, lda, n
      character*1 side

      double precision a(lda,*), s(*)

      double precision aij, temp
      integer i, j

      if ((min(n,k1).lt.1).or.(k2.le.k1).or.(k2.gt.n)) return
      if (side.eq.'l') then
c
c        apply the permutations to columns n back to k1.
c
         do 40 j = n, k1, -1
            if (j.ge.k2) then
               aij = a(k2,j)
            else
c
c              form  the  additional sub-diagonal element  h(j + 1, j)
c              and store it in s(j).
c
               aij = 0d0
               s(j) = a(j,j)
            end if
            do 20 i = min(k2,j) - 1, k1, -1
               temp = a(i,j)
               a(i+1,j) = temp
               aij = aij
   20       continue
            a(k1,j) = aij
   40    continue
      else if (side.eq.'r') then
c
c        apply  the  plane interchanges to  columns  k1  up to
c        (k2 - 1) and  form   the   additional  sub-diagonal
c        elements,   storing  h(j + 1, j) in s(j).
c
         do 80 j = k1, k2 - 1
            do 60 i = 1, j
               temp = a(i,j+1)
               a(i,j+1) = a(i,j)
               a(i,j) = temp
   60       continue
            s(j) = a(j+1,j+1)
            a(j+1,j+1) = 0d0
   80    continue
      end if

      end

      subroutine smcopy (matrix, m, n, a, lda, b, ldb)
c  smcopy  copies  the  m by n  matrix  a  into  the  m by n  matrix  b.
c
c  if   matrix = 'g'    then  a  and  b  are  regarded as  general
c                             matrices,
c  if   matrix = 'u'    then  a  and  b  are  regarded  as   upper
c                             triangular,  and only  elements  for which
c                             i.le.j  are referenced,
c  if   matrix = 'l'    then  a  and  b  are  regarded  as   lower
c                             triangular,  and only  elements  for which
c                             i.ge.j  are referenced.

      implicit none

      character*1 matrix
      integer m, n, lda, ldb
      double precision a(lda, *), b(ldb, *)
      integer i, j


      if (matrix.eq.'g') then
         do 20 j = 1, n
            do 10 i = 1, m
               b(i, j) = a(i, j)
   10       continue
   20    continue
      else if (matrix.eq.'u') then
         do 40 j = 1, n
            do 30 i = 1, min(m, j)
               b(i, j) = a(i, j)
   30       continue
   40    continue
      else if (matrix.eq.'l') then
         do 60 j = 1, min(m, n)
            do 50 i = j, m
               b(i, j) = a(i, j)
   50       continue
   60    continue
      end if

      end

      subroutine icopy (n, x, incx, y, incy)

c  icopy does the operation
c     y := x

      implicit none

      integer incx, incy, n
      integer x(*), y(*)
      integer i, ix, iy

      if (n.gt.0) then
         if ((incx.eq.incy).and.(incy.gt.0)) then
            do 10, iy = 1, 1 + (n - 1)*incy, incy
               y(iy) = x(iy)
   10       continue
         else
            if (incx.ge.0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incy.gt.0) then
               do 20, iy = 1, 1 + (n - 1)*incy, incy
                  y(iy) = x(ix)
                  ix      = ix      + incx
   20          continue
            else
               iy = 1 - (n - 1)*incy
               do 30, i = 1, n
                  y(iy) = x(ix)
                  iy      = iy      + incy
                  ix      = ix      + incx
   30          continue
            end if
         end if
      end if

      end

      subroutine sgeqr(m,n,a,lda,zeta,ifail)

c  sgeqr  finds  the  qr factorization  of the real  m by n,  m .ge. n,
c  matrix a,  so that  a is reduced to upper triangular form by means of
c  orthogonal transformations.

c  the m by n matrix a is factorized as
c
c     a = q*(r)   when   m.gt.n,
c           (0)
c
c     a = q*r       when   m = n,
c
c  where  q  is an  m by m orthogonal matrix and  r  is an  n by n upper
c  triangular matrix.
c
c  the  factorization  is  obtained  by  householder's  method. the  kth
c  transformation matrix, q(k), which is used  to introduce zeros into
c  the kth column of a is given in the form
c
c     q(k) = (i     0  ),
c              (0  t(k))
c
c  where
c
c     t(k) = i - u(k)*u(k)',
c
c     u(k) = (zeta(k)),
c              (   z(k))
c
c  zeta(k)  is a scalar and  z(k)  is an  (m - k)  element vector.
c  zeta(k) and z(k)  are chosen to annhilate the elements  below the
c  triangular part of  a.
c
c  the vector  u(k) is returned in the kth element of  zeta and in the
c  kth column of a, such that zeta(k) is in zeta(k) and the elements
c  of  z(k) are in  a(k + 1, k), ..., a(m, k).  the elements of  r
c  are returned in the upper triangular part of  a.
c
c  q is given by
c
c     q = (q(n)*q(n - 1)*...*q(1))'.

c  m      - integer.
c
c           on entry, m must specify the number of rows of  a. m must be
c           at least  n.
c
c           unchanged on exit.
c
c  n      - integer.
c
c           on entry, n must specify the number of columns of  a. n must
c           be  at  least zero. when  n = 0  then an immediate return is
c           effected.
c
c           unchanged on exit.
c
c  a      - real             array of dimension (lda, n).
c
c           before entry, the leading  m by n  part of the array  a must
c           contain the matrix to be factorized.
c
c           on exit, the  n by n upper triangular part of a will contain
c           the upper triangular matrix r and the  m by n strictly lower
c           triangular  part   of   a   will  contain  details   of  the
c           factorization as described above.
c
c  lda    - integer.
c
c           on entry, lda  must  specify  the  leading dimension of  the
c           array  a  as declared in the calling (sub) program. lda must
c           be at least  m.
c
c           unchanged on exit.
c
c  zeta   - real             array of dimension at least (n).
c
c           on exit,  zeta(k)  contains the scalar  zeta(k)  for the
c           kth  transformation.  if  t(k) = i  then  zeta(k) = 0.0,
c           otherwise  zeta(k)  contains  zeta(k) as described above
c           and  zeta(k) is always in the range  (1.0, sqrt(2.0)).
c
c  ifail  - integer.
c
c           before entry,  ifail  must contain one of the values -1 or 0
c           or 1 to specify noisy soft failure or noisy hard failure  or
c           silent soft failure. (see chapter p01 for further details.)
c
c           on successful  exit  ifail  will be  zero,  otherwise  ifail
c           will  be set to  -1  indicating that an  input parameter has
c           been  incorrectly  set. see  the  next section  for  further
c           details.

c  ifail = -1
c
c     one or more of the following conditions holds:
c
c        m   .lt. n
c        n   .lt. 0
c        lda .lt. m
c
c  if  on  entry,  ifail  was  either  -1 or 0  then  further diagnostic
c  information  will  be  output  on  the  error message  channel.

      character*6 srname
      parameter (srname='sgeqr ')

      integer ifail, lda, m, n

      double precision a(lda,*), zeta(*)

      double precision temp
      integer k, la

c     do the factorization.

      if (n.eq.0) then
         ifail = 0
         return
      end if
      la = lda
      do 20 k = 1, min(m-1,n)

c        use a  householder reflection  to  zero the  kth column  of  a.
c        first set up the reflection.

         call sgrfg (m-k,a(k,k),a(k+1,k),1,0d0,zeta(k))
         if ((zeta(k).gt.0d0).and.(k.lt.n)) then
            if ((k+1).eq.n) la = m - k + 1

c           temporarily  store  beta and  put  zeta(k)  in  a(k, k).

            temp = a(k,k)
            a(k,k) = zeta(k)

c           now do the operation  a := q(k)*a.

c           let  b  denote  the bottom  (m - k + 1) by (n - k)  part
c           of  a.

c           first form   work = b'*u.  (work  is stored in the elements
c           zeta(k + 1), ..., zeta(n).)

            call dgemv('t',m-k+1,n-k,1d0,a(k,k+1),la,a(k,k),
     *                 0d0,zeta(k+1))

c           now form  b := b - u*work'.

            call dger1 (m-k+1,n-k,a(k,k),zeta(k+1),a(k,k+1),la)

c           restore beta.

            a(k,k) = temp
         end if
   20 continue

c     set the final  zeta  when  m.eq.n.

      if (m.eq.n) zeta(n) = 0d0

      ifail = 0

99999 format ('    the input parameters contained ',i2,' error(s)')
      end

      subroutine iload (n, const, x, incx)

c  iload does the operation

c     x = const*e,   e' = (1  1 ... 1).

      integer const, incx, n
      integer x(*)
      integer ix

      if (n.gt.0) then
         if (const.ne.0) then
            do 10, ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = const
   10       continue
         else
            do 20, ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = 0
   20       continue
         end if
      end if

      end

      subroutine smlod1 (m, n, a, lda)
c----------------------------------------------------------------------
c modified blas routine smload (smload)
c----------------------------------------------------------------------
      implicit none

      integer i, j, lda, m, n

      double precision a(lda, *)
c----------------------------------------------------------------------
         do j = 1, n
            do i = 1, m
               a(i,j) = 0d0
            end do
         end do

         do i = 1, min(m,n)
            a(i,i) = 1d0
         end do

      end

      subroutine smload (matrix, m, n, const, diag, a, lda)
c----------------------------------------------------------------------
c  smload forms the m by n matrix a given by

c     a(i, j) = ( diag  i.eq.j,
c                 (
c                 (const  i.ne.j.

c  if   matrix = 'g' or 'g'   then  a  is regarded  as a general matrix,
c  if   matrix = 'u' or 'u'   then  a  is regarded  as upper triangular,
c                             and only  elements  for which  i.le.j  are
c                             referenced,
c  if   matrix = 'l' or 'l'   then  a  is regarded  as lower triangular,
c                             and only  elements  for which  i.ge.j  are
c                             referenced.
c----------------------------------------------------------------------
      character*1 matrix
      double precision const, diag
      integer lda, m, n
      double precision a(lda, *)
      integer i, j

      if (matrix.eq.'g') then
         do 20 j = 1, n
            do 10 i = 1, m
               a(i, j) = const
   10       continue
   20    continue
         if (const.ne.diag) then
            do 30 i = 1, min(m, n)
               a(i, i) = diag
   30       continue
         end if
      else if (matrix.eq.'u') then
         do 50 j = 1, n
            do 40 i = 1, min(m, j)
               a(i, j) = const
   40       continue
   50    continue
         if (const.ne.diag) then
            do 60 i = 1, min(m, n)
               a(i, i) = diag
   60       continue
         end if
      else if (matrix.eq.'l') then
         do 80 j = 1, min(m, n)
            do 70 i = j, m
               a(i, j) = const
   70       continue
   80    continue
         if (const.ne.diag) then
            do 90 i = 1, min(m, n)
               a(i, i) = diag
   90       continue
         end if
      end if

      end

      double precision function dlantr(norm,uplo,diag,m,n,a,lda,work)

c  dlantr  returns the value of the one norm,  or the frobenius norm, or
c  the  infinity norm,  or the  element of  largest absolute value  of a
c  trapezoidal or triangular matrix a.

c  dlantr returns the value
c
c     dlantr = (max(abs(a(i,j))), norm = 'm' or 'm'
c              (
c              (norm1(a),         norm = '1', 'o' or 'o'
c              (
c              (normi(a),         norm = 'i' or 'i'
c              (
c              (normf(a),         norm = 'f', 'f', 'e' or 'e'
c
c  where  norm1  denotes the  one norm of a matrix (maximum column sum),
c  normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
c  normf  denotes the  frobenius norm of a matrix (square root of sum of
c  squares).  note that  max(abs(a(i,j)))  is not a  matrix norm.

c  norm    (input) character*1
c          specifies the value to be returned in dlantr as described
c          above.
c
c  uplo    (input) character*1
c          specifies whether the matrix a is upper or lower trapezoidal.
c          = 'u':  upper trapezoidal
c          = 'l':  lower trapezoidal
c          note that a is triangular instead of trapezoidal if m = n.
c
c  diag    (input) character*1
c          specifies whether or not the matrix a has unit diagonal.
c          = 'n':  non-unit diagonal
c          = 'u':  unit diagonal
c
c  m       (input) integer
c          the number of rows of the matrix a.  m >= 0, and if
c          uplo = 'u', m <= n.  when m = 0, dlantr is set to zero.
c
c  n       (input) integer
c          the number of columns of the matrix a.  n >= 0, and if
c          uplo = 'l', n <= m.  when n = 0, dlantr is set to zero.
c
c  a       (input) real array, dimension (lda,n)
c          the trapezoidal matrix a (a is triangular if m = n).
c          if uplo = 'u', the leading m by n upper trapezoidal part of
c          the array a contains the upper trapezoidal matrix, and the
c          strictly lower triangular part of a is not referenced.
c          if uplo = 'l', the leading m by n lower trapezoidal part of
c          the array a contains the lower trapezoidal matrix, and the
c          strictly upper triangular part of a is not referenced.  note
c          that when diag = 'u', the diagonal elements of a are not
c          referenced and are assumed to be one.
c
c  lda     (input) integer
c          the leading dimension of the array a.  lda >= max(m,1).
c
c  work    (workspace) real array, dimension (lwork),
c          where lwork >= m when norm = 'i'; otherwise, work is not
c          referenced.
      implicit none 

      integer             lda, m, n
      character                        diag, norm, uplo
      double precision a(lda,*), work(*)
      double precision scale, sum, value
      integer             i, j
      logical                          udiag

      if (min(m,n).eq.0) then
         value = 0d0
      else if (norm.eq.'m') then
c
c        find max(abs(a(i,j))).
c
         if (diag.eq.'u') then
            value = 1d0
            if (uplo.eq.'u') then
               do 40 j = 1, n
                  do 20 i = 1, min(m,j-1)
                     value = max(value,dabs(a(i,j)))
   20             continue
   40          continue
            else
               do 80 j = 1, n
                  do 60 i = j + 1, m
                     value = max(value,dabs(a(i,j)))
   60             continue
   80          continue
            end if
         else
            value = 0d0
            if (uplo.eq.'u') then
               do 120 j = 1, n
                  do 100 i = 1, min(m,j)
                     value = max(value,abs(a(i,j)))
  100             continue
  120          continue
            else
               do 160 j = 1, n
                  do 140 i = j, m
                     value = max(value,abs(a(i,j)))
  140             continue
  160          continue
            end if
         end if
      else if (norm.eq.'o'.or.norm.eq.'1') then
c
c        find norm1(a).
c
         value = 0d0

         udiag = (diag.eq.'u')

         if (uplo.eq.'u') then
            do 220 j = 1, n
               if (udiag.and.(j.le.m)) then
                  sum = 1d0
                  do 180 i = 1, j - 1
                     sum = sum + dabs(a(i,j))
  180             continue
               else
                  sum = 0d0
                  do 200 i = 1, min(m,j)
                     sum = sum + dabs(a(i,j))
  200             continue
               end if
               value = max(value,sum)
  220       continue
         else
            do 280 j = 1, n
               if (udiag) then
                  sum = 1d0
                  do 240 i = j + 1, m
                     sum = sum + dabs(a(i,j))
  240             continue
               else
                  sum = 0d0
                  do 260 i = j, m
                     sum = sum + dabs(a(i,j))
  260             continue
               end if
               value = max(value,sum)
  280       continue
         end if
      else if (norm.eq.'i') then
c
c        find normi(a).
c
         if (uplo.eq.'u') then
            if (diag.eq.'u') then
               do 300 i = 1, m
                  work(i) = 1d0
  300          continue
               do 340 j = 1, n
                  do 320 i = 1, min(m,j-1)
                     work(i) = work(i) + dabs(a(i,j))
  320             continue
  340          continue
            else
               do 360 i = 1, m
                  work(i) = 0d0
  360          continue
               do 400 j = 1, n
                  do 380 i = 1, min(m,j)
                     work(i) = work(i) + dabs(a(i,j))
  380             continue
  400          continue
            end if
         else
            if (diag.eq.'u') then
               do 420 i = 1, n
                  work(i) = 1d0
  420          continue
               do 440 i = n + 1, m
                  work(i) = 0d0
  440          continue
               do 480 j = 1, n
                  do 460 i = j + 1, m
                     work(i) = work(i) + dabs(a(i,j))
  460             continue
  480          continue
            else
               do 500 i = 1, m
                  work(i) = 0d0
  500          continue
               do 540 j = 1, n
                  do 520 i = j, m
                     work(i) = work(i) + dabs(a(i,j))
  520             continue
  540          continue
            end if
         end if
         value = 0d0
         do 560 i = 1, m
            value = max(value,work(i))
  560    continue

      else if (norm.eq.'f'.or.norm.eq.'e') then
c
c        find normf(a).
c
         if (uplo.eq.'u') then
            if (diag.eq.'u') then
               scale = 1d0
               sum = min(m,n)
               do 580 j = 2, n
                  call sssq (min(m,j-1),a(1,j),1,scale,sum)
  580          continue
            else
               scale = 0d0
               sum = 1d0
               do 600 j = 1, n
                  call sssq (min(m,j),a(1,j),1,scale,sum)
  600          continue
            end if
         else
            if (diag.eq.'u') then
               scale = 1d0
               sum = min(m,n)
               do 620 j = 1, n
                  call sssq (m-j,a(min(m,j+1),j),1,scale,sum)
  620          continue
            else
               scale = 0d0
               sum = 1d0
               do 640 j = 1, n
                  call sssq (m-j+1,a(j,j),1,scale,sum)
  640          continue
            end if
         end if
         value = scale*dsqrt(sum)
      end if
c
      dlantr = value

      end

      subroutine cmqmul (mode,n,nz,nfree,nq,unitq,kx,v,zy,wrk)
c----------------------------------------------------------------------
c     transform the vector  v  in various ways using the
C     matrix  Q = ( Z  Y )  defined by the input parameters.
C
C        MODE               result
C        ----               ------
C
C          1                v = Z v
C          2                v = Y v
C          3                v = Q v
C
C     On input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
C     on output, v  is a full n-vector.
C
C
C          4                v = Z'v
C          5                v = Y'v
C          6                v = Q'v
C
C     On input,  v  is a full n-vector.
C     On output, v  is ordered as  ( v(free)  v(fixed) ).
C
C          7                v = Y'v
C          8                v = Q'v
C
C     On input,  v  is a full n-vector.
C     On output, v  is as in modes 5 and 6 except that v(fixed) is not
C     set.
c----------------------------------------------------------------------
      implicit none

      integer mode, n, nfree, nq, nz
      logical           unitq
      double precision v(n), wrk(n), zy(nq,*)
      integer kx(n)
      integer j, j1, j2, k, l, lenv, nfixed
c----------------------------------------------------------------------
      nfixed = n - nfree
      j1 = 1
      j2 = nfree
      if (mode.eq.1.or.mode.eq.4) j2 = nz
      if (mode.eq.2.or.mode.eq.5.or.mode.eq.7) j1 = nz + 1
      lenv = j2 - j1 + 1

      if (mode.le.3) then

         if (nfree.gt.0) call sload (nfree,0d0,wrk,1)

c        copy  v(fixed)  into the end of  wrk.

         if (mode.ge.2.and.nfixed.gt.0) call dcopy (nfixed,v(nfree+1),
     *       1,wrk(nfree+1),1)

c        set  wrk  =  relevant part of  zy * v.

         if (lenv.gt.0) then
            if (unitq) then
               call dcopy (lenv,v(j1),1,wrk(j1),1)
            else

               call dgemv ('n',nfree,j2-j1+1,1d0,zy(1,j1),nq,v(j1),1d0,
     *                    wrk)

            end if
         end if

c        expand  wrk  into  v  as a full n-vector.

         call sload (n,0d0,v,1)
         do k = 1, nfree
            j = kx(k)
            v(j) = wrk(k)
         end do 

c        copy  wrk(fixed)  into the appropriate parts of  v.

         if (mode.gt.1) then
            do 40 l = 1, nfixed
               j = kx(nfree+l)
               v(j) = wrk(nfree+l)
   40       continue
         end if

      else

c        mode = 4, 5, 6, 7  or  8.

c        put the fixed components of  v  into the end of  wrk.

         if (mode.eq.5.or.mode.eq.6) then
            do 60 l = 1, nfixed
               j = kx(nfree+l)
               wrk(nfree+l) = v(j)
   60       continue
         end if

c        put the free  components of  v  into the beginning of  wrk.

         if (nfree.gt.0) then
            do 80 k = 1, nfree
               j = kx(k)
               wrk(k) = v(j)
   80       continue

c           set  v  =  relevant part of  zy' * wrk.

            if (lenv.gt.0) then

               if (unitq) then
                  call dcopy (lenv,wrk(j1),1,v(j1),1)
               else
                  call dgemv ('t',nfree,j2-j1+1,1d0,zy(1,j1),nq,wrk,
     *                       0d0,v(j1))
               end if

            end if
         end if

c        copy the fixed components of  wrk  into the end of  v.

         if (nfixed.gt.0.and.(mode.eq.5.or.mode.eq.6))
     *       call dcopy (nfixed,wrk(nfree+1),1,v(nfree+1),1)

      end if

      end
