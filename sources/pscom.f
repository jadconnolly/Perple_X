      subroutine psytic (x0, y0, dy, tic, tic1, tic2)

      implicit none

      include 'perplex_parameters.h'

      double precision x0, y0, dy, tic, tic1, tic2, y, dy2

      integer i

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

c psytic - subroutine to draw y-axis ticks.

      y = y0 

      call psmove (x0,y0)

      if (.not.half.and..not.tenth) then 

         do while (y.lt.ymax) 
            call psrlin (tic,0d0,1d0,width)
            call psrmov (-tic,dy)
            y = y + dy 
         end do 

      else if (.not.tenth) then

         dy2 = dy / 2d0

         do while (y.lt.ymax) 

            call psrlin (tic,0d0,1d0,width)
            call psrmov (-tic,dy2)
            y = y + dy2
            if (y.ge.ymax) exit 
            call psrlin (tic1, 0d0, 1d0, width)
            call psrmov (-tic1,dy2)
            y = y + dy2
         end do 

         if (y0-dy2.gt.ymin)
     *      call psline (x0,y0-dy2,x0+tic1,y0-dy2,1d0,width)

      else

         dy2 = dy / 10d0

         do while (y.lt.ymax)

            call psrlin (tic, 0d0, 1d0, width)
            call psrmov (-tic, dy2)

            y = y + dy2

            do i = 1, 4
               if (y.ge.ymax) exit
               call psrlin (tic2, 0d0, 1d0, width)
               call psrmov (-tic2, dy2)
               y = y + dy2
            end do 

            if (y.ge.ymax) exit 

            call psrlin (tic1, 0d0, 1d0, width)
            call psrmov (-tic1, dy2)
            y = y + dy2

            do i = 1, 4
               if (y.ge.ymax) exit 
               call psrlin (tic2, 0d0, 1d0, width)
               call psrmov (-tic2, dy2)
               y = y + dy2
            end do 

         end do 

         if (y0-dy2.lt.ymin) return

         y = y0 - dy2
         call psmove (x0,y)

         do i = 1, 4
            if (y.le.ymin) return
            call psrlin (tic2, 0d0, 1d0, width)
            call psrmov (-tic2, -dy2)
            y = y - dy2
         end do 

         if (y.le.ymin) return

         call psrlin (tic1, 0d0, 1d0, width)
         call psrmov (-tic1, -dy2)

         y = y - dy2 

         do i = 1, 4
            if (y.le.ymin) exit
            call psrlin (tic2, 0d0, 1d0, width)
            call psrmov (-tic2, -dy2)
            y = y - dy2
         end do 

      end if 

      end

c------------------------------------------------------------------
      subroutine psxtic (y0, x0, dx, tic, tic1, tic2)

      implicit none 

      include 'perplex_parameters.h'

      double precision x0, y0, dx, tic, tic1, tic2, x, dx2

      integer i

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

c psxtic - subroutine to draw x-axis ticks.

      x = x0 
      call psmove (x0,y0)

      if (.not.half.and..not.tenth) then 

         do while (x.lt.xmax) 
            call psrlin (0d0, tic, 1d0, width)
            call psrmov (dx, -tic)
            x = x + dx 
         end do 

      else if (.not.tenth) then

         dx2 = dx / 2d0

         do while (x.lt.xmax) 
            call psrlin (0d0, tic, 1d0, width)
            call psrmov (dx2, -tic)
            x = x + dx2
            if (x.ge.xmax) exit
            call psrlin (0d0, tic1, 1d0, width)
            call psrmov (dx2, -tic1)
            x = x + dx2
         end do 

         if (x0-dx2.gt.xmin)
     *      call psline (x0-dx2, y0, x0-dx2, y0+tic1, 1d0, width)

      else

         dx2 = dx / 10d0

         do while (x.le.xmax) 

            call psrlin (0d0, tic, 1d0, width)
            call psrmov (dx2, -tic)
            x = x + dx2

            do i = 1, 4
               if (x.ge.xmax) exit
               call psrlin (0d0, tic2, 1d0, width)
               call psrmov (dx2, -tic2)
               x = x + dx2
            end do 

            if (x.ge.xmax) exit 
            call psrlin (0d0, tic1, 1d0, width)
            call psrmov (dx2, -tic1)
            x = x + dx2

            do i = 1, 4
               if (x.ge.xmax) exit
               call psrlin (0d0, tic2, 1d0, width)
               call psrmov (dx2, -tic2)
               x = x + dx2
            end do 

         end do 

         if (x0-dx2.lt.xmin) return

         x = x0 - dx2

         call psmove (x, y0)

         do i = 1, 4
            if (x.le.xmin) return
            call psrlin (0d0, tic2, 1d0, width)
            call psrmov (-dx2, -tic2)
            x = x - dx2
         end do

         if (x.le.xmin) return
         call psrlin (0d0, tic1, 1d0, width)
         call psrmov (dx2, -tic1) 
         x = x - dx2 

         do i = 1, 4
            if (x.le.xmin) exit
            call psrlin (0d0, tic2, 1d0, width)
            call psrmov (-dx2, -tic2)
            x = x - dx2
         end do 

      end if 

      end
 
      subroutine psylbl (y0, dy, tmin)
 
c psylbl - subroutine to put on y-axis labels.

      implicit none

      include 'perplex_parameters.h'

      double precision y0,dy,tmin,x,y,xdc,ydc

      integer i,iy1
 
      character*12 numbs(40)
 
      integer ic(40)

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      xdc = 1.17d0 * dcx * nscale 
      ydc = .667d0 * dcy * nscale
      tmin = 1d30
 
      call psnum (y0,ymax,dy,ic,iy1,numbs)
 
      y = y0
 
      do i = 1, iy1
         x = xmin - (dfloat(ic(i)+1) * xdc) 
         if (x.lt.tmin) tmin = x
         call pstext (x, y + ydc, numbs(i), ic(i))
         if (lgrid) call psline (xmin, y, xmax, y, 10d0, 0d0)
         y = y + dy
      end do 
 
      end
 
      subroutine psnum (rmin,rmax,dr,ic,i1,numbs)

      implicit none

      double precision rmin,rmax,dr,r,big,small
 
      character*1 text(12), numbs(*)*12, next(12)
 
      integer ic(*),i,k,i1,icase,int,j

      i1 = 1 + int ((rmax - rmin) / dr)

      if (rmin.gt.rmax) then 
         big = rmin
         small = rmax
      else 
         big = rmax
         small = rmin
      end if 
 
      if (big.gt.9999.9d0.and.big.lt.99999.9d0.and.small.gt.-big) then
         icase = 1
      else if (big.gt.999.9d0.and.big.le.9999.9d0
     *                        .and.small.gt.-big) then
         icase = 2
      else if (big.gt.99.9d0.and.big.le.999.9d0
     *                       .and.small.gt.-big) then
         icase = 3
      else
         icase = 4
      end if
 
      r = rmin
 
      do i = 1, i1

         if (icase.eq.2) then 
             write (numbs(i),'(i4)') int(r)
         else if (icase.eq.3) then 
            write (numbs(i),'(i3)') int(r)
         else if (icase.eq.4) then
c                                 GH, 12/23/21
            write (numbs(i),'(1pg10.3)') r
         else 
            write (numbs(i),'(i5)') int(r)
         end if 

         read (numbs(i),'(12a)') text

         k = 0

         do j = 1, 12
            if (text(j).eq.' ') cycle
            k = k + 1
            next(k) = text(j)
         end do 

         ic(i) = k
         write (numbs(i),'(12a)') (next(j),j=1,k)
         r = r + dr

      end do 
 
      end
 
      subroutine psxlbl (x0, dx)
 
c psxlbl - subroutine to put on x-axis labels.

      implicit none

      include 'perplex_parameters.h'

      integer ic(40),i,ix1

      double precision x0, dx, x, xdc, y, t
 
      character*12 numbs(40)

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      y = ymin - 1.4d0*nscale*dcy
 
      x = x0
 
      xdc = dcx*nscale/1.75d0
 
      call psnum (x0,xmax,dx,ic,ix1,numbs)
 
      do i = 1, ix1

         if (x.ne.xmin) then 
            t = x - dfloat(ic(i)) * xdc 
            call pstext (t, y, numbs(i), ic(i))
            if (lgrid) call psline (x, ymin, x, ymax, 10d0, 0d0)
         end if 

         x = x + dx

      end do 
 
      end
c----------------------------------------------------------------------
      subroutine psaxes (jop0)
 
c psaxes - subroutine to output (sloppy) axes.

      implicit none

      include 'perplex_parameters.h'

      double precision x0,y0,dx,dy,xtic,ytic,xtic1,ytic1,
     *                 xtic2,ytic2,tmin,x,y,vlo,vhi,fpoly

      integer jop0, i, nchar, ic(2), nblen

      character*8 record*64, yes*1, nums(2)*12

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      character vname*8, xname*8
      common/ csta2 /xname(k5),vname(l2)

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
c----------------------------------------------------------------------
c                                 statement function to evaluate T(P) or P(T)
      fpoly(X) = c0 + X*(c1 + X*(c2 + X*(c3 + X*(c4 + X*c5))))
c----------------------------------------------------------------------
      x0 = xmin
      y0 = ymin
      dx = xlen / 5d0
      dy = ylen / 5d0

      xtic  = xlen/45d0/xfac
      xtic1 = xtic*.67d0
      xtic2 = xtic1*.67d0
 
      ytic  = ylen/45d0
      ytic1 = ytic*.67d0
      ytic2 = ytic1*.67d0 

      if (jop0.eq.1) then

         write (*,'(/,a)') 'Modify default axes numbering (y/n)?'
         read (*,'(a)') yes

         if (yes.eq.'y'.or.yes.eq.'Y') then
            write (*,1030) 'X', x0, dx
            read (*,*) x0, dx
            write (*,1030) 'Y', y0, dy
            read (*,*) y0, dy
         end if 

      end if 
c                                       draw axes box
      call psrect (xmin,xmax,ymin,ymax,1d0,width,0)
c                                       draw left vertical axis tics
      call psytic (xmin, y0, dy, xtic, xtic1, xtic2)
c                                       draw right vertical axis
      call psytic (xmax, y0, dy, -xtic, -xtic1, -xtic2)
c                                       draw bottom horizontal axis
      call psxtic (ymin, x0, dx, ytic, ytic1, ytic2)
c                                       draw top horizontal axis
      call psxtic (ymax, x0, dx, -ytic, -ytic1, -ytic2)
 
      call pssctr (ifont,nscale,nscale,0d0)
c                                       numeric axis labels:
      call psylbl (y0, dy, tmin) 
 
      call psxlbl (x0, dx)
c                                     x-axis name
      call pssctr (ifont,nscale,nscale,0d0)  

      x = xmin + 0.5d0 * xlen - 2d0*dcx*nscale
      y = ymin - 4d0*dcy*nscale

      call pstext (x,y,vnm(1),0)
c                                       y-axis name
      call pssctr (ifont,nscale,nscale,90d0)

      x = tmin - 3.33d0*dcx*nscale
      y = ymin + 0.5d0 * ylen - 2.5d0*dcy*nscale
 
      call pstext (x,y,vnm(2),0)
c                                       sectioning constraints
      if (jvar.gt.2) then

         call pssctr (ifont,nscale,nscale,0d0)
         y = ymax + 1.2d1*dcy*nscale

         do i = 3, jvar
c                                       modfied GH, 12/23/21
            if (i.eq.3 .and. idep.gt.0) then

               vlo = fpoly(vmin(iind)) 
               vhi = fpoly(vmax(iind)) 

               call psnum (vlo,vhi,vhi-vlo,ic,nchar,nums)

               write (record,'(a)') 
     *               vnm(i)(1:nblen(vnm(i)))//
     *               ' = f('//vname(iind)(1:1)//' = '//
     *               nums(1)(1:nblen(nums(1)))//'-'//
     *               nums(2)(1:nblen(nums(2)))//')'

            else

               write (record,1000) vnm(i),vmn(i)

            end if

            nchar = nblen(record)
            call psublk (record,nchar)
            call pstext (xmin,y,record,nchar)
            y = y - 2.4*dcy*nscale
         end do 
 
      end if
 
1000  format (a,'=',1pg9.3)
1030  format (/,'Enter the starting value and interval for',
     *          ' major tick marks on',/,'the ',a,'-axis (',
     *          ' current values are:',2(1x,g9.3),')',/, 
     *          'Enter the new values:')

      end

      subroutine psaxop (jcopt,jop0,iop1)
c----------------------------------------------------------------------
c psaxop - subroutine to make graphics transformation and get some options

      implicit none

      include 'perplex_parameters.h'

      character yes*1

      integer jop0,iop1

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      integer  iop0 
      common / basic /iop0

      integer jcopt
c----------------------------------------------------------------------

      jop0 = 0

      
      if (jcopt.eq.3) then 
c                                     Should just prompt for axis
c                                     numeration.
         jop0 = iop0

      else if (iop0.eq.1) then

         write (*,1090) 
         read (*,1030) yes
         if (yes.eq.'y'.or.yes.eq.'Y') jop0 = 1

      end if 

      if (jop0.eq.1.and.jcopt.ne.3) then

         write (*,1070) 
         read (*,1030) yes
         iop1 = 0
 
         if (yes.eq.'y'.or.yes.eq.'Y') then
            write (*,1040) vnm(1),vmn(1),vmx(1)
            read (*,*) vmn(1),vmx(1)
            write (*,1040) vnm(2),vmn(2),vmx(2)
            read (*,*) vmn(2),vmx(2)
            iop1 = 1
            write (*,1080) 
         end if

      end if

      xmax = vmx(1)
      xmin = vmn(1)
      ymax = vmx(2)
      ymin = vmn(2)
      xlen = xmax - xmin
      ylen = ymax - ymin
c                                     set default char size, based
c                                     on font 7 (helv 14) the axes
c                                     length is about 85 characters
      dcx = xlen / 85d0 * cscale / xfac
      dcy = ylen / 85d0 * cscale 
c                                     set axes scales
      call psssc2 (xmin,xmax,ymin,ymax)

1030  format (a)
1040  format (/,'Enter new min and max for ',a8,' old values ',
     *        ' were: ',2(g11.5,1x))
1070  format (/,'Modify x-y limits (y/n)? ')
1080  format ('This may be sloppy. ')
1090  format (/,'Modify drafting options (y/n)?',/,
     *          '  answer yes to modify:',/,
     *          '   - field labeling',/,
     *          '   - x-y plotting limits',/,
     *          '   - axes numbering')
     
      end 

      subroutine psblrb (nline)
c----------------------------------------------------------------------
c psblrb - place an nline line text blurb (title, sectioning conditions) in the
c upper left corner of a diagram, text is in tname (csta8).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      integer i, nline, nchar

      double precision x,y

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      character*162 title
      common/ csta8 /title(4)
c---------------------------------------------------------------
      call pssctr (ifont,nscale,nscale,0d0)

      x = xmin - 1d1*dcx*nscale
      y = ymax + 2.3d1*dcy*nscale

      do i = 1, nline
         nchar = 162 
         call psublk (title(i),nchar)
         call pstext (x,y,title(i),nchar)
         y = y - 2.4*dcy*nscale
      end do 

      end 

      subroutine pslbtx 
c---------------------------------------------------------------------- 
c pslbtx - subprogram to read and draw text label data if present 
c at the end of a plot file, format:
c           x - y position 
c           text
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 text

      integer ier 

      double precision x,y

      call pssctr (ifont,nscale,nscale,0d0)

      do 

         read (n4,*,iostat=ier) x, y
         if (ier.ne.0) exit 
         read (n4,'(a)') text

         call pstext (x,y,text,10)

      end do 

      end 


      subroutine getfil (fname,lun,ier)

      implicit none

      include 'perplex_parameters.h'

      character*100 fname, yes*1

      integer ier, lun
c----------------------------------------------------------------------         
      open (lun,iostat=ier,file=fname,status='old')

      if (ier.ne.0) then
       
         write (*,1010) fname
         read (*,'(a)') yes

         if (yes.eq.'Y'.or.yes.eq.'y') return

         stop

      end if 

1010  format (/,'No such file as:',/,a,/,'Try again (y/n)?',/)

      end 

      subroutine redtab (lun)
c----------------------------------------------------------------------
c redtab - reads a Perple_X tab format file from logical unit number lun.
c chooses variables if necessary.
c see www.perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format.txt
c for details of the format.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ratio, warn1, eof

      integer i, j, dvar, dvar1, ier, lun, inc(l3)

      character tag*5, y

      double precision row(i11)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      integer ix,iy,mvar
      double precision z
      common/ dim   /z(nx,ny),ix,iy,mvar

      save warn1
      data warn1/.true./
c----------------------------------------------------------------------
      ratio = .false.

      read (lun,'(1x,a)') tag

      if (tag.ne.'6.6.6') then
         write (*,1000) tag
         stop
      end if 
c                                 text title
      read (lun,'(a)') title
c                                 number of independent variables
      read (lun,*) jvar

      if (jvar.gt.2) then 
         write (*,1010) jvar
         stop
      end if 
c                                 for each independent variable
      do i = 1, jvar
c                                 name
         read (lun,*) vnm(i)
c                                 minumum value 
         read (lun,*) vmn(i)
c                                 increment vlaue
         read (lun,*) dvr(i)
c                                 number of values
         read (lun,*) inc(i)

         vmx(i) = vmn(i) + (inc(i)-1)*dvr(i)

      end do 
c                                 number of dependent variables
      read (lun,*) mvar

      if (mvar.gt.i11) then
         write (*,1020) mvar,i11
         stop
      end if 
c                                 read names of dependent properties
      read (lun,'(80(a14,1x))') (dname(i),i=1,mvar)
c                                 read data 
      if (jvar.eq.2) then 

         ix = inc(1)
         iy = inc(2)

         if (ix.gt.nx) call error (1,dvr(1),ix,'NX, REDTAB')
         if (iy.gt.ny) call error (1,dvr(2),iy,'NY, REDTAB')

         if (mvar.gt.1) then
c                                 get dependent variable choice, single 
c                                 variable or ratio
            write (*,1030) 
            read (*,'(a)') y

            if (y.eq.'Y'.or.y.eq.'y') then

               ratio = .true.

               do 

                  write (*,1040) 'numerator'
                  write (*,1050) (i,dname(i),i=1,mvar)
                  read (*,*,iostat=ier) dvar 

                  if (ier.ne.0.or.dvar.lt.1.or.dvar.gt.mvar) then 
                     call rerr
                     cycle
                  end if 

                  exit 

               end do 
              
               do 

                  write (*,1040) 'denominator'
                  write (*,1050) (i,dname(i),i=1,mvar)
                  read (*,*,iostat=ier) dvar1 

                  if (ier.ne.0.or.dvar.lt.1.or.dvar.gt.mvar) then 
                     call rerr
                     cycle
                  end if 

                  exit 

               end do 

            else 

               do 

                  write (*,1060) 
                  write (*,1050) (i,dname(i),i=1,mvar)
                  read (*,*,iostat=ier) dvar
 
                  if (ier.ne.0.or.dvar.lt.1.or.dvar.gt.mvar) then 
                     call rerr
                     cycle
                  end if 

                  exit 

               end do 

            end if 

         else 

            dvar = 1

         end if      
c                                 make a title
         call mertxt (title,dname(dvar),title,1)
c                                 read the data 
         do j = 1, iy

            do i = 1, ix
c                                 read data with nan-check

               call redrow (row,lun,eof)

               if (ratio) then 

                  if (row(dvar1).ne.0d0) then

                     z(i,j) = row(dvar)/row(dvar1)

                  else 

                     if (warn1) then 
                        write (*,1100)
                        warn1 = .false.
                     end if 
 
                     if (isnan(nopt(7))) then 
                        z(i,j) = 0d0
                     else
                        z(i,j) = nopt(7)
                     end if 

                  end if 

               else 

                  z(i,j) = row(dvar)

               end if 

            end do 

         end do 

      else
c                                 a 1d table, determine number of
c                                 rows by reading to end of file
         iy = 1

         do
c                                 read data with eof and nan-check
            call redrow (row,lun,eof)

            if (eof) then 
               iy = iy - 1
               exit
            end if

            do i = 1, mvar
               z(iy,i) = row(i)
            end do 

            iy = iy + 1

            if (iy.gt.nx) call error (1,dvr(1),iy+100,'NX, REDTAB 2')

         end do
c                                 select variables
         if (mvar.eq.2) then 
c                                 no choices
            inv(1) = 1 
            inv(2) = 2

         else 
c                                 get independent variable (inv(1))
            write (*,1070) 
            write (*,1050) (i,dname(i),i=1,mvar)
            call rdnumb (row(1),0d0,inv(1),1,.false.)
            if (inv(1).lt.0.or.inv(1).gt.mvar) inv(1) = 1
c                                 get dependent variables (inv(2..dvar))
            write (*,1080)

            dvar = 1

            do 

               read (*,*,iostat=ier) dvar1 

               if (ier.ne.0.or.dvar1.gt.mvar.or.dvar1.lt.0) then 
                  call rerr
                  cycle
               end if    
               
               if (dvar1.eq.0) exit 

               dvar = dvar + 1
               inv(dvar) = dvar1

            end do 

            mvar = dvar

            if (dvar.eq.1) then 
               write (*,1090) 
               stop
            end if 

            read (dname(inv(1)),'(a8)') vnm(1) 

            if (mvar.eq.2) then 
               read (dname(inv(2)),'(a8)') vnm(2) 
               dname(inv(2)) = ' '
            else
               vnm(2) = ' '
            end if 

         end if 

      end if      

1000  format (/,'**error ver666** the version tag (',a,') in the input '
     *       ,'data file is inconsistent',/,'with this version of '
     *       ,'Perple_X, update from www.perplex.ethz.ch or modify the',
     *      /,'file format to be consistent with the description at:',/,
     *        'perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format',
     *        '.txt',/)
1010  format (/,'**error ver667** the data table is a function of ',i2,
     *       ' independent variables',/,'Perple_X plotting programs ',
     *       'cannot plot tables as a function of 2 variables.',/)
1020  format (/,'**error ver668** too many dependent variables ',i3,
     *       ' increase dimension i11 (',i3,')',/,
     *       'and recompile Perple_X',/)
1030  format (/,'Plot the ratio of two dependent variables (Y/N)?')
1040  format (/,'Select the ',a,' variable:',/)
1050  format (4x,i2,' - ',a)
1060  format (/,'Select the dependent variable to be contoured:')
1070  format (/,'Select x-axis variable [default variable 1]:',/)
1080  format (/,'Select y-axis variables from the list above',/,
     *          'one per line, enter zero to finish:',/)
1090  format (/,'You did not choose any dependent variables, I quit!',/)
1100  format (/,'**warning ver670** the denominator of a ratio is zero,'
     *      ,' infinite ratios will replaced',/,
     *       'by the bad_number value if bad_number is a number, ',
     *       'otherwise the ratio is set to 0',/)
      end 

      subroutine rdopt 
c----------------------------------------------------------------------
c rdopt - looks for the perplex_plot_option.dat file, if it finds
c the option file it scans the file for keywords and sets options
c accordingly.

c option variables - keyword associations
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier

      double precision dsx,dsy,dtx,dty,drot
      
      character*3 key*22, val, nval1*12, nval2*12,
     *            nval3*12,opname*40, strg*40,strg1*40
      
      character font*40
      common/ myfont /font
c----------------------------------------------------------------------
c                                 default option values:
c                                 -------------------------------------
c                                 picture transformation
      dsx = 0.18d0
      dsy = 0.18d0
      dtx = 160d0
      dty = 220d0
      drot = 0d0
c                                 bounding box
      bbox(1) = 0
      bbox(2) = 0
      bbox(3) = 800
      bbox(4) = 800
c                                 character scaling
      cscale = 1d0
c                                 axis number/label scale 
      nscale = 1.2d0
c                                 field label scale
      ascale = 0.72d0
c                                 b-splines to fit data curves
      spline = .true.
c                                 half interval axis tick marks
      half = .true.
c                                 tenth interval axis tick marks
      tenth = .false.
c                                 superpose grid
      lgrid = .false. 
c                                 field_fill
      fill = .true.
c                                 field_fill_scale => adjust gray scale increments
      plopt(1) = .true.
c                                 field_fill_zero
      plopt(2) = .true.
c                                 numeric_field_label
      plopt(3) = .false.
c                                 george's bracket/data plots
      plopt(4) = .false.
c                                 field labels
      label = .true. 
c                                 font - default font, the ifont choice
c                                 is for idraw postscript and is probably
c                                 archaic
      ifont = 7
      font = 'Helvetica'
c                                 line_width - default line width (points, integer)
      width = 1d0
c                                 plot_aspect_ratio - (x_axis_length/y-axis_length)
      xfac = 1d0 
c                                 replicate_label - minimum separation before writing
c                                 a replicate label
      rlabel = 0.025d0
c                                 -------------------------------------
c                                 look for file
      opname = 'perplex_plot_option.dat'
      
      open (n8, file = opname, iostat = ier, status = 'old')
c                                 no option file
      if (ier.ne.0) write (*,1120) opname
c                                 read cards to end of option file
      do while (ier.eq.0) 

         call redcd1 (n8,ier,key,val,nval1,nval2,nval3,strg,strg1)
c                                 ier ne 0 bad record or eof
         if (ier.ne.0) exit 
c                                 if here we have a keyword and value
         if (key.eq.'font') then
         
            font = strg 

         else if (key.eq.'axis_label_scale') then 
c                                 eps bounding box.   
            read (strg,*) nscale

         else if (key.eq.'field_label_scale') then 
c                                 eps bounding box.   
            read (strg,*) ascale

         else if (key.eq.'text_scale') then 
c                                 eps bounding box.   
            read (strg,*) cscale

         else if (key.eq.'picture_transformation') then
c                                 picture transformation.   
            read (strg1,*) dsx,dsy,dtx,dty,drot  

         else if (key.eq.'half_ticks') then 
c                                 half interval tick marks on axes  
            read (strg,*) half
          
         else if (key.eq.'tenth_ticks') then 
c                                 tenth interval tick marks on axes   
            read (strg,*) tenth
                
         else if (key.eq.'grid') then 
c                                 gridding   
            read (strg,*) lgrid

         else if (key.eq.'field_fill') then 
c                                 fill phase fields (gridded minimization).   
            read (strg,*) fill

         else if (key.eq.'field_fill_scale') then 
c                                 scale fills according to variance range.   
            read (strg,*) plopt(1)

         else if (key.eq.'field_fill_zero') then 
c                                 scale from minimum variance rather than 
c                                 divariant.   
            read (strg,*) plopt(2)

         else if (key.eq.'field_label') then 
c                                 label phase fields (gridded minimization).   
            read (strg,*) label

         else if (key.eq.'numeric_field_label') then 
         
            read (strg,*) plopt(3)

         else if (key.eq.'plot_extra_data') then 
         
            read (strg,*) plopt(4)

         else if (key.eq.'splines') then 
c                                 eps bounding box.   
            read (strg,*) spline
                    
         else if (key.eq.'bounding_box') then 
c                                 eps bounding box.   
            read (strg1,*) bbox
            
         else if (key.eq.'line_width') then 
c                                 default line width in points (integer).   
            read (strg,*) width
            
         else if (key.eq.'plot_aspect_ratio') then 
c                                 x_axis_lenght/y-axis_length.   
            read (strg,*) xfac
            
         else if (key.eq.'replicate_label') then 
c                                 replicate lable tolerance for PSSECT
c                                 replicate fields are labled only if they
c                                 are further apart than the normalized distance 
c                                 specified by this keyword.   
            read (strg,*) rlabel
      
         else if (key.eq.'page_size') then 
c                                 perplex_pdf, do nothing
         else if (key.eq.'new_font') then 
c                                 perplex_pdf, do nothing
         else if (key.eq.'plot_output_type') then 
c                                 perplex_pdf, do nothing
         else if (key.ne.'|') then 

            write (*,1110) key

         end if 

      end do 
                
      close (n8)
c                                 --------------------------------------    
      dsx = dsx * xfac
 
      call psstrn (dsx,dsy,dtx,dty,drot)
c                                 --------------------------------------
c                                 output 
      write (*,1000) 

      write (*,1010) nscale, bbox, fill, label, plopt(3), rlabel, 
     *               ascale, font, lgrid, half, width, dsx, 
     *               dsy, dtx, dty, drot, xfac, spline, tenth, 
     *               cscale

      write (*,1020) 
c                                 -------------------------------------
1000  format (/,'Perple_X plot options are currently set as:',//,
     *      '    Keyword:               Value:     Permitted values ',
     *          '[default]:')
1010  format (4x,'axis_label_scale       ',f4.2,7x,'[1.2] (rel)',/,
     *        4x,'bounding_box :',/,
     *        28x,i4,6x,'[0] x-min (pts)',/,
     *        28x,i4,6x,'[0] y-min (pts)',/,
     *        28x,i4,6x,'[800] x-length (pts)',/,
     *        28x,i4,6x,'[800] y-length (pts)',/,
     *        4x,'field_fill             ',l1,10x,'[T] F',/,
     *        4x,'field_label            ',l1,10x,'[T] F',/,
     *        4x,'numeric_field_label    ',l1,10x,'[F] T, if T ',
     *           'PSSECT writes list to *_assemblages.txt',/,
     *        4x,'replicate_label        ',f5.3,6x,'0->1 [0.025]',/,
     *        4x,'field_label_scale      ',f4.2,7x,'[0.72] (rel)',/,
     *        4x,'font                   ',a,/,
     *        4x,'grid                   ',l1,10x,'[F] T',/,
     *        4x,'half_ticks             ',l1,10x,'[T] F',/,
     *        4x,'line_width             ',f4.2,7x,'0-99 [1.] (pts)',/,
     *        4x,'picture_transformation :',/,
     *        28x,g9.3,1x,'[0.18] x-scale (rel)',/,
     *        28x,g9.3,1x,'[0.18] y-scale (rel)',/,
     *        28x,g9.3,1x,'[0.18] x-translation (pts)',/,
     *        28x,g9.3,1x,'[0.18] y-translation (pts)',/,
     *        28x,g9.3,1x,'[0.0]  rotation (deg)',/,
     *        4x,'plot_aspect_ratio      ',f5.3,6x,
     *           '[1.0] x_axis_length/y_axis_length',/,
     *        4x,'splines                ',l1,10x,'[T] F',/,
     *        4x,'tenth_ticks            ',l1,10x,'[F] T',/,
     *        4x,'text_scale             ',f5.3,6x,'[1.] (rel)',/,
     *        4x,'plot_extra_data        ',l1,10x,
     *           '[T] F, to plot, e.g., experimental observations',/)

1020  format (/,'To change these options edit or create ',
     *        'the plot option file',/,'See: ',
     *        'www.perplex.ethz.ch/perplex_plot_options.html',/)
1110  format (/,'Warning: unrecognized option text: ',a,/,
     *       'If the text is intentional, check spelling and case.',/) 
1120  format (/,'Warning: the plot option file: ',a,/,
     *       'was not found, default option values will be used.',/) 

      end

      subroutine redrow (row,lun,eof)
c----------------------------------------------------------------------
c redrow - reads a row of table entries and checks for bad data.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical eof, warn1

      integer k, ier, lun

      character numbs(i11)*14

      double precision row(*)

      integer ix,iy,mvar
      double precision z
      common/ dim   /z(nx,ny),ix,iy,mvar

      save warn1
      data warn1/.true./
c----------------------------------------------------------------------
c                                 format/char read necessary for NaN check
c                                 with compaq fortran:
      read (lun,'(80(a14,1x))',iostat=ier) (numbs(k),k=1,mvar)

      if (ier.ne.0) then 
c                                 presumably eof
         eof = .true.
         return 

      else 

         eof = .false.

      end if 
c                                 convert strings to numbers:
      do k = 1, mvar

         read (numbs(k),'(g14.7)',iostat=ier) row(k)

         if (ier.ne.0.or.isnan(row(k))) then 
    
             if (warn1) then 
                call warn (4,row(1),k,numbs(k))
                warn1 = .false.
             end if 

             row(k) = 0d0

          end if 

      end do

      end 

      subroutine psdat
c----------------------------------------------------------------
c psdat - subroutine to add data points and brackets to a plot

c contributed by George Helffrich, ELSI, 9/23/2021.

c lines starting with * are comments, and text on any line after # or | is
c ignored.

c Input line formats:
c 1) this draws a line on the plot
c > L      ## header indicating group of lines; ends with > (or new > L)
c   x1 y1
c   x2 y2
c   .. ..
c   xn yn
c >
c 2) points
c   xp yp symb size gfill        ## 5 blank separated fields
c 3) error bars
c   xp yp ex ey symb size gfill  ## 7 blank separated fields

c For all types: xi, yi - X and Y coordinates on plot.

c For types 2) & 3): symb = symbol code (see comments below in code for shapes),
c   size = real number (1 is normal), gfill = gfill code (0 = none, 1-7 gray
c   scale (lightest->darkest), 8-15 patterns.

c For type 3): ex, ey - horiz. & vert. half-width of error bar.

c IN FUTURE:
c extend > L to add optional line feature keywords, e.g.
c    >L style=n width=m symbol=j
c to set a line style (dashed, dotted, etc.), line width and whether to plot
c a symbol at each point too.  Default would be solid line, std. width, no
c symbols.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision SMALL, RADIUS, SQRC, DMDC, TRC0, TRC1, TRC2,
     *   SQRT2
      parameter (SMALL=2*0.25d0, RADIUS=2*0.375d0,
     * SQRC=0.88622692545275801364d0, ! /* sqrt(pi / 4) */
     * DMDC=1.25331413731550025119d0, !	/* sqrt(pi / 4) * sqrt(2) */
     * TRC0=1.55512030155621416073d0, !	/* sqrt(4 * pi/(3 * sqrt(3))) */
     * TRC1=1.34677368708859836060d0, !	/* TRC0 * sqrt(3) / 2 */
     * TRC2=0.77756015077810708036d0, !	/* TRC0 / 2 */
     * SQRT2=1.414213562373095048801688724210d0 ! sqrt(2)
     * )

      character file*72

      integer type, ier, symb, gfill, ix, npts, ifg, ibg

      double precision x, y, size, sigx, sigy, xy(2), sig(2),
     *   lx(nx), ly(nx), r, xc, yc, rline, cwidth, xx(4), yy(4)

      character line*128

      integer nblen
      external nblen

      equivalence (x,xy(1)),(y,xy(2))
      equivalence (sigx,sig(1)),(sigy,sig(2))

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      write (*,1000) 
      read (*,'(a)') file
1000  format ('Enter file name for the plot_extra_data option:')

      open (n8,file=file,status='old',iostat=ier)
      if (ier.ne.0) then
         write (*,*) '**Bad plot annotation file: ',
     *         file(1:nblen(file))
         return
      endif

      do
100      continue
         read(n8,'(a)',iostat=ier) line
         if (ier.ne.0) exit
c                                 skip blank or comment lines
         if (0.ne.index('*#|',line(1:1))) cycle
         if (line.eq.' ') cycle
c                                 omit comment text
         ix = index(line, '#')
         if (ix.eq.0) ix = index(line, '|')
         if (ix.eq.0) ix = nblen(line)

c                                 check for line data
c                                 in future, parse '> L' line for parameters
c                                 like lty=, lwd=
         if (ix.ge.3 .and. line(1:3) .eq. '> L') then

            do npts=1,nx

               do 
                  read(n8,'(a)',iostat=ier) line
                  if (ier.ne.0) exit
                  if (line(1:1) .eq. '>') exit
                  if (0.ne.index('*#|',line(1:1))) cycle
                  ix = index(line, '#')
                  if (ix.eq.0) ix = index(line, '|')
                  if (ix.eq.0) ix = nblen(line)
                  read(line(1:ix),*,iostat=ier) lx(npts),ly(npts)
                  if (ier.eq.0) exit
                   write(*,*) '**Bad line point: ',line(1:nblen(line))
               end do

               if (ier.ne.0) exit

               if (line(1:1).eq.'>') goto 100

            end do

            if (npts.gt.2) then
               do ix=1,npts-2
                  call psline(lx(ix),ly(ix),lx(ix+1),ly(ix+1),1d0,1d0)
               end do
            end if

            if (ier.ne.0) exit 

            if (line(1:1).eq.'>') go to 100

            cycle

         end if

         type = 1
         read(line(1:ix),*,iostat=ier)
     *      xy(iv1), xy(iv2), sig(iv1), sig(iv2), symb, size, gfill
         if (ier.ne.0) then
            type = 0
            read(line(1:ix),*,iostat=ier)
     *         xy(iv1), xy(iv2), symb, size, gfill
         end if
         if (ier.ne.0) then
            write(*,*) '**Bad point file line: ',line(1:nblen(line))
            cycle
         end if
         if (gfill.lt.0 .or. gfill.gt.15) then
            write(*,*) '**Bad fill in line: ',line(1:nblen(line))
            cycle
         end if

c                                 various plot symbols
         rline = 1d0
         cwidth = 1d0
         if (symb.eq.0) then
c           0: square
            xc = dcx*size*RADIUS
            yc = dcy*size*RADIUS
            call psrect (x-xc, x+xc, y-yc, y+yc, rline, cwidth, gfill)
         else if (symb.eq.1) then
c           1: circle
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call pselip (x,y, xc, yc, 1d0, 0d0, gfill, ifg, ibg)
         else if (symb.eq.2) then
c           2: triangle - point up
            r =  TRC0 * RADIUS * dcy * size
            yc = TRC2 * RADIUS * dcy * size
            xc = TRC1 * RADIUS * dcx * size
            xx(1) = x
            yy(1) = y+r
            xx(2) = x+xc
            yy(2) = y-yc
            xx(3) = x-xc
            yy(3) = y-yc
            call pspygn (xx,yy,3,rline,cwidth,gfill)
         else if (symb.eq.3) then
c           3: plus
            xc = SQRT2*RADIUS*dcx*size
            yc = SQRT2*RADIUS*dcy*size
            call psline (x-xc,y,x+xc,y,rline,cwidth)
            call psline (x,y-yc,x,y+yc,rline,cwidth)
         else if (symb.eq.4) then
c           4: times
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call psline (x-xc,y-yc,x+xc,y+yc,rline,cwidth)
            call psline (x-xc,y+yc,x+xc,y-yc,rline,cwidth)
         else if (symb.eq.5) then
c           5: diamond
            xc = SQRT2 * RADIUS * dcx * size
            yc = SQRT2 * RADIUS * dcy * size
            xx(1) = x-xc
            yy(1) = y
            xx(2) = x
            yy(2) = y+yc
            xx(3) = x+xc
            yy(3) = y
            xx(4) = x
            yy(4) = y-yc
            call pspygn (xx, yy, 4, rline, cwidth, gfill)
         else if (symb.eq.6) then
c           6: triangle - point down
            r =  TRC0 * RADIUS * dcy * size
            yc = TRC2 * RADIUS * dcy * size
            xc = TRC1 * RADIUS * dcx * size
            xx(1) = x
            yy(1) = y-r
            xx(2) = x+xc
            yy(2) = y+yc
            xx(3) = x-xc
            yy(3) = y+yc
            call pspygn (xx, yy, 3, rline, cwidth, gfill)
         else if (symb.eq.7) then
c           7: square and times superimposed
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call psrect (x-xc, x+xc, y-yc, y+yc, rline, cwidth, gfill)
            call psline (x-xc,y-yc,x+xc,y+yc,rline,cwidth)
            call psline (x-xc,y+yc,x+xc,y-yc,rline,cwidth)
         else if (symb.eq.8) then
c           8: plus and times superimposed
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call psline (x-xc,y-yc,x+xc,y+yc,rline,cwidth)
            call psline (x-xc,y+yc,x+xc,y-yc,rline,cwidth)
            xc = SQRT2*RADIUS*dcx*size
            yc = SQRT2*RADIUS*dcy*size
            call psline (x-xc,y,x+xc,y,rline,cwidth)
            call psline (x,y-yc,x,y+yc,rline,cwidth)
         else if (symb.eq.9) then
c           9: diamond and plus superimposed
            xc = SQRT2 * RADIUS * dcx * size
            yc = SQRT2 * RADIUS * dcy * size
            call psline (x-xc,y,x+xc,y,rline,cwidth)
            call psline (x,y-yc,x,y+yc,rline,cwidth)
            xx(1) = x-xc
            yy(1) = y
            xx(2) = x
            yy(2) = y+yc
            xx(3) = x+xc
            yy(3) = y
            xx(4) = x
            yy(4) = y-yc
            call pspygn (xx, yy, 4, rline, cwidth, gfill)
         else if (symb.eq.10) then
c           10: circle and plus superimposed
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call pselip (x,y, xc, yc, 1d0, 0d0, gfill, ifg, ibg)
            call psline (x-xc,y,x+xc,y,rline,cwidth)
            call psline (x,y-yc,x,y+yc,rline,cwidth)
         else if (symb.eq.11) then
c           11: superimposed triangles
            r =  TRC0 * RADIUS * dcy * size
            yc = TRC2 * RADIUS * dcy * size
            yc = 0.5 * (yc + r)
            xc = TRC1 * RADIUS * dcx * size
            xx(1) = x
            yy(1) = y-r
            xx(2) = x+xc
            yy(2) = y+yc
            xx(3) = x-xc
            yy(3) = y+yc
            call pspygn (xx, yy, 3, rline, cwidth, gfill)
            xx(1) = x
            yy(1) = y+r
            xx(2) = x+xc
            yy(2) = y-yc
            xx(3) = x-xc
            yy(3) = y-yc
            call pspygn (xx, yy, 3, rline, cwidth, gfill)
         else if (symb.eq.12) then
c           12: square and plus superimposed
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call psrect (x-xc, x+xc, y-yc, y+yc, rline, cwidth, gfill)
            call psline (x-xc,y,x+xc,y,rline,cwidth)
            call psline (x,y-yc,x,y+yc,rline,cwidth)
         else if (symb.eq.13) then
c           13: circle and times superimposed
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call pselip (x,y, xc, yc, 1d0, 0d0, gfill, ifg, ibg)
            call psline (x-xc,y-yc,x+xc,y+yc,rline,cwidth)
            call psline (x-xc,y+yc,x-xc,y-yc,rline,cwidth)
         else if (symb.eq.14) then
c           14: square and point-up triangle superimposed
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            xx(1) = x
            yy(1) = y+yc
            xx(2) = x+xc
            yy(2) = y-yc
            xx(3) = x-xc
            yy(3) = y-yc
            call psrect (x-xc, x+xc, y-yc, y+yc, rline, cwidth, gfill)
            call pspygn (xx, yy, 3, rline, cwidth, 0)
         else if (symb.eq.15) then
c           15: filled square
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call psrect (x-xc, x+xc, y-yc, y+yc, rline, cwidth, gfill)
         else if (symb.eq.16) then
c           16: filled circle
            xc = RADIUS * dcx * size;
            yc = RADIUS * dcy * size;
            call pselip (x,y, xc, yc, rline, cwidth, gfill, ifg, ibg)
         else if (symb.eq.17) then
c           17: filled point-up triangle
            r =  TRC0 * RADIUS * dcy * size
            yc = TRC2 * RADIUS * dcy * size
            xc = TRC1 * RADIUS * dcx * size
            xx(1) = x
            yy(1) = y+r
            xx(2) = x+xc
            yy(2) = y-yc
            xx(3) = x-xc
            yy(3) = y-yc
            call pspygn (xx, yy, 3, rline, cwidth, gfill)
         else if (symb.eq.18) then
c           18: filled diamond
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            xx(1) = x-xc
            yy(1) = y
            xx(2) = x
            yy(2) = y+yc
            xx(3) = x+xc
            yy(3) = y
            xx(4) = x
            yy(4) = y-yc
            call pspygn (xx, yy, 4, rline, cwidth, gfill)
         else if (symb.eq.19) then
c           19: filled circle
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call pselip (x,y, xc, yc, rline, cwidth, gfill, ifg, ibg)
         else if (symb.eq.20) then
c           20: small circle
            xc = SMALL * dcx * size
            yc = SMALL * dcy * size
            call pselip (x,y, xc, yc, rline, cwidth, gfill, ifg, ibg)
         else if (symb.eq.21) then
c           21: circles
            xc = RADIUS * dcx * size
            yc = RADIUS * dcy * size
            call pselip (x,y, xc, yc, rline, cwidth, gfill, ifg, ibg)
         else if (symb.eq.22) then
c           22: squares
            xc = RADIUS * SQRC * dcx * size
            yc = RADIUS * SQRC * dcy * size
            call psrect (x-xc, x+xc, y-yc, y+yc, rline, cwidth, gfill)
         else if (symb.eq.23) then
c           23: diamonds
            xc = RADIUS * DMDC * dcx * size
            yc = RADIUS * DMDC * dcy * size
            xx(1) = x          
            yy(1) = y-yc 
            xx(2) = x+xc
            yy(2) = y 
            xx(3) = x          
            yy(3) = y+yc 
            xx(4) = x-xc
            yy(4) = y 
            call pspygn (xx, yy, 4, rline, cwidth, gfill)
         else if (symb.eq.24) then
c           24: triangle (point up)
            r =  TRC0 * RADIUS * dcy * size
            yc = TRC2 * RADIUS * dcy * size
            xc = TRC1 * RADIUS * dcx * size
            xx(1) = x
            yy(1) = y+r
            xx(2) = x+xc
            yy(2) = y-yc
            xx(3) = x-xc
            yy(3) = y-yc
            call pspygn (xx, yy, 3, rline, cwidth, gfill)
         else if (symb.eq.25) then
c           25: triangle (point down)
            r =  TRC0 * RADIUS * dcy * size
            yc = TRC2 * RADIUS * dcy * size
            xc = TRC1 * RADIUS * dcx * size
            xx(1) = x
            yy(1) = y-r
            xx(2) = x+xc
            yy(2) = y+yc
            xx(3) = x-xc
            yy(3) = y+yc
            call pspygn (xx, yy, 3, rline, cwidth, gfill)
         else
            write(*,*) '**Bad symbol in line: ',line(1:nblen(line))
         endif

c                                 add error bars if requested
         if (type.eq.1) then
            call psmove (x,y)
            call psrlin (0d0, sigy, 1d0,1d0)
            call psmove (x,y)
            call psrlin (0d0,-sigy, 1d0,1d0)
            call psmove (x,y)
            call psrlin (sigx, 0d0, 1d0,1d0)
            call psmove (x,y)
            call psrlin (-sigx,0d0, 1d0,1d0)
         end if

      end do

      close(n8)

      end

      integer function nblen(str)
c--------------------------------------------------------------- 
c george's function for psgrid, 
c replicates a function already in perple_x
c nblen - function to return nonblank length of a string
      character str*(*)
      integer i

      do i=len(str),1,-1
         if (str(i:i) .ne. ' ') exit
      end do
      nblen = i
      end