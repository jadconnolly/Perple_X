      PROGRAM PSTABL 
c----------------------------------------------------------------------
c A program to plot Perple_X "tab" (table) format files.
c                                              JADC March 15, 2011.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier

      logical ratio
 
      character y*1

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer ix,iy,mvar
      double precision z
      common/ dim   /z(nx,ny),ix,iy,mvar

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer  iop0 
      common / basic /iop0
c----------------------------------------------------------------------
c                                 version info
      call vrsion (6)
c                                 read plot option file, set
c                                 default transformation
      call rdopt 

      ratio = .false.
c                         
      do 
c                                 get input file 
         write (*,'(/,a)') 
     *         'Enter the tab file name [without the .tab suffix]:'

         call readrt
         call mertxt (tfname,prject,'.tab',0)

         open (n4,iostat=ier,file=tfname,status='old')

         if (ier.eq.0) exit
       
         write (*,1010) tfname
         read (*,'(a)') y
         if (y.eq.'Y'.or.y.eq.'y') cycle 

         stop

      end do 
c                                 read the table file, get plot variable
c                                 choices
      call redtab (n4)                    
c                                 jvar is the dimension of the table (1- or 2-d)
      if (jvar.eq.2) then
c                                 for 2d case query for ratio plots
         write (*,1020) 
         read (*,'(a)') y

         if (y.eq.'Y'.or.y.eq.'y') then

            ratio = .true.
c                                 get second table file name 
            do 

               write (*,1040)               
               call readrt
               call mertxt (tfname,prject,'.tab',0)              
               open (n5,iostat=ier,file=tfname,status='old')

               if (ier.eq.0) exit
       
               write (*,1010) tfname
               read (*,'(a)') y

               if (y.eq.'Y'.or.y.eq.'y') cycle 

               stop

            end do  

         end if

      end if 
c                                 open output file 
      call psopen
c                                 allow drafting options prompt
      iop0 = 0

      write (*,'(/,a)') 'Modify the default plot (y/n)?'
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') iop0 = 1

      if (jvar.eq.2) then 
c                                 contour plotting
         call pstab2 (ratio)

      else 
c                                 x-y plotting
         call pstab1

      end if

      if (plopt(4)) call psdat
 
      call psclos
 
      close (n4)

c1000  format (/,'Enter the complete plot file name [e.g., ',
c     *       'my_project.tab or my_project.ctr]:')


1010  format (/,'**warning ver191** cannot find file',/,a,/,
     *       'run WERAMI/FRENDLY to generate the ',
     *       'file or try a different name (y/n)?')
1020  format (/,'Contour the ratio of values in separate tab ',
     *       'files (y/n)?',//,'If you answer yes the data from the ',
     *       'file just read will define the',/,'numerator of the '
     *       'ratio and you will be prompted next for a file',/,
     *       'containing the data for the denominator.')
1040  format (/,'Enter the name of the tab file that ',
     *          'contains the denominator data:')

      end

      subroutine pstab2 (ratio)
c--------------------------------------------------------------------- 
c pstab2 - subroutine to plot (contour) 2d-table data
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ratio

      character y*1

      integer i,j,jmn,imn,imx,jop0,ncon,jmx,iop1,jy,jx

      double precision xpmn,xpmx,cmin,cmax,dcon,ypmx,ypmn,z0min,z0max

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer  iop0 
      common / basic /iop0

      integer ix,iy,mvar
      double precision z,zt 
      common/ dim   /z(nx,ny),ix,iy,mvar
      common/ dim1  /zt(nx,ny)

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      double precision zmin,zmax
      common/ stuff /zmax,zmin

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
c----------------------------------------------------------------------
c                                 if ratio, making a plot of data 
c                                 in two files.
      if (ratio) then 
c                                 copy the z-data in the first file to
c                                 zt. 
         jx = ix
         jy = jy 
         do i = 1, ix
            do j = 1, iy
               zt(i,j) = z(i,j)
            end do
         end do 

         call redtab (n5)

         if (jx.ne.ix.or.jy.ne.iy) then 
            write (*,'(a)') 'the plots do not have consistent ',
     *                      'dimensions'
            stop
         end if 

         do i = 1, ix
            do j = 1, iy
               z(i,j) = zt(i,j)/z(i,j)
            end do
         end do

      end if 

      if (iop0.eq.1) then 
c                                 modify axes option, data limits, etc
         write (*,'(a)') 
     *         'Contour the log10 of the dependent variable (y/n)?'
         read (*,'(a)') y

         if (y.eq.'y'.or.y.eq.'Y') then 
            do j = 1, iy
               do i = 1, ix
                  if (z(i,j).ne.0d0) z(i,j) = dlog10(dabs(z(i,j)))
               end do 
            end do  
         end if 
c                                 plot limits
         write (*,'(/,a)') 'Reset plot limits (y/n)?' 
         read (*,'(a)') y

         if (y.eq.'y'.or.y.eq.'Y') then 

            write (*,1070) vmx(2),vmn(2),vmx(1),vmn(1)
            read (*,*) ypmx,ypmn,xpmx,xpmn
           
            imn = int(xpmn/dvr(1)) + 1
            imx = int(xpmx/dvr(1)) + 1
            jmn = int(ypmn/dvr(2)) + 1
            jmx = int(ypmx/dvr(2)) + 1

            ix = (imx-imn+1)
            iy = (jmx-jmn+1)
            vmx(1) = xpmn + (ix-1)*dvr(1)
            vmx(2) = ypmn + (iy-1)*dvr(2)
            vmn(1) = xpmn
            vmn(2) = ypmn 
c                                      reload mini matrix:
            do i = 1, ix
               do j = 1, iy
                  z(i,j) = z(i+imn-1,j+jmn-1)
               end do
            end do

         end if
 
      end if 
c                                 get some options and
c                                 set up transformations
      call psaxop (3,jop0,iop1)
        
      zmin = 1d9
      zmax = -1d9
      z0min = 1d30
      z0max = -1d30
c                                      set up contour intervals                                      
      do i = 1, ix
         do j = 1, iy 
            if (z(i,j).lt.zmin) zmin = z(i,j)
            if (z(i,j).gt.zmax) zmax = z(i,j)
            if (z(i,j).lt.z0min.and.z(i,j).ne.0d0) z0min = z(i,j)
            if (z(i,j).gt.z0max.and.z(i,j).ne.0d0) z0max = z(i,j)
         end do 
      end do 
c                                      set up contour intervals
      write (*,1020) zmin, zmax, z0min, z0max 
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') then

         write (*,'(a)') 'Enter min, max and interval for contours:'
         read (*,*) cmin, cmax, dcon
         ncon = int((cmax-cmin)*1.0001d0/dcon) + 1

      else 

         dcon = (zmax-zmin)/11.
         cmin = zmin + 0.5d0 * dcon
         ncon = 11

      end if 

      call pscontor (cmin,ncon,dcon)
 
      call psaxes (jop0)

1020  format ('Contoured variable range: ',g14.6,'->',g14.6,/,
     *        'Range excluding zero values: ',g14.6,'->',g14.6,//,
     *        'Modify default contour interval (y/n)?')
1070  format (/,'Old values were: ',4(g12.4),/,'Enter new values:')

      end

      subroutine pstab1 
c--------------------------------------------------------------------- 
c pstab1 - subroutine to plot 1d-tables.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character yes*1,short*10

      integer i,j,k,jop0,npts,iop1,i1,i0

      double precision x(l5),y(l5),rline,w,x0min,x0max,y0min,y0max

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer  iop0 
      common / basic /iop0

      integer ix,iy,mvar
      double precision z
      common/ dim   /z(nx,ny),ix,iy,mvar

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title
c----------------------------------------------------------------------
      if (iop0.eq.1) then 
c                                 log transformations
         write (*,1040) 
         read (*,'(a)') yes

         if (yes.eq.'y'.or.yes.eq.'Y') then

            do i = 1, mvar 

               write (*,1050) dname(inv(i))
               read (*,'(a)') yes

               if (yes.eq.'y'.or.yes.eq.'Y') then

                  read (dname(inv(i)),'(a10)') short
                  write (dname(inv(i)),'(a4,a10)') 'log_',short

                  do j = 1, iy
                     
                     z(j,inv(i)) = dlog10(z(j,inv(i)))
                    
                  end do 
           
               end if 

            end do 
     
         end if 

      end if 

      xmin = 1d99
      xmax = -1d99
      ymin = 1d99
      ymax = -1d99
      y0min = 1d99
      y0max = -1d99
      x0min = 1d99
      x0max = -1d99
c                                 find x-data limits        
      i = inv(1)

      do j = 1, iy 
            
         w = z(j,i) 

         if (w.lt.xmin) xmin = w
         if (w.gt.xmax) xmax = w
         if (w.lt.x0min.and.w.ne.0d0) x0min = w         
         if (w.gt.x0max.and.w.ne.0d0) x0max = w   

      end do 
c                                 find y-data limits
      do k = 2, mvar
         
         i = inv(k)

         do j = 1, iy 
            
            w = z(j,i) 

            if (w.lt.ymin) ymin = w
            if (w.gt.ymax) ymax = w
            if (w.lt.y0min.and.w.ne.0d0) y0min = w         
            if (w.gt.y0max.and.w.ne.0d0) y0max = w   

         end do 

      end do 

      write (*,1020) dname(inv(1)),xmin,xmax,x0min,x0max
      write (*,1020) 'y-axis variables',ymin,ymax,y0min,y0max

      if (iop0.eq.1) then 
c                                 reset plot limits
         write (*,1060) 
         read (*,'(a)') yes

         if (yes.eq.'y'.or.yes.eq.'Y') then

             write (*,1070) 'lower limit',dname(inv(1)),xmin    
             call rdnumb (xmin,xmin,inv(1),inv(1),.true.)
     
             write (*,1070) 'upper limit',dname(inv(1)),xmax    
             call rdnumb (xmax,xmax,inv(1),inv(1),.true.) 
  
             write (*,1070) 'lower limit','y-axis',ymin    
             call rdnumb (ymin,ymin,inv(1),inv(1),.true.)
     
             write (*,1070) 'upper limit','y-axis',ymax    
             call rdnumb (ymax,ymax,inv(1),inv(1),.true.)  

          end if 

      end if   

      vmx(1) = xmax
      vmx(2) = ymax
      vmn(1) = xmin
      vmn(2) = ymin

      rline = 1d0
c                                 get some options and
c                                 set up transformations
      call psaxop (3,jop0,iop1)
c                                 set character transformation? 
      call pssctr (ifont,nscale,nscale,0d0)
c                                 plot loop: 
      i0 = inv(1)

      do i = 2, mvar

         npts = 0 
         i1 = inv(i)
c                                 filter data 
         do j = 1, iy

            if (z(j,i0).lt.xmin.or.z(j,i0).gt.xmax.or.
     *          z(j,i1).lt.ymin.or.z(j,i1).gt.ymax) cycle

            npts = npts + 1

            if (npts.gt.l5) call error (1,x(1),iy,'L5')
            
            x(npts) = z(j,i0)
            y(npts) = z(j,i1)

         end do 

         if (npts.lt.2) cycle
c                                 draw data 
         if (spline) then
            call psbspl (x,y,npts,rline,width,0)
         else 
            call pspyln (x,y,npts,rline,width,0)
         end if
c                                 label curve
         k = npts/2
         call pstext (x(k)+dcx,y(k)+4d0*dcy,dname(inv(i)),14)

      end do  
c                                 write a title
      call pstext (xmin,ymax+1d1*dcy,title,162)

      call psaxes (jop0)
 
1020  format ('Range of ',a,' is:',g14.6,'->',g14.6,/,
     *        'Range excluding zero values is:',g14.6,'->',g14.6,/)
1040  format (/,'Convert to logarithmic variables (y/n)? ',
     *          'NOTE: if a variable is',/,'already logarithmic, then ',
     *          'answering yes will plot its log-log value.')
1050  format (/,'Plot base 10 logarithm of ',a,'  (y/n)?')
1060  format (/,'Reset plot limits (y/n)?')
1070  format (/,'Enter ',a,' for ',a,' [default=',g14.7,']:')

      end
