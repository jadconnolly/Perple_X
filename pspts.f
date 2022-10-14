 
c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy and
c Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c PSVDRAW - a program to generate postscript code for x-y plots and chemographic
c phase diagrams.  The input data format is consistent with the files generated
c by the programs VERTEX and FRENDLY.
 
c Please do not distribute any part of this source.
 
      PROGRAM PSPTS

      implicit none

      include 'perplex_parameters.h'

      integer ier

      character*1 yes

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer  iop0 
      common / basic /iop0
c----------------------------------------------------------------------
c                                 version info
      call vrsion (6)
c                                 set iop0 to 1 to allow 
c                                 drafting ptompts
      iop0 = 1

      do 
c                                 get input file 
         write (*,1000) 
          
         call readrt

         call mertxt (tfname,prject,'.pts',0)
         
         open (n4,iostat=ier,file=tfname,status='old')

         if (ier.ne.0) then
       
            write (*,1010) tfname
            read (*,'(a)') yes

            if (yes.eq.'Y'.or.yes.eq.'y') cycle 

            stop

         end if

         exit  

      end do 
c                                 read plot option file, set
c                                 default transformation
      call rdopt 
c                                 open output file 
      call psopen 

      call psxypl 
 
      call psclos
 
      close (n4)
 
1000  format (/,'Enter the POINT plot file name [',
     *       'without the .pts suffix]:')
1010  format (/,'**warning ver191** cannot find file:',/,a,/,
     *       'run WERAMI to generate the ',
     *       'file or try a different name (y/n)?')
 
      end

      subroutine psxypl
c--------------------------------------------------------------------- 
c psxypl - subroutine to output x-y plot.
 
      implicit none
 
      include 'perplex_parameters.h'

      integer iop1, jop0, ier, isym

      double precision x, y
 
      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)   
c--------------------------------------------------------------------- 
      jvar = 2
      vmn(1) = 1d30
      vmx(1) = -1d30
      vmn(2) = 1d30
      vmx(2) = -1d30
      vnm(1) = 'x axis'
      vnm(2) = 'y axis'
c                                 read the data to get the range  
      do                             
         read (n4,*,iostat=ier) isym, x, y
         if (ier.ne.0) exit 
         if (x.lt.vmn(1)) vmn(1) = x
         if (x.gt.vmx(1)) vmx(1) = x
         if (y.gt.vmx(2)) vmx(2) = y
         if (y.lt.vmn(2)) vmn(2) = y
      end do 
c                                 get some options and
c                                 set up transformations
      call psaxop (1,jop0,iop1)
 
      call psipts 
 
      call psaxes (jop0)

      end

      subroutine psipts 
c---------------------------------------------------------------
c psipts - subprogram to x-y points.

      implicit none

      include 'perplex_parameters.h'

      integer ier, isym, jsym, icfg, icbg

      double precision r,x,y

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen 
c---------------------------------------------------------------  
      rewind (n4)

      do 

         read (n4,*,iostat=ier) isym,x,y
         if (ier.ne.0) exit

         r = 0.78d0

 
         if (isym.lt.4) then

            if (isym.eq.0) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,7,0,1)

            else if (isym.eq.1) then
 
               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,7)

            else if (isym.eq.2) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,1,0,1)

            else if (isym.eq.3) then

               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,1)

            end if 

         else 
 
            r = 0.38d0

            if (isym.eq.4) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,7,0,1)

            else if (isym.eq.5) then
 
               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,7)

            else if (isym.eq.6) then

               call pselip (x,y,r*dcx,r*dcy,1d0,0d0,1,0,1)

            else if (isym.eq.7) then

               call psrect (x-r*dcx,x+r*dcx,y-r*dcy,y+r*dcy,1d0,0d0,1)

            else 

               if (isym.gt.1000) then
                   r = 1.
                   jsym = isym - 1000
               else
                   r = 0.5
                   jsym = isym
               end if

               if (jsym.gt.24) jsym = 24

               if (jsym.le.12) then 
                  icfg = jsym
                  icbg = 0
               else 
                  icfg = jsym - 12
                  icbg = 1
               end if 

               call pselip (x,y,r*dcx,r*dcy,1d0,1d0,7,icfg,icbg)

            end if 

         end if 

      end do 
 
      end
