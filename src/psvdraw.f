 
c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy and
c Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c PSVDRAW - a program to generate postscript code for x-y plots and chemographic
c phase diagrams.  The input data format is consistent with the files generated
c by the programs VERTEX and FRENDLY.
 
c Please do not distribute any part of this source.
 
      PROGRAM PVDRAW
 
      implicit none

      include 'perplex_parameters.h'

      character yes*1

      integer ier 

      integer  iop0 
      common / basic /iop0

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iphct
      common/ ln4 /iphct

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 8
c                                 version info
      call vrsion (6)
c                                 default no modification prompts
      iop0 = 0
 
      do 
c                                 get input file 
         write (*,1000) 
          
         call readrt

         call mertxt (tfname,prject,'.plt',0)
         
         open (n4,iostat=ier,file=tfname,status='old')

         if (ier.ne.0) then
       
            write (*,1010) tfname
            read (*,'(a)') yes

            if (yes.eq.'Y'.or.yes.eq.'y') cycle 

            stop

         end if

         exit  

      end do        

      read (n4,*,iostat=ier) icopt
      if (ier.ne.0) call error (67,1d0,iop0,tfname)
      if (icopt.gt.3) call error (66,1d0,iop0,'PSVDRW')
c                                 read plot option file, set
c                                 default transformation
      call rdopt 
c                                 open output file 
      call psopen
  
      if (icopt.ne.0) then 
         write (*,1020) 
         read (*,'(a)') yes
         if (yes.eq.'y'.or.yes.eq.'Y') iop0 = 1
      end if 
  
      if (icopt.eq.1) then

         call psxypl 
 
      else if (icopt.eq.0) then
 
         call pschem 
 
      else if (icopt.eq.3) then

         call psmixd 

      else  
 
         call error (67,1d0,iop0,tfname)
 
      end if
 
      call psclos
 
      close (n4)

1000  format (/,'Enter the project or plot file name [i.e., without ',
     *       'the .plt suffix]:')
1010  format (/,'**warning ver191** cannot find file:',/,a,/,
     *       'run VERTEX, FRENDLY, SPECIES or PT2CURV to generate the ',
     *       'file or try a different name (y/n)?')
1020  format (/,'Modify the default plot (y/n)?')

      end
c---------------------------------------------------------------------
      subroutine psxypl 
 
c psxypl - subroutine to output x-y plot.

      implicit none

      logical err
 
      include 'perplex_parameters.h'

      integer jop0
c                                 rewind n4, cause icopt was already read
      rewind (n4)
c                                 read plot file header
      call plinp (err)
c                                 get user options and read
c                                 rest of plot file, draw data
      call psdplt (jop0)         
c                                 draw axes
      call psaxes (jop0)

      end

      subroutine psdplt (jop0)
c----------------------------------------------------------------------
c psdplt - draw xy-plot
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,lop(15),iop2,iop3,iop5,iop6,iop7,iop4c,iop4p,jop0,
     *        iend,iop1,jop4,iop4
 
      character yes*1, prompt*14

      double precision fac2,fac3

      integer  iop0 
      common / basic /iop0

      integer iphct
      common/ ln4 /iphct

      integer isat
      common/ wee /isat
c---------------------------------------------------------------
      iop2 = 0
      iop3 = 0
      iop5 = 0
      iop6 = 0
      iop7 = 0
c                                 jop4 is the mode of variance restriction, 0 = off, 1 < threshold, 2 = threshold
      jop4 = 0
c                                 iop4 is the threshold for variance restrictions
      iop4 = 99 
c                                 the threshold explicitly for curves and points
      iop4c = iop4
      iop4p = iop4 
 
      fac2 = 2d-1
      fac3 = 5d-2
c                                 iop2 = 1 suppress curve labels
c                                 iop3 = 1 suppress point labels                                 
      do i = 1, 15
         lop(i) = 0
      end do 
c                                 get some options and
c                                 set up transformations
      call psaxop (icopt,jop0,iop1)

      if (iop0.eq.0) goto 30
c                                 phase field restrictions by variance 
      write (*,1160)       
      read (*,'(a)') yes       
c                                 variance restrictions for projections
      if (yes.eq.'y'.or.yes.eq.'Y') then
c                                 get type of variance restriction (later called jop4)
         write (*,1090)
c                                 get choice
         call rdnumb (fac2,0d0,iop4,1,.false.)
         if (iop4.ne.2) iop4 = 1
c                                 get the variance threshold
         write (*,1100) 
c                                 get choice
         call rdnumb (fac2,0d0,jop4,2,.false.)
         if (jop4.lt.1.or.jop4.gt.99) jop4 = 2 

         iop4c = jop4
         iop4p = jop4

         if (iop4.eq.1.and.jop4.gt.2) then
c                                 allow suppression of pseudounivariant curves
            write (*,1170) 
            read (*,'(a)') yes 
c                                 general restriction on pseudounivariant curves
            if (yes.eq.'y'.or.yes.eq.'Y') iop4c = 2
         end if 

         if (iop4.eq.1.and.jop4.gt.1) then

            write (*,1180) 
            read (*,'(a)') yes 
c                                 general restriction on pseudoinvariant points
            if (yes.eq.'y'.or.yes.eq.'Y') iop4p = 1
         end if 
      end if  
c                                 restrict by assemblage
      write (*,1150) 
      read (*,'(a)') yes

      if (yes.eq.'y'.or.yes.eq.'Y') then
c                               if saturated components write warning:
         if (isat.ne.0) write (*,1400) 
         write (*,1140) 
         read (*,'(a)') yes
         if (yes.eq.'y'.or.yes.eq.'Y') then
            iop5 = 1
            prompt='present in the'
            call rname (1,prompt)
         end if
c                               restrict by phase presence
         write (*,1120)  
         read (*,'(a)') yes
         if (yes.eq.'y'.or.yes.eq.'Y') then 
            iop6 = 1
            prompt=' absent in all'
            call rname (2,prompt)
         end if 
c                               restrict by phase absence
         write (*,1130)
         read (*,'(a)') yes
         if (yes.eq.'y'.or.yes.eq.'Y') then 
            iop7 = 1
            prompt='present in all'
            call rname (3,prompt)
         end if 
      end if
c                                       equilibrium labeling
c                                       suppress/modify labels  
      write (*,1225)
      write (*,1230)
c                                       modify eq labeling
      read (*,'(a)') yes

      if (yes.eq.'y'.or.yes.eq.'Y') then  

            if (icopt.eq.1) then 
c                                       suppress curve labels?
               write (*,1220) 
            else  
               write (*,1215) 
            end if 
            read (*,'(a)') yes
            if (yes.eq.'y'.or.yes.eq.'Y') iop2 = iabs(iop2-1)
 
            if (iop2.eq.0) then
c                                       curve label modifications
               write (*,1070)
               read (*,'(a)') yes
 
               if (yes.eq.'y'.or.yes.eq.'Y') then
                  write (*,1080)
                  read (*,'(a)') yes
                  if (yes.eq.'y'.or.yes.eq.'Y') lop(15) = 1 

10                write (*,1030)
                  read (*,*) fac2
                  if (fac2.lt.0d0.or.fac2.gt.1d0) then
                     write (*,1200) fac2
                     goto 10
                  end if
20                write (*,1110) fac2
                  read (*,*) fac3
                  if (fac3.lt.0d0.or.fac3.gt.fac2) then
                     write (*,1200) fac3
                     goto 20
                  end if
               end if
            end if 
c                                      suppress ip point labels
            if (icopt.eq.1) then 
               write (*,1190)
            else
               write (*,1195)
            end if 
            read (*,'(a)') yes
            if (yes.eq.'y'.or.yes.eq.'Y') iop3 = iabs(iop3-1)

      end if     
c                                 draw univariant curves
30    call pscurv (iop1,iop2,iop4,iop4c,iop5,iop6,iop7,
     *             lop,fac2,fac3,iend)

      if (iend.eq.0)  
     *    call psipts (iop1,iop3,iop4,iop4p,iop5,iop6,iop7)
c                                 look for text labels
      call pslbtx 

      call psblrb (4)

1030  format (/,'Enter minimum fraction of the axes length that a',/
     *        'curve must be to receive a text label (0-1): ')
1070  format (/,'Change default labeling of curve segments (y/n)?')
1080  format (/,'Suppress labels of pseudounivariant curves (y/n)?')
1090  format (/,'Select true variance restriction to be applied to ',
     *          'pseudo-invariant/univariant equilibria:',/,
     *        '  1 - show all fields with true variance < than a ',
     *        'specified value [default]',/,
     *        '  2 - show all fields with a specified true variance ')
1100  format (/,'Enter the true variance to be used for this ',
     *          'restriction [1-99, default = 2]:')
1110  format (/,'Enter minimum fraction of the axes length that a',/
     *        'curve must be to receive a numeric label (0-',f5.3,'):')
1120  format (/,'Show only without phases (y/n)? ')
1130  format (/,'Show only with phases (y/n)? ')
1140  format (/,'Show only with assemblage (y/n)? ')
1150  format (/,'Restrict phase fields by phase identities (y/n)?',/,
     *        '  answer yes to:',/,
     *        '   - show fields that contain a specific assemblage',/,
     *        '   - show fields that do not contain specified phases',/,
     *        '   - show fields that contain any of a set of specified',
     *            ' phases ')
1160  format (/,'Restrict phase fields by variance (y/n)?'/,
     *          '  answer yes to:',/,
     *       '   - suppress pseudounivariant curves and/or',
     *       ' pseudoinvariant points of',/,
     *       '     a specified true variance.')
1170  format (/,'Suppress pseudounivariant curves (y/n)? ')
1180  format (/,'Suppress pseudoinvariant points (y/n)? ')
1190  format (/,'Suppress point labels (y/n)? ')
1195  format (/,'Show point labels (y/n)? ')
1200  format (/,g13.6,' is an invalid value.',/)
1215  format (/,'Show curve labels (y/n)? ')
1220  format (/,'Suppress curve labels (y/n)? ')
1225  format (/,'Modify default equilibrium labeling (y/n)?',/,
     *       '  answer yes to:')
1230  format ('   - modify/suppress [pseudo-] univariant curve labels',/
     *       ,'   - suppress [pseudo-] invariant point labels')
1400  format (/,'WARNING: You can not specify saturated phases or',
     *          ' phases determined by',/,'component saturation',
     *          ' constraints in these restrictions.',/)       

      end 
c----------------------------------------------------------------------
      subroutine pschem
 
c pschem - subroutine to output ternary chemographies.

      implicit none

      include 'perplex_parameters.h'
 
      character title*162, record*72, yes*1,xname(k5)*8
 
      double precision x3(3),y3(3),xx(j9),yy(j9),style,x1,y1,y,dyt,yt,xt
 
      integer iperm(2,3),i,j,iflag,isat,iop1,kvert,id,nchar

      integer icp,istct,ipoint,ifct,ipot,jas,jd

      logical vline, tlbl

      save x3,y3,dyt,yt,xt,iperm 

      data x3,y3,dyt,yt,xt,iperm/0d0,0.5d0,1d0,0d0,0.866025d0,0d0,
     *     3d2,98d1,1d2,1,2,1,3,2,3/

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer ikp
      common/ phase /ikp(k1)
      
      integer idf,ib,iasmbl,ivchk
      double precision x
      common/ asmbl /x(2,k1),idf(3,k2),ib,iasmbl(k2),ivchk(k1)
 
      integer iphct
      common/ ln4 /iphct
c----------------------------------------------------------------------
      iflag = 0
c                                  read header information 
      read (n4,*) icp,istct,iphct,ipoint,ifct,isat,ipot,isoct

      if (icp.ne.3) call error (64,xfac,iphct,'PSCHEM')
 
      if (iphct.gt.k1.or.jvar.gt.5.or.ib.gt.k2.or.ipot.gt.5) 
     *   call error (65,xfac,iphct,'PSCHEM') 

      read (n4,'(a)') (vnm(i),i=1,ipot)
      read (n4,'(a)') title
c                                  names
      read (n4,'(10a)') (names(i), i = 1, iphct)
c                                  phase coordinates
      read (n4,*) ((x(j,i),j = 1, 2), i=istct, iphct)
      read (n4,*) (ikp(i), i = 1, iphct)
c                                  solution names
      if (isoct.ne.0) read (n4,'(8a)') (fname(i),i = 1, isoct)
c                                  write graphics code names for x variables:
      read (n4,'(10a)') (xname(i), i = 1, icp)
c                                  -----------------------------------------------
c                                  convert to equilateral coordinates:
      do i = istct, iphct
         x(1,i) = x(1,i) + 0.5d0 * x(2,i)
         x(2,i) = x(2,i) * 0.866025d0
      end do 
c                                 default settings
c                                 suppress tie lines
      iop1 = 0
c                                 black fill for homogeneous fields
      fill = .true.
c                                 variable line thickness
      vline = .true.
c                                 label phases
      tlbl = .true.
c                                 modify default plot?
      write (*,1010) 
      read (*,'(a)') yes

      if (yes.eq.'y'.or.yes.eq.'Y') then
c                                 draw tie lines?
         write (*,1060)
         read (*,'(a)') yes

         if (yes.eq.'y'.or.yes.eq.'Y') then

            iop1 = 1
c                                 suppress homogeneous fields fill
            write (*,1030)
            read (*,'(a)') yes

            if (yes.eq.'y'.or.yes.eq.'Y') fill = .false. 
c                                 use one line thickness
            write (*,1070)
            read (*,'(a)') yes

            if (yes.eq.'y'.or.yes.eq.'Y') vline = .false.

         end if 

c                                 suppress phase labels?
         write (*,1080)
         read (*,'(a)') yes

         if (yes.eq.'y'.or.yes.eq.'Y') tlbl = .false.

      end if 

      do 
c                                  read simple variables
         read (n4,*,end=99) (vmn(i), i = 1, ipot)
         read (n4,*) ib
c                                  read phase configurations   
         read (n4,*) ((idf(j,i),j = 1, icp), i = 1, ib)
 
         do i = 1, iphct
            ivchk(i) = 0 
         end do 
c                                  read assemblage flags
         read (n4,*) (iasmbl(i), i = 1, ib)
c                                  read saturated phase id's
         if (isat.gt.0) then 
            read (n4,*) isat
            read (n4,*) (idss(i), i = 1, isat)
         end if 
 
         iflag = (iflag-1)**2
 
         if (iflag.eq.0) then
            xt = 375d0
         else
            yt = yt - dyt
            xt = 75d0
         end if
 
         call psssc2 (0d0,1d0,0d0,1d0)
 
         call psstrn (8d-2,8d-2,xt,yt,0d0)

         if (iop1.eq.0) then 
c                                  don't draw tielines
c                                  first divariant 
            do i = 1, ib                      
               if (iasmbl(i).eq.2) then
                  call merger (i,iperm,kvert,ipoint,xx,yy)
                  call pspygn (xx,yy,kvert,0d0,0d0,7)
               end if
            end do
c                                  third univariant
            do i = 1, ib                      
               if (iasmbl(i).eq.1) then
                  call merger (i,iperm,kvert,ipoint,xx,yy)
                  call pspygn (xx,yy,kvert,1d0,0d0,3)
               end if
            end do
c                                  last invariant
            do i = 1, ib                      
               if (iasmbl(i).eq.0) then
                  do j = 1, 3
                     id = idf(j,i)
                     ivchk(id) = 1
                     xx(j) = x(1,id) 
                     yy(j) = x(2,id) 
                  end do
                  call pspygn (xx,yy,3,1d0,0d0,0)
               end if
            end do   

         else 
c                                  old mode, draw tielines:
            do i = 1, ib
 
               jas = iasmbl(i)
 
               do j = 1, 3
                  id = idf(j,i)
                  if (jas.eq.0.or.id.le.ipoint) ivchk(id) = 1
                  xx(j) = x(1,id) 
                  yy(j) = x(2,id) 
               end do 
c                                           draw tie lines:
               if (jas.eq.0) then
c                                           invariant
                  if (vline) then
                     call pspygn  (xx,yy,3,1d0,1d0,0)
                  else 
                     call pspygn  (xx,yy,3,1d0,0d0,0)
                  end if 

               else if (jas.eq.1) then
c                                           univariant
                  if (vline) then
                     style = 1
                  else 
                     style = 0
                  end if

                  if (fill) then 

                     call pspygn  (xx,yy,3,0d0,0d0,3)
c                                 draw a thin line between same phase
                     do j = 1, 3
                        id = idf(iperm(1,j),i)
                        jd = idf(iperm(2,j),i)
c                                           draw the line:
                        if (ikp(jd).eq.ikp(id)) call psline 
     *                     (xx(iperm(1,j)),yy(iperm(1,j)),
     *                      xx(iperm(2,j)),yy(iperm(2,j)),style,0d0)
                     end do
                  else
                     call pspygn  (xx,yy,3,1d0,0d0,0)
                  end if 

               else if (jas.eq.2) then
c                                           divariant
                  if (fill) then 
                     call pspygn  (xx,yy,3,0d0,0d0,7)
                  else
                     call pspygn  (xx,yy,3,1d0,0d0,0)
                  end if 
               end if 

            end do 

         end if
 
         if (tlbl) then 
            do i = istct, iphct
               if (ivchk(i).ne.0) then
                  call pssctr (ifont,nscale,nscale,0d0)
                  x1 = x(1,i) + 15d-3
                  y1 = x(2,i) 
                  call pstext (x1,y1,names(i),8)
                  call pselip (x1-0.015d0,y1,7d-3,7d-3,0d0,0d0,7,0,1)
               end if
            end do
         end if 
c                                       draw bounding triangle
         call pspygn  (x3,y3,3,1d0,0d0,0)
c                                       sectioning constraints
         call pssctr (ifont,nscale,nscale,0d0)

         y = 1d0

         do i = 1, ipot
            write (record,1040) vnm(i),vmn(i)
            nchar = 0 
            call psublk (record,nchar)
            call pstext (0.75d0,y,record,nchar)
            y = y - 0.06d0*nscale
         end do 
 
         if (ifct.ne.0) then 
            call pstext (0.75d0,y,'(fluid saturated)',17)
            y = y - 0.06d0*nscale
         end if 

         if (isat.gt.0) then
            write (record,1050) (names(idss(i)),i=1,isat)
            nchar = 72 
            call psublk (record,nchar)
            call pstext (0.75d0,y,record,nchar)
         end if 
 
      end do 

1010  format ('Modify the default plot (y/n)?')
1030  format ('Suppress phase field fills (y/n)?')
1040  format (a8,'=',g9.3)
1050  format ('+ ',6(a8,' '))
1060  format ('Draw tielines (y/n)?')
1070  format ('Suppress variable line thickness (y/n)?')
1080  format ('Suppress phase labels (y/n)?')

99    end

      subroutine merger (i,iperm,kvert,ipoint,xx,yy)

      implicit none
 
      include 'perplex_parameters.h'
 
      double precision xx(j9), yy(j9)
 
      integer iperm(2,3),idv(j9),idp(3),jdv(3),jas,i,kvert,j,id,
     *                   kas,k,l,jvert,kp,kp2,ipoint,ll1,ll2,ll3

      integer idf,ib,iasmbl,ivchk
      double precision x
      common/ asmbl /x(2,k1),idf(3,k2),ib,iasmbl(k2),ivchk(k1)

      integer ikp
      common/ phase /ikp(k1)
c                                    merge high variance fields
c                                    load first part into polygon:
      jas = iasmbl(i)
      kvert = 3

      do j = 1, 3

         id = idf(j,i)
         idv(j) = id
c                                    identify the phases
         if (ikp(id).ne.0) then 

            idp(j) = -ikp(id)
c                                    tag used endmember compounds:
            if (id.le.ipoint) then
               ivchk(id) = 1
            end if 

         else 

            idp(j) = id
c                                    tag used compounds:
            ivchk(id) = 1

         end if 
      end do

      do 
c                                    now search for the common simplexes:
         do 20 j = 1, ib

            kas = iasmbl(j)
c                                    no match, reject:
            if (kas.ne.jas.or.j.eq.i) cycle

            do 30 k = 1, 3
               id = idf(k,j)
               jdv(k) = id
c                                    check if the phases match: 
               if (ikp(id).ne.0) then 
                  if (id.lt.ipoint) ivchk(id) = 1
                  do l = 1, 3
                     if (-ikp(id).eq.idp(l)) goto 30
                  end do 
c                                    no match, reject:
                  goto 20                    
               else 
                  do l = 1, 3
                     if (id.eq.idp(l)) goto 30
                  end do
c                                    no match, reject:
                  goto 20
               end if 
30          continue
c                                    ok, the simplex has the same assemblage
c                                    now match vertices:
            jvert = kvert 

            do k = 1, kvert
               if (k.lt.kvert) then 
                  kp = k + 1
               else
                  kp = 1
               end if 

               do l = 1, 3
                  if (jdv(iperm(1,l)).eq.idv(k).and.
     *                jdv(iperm(2,l)).eq.idv(kp).or.
     *                jdv(iperm(2,l)).eq.idv(k).and.
     *                jdv(iperm(1,l)).eq.idv(kp)) then
                     if (kp.eq.1) then
c                                    tag the common join on the end
                        idv(kvert+1) = jdv(4-l)
                     else 
c                                    the simplex might share the next
c                                    segment as well:
                        if (jas.gt.1) then 
                           if (kp.lt.kvert) then
                              kp2 = kp + 1
                           else 
                              kp2 = 1
                           end if 

                           do ll2 = 1, 3
                              if (jdv(iperm(1,ll2)).eq.idv(kp).and.
     *                            jdv(iperm(2,ll2)).eq.idv(kp2).or.
     *                            jdv(iperm(2,ll2)).eq.idv(kp).and.
     *                            jdv(iperm(1,ll2)).eq.idv(kp2)) then
c                                     the simplex matches two segments
c                                     remove the matched point
                                 do ll3 = kp+1, kvert
                                    idv(ll3-1) = idv(ll3) 
                                 end do 
                                 iasmbl(j) = -1
                                 jvert = jvert - 1
                                 goto 40 
                              end if 
                           end do 
                        end if 
c                                    found a common join, displace the old vertices
                        do ll1 = kvert + 1, kp, -1
                           if (ll1.gt.j9) call error (1,0d0,ll1,'J9')
                           idv(ll1) = idv(ll1-1)
                        end do
c                                    insert the new vertex
                        idv(kp) = jdv(4-l)
                     end if 
c                                    reset the assemblage flag:
                     iasmbl(j) = -1 
                     jvert = kvert + 1
                     goto 40 
                  end if 
               end do 
            end do 
20       continue
c                                    if got to here then there were no
c                                    no matches, save the polygon
         do j = 1, kvert
c                                    the polygon now has kvert vertices
c                                    load the vertices
            if (idv(j).eq.0) cycle

            xx(j) = x(1,idv(j))
            yy(j) = x(2,idv(j))

         end do 

         exit 
c                                    if here found a match, check the 
c                                    list for more
40       kvert = jvert 

      end do 

      end 
c----------------------------------------------------------------------
      subroutine pscurv (iop1,iop2,jop4,iop4,iop5,iop6,iop7,
     *                   lop,fac2,fac3,iend)

c pscurv - subroutine to output curves.

      implicit none
 
      include 'perplex_parameters.h'

      integer k,lchk(k2),lop(15),idv(k7),iend,ipt,ird,i,
     *        ivar,ipr,iop1,iop2,iop4,iop5,iop6,iop7,jvct,imatch,jpr,
     *        kvar,jop4
 
      double precision x(l5),y(l5),rline,yfrac,xfrac,fac2,fac3,cwidth

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      integer ivct,iplus,iminus,idr
      double precision vnu
      common/ rxn /vnu(k7),idr(k7),ivct,iplus(k7),iminus(k7)

      integer ixct,iex,jex,ict
      common/ excl1 /ixct(3),iex(50,3),jex(50,3),ict(3)

      character*10 xnams
      common/ excl4 /xnams(50,3)
c----------------------------------------------------------------------
      iend = 0

      do i = 1, 3
         ict(i) = 0
      end do 

      do 

         read (n4,*,end=100) ipt,ird,ivar,ivct,(idr(i),i=1,ivct)
         read (n4,*) (vnu(i),i=1,ivct)
c                                 the end of curve data is indicated by 
c                                 a ordinate counter (ipt) of 1. 
         if (ipt.eq.1) then
            goto 99
         else if (ipt.eq.0) then
            cycle 
         end if

         ipr = ipt/2

         if (ird.gt.k2) call error (1,xfac,ird,'K2 (PSCURV)') 
         if (ipr.gt.l5*2) call error (1,xfac,l5,'L5 (PSCURV)')
 
         read (n4,*) (x(i),y(i),i=1,ipr)
c                                 do NaNchk
         do i = 1, ipr
            call nanchk (x(i),y(i),'PSCURV/PSVDRAW')
         end do
c                                 variance restrictions:
         if (jop4.ne.0) then 
            if (jop4.eq.1.and.ivar.gt.iop4-1) then 
               cycle 
            else if (jop4.eq.2.and.ivar.ne.iop4) then 
               cycle 
            end if 
         end if 
c                                 copy phase id's into idv
c                                 array for phase restrictions:
         do i = 1, ivct 
            idv(i) = idr(i)
         end do 

         jvct = ivct  

         call checkr (iop5,iop6,iop7,idv,k7,jvct,imatch)

         if (imatch.eq.1) cycle 

         if (iop6.eq.1) ict(2) = ict(2) + 1
 
         if (iop1.eq.1) then
            jpr = 0
            do i = 1, ipr
               if (x(i).gt.xmax.or.x(i).lt.xmin.or.
     *             y(i).gt.ymax.or.y(i).lt.ymin) cycle
               jpr = jpr + 1
               x(jpr) = x(i)
               y(jpr) = y(i)
            end do 
            if (jpr.eq.0) cycle 
            ipr = jpr
         end if

         kvar = ivar 
         cwidth = 0
 
         if (ivar.eq.1) then
            cwidth = 2
         else if (ivar.eq.98) then
            ivar = 1
            cwidth = 2
         else if (ivar.eq.99) then
            ivar = 1
         end if
 
         rline = dfloat(ivar)
         cwidth = cwidth * width
 
         if (spline) then
            call psbspl (x,y,ipr,rline,cwidth,0)
         else 
            call pspyln (x,y,ipr,rline,cwidth,0)
         end if
 
         if (ird.eq.0) cycle 
 
         if (iop2.eq.0.and.lchk(ird).eq.0) then

            if (lop(15).eq.1.and.kvar.ne.1) cycle 
 
            yfrac = dabs(y(ipr)-y(1)) / ylen
            xfrac = dabs(x(ipr)-x(1)) / xlen

            if (xfrac.lt.fac3.and.yfrac.lt.fac3) then
c                                     curve is too short
               cycle 
            else if (xfrac.lt.fac2.and.yfrac.lt.fac2) then
c                                     use a numeric label
               call psalbl (x,y,ipr,ivar,ird,0,1d0,0)
            else
c                                     use a text label
               if (ivct.ne.0) 
     *            call psalbl (x,y,ipr,ivar,ird,0,1d0,1)
            end if

            lchk(ird) = 1
 
         end if
 
      end do 

100   write (*,1010) 
      iend = 1
 
99    if (iop5.eq.1) write (6,*) ict(1),
     *               ' curves have the assemblage: ',
     *               (xnams(k,1),' ',k = 1, ixct(1))

      if (iop6.eq.1) write (6,*) ict(2),
     *               ' curves have none of the phases: ',
     *               (xnams(k,2),' ',k = 1, ixct(2))

      if (iop7.eq.1) write (6,*) ict(3),
     *               ' curves have one of the phases: ',
     *               (xnams(k,3),' ',k = 1, ixct(3))

1010  format ('End-of-file! if this is a VERTEX plot file then',
     *        ' VERTEX terminated incorrectly',/)
      end

c---------------------------------------------------------------
      subroutine plumin (ip,im)

c subprogram to order a reactants according to sign)

      implicit none

      include 'perplex_parameters.h'

      integer ip,im,i

      double precision vp(k7), vm(k7)

      integer ivct,iplus,iminus,idr
      double precision vnu
      common/ rxn /vnu(k7),idr(k7),ivct,iplus(k7),iminus(k7)

      ip = 0
      im = 0
c                             load id's of phases with + or - coeffs
      do i = 1, ivct
         if (vnu(i).gt.0d0) then
            ip = ip + 1
            vp(ip) = vnu(i)
            iplus(ip) = idr(i)
         else
            im = im + 1
            vm(im) = vnu(i)
            iminus(im) = idr(i)
         end if
      end do 

      do i = 1, im
         vnu(i) = vm(i)
         idr(i) = iminus(i)
      end do 

      do i = 1, ip
         vnu(i+im) = vp(i)
         idr(i+im) = iplus(i)
      end do 

      end
c---------------------------------------------------------------
      subroutine rxntxt (string,iend)
 
c subprogram to write a text label for a reaction
 
      implicit none
 
      include 'perplex_parameters.h'

      integer ip,im,i,j,is,ist,id,iend
 
      character text(400)*1, string*(*)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer ikp
      common/ phase /ikp(k1)

      integer ivct,iplus,iminus,idr
      double precision vnu
      common/ rxn /vnu(k7),idr(k7),ivct,iplus(k7),iminus(k7)

c                             load id's of phases with + or - coeffs
      call plumin (ip,im)
c                             now dump the names into the array text
      ist = 1
      is = 1

80    do i = is, im
         id = idr(i)
         if (ikp(id).eq.0) then
c                             simple compound:
            read (names(id),'(8a)') (text(j), j = ist, ist + 7)
            text(ist+8) = ' '
            ist = ist + 9
         else
c                             solution phases:
            read (fname(ikp(id)),'(10a)') (text(j), j = ist, ist + 9)
            text(ist+10) = '('
            read (names(id),'(8a)') (text(j), j = ist + 11, ist + 18)
            text(ist+19) = ')'
            text(ist+20) = ' '
            ist = ist + 21
         end if
      end do 
 
      if (is.eq.1) then
         text(ist) = '='
         text(ist+1) = ' '
         ist = ist + 2
         is = im + 1
         im = ivct
         goto 80
      else
         text(ist) = ' '
      end if
c                             now filter out blanks
      iend = 1

      do i = 2, ist
         if (text(i).eq.' '.and.text(i+1).eq.' ') then
            cycle
         else if (text(i).eq.' '.and.text(i+1).eq.')') then
            cycle
         else if (text(i).eq.' '.and.text(i+1).eq.'(') then
            cycle
         end if
         iend = iend + 1
         text(iend) = text(i)
      end do 

      if (iend.gt.400) iend = 400
      write (string,'(400a)') (text(i),i=1,400)
 
      end

      subroutine psipts (iop1,iop3,jop4,iop4,iop5,iop6,iop7)
c---------------------------------------------------------------------- 
c psipts - subprogram to output invariant points.

      implicit none

      include 'perplex_parameters.h'

      double precision x(1),y(1),r,rline
 
      integer ipds(k8),k,iop1,iop3,iop4,iop5,iop6,iop7,iop9,i,j,ipct,
     *        incts,ipid,ipvar,imatch,jop4,ier

      integer ixct,iex,jex,ict
      common/ excl1 /ixct(3),iex(50,3),jex(50,3),ict(3)

      character*10 xnams
      common/ excl4 /xnams(50,3)

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      data iop9/0/
c----------------------------------------------------------------------
      do i = 1, 3
        ict(i) = 0
      end do 
 
      read (n4,*) ipct, incts
 
      if (ipct.eq.0) goto 99
 
      do i = 1, ipct
 
         read (n4,*,iostat=ier) ipid,ipvar,(ipds(j), j=1,incts),
     *                          (var(j), j=1,jvar)

         if (ier.ne.0) then
c                                 who know's how this condition comes about,
c                                 likely the user edited the file.
            call warn (99,0d0,0,'the list of invariant points in the '
     *             //'plt file is incomplete.')

            exit

         end if 

         x(1) = var(1)
         y(1) = var(2)

         call nanchk (x(1),y(1),'PSCURV/PSVDRAW')

c                                 variance restrictions:
         if (jop4.ne.0) then 
            if (jop4.eq.1.and.ipvar.gt.iop4-1) then 
               cycle
            else if (jop4.eq.2.and.ipvar.ne.iop4) then 
               cycle
            end if 
         end if 
 
         if (iop1.eq.1) then
            if (x(1).gt.xmax.or.x(1).lt.xmin.or.
     *          y(1).gt.ymax.or.y(1).lt.ymin) cycle
          end if
 
         call checkr (iop5,iop6,iop7,ipds,k8,incts,imatch)
         if (imatch.eq.1) cycle 

         if (iop6.eq.1) ict(2) = ict(2) + 1 
c                                 make radius proportional to variance 
         r = .78d0/(ipvar+1)
 
         call pselip (x(1),y(1),r*dcx,r*dcy,0d0,0d0,7,0,1)
 
         if (iop3.eq.0) call psalbl (x,y,1,ipvar,ipid,1,rline,iop9)
 
      end do 
 
      if (iop5.eq.1) write (6,*) ict(1),
     *               ' points have the assemblage: ',
     *               (xnams(k,1),' ',k = 1, ixct(1))

      if (iop6.eq.1) write (6,*) ict(2),
     *               ' points do not have any of the phases: ',
     *               (xnams(k,2),' ',k = 1, ixct(2))

      if (iop7.eq.1) write (6,*) ict(3),
     *               ' points have one of the phases: ',
     *               (xnams(k,3),' ',k = 1, ixct(3))
 
99    end

c-------------------------------------------------------------------------
      subroutine psmixd 

      implicit none
 
c psmixd - subroutine to draw binary mixed variable diagrams

      include 'perplex_parameters.h'

      integer icp,ipoint,ifct,isat,ipot,i,ird,ivar
  
      character*8 title*162, string*(lchar), y*1, tname(5),xname(k5)
  
      integer idf(3),jphi(k1),igo,jop0,iop1,iop2,iop3,jb,
     *        jplus,jminus,isum,idif,j,i1,id1,it,i2,id2,jt,itot,
     *        i00,imis,itoc,j1,iend

      double precision v(l2),rx,ry,t0,tlast,t1,rline,ttext

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer ivct,iplus,iminus,idr
      double precision vnu
      common/ rxn /vnu(k7),idr(k7),ivct,iplus(k7),iminus(k7)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer ikp
      common/ phase /ikp(k1)

      integer iphi,ivph,ib,istct
      double precision x
      common/ tx /x(k1),iphi(k1),ivph(k1),ib,istct

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen

      integer  iop0 
      common / basic /iop0

      integer iphct
      common/ ln4 /iphct

      data idf,igo/4*0/
c                                  start-of-header
c                                  ------------------------------
c                                  read simple variables
      read (n4,*) icp,istct,iphct,ipoint,ifct,isat,ipot,isoct

      if (icp.ne.2) call error (64,xfac,ifont,'PSMIXD')
  
      if (ib.gt.k1.or.ipot.gt.5) call error (65,xfac,ifont,'PSMIXD')

      if (isoct.ne.0) write (*,1060)
 
      read (n4,'(a)') (tname(i), i = 1, ipot)
      read (n4,'(a)') title
c                                 names
      read (n4,'(10a)') (names(i),i=1,iphct)
c                                 phase coordinates
      read (n4,*) (x(i), i = istct, iphct)
c                                 read solution identifier flags
      read (n4,*) (ikp(i),i=1,iphct)
c                                 read solution names
      if (isoct.ne.0) read (n4,'(8a)') (fname(i), i=1, isoct)

      read (n4,'(10a)') (xname(i),i = 1, icp)
c                                 end-of-header
c                                 ------------------------------
      read (n4,*) (v(i), i = 1, ipot)
c                                 read stable phase ids
      read (n4,*) ib
      read (n4,*) (iphi(i), i = 1, ib)
c                                 read saturated phase id's
      if (isat.gt.0) then 
         read (n4,*) isat
         read (n4,*) (idss(i),i=1,isat)
      end if 
c                                 read index of y variable
      read (tname(1),'(a)') vnm(2)
      vnm(1) = xname(2)

      read (n4,*) vmn(2), vmx(2) 

      vmn(1) = 0d0
      vmx(1) = 1d0
c                                 move sectioning variables
c                                 inot vnm array
      do i = 2, ipot
         vnm(i+1) = tname(i)
         vmn(i+1) = v(i)
      end do 

      jvar = ipot + 1

      call psaxop (icopt,jop0,iop1)
c                              get some other options
      iop1 = 1
      iop2 = 0
      iop3 = 0

      rx = dcx/4d0/xfac
      ry = dcy/4d0

      if  (iop0.eq.1) then 
         write (*,2050)
         read (*,'(a)') y
         if (y.eq.'y'.or.y.eq.'Y') iop2 = 1

         write (*,2030)
         read (*,'(a)') y
         if (y.eq.'y'.or.y.eq.'Y') then
            iop1 = 0
         else 
            write (*,2040)
            read (*,'(a)') y
            if (y.eq.'y'.or.y.eq.'Y') iop1 = 2
         end if 
         if (isoct.gt.0) then
            write (*,2060)
            read (*,'(a)') y
            if (y.eq.'y'.or.y.eq.'Y') iop3 = 1
         end if 
      end if 

      t0 = vmn(2)
      tlast = 0d0    
c                              get the variance of each compound:
      do i = 1, ib
         call getva1 (i,iop3)
      end do 

      do 
         read (n4,*) ird, ivct, ivar, t1, (idr(i), i = 1, ivct),
     *                                 (vnu(i), i = 1, ivct)

         if (ivct.eq.1) goto 90

         do i = 1, ib 
            jphi(i) = iphi(i)
         end do 

         jb = ib 

         call plumin (jplus,jminus)

         do i = 2, ib - 1
            rline = ivph(i)
            if (rline.ne.0d0) call psline 
     *                     (x(iphi(i)),t0,x(iphi(i)),t1,rline,0d0)
         end do 

         isum = 0
         idif = 0

         do 50 i = 1, ivct
            if (ikp(idr(i)).ne.0) then
               isum = isum + 1
               do j = 1, idif
                  if (idf(j).eq.ikp(idr(i))) goto 50
               end do 
               idif = idif + 1
               idf(idif) = ikp(idr(i))
            end if 
50       continue 
c                              consider different cases:
         if (ivct.eq.2) then
c                              get index (id) of reacted phase:
            do i = 1, ib
               if (iphi(i).eq.idr(1)) then
                  i1 = i - 1 
                  if (i1.lt.1) i1 = 1
                  id1 = iphi(i1)
                  it = i
                  i2 = i + 1
                  if (i2.gt.ib) i2 = ib 
                  id2 = iphi(i2)   
                  goto 55
               end if
            end do 

            write (*,*) 'nope'
            write (*,*) (iphi(i),i=1,ib)
            write (*,*) idr(1),ird,t1

            cycle 
c                              check if only a compositional
c                              change in a solution:
55          if (isum.eq.2.and.ikp(idr(1)).eq.ikp(idr(2))) then
c                              draw an open dot, but don't
c                              do anything else:         
               call pselip (x(idr(1)),t1,rx,ry,1d0,0d0,0,0,1)
               iphi(it) = idr(2)          
            
               cycle 

            end if  
c                              draw dot         
            call pselip (x(idr(1)),t1,rx,ry,1d0,0d0,7,0,1)
c                              draw line
            call psline (x(id1),t1,x(id2),t1,1d0,0d0)
c                              replace reacted phase:
            iphi(it) = idr(2)
            call getva1 (it,iop3)
            if (it.gt.1) call getva1 (it-1,iop3)
            if (it.lt.ib) call getva1 (it+1,iop3)

         else  
c                              classify reaction, jminus = 2 
c                              eutectoid, else peritectoid.
            if (jminus.eq.2) then
    
               i1 = idr(1)
               i2 = idr(2)
               i00 = idr(3)
c                              find indices of old phases:
               do i = 1, ib
                  if (iphi(i).eq.i1) id1 = i
                  if (iphi(i).eq.i2) id2 = i 
               end do 
                  
            else 

               i00 = idr(1)
               i1 = idr(2)
               i2 = idr(3)
               id1 = 0
               id2 = 0 
c                                 find indices of old phases:
               do i = 1, ib
                  if (iphi(i).eq.i1) id1 = i
                  if (iphi(i).eq.i2) id2 = i
               end do 

               if (id1.eq.0.or.id2.eq.0) then 
c                                 the phase with vnu < 0 is NOT
c                                 the new phase.
                  write (*,*) 'nope'
                  write (*,*) (iphi(i),i=1,ib)
                  write (*,*) idr(1), ird,t1
                  
                  cycle 

               end if 
            end if 

            if (id1.gt.id2) then
               it = id1
               jt = i1
               i1 = i2
               i2 = jt
               id1 = id2
               id2 = it
               if (jminus.eq.2) then
                  idr(1) = i1
                  idr(2) = i2
               else 
                  idr(2) = i1
                  idr(3) = i2
               end if
            end if

            if (isum.lt.2.or.idif.eq.3) then 
c                              all phases different
c                              draw line
               call psline (x(i1),t1,x(i2),t1,1d0,0d0)
               call pselip (x(i00),t1,rx,ry,1d0,0d0,7,0,1)

            else if (isum.eq.2) then
               if (idif.eq.2) then 
c                              all phases different
                  call psline (x(i1),t1,x(i2),t1,1d0,0d0)
                  call pselip (x(i00),t1,rx,ry,1d0,0d0,7,0,1)

               else 
c                              one solution + one compound
                  if (ikp(i00).eq.0) then
c                              compound in the middle
c                              change variance of reaction
                     ivar = 1

                     call pselip (x(i00),t1,rx,ry,1d0,0d0,7,0,1)

                     call psline (x(i1),t1,x(i2),t1,1d0,0d0)          

                  else 
c                              pseudoinvariant
                     if (ikp(i00).eq.0) then
                        call psline (x(i1),t1,x(i2),t1,1d0,0d0)
                     else if (ikp(i1).eq.0) then
                        call psline (x(i00),t1,x(i2),t1,1d0,0d0)
                     else 
                        call psline (x(i00),t1,x(i1),t1,1d0,0d0)
                     end if 
                  end if
               end if

            else if (isum.eq.3.and.idif.eq.1.and.igo.eq.0) then
c                              reaction on solvus
c                              check if on left line
               itot = 0
               call miscib (x(i1),x(i00),ikp(i00),imis,igo)
               if (imis.eq.0) call psline (x(i00),t1,x(i1),t1,1d0,0d0)
               itot = itot + imis
               call miscib (x(i00),x(i2),ikp(i00),imis,igo)
               itot = itot + imis
               if (imis.eq.0) call psline (x(i00),t1,x(i2),t1,1d0,0d0)
               if (itot.eq.0) then
                  ivar = 1
                  call pselip (x(i00),t1,rx,ry,1d0,0d0,7,0,1)
               end if 


            else 
c                              two solutions present.
c                              solvus possible, find
c                              the pseudocompounds of
c                              the same solution:
               if (ikp(i00).ne.ikp(i1).and.ikp(i00).ne.ikp(i2)) then
c                              i1 and i2 are the same:
                  call miscib (x(i1),x(i2),ikp(i1),imis,igo)
               else if (ikp(i00).eq.ikp(i1)) then
c                              i00 and i1 are the same:
                  call miscib (x(i1),x(i00),ikp(i00),imis,igo)
               else 
c                              i00 and i2 are the same:
                  call miscib (x(i00),x(i2),ikp(i00),imis,igo)
               end if 

               if (imis.eq.1) then
                     ivar = 1
c                              eutectic/peritectic with solvus
                     call psline (x(i1),t1,x(i2),t1,1d0,0d0)
                     call pselip (x(i00),t1,rx,ry,1d0,0d0,7,0,1)
               else if (ikp(i1).eq.ikp(i00)) then
c                              pseudounivariant, left
                  call psline (x(i1),t1,x(i00),t1,1d0,0d0)

               else if (ikp(i2).eq.ikp(i00)) then
c                              pseudounivariant, right
                  call psline (x(i2),t1,x(i00),t1,1d0,0d0)

               else
c                              else center
c                              compound in the middle
c                              change variance:
                  itoc = 1
                  do j = 1, ib
                     if (iphi(j).eq.i00) cycle
                     if (ikp(iphi(j)).eq.ikp(i00)) itoc = 0
                  end do 

                  if (itoc.eq.1) ivar = 1

                  if (ivar.eq.1) call pselip 
     *                           (x(idr(1)),t1,rx,ry,1d0,0d0,7,0,1)

                  call psline (x(iphi(id1)),t1,
     *                         x(iphi(id2)),t1,1d0,0d0)

               end if 
            end if 
c                              get new chemography:
            if (jminus.eq.2) then
c                              eutectoid, number of phases
c                              increases.
               do i = ib, id2, -1
                  iphi(i+1) = iphi(i)
                  ivph(i+1) = ivph(i)
               end do 
               ib = ib + 1
               iphi(id1+1) = i00
               do i = id1, id1 + 2
                  call getva1 (i,iop3)
               end do 
            else
c                              peritectoid, number of phases
c                              decreases.
               do i = id2, ib
                  iphi(i-1) = iphi(i)
                  ivph(i-1) = ivph(i)
               end do 
               ib = ib - 1
               do i = id1, id1 + 1
                  call getva1 (i,iop3)
               end do 
            end if            
         end if

         if (iop2.eq.1) then 
            do 170 i = 1, ib - 1
               if (ikp(iphi(i)).eq.0) cycle
               i1 = i + 1
               if (ikp(iphi(i)).ne.ikp(iphi(i1))) cycle
               do j = 1, jb - 1
                  if (x(jphi(j)).gt.x(iphi(i))) then
                     goto 170
                  else if (x(jphi(j)).eq.x(iphi(i))) then
                     j1 = j + 1
                     if (ikp(jphi(j)).eq.ikp(iphi(i)).and.
     *                   ikp(jphi(j)).eq.ikp(jphi(j1))) then
c                                  draw a rectangle from
c                                  i-i1-j-j1
                         call psrect (x(iphi(i)),x(iphi(i1)),
     *                                t0,t1,0d0,width,2)
                     end if
                  end if
               end do 
170         continue 
         end if 

         t0 = t1 

         if (iop1.eq.0) then 
            cycle
         else if (iop1.eq.1) then
            if (ivar.ne.1) cycle  
            call rxntxt (string,iend)
         else 
            call rxntxt (string,iend)
         end if 

         call pssctr (ifont,ascale,ascale,0d0)
         ttext = t1 + 1.25d0*dcy*ascale 
         if (ttext.lt.tlast+1.25d0*dcy*ascale) 
     *       ttext = tlast + 1.25d0*dcy*ascale
         call pstext (xmax+2d0*dcx,ttext,string,iend)
         tlast = ttext

      end do          

90    t1 = vmx(2)
c                                   fill in last part of one phase regions
      if (iop2.eq.1) then 
         do i = 1, ib - 1
            if (ikp(iphi(i)).eq.0) cycle
            i1 = i + 1
            if (ikp(iphi(i)).ne.ikp(iphi(i1))) cycle
c                                  draw a rectangle from
c                                  i-i1-j-j1
                      call psrect (x(iphi(i)),x(iphi(i1)),
     *                             t0,t1,0d0,width,2)
         end do 
      end if 

      do i = 2, ib - 1
         rline = ivph(i)
         if (rline.eq.0d0) cycle
         call psline (x(iphi(i)),t0,x(iphi(i)),t1,rline,0d0)
      end do 

      call psaxes (jop0)

1060  format (/,'WARNING: Psvdraw may not draw t-x diagrams',
     *          ' correctly if immiscibility',/,'occurs in a',
     *          ' projected ternary or higher order solution, to',
     *          ' avoid',/,'problems turn the miscibility test',
     *          ' off (modify default options).',/)
2030  format (/,'Suppress reaction labels (y/n)?')
2040  format (/,'Show pseudounivariant reaction labels (y/n)')
2050  format (/,'Try to fill one phase fields (y/n)?')
2060  format (/,'Skip tests for immiscibility (y/n)?')
 
      end
c-------------------------------------------------------------------
      subroutine getva1 (i,igo)
 
c getva1 - subroutine to get variance of compound i in a binary row

c     ivphi set 1 if true compound
c     ivphi set 0 if pseudocompound
c     ivphi set 2 if solvus limb and igo = 0.
c     if igo = 1 skip miscibility test

      implicit none
  
      include 'perplex_parameters.h'

      integer i,igo,i1,i2,id,id1,id2,imis
  
      integer ikp
      common/ phase /ikp(k1)      

      integer iphi,ivph,ib,istct
      double precision x
      common/ tx /x(k1),iphi(k1),ivph(k1),ib,istct

      if (i.eq.1.or.i.eq.ib) then 
         ivph(i) = 1
         return
      end if 

         i1 = i - 1
         i2 = i + 1

         id1 = iphi(i1)        
         id = iphi(i)
         id2 = iphi(i2)

         if (ikp(id).eq.0) then 
c                              true compound
            ivph(i) = 1

         else if (ikp(id1).eq.ikp(id).and.ikp(id2).eq.ikp(id)) then
c                              only one solution:
            ivph(i) = 0
            call miscib (x(id1), x(id), ikp(id), imis, igo)
            if (imis.eq.1) ivph(i) = 1
            call miscib (x(id), x(id2), ikp(id), imis, igo)
            if (imis.eq.1) ivph(i) = 1
           
         else 

            ivph(i) = 1

         end if 

      end 
c-------------------------------------------------------------------
      subroutine miscib (x1, x2, ids, imis, igo)

c miscib - subroutine check for immiscibility in a binary phase

c imis = 1 => solvus between solution(ids) compositions x1 and x2

      implicit none
  
      include 'perplex_parameters.h'

      integer imis,igo,i,ids

      double precision x1,x2
  
      integer ikp
      common/ phase /ikp(k1)

      integer iphi,ivph,ib,istct
      double precision x
      common/ tx /x(k1),iphi(k1),ivph(k1),ib,istct

      integer iphct
      common/ ln4 /iphct

      imis = 0
c                                igo=1, disabled cause of projections
      if (igo.ne.1) then 

         do i = istct, iphct

            if (ikp(i).ne.ids) cycle

            if (x(i).gt.x1.and.x(i).lt.x2) then 
               imis = 1
               exit 
            end if 

         end do 

      end if 

      end 

      subroutine checkr (iop5,iop6,iop7,id,n,ivct,imatch)

c checkr  - subroutine to check exclusions by phase identity
c           if imatch = 0 the assemblage id matches the criteria

      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,iop5,iop6,iop7,n,ivct,imatch,itic,itis

      integer id(n), kex(k8)

      integer ixct,iex,jex,ict
      common/ excl1 /ixct(3),iex(50,3),jex(50,3),ict(3)

      integer iphct
      common/ ln4 /iphct

      imatch = 0

      if (iop5.eq.1.or.iop6.eq.1.or.iop7.eq.1) then  

         imatch = 1

         if (iop5.eq.1) then       
            if (ivct.lt.ixct(1)) return  
c                            show only with assemblage:
            itic = 0
            do 2 i = 1, ivct
c                            if you get a match then look at
c                            look for the next phase
               call checki (1,id(i),kex(i))
                  if (kex(i).ge.0) then
 
                     if (isoct.gt.0.and.kex(i).ne.0.and.
     *                              i.gt.1.and.itic.gt.0) then
                        do k = 1, i - 1
                           if (kex(k).eq.kex(i)) goto 2
                        end do 
                     end if
                     itic = itic + 1
                  end if
2           continue
c                             itic must = ixct if all phases match
            if (itic.lt.ixct(1))  return
c                             if here the assemblage contains
c                             all the phases requested:
            ict(1) = ict(1) + 1
         end if

         if (iop6.eq.1) then
c                             reject fields that contain a phase:
            do j = 1, ivct
               call checki (2,id(j),itis)
c                            name matches a phase in the field
               if (itis.ge.0) return
 
            end do
c                             no match
         end if

         if (iop7.eq.1) then
            do j = 1, ivct
               call checki (3,id(j),itis)
               if (itis.ge.0) then
c                            name matches a phase in the field
                  ict(3) = ict(3) + 1
                  goto 90
               end if
            end do
c                            no match
            return
         end if

90       imatch = 0 

      end if

      end 

      subroutine plinp (err)
c---------------------------------------------------------------------- 
c plinp - subroutine to read x-y plot file header.

      implicit none
 
      include 'perplex_parameters.h'

      integer i,ier,jpt(2),jv(l3),iind,idep

      logical err

      double precision c0,c1,c2,c3,c4

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer ikp
      common/ phase /ikp(k1)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character*162 title
      common/ csta8 /title(4)

      integer iphct
      common/ ln4 /iphct

      integer isat
      common/ wee /isat
c----------------------------------------------------------------------
      read (n4,*,iostat=ier) icopt
      if (ier.ne.0) call error (67,c0,i,'PLINP')
      if (icopt.gt.3) call error (66,c0,i,'PLINP')

      read (n4,*) iphct, isoct
      read (n4,*) isat

      if (iphct.gt.0) then
         if (iphct.gt.k1) call error (65,c0,i,'PLINP')
         read (n4,'(10a)') (names(i),i = 1, iphct)
         read (n4,*) (ikp(i), i = 1, iphct)
      end if

      if (isoct.ne.0) read (n4,'(8a)') (fname(i), i = 1, isoct) 

      read (n4,'(a)') title
      read (n4,*) jvar,(jv(i),i=1,jvar),jpt
      read (n4,*) iind,idep,c0,c1,c2,c3,c4
      read (n4,*) (vmx(i),vmn(i),i=1,jvar)
      read (n4,'(a)') (vnm(i),i=1,jvar)

      end
c----------------------------------------------------------------------
      subroutine psalbl (x,y,ipr,ivar,ird,ityp,rline,iop9)
 
c psalbl - subroutine to output reaction(ityp=0) or ip(ityp=1)
c          name label.

      implicit none

      include 'perplex_parameters.h'

      character*6 lnms, pnms, string*(lchar)

      integer ipr,ivar,ird,ityp,iop9,imid,i,iend

      double precision fac,dx,dy,delta,xmid,ymid,theta,clx,cly,
     *                 rline,x(ipr),y(ipr)

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen 
c----------------------------------------------------------------------
 
      if (ivar.eq.1) then
            fac = ascale
      else
            fac = 0.75d0 * ascale
      end if
 
      if (ityp.eq.0) then
c                                       get curve midpoint
         dx = (x(ipr) - x(1))
         dy = (y(ipr) - y(1))
         delta = 1d30
 
         if (abs(dx)/xlen.gt.abs(dy)/ylen) then
            xmid = x(1) + dx / 2d0
 
            do i = 2, ipr
               if (abs(x(i)-xmid).lt.delta) then
                  delta = abs(x(i)-xmid)
                  imid = i
               end if
            end do 
 
         else
            ymid = y(1) + dy / 2d0
            do i = 2, ipr
               if (abs(y(i)-ymid).lt.delta) then
                  delta = abs(y(i)-ymid)
                  imid = i
               end if
            end do 
         end if
 
         if (iop9.eq.0) then
c                                   numeric label:
             if (ird.lt.100000) then
                write (lnms,'(i5)') ird
             else 
                write (lnms,'(i6)') ird
             end if

             call pssctr (ifont,fac,fac,0d0)

             if (ird.lt.10) then
                call pselip (x(imid),y(imid),
     *                      fac*dcx,fac*1.2d0*dcy,rline,0d0,1,0,1)
                call pstext (x(imid)-2.45d0*dcx*fac,
     *                      y(imid)+0.8d0*dcy*fac,lnms,6)
             else if (ird.lt.100) then
                call pselip (x(imid),y(imid),fac*1.74d0*dcx,
     *                       fac*1.6d0*dcy,rline,0d0,1,0,1)
                call pstext (x(imid)-2.71d0*dcx*fac,
     *                      y(imid)+0.8d0*dcy*fac,lnms,6)
             else if (ird.lt.1000) then
                call pselip (x(imid),y(imid),
     *               fac*2.75d0*dcx,fac*1.83d0*dcy,rline,0d0,1,0,1)
                call pstext (x(imid)-3.04d0*dcx*fac,
     *                       y(imid)+0.8d0*dcy*fac,lnms,6)
             else if (ird.lt.10000) then 
                call pselip (x(imid),y(imid),
     *               fac*3.75d0*dcx,fac*2d0*dcy,rline,0d0,1,0,1)
                call pstext (x(imid)-4d0*dcx*fac,
     *                       y(imid)+0.8d0*dcy*fac,lnms,6)
             else 
                call pselip (x(imid),y(imid),
     *               fac*4.75d0*dcx,fac*2d0*dcy,rline,0d0,1,0,1)
                call pstext (x(imid)-4d0*dcx*fac,
     *                       y(imid)+0.8d0*dcy*fac,lnms,6)
             end if

         else
            if (x(imid).eq.x(imid-1)) then
               theta = 1.5708d0
            else if (y(imid).eq.y(imid-1)) then
               theta = 0d0
            else 
c                                   text label, get angle:
               theta = atan ( ((y(imid)-y(imid-1)) / ylen)
     *                      / ((x(imid)-x(imid-1)) / xlen) )
            end if 
 
            call rxntxt (string,iend)
c                                   normalized x text length:
            clx = iend * dcx * fac  / xlen / 2d0
c                                   normalized y text displacement:
            cly = .5d0 * dcy / ylen
c                                   get displacements from midpoint:
            dy = (sin(theta) * clx + cos(theta) * cly) * ylen
 
            dx = (cos(theta) * clx - sin(theta) * cly) * xlen
 
            theta = 57.29578d0 * theta
 
            call pssctr (ifont,fac,fac,theta)
            call pstext (x(imid)-dx,y(imid)-dy,string,iend)
 
         end if

      else
c                                  invariant point labels:
         if (ird.gt.999) then
            write (pnms,1040) ird
         else if (ird.gt.99) then 
            write (pnms,1010) ird
         else if (ird.gt.9) then
            write (pnms,1020) ird
         else
            write (pnms,1030) ird
         end if
 
         call pssctr (ifont,fac*nscale,fac*nscale,0d0)
         call pstext (x(1)+.7d0*dcx,y(1)+.5d0*dcy,pnms,6)
 
      end if

1010  format ('(',i3,') ')
1020  format ('(',i2,')  ')
1030  format ('(',i1,')   ')
1040  format ('(',i4,')')
 
      end

c---------------------------------------------------------------
      subroutine rname (iw,prompt)

      implicit  none

      integer iw,jxct

      character*10 xnam, prompt*14

      integer ixct,iex,jex,ict
      common/ excl1 /ixct(3),iex(50,3),jex(50,3),ict(3)

      character*10 xnams
      common/ excl4 /xnams(50,3)

      jxct = 1

      do 
         write (*,1040) prompt
         read (*,'(a)') xnam

         if (xnam.eq.' ') exit

         call matchi (xnam,iex(jxct,iw),jex(jxct,iw))

         if (iex(jxct,iw).eq.-1) then
            write (*,1100) xnam
            cycle
         end if

         xnams(jxct,iw) = xnam
         jxct = jxct + 1

      end do 
 
      ixct(iw) = jxct - 1

1040  format (/,'Enter the name of a phase ',a,' fields',
     *        ' (left justified, <cr> to finish): ')
1100  format (/,'No such entity as ',a,', try again: ')

      end 

      subroutine matchi (unnown,itis,icpd)
c----------------------------------------------------------------
 
c matchi - subroutine to determine if the string unnown is a valid
c          solution or compound name.
 
c   itis = 0 if compound, = ikp if solution, and =-1 if invalid
c   icpd = the id of the compound if itis = 0

      implicit none
 
      character*10 unnown
 
      include 'perplex_parameters.h'

      integer i,itis,icpd

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer iphct
      common/ ln4 /iphct
c----------------------------------------------------------------------
      do i = 1, isoct
         if (unnown.eq.fname(i)) then
             itis = i
             return
         end if
      end do
 
      do i = 1, iphct
         if (unnown.eq.names(i)) then
            itis = 0
            icpd = i
            return
         end if
      end do 
 
      itis = -1

      end

      subroutine checki (iw,iun,itis)
c----------------------------------------------------------------------
c checki - subroutine to determine if the phase iun is in the
c          exclude lists (iex, jex), two versions one in psect checks
c          only one index list, the other on psvdraw checks two
c          lists.

c   itis = -1 if not, itis>=0 if it is.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, iw, id, iun, itis
 
      integer ixct, iex, ict, jex
      common/ excl1 /ixct(3),iex(50,3),jex(50,3),ict(3)

      integer ikp
      common/ phase /ikp(k1)
c----------------------------------------------------------------------

      if (ikp(iun).ne.0) then
         id = ikp(iun)
         do i = 1, ixct(iw)
            if (id.eq.iex(i,iw)) then
               itis = iex(i,iw)
               return
             end if
         end do 
      end if
 
      do i = 1, ixct(iw)
         if (iun.eq.jex(i,iw)) then
            itis = 0
            return
         end if
      end do 
 
      itis = -1

      end

