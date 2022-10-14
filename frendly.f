 
c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy and
c Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c Please do not distribute this source.
 
c------------------------------------------------------------------------
 
c FRENDLY - a thermodynamic calculator for petrologic problems.
 
c------------------------------------------------------------------------

      program FRNDLY
 
      implicit none
 
      include 'perplex_parameters.h'
 
      logical nonlin

      character uname*8, y*1, rxny*1, opname*100

      integer i, j, k, l, idiag, ier, iord

      double precision coef(0:10)

      double precision gcpd
      external gcpd

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer inc,jpot
      common/ cst101 /inc(l2),jpot

      integer iam
      common/ cst4 /iam

      save idiag
c----------------------------------------------------------------------- 
c                                 iam flag indicates the Perple_X program
      iam = 5
c                                 version info
      call vrsion (6)
c                                 assign data files
      call fopen2 (2,opname)
c                                 read options
      opname = 'perplex_option.dat'
      call redop1 (.false.,opname)
c                                 override T_melt option
      nopt(20) = 0d0
c                                 flag to indicate whether a plot file is open
      idiag = 0
c                                 harass the user for no reason
      write (*,1040)
      read (*,'(a)') uname

      if (uname.eq.' ') uname = ' Nimrod '

      write (*,1050) uname
      read (*,'(a)') y

      if (y.ne.'y'.and.y.ne.'Y') then
         write (*,1060) uname
         stop
      end if

      do 

         write (*,1030)
c                                 read icopt, default icopt = 2.
         call rdnumb (r,0d0,icopt,2,.false.)

         if (icopt.eq.4) exit 

         call jnput2 (rxny,uname)
 
         nonlin = .false.

         if (icopt.eq.1) then
c                                 calculating equilibrium properties for a
c                                 reaction:
c                                 select variables and set up plot file:
            if (idiag.eq.0) call setplt (.false.,nonlin,coef,iord)
 
            idiag = 1
 
            do 

               call eqrxn 
 
               write (*,1070)
               read (*,'(a)') rxny

               if (rxny.ne.'y'.and.rxny.ne.'Y') exit 
  
               call jnput2 (rxny,uname)
                
            end do
 
         else if (icopt.eq.2) then
c                                 calculating properties at arbitrary conditions:
            write (*,1180) 
            read (*,'(a)') y

            if (y.eq.'y'.or.y.eq.'Y') then
c                                 tabulated properties
               write (*,1190)         
               read (*,'(a)') y

               if (y.eq.'y'.or.y.eq.'Y') then 
c                                 non-linear 1d table
                  nonlin = .true.
                  call setplt (.true.,nonlin,coef,iord)
               else 
c                                 linear 1d or 2d table
                  call setplt (.true.,nonlin,coef,iord)       
               end if 

               do k = 1, inc(iv(3))

                  v(iv(3)) = vmin(iv(3)) + dfloat(k-1)*dv(iv(3))

                  do j = 1, inc(iv(2))

                     v(iv(2)) = vmin(iv(2)) + dfloat(j-1)*dv(iv(2))

                     do i = 1, inc(iv(1))

                        v(iv(1)) = vmin(iv(1)) + dfloat(i-1)*dv(iv(1))

                        if (nonlin) then 

                           v(iv(2)) = 0d0

                           do l = 0, iord
                              v(iv(2)) = v(iv(2)) + coef(l)*v(iv(1))**l
                           end do

                        end if 

                        call calphp 

                        call outphp (.true.)

                     end do 

                  end do 

               end do 
    
               close (n4)

               write (*,1120)

               cycle 

            else 

               do 
c                                 interactively entered conditions
                  write (*,1100)
                  read (*,*,iostat=ier) v(1),v(2)
                  if (ier.ne.0) then 
                     call rerr
                     cycle
                  end if 

                  if (v(1).eq.0d0) exit 

                  if (ifct.gt.0) then

                     do 
                        write (*,1110)
                        read (*,*,iostat=ier) v(3)
                        if (ier.eq.0) exit
                        call rerr
                     end do 

                  end if

                  call calphp 

                  call outphp (.false.)

                  write (*,1090)
                  read (*,'(a)') y

                  if (y.eq.'y'.or.y.eq.'Y') call change 
c
               end do 

            end if 

         else if (icopt.eq.5) then 
c                                 special file output
             do 
                write (*,*) 'Enter P-T file name '
                read (*,*) opname
                open (666,file=opname,status='old',iostat=ier)
                if (ier.ne.0) cycle
                exit
             end do

             do 
                write (*,*) 'Enter output file name '
                read (*,*) opname
                open (667,file=opname,status='new',iostat=ier)
                if (ier.ne.0) cycle
                exit
             end do 

             do 
                read (666,*,iostat=ier) v(1), v(2)
                if (ier.ne.0) exit
                write (667,'(3(g14.6,1x))') gcpd(1,.false.), v(1), v(2)
             end do

             close (666)
             close (667)

         else if (icopt.eq.3) then
c                                 create a new data base entry
            call nentry
 
         end if
 
      end do 
 
      write (*,1130) uname
 
      if (idiag.eq.1) write (n4,1010) 1,1,1,1,1,0,0,0,0,1d0,0,0

1000  format (80(g14.7,1x))
1010  format (9(1x,i1),/,f3.1,/,2(1x,i1))
1030  format (/,'Choose from the following options:',
     *  //,2x,'1 - equilibrium coordinates for a reaction.',
     *   /,2x,'2 - [default] thermodynamic properties for a phase or',
     *          ' reaction relative to',/,6x,'the reference state.',
     *   /,2x,'3 - create new thermodynamic data file entries.',
     *   /,2x,'4 - quit.',
     *  //,'With options 1-2 you may also modify',
     *     ' thermodynamic data, the modified data',/,'can then',
     *     ' be stored as a new entry in the thermodynamic data',
     *     ' file.',/)
1040  format (/,'Hi! I am a user freindly program.',/,
     *          'What is your name? ')
1050  format (/,'Hi ',a,'!, I like you and I hope we have fun!',/,
     *          'Can you say "therm-o-dy-nam-ics" (y/n)?')
1060  format (/,'Weeell ', a,' why dont you go read Gibbs and try me',
     *          ' again later?')
1070  format ('Calculate a different equilibrium (y/n)?')
1090  format ('Modify or output thermodynamic parameters (y/n)? ')
1100  format ('Enter P(bars) and T(K) [zeroes to quit]:')
1110  format ('Enter X(CO2/O) in fluid phase:')
1120  format (/,'The table has been written.',/)
1130  format (/,'Have a nice day ',a,'!',/)
1180  format ('Write a properties table (Y/N)?')
1190  format ('Compute properties along a path (Y/N)?')
 
      end

      subroutine chptx
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer ier, i, j

      character*8 vname,xname
      common/ csta2 /xname(k5),vname(l2)

      double precision delv
      common/ cst63 /delv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
 
      write (*,1020)

      do i = 1, ipot

         j = iv(i) 

         do 
            write (*,1000) vname(j),vmin(j),vmax(j)
            read (*,*,iostat=ier) vmin(j),vmax(j)

            if (j.eq.3.and.vmin(j).lt.0d0.or.j.eq.3.and.vmax(j).gt.1d0
     *          .or.j.ne.3.and.vmin(j).ge.vmax(j).or.ier.ne.0) then

               write (*,1010) 
               cycle

            end if

            exit 

         end do

         v(j) = vmin(j)
         delv(j) = vmax(j) - vmin(j) 
         dv(j) = delv(j) / 4d1

      end do 
 
      call concrt

1000  format (/,'Enter new min/max values for ',a,' (',
     *           'old values were ',g12.5,',',g12.5,')',/)
1010  format (/,'Try again.',/)
1020  format (/,'This option does not change plot limits!'
     *         ,'To do this, modify default plot options',
     *        /,'while running PSVDRAW.',/)
      end 
 
      subroutine eqrxn 
c------------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i
   
      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer inc,jpot
      common/ cst101 /inc(l2),jpot
c------------------------------------------------------------------------
c                                search for an equilibrium point
c                                on the x-y coordinate frame.
      do i = 1, inc(iv(3))

         v(iv(3)) = vmin(iv(3)) + dfloat(i-1)*dv(iv(3))
         call newhld

      end do 
 
      end
 
      subroutine setplt (table,nonlin,coef,iord)
c----------------------------------------------------------------------
c select variables for a plot or table
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character y*1,n4name*100,title*100,text*200

      integer i,j,ier,ic,ix,iord 

      logical table, nonlin

      character*14 tags(27)

      double precision coef(0:10)

      integer inc,jpot
      common/ cst101 /inc(l2),jpot

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      double precision delv
      common/ cst63 /delv(l2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer io3,io4,io9
      common/ cst41 /io3,io4,io9

      logical oned
      common/ cst82 /oned

      save tags

      data tags/      'g(J/mol)     ','h(J/mol)     ','log10[Keq]   ',
     *'s(J/mol/K)   ','cp(J/mol/K)  ','v(J/mol/bar) ','alpha(1/K)   ',
     *'beta(1/bar)  ','N(g/mol)     ','rho(kg/m3)   ','Gruneisen_T  ',
     *'Ks(bar)      ','KsT(bar/K)   ','KsP          ','Gs(bar)      ',
     *'GsT(bar/K)   ','GsP          ',
     *'v0(km/s)     ','v0T(km/s/K)  ','v0P(km/s/bar)',
     *'vp(km/s)     ','vpT(km/s/K)  ','vpP(km/s/bar)','vs(km/s)     ',
     *'vsT(km/s/K)  ','vsP(km/s/bar)','cp/cv        '/
c----------------------------------------------------------------------
 
      do i = 1, l2
         iv(i) = i
         vmin(i) = 0d0
         vmax(i) = vmin(i)
         dv(i) = 1d0
         inc(i) = 1
      end do 

      oned = .false.
c                                 query for 1d table
      if (table) then
         if (nonlin) then
            oned = .true.
         else
            write (*,1110) 
            read (*,'(a)') y
            if (y.eq.'y'.or.y.eq.'Y') oned = .true.
         end if 
      end if 

      do 
c                                 select the x variable (iv(1)):
         if (oned) then 
            write (*,1120)
         else 
            write (*,2130)
         end if 

         write (*,2140) (j,vname(iv(j)), j = 1, ipot)
         write (*,*)
         read (*,*,iostat=ier) ic
         if (ier.eq.0) exit
         call rerr 

      end do 

      ix = iv(1)
      iv(1) = iv(ic)
      iv(ic) = ix
c                                 get x variable limits and increment
      do 

         if (table) then
            write (*,1000) vname(iv(1))
            read (*,*,iostat=ier) vmin(iv(1)),vmax(iv(1)),dv(iv(1)) 
         else 
            write (*,2150) vname(iv(1))
            read (*,*,iostat=ier) vmin(iv(1)),vmax(iv(1))
         end if 

         if (ier.eq.0) exit
         call rerr

      end do 

      if (oned.and.nonlin) then 
c                                 select the dependent variable
         jpot = ipot

         if (ipot.eq.3) then

            do 

                write (*,1030)
                write (*,2140) (j,vname(iv(j)), j = 2, ipot)
                write (*,*)
                read (*,*,iostat=ier) ic

                if (ier.eq.0) exit 
                call rerr

             end do 

         else

            ic = 2

         end if

         ix = iv(2)
         iv(2) = iv(ic)
         iv(ic) = ix

         vmax(iv(2)) = vmin(iv(2))

         if (ipot.eq.3) then 
c                                 specify the value for the 3rd variable
            write (*,2180) vname(iv(3))
            read (*,*,iostat=ier) vmin(iv(3))
            vmax(iv(3)) = vmin(iv(3))

         end if 
c                                 get the path function
         write (*,1040) vname(iv(2)),vname(iv(1))
         read (*,*) iord

         do i = 0, iord
            write (*,1050) i
            read (*,*) coef(i)
         end do

         write (text,1070) vname(iv(2)), 
     *                     (coef(i),vname(iv(1)),i,i=0,iord)
         call unblnk (text)
         write (*,1070) text

      else if (oned) then 
c                                 specify the sectioning variables
         jpot = ipot

          do i = 2, ipot

             do 
                write (*,2180) vname(iv(i))
                read (*,*,iostat=ier) vmin(iv(i))
                if (ier.eq.0) exit 
                call rerr
             end do

             vmax(iv(i)) = vmin(iv(i))

          end do 

      else
c                                 specify the additional variables
         jpot = 2

         if (ipot.eq.3) then

            do 

                write (*,2110)
                write (*,2140) (j,vname(iv(j)), j = 2, ipot)
                write (*,*)
                read (*,*,iostat=ier) ic

                if (ier.eq.0) exit 
                call rerr

             end do 

         else

            ic = 2

         end if

         ix = iv(2)
         iv(2) = iv(ic)
         iv(ic) = ix
 
         do 

            if (table) then
               write (*,1000) vname(iv(2))
               read (*,*,iostat=ier) vmin(iv(2)),vmax(iv(2)),dv(iv(2)) 
            else
               write (*,2150) vname(iv(2))
               read (*,*,iostat=ier) vmin(iv(2)),vmax(iv(2))
            end if 

            if (ier.eq.0) exit 
            call rerr

         end do 
c                                 third variable?
         if (ipot.eq.3) then 

            if (table) then 

               write (*,1010) vname(iv(3))
               read (*,'(a)') y
c
               if (y.eq.'y'.or.y.eq.'Y') then
c                                 setting jpot = 3 will cause plotting programs 
c                                 to crash. 
                  jpot = 3

                  write (*,1000) vname(iv(3))
                  read (*,*,iostat=ier) vmin(iv(3)),vmax(iv(3)),
     *                                  dv(iv(3)) 

               end if 
  
            else 

               write (*,2160) vname(iv(3))
               read (*,'(a)') y

               if (y.eq.'y'.or.y.eq.'Y') then

                  jpot = 3

                  do 
                     write (*,2150) vname(iv(3))
                     read (*,*,iostat=ier) vmin(iv(3)),vmax(iv(3))
                     if (ier.eq.0) exit               
                     call rerr
                  end do

                  do  
                     write (*,1060)
                     read (*,*,iostat=ier) inc(iv(3))
                     if (ier.eq.0) exit 
                     call rerr
                  end do 

                  if (inc(iv(3)).lt.1) inc(iv(3)) = 1

               end if 

            end if 
c                                 set third variable if unused:
            if (jpot.eq.2) then 

               do 

                  write (*,2180) vname(iv(3))
                  read (*,*,iostat=ier) vmin(iv(3))
                  if (ier.eq.0) exit 
                  call rerr

               end do 

               vmax(iv(3)) = vmin(iv(3))
               dv(iv(3)) = 0d0
               inc(iv(3)) = 1

            end if 

         end if 

      end if 
c                                 increments or counters:
      if (table) then
 
         do i = 1, jpot

            inc(iv(i)) = idint(dabs(vmax(iv(i))-vmin(iv(i)))/dv(iv(i))) 
     *                   + 1

         end do 

      else 

         do i = 1, jpot
            inc(iv(i)) = 100
            delv(iv(i)) = vmax(iv(i))-vmin(iv(i)) 
            dv(iv(i)) = delv(iv(i)) / dfloat(inc(iv(i)))
         end do 
c                                 set convergence criteria for univeq:
         call concrt
c                                 plot file output?
         write (*,1020) 
         read (*,'(a)') y

         if (y.eq.'y'.or.y.eq.'Y') then
            io4 = 0
         else 
            io4 = 1
         end if       

      end if 
c                                 setup output files:
   
      if (table.or.io4.eq.0) then
c                                 query output file name
         write (*,1080)
c                                 readrt loads the root into prject
         call readrt 

         if (table) then 
            call mertxt (n4name,prject,'.tab',0)
         else
            call mertxt (n4name,prject,'.plt',0)
         end if   

         open (n4,file=n4name)       
c                                  plot blurb
         if (oned) then 

            call plblrb (4)

         else if (table) then 

            call plblrb (1)

         else

            call plblrb (2)

         end if   
c                                 query for title
         write (*,'(/,a)') 'Enter calculation title:'
         read (*,'(a)') title
 
      end if
c                                 write file headers
      if (table) then 
c                                 terminal info on variables
c        write (*,1090)
c        write (*,'(80(a14,1x))') (vname(iv(i)),i=1,jpot),tags
c                                 write version flag       
         write (n4,'(a,/,a)') '|6.6.6',title 
         write (n4,*) jpot

         do i = 1, jpot
            write (n4,*) vname(iv(i))
            write (n4,*) vmin(iv(i))
            write (n4,*) dv(iv(i))
            write (n4,*) inc(iv(i))
         end do 

         write (n4,*) 27+jpot
         write (n4,'(80(a14,1x))') (vname(iv(i)),i=1,jpot),tags

      else if (io4.eq.0) then 

         write (n4,*) 1
         write (n4,*) 0, 0, 0
         write (n4,*) 0, 0, 0, 0, 0, 0
         write (n4,'(a)') title,' ',' ',' '
         write (n4,*) ipot,(iv(i),i=1,ipot),1,2
         write (n4,'(a)') '0 0 0 0. 0. 0. 0. 0.'
         write (n4,'(6(g11.5,1x))') (vmax(iv(i)),vmin(iv(i)),i=1,ipot)
         write (n4,'(a)') (vname(iv(i)),i=1,ipot)

      end if 

1000  format ('Enter minimum, maximum, and increment for ',a,':')
1010  format (/,'Make the table also a function of ',a,' (y/n)?',//,
     *       'WARNING: if you answer yes, then the resulting 3d table ',
     *       'cannot be plotted ',/,
     *       'with current Perple_X programs or scripts.',/)
1020  format (/,'Generate a plot file (y/n)?')
1030  format (/,'Select the dependent path variable:',/)
1040  format (/,'Profile must be described by the function',/,a,
     *        ' = Sum ( c(i) * ',a,' ^i, i = 0..n)',/,'Enter n (<10)')
1050  format (/,'Enter c(',i2,')')
1060  format (/,'Enter number of sections:')
1070  format (a,'=',5('+(',g12.6,')','*',a,'^',i1))
1080  format (/,'Enter the output file name [without the ',
     *          '.plt/.tab suffix, default is my_project]:')
1090  format (/,'Table columns will be:',/)
1100  format (/,'Your polynomial is:',/,a)
1110  format (/,'Make a 1-dimensional (e.g., isobaric) table (y/n)?')
1120  format (/,'Select the independent (x) variable:',/)

2130  format (/,'Select the first independent (x) variable:',/)
2140  format (10x,i1,' - ',a)
2150  format (/,'Enter minimum and maximum values for ',a,':')
2110  format (/,'Select the second independent (y) variable:',/)
2180  format (/,'Specify the value for: ',a)
2160  format (/,'Calculate multiple sections as a function of ',a,
     *          ' (y/n)?')
 
      end 

      subroutine newhld
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character y*1

      integer ivi,ivd,ier,igo

      double precision div

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*8 vname,xname
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5
c-----------------------------------------------------------------------
c                                  initialization:
c       ifuk = 0
10    write (*,1160)
      write (*,1170) vname(iv1),vname(iv2)
c                                  write potential variable sectioning
c                                  constraints:
      if (ipot.gt.2) write (*,1180) vname(iv3), v(iv3)
c                                  set starting values for search
      v(iv1)=vmin(iv1)
      v(iv2)=vmin(iv2)
c                                  test stability along an edge of the
c                                  diagrams coordinate frame:
      call search (ivi,ivd,div,ier)

      if (ier.eq.1) then
         write (*,1010)
         goto 20
      end if
 
      call trace (ivd,ivi,div,igo)
c                                  this would permit continuation of the
c                                  search on failures from trace if iste
c                                  were not initialized in search (but here).      
c      if (igo.eq.1.and.ifuk.eq.0) then
c         ifuk = 1
c         goto 30
c      else if (igo.eq.1.and.ifuk.eq.1) then
c         write (*,1030)
c      end if         
 
20    write (*,1040)
      read (*,1000) y
      if (y.eq.'y'.or.y.eq.'Y') then
         call chptx
         goto 10
      end if 
      write (*,1020)
      read (*,1000) y
      if (y.ne.'y'.and.y.ne.'Y') return
      call change 
      goto 10
 
1000  format (a)
1010  format (/,'Equilibrium is not in specified',
     *          ' coordinate frame.',/)
1020  format (/,'Modify data and',
     *        ' recalculate the equilibrium (y/n)? ')
1040  format (/,'Change PTX limits (y/n)?',/)
1160  format (/,'-------------------------------------------------'
     *          ,'---------------',/)
1170  format ('The ',a,'-',a,' loci of the univariant field'
     *        ,' follows:')
1180  format ('(subject to the constraint ',a,'=',g12.6,')',/)


99    end
 
      subroutine search (ivi,ivd,div,ier)
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ivi,ivd,ier,i
    
      double precision div,ddv,gval,gst
 
      double precision a,b
      integer ipvt,idv,iophi,idphi,iiphi,iflg1
      common/ cst23 /a(k8,k8),b(k8),ipvt(k8),idv(k8),
     *               iophi,idphi,iiphi,iflg1

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c-----------------------------------------------------------------------
c                                 initialization
      ier = 0

      v(iv1) = vmin(iv1)
      v(iv2) = vmin(iv2)

      call grxn (gst)
 
      do i = 1, 4
c                                 set default dependent and independent
c                                 variables and increments for each edge
c                                 to be searched.
         if (i.eq.1) then 
c                                 traverse 1.
            ivi = iv2
            ivd = iv1
            ddv = dv(ivd)
            div = dv(ivi)
            v(ivi) = vmin(ivi)

         else if (i.eq.2) then 
c                                 traverse 2.
            ivi = iv1
            ivd = iv2
            ddv = dv(ivd)
            div = -dv(ivi)
            v(ivi) = vmax(ivi)

         else if (i.eq.3) then 
c                                 traverse 3.
            ivi = iv2
            ivd = iv1
            ddv = -dv(ivd)
            div = -dv(ivi)
            v(ivi) = vmax(ivi)

         else if (i.eq.4) then 
c                                 traverse 4.
            ivi = iv1
            ivd = iv2
            ddv = -dv(ivd)
            div = dv(ivi)
            v(ivi) = vmin(ivi)

         end if 

         do 
c                                 begin search:
            v(ivd) = v(ivd) + ddv
c                                 out of range?:
            if (i.le.2) then 

               if (v(ivd).gt.vmax(ivd)) v(ivd) = vmax(ivd)

            else 

               if (v(ivd).lt.vmin(ivd)) v(ivd) = vmin(ivd)

            end if 
c                                 calculate phase energies:
            call grxn (gval)

            if (gval*gst.lt.0d0) return
c                                 check if search is in range:
            if (i.le.2) then 
               if (v(ivd).ge.vmax(ivd)) exit 
            else 
               if (v(ivd).le.vmin(ivd)) exit 
            end if 

         end do  
c                                 next traverse.
      end do 
c                                 set this ier flag to indicate
c                                 the reaction wasn't found
      ier=1
c                                 done:
      end
 
      subroutine trace (iovd,iovi,odiv,igo)
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer iovd,iovi,ivd,igo,ier,jer,icter,ird,ivi

      double precision odiv,div

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
     
      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2
c-----------------------------------------------------------------------
      ivi = iovi
      ivd = iovd
      igo = 0

50    call univeq (ivd,ier)
c                                 if univeq fails on a bounding edge
c                                 write error message and return:
      goto (9000,9000),ier
c                                 set the increment for the iv
      div = odiv
c                                 initialize counters
      ipt2 = 0
      icter = 0
c                                 assign the 1st point
      call assptx
c                                 follow the equilibrium
60    call sfol1 (ivd,ivi,ier,div)
c                                 sfol1 returns ier=1
      goto (70,70),ier
      ivi = iovi
      ivd = iovd
      goto 9999

70    call switch (div,ivi,ivd,jer)
      goto (75),jer
      icter=icter+1
      if (icter.lt.4) goto 60
 
75    call warn (10, v(ivi), igo, 'TRACE')
 
      call outrxn
      ivi = iovi
      ivd = iovd
c                                 return on error
      goto 9999
c                                 error in univeq:
9000  call warn (79, v(ivi), ird, 'TRACE')
      write (*,*) ' failed at P=',v(1),' T=',v(2),' XCO2 =',v(3)
      goto (9999), igo
      ivi = iovd
      ivd = iovi
      igo = 1
      goto 50
 
9999  end
 
      subroutine sfol1 (ivd,ivi,ier,dv)
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ivi,ivd,ier

      double precision dv
 
      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,odv
      common/ cst9  /vmax(l2),vmin(l2),odv(l2)

      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2
c-----------------------------------------------------------------------
      do 
c                                 begin traverse:
         v(ivi)=v(ivi)+dv
c                                 is search in range?
         if (v(ivi).gt.vmax(ivi)) then
            v(ivi)=vmax(ivi)
         else if (v(ivi).lt.vmin(ivi)) then
            v(ivi)=vmin(ivi)
         end if
c                                 solve for the equilibrium conditions:
         call univeq (ivd,ier)
c                                 on error return:
c                                 calling routine will switch variables.
         if (ier.ne.0) return 

         if (ipt2.gt.449) goto 9000
c                                 dependent v in range? if
c                                 greater than the maximum value for v
c                                 or less than the minimum value for v
c                                 reset conditions, and
c                                 switch independent/dependent variables
         if (v(ivd).gt.vmax(ivd)) then

            v(ivd)=vmax(ivd)

         else if (v(ivd).lt.vmin(ivd)) then

            v(ivd)=vmin(ivd)

         else

            call assptx

            if ((v(ivi).eq.vmax(ivi)).or.(v(ivi).eq.vmin(ivi)))
     *         goto 9000
            
            cycle 

         end if

         exit 

      end do 
c                                 solve for the equilibrium with
c                                 the switched variables:
      call univeq (ivi,ier)
c                                 assign final value
      if (ier.eq.0) call assptx
c                                 output the traversed equilibrium:
9000  call outrxn

      ier = 0

      end
 
      subroutine outrxn
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,l
 
      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ipt2
      double precision ptx
      common/ cst32 /ptx(l5),ipt2

      integer io3,io4,io9
      common/ cst41 /io3,io4,io9
c-----------------------------------------------------------------------
      if (iphct.gt.4) goto 20
      write (*,1050) (vnu(l),names(l),l=1,iphct)
      goto 30
20    write (*,1050) (vnu(l),names(l),l=1,4)
      write (*,1060) (vnu(l),names(l),l=5,iphct)
30    write (*,*)
      write (*,1000) (ptx(i),i=1,ipt2)
      write (*,*)
      goto 10
 
10    goto (99), io4
      if (ipt2.eq.0) goto 99
 
      write (n4,1010) ipt2,0,1,iphct,(i,i=1,iphct),0,0,0,0
      write (n4,1020) (vnu(l),l=1,iphct)
      write (n4,1000) (ptx(i),i=1,ipt2)
 
1000  format (3(1x,g10.4,1x,g10.4,3x))
1010  format (20(i5,1x))
1020  format (10(g9.3,1x))
1050  format (/,4(1x,g9.3,1x,a))
1060  format (6x,4(1x,g9.3,1x,a),/,6x,4(1x,g9.3,1x,a))
99    end

      subroutine change 
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,ier,id,jdis,imurg,kv,ichk,h,k,kd,jd
 
      character y*1

      double precision vsum

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(11),wstrg(m16),
     *               e16st(12)

      integer eos
      common/ cst303 /eos(k10)
c-----------------------------------------------------------------------
      write (*,1110)
      read (*,1050) y
 
      if (y.eq.'y'.or.y.eq.'Y') then 
         
         if (iphct.gt.1.or.vnu(1).ne.1d0) then
            write (*,1120)
            read (*,1050) y
         else 
            y = 'n'
         end if 
 
         if (y.ne.'y'.and.y.ne.'Y') then

            if (iphct.gt.1) then 
 
               do 
                  write (*,1000) (i,names(i),i=1,iphct)
                  read (*,*,iostat=ier) id
                  if (ier.eq.0) exit 
                  call rerr
               end do 

            else
               id = 1
            end if 

            write (*,1040) names(id)
            read (*,1050) names(id)

            call unlam (tm,id)

            call unver (
c                                 g0, s0, v0
     *            thermo(1,id),thermo(2,id),thermo(3,id),
c                                 c1-c8
     *            thermo(4,id),thermo(5,id),thermo(6,id),thermo(7,id),
     *            thermo(8,id),thermo(9,id),thermo(10,id),thermo(24,id),
c                                 b1-b11 
     *            thermo(11,id),thermo(12,id),thermo(13,id),
     *            thermo(14,id),thermo(15,id),thermo(16,id),
     *            thermo(17,id),thermo(18,id),thermo(19,id),
     *            thermo(20,id),thermo(21,id),
c                                 ref stuff
     *            tr,pr,eos(id))
c                                 add in activity correction
            thermo(1,k10) = thermo(1,id)
            thermo(2,k10) = thermo(2,id)
            thermo(1,id) = thermo(1,id) + tr * r * dlog (act(id))
            thermo(2,id) = thermo(2,id) - r * dlog (act(id))
 
            call append (n2)
            call outdat (n2,id,2)
c                                 reset data
            thermo(1,id) = thermo(1,k10) 
            thermo(2,id) = thermo(2,k10) 

         else

            id = k5
            idis(id) = 0
            ltyp(id) = 0
            lct(id)  = 0
            eos(id)  = eos(1)
            jdis = 0

            write (*,1130)
            read (*,1050) names(id)
 
            do i = 1, k4
               thermo(i,id) = 0d0
            end do 

            do i = 1, k0
               cp0(i,id) = 0d0
            end do 
 
            imurg = 0
            vsum = 0d0

            do j = 1, iphct

               if (thermo(18,j).ne.0d0.and.iphct.gt.1) 
     *            call warn (45,r,iphct,names(id))
               if (eos(j).ne.eos(id)) call warn (99,r,j,'combining di'//
     *            'fferent EoS may not function correctly')

               if (eos(j).ne.13) then

                  do i = 1, k4
                     thermo(i,id) = thermo(i,id) + vnu(j) * thermo(i,j)
                  end do

               else
c                                 komabayashi, caloric molar weighted
                  do i = 1, 10
                     thermo(i,id) = thermo(i,id) + vnu(j) * thermo(i,j)
                  end do

                  vsum = vsum + vnu(j) * thermo(3,j)
c                                 compressibility and expansivity volume weighted
c                                 and k'(T) for lack of a better choice. 
                  do i = 11, 20
                     thermo(i,id) = thermo(i,id) 
     *                            + vnu(j) * thermo(i,j) * thermo(3,j) 
                  end do

               end if

               do i = 1, k0
                  cp0(i,id) = cp0(i,id) + vnu(j) * cp0(i,j)
               end do 

               if (idis(j).ne.0) then

                  idis(id) = m9
                  jdis = jdis + 1
                  jd = idis(j) 
                  if (jdis.gt.1) goto 91

                  do i = 1, 7      
                     therdi(i,m9) = vnu(j) * therdi(i,jd)
                  end do 

                  therdi(8,m9) = therdi(8,jd)
                  therdi(9,m9) = therdi(9,jd)

               end if

               if (ltyp(j).ne.0) then

                  if (ltyp(id).ne.0) goto 91

                  ltyp(id) = ltyp(j)
                  kd = iphct + 1
                  lmda(id) = kd
                  jd = lmda(j)
                  lct(id) = lct(j) 

                  if (ltyp(j).eq.4) then

                     therlm(1,lct(id),kd) = therlm(1,1,jd)
                     therlm(2,lct(id),kd) = vnu(j) * therlm(2,1,jd)
                     therlm(3,lct(id),kd) = vnu(j) * therlm(3,1,jd)

                  else if (ltyp(j).eq.1) then

                     do i = 1, lct(j)
                        therlm(1,i,kd) = vnu(j) * therlm(1,i,jd)
                        therlm(2,i,kd) = vnu(j) * therlm(2,i,jd)
                        therlm(5,i,kd) = vnu(j) * therlm(5,i,jd)
                        therlm(6,i,kd) = vnu(j) * therlm(6,i,jd)
                        therlm(3,i,kd) = therlm(3,i,jd)
                        therlm(4,i,kd) = therlm(4,i,jd)
                        therlm(7,i,kd) = therlm(7,i,jd)
                        therlm(8,i,kd) = therlm(8,i,jd)
                     end do

                  else if (ltyp(j).le.3) then

                     do k = 1, lct(j) 

                        therlm(1,k,kd) = therlm(1,k,jd)
                        therlm(2,k,kd) = therlm(2,k,jd)

                        do h = 3, 12
                           therlm(h,k,kd) = therlm(h,k,jd)*vnu(j)
                        end do

                     end do 

                  end if

               end if

            end do
c                                renormalize volumetric averages:
            if (eos(id).eq.13) then
c                                renormalize volumetric averages:
               do i = 11, 20
                  thermo(i,id) = thermo(i,id)/vsum
               end do

            end if

            call unlam (tm,id)

            call unver (
c                                 g0,s0, v0
     *            thermo(1,id),thermo(2,id),thermo(3,id),
c                                 c1-c8
     *            thermo(4,id),thermo(5,id),thermo(6,id),
     *            thermo(7,id),thermo(8,id),thermo(9,id),
     *            thermo(10,id),thermo(24,id),
c                                 b1-b11 
     *            thermo(11,id),thermo(12,id),thermo(13,id),
     *            thermo(14,id),thermo(15,id),thermo(16,id),
     *            thermo(17,id),thermo(18,id),thermo(19,id),
     *            thermo(20,id),thermo(21,id),
c                                 ref stuff
     *            tr,pr,eos(id))

c                                 add in activity correction
            thermo(1,k10) = thermo(1,id)
            thermo(2,k10) = thermo(2,id)

            do i = 1, iphct

               thermo(1,id) = thermo(1,id) 
     *                      + vnu(i) * tr * r * dlog(act(i))
               thermo(2,id) = thermo(2,id) 
     *                      - vnu(i) * r * dlog(act(i))

            end do 
c                                 output the data 
            call append (n2)
            call outdat (n2,id,2)
c                                 reset data
            thermo(1,id) = thermo(1,k10) 
            thermo(2,id) = thermo(2,k10)  

         end if

         return 

      else 
 
         do 
c                                 phase loop
            if (iphct.gt.1) then

               do 
                  write (*,1000) (i,names(i),i=1,iphct)
                  read (*,*,iostat=ier) id
                  if (ier.eq.0) exit 
                  call rerr
               end do 

            else

               id = 1

            end if 

            ichk = 0

            write (*,1010)

            do
c                                 property loop 
               do 
                  write (*,1020) (i,strgs(i),i=1,18)
                  read (*,*,iostat=ier) kv
                  if (ier.eq.0) exit 
                  call rerr
               end do 

               if (kv.eq.0) then

                  if (ichk.eq.0) exit
c                                 write entry to permanant file:
                  write (*,1070) names(id)
                  read (*,1050) y
c
                  if (y.eq.'y'.or.y.eq.'Y') then
c                                 add in activity correction
                     thermo(1,k10) = thermo(1,id)
                     thermo(2,k10) = thermo(2,id)
                     thermo(1,id) = thermo(1,id) 
     *                            + tr * r * dlog (act(id))
                     thermo(2,id) = thermo(2,id) 
     *                            - r * dlog (act(id))
 
                     call append (n2)
                     call outdat (n2,id,2)

                     thermo(1,id) = thermo(1,k10) 
                     thermo(2,id) = thermo(2,k10) 

                  end if
c
                  call conver (
c                                 g0, s0, v0:
     *            thermo(1,id), thermo(2,id), thermo(3,id),
c                                 c1-c8 
     *            thermo(4,id),thermo(5,id), thermo(6,id), thermo(7,id),
     *            thermo(8,id),thermo(9,id), thermo(10,id),
     *            thermo(24,id),
c                                 b1-b12
     *            thermo(11,id),thermo(12,id),thermo(13,id),
     *            thermo(14,id),thermo(15,id),thermo(16,id),
     *            thermo(17,id),thermo(18,id),thermo(19,id),
     *            thermo(20,id),thermo(21,id),thermo(22,id),
c                                 b13 on return
     *            thermo(23,id),
c                                 ref stuff
     *            tr,pr,r,0)

                  exit 

               end if
c
               if (ichk.eq.0) then

                  call unlam (tm,id)

                  call unver (
c                                 g0,s0, v0
     *            thermo(1,id),thermo(2,id),thermo(3,id),
c                                 c1-c8
     *            thermo(4,id),thermo(5,id),thermo(6,id),
     *            thermo(7,id),thermo(8,id),thermo(9,id),
     *            thermo(10,id),thermo(24,id),
c                                 b1-b11 
     *            thermo(11,id),thermo(12,id),thermo(13,id),
     *            thermo(14,id),thermo(15,id),thermo(16,id),
     *            thermo(17,id),thermo(18,id),thermo(19,id),
     *            thermo(20,id),thermo(21,id),
c                                 ref stuff
     *            tr,pr,eos(id))

                  write (*,1040) names(id)
                  read (*,1050) names(id)
                  ichk = 1

               end if

               if (kv.eq.19) then 

                  do 
                     write (*,1030) 'thermodynamic activity',
     *                              names(id),act(id)
                     read (*,*,iostat=ier) act(id)
                     if (ier.eq.0) exit 
                     call rerr
                  end do 

               else if (kv.eq.20) then 

                  do 
                     write (*,1030) 'stoichiometric coefficient',
     *                              names(id),vnu(id)
                     read (*,*,iostat=ier) vnu(id)
                     if (ier.eq.0) exit 
                     call rerr
                  end do 
 
                  if (id.eq.idf(1)) then
                     vuf(1) = vnu(id)
                  else if (id.eq.idf(2)) then
                     vuf(2) = vnu(id)
                  end if

               else 

                  do 
                     write (*,1030) strgs(kv),names(id),thermo(kv,id)
                     read (*,*,iostat=ier) thermo(kv,id)
                     if (ier.eq.0) exit 
                     call rerr
                  end do 
         
               end if

               cycle
c                                 end property loop
            end do 

            if (iphct.eq.1) exit

            write (*,1150)
            read (*,1050) y

            if (y.ne.'y'.and.y.ne.'Y') exit 
c                                end of phase loop
         end do 
 
      end if 

      if (ifct.eq.0) return 

      write (*,1160)
      read (*,1050) y
      if (y.ne.'y'.and.y.ne.'Y') return 

      call rfluid (1)
c                                 for multispecies fluids set
c                                 up species indices and name
c                                 of independent variable
      call setins (ifug)

      return 

91    call warn (9,t,ifct,names(id))

1000  format (/,'Select phase to modify or output:',9(/,6x,i2,') ',a))
1010  format (/,'For definitions of the EoS parameters refer to:',/
     *        /,'  www.perplex.ethz.ch/perplex_thermodynamic_data_file',
     *          '.html',//,'Units: J, bar, K, mole')
1020  format (/,'Enter the index of the parameter to be modified:',//,
     *          9(2(1x,i2,') ',a,23x),/),
     *         ' 19) thermodynamic activity',3x,
     *         ' 20) reaction coefficient',//,
     *         'Enter zero when you are finished: ')
1030  format (/,'Old value for ',a,' of ',a,' was ',g15.8,/,
     *          'Enter new value: ')
1040  format (/,'Enter a name (<9 characters left justified) to',
     *          ' distinguish the modified',/,'version of ',a,'.',
     *          ' WARNING: if you intend to store modified data',/,
     *          'in the data file this name must be unique.',/)
1050  format (a)
1070  format (/,'Store ',a,' as an entry',
     *          ' in the thermodynamic data file (y/n)?')
1110  format (/,'Do you only want to output data (y/n)? ')
1120  format (/,'Output reaction properties (y/n)? ')
1130  format (/,'Enter an 8 character name for the reaction: ')
1150  format (/,'Modify properties of another phase (y/n)? ')
1160  format (/,'Change fluid equation of state (y/n)? ')

      end
 
      subroutine nentry
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, ier 
 
      character y*1

      integer eos
      common/ cst303 /eos(k10)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(11),wstrg(m16),
     *               e16st(12)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      ier = 0

      write (*,1000) tr,pr
 
      do 
  
         write (*,1010)
         read (*,'(a)') names(k10)
c                                 composition:
         write (*,1020) names(k10), (cmpnt(i),i=1,icmpn)
         write (*,1030) 
c                                 read the text string:
         call formul (5)

         write (*,1050)

         do i = 1, 18
5019        write (*,1040) strgs(i),names(k10)
            read (*,*,iostat=ier) thermo(i,k10)
            call rerror (ier,*5019)
         end do 

c                                 classify eos:
         if (thermo(3,k10).lt.0d0) then 
c                                 negative "volume" signals one of the 
c                                 stixrude EoS's, for these EoS's "s" is
c                                 +/-n - number of atoms pfu.
            if (thermo(2,k10).gt.0d0) then

               eos(k10) = 5
c                                 stixrude & bukowinski JGR '93
            else

              eos(k10) = 6
c                                 stixrude & lithgow-bertelloni GJI '05
            end if 

         else if (thermo(18,k10).eq.0d0) then

            eos(k10) = 1
c                                 normal polynomial vdp term:
         else

            if (thermo(16,k10).eq.0d0) then
               eos(k10) = 3
            else if (thermo(18,k10).ge.3d0) then
               eos(k10) = 2
            else if (thermo(18,k10).le.3d0) then
               eos(k10) = 4
            else 
               eos(k10) = 7
            end if

         end if 
 
         call append (n2)

         call outdat (n2,k10,0)
 
         write (*,1110)
         read (*,'(a)') y
         if (y.ne.'y'.and.y.ne.'Y') exit

      end do 
 
1000  format (/,'This entry will be for a T = ',g13.6,'(K) P=',
     *        g13.6,'(bar)',/,'reference state (Units: J, bar, K).',/)
1010  format ('Enter name for your entry, <8 characters, left',
     *        ' justified.',/,'WARNING: this name must not duplicate',
     *        ' an entry already',/,'in the data file!')
1020  format (/,'Enter the molar formula of ',a,' in terms of the ',
     *          'following components:',(12(1x,a)))
1030  format (/,'Indicate the component stoichiometry by an number ',
     *          'enclosed in parentheses',/,'following the CASE ',
     *          'SENSITIVE component name, no blanks, e.g.:',/,20x,
     *          'CAO(1)AL2O3(1)SIO2(2)',/)
1040  format ('Enter parameter ',a,' for ',a,':')
1050  format (/,'For definitions of the following parameters refer to:',
     *       //,'  www.perplex.ethz.ch/perplex_thermodynamic_data_file',
     *          '.html',/)
1110  format (/,'Make another entry (y/n)?')

      end

      subroutine append (lun)
c---------------------------------------------------------
c routine to allow frendly to append data to the thermo
c data file, necessary for microsoft NT compiler.
c---------------------------------------------------------
      implicit none

      integer lun, ier

      character*1 record
c                          find end of file:
      do 
         read (lun,*,iostat=ier) record
         if (ier.ne.0) exit
      end do 

      backspace (lun)
c                          start new record:
      write (lun,*)
c                          reposition at new
c                          record:
      backspace (lun)

      end

      subroutine jnput2 (rxny,uname)
c----------------------------------------------------------------------
c jnput2 reads all data in the thermodynamic data file, makes a list of
c possible phases and, depending on icopt, allows users to select 
c specific phases.

c icopt = 1-4 frendly
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*1 uname*8, y, rxny,mnames(k16*k17)*8
 
      double precision vvv

      logical eof, first, match

      integer inames, jcmpn, i, j, k, l, ier,
     *        isct, jj, itic, jphct

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      character*8 name
      common/ csta6 /name

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ic
      common/ cst42 /ic(k0)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer make, mkst, mkend
      common / cst335 /make(k10)

      integer jcv,jvct,jpv,jtv,jf
      common / mu2f1 /jcv(2),jvct,jpv,jtv,jf(2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character tname
      common/ csta10 /tname(2)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision atwt
      common/ cst45 /atwt(k0)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      save first, inames, mnames
      data first/.true./
c-----------------------------------------------------------------------

      call topn2 (1)

      if (first) then

         call eohead (n2)

         call smakcp (inames,mnames)

         first = .false.

         call eohead (n2)

         do i = 1, l2
            iv(i) = i
         end do

c         if (icmpn.gt.k5) write (*,3020) k5

      end if 
c                              this check is here to avoid
c                              redimensioning comps, but cause
c                              the maximum number of components
c                              counted by frendly to be k5, if you
c                              want as many components in frendly as
c                              in the data base you must set k0 = k5.
      jcmpn = icmpn
      icomp = icmpn
c                                 component pointers
      do i = 1, icmpn
         ic(i) = i
      end do

      if (icmpn.gt.k5) then
         icomp = k5
         jcmpn = k5
      end if 
c                               list database phases:
      write (*,1030)
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') then 

         write (*,3000) (cmpnt(i),i=1,jcmpn)

         do 

            call getphi (name,.true.,eof)

            if (eof) then
c                               write make list:
               do i = 1, nmak

                  write (*,3010) mknam(i,mknum(i)+1),
     *                           (mcomp(i,j),j=1,jcmpn)
               end do 

               exit 

            end if
c                               acceptable data, count the phase
            write (*,3010) name, (comp(i),i=1,jcmpn)
 
         end do

      end if 

      if (icopt.eq.3) goto 99
c                               initialize: 
      do k = 1, k0
         cp0(k,1) = 0d0
      end do  
c                               initialization for k10 endnmembers
      do i = 1, k10
         make(i) = 0 
         names(i) = ' '
      end do
 
      do k = 1, k5+1
         lmda(k) = 0
         idis(k) = 0
         vnu(k) = 0d0
         act(k) = 0d0
         ltyp(k) = 0
         lct(k) = 0
      end do
 
      do k = 1, m7
          do l = 1, m6
            tm(k,l) = 0d0
          end do 
      end do 
 
      iphct = 0
      lamin = 0 
      ipot = 2
      ifct = 0
      vuf(1) = 0d0
      vuf(2) = 0d0
      idf(1) = 0
      idf(2) = 0
      idf(3) = 0 
 
      if (icopt.eq.1) then
 
         rxny='y'
5002     write (*,4010)
         read (*,*,iostat=ier) isct
         call rerror (ier,*5002)
 
      else if (icopt.lt.8) then
 
         write (*,4000)
         read (*,'(a)') rxny
 
         if (rxny.eq.'y'.or.rxny.eq.'y') then
5007        write (*,4010)
            read (*,*,iostat=ier) isct
            call rerror (ier,*5007)
         else
            isct = 1
         end if

      end if
c                               jj is the pointer to the reaction props
      jj = isct+1
c                               initialize fluid counter outside loop, 9/18/20
      ifct = 0 
c                               get composition vectors for entities
c                               defined by a make definition:
      do i = 1, isct

30       match = .false.

         if (rxny.ne.'y'.and.rxny.ne.'Y') then
            write (*,4020)
         else
            write (*,4030) i
         end if
 
         read (*,'(a)') exname(1)
c                               general input data for main program
         call eohead (n2)
c                               first look in the real data:
         do  

           call getphi (name,.true.,eof)
c                               looked at all real data
           if (eof) exit

           if (exname(1).eq.name) exit 

         end do 

         if (eof) then 

            ieos = 0

            match = .false.
c                                look in the make list
            do k = 1, nmak

               name = mknam(k,mknum(k)+1)

               if (mksat(k).or.exname(1).ne.name) cycle
c                                load make data 
               do j = 1, icmpn
                  comp(j) = mcomp(k,j)
               end do 

               match = .true.

               make(iphct+1) = k

               exit 

            end do 

            if (.not.match) then 
c                                 no match with named phase
               write (*,4050) uname
               goto 30

            end if 

         end if 
c                                 set special flag if O2
         iphct = iphct + 1

         if (exname(1).eq.'O2      ') idf(3) = iphct 

         if (rxny.eq.'Y'.or.rxny.eq.'y') then

            do 
               write (*,1040) exname(1)
               read (*,*,iostat=ier) vvv
               if (ier.eq.0) exit
               call rerr
            end do 

         else
            vvv = 1d0
         end if
 
         if (lopt(7)) then 

            do k = 1, ispec
               if (exname(1).ne.cmpnt(idspe(k))) cycle 
               vuf(k) = vvv
               idf(k) = iphct
               ifct = ifct + 1
               exit 
            end do

         end if
c                                 get activity:
         if (ifct.gt.0) then 

            do k = 1, ispec
               if (exname(1).ne.cmpnt(idspe(k))) cycle 
               act(iphct) = 1d0
               goto 40
            end do 

         end if 

         write (*,4070) exname(1)
c                                 read icopt, default icopt = 2.
         call rdnumb (act(iphct),1d0,icopt,0,.true.)
c                               reaction coefficient:
40       vnu(iphct)= vvv
c                               store the data 
         call loadit (iphct,.false.,.true.)

      end do 

      if (ifct.gt.0) ipot = 3
c                                load the make dependencies               
      jphct = 21
c                                read header
      call eohead (n2)

      do 

         call getphi (name,.false.,eof)

         if (eof) exit

         do i = 1, inames

            if (name.ne.mnames(i)) cycle
c                                matched a name
            jphct = jphct + 1
c                                store thermodynamic parameters:
            call loadit (jphct,.false.,.false.)

         end do 

      end do 
c                                set counters
      mkend = jphct
      mkst  = 21

      do i = 1, nmak
c                                remake pointer array for makes 
         do j = 1, mknum(i)
            do k = mkst, mkend
               if (names(k).ne.mknam(i,j)) cycle
               mkind(i,j) = k
            end do
         end do 
      end do  
c                                 select equation of state for the
c                                 saturated phase.
      if (ifct.gt.0) call rfluid (1)
c                                 compute formula weights
      do l = 1,iphct

         props(17,l) = 0d0 

         do k= 1, k0
            props(17,l) =  props(17,l) + cp0(k,l)*atwt(k)
         end do

      end do 
c
      if (rxny.eq.'y'.or.rxny.eq.'Y') then
c                                 check reaction stoichiometries
55       do k = 1, k0
            cp0(k,jj) = 0d0
         end do 

         do l = 1,iphct
            do k = 1, k0
               cp0(k,jj) = cp0(k,jj) + vnu(l)*cp0(k,l)
            end do 
         end do 

         itic = 0
         do k = 1, k0
            if (cp0(k,jj).eq.0d0) cycle
            itic = itic+1
         end do 

         if (itic.eq.0) goto 99
 
         do k = 1, k0
            if (cp0(k,jj).eq.0d0) cycle
            write (*,2460) cp0(k,jj),cmpnt(k)
         end do 
 
         write (*,1070)
         read (*,'(a)') y
         if (y.ne.'y'.and.y.ne.'Y') goto 99
         call stoich
         goto 55
 
      end if
c                               position pointer to last record of
c                               data file to allow writing:
99    do
         read (n2,'(a)',iostat=i) name
         if (i.ne.0) exit
      end do

1030  format (/,'List database phases (y/n)? ')
1040  format ('Enter reaction coefficient for: ',a,
     *        ' products (+), reactants (-): ')
1041  format (/,'Select a phase or species to be used to define ',
     *        a,'(default =',a,')')
1051  format (/,'The phase or reaction has a T-dependent disordering ',
     *          'function.',/)
1061  format (/,'The phase or reaction has a transition-dependent ',
     *          'function.',/)
1070  format (/,'Change reaction coefficients (y/n)? ')
2460  format ('Warning ** reaction does not balance by ',g13.6,
     *        ' moles of ',a5)
3000  format (t30,'composition',/,' phase',4x,12(a5,1x))
3010  format (1x,a8,1x,12(f5.2,1x))
3020  format (/,'too many components, only first ',i2,
     *          ' are listed.',/)
4020  format ('Calculate thermodynamic properties for phase: ')
4030  format ('Enter phase or species number ',i2,
     *       ' in your reaction: ')
4050  format ('Sorry ',a,', that name is invalid, try again.',/)
4070  format ('Enter activity [default is 1] of ',a)
4000  format (/,
     *    'Calculate thermodynamic properties for a reaction (y/n)?')
4010  format (/,'How many phases or species in the reaction? ')

      end

      subroutine stoich
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character y*1

      integer i,ier,id

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr
c-----------------------------------------------------------------------
      ier = 0
20    write (*,1000) (i,names(i),vnu(i),i=1,iphct)
      write (*,*)
      read (*,*,iostat=ier) id
      call rerror (ier,*20)
 
5019  write (*,1030) names(id),vnu(id)
      read (*,*,iostat=ier) vnu(id)
      call rerror (ier,*5019)
 
      if (id.eq.idf(1)) then
         vuf(1) = vnu(id)
      else if (id.eq.idf(2)) then
         vuf(2) = vnu(id)
      end if
 
      write (*,1020)
      read (*,1010) y
      if (y.ne.'Y'.and.y.ne.'y') goto 99
      goto 20
 
1000  format (/,'Enter number of phase to be modified:',
     *        9(/,6x,i2,') ',a,' reaction coeff.=',f8.4))
1010  format (a1)
1020  format (/,'Modify coefficient of another phase (y/n)? ')
1030  format (/,'Old coefficient for ',a,' was ',f8.4,
     *          ' enter new value: ')
99    end

      subroutine calphp 
c-----------------------------------------------------------------------
c an interface to getphp to get properties of phases or reaction
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      logical sick(i8), ssick, ppois, bulkg, bsick

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c----------------------------------------------------------------------
c                                 initialize
      ntot = iphct
c                                 initialize system properites
      call insysp (ssick,ppois,bulkg,bsick)

      rxn = .false.

      do i = 1, ntot
c                                 check for reactions
         if (vnu(i).lt.0d0) rxn = .true.

      end do 

      do i = 1, ntot
c                                 set molar amount of phase 
         props(16,i) = 1d0
c                                 getphp uses the sign of the phase
c                                 pointer to discrimate between pure
c                                 compounds and solutions
         call getphp (-i,i,sick,ssick,ppois,bulkg,bsick)

      end do 
c                                 compute aggregate properties:
      call gtsysp (sick,ssick,bulkg,bsick)

      end 

      subroutine outphp (table)
c----------------------------------------------------------------------
c outphp - output properties of a phase, reaction or phase aggregat, 
c called by frendly. if meemum, prints chemical potentials.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer pt2prp(26),i

      double precision lgk

      logical table 

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision v,tr,pr,r,ps
      common/ cst5 /v(l2),tr,pr,r,ps

      integer inc,jpot
      common/ cst101 /inc(l2),jpot

      character*8 vname,xname
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      save pt2prp

      data pt2prp/11, 2,15,12, 1,13,14,17,10, 3, 4,18,20, 5,19,21,
     *             6,22,25, 7,23,26, 8,24,27, 28/
c----------------------------------------------------------------------

      lgk = -(psys(11)/r/v(2))/dlog(1d1)

      if (table) then 

         write (n4,1000) (v(iv(i)), i = 1, jpot),
     *                   psys(pt2prp(1)),psys(pt2prp(2)),lgk,
     *                   (psys(pt2prp(i)),i=3,26)

      else 

         write (*,1050)
         write (*,1030) (vname(iv(i)),v(iv(i)), i = 1, ipot)
         write (*,1010) psys(pt2prp(1))/1d3,psys(pt2prp(2))/1d3,lgk,
     *                  (psys(pt2prp(i)),i=3,7)

         if (.not.rxn) then 
            write (*,1020) psys(pt2prp(26)),(psys(pt2prp(i)),i=8,25)
         else 
            write (*,1040) 
         end if

      end if 

1000  format (40(g14.7,1x))
1010  format (/,'apparent Gibbs energy (kJ/mol) = ',g14.7,
     *        /,'apparent enthalpy (kJ/mol) ',t32,'= ',g14.7,
     *        /,'log10[Keq] ',t32,'= ',g14.7,/,
     *        /,'entropy (J/mol/K) ',t32,'= ',g14.7,
     *        /,'heat capacity (J/mol/K) ',t32,'= ',g14.7,/,
     *        /,'volume (J/mol/bar) ',t32,'= ',g14.7,
     *        /,'expansivity (1/K) ',t32,'= ',g14.7,
     *        /,'compressibility (1/bar) ',t32,'= ',g14.7)
1020  format (/,'heat capacity ratio (cp/cv)',t32,'= ',g14.7,
     *        /,'formula weight (g/mol) ',t32,'= ',g14.7,
     *        /,'density (kg/m3) ',t32,'= ',g14.7,/,
     *        /,'Gruneisen_T ',t32,'= ',g14.7,//,
     *        'Adiabatic elastic moduli:',/,
     *        t30,' T derivative',t45,' P derivative',/
     *        2x,'Ks(bar) = ',g14.7,t30,g14.7,t45,g14.7,/,
     *        2x,'Gs(bar) = ',g14.7,t30,g14.7,t45,g14.7,//,
     *        'Seismic velocities:',/,
     *        t30,' T derivative',t45,' P derivative',/
     *        2x,'v0(km/s) = ',g14.7,t30,g14.7,t45,g14.7,/,
     *        2x,'vp(km/s) = ',g14.7,t30,g14.7,t45,g14.7,/,
     *        2x,'vs(km/s) = ',g14.7,t30,g14.7,t45,g14.7,//,40('-'),/)
1030  format (29x,a,' = ',g12.6)
1040  format (/,40('-'),/)
1050  format (/,40('-'),//,'Thermodynamic properties at:',/)

      end 
