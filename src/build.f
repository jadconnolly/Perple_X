      program build 
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *  build.mar.10.1991   *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------
c build is a fortran program which interactively creates the input
c file read from unit n1 by the vertex program.  build reads data
c from three sources n2, *, and n9.  the output file is written to
c unit n1
c-----------------------------------------------------------------------
c files (see vertex program documentation for additional information):
c
c  l.u.n  i/o       purpose
c  -----  ---   --------------------------------------------------------
c   n2     i    data file containing the names, compositional vectors,
c               and standard state thermodynamic data for potentially
c               stable condensed phases and fluid species.
c   n9     i    optional data file which contains data on the
c               compositional equations of state and subdivision schemes
c               to be used for solution phases.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision c(0:4)

      character mnames(k16*k17)*8, n9name*100, dsol*100, opname*100,
     *          mname(k5)*5, oname(k5)*5, pname(k5)*5, cfname*100,
     *          uname(k0)*5, group*28, aname(i9)*6, amount*5, new*3, 
     *          text*256, dtext*200, title*162, y*1, lname(i9)*22,
     *          sname(i9)*10, tn1*6, tn2*22, fname(i9)*10

      integer i, k, l, iind, im, idum, ivct, jcth, j, ier, idep, 
     *        gct(i9), gid(i9,i9), ict, idsol, inames

      logical eof, good, oned, findph, first, feos, chksol

      external chksol, findph

      character prject*100,tfname*100
      common/ cst228 /prject,tfname

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      character*8 name
      common/ csta6 /name

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)
 
      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k0),ictr(k0),itrans

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer eos
      common/ cst303 /eos(k10)

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      double precision buf
      common/ cst112 /buf(5)

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx
  
      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ic
      common/ cst42 /ic(k0)

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m1),iterm,iord,istot,jstot,kstot

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer iam
      common/ cst4 /iam

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      data dsol/'solution_model.dat'/
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 4
c                                 version info
      call vrsion (6)
c                                 initialize strings (necessary for some OS). 
      opname = ' '
      n9name = ' '
      cfname = ' '
      title = ' '

      write (*,7020)
c                                 name and open computational option file (unit n1)      
      call fopen1 
c                                 choose and open thermo data file (unit n2)
      call fopen2 (1,opname)
c                                 get computational option file name
      write (*,1170) 
      read (*,'(a)') opname
 
      if (opname.eq.' ') opname = 'perplex_option.dat'
c                                 initialization:
      do i = 1, l2
         iv(i) = i
         vmin(i) = 0d0
         vmax(i) = 0d0
         dv(i) = 0d0 
      end do 

      iv(5) = 3
      idum = 0
      aqst = -1
      ixct = 0
      isat = 0
      jmct = 0
      ifct = 0
      io3 = 1
      ivct = 2
      icont = 0
      iphct = 0
      isoct = 0 
      jmuct = 0
      first = .true. 

      call redop1 (.false.,opname)

      do i = 0, 4
         c(i) = 0d0
      end do 
c                                 Read THERMODYNAMIC DATA file header
      call topn2 (3)
c                                 ====================================
c                                 choose problem class
      write (*,1490)

      call rdnumb (c(0),0d0,icopt,2,.false.)

      if (icopt.eq.7) then
c                                 warn for frac2d calculations about the
c                                 aux file and set icopt, icopt may be 
c                                 changed to 10 in varich if fileio.
         icopt = 9
         call mertxt (tfname,prject,'.aux',0)
         write (*,1500) tfname

      end if
c                                 choose chemical components
      call compch (ivct,feos,mname,pname,oname,uname)
c                                 physical variable choices and ranges
c                                 set icopt to its internal value:
      call varich (c,ivct,iind,oned,idep,iord,jcth,
     *             amount,dtext,opname,pname,cfname)
c                                 warn about the use of chemical potentials
c                                 in different types of calculations
      if (jmct.gt.0) then 
c                                 count variables
         j = 0

         do i = 1, 2
            if (iv(i).gt.3) j = j + 1
         end do

         if (icopt.eq.0.or.icopt.eq.5.or.
     *       (j.eq.0.and.icopt.eq.1.or.icopt.eq.3)) then 

            write (*,3110)

         else if (icopt.eq.1.or.icopt.eq.3) then 
c                                 an spd or mixed variable diagram as
c                                 a function of chemical potentials
            write (*,3111)

         end if

      end if
c                                  ================================
c                                  ask if the user wants print output,
c                                  graphical output is automatic. 
      write (*,2000) 
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') then 
         io3 = 0
      else 
         io3 = 1
      end if 
c                                  warn about non-plottable results:
      if (icopt.eq.0.and.icp.ne.3) then
         call warn (37,0d0,i,'BUILD')
      else if (icopt.eq.3.and.icp.ne.2) then
         call warn (38,0d0,i,'BUILD')
      else if (jcth.gt.0.and.jcth.lt.icp) then
         call warn (39,0d0,i,'BUILD')
      end if 
c                                 ======================================
c                                 now start phase data:
c                                 get composition vectors for entities
c                                 defined by a make definition:
      call makecp (inames,mnames,first)

      call eohead (n2)

      write (*,*) ' '
c                                 next get consistent real phases
      do 
 
         call getphi (name,.true.,eof)

         if (eof) exit
 
         call chkphi (0,name,good)

         if (.not.good) cycle
 
         iphct = iphct + 1

         call loadit (iphct,.false.,.true.)
 
      end do 
c                                 get all made phases
      do i = 1, nmak
         iphct = iphct + 1
         names(iphct) = mknam(i,mknum(i)+1)
      end do 
c                                 Excluded phases:
      write (*,2080)
      read (*,'(a)') y
 
      if (y.eq.'Y'.or.y.eq.'y') then
 
         write (*,2034)
         read (*,'(a)') y
 
         if (y.ne.'y'.and.y.ne.'Y') then

            name = 'b'

            write (*,2021)
 
            do while (ixct.lt.h8.and.name.ne.' ') 

               read (*,'(a)') name

               if (name.eq.' ') exit 

               good = .false.
 
               do i = 1, iphct
                  if (names(i).eq.name) then 
                     ixct = ixct + 1
                     exname(ixct) = name
                     good = .true.
                     exit 
                  end if 
               end do 

               if (ixct.eq.h8) then 
                  call warn (8,0d0,h8,'BUILD')
                  exit 
               end if 

               if (good) cycle

               write (*,1020) name
 
            end do
            
         else

            do i = 1, iphct
 
               write (*,1130) names(i)
               read (*,'(a)') y
 
               if (y.eq.'Y'.or.y.eq.'y') then
                  ixct = ixct + 1
                  if (ixct.gt.h8) then 
                     call warn (8,0d0,h8,'BUILD')
                     ixct = ixct - 1
                     exit
                  end if 
                  exname(ixct) = names(i)
              end if
 
           end do 
 
        end if
 
      end if
c                                check for fluid species EoS if no special components:
c                                as of 3/16/2016 this was incorrect, so far as i can 
c                                an EoS is automatically associated with fluid species
c                                if eos(i) = 10 (ideal) gas, or 117 > eos(i) > 103, therefore
c                                101 (H2O) and 102 (CO2) are the only special cases that
c                                require the user to specify an EoS

c                                for true ternary c-o-h fluid calculations H2O and CO2
c                                should be excluded and this loop should not request that
c                                the user specify an EoS
      good = .true.

      if (.not.feos) then
 
         do i = 1, iphct

            if (eos(i).eq.201.or.eos(i).eq.202) then 

               do j = 1, ixct

                  if (exname(j).eq.names(i)) then
                     good = .false.
                     exit 
                  end if  

               end do 

               if (.not.good) cycle
c                               got one, get the EoS
               write (*,'(/,a,/)') 
     *               'One or more endmembers require a fluid EoS'

               call rfluid (1) 

               exit 

            end if 

         end do 

      end if 
c                                read solution phases.
      write (*,2500)
      read (*,'(a)') y
 
      if (y.eq.'y'.or.y.eq.'Y') then
c                                 get the file containing the solution models
         do
 
            write (*,3010)
            read (*,'(a)') n9name
            if (n9name.eq.' ') n9name = dsol
            write (*,*) 
            open (n9,file=n9name,iostat=ier,status='old')
 
            if (ier.ne.0) then
c                                 system could not find the file
               write (*,3020) n9name
               write (*,7050)
               read (*,'(a)') y
 
               if (y.ne.'Y'.and.y.ne.'y') goto 999
               cycle
c                                 try again
            end if
 
            ict = 0
            ipoint = iphct
c                                 format test line
            read (n9,*) new
c                                 check version compatability
            if (chksol(new)) exit 

            call warn (3,0d0,i,new)

            cycle 

         end do 

         do 
c                                 read candidates:
            call rmodel (tn1,tn2)
c                                 istot = 0 = eof
            if (istot.eq.0) exit 
c                                 don't allow fluid models if 
c                                 the system is fluid saturated:
            if (jsmod.eq.0.and.ifct.gt.0) cycle
c                                 check for endmembers:
            call cmodel (im,idsol,first,good)

            if (jsmod.eq.39) then
               if (jstot.eq.0) cycle
            else if (jstot.lt.2) then 
               cycle
            end if

            ict = ict + 1
            if (ict.gt.i9) call error (24,0d0,i9,'build')

            fname(ict) = tname
            aname(ict) = tn1
            lname(ict) = tn2

         end do 
c                                 group models according to fullname
         j = 0

         do i = 1, ict
            gct(i) = 0 
         end do 

         do i = 1, ict 

            good = .false.

            do k = 1, j

               if (lname(i).eq.lname(gid(k,1))) then 
c                                 in group j
                  gct(k) = gct(k) + 1
                  gid(k,gct(k)) = i
                  good = .true.
                  exit

               end if 

            end do 

            if (good) cycle
c                                 in a new group
            j = j + 1
            gct(j) = 1
            gid(j,gct(j)) = i

         end do 
c                                 we have the list, ask user for choices
         if (ict.eq.0) then
c                                 no valid models
            write (*,7040)
 
         else 

            write (*,2510)

            if (j.eq.1) then 
c                                 just one group
               write (*,'(6(2x,a))') (fname(i), i = 1, ict)

            else

               do i = 1, j

                  call mertxt (group,lname(gid(i,1)),' models:',1)

                  l = 5
                  if (gct(i).lt.6) l = gct(i)

                  write (*,'(a25,5(1x,a))') group,
     *                                 (fname(gid(i,k)), k = 1, l)

                  if (gct(i).lt.6) cycle 

                  do k = 6, gct(i)-4, 5
                     write (*,'(t26,5(1x,a))') 
     *                                 (fname(gid(i,l)), l = k, k+4)
                  end do
 
                  if (k-5.lt.gct(i)) write (*,'(t26,5(1x,a))')
     *                                 (fname(gid(i,l)), l = k, gct(i))

               end do 

            end if 

            write (*,1180) 
c                                 initialize tname, because reading a 
c                                 null record may not do that. 
            tname  = ' '

            do 

               read (*,'(a)') tname
               if (tname.eq.' ') exit
c                                 check if same name entered twice
               do i = 1, isoct
                  if (tname.eq.sname(i)) cycle
               end do 
c                                 check if name in list
               good = .false.

               do i = 1, ict
                  if (tname.eq.fname(i)) then
                     isoct = isoct + 1
                     if (isoct.gt.h9) call error (25,0d0,h9,'BUILD') 
                     sname(isoct) = tname
                     good = .true.
                     exit 
                  end if 
               end do 

               if (good) cycle

               write (*,2310) tname

            end do  
         end if

      else 

         n9name = ' '

      end if 
c                                 get title
      write (*,7070) 
      read (*,'(a)') title
c                                 output options etc:
c                                 thermodynamic data file is output to n1 by fopen
c                                 print output request:
      if (io3.eq.0) then 
         write (n1,'(a)') 'print    | no_print suppresses print output'
      else 
         write (n1,'(a)') 'no_print | print generates print output'
      end if
c                                 plot output request  
      write (n1,'(a)') 'plot     | obsolete 6.8.4+'
c                                 solution model file request
      call mertxt (text,n9name,
     *            '| solution model file, blank = none',5)
      write (n1,'(a)') text
c                               if special dependencies, put them in the
c                               title
      if (title.eq.' '.and.idep.ne.0) then 
         text = dtext
      else if (title.eq.' ') then
         text = 'Title Text'
      else if (idep.eq.0) then 
         text = title
      else
         call mertxt (text,title,dtext,5)
      end if 
c                                 title
      write (n1,'(a)') text
c                                 option file name
      call mertxt (text,opname,'| Perple_X option file',5)
      write (n1,'(a)') text 
c                                 computational mode:
      write (n1,1010) icopt,'calculation type: 0- composition,'//
     *                ' 1- Schreinemakers, 3- Mixed, 4- swash,'//
     *                ' 5- gridded min, 7- 1d fract, 8- gwash,'//
     *                ' 9- 2d fract, 10- 7 w/file input,'//
     *                ' 11- 9 w/file input, 12- 0d infiltration'
c                                 coordinate file name if necessary
      if (icopt.eq.10.or.icopt.eq.11) then 
         call mertxt (text,cfname,'| coordinate file',5)
         write (n1,'(a)') text
      end if 

      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06' 
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) itrans,'number component transformations'
      write (n1,1010) icmpn,'number of components in the data base'

      if (itrans.gt.0) then
         do i = 1, itrans
            write (n1,1120) uname(ictr(i)), ictr(i)
            write (n1,1125) (ctrans(j,i), j = 1, icmpn)
         end do 
      end if 

      write (n1,1010) iwt,'component amounts, 0 - mole, 1 mass'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 06'
      write (n1,1010) idum,'unused place holder, post 05'
c                                 output saturated phase eos choice:
      write (n1,1010) ifug,'ifug EoS for saturated phase'
      if (ifug.eq.8 .or.ifug.eq.10.or.ifug.eq.12.or.ifug.eq.16.or.
     *    ifug.eq.17.or.ifug.eq.19.or.ifug.eq.20.or.ifug.eq.24.or.
     *    ifug.eq.25)  write (n1,1330) ibuf, hu, dlnfo2, elag,
     *              'ibuf, ipro, choice dependent parameter, ln(ag)'
      if (ibuf.eq.5) write (n1,4020) buf,'a-e'

      if (oned) then 
         i = 1
      else 
         i = 2
      end if 

      write (n1,1010) i,'gridded minimization dimension (1 or 2)'
c                                 potential variable dependencies
      write (n1,1010) idep,'special dependencies: '//
     *                '0 - P and T independent, 1 - P(T), 2 - T(P)'
      write (n1,1560) c
c                                 output component data:
      write (n1,3060) 'thermodynamic component list'

      do i = 1, icp 
         if (i.gt.jcth) then 
            write (n1,3000) pname(i),0,0.,0.,0.,'unconstrained'
         else 
            write (n1,3000) pname(i),icont,(dblk(j,i),j=1,3),amount
         end if 
      end do 

      write (n1,3050) 'thermodynamic component list'

      write (n1,3060) 'saturated component list'

      do i = icp + 1, icp + isat

         if (i.gt.jcth) then 
            write (n1,3000) pname(i),0,0.,0.,0.,'unconstrained'
         else 
            write (n1,3000) pname(i),icont,(dblk(j,i),j=1,3),amount
         end if

      end do 

      write (n1,3050) 'saturated component list'

      write (n1,3060) 'saturated phase component list'
      do i = 1, ifct
         write (n1,'(a)') mname(i)
      end do 
      write (n1,3050) 'saturated phase component list'

      write (n1,3060) 'independent potential/fugacity/activity list'
      do i = 1, jmct
         write (n1,'(3(a,1x))') oname(i),vname(3+i),afname(i)
      end do 
      write (n1,3050) 'independent potential list'

      write (n1,3060) 'excluded phase list'
      do i = 1, ixct
         write (n1,'(a)') exname(i)
      end do 
      write (n1,3050) 'excluded phase list'

      write (n1,3060) 'solution phase list'
      do i = 1, isoct
         write (n1,'(a)') sname(i)
      end do 
      write (n1,3050) 'solution phase list'
c                                 output variable choices and values:
      write (n1,4020) (vmax(i), i= 1, l2),'max p, t, xco2, mu_1, mu_2'
      write (n1,4020) (vmin(i), i= 1, l2),'min p, t, xco2, mu_1, mu_2'
      write (n1,4020) (dv(i), i= 1, l2),'unused place holder post 06'

      if (oned) then 
         write (n1,1310) (iv(i), i = 1, l2),
     *  'index of the independent & sectioning variables'
      else if (icopt.ne.0) then 
         write (n1,1310) (iv(i), i = 1, l2),
     *  'indices of 1st & 2nd independent & sectioning variables'
      else
         write (n1,1310) (iv(i), i = 1, l2),
     *  'independent variables indices'
      end if 
c                                 get conditions for composition
c                                 diagrams:
      if (icopt.eq.0) then

         if (ifct.eq.1) write (*,7150) 'fluid', vname(3) 

         i = 0

         do 

            i = i + 1
 
            write (*,6020) (vname(iv(j)),j=1,ivct)

            do 
               write (*,6010) i
               read (*,*,iostat=ier) (vmin(iv(j)),j=1,ivct)
               if (ier.eq.0) exit
               call rerr
            end do

            if (vmin(iv(1))+vmin(iv(2)).eq.0d0) exit 
            write (n1,1340) vmin

         end do 

      else if (icopt.eq.12) then 

         write (*,'(a)') 'Enter the number of water aliquots to be '//
     *                   'added to the system:'
         read (*,*) i

         write (*,'(a)') 'Enter the aliquot size (mole):'
         read (*,*) v(1)

         write (n1,'(i4,1x,g14.6,1x,a)')  i, v(1),
     *                         'Number and size (mole) of fluid '//
     *                         'aliquots for 0d infilitration model'

      end if
 
99    endfile (n1)
 
      close (n1)
 
999   stop

1010  format (i5,1x,a,a,a,a,a)
1020  format (a,' does not exist in the selected data base')
1120  format (a,1x,i2,' component transformation')
1125  format (13(f6.2,1x))
1130  format ('Exclude ',a,' (Y/N)?')
1170  format (/,'Enter the computational option file name ',
     *       '[default = perplex_option.dat]:',/,
     *       'See: www.perplex.ethz.ch/perplex_options.html')
1180  format (/,'For details on these models see: ',
     *        'www.perplex.ethz.ch/perplex_solution_model_glossary.html'
     *      /,'or read the commentary in the solution model file.',/)
1310  format (/,5(i2,1x),2x,a,/)
1330  format (i2,1x,i2,1x,g13.6,1x,g13.6,1x,a)
1340  format (5(g13.6,1x))
1490  format ('Specify computational mode:',//,
     *    5x,'1 - Convex-Hull minimization',/, 
     *    5x,'2 - Constrained minimization on a 2d grid [default]',/,
     *    5x,'3 - Constrained minimization on a 1d grid',/,
     *    5x,'4 - Output pseudocompound data',/,
     *    5x,'5 - 1-d Phase fractionation',/,
     *    5x,'6 - 0-d Infiltration-reaction-fractionation',/,
     *    5x,'7 - 2-d Phase fractionation (FRAC2D and TITRATE ',
     *               'reactive transport models)',//,
     *        'Use Convex-Hull minimization for Schreinemakers ',
     *        'projections or phase diagrams',/,
     *        'with > 2 independent variables. ',
     *        'Use constrained minimization for phase diagrams',/,
     *        'or phase diagram sections ',
     *        'with < 3 independent variables.')
1500  format (/,'Reactive transport models require an auxilliary input '
     *       ,'file ',a30,/,'that you must create yourself after ',
     *        'running BUILD. See:',//,5X,
     *        'www.perplex.ethz.ch/perplex/Perple_X_FRAC2D.html')
1560  format (5(g12.6,1x),'Geothermal gradient polynomial coeffs.')
2000  format (/,'Output a print file (Y/N)?')
2021  format ('Enter names, 1 per line, press <enter> to finish:')
2034  format ('Do you want to be prompted for phases (Y/N)?')
2080  format (/,'Exclude pure and/or endmember phases (Y/N)?')
2110  format (/,'Calculations with saturated components (Y/N)?')
2310  format (/,a,' is invalid. Check spelling, upper/lower case match',
     *        ', and do not use leading blanks. Try again:',/)
2500  format (/,'Include solution models (Y/N)?')
2510  format (/,'Select models from the following list, enter',
     *        ' 1 per line, press <enter> to finish',/)
3000  format (a,1x,i1,1x,3(g12.6,1x),a,' amount')
3010  format ('Enter the solution model file name [default = ',
     *        'solution_model.dat]: ')
3020  format (/,'**error ver191** FOPEN cannot find file:',/,a,/)
3050  format ('end ',a,/)
3060  format (/,'begin ',a)
3110  format (/,'**warning ver110** in this mode Perple_X  will not che'
     *       ,'ck whether conditions are',/,'supersaturated with resp',
     *        'ect to mobile components.',//,
     *        'To explicitly compute saturation surfaces:',/,
     *        '  1) use convex hull optimization and do not use ',
     *        'activities or fugacities as',/,
     *        '     independent variables',/,
     *        '  2) or compute the saturation surface with FLUIDS or ',
     *        'FRENDLY.',/)
3111  format (/,'**warning ver111** before running this calculation ',
     *        'understand the null_phase',/,
     *        'perplex option. If null_phase is true and the initial ',
     *        'condition for a calculation',/,
     *        'is supersaturated, then CONVEX will fail.',/)
4020  format (5(g14.6,1x),a)
6020  format (/,'Specify values for:',/,(10x,5(a,2x)))
6010  format ('For calculation ',i2,', enter zeros to finish.')
7020  format (//,'NO is the default (blank) answer to all Y/N prompts',
     *        /)
7040  format (/,'The solution model file contains no relevant models.',/
     *       )
7050  format ('Try again (Y/N)?')
7070  format (/,'Enter calculation title:')
7150  format (/,'*Although only one component is specified for the ',a,
     *       ' phase, its equation of state ',/,
     *       'permits use of its compositional variable: ',a,'.',/) 

      end
 
      subroutine depend (ivct,idep,iind,iord,c,dtext)
c---------------------------------------------------------------------------
c subroutine to reset variable flags and counters if one variable is
c made dependent on another; the routine also prompts and reads the 
c functional dependence
c---------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character xname*8, dtext*200, vname*8

      integer i, j1, j2, ivct, iind, idep, iord, ier

      double precision c(0:4)
 
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c                                 reset the variable counters and flags:
      ivct = ivct - 1

      if (idep.eq.1) then 
         iind = 2
         j1 = 1
         j2 = 2
         do i = 1, ivct 
            iv(i) = iv(i+1)
         end do 
      else 
         iind = 1
         j1 = 2
         j2 = 1
         do i = 2, ivct 
            iv(i) = iv(i+1)
         end do 
      end if 

      iv(ivct+1) = idep
c                                  get the functional dependence
      do

         write (*,1070) vname(j1),vname(j2)
         read (*,*,iostat=ier) iord

         if (ier.ne.0) then
            call rerr
            cycle 
         end if

         if (iord.lt.0.or.iord.gt.5) cycle 
         
         exit

      end do

      do i = 0, iord

         do 
            write (*,1080) i
            read (*,*,iostat=ier) c(i)
            if (ier.eq.0) exit
            call rerr
         end do

      end do 

      write (*,'(/)')
c                                  write a text version for the title
      write (dtext,1580) vname(idep),c(0),
     *                      (c(i),vname(iind),i,i=1,iord)
      call deblnk (dtext)

1080  format ('Enter c(',i2,')')
1070  format (/,'The dependence must be described by the polynomial',//,
     *        a,' = Sum ( c(i) * [',a,']^i, i = 0..n)',//,
     *       'Enter n (<5)')
1580  format (a,' = ',g12.6,4(' + ',g12.6,' * ',a,'^',i1))

      end 

      subroutine redvar (ind,iprompt)
c----------------------------------------------------------------------
c redvar interactively reads and checks values for the primary variable 
c indexed by ind.
c
c if iprompt = 1 reads vmin(ind) and vmax(ind) are read.
c else reads just vmin(ind)
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ind, iprompt, icount, ier

      logical numbad

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------
      do 

         if (iprompt.eq.1) then 
            icount = 2
            write (*,1010) vname(iv(ind))
         else if (iprompt.eq.2) then
            icount = 1
            write (*,1020) vname(iv(ind))
         else if (iprompt.eq.3) then 
            icount = 1
            write (*,1030) vname(iv(ind))
         end if 

         if (icount.eq.1) then 
            read (*,*,iostat=ier) vmin(iv(ind))
         else 
            read (*,*,iostat=ier) vmin(iv(ind)),vmax(iv(ind))
         end if 

         if (ier.ne.0) then 
            write (*,1000)
            cycle 
         end if 

         if (numbad(1,ind)) cycle 

         if (icount.eq.2) then 
            if (numbad(2,ind)) cycle  
         end if 

         exit 
    
      end do 

1000  format (/,' Your input is incorrect, probably you are using ',
     *        'a character where',/,' you should be using a number ',
     *        'or vice versa, try again...',/)
1010  format (/,'Enter minimum and maximum values, respectively,',
     *        ' for: ',a)
1020  format (/,'Specify sectioning value for: ',a)
1030  format (/,'Specify the value for: ',a)

      end 

      logical function numbad (num,ind)
c----------------------------------------------------------------------
c numbad checks if a primary variable limit is reasonable, the variable
c is indexed by ind and is the lower (vmin) or upper (vmax) bound 
c depending on the value of num (1 or 2, respectively). 
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character*1 y

      integer num, ind, jnd

      double precision value 

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------

      numbad = .false.

      jnd = iv(ind) 

      if (num.eq.1) then 
         value = vmin(jnd)
      else
         value = vmax(jnd)
      end if 

      if (jnd.eq.1.or.jnd.eq.2) then 
c                                 pressure-temperauture
         if (value.le.0d0) then 
            call warn (57,value,jnd,vname(jnd))
            numbad = .true.
         end if 

      else if (jnd.eq.3) then 
c                                 phase composition
         if (value.lt.0d0.or.value.gt.1d0) then
            call warn (44,value,jnd,vname(jnd))
            numbad = .true.
         end if 
 
      else 
c                                 potential/fugacity/activity
         if (imaf(jnd-3).eq.3.and.value.gt.0d0) then
c                                 imaf = 3, warn on activity > 1.
            call warn (43,value,jnd,vname(jnd))
            numbad = .true.
         end if 

      end if 

      if (numbad) then 
         read (*,'(a)') y
         if (y.ne.'y'.and.y.ne.'Y') numbad = .false.
      end if  

      end 

      subroutine chknam (igood,jcmpn,iflu,good,char5,qname,uname)
c----------------------------------------------------------------------
c chknam - looks for a match between string char5 and the jcmpn strings in 
c array qname. if the match is found good = .true, index is the index of
c the string in uname, and the string is eliminated from uname and jcmpn
c decremented. igood is the index of the equivalent string in the uname
c array.

c if iflu = 0, then chknam first matches fluid components in uname
c before eliminating the component from qname.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      character*5 qname(k0), uname(k0), char5

      logical good, go

      integer igood, jcmpn, iflu, i, j

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
c----------------------------------------------------------------------

      good = .false.

      if (iflu.eq.0) then

         go = .true.

         do i = 1, ispec
            if (char5.ne.uname(idspe(i))) cycle
            go = .false.
         end do 
c                                 special check for saturated phase
c                                 components. 
         if (go) then 
            write (*,1000) char5
            return
         end if

      end if 

      do i = 1, jcmpn
         if (qname(i).eq.char5) then 
c                                  eliminate used components
            do j = i+1, jcmpn
               qname(j-1) = qname(j)
            end do 
            good = .true.
            exit 
         end if 
      end do

      if (good) then 
c                                 decrement qname counter
         jcmpn = jcmpn - 1
c                                 find index in uname array
         do i = 1, icmpn
            if (char5.ne.uname(i)) cycle 
            igood = i
            exit
         end do 

      else 

         write (*,1000) char5

      end if 

1000  format (/,a,' is invalid. Check spelling, upper/lower case match',
     *        ', and do not use leading blanks. Try again:',/)
      end

      subroutine compch (ivct,feos,mname,pname,oname,uname)
c---------------------------------------------------------------------------
c interactively choose components for build.
c---------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, jcmpn, ivct, igood, ier, jsat(h5), ima, jspec, jc(2)

      logical satflu, mobflu, good, findph, eof, quit, feos

      character mname(*)*5, qname(k0)*5, char5*5, uname(*)*5,
     *          y*1, oname(*)*5, nname(k5)*5, char6*6, 
     *          pname(*)*5, fugact(3)*8

      external findph

      integer ic
      common/ cst42 /ic(k0)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      character*8 name
      common/ csta6 /name

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      data fugact/'chem_pot','fugacity','activity'/

c---------------------------------------------------------------------------
      jcmpn = icmpn
      satflu = .false.
      mobflu = .false.
      feos = .false.
      char5 = ' '
c                                 Component stuff first:
      do i = 1, icmpn
         uname(i) = cmpnt(i)
         qname(i) = cmpnt(i)
         pname(i) = cmpnt(i)
      end do

      if (lopt(7).and.icopt.ne.9) then 
c                                 components of saturated phase:
         write (*,2030) 
         read (*,'(a)') y

         if (y.eq.'y'.or.y.eq.'Y') then
c                                   write prompt
            write (*,2031) 
            write (*,'(12(1x,a))') (uname(idspe(i)), i = 1, ispec)
            write (*,2021)
c                                   write blurb
            write (*,1030)

            do 

               read (*,3000) char5 

               if (char5.ne.' ') then
c                                 check choice 
                  call chknam (igood,jcmpn,1,good,char5,qname,uname)

                  if (.not.good) cycle
c                                 count component
                  ifct = ifct + 1
                  jc(ifct) = igood
                  mname(ifct) = char5

                  if (ifct.eq.2) then 
                     write (*,*)
                     exit 
                  end if   

               else
c                                 blank input 
                   exit 

               end if 

            end do 

            if (ifct.ne.0) then
               iv(3) = 3
               ivct = 3
            end if

         end if

      end if
c                                 number of used special components
      jspec = ifct 

      if (icopt.ne.9) then 
c                                 saturated components
         write (*,2110)
         read (*,'(a)') y
 
         if (y.eq.'Y'.or.y.eq.'y') then

            call warn (15,nopt(1),i,'BUILD')
 
            write (*,2032) h5+1
            write (*,'(12(1x,a))') (qname(i),i=1,jcmpn)
            write (*,2021)

            do

               read (*,3000) char5 

               if (char5.ne.' ') then
c                                 check if it's a saturated phase component:
                  good = .false.

                  do i = 1, ispec
                     if (char5.ne.uname(idspe(i))) cycle
                     good = .true.
                  end do 

                  if (good) then 

                     if (ifct.gt.0) then 
                        call warn (7,nopt(1),i,'BUILD')
                        write (*,1000) 
                        read (*,'(a)') y
                        if (y.eq.'y'.or.y.eq.'Y') then
                           satflu = .true.
                           write (*,1070)
                        else 
                           write (*,1080) 
                           cycle
                        end if
                     else 
                        satflu = .true.
                     end if  

                  end if 

                  call chknam (igood,jcmpn,1,good,char5,qname,uname)

                  if (.not.good) cycle
c                                 count component
                  isat = isat + 1    
                  jsat(isat) = igood
                  nname(isat) = char5

                  if (satflu) then
                     jspec = jspec + 1
                     mname(jspec) = char5
                  end if 

               else 
c                                 blank input
                  exit   
             
               end if 

            end do 
 
         end if
c                                 mobile components:
         write (*,2040)
         read (*,'(a)') y
 
         if (y.eq.'y'.or.y.eq.'Y') then
         
            do 
c                                 write prompt
               write (*,2050) 
               if (jmct.eq.1) write (*,2055)
               write (*,'(12(1x,a))') (qname(i),i=1,jcmpn)
               read (*,3000) char5 
               if (char5.eq.' ') exit 
c                                 check if it's a saturated phase component,
c                                 or a saturated component:
               good = .false.

               do i = 1, ispec
                  if (char5.ne.uname(idspe(i))) cycle
                  good = .true.
               end do 

               if (good) then

                  if (ifct.gt.0.or.satflu) then  
                     call warn (5,nopt(1),i,'BUILD')
                     write (*,1000) 
                     read (*,'(a)') y
                     if (y.eq.'y'.or.y.eq.'Y') then
                        mobflu = .true.
                        write (*,1070)
                     else 
                        write (*,1080) 
                        cycle
                     end if
                  else
                     mobflu = .true.
                  end if 

               end if 

               call chknam (igood,jcmpn,1,good,char5,qname,uname)

               if (.not.good) cycle
c                                 count component
               jmct = jmct + 1    
               ivct = ivct + 1
               oname(jmct) = char5
c                                 ask if chemical potential, activity
c                                 or fugacity
               write (*,2060) char5
               do 
                  read (*,*,iostat=ier) ima
                  if (ier.eq.0) exit
               end do 
c                                 if fugacity/activity get the reference
c                                 phase
               if (ima.eq.2.or.ima.eq.3) then
                     
                  call topn2 (2)

                  do 
 
                     call getphi (name,.false.,eof)

                     if (eof) exit

                     if (findph(igood)) then
                        iphct = iphct + 1
                        names(iphct) = name
                     end if 

                  end do 
c                                 names contains the list of candidates
                  if (iphct.eq.0) then

                     write (*,2061) fugact(ima),char5
                     ima = 1
                     imaf(jmct) = ima
                     jmuct = jmuct + 1

                  else if (iphct.eq.1) then 

                     write (*,2062) names(1),char5,fugact(ima)
                     afname(jmct) = names(1)
                     icp = 1

                  else 

                     write (*,2063) char5,fugact(ima)
                     write (*,2064) (names(i),i=1,iphct)
                     good = .false.
                     do 
                        read (*,'(a)') name
                        do i = 1, iphct
                           if (name.eq.names(i)) then
                              good = .true.
                              icp = i  
                              afname(jmct) = name
                              exit 
                           end if
                        end do 
                        if (good) exit
                        write (*,2310) name
                     end do

                     write (*,1190) 

                  end if

               else

                  ima = 1
                  afname(jmct) = ' '

               end if
c                                 now make the variable name
               if (ima.eq.2.or.ima.eq.3) then

                  read (names(icp),*) char6
                  if (ima.eq.2) then 
                     write (vname(3+jmct),'(a,a)') 'f_',char6
                  else 
                     write (vname(3+jmct),'(a,a)') 'a_',char6
                  end if

                  write (*,2065) char5,fugact(ima),vname(3+jmct)
c                                 if name is a special phase component 
c                                 set flag for eos request
                  do i = 1, ispec
                     if (names(icp).ne.uname(idspe(i))) cycle
                     jspec = jspec + 1
                     mname(jspec) = uname(idspe(i))
                  end do 

               else 

                  write (vname(3+jmct),'(a,a)') 'mu_',char5
                  write (*,2066) char5,vname(3+jmct)
c                                 if component is a special phase component 
c                                 set flag for eos request
                  do i = 1, ispec
                     if (char5.ne.uname(idspe(i))) cycle
                     jspec = jspec + 1
                     mname(jspec) = uname(idspe(i))
                  end do 

               end if 
                
               imaf(jmct) = ima
               iphct = 0 
               if (jmct.eq.i6) exit  

            end do 
c                                 reset n2 
            call topn2 (2)
c                                 end of mobile component loop
         end if
 
         if (ifct.eq.0) then 
c                                 no saturated fluid phase 
            iv(3)=4
            iv(4)=5

         else 
c                                 a saturated phase is present
c                                 but may have only one component
            iv(4)=4
            iv(5)=5

         end if
      end if 
c                                 Thermodynamic components:
      write (*,2070) 
      write (*,'(12(1x,a))') (qname(i),i=1,jcmpn)
      write (*,2021)

      icp = 0 

      do 

         read (*,3000) char5

         if (char5.ne.' ') then 
c                                 check if it's a saturated phase component,
c                                 mobile or a saturated component:
            if (ifct.gt.0.or.satflu.or.mobflu) then 

               good = .false.

               do i = 1, ispec
                  if (char5.ne.uname(idspe(i))) cycle
                  good = .true.
               end do 
               
               if (good) then 
 
                  call warn (5,0d0,i,'BUILD')
                  write (*,1000) 
                  read (*,'(a)') y
                  if (y.eq.'y'.or.y.eq.'Y') then
                     write (*,1070)
                  else 
                     write (*,1080) 
                     cycle
                  end if
 
               end if 

            end if 

            call chknam (igood,jcmpn,1,good,char5,qname,uname)

            if (.not.good) cycle
c                                  good name, count and save index
            icp = icp + 1

            if (icp.gt.k5) call error (197,0d0,icp,'BUILD')

            ic(icp) = igood
            pname(icp) = char5

         else 
c                                 blank input, check counter
            if (icp.eq.0) then 
               call warn (19,0d0,k5,'build')
               cycle
            end if

            exit 

         end if 

      end do
c                                 check for fluid components 
c                                 in thermodynamic composition space.
      do i = 1, icp

         do j = 1, ispec
            if (pname(i).ne.uname(idspe(j))) cycle
            jspec = jspec + 1
            mname(jspec) = pname(i)
         end do

      end do

      if (jspec.gt.0) then
c                                 get fluid equation of state
c                                 write explanation
         write (*,1010) (mname(i), i = 1, jspec)
         write (*,1020)

         if (ifct.ne.0) write (*,1025) 

         feos = .true.
         call rfluid (1)
         write (*,'(/)')
c                                 eliminate composition variable
c                                 for saturated fluid if constrained
c                                 fugacity EoS is used:
c                                 probably should change iv(3)?
         if ((ifug.eq.9.or.ifug.eq.24).and.ifct.gt.1) 
     *      ivct = ivct - 1

      end if
c                                 make pname into a list of retained (non-mobile) components
      do i = icp + 1, icp + isat
         pname(i) = nname(i-icp)
      end do
c                                 jcmpn unused components are in qname
c                                 get the indices for chkphi
      do i = 1, icmpn

         quit = .false.
         icout(i) = 0

         do j = 1, jcmpn

            if (qname(j).ne.uname(i)) cycle

            quit = .true.
            exit

         end do 

         if (quit) cycle 

         icout(i) = 1

      end do
c                                 also pad out ic array with saturated
c                                 component pointers for chkphi
      do i = 1, isat 
         ic(icp+i) = jsat(i)
      end do

      do i = 1, ifct 
         ic(icp+isat+i) = jc(i)
      end do 

1000  format ('Are you sure you want to do this (y/n)?')
1010  format ('Because the thermodynamic data file identifies:',
     *        5(1x,a5))
1020  format ('as special components, you will be prompted next for ',
     *        'the EoS to be used for ',/,'the ',
     *        'corresponding composants and mixtures thereof. To over',
     *        'ride this behavior, ',/,'e.g., to use a generic hybrid ',
     *        'fluid EoS, delete the special_component section ',/,
     *        'from the thermodynamic data file header.',/)
1025  format ('NOTE: the EoS choice specified here will override the ',
     *        'EoS choice specified by',/,'the hybrid_EoS option, to',
     *        'override this behavior delete the special_component',/,
     *        'section from the the thermodynamic data file header.',/)
1030  format (/,'For C-O-H fluids it is only necessary to select ',
     *       'volatile species present in',/,'the solids of interest. ',
     *       'If the species listed here are H2O and CO2, then to',/,
     *       'constrain O2 chemical potential to be consistent with ',
     *       'C-O-H fluid speciation',/,'treat O2 as a saturated ',
     *       'component. Refer to the Perple_X Tutorial for details.',/)
1070  format ('Ok, but don''t say i didn''t warn you.')
1080  format ('Wise move, choose another component.')
1190  format (/,'**warning ver064** in general it is wise to exclude ',
     *        'the unused phases from',/,'the above list from your ',
     *        'calcultion. NOTE: you will not be prompted for these ',/,
     *        'phases if you select the automatic phase exclusion ',
     *        'prompt. To exclude unused',/,'reference phases either',
     *        'do not select the prompt option or edit the problem',/,
     *        'definition after running BUILD',/)
2021  format ('Enter names, 1 per line, press <enter> to finish:')
2030  format (/,'Calculations with a saturated fluid (Y/N)?')
2031  format (/,'Select the independent saturated fluid components:')
2032  format (/,'Select < ',i1,' saturated components from the set:')
2040  format (/,'Use chemical potentials, activities or fugacities as',
     *        ' independent',/,'variables (Y/N)?')
2050  format (/,'Specify a component whose chemical potential, activi',
     *       'ty or fugacity is',/,'to be independent, ',
     *       'press <enter> to finish:')
2055  format (/,'In other words, if you only want one mobile component',
     *       ', then press <enter>',/)
2060  format (/,'Component ',a,' is to be characterized by:',//,
     *       5x,'1 - chemical potential [default]',/,
     *       5x,'2 - log10(fugacity)',/,
     *       5x,'3 - log10(activity)')
2061  format (/,'The data file contains no phase suitable to define ',
     *        'the',a,' of ',a,/,'This component will be ',
     *        'characterized by its chemical potential',/)
2062  format (/,'Phase ',a,' will be used to define ',a,' ',a)
2063  format (/,'Select the phase to be used to define ',a,' ',
     *        a,' from the following:',/)
2064  format (4(5x,9(a,1x),/))
2065  format (/,'The log10(',a,' ',a,') variable is named: ',a,/)
2066  format (/,'The chemical potential of ',a,' is named: ',a,/)
2070  format (/,'Select thermodynamic components from the set:')
2110  format (/,'Calculations with saturated components (Y/N)?')
2310  format (/,a,' is invalid. Check spelling, upper/lower case match',
     *        ', and do not use leading blanks. Try again:',/)
3000  format (a,1x,i1,1x,3(g12.6,1x),a,' amount')

      end 

      subroutine varich (c,ivct,iind,oned,idep,iord,jcth,
     *                   amount,dtext,opname,pname,cfname)
c---------------------------------------------------------------------------
c interatctively choose physical variables for build.
c---------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, ivct, ier, iind, idep, iord, kvct, jc, icth,
     *        jcth, loopx, loopy, ind, ix, jst, jvct

      logical oned

      character y*1, dtext*(*), amount*5, stext*11, nc(3)*2, 
     *          opname*(*), pname(*)*5, cfname*100

      double precision c(0:4)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      character prject*100,tfname*100
      common/ cst228 /prject,tfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)
c-----------------------------------------------------------------------
      nc(1) = 'C0'
      nc(2) = 'C1'
      nc(3) = 'C2'
      iwt = 0
      jcth = 0
      amount = 'molar '
      dtext = ' '
      cfname = ' '
c                                 coming in icopt has the following values

c                                 1 - convex
c                                 2 - 2d grid
c                                 3 - 1d grid
c                                 4 - swash
c                                 5 - 1d frac
c                                 6 - 0d frac
c                                 9 - frac2d

c                                 reorder for oned flag
      if (icopt.eq.3.or.icopt.eq.5.or.icopt.eq.9) then
 
         if (icopt.ne.9) oned = .true.

         if (icopt.ne.3.and.icp.eq.1) call error (53,0d0,i,'BUILD')
c                                 fractionation from a file
         if (icopt.eq.9) then
            write (*,1080)
         else 
            write (*,1090)
         end if 

         read (*,'(a)') y

         if (y.eq.'y'.or.y.eq.'Y') then

            write (*,1100) (vname(iv(i)),i=1,ivct)
            write (*,'(/)')
            write (*,2010) 'coordinate','coor.dat'
            read (*,'(a)') cfname
            if (cfname.eq.' ') cfname = 'coor.dat'

            open (n8,file=cfname,iostat=ier,status='old')
            if (ier.ne.0) write (*,1140) cfname

            close (n8)

            if (icopt.eq.9) then
               icopt = 11
            else
               icopt = 10
            end if 

            oned = .false.
            icont = 1

         end if

      else

         oned = .false.

      end if

c                                 icopt now has the following values

c                                 1 - convex
c                                 2 - 2d grid
c                                 3 - 1d grid - no fileio (oned true)
c                                 4 - swash
c                                 5 - 1d frac - no fileio (oned true)
c                                 6 - 0d frac
c                                 9 - frac2d - no fileo
c                                10 - 1d - with fileio (oned false)
c                                11 - frac2d - with fileio (oned false)

      if (icopt.gt.2.and.icopt.lt.10.and.icopt.ne.6.and.icopt.ne.9) 
     *                                               icopt = icopt - 1

c                                 icopt now has the following values

c                                 1 - convex
c                                 2 - 2d grid
c                                 2 - 1d grid - no fileio (oned true)
c                                 3 - swash
c                                 4 - 1d frac - no fileio (oned true)
c                                 6 - 0d frac
c                                 9 - frac2d - no fileo
c                                10 - 1d - with fileio (oned false)
c                                11 - frac2d - with fileio (oned false)

      if (icopt.ne.9) then 
c                                 ====================================
c                                 ask if p = f(T) or vice versa for all
c                                 phase diagram calculations:
         do

            idep = 0

            if (icopt.ne.3.and.icopt.ne.10.and.icopt.ne.6) then 

               write (*,1050) vname(1),vname(2)
               if (icopt.eq.4) write (*,1160) 

               read (*,'(a)') y

               if (y.eq.'Y'.or.y.eq.'y') then

                  do 
                     write (*,1060) vname(1),vname(2),vname(2),vname(1)
                     read (*,*,iostat=ier) idep
                     if (ier.eq.0) exit
                     call rerr
                  end do
c                                 cycle on bad choice
                  if (idep.lt.1.or.idep.gt.2) cycle
c                                 reset the variable counters and flags:
c                                 depend changes the dependent variable
c                                 pointer to iv(ivct), iind points to the
c                                 position of the independent variable in 
c                                 the array v, idep to the dependent v.  
                  call depend (ivct,idep,iind,iord,c,dtext)
c                                 kvct is the number of independent potentials
                  kvct = ivct - 1

                  if (icp.eq.1.and.kvct.eq.1) oned = .true.

               end if

            end if 

            exit 

         end do

      else 

         idep = 0

      end if 

      if (oned) then 

         icont = 1
c                                 ======================================
c                                 for 1d calculations get the independent
c                                 variable, currently composition is not
c                                 allowed. 
         call getxvr (ivct,jvct,icont,jc,oned,'the independent')
c                                 get sectioning variables values:
         do j = 1, ivct
            if (iv(j).eq.jc) cycle 
            call redvar (j,2)
            vmax(iv(j)) = vmin(iv(j))
         end do     
c                                  set the icopt flag to its final value
         if (icopt.eq.2) then
c                                  1-d minimization
            icopt = 5

         else 
c                                  fractionation, also write the grid blurb
            icopt = 7
            write (*,3100) grid(4,1),grid(4,2)

         end if

      else if (icopt.eq.10) then
c                                  1-d "gridded min" from a file, only allow p-t
         icont = 1

         call getxvr (ivct,jvct,icont,jc,oned,
     *                                'the independent path')

      else if (icopt.eq.9.or.icopt.eq.11) then
c                                  2-d fractionation, only allow p and t
         ivct = 2
         icont = 0 
c                                  choose the primary variable (IV(1)):
         call getxvr (ivct,jvct,icont,jc,oned,
     *                'the primary (usually pressure)')

      else if (icopt.eq.6) then 
c                                  0-d infiltration/reaction/fractionation
         icopt = 12

         icont = 1

         oned = .true.

         write (*,3120)

         do i = 1, 2
            call redvar (i,3)
            vmax(iv(i)) = vmin(iv(i))
         end do

      else if (icopt.eq.3) then 
c                                  ==================================-
c                                  swash, not a lot to do
         icopt = 4   


      else if (icopt.eq.2) then
c                                  =====================
c                                  gridded minimization:
         icopt = 5
         jvct = 0
         icont = 1
c                                  Select the x variable
         call getxvr (ivct,jvct,icont,jc,oned,'x-axis')


         if (ivct.eq.2.and.icont.eq.1) then 
c                                 there is no C variable and there 
c                                 are only 2 potentials, 
c                                 the y variable must be iv(2)     
            call redvar (2,1)
            jvct = ivct   

         else 
c                                 select the y variable 
            if (ivct.gt.1.or.icont.eq.2.and.icp.gt.2) then
               jst = 2
               if (icont.eq.2) jst = 1

               do 

                  write (*,2130)

                  do 
                     write (*,2140) (j,vname(iv(j)), j = jst, ivct)
                     if (icp.gt.2.and.icont.eq.2) write (*,1480) j
                     write (*,*) ' '
                     read (*,*,iostat=ier) jc
                     if (ier.eq.0) exit
                     call rerr
                  end do
 
                  if (jc.gt.ivct+1.or.jc.lt.jst) then
                     write (*,1150)
                     cycle 
                  else if (jc.eq.ivct+1) then
                     icont = 3
                  end if

                  exit 

               end do 

            else if (icont.eq.2) then
 
               jc = 1 

            else 

               jc = 2
 
            end if

            if (icont.lt.3) then 
               ind = 2
               if (icont.eq.2) ind = 1
               ix = iv(ind)
               iv(ind) = iv(jc)
               iv(jc) = ix
               jvct = jvct + 1
               call redvar (ind,1)
            end if  
         end if
c                                 get sectioning variables values:
         do j = jvct+1, ivct
            call redvar (j,2) 
            vmax(iv(j)) = vmin(iv(j))
         end do 

         write (*,3070)
c                                  inform the user of the grid settings:
         do j = 1, 2
            if (j.eq.1) then 
               stext = 'exploratory'
            else
               stext = 'auto-refine'
            end if 

            loopy = (grid(2,j)-1) * 2**(grid(3,j)-1) + 1
            loopx = (grid(1,j)-1) * 2**(grid(3,j)-1) + 1

            write (*,3080) stext,grid(3,j),grid(1,j),grid(2,j),loopx,
     *                     loopy

         end do 

         write (*,3090) opname

      else if (icopt.eq.1) then
c                                  =========================
c                                  Normal computational mode
         write (*,1500)
         if (ivct.gt.1) write (*,1590)
c                                  get choice
         call rdnumb (c(0),0d0,icopt,0,.false.)

         if (icopt.lt.1.or.icopt.gt.2) icopt = 0
c                                  convert to internal values
         if (icopt.eq.2) then 
            icopt = 1
         else if (icopt.eq.1) then 
            icopt = 3
         end if 

         if (icopt.eq.1) then
c                                  Select the x variable (IV(1)):
            call getxvr (ivct,jvct,icont,jc,oned,'x-axis')
c                                  select the y variable (iv(2)):
            if (ivct.gt.2) then
 
               do 

                  write (*,2130)

                  do 
                     write (*,2140) (j,vname(iv(j)), j = 2, ivct)
                     read (*,*,iostat=ier) jc
                     if (ier.eq.0) exit
                     call rerr
                  end do

 
                  if (jc.gt.ivct.or.jc.lt.2) then
                     write (*,1150)
                     cycle
                  end if

                  exit
         
               end do 
 
            else
 
               jc=2
 
            end if
 
            ix = iv(2)
            iv(2) = iv(jc)
            iv(jc) = ix
 
            call redvar (2,1)
c                                  specify sectioning variables (iv(3)):
            do j = 3, ivct
               call redvar (j,2) 
               vmax(iv(j)) = vmin(iv(j))
            end do 

         else if (icopt.eq.3) then
c                                  select the y variable (iv(1)):
            call getxvr (ivct,jvct,icont,jc,oned,'y-axis')
c                                  specify sectioning variable (iv(2)):
            do j = 2, ivct
               call redvar (j,2) 
               vmax(iv(j)) = vmin(iv(j))
            end do 
         end if

      else if (icopt.eq.11) then 
         icont = 1
      end if

c                                 icopt now has the following final values for output

c                                 0 - convex composition
c                                 1 - convex spd
c                                 3 - convex mixed
c                                 5 - 2d grid
c                                 5 - 1d grid - no fileio (oned true)
c                                 4 - swash
c                                 7 - 1d frac - no fileio (oned true)
c                                12 - 0d frac (oned true)
c                                 9 - frac2d - no fileo
c                                10 - 1d - with fileio (oned false)
c                                11 - frac2d - with fileio (oned false)

      icth = icp + isat
      
      if (icopt.gt.4) then 
c                                 =========================
c                                 Compositional constraints, only
c                                 for constrained minimization
         if (isat.ne.0) then 
            write (*,1450)
            read (*,'(a)') y 
            if (y.ne.'y'.and.y.ne.'Y') then 
               jcth = icp
            else
               write (*,1460)
            end if 
         else 
            jcth = icp 
         end if 

         if (jcth.ne.icp)  then 

            jcth = icp

            do i = icp + 1, icth
               write (*,1430) pname(i)
               read (*,'(a)') y
               if (y.eq.'Y'.or.y.eq.'y') then
                  jcth = jcth + 1
               else 
                  exit
               end if 
            end do 

         end if 

         if (icopt.ne.9.and.icopt.ne.11) then 
c                                ask if weight amounts, otherwise
c                                use molar amounts.
            if (icopt.ne.12) then
               write (*,1420) 
               read (*,'(a)') y
               if (y.eq.'Y'.or.y.eq.'y') then 
                  iwt = 1
                  amount = 'mass '
               end if
            end if 

            write (*,1380)
c                                get the bulk composition:
            if (icopt.eq.12) then
c                                0-d infiltration special case
               do i = 1, 2

                  do 
c                                initial bulk, then infiltrant:
                     if (i.eq.1) then 
                        write (*,1370) amount
                     else 
                        write (*,1360) amount
                     end if 

                     write (*,'(12(1x,a))') (pname(j), j = 1, jcth)
                     write (*,1410) 
                     read (*,*,iostat=ier) (dblk(i,j), j = 1, jcth)
                     if (ier.eq.0) exit
                     call rerr

                  end do

               end do 

            else if (icont.eq.1) then 
c                                fixed bulk composition
               do 
                  write (*,1390) amount
                  write (*,'(12(1x,a))') (pname(j), j = 1, jcth)
                  write (*,1410) 
                  read (*,*,iostat=ier) (dblk(1,j), j = 1, jcth)
                  if (ier.eq.0) exit
                  call rerr
               end do

            else 
c                                user must define compositional variables
               if (lopt(1)) then
c                                closed c-space
                  if (icont.eq.2) then
                     write (*,1520) '   C = C0*(1 - X_C1) + C1*X_C1'
                     write (*,1510) '   C = C0 + C1*X_C1'
                  else 
                     write (*,1540) 
     *            '   C = C0*(1 - X_C1 - X_C2) + C1*X_C1 + C2*X_C2'
                     write (*,1510) '   C = C0 + C1*X_C1 + C2*X_C2'
                  end if
 
               else 
c                                 open c-space
                  if (icont.eq.2) then
                     write (*,1520) '   C = C0 + C1*X_C1'
                     write (*,1510) '   C = C0*(1 - X_C1) + C1*X_C1'
                  else 
                     write (*,1540) '   C = C0 + C1*X_C1 + C2*X_C2'
                     write (*,1510) 
     *            '   C = C0*(1 - X_C1 - X_C2) + C1*X_C1 + C2*X_C2'
                  end if

               end if 
 
               do i = 1, icont

                  do 
                     write (*,1390) amount
                     write (*,'(12(1x,a))') (pname(j), j = 1, jcth)
                     write (*,1530) nc(i)
                     read (*,*,iostat=ier) (dblk(i,j), j = 1, jcth)
                     if (ier.eq.0) exit
                     call rerr
                  end do

               end do
 
            end if

         end if

      end if 

1050  format (/,'The data base has ',a,' and ',a,' as default ',
     *       'independent potentials.',/,'Make one dependent on the ',
     *       'other, e.g., as along a geothermal gradient (y/n)? ')
1060  format (/,'Select dependent variable:',//,'  1 - ',a,' = f(',a,')'
     *       ,/,'  2 - ',a,' = f(',a,')',/)
1080  format (/,'Enter P(Z0,DZ), T(Z0,DZ) from a file (Y/N)?')
1090  format (/,'Enter path coordinates from a file (Y/N)?')
1100  format (/,'In this mode VERTEX/WERAMI read path coordinates',
     *        'from a file',/,'the file must be formatted so that',
     *        ' the coordinates of each point',/,'are on a separate',
     *        ' line, the coordinates are, in order:',4x,5(a,2x))
1140  format (/,'File: ',a,/,'Does not exist, you must create it',
     *        ' before running VERTEX.',/)
1150  format (/,'hunh?',/)
1160  format (/,'Answer yes to specify a P-T path for phase ',
     *          'fractionation calculations.',/)
1210  format ('Select the path variable for the calculation:',/)
1360  format (/,'Enter the ',a,' amounts of the components per mole',
     *        ' of infiltrating fluid',/,'these amounts define a mole',
     *        ' of the fluid:')
1370  format (/,'Enter the ',a,' amounts of the components ',
     *        'present initially in the system:')
1380  format (/,'The amounts you enter next need not be normalized;',
     *        ' regardless of units, they',/,
     *        'define the molar amount of the system',/)
1390  format (/,'Enter the ',a,' amounts of the components:')
1410  format ('for the bulk composition of interest:')
1420  format (/,'Specify component amounts by mass (Y/N)?')
1430  format ('Constrain component ',a,' (Y/N)?')
1450  format (/,'All thermodynamic components must be ',
     *         'constrained,',/'specify ',
     *         'saturated components also (Y/N)?')
1460  format (/,'The next prompts are for the saturated',
     *        ' component constraints.',/,
     *        'Answering no at any point',
     *        ' completes the set of constraints.',/)
1470  format (5x,i1,' - Composition X_C1* (user defined)')
1480  format (5x,i1,' - Composition X_C2 (user defined)')
1500  format (//,'Specify number of independent potential variables:',
     *         /,5x,'0 - Composition diagram [default]',/,
     *           5x,'1 - Mixed-variable diagram')
1510  format (/,'To compute bulk compositions as:',a,/,'change the ',
     *       'computational option keyword closed_c_space.')
1520  format (/,'The bulk composition of the system will be computed ',
     *       'as:',/,a,/,
     *       'where X_C1 varies between 0 and 1, and C0 and C1 are ',
     *       'the compositions',/,'specified next.')
1530  format ('to define the composition ',a)
1540  format (/,'The bulk composition of the system will be computed ',
     *       'as:',/,a,/,
     *       'where X_C1 and X_C2 vary between 0 and 1, and C0, C1 ',
     *       'and C2 are the',/,'compositions specified next.')
1570  format (/,'*X_C1 can not be selected as the y-axis variable',/)
1590  format (5x,'2 - Sections and Schreinemakers-type diagrams')
2010  format ('Enter the ',a,' file name [default = ',a,']:')
2111  format (/,'Select x-axis variable:')
2130  format (/,'Select y-axis variable:')
2140  format (5x,I1,' - ',a)
2210  format (/,'Select vertical axis variable:')
3070  format (/,'For gridded minimization, grid resolution is ',
     *        'determined by the number of levels',/,
     *        '(grid_levels) and the resolution at the lowest level ',
     *        'in the X- and Y-directions (x_nodes',/,
     *        'and y_nodes) these parameters are currently set for ',
     *        'the exploratory and autorefine cycles',/,
     *        'as follows:',//,'stage        grid_levels  xnodes  ',
     *        'ynodes    effective resolution')
3080  format (a,7x,i1,8x,i4,4x,i4,6x,i4,' x',i4,' nodes')
3090  format (/,'To change these options edit or create',
     *        ' the file ',a,/,'See: ',
     *        'www.perplex.ethz.ch/perplex_options.html#grid_parameters'
     *        ,/)
3100  format (/,'For phase fractionation calculations the number of ',
     *        'points computed along the path',/,
     *        'is determined by the 1d_path parameter. The values ',
     *        'for this parameter are currently',/,
     *        'set to ',i3,' and ',i3,' points for the exploratory ',
     *        'and autorefine cycles.')
3120  format (/,'For zero-dimensional infiltration calculations:',/,
     *        3x,' - the fluid must be described by a hybrid EoS ',
     *           'solution model',/,
     *        3x,' - the fluid components must be specified as ',
     *           'thermodynamic components',/,
     *          'Additionally, to account for solute ',
     *           'chemistry:',/,
     *         3x,' - select a thermodynamic data file that includes ',
     *            'electrolyte data',/,
     *         3x,' - set aq_lagged_speciation and refine_endmembers to'
     *           ,' true (T)')
7150  format (/,'*Although only one component is specified for the ',a,
     *       ' phase, its equation of state ',/,
     *       'permits use of its compositional variable: ',a,'.',/) 

      end 


      subroutine getxvr (ivct,jvct,icont,jc,oned,text)
c----------------------------------------------------------------------
c read primary variable 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character text*(*)

      logical oned

      integer j,ivct, ier, ix, jc, icont, jvct, jcont

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*8 vname, xname
      common/ csta2 /xname(k5),vname(l2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------
c                                 jcont = 0, no bulk composition allowed
      jcont = 0
c                                 at this point:
c                                 1 - convex
c                                 2 - 2d grid (oned false)
c                                 2 - 1d grid - no fileio (oned true)
c                                 3 - swash
c                                 4 - 1d frac - no fileio (oned true)
c                                 6 - 0d frac
c                                 9 - frac2d - no fileo
c                                10 - 1d - with fileio (oned false)
c                                11 - frac2d - with fileio (oned false)
      if (icopt.eq.1) then 
c                                 2d projection, no bulk composition allowed
      else if (icopt.eq.3) then 
c                                 1d projection, no bulk composition allowed
      else if (icopt.eq.9.or.icopt.eq.11) then 
c                                 2-d fractionation
      else if ((icopt.eq.2.and.oned).or.icopt.eq.10) then
c                                 1d-gridded min file I/O
      else if (icopt.eq.4.or.icopt.eq.2.or.icopt.eq.5) then
c                                 icopt = 4 -> 1d with P(T) or T(P)
c                                 1d- & 2d-gridded min
c                                 bulk composition allowed
         if (icp.gt.1) jcont = 1
      else 
         call errdbg ('unanticipated icopt value in getxvr')
      end if

      do 

         write (*,2111) text

         do

            write (*,2140) (j,vname(iv(j)), j = 1, ivct)
c                                 bulk composition variable
            if (jcont.eq.1) write (*,1470) j
c                                 saturated phase composition explanation
            if (ifct.eq.1) write (*,7150) vname(3)
c                                 warn that the bulk composition
c                                 variable can't be used on y
            if (.not.oned.and.jcont.eq.1) write (*,1570)

            read (*,*,iostat=ier) jc
            if (ier.eq.0) exit
            call rerr

         end do

         if (jcont.eq.1.and.jc.gt.ivct+1.or.
     *       jcont.eq.0.and.jc.gt.ivct.or.
     *                                         jc.lt.1) then
            write (*,'(/,''hunh?'',/)')
            cycle

         else if (jc.eq.ivct+1) then 
c                                 bulk composition variable, limits
c                                 and position are automatic
            icont = 2

         else 
c                                 normal variable, swap positions
            jvct = 1
            ix = iv(1)
            iv(1) = iv(jc)
            iv(jc) = ix
c                                 read limits
            call redvar (1,1)

         end if

         exit 

      end do

1470  format (5x,i1,' - Composition X_C1* (user defined)')
1570  format (/,'*X_C1 can not be selected as the y-axis variable',/)
2111  format (/,'Select ',a,' variable:')
2140  format (5x,I1,' - ',a)
7150  format (/,'*Although only one component is specified for the ',a,
     *       ' phase, its equation of state',/,
     *       ' permits use of its compositional variable: ',a,'.',/)

      end