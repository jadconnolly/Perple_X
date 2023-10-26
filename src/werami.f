
      program werami
c----------------------------------------------------------------------
c       lun restrictions
c       n8 - echo 
c       n7 - geotherm
c       n5 - out
c       n6 - special out
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, fake, err

      integer imode, ierr, i

      character*100 n5name,n6name

      logical oned
      common/ cst82 /oned

      integer ivar,ind
      common/ cst83 /ivar,ind

      integer iam
      common/ cst4 /iam

      logical fileio, flsh, anneal, verbos, siphon, colcmp, usecmp
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,verbos,siphon,
     *                usecmp, colcmp

      integer icps, jcx, jcx1, kds
      logical stol, savg, spec
      double precision rcps, a0
      common/ comps /rcps(2*k5,m13),a0(m13,2),icps(2*k5,m13),jcx(m13),
     *               jcx1(m13),kds(m13),stol(i11),savg(i11),spec(m13)
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 3
c                                 version info
      call vrsion (6)
c                                 initialize some flags
      first = .true.
      outprt = .false.
      fake   = .false.
      rxn = .false.
c                                 this could be eliminated by passing first 
c                                 to chsprp.
      do i = 1, h9
         stol(i) = .false.
      end do 
c                                 initialize ind only used for outmod, but
c                                 only set for 1d calcs in getind
      ind = 1
c                                 read input from unit n1 (terminal/disk).
c                                 input1 also initializes:
c                                 equilibrium counters; units n2 n4 and n6;
c                                 and the limits for numerical results.
      call input1 (first,err)
c                                 set ivar flag, this indicates the number
c                                 of possible independent plotting variables, jvar
c                                 indicates the number of thermodynamic variables
      ivar = 2
      if (icopt.eq.7.and.fileio) ivar = 3
c                                 don't allow users to do anything
c                                 other than gridded min
      if (icopt.lt.5.and.icopt.ne.2) call error (4,1d0,icopt,'PSVDRAW')
c                                 titration/0-d fractionation
      if (icopt.eq.12) oned = .true.
c                                 read thermodynamic data on unit n2:
      call input2 (fake)
c                                 read autorefine lists
      call setau1 
c                                 read data for solution phases on n9:
      call input9 (fake)
c                                 seismic data summary file
      if (lopt(50)) call outsei

      call setau2 
c                                 read the plt/blk files
      call interm (.false.,err)
c                                 organize variables 
      call getvar
c                                 initialize variables
      call setvar

      do 

         write (*,1000) 
         if (.not.oned) write (*,1010)
         write (*,1020)
         if (.not.oned) write (*,1025)
         write (*,1026)

         read (*,*,iostat=ierr) imode
         if (ierr.ne.0) cycle 

         if (first.and.imode.eq.1) then 
c                                 make console output echo to rpl file
            call fopenn (n8,0,n5name,n6name)
            first = .false.

         end if 

         if (imode.eq.1) then          
c                                 mode 1 - print assemblage data at a 
c                                 specific p-t condition
            call mode1 

         else if (imode.eq.2) then 

            if (oned) then 
               write (*,1030) 
               cycle
            end if 

            call mode2 

         else if (imode.eq.3) then 
c                                extract the data
            if (oned) then 
               call mode31 
            else 
               call mode3 
            end if 

         else if (imode.eq.4) then 

            call mode4 

         else if (imode.eq.0) then 

            exit 

         end if 

      end do 
c                                 close "echo" file
      if (imode.eq.1.or.imode.eq.3) close (n8)

1000  format ('Select operational mode:',/,
     *        4x,'1 - properties at specified conditions')
1010  format (4x,'2 - properties on a 2d grid')
1020  format (4x,'3 - properties along a 1d path')        
1025  format (4x,'4 - as in 3, but input from file')
1026  format (4x,'0 - EXIT')
1030  format (/,'Invalid choice for 1d grids',/)

      end 

      subroutine mode2 
c----------------------------------------------------------------------
c sample data on an x-y grid 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node, change, readyn

      integer i, j, nxy(2), dim

      double precision tmin(2), tmax(2), dx(2)

      character n6name*100, n5name*100

      external readyn

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      integer ivar,ind
      common/ cst83 /ivar,ind

      character vnm*8
      common/ cxt18a /vnm(l3)
c----------------------------------------------------------------------
      node = .false. 
      dim = 2
c                                 select the property
      call chsprp
c                                 set up coordinates etc
c                                 allow restricted plot limits
      write(*,1040)

      if (readyn()) then 

         change = .false.

         do i = 1, 2
30          write (*,1060) vnm(i),vmn(i),vmx(i)
            read (*,*,err=30) tmin(i),tmax(i)
            if (tmin(i).ne.vmn(i).or.tmax(i).ne.vmx(i)) change = .true.
         end do

         if (lopt(48)) then 
            write (*,'(/,a,/)') '**warning ver084** sample_on_grid is '
     *                        //'disabled for arbitrary limits'
            lopt(48) = .false.
         end if

      else

         tmin(1) = vmn(1)
         tmin(2) = vmn(2)
         tmax(1) = vmx(1)
         tmax(2) = vmx(2)

      end if

      if (lopt(48).and.(icopt.eq.5.or.icopt.eq.2)) then

         if (lopt(47)) then
c                                 unsplt, could work out sample on grid with 
c                                 the spt file, but i'm too lazy.
           nxy(1) = loopx
           nxy(2) = loopy

           write (*,1000) loopx, loopy

         else 

            call getlvl (i)

            nxy(1) = (loopx - 1)/ 2**(jlev-i) + 1
            nxy(2) = (loopy - 1)/ 2**(jlev-i) + 1

         end if 

      else if (lopt(48)) then 
c                                 frac2d calculations
         nxy(1) = loopx
         nxy(2) = loopy

      else
c                                 arbitrary number of grid points
         write (*,1080) vnm(1),vnm(2)
         read (*,*) nxy

      end if

      do i = 1, 2
         dx(i) = (tmax(i)-tmin(i))/dfloat(nxy(i)-1)
      end do 

      call tabhed (n5,tmin,dx,nxy,dim,n5name,n6name)

      do j = 1, nxy(2)

         var(2) = tmin(2) + dx(2)*dfloat(j-1)
c                                 round off tests:
         if (var(2).gt.tmax(2)) then
            var(2) = tmax(2)
         else if (var(2).lt.tmin(2)) then 
            var(2) = tmin(2)
         end if 

         do i = 1, nxy(1) 

            var(1) = tmin(1) + dx(1)*dfloat(i-1)
c                                 round off tests:
            if (var(1).gt.tmax(1)) then
               var(1) = tmax(1)
            else if (var(1).lt.tmin(1)) then 
               var(1) = tmin(1)
            end if

            call polprp (dim)
 
         end do
 
      end do 
c                                 wrap up the calculation
      call finprp (dim,n5name,n6name,node) 

1000  format (/,'**warning ver638** for output files generated by UNS',
     *       'PLT, sample_on_grid uses the',/,'highest resolution pos',
     *       'sible (',i4,'x',i4,'), if this is excessive set ',/,
     *       'sample_on_grid to false and restart WERAMI',/)
1040  format (/,'Change default variable range (y/n)?')
1060  format (/,'Current limits on ',a,' are: ',g14.7,'->',g14.7,/,
     *          'Enter new values:')
1080  format (/,'Enter number of nodes in the ',a,' and ',a,
     *          ' directions:')

      end 

      subroutine mode1 
c----------------------------------------------------------------------
c read x-y coordinates for a 2-d section from the console. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical quit, nodata

      integer itri(4),jtri(4),ijpt

      double precision wt(3)
c----------------------------------------------------------------------
      do 

         call readxy (quit)

         if (quit) exit

         call triang (itri,jtri,ijpt,wt)

         if (ijpt.eq.0) then 
            nodata = .true.
         else 
            call getloc (itri,jtri,ijpt,wt,nodata)
         end if 

         if (nodata) then 
            write (*,1000) 
         else 
            call calpr0 (6)
            call calpr0 (n8)
         end if

      end do 

1000  format (/,'No data at this condition, presumably because',
     *          ' minimization failed.',/)

      end 

      subroutine setval 
c--------------------------------------------------------------------
c setval sets the values of the thermodynamic variables after a call
c to readxy gets the section coordinates in array var, five cases:

c icopt = 10 -> var(1) is the index of coordinate array
c icopt = 9  -> frac2d calculations x-y indirectly related to p-t
c icont = 1  -> independent variables are the 1st and 2nd potentials
c icont = 2  -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3  -> independent variables are compositional variables
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ind, j

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      logical fileio, flsh, anneal, verbos, siphon, colcmp, usecmp
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,verbos,siphon,
     *                usecmp, colcmp
c----------------------------------------------------------------------
      if (icopt.eq.7.and.fileio) then
c                                 fractionation from a coodinate file
c                                 only exact nodal coordinates allowed
          ind = idint(var(1))

          do j = 1, ipot
             v(jv(j)) = vip(j,ind)
          end do 

          var(2) = v(jv(1))
          var(3) = v(jv(2))

      else if (icopt.eq.12) then 
c                                 0-d infiltration, var(2) is the
c                                 number of aliquots added. var(1)
c                                 is the cumulative moles of aliquot.
c         var(1) = (var(2)-1d0) * nopt(36)
          var(2) = var(1)/nopt(36) + 1d0

      else if (icopt.eq.9) then 
c                                 change sign on dz because of downward
c                                 directed depth coordinate.
         call fr2dpt (var(1),var(2))

         var(3) = v(1)
         var(4) = v(2)

      else if (icont.eq.1) then 

         v(iv1) = var(1)
         v(iv2) = var(2)
         call incdp0
         if (idep.ne.0) var(jvar) = v(idep) 

      else if (icont.eq.2) then 

         cx(1) =  var(1)
         call setblk 

         v(iv1) = var(2)
         call incdep (iv1)
         if (idep.ne.0) var(jvar) = v(idep) 

      else 

         cx(1) = var(1)
         cx(2) = var(2)
         call setblk

      end if 

      end

      subroutine readxy (quit)
c----------------------------------------------------------------------
c read x-y coordinates for a 2-d section from the users console, these
c are then assigned to the thermodynamic variables by calling setval.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ier
      logical quit

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3) 

      logical oned
      common/ cst82 /oned
c----------------------------------------------------------------------
      if (.not.oned) then

         do 

            quit = .false.

            write (*,1000) vnm(1),vnm(2)
            read (*,*,iostat=ier) var(1),var(2)
            if (ier.ne.0) cycle 

            if (var(1)+var(2).eq.198d0) quit = .true.

            if (.not.quit) then 
                
               quit = .false.

               do i = 1, 2
                  if (vmn(i).lt.vmx(i)) then
                     if (var(i).lt.vmn(i).or.var(i).gt.vmx(i)) then  
                        write (*,1010) vnm(i),vmn(i),vmx(i)
                        quit = .true.
                     end if 
                  else 
                     if (var(i).lt.vmx(i).or.var(i).gt.vmn(i)) then  
                        write (*,1010) vnm(i),vmn(i),vmx(i)
                        quit = .true.  
                     end if 
                  end if 
               end do 

               if (quit) cycle

            end if 

            exit
 
         end do 

      else 

         do 

            quit = .false.

            write (*,1020) vnm(1)
            read (*,*) var(1)


            if (var(1).eq.999d0) quit = .true.

            if (.not.quit) then 
               if (vmn(1).lt.vmx(1)) then
                  if (var(1).lt.vmn(1).or.var(1).gt.vmx(1)) then  
                     write (*,1010) vnm(1),vmn(1),vmx(1)
                     cycle  
                  end if 
               else 
                  if (var(1).lt.vmx(1).or.var(1).gt.vmn(1)) then  
                     write (*,1010) vnm(1),vmn(1),vmx(1)
                     cycle  
                  end if 
               end if 
            end if
            
            exit
 
         end do 

      end if 

      if (.not.quit) call setval

1000  format (/,'Enter ',a,' and ',a,' (99 and 99 to quit):')
1010  format (/,'The plot file range for ',a,' is ',g12.4,' - ',g12.4,
     *        /,'Try again:',/)
1020  format (/,'Enter ',a,' (999 to quit):')

      end 

      subroutine polprp (dim)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical nodata, ongrid

      integer itri(4), jtri(4), ijpt, lop, icx, komp, i, j, dim

      double precision wt(3), mode

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname
c----------------------------------------------------------------------
c                                 set variables to x-y value
      call setval

      do i = 1, iprop

         lop = kop(i)
         icx = kcx(i)
         komp = k2c(i)
         lflu = kfl(i)
c                                 initialize, necessary?
         prop(i) = nopt(7)
c                                 get node(s) to extract value
         call triang (itri,jtri,ijpt,wt)

         if (ijpt.eq.0) then 
c                                 missing data at the node
            call badnum (dim)

            return 

         else 
c                                 compute all properties
            call getloc (itri,jtri,ijpt,wt,nodata)

            if (nodata) then 

               call badnum (dim)

               return  

            else 
c                                 get the properties of interest
               if (lop.eq.25) then 
c                                 get all modes
                  call allmod

                  exit 

               else if (lop.eq.36.or.lop.eq.38) then 
c                                 multiple property lists, PHEMGP
                  call allprp (dim)

                  exit 

               else if (lop.eq.24) then 
c                                 assemblage index request, lots
c                                 of redundant calc, but should be 
c                                 moved to getptp.
c                                 no need to call triang/getlow
                  call xy2ij (itri(1),jtri(1),nodata,ongrid)

                  prop(i) = iap(igrd(itri(1),jtri(1)))
           
               else if (lop.eq.40) then 
c                                 solute chemistry, call getprp (mode computation) 
c                                 to check if the solvent is stable, and to decide 
c                                 which phase (komp) is to be used in the case of 
c                                 immiscibility.
                  call getprp (mode,7,icx,komp,.false.) 
c                                 decide whether to output lagged or back-calculated
c                                 speciation into props:
                  if (komp.eq.0) then
c                                 phase not stable

                     do j = 1, iprop
                        prop(j) = nopt(7)
                     end do

                  else if (lopt(32).and..not.kfl(1)) then 
c                                  do the lagged speciation calculation and output
c                                  the result to prop.
                     call lagprp (komp) 

                  else
c                                  output back-calculated result to props, -1 signals
c                                  tab file output.  
                     call aqrxdo (komp,-1)

                  end if 

                  exit 

               else 
c                                 general properties
                  call getprp (prop(i),lop,icx,komp,.false.)

               end if 

            end if 

         end if 

      end do 

      if (lop.ne.36.and.lop.ne.38) call outprp (dim)

      end 

      subroutine badnum (dim)  
c----------------------------------------------------------------
c badnum: writes missing data message; assigns and outputs badnumber
c data.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer dim, i 

      character vnm*8
      common/ cxt18a /vnm(l3)   

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c----------------------------------------------------------------------
      write (*,1000) vnm(1),var(1),vnm(2),var(2),nopt(7)
c                                 for phemgp format:
      ntot = 0
      tname = 'Missing data  '
c                                 in general
      do i = 1, iprop
         prop(i) = nopt(7)
      end do 

      call outprp (dim)

1000  format ('Missing data at ',a,'=',g12.5,', ',a,'=',g12.5,
     *        ' assigned ',g12.5,' to all properties')

      end

      logical function bdnum (num) 
c----------------------------------------------------------------------
c test if num eq bad_number with check to avoid NaN fp error
c----------------------------------------------------------------------
      implicit none

      double precision num

      include 'perplex_parameters.h'
c----------------------------------------------------------------------
      bdnum = .false.

      if (isnan(nopt(7))) then 
         if (isnan(num)) bdnum = .true.
      else
         if (num.eq.nopt(7)) bdnum = .true.
      end if

      end 

      subroutine outprp (dim)
c----------------------------------------------------------------------
c outprp outputs properties computed by chsprp for dim-dimensional tables: 
c   dim    = 1 tables include independent variables;
c   dim    = 2, kcx(1) ne 999, tables include independent variables if lopt(15) = T;
c   kcx(1) = 999, write phemgp format that includes a name, counter, and 
c            independent variables. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer dim, i

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer ivar,ind
      common/ cst83 /ivar,ind

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c----------------------------------------------------------------------
      if (kcx(1).eq.999) then 
c                                 write phemgp format
         write (n5,'(a14,1x,7x,i2,6x,200(g16.6e3,1x))') tname,ntot,
     *                                               (var(i),i=1,ivar), 
     *                                               (prop(i),i=1,iprop)

      else if (lopt(15).or.dim.eq.1) then 
c                                 write spreadsheet tab format
c         prop(5) = prop(6)+prop(8)*2d0+ prop(9)
c         prop(5) = prop(12)+ prop(13)

         write (n5,'(200(g14.6e3,1x))') (var(i),i=1,ivar), 
     *                                  (prop(i),i=1,iprop)

      else 
c                                 write compact tab format
         write (n5,'(200(g14.6e3,1x))') (prop(i),i=1,iprop)
      end if 
c                                 check property ranges       
      do i = 1, iprop
c                                
         if (isnan(prop(i))) cycle
c                                 first check eliminates logical 
c                                 comparisons with NaN's that compaq
c                                 fortran doesn't like.
         if (.not.isnan(nopt(7))) then 
            if (prop(i).eq.nopt(7)) cycle 
         end if 
         
         if (prop(i).gt.prmx(i)) prmx(i) = prop(i)
         if (prop(i).lt.prmn(i)) prmn(i) = prop(i)
 
      end do

      end 

      subroutine getprp (prop,lop,icx,komp,aprp)
c----------------------------------------------------------------
c getprp gets properties:

c   aprp - true if called by allprp, false otherwise

c   jd   - is the pointer to the assemblage
c   jflu - if 1, bulk properties include fluid phase; else 0

c   lop  - flag indicating the property chosen
c   icx  - if lop = 6, the component chosen
c   icx  - if the identity of the phase chosen
c          icx = 0 if bulk properties requested.
c   komp - if lop = 7, output pointer to the phase in assemblage
c          if lop = 8, input pointer to component

c   icps - if lop = 8, the indices of the components. 
c   
c requestable properties (indicated by lop)

c 1                 Specific enthalpy (J/m3)',
c 2                 Density (kg/m3)',
c 3                 Specific Heat capacity (J/K/m3)',
c 4                 Expansivity (1/K, for volume)',
c 5                 Compressibility (1/bar, for volume)
c 6                 Weight percent of a component
c 7                 Mode (Vol %) of a compound or solution',
c 8                 Composition of a solution'
c 9                 Grueneisen thermal ratio',
c 10                Adiabatic bulk modulus (bar)',
c 11                Sound velocity (km/s)
c 12                Shear modulus (bar)',
c 13                P-wave velocity (km/s)',
c 14                S-wave velocity (km/s)',
c 15                Vp/Vs
c 16                Specific Entropy (J/K/m3)'
c 17                Entropy (J/K/kg)'
c 18                Enthalpy (J/kg)'
c 19                Heat Capacity (J/K/kg)'
c 20                Specific mass of phase (kg/m3)
c 21                Poisson's Ratio 
c 22                Molar Volume (J/bar) 
c 23                Chemical potentials (J/mol)
c 24                not passed to getprp
c 25                output all modes
c 26                Sound velocity temperature derivative (km/s/K)
c 27                P-wave velocity temperature derivative (km/s/K)
c 28                S-wave velocity temperature derivative (km/s/K)
c 29                Adiabatic bulk modulus temperature derivative (bar/K)
c 30                Shear modulus temperature derivative (bar/K)
c 31                Sound velocity pressure derivative (km/s/bar)
c 32                P-wave velocity pressure derivative (km/s/bar)
c 33                S-wave velocity pressure derivative (km/s/bar)
c 34                Adiabatic bulk modulus pressure derivative (unitless)
c 35                Shear modulus pressure derivative (unitless)
c 36                All properties of a phase or the system
c 37                Amount of a phase per unit system

c system properties from seismo:

c 1 vsys       - sys molar volume 
c 2 esys       - sys molar enthalpy 
c 3 grun       - sys gruenesien ratio
c 4 bulkm      - sys adiabatic bulk modulus
c 5 mu         - sys shear modulus
c 6 vel0       - sys sound velocity 
c 7 velp       - sys p-wave velocity
c 8 vels       - sys s-wave velocity
c 9 vp/vs      
c 10 rsys      - sys density
c 11 csys      - sys specific heat capacity
c 12 cp
c 13 alpha
c 14 beta 
c 15 entropy
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lop, icx, id, komp

      logical aprp

      double precision prop, r, gtcomp, mode(3)

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus
c----------------------------------------------------------------------
      if (lop.eq.6) then 
c                                 wt % of component 
         if (aflu.and.lflu.or.(.not.aflu)) then 
c                                 include fluid
            if (lopt(23)) then 
               prop = fbulk(icx)*atwt(icx)/psys(17)*1d2
            else 
               prop = fbulk(icx)/gtot*1d2
            end if 
         else 
c                                 exclude fluid
            if (lopt(23)) then
               prop = fbulk1(icx)*atwt(icx)/psys1(17)*1d2
            else 
               prop = fbulk1(icx)/gtot1*1d2
            end if 

         end if 
      
      else if (lop.eq.23) then 
c                                 chemical potential
         prop = mu(icx)

      else 

         if (icx.eq.0) then 
c                                 a system property is  requested:
c                                 if psys1 = 0, the system is just fluid.
            if (aflu.and.lflu.or.(.not.aflu).or.psys1(1).eq.0d0) then 
c                                 if lflu the property is to include 
c                                 fluid (if present), i.e, use psys/ptot array:
               if (lop.eq.1) then 
c 1                               Specific enthalpy (J/m3)
c                                 psys(2)  is J/mole psys(1) is J/bar/mole (volume)
                  prop = psys(2)/psys(1)*1d5
               else if (lop.eq.2) then 
c 2                               Density (kg/m3)
                  prop = psys(10)
               else if (lop.eq.3) then 
c 3                               Specific Heat capacity (J/K/m3)
                 prop = psys(12)/psys(1)*1d5
               else if (lop.eq.4) then
c 4                               Expansivity (1/K, for volume)
                  prop =  psys(13)
               else if (lop.eq.5) then 
c 5                               Compressibility (1/bar, for volume)
                  prop =  psys(14)
               else if (lop.eq.7) then 
c                                 Mode %
                  prop = 1d2
               else if (lop.ge.9.and.lop.le.15) then 
c                                 gruneisen T, K, mu, Vphi, vp, vs, vp/vs
                  prop = psys(lop-6) 
               else if (lop.eq.16) then 
c                                 specific s (j/k/m3)
                  prop = psys(15)/psys(1)*1d5
               else if (lop.eq.17) then 
c                                 S (J/K/kg)
                  prop = psys(15)/psys(1)*1d5/psys(10)
               else if (lop.eq.18) then 
c                                 H (J/kg)
                  prop = psys(2)/psys(1)*1d5/psys(10)
               else if (lop.eq.19) then
c                                 Cp (J/K/kg)
                  prop = psys(12)/psys(1)*1d5/psys(10)
               else if (lop.eq.21) then 
c                                 Poisson's ratio
                  if (psys(8).eq.0d0) then 
                     prop = 0.5d0
                  else 
                     r = (psys(7)/psys(8))**2
                     prop = 0.5d0*(r-2d0)/(r-1d0)
                  end if 
               else if (lop.eq.22) then
c                                 molar volume
                  prop = psys(1)
               else if (lop.eq.29.or.lop.eq.30) then
c                                 (Ks or Mu)_T
                  prop = psys(lop-11)
               else if (lop.eq.34.or.lop.eq.35) then
c                                 (Ks or Mu)_P
                  prop = psys(lop-14)
               else if (lop.ge.26.and.lop.le.28) then
c                                 (Vphi, Vp, Vs)_T
                  prop = psys(lop-4)
               else if (lop.ge.31.and.lop.le.33) then
c                                 (Vphi, Vp, Vs)_P
                  prop = psys(lop-6)
               else if (lop.eq.39) then 
c                                 cp/cv
                  prop = psys(28)
               end if 
            else 
c                                 fluid absent system properties:
               if (lop.eq.1) then 
c 1                               Specific enthalpy (J/m3)
                  prop = psys1(2)/psys1(1)*1d5
               else if (lop.eq.2) then 
c 2                               Density (kg/m3)
                  prop = psys1(10)
               else if (lop.eq.3) then 
c 3                               Specific Heat capacity (J/K/m3)
                  prop = psys1(12)/psys1(1)*1d5
               else if (lop.eq.4) then
c 4                               Expansivity (1/K, for volume)
                  prop =  psys1(13)
               else if (lop.eq.5) then 
c 5                               Compressibility (1/bar, for volume)
                  prop =  psys1(14)
               else if (lop.eq.7) then 
c                                 Mode %
                  prop = 1d2
               else if (lop.ge.9.and.lop.le.15) then 
c                                 m  gruneisen T, K, mu, Vphi, vp, vs, vp/vs
                  prop = psys1(lop-6) 
               else if (lop.eq.16) then
c                                 m  specific s (j/k/m3)
                  prop = psys1(15)/psys1(1)*1d5
               else if (lop.eq.17) then
c                                 S (J/K/kg)
                  prop = psys1(15)/psys1(1)*1d5/psys1(10)
               else if (lop.eq.18) then 
c                                 H (J/kg)
                  prop = psys1(2)/psys1(1)*1d5/psys1(10)
               else if (lop.eq.19) then
c                                 Cp (J/K/kg)
                  prop = psys1(12)/psys1(1)*1d5/psys1(10)
               else if (lop.eq.21) then 
c                                 Poisson's ratio
                  if (psys1(8).eq.0d0) then 
                     prop = 0.5d0
                  else 
                     r = (psys1(7)/psys1(8))**2
                     prop = 0.5d0*(r-2d0)/(r-1d0)
                  end if 
               else if (lop.eq.22) then
c                                 molar volume
                  prop = psys1(1)
               else if (lop.eq.29.or.lop.eq.30) then
c                                 (Ks or Mu)_T
                  prop = psys1(lop-11)
               else if (lop.eq.34.or.lop.eq.35) then
c                                 (Ks or Mu)_P
                  prop = psys1(lop-14)
               else if (lop.ge.26.and.lop.le.28) then
c                                 (Vphi, Vp, Vs)_T
                  prop = psys1(lop-4)
               else if (lop.ge.31.and.lop.le.33) then
c                                 (Vphi, Vp, Vs)_P
                  prop = psys1(lop-6)
               else if (lop.eq.38) then 
c                                 cp/cv
                  prop = psys1(28)
               end if 
            end if 

         else 

            if (aprp) then
c                                 call from allprp
               id = icx

            else 
c                                 normal call by phase
c                                 find the phase index
               call soltst (id,icx)

            end if 
                                 
            if (id.eq.0.and.lop.ne.7) then 
c                                 the phase is absent, set 
c                                 the property to bad_number
               prop = nopt(7)

            else if (id.ne.0) then 

               if (lop.eq.1) then
c                                 specific enthalpy 
                   prop = props(2,id)/props(1,id)*1d5
               else if (lop.eq.2) then 
c                                 density (kg/m3)
                   prop = props(10,id)
               else if (lop.eq.3) then 
c                                 specific cp 
                   prop = props(12,id)/props(1,id)*1d5
               else if (lop.eq.4) then 
c                                 expansivity
                   prop = props(13,id)
               else if (lop.eq.5) then 
c                                 compressibility
                   prop = props(14,id)
               else if (lop.eq.7) then                           
c                                 mode (%)
                   call gtmode (mode,id)
                   prop = mode(iopt(3)+1)
                   komp = id 
               else if (lop.eq.8) then 
c                                 composition (external function)
                  prop = gtcomp(id,icx,komp)
               else if (lop.ge.9.and.lop.le.15) then 
c                                 gruneisen T, K, mu, Vphi, vp, vs, vp/vs
                  prop = props(lop-6,id) 
               else if (lop.eq.16) then
c                                 specific s (j/k/m3)
                  prop = props(15,id)/props(1,id)*1d5
               else if (lop.eq.17) then 
c                                 S (J/K/kg)
                  prop = props(15,id)/props(1,id)*1d5/props(10,id)
               else if (lop.eq.18) then 
c                                 H (J/kg)
                  prop = props(2,id)/props(1,id)*1d5/props(10,id)
               else if (lop.eq.19) then
c                                 Cp (J/K/kg)
                  prop = props(12,id)/props(1,id)*1d5/props(10,id)
               else if (lop.eq.20) then                           
c                                 specific weight of a phase is the mass per
c                                 m3 of solid+melt

c                                 the number of moles of system/m3 is 1d5/psys(1)
c                                 twt is g/mol phase                                        
                  if (aflu.and.lflu.or.(.not.aflu)) then
c                                 total mode:
                     prop = props(16,id)*props(17,id)*1d2/psys(1)
                  else 
c                                 solid only mode:
                     prop = props(16,id)*props(17,id)*1d2/psys1(1)
                  end if 
               else if (lop.eq.21) then 
c                                 Poisson's ratio
                  if (props(8,id).eq.0d0) then 
                     prop = 0.5d0
                  else 
                     r = (props(7,id)/props(8,id))**2
                     prop = 0.5d0*(r-2d0)/(r-1d0)
                  end if 
               else if (lop.eq.22) then
c                                 molar volume
                  prop = props(1,id)
               else if (lop.eq.29.or.lop.eq.30) then
c                                 (Ks or Mu)_T
                  prop = props(lop-11,id)
               else if (lop.eq.34.or.lop.eq.35) then
c                                 (Ks or Mu)_P
                  prop = props(lop-14,id)
               else if (lop.ge.26.and.lop.le.28) then
c                                 (Vphi, Vp, Vs)_T
                  prop = props(lop-4,id)
               else if (lop.ge.31.and.lop.le.33) then
c                                 (Vphi, Vp, Vs)_P
                  prop = props(lop-6,id)

               else if (lop.eq.37) then                           
c                                 absolute amount of phase
c                                 per unit system
                  if (iopt(3).eq.0) then 
c                                 volume J/bar multiply by 1d-5 to get (m3)
                     prop = props(1,id)*props(16,id)
                  else if (iopt(3).eq.1) then   
c                                 mass (kg) 
                     prop = props(16,id)*props(17,id)/1d3
                  else if (iopt(3).eq.2) then 
c                                 mol 
                     prop = props(16,id)
                  end if 

               else if (lop.eq.39) then 
c                                 cp/cv
                     prop = props(28,id)

               end if 

            end if

         end if 

      end if                
 
      end 


      double precision function gtcomp (id,icx,jcomp)
c-------------------------------------------------------------------
c function comp returns icomp'th composition of phase id
c in an assemblage whose properties have been defined in routine
c seismo, the composition is defined in routine mkcomp. 

c returns -1d99 if composition jcomp is not for phase id.
c ------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision totden, comp

      integer jcomp, j, id, icx

      integer icps, jcx, jcx1, kds
      logical stol, savg, spec
      double precision rcps, a0
      common/ comps /rcps(2*k5,m13),a0(m13,2),icps(2*k5,m13),jcx(m13),
     *               jcx1(m13),kds(m13),stol(i11),savg(i11),spec(m13)
c----------------------------------------------------------------------
      if (icx.eq.kds(jcomp)) then 

         comp = a0(jcomp,1)
         totden = a0(jcomp,2)
c                                 now compute the composition:
         if (spec(jcomp)) then 
c                                 composition in terms of species:
c                                 numerator:
            do j = 1, jcx(jcomp)
               comp = comp + rcps(j,jcomp)*ysp(icps(j,jcomp),id)
            end do
c                                 denominator:
            do j = jcx(jcomp)+1, jcx1(jcomp)
               totden = totden + rcps(j,jcomp)*ysp(icps(j,jcomp),id)
            end do     

         else 
c                                 compoisiton in terms of components:
c                                 numerator:
            do j = 1, jcx(jcomp)
               comp = comp + rcps(j,jcomp)*pcomp(icps(j,jcomp),id)
            end do
c                                 denominator:
            do j = jcx(jcomp)+1, jcx1(jcomp)
               totden = totden + rcps(j,jcomp)*pcomp(icps(j,jcomp),id)
            end do     

         end if 
c                                 numerator/denominator:        
         if (totden.ne.0d0) comp = comp / totden

      else 

         comp = -1d99
 
      end if 

      gtcomp = comp

      end

      subroutine soltst (index,icx)  
c-------------------------------------------------------------------
c soltst checks for solvi, if solution icx has a solvus in 
c assemblage jd, soltsts asks the user to define a compositional
c criterion to decide which phase is relevant for property 
c contours. the index of the chosen phase is returned as index
c-------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h' 

      character cprop*6

      logical max, readyn

      integer choice, index, ksol(k5), isol, i, j, icx, 
     *        jsol, ier, phase, mode

      double precision cmin(k5) ,cmax(k5), tcomp, gtcomp

      external readyn

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character*5 cname
      common/ csta4 /cname(k5)

      integer javg,jdsol
      common/ cxt5 /javg,jdsol(k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer icps, jcx, jcx1, kds
      logical stol, savg, spec
      double precision rcps, a0
      common/ comps /rcps(2*k5,m13),a0(m13,2),icps(2*k5,m13),jcx(m13),
     *               jcx1(m13),kds(m13),stol(i11),savg(i11),spec(m13)

      save cprop, cmin, cmax, mode, max
c----------------------------------------------------------------------
      isol = 0
c                                 flag indicating whether phase speciation
c                                 is to be averaged:
      javg = 0 
c                                 how many times does the phase occur?
      do i = 1, iavar(3,ias)
         if (icx.eq.idasls(i,ias)) then
             isol = isol + 1
             jdsol(isol) = i
         end if 
      end do
 
      if (isol.ne.0) then
         index = jdsol(1)
      else 
         index = 0
      end if 
c                                 the phase doesn't occur or occurs once
      if (isol.lt.2) return 

      phase = idasls(index,ias)

c      if (lopt(34)) then
c                                 DEBUG DEBUG DANGER 
c                                 debug dump
c         write (*,*) var(1),var(2)

c         write (*,*) p,t
c         write (*,1000) isol, fname(phase), cprop
c         write (*,1040) (cname(i), i = 1, icomp)
c         do i = 1, isol 
c            write (*,1050) (pcomp(j,jdsol(i)), j = 1, icomp)
c         end do

c      end if 
c                                 if here, the phase must be a solution
10    if (.not.stol(phase)) then
c                                 immisicible phases are present (isol>1)
c                                 but there is no criterion (stol = .false.)
         stol(phase) = .true.

         if (iopt(2).eq.0) then 
            cprop = 'molar '
         else
            cprop = 'weight'
         end if
c                                tell the user
         write (*,1000) isol, fname(phase), cprop
         write (*,1040) (cname(i), i = 1, icomp)
         do i = 1, isol 
            write (*,1050) (pcomp(j,jdsol(i)), j = 1, icomp)
         end do 

         call prtptx
c                                choose action
         write (*,1001) 
         call rdnumb (nopt(8),nopt(8),choice,1,.false.)

         if (choice.eq.5) then
c                                 average the compositions, turn 
c                                 solvus testing off
            mode = 0

            call warn (99,nopt(8),i,'averaging is molar weighted, '//
     *              'this is incorrect for volumetric properties, '//
     *              'e.g., density and moduli')

            call warn (99,nopt(8),i,'phase speciation is averaged '//
     *                            'incorrectly except in the case of '//
     *                            'lagged-aqueous phase speciation')

         else
c                                 savg is a flag used only by 
c                                 werami to indicate whether 
c                                 solutions should be averaged 
c                                 within an exisiting set of solvus
c                                 criteria, initialized here.
            savg(phase) = .false.

            if (choice.lt.3) then
c                                 using a max/min criteria
               mode = 1

               if (choice.eq.1) then 
                  max = .true.
               else 
                  max = .false.
               end if 
c                                 prompt for a single composition
               write (*,1005)
               write (*,1006)

               call mkcomp (1+k5,phase)

            else if (choice.eq.3) then 
c                                 just range criteria
               mode = 2 

               write (*,1010) isol,isol-1
               write (*,1006)
            
               do i = 1, isol-1

                  call mkcomp (i+k5,phase)
c                                 get the range for the compositional
c                                 variable:
                  do 
                     write (*,1020) i
                     read (*,*,iostat=ier) cmin(i), cmax(i)
                     if (ier.eq.0) exit
                     call rerr
                  end do 


               end do 

            else
c                                 a combination of a range and a min/max
               mode = 3

               write (*,1015)
               write (*,1006) 
               write (*,1008)
c                                 the range variable 
               call mkcomp (1+k5,phase)
               
               do 
                  write (*,1020) 1
                  read (*,*,iostat=ier) cmin(1), cmax(1)
                  if (ier.eq.0) exit
                  call rerr
               end do 
c                                 the extremal variable
               write (*,1007) 
               call mkcomp (2+k5,phase)
               write (*,1025)

               if (readyn()) then 
                  max = .true.
               else 
                  max = .false.
               end if
 
            end if 

         end if 

      end if 
           

      if (mode.eq.0) then 
c                                 average immiscible compositions/props
         call avgcmp (isol,jdsol)
c                                 set javg = the number of occurences
         javg = isol

      else if (mode.eq.1) then 
c                                 min/max
         if (max) then 
            cmin(1) = -1d99
         else 
            cmin(1) = 1d99
         end if 

         do i = 1, isol

            tcomp = gtcomp (jdsol(i),idasls(jdsol(i),ias),k5+1)

            if (max) then 
               if (tcomp.ge.cmin(1)) then 
                  index = jdsol(i)
                  cmin(1) = tcomp
               end if 
            else 
               if (tcomp.lt.cmin(1)) then 
                  index = jdsol(i)
                  cmin(1) = tcomp
               end if  
            end if
 
         end do  

      else if (mode.eq.2) then 
c                                 range
c                                 identify the immiscible phase of interest     
c                                 test which phase (if any) satisfy
c                                 the compositional criteria:
         jsol = 0

         do 20 i = 1, isol
c                                 for each phase, test the isol-1 
c                                 conditions:
            do j = 1, isol-1
c                                 comp is a function that returns
c                                 the j+1th composition 
               tcomp = gtcomp (jdsol(i),idasls(jdsol(i),ias),k5+j)
c                                 the composition is not relevant
               if (tcomp.eq.-1d99) cycle
c                                 the composition is out of bounds
               if (tcomp.lt.cmin(j).or.tcomp.gt.cmax(j)) goto 20

            end do
c                                 the solution past all tests
            jsol = jsol + 1
            ksol(jsol) = jdsol(i)

20       continue 

         if (jsol.gt.1.and..not.savg(phase)) then 
c                                 two or more phases satisfy the
c                                 existing criteria
            write (*,1060) jsol,fname(phase)
            write (*,1040) (cname(i), i = 1, icomp)
            do i = 1, jsol 
               write (*,1050) (pcomp(j,ksol(i)), j = 1, icomp)
            end do 

            write (*,1070)
c                                 get choice
            call rdnumb (tcomp,0d0,choice,1,.false.)

            if (choice.eq.2) then 
c                                 2 - average within existing criterion 
               savg(phase) = .true.

            else if (choice.eq.3) then 
c                                 3 - ignore and hope for the best
               call avgcmp (jsol,ksol)

            else 
c                                 not 2 or 3, redefine the criterion
               stol(phase) = .false.
               goto 10             

            end if 

         end if 
c                                 user has elected to average within
c                                 existing criterion
         if (jsol.gt.1.and.savg(phase)) call avgcmp (jsol,ksol) 

         if (jsol.eq.0) then 
            index = jdsol(1)
         else 
            index = ksol(1)
         end if 

      else 
c                                 min/max + range, first check range:
         jsol = 0

         do i = 1, isol
c                                 for each phase, test the range
            tcomp = gtcomp (jdsol(i),idasls(jdsol(i),ias),k5+1)
c                                 the composition is not relevant or out
c                                 of bounds, cycle
            if (tcomp.eq.-1d99.or.tcomp.lt.cmin(1).or.tcomp.gt.cmax(1))
     *         cycle 

c                                 the phase passed all tests
            jsol = jsol + 1
            ksol(jsol) = jdsol(i)

         end do  
c                                 if no phase passed range criterion
c                                 load them all in to ksol and use the
c                                 min/max, could add a 1 time warning.
         if (jsol.eq.0) then 
            jsol = isol
            do i = 1, jsol
               ksol(i) = jdsol(i)
            end do 
         end if 
c                                 min/max
         if (max) then 
            cmin(2) = -1d99
         else 
            cmin(2) = 1d99
         end if 

         do i = 1, jsol

            tcomp = gtcomp (ksol(i),idasls(ksol(i),ias),k5+2)

            if (max) then 
               if (tcomp.ge.cmin(2)) then 
                  index = jdsol(i)
                  cmin(2) = tcomp
               end if 
            else 
               if (tcomp.lt.cmin(2)) then 
                  index = jdsol(i)
                  cmin(2) = tcomp
               end if  
            end if
 
         end do  

      end if 

1000  format (/,i1,' immiscible phases of ',a,/,'coexist with the ',
     *        'following ',a,' compositions:',/)
1001  format (/,'Identify the phase of interest by:',//,3x,
     *        '1 - the maximum value of a composition [default].',/,3x,
     *        '2 - the minimum value of a composition.',/,3x,
     *        '3 - the range of one or more compositions.',/,3x,
     *        '4 - a combination of the above.',/,3x,
     *        '5 - average the compositions of immiscible phases.')
1005  format (/,'The following prompts define the composition C[1] to ',
     *         'be used to identify the',/,'phase of interest.')
1006  format (/,'NOTE: discriminatory criteria are only applied when ',
     *        'immiscible phases coexist.',/,'If only one phase of a',
     *        'solution is stable, then data for this phase is output',/
     *       ,'regardless of whether the phase meets the criteria ',
     *        'specified here.')
1007  format (/,
     *     'The following prompts refer to the extremal variable C[2]')
1008  format (/,
     *        'The following prompts refer to the range variable C[1]')
1010  format (/,'The following prompts define the compositional ',
     *        'variable(s) (C[i]) to be used',/,
     *        'to identify the phase of interest.',//,
     *        'As there are ',i1,' coexisting phases',
     *        ' you will be prompted',/,'for ',i1,
     *        ' compositional variable(s).')
1015  format (/,'In this mode you will be prompted to specify: ',//,3x,
     *        '- a compositional variable C[1] that is used to define a'
     *       ,' range criterion',/,3x,
     *        '- a compositional variable C[2] that is used to define a'
     *       ,'minimum/maximum criterion',//,
     *        'The range criterion is applied first to retrict the ',
     *        'identity of the phase of interest.',//,'The minimum/',
     *        'maxumum criterion is then used to select the phase of',
     *        'interest from',/,'those phases that satisfy the range ',
     *        'criterion. NOTE: if no phases satisfy',/,'the range ',
     *        'criterion, then the phase of interest will be ',
     *        'identified by the',/,'minimum/maximum criterion; this',
     *        'may lead to surprising results.')
1020  format (/,'Enter the range (minimum, maximum) of C[',i1,'] that ',
     *        'defines the phase of interest:')
1025  format (/,'Identify the phase of interest as the phase with the ',
     *        'highest value of C[2] (Y/N)?')
1040  format (/,4x,20(a,4x))
1050  format (3x,20(f7.3,2x))
1060  format (/,i1,' coexisting phases of ',a,' satisfy your ',
     *        ' compositional criteria',/,'with the compositions:')
1070  format (/,'Choose an option:',/,3x,
     *        '1 - redefine the compositional criteria [default].',/,3x,
     *        '2 - average the compositions of all phases that',
     *        ' meet the existing criterion.',/,3x,
     *        '3 - ignore this instance.')
1080  format (/,'Identify the phase of interest as:',/,3x,
     *        '1 - the phase with the maximum value of C[',i1,'].',/,3x,
     *        '2 - the phase with the minimum value of C[',i1,'].',/,3x,
     *        '3 - the phase with C[',i1,'] in a specific range ',
     *        '[default].',/)

      end 

      subroutine mode3 
c----------------------------------------------------------------------
c sample data on an x-y path defined by interactive user input 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node, ok, readyn

      integer i, j, icurve, ivi, ivd, iord, ipts, jpts, ier, k(2), dim

      double precision coef(0:10), dxy(2), xyp(2,2), s, d

      character*100 n5name, n6name, text*320

      external readyn

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer ivar,ind
      common/ cst83 /ivar,ind
c----------------------------------------------------------------------
      node = .false.
      dim = 1
c                                 set up path information
      icurve = 0 
c                                 ask if non-linear path
10    write (*,1200)

      if (readyn()) then 

         icurve = 1
c                                 select independent variable:
5        write (*,1160) (i,vnm(i),i= 1, 2)
         read (*,*,err=5) ivi

         if (ivi.eq.1) then 
            ivd = 2
         else
            ivi = 2
            ivd = 1
         end if

         write (*,1210) vnm(ivd),vnm(ivi)
         read (*,*) iord

         do i = 0, iord
            write (*,1220) i
            read (*,*) coef(i)
         end do

         write (text,1350) vnm(ivd),(coef(i),vnm(ivi),i,i=0,iord)
         call unblnk (text)
         write (*,1340) text
c                                 ask if ok.
         write (*,1320)

         if (readyn()) goto 10
c                                 it's ok.
         dxy(ivi) = vmx(ivi)-vmn(ivi)
         dxy(ivd) = vmx(ivd)-vmn(ivd)
         xyp(ivi,1) = vmn(ivi)
         xyp(ivd,1) = vmn(ivd)

      else 
c                                 linear path
         ivi = 1
         ivd = 2
         i = 1 

         do 

            write (*,1140) i, vnm(1), vnm(2)
            read (*,*,iostat=ier) xyp(1,i), xyp(2,i)
            if (ier.ne.0) cycle

            ok = .true.

            do j = 1, 2

               if (vmn(j).lt.vmx(j)) then 

                  if (xyp(j,i).ge.vmn(j).and.xyp(j,i).le.vmx(j)) cycle 
 
                  write (*,1010) vnm(j),vmn(j),vmx(j)
                  ok = .false. 

               else

                  if (xyp(j,i).ge.vmx(j).and.xyp(j,i).le.vmn(j)) cycle  
                  write (*,1010) vnm(j),vmn(j),vmx(j)
                  ok = .false.

               end if
 
            end do 

            if (.not.ok) cycle

            i = i + 1

            if (i.eq.3) then 

               do j = 1, 2 
                  dxy(j) = xyp(j,2) - xyp(j,1)
               end do 

               if (dxy(1).eq.0d0.and.dxy(2).eq.0d0) then

                  write (*,*) 
     *               'initial and final coordinates cannot be identical'
                  i = 1
                  cycle

               else if (dxy(1).eq.0d0) then 
c                                 parallel to the x-axis
                     ivi = 2
                     ivd = 1
                     s = 0d0

               else if (dxy(2).eq.0d0) then 
c                                 parallel to y axis
                  s = 0d0 

               else 
c                                 compute slope
                  s = dxy(ivd)/dxy(ivi)

               end if 

               exit 

            end if 

         end do 

      end if 
c                                 set up counters, pointers:
      do 

         write (*,1150) 
         read (*,*,iostat=ier) ipts
         if (ipts.lt.2) ipts = 2
         if (ier.eq.0) exit 

      end do 

      d = dxy(ivi)/dfloat(ipts - 1)
c                                 ind could be set by getind              
      ind = ivi

      jpts = 0 
c                                 select properties
      call chsprp
c                                 write file header
      call tabhed (n5,dxy,dxy,k,dim,n5name,n6name)
c                                 compute properties
      do i = 1, ipts

         var(ivi) = xyp(ivi,1) + dfloat(i-1)*d

         if (icurve.eq.0) then 

            var(ivd) = xyp(ivd,1) + s*(var(ivi)-xyp(ivi,1))
            jpts = 1

         else
 
            var(ivd) = 0d0

            do j = 0, iord
               var(ivd) = var(ivd) + coef(j)*var(ivi)**j
            end do

            if (var(ivd).le.vmx(ivd).and.var(ivd).ge.vmn(ivd)) then 
               jpts = jpts + 1
            else 
               cycle
            end if
 
         end if 

         call polprp (dim)

      end do 

      if (jpts.eq.0) call error (72,s,i,'your polynomial yields no cond'
     *//'itions with the coordinate frame of the calculation')

      call finprp (dim,n5name,n6name,node) 

1010  format (/,'The plot file range for ',a,' is ',g12.4,' - ',g12.4,
     *        /,'Try again:',/)
1140  format (/,'Enter endpoint ',i1,' (',a,'-',a,') coordinates:')
1150  format (/,'How many points along the profile?')
1160  format (/,'Select independent variable: ',2(/,1x,i1,' - ',a))
1200  format (/,'Construct a non-linear profile (y/n)?')
1210  format (/,'Profile must be described by the function',/,a,
     *        ' = Sum ( c(i) * ',a,' ^i, i = 0..n)',/,'Enter n (<10)')
1220  format (/,'Enter c(',i2,')')
1320  format (/,'Change the profile (Y/N)?')
1330  format (/,'Your polynomial does not yield conditions within',
     *          ' the computational coordinate frame.',/)
1340  format (/,'Your polynomial is:',/,a)
1350  format (a,'=',11('+(',g12.6,')*',a,'^',i1))

      end 

      subroutine mode31
c----------------------------------------------------------------------
c sample data in a 1d section 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node

      integer i, ipts, k(2), dim 

      character*100 n5name, n6name

      double precision xy(2)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer ivar,ind
      common/ cst83 /ivar,ind



c----------------------------------------------------------------------
      node = .false.
      dim = 1
c                                 select properties:
      call chsprp
c                                 see pre-6.8.4 for limited range 
c                                 plotting. here assume only two 
c                                 cases icopt = 12 = nodal coordinate
c                                 is the second variable, and other-
c                                 wise the nodal coordinate is 
c                                 set the independent output variable
      ind = 1

      if (icopt.eq.12) then 
         ipts = iopt(36) + 1
      else 
         ipts = int( (vmx(ind)-vmn(ind))/dvr(ind) + 1)
      end if 
c                                 name and open plot file, write header 
      call tabhed (n5,xy,xy,k,dim,n5name,n6name)

      do i = 1, ipts

         var(ind) = vmn(ind) + dfloat(i-1) * dvr(ind)

         call polprp (dim)

      end do 

      call finprp (dim,n5name,n6name,node) 

      end 

      subroutine mode4
c----------------------------------------------------------------------
c sample data on an x-y path defined by input from a file
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node

      integer i, icoors, jmode, ixy, inc, ierr, k(2), dim

      double precision r(2), tmin1, tmax1,a1, b1, c1, d1, x0, x1,  
     *                 dt1, x, y, xx(5*l5), yy(5*l5), pmin, pmax

      character*100 n5name, n6name, dname

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar 

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      logical oned
      common/ cst82 /oned
c----------------------------------------------------------------------
      node = .false.
      dim = 1

      do 

         if (oned) then 

            write (*,1280) 
            jmode = 2 

         else 
c                                 select input type
            write (*,1290) 
            read (*,*) jmode   
            write (*,1300) 

         end if 

         read (*,'(a)') dname
         open (n7,file=dname,iostat=ierr)
         if (ierr.eq.0) exit 
         write (*,*) 'No such data file as: ',dname

      end do 

      if (jmode.eq.1) then 
c                                 points from a polynomial
c                                 the association of "x" and "y"
c                                 with the actual variables of the
c                                 diagram is determined by the flag
c                                 ixy 
         read (n7,*) pmin,pmax,ixy
         read (n7,*,iostat=ierr) tmin1,tmax1,dt1,a1,b1,c1,d1
         if (ierr.ne.0) return

         if (dt1.lt.0d0) then 
            x0 = tmax1
            x1 = tmin1
         else 
            x0 = tmin1
            x1 = tmax1
         end if 

         x = x0
c                                 select properties
         call chsprp 
c                                 write plot file header
         call tabhed (n5,r,r,k,dim,n5name,n6name)
c                                 set counter for cumulative modes
         i = 1

         do 

            y = a1 + b1*x + c1*x**2 + d1*x**3

            if (y.le.pmax.and.y.ge.pmin) then 
               if ((dt1.gt.0d0.and.x.le.x1).or.
     *             (dt1.lt.0d0.and.x.ge.x1)) then
c                                 condition is in bounds
                  if (ixy.eq.0) then 
                     var(1) = x
                     var(2) = y 
                  else
                     var(1) = y
                     var(2) = x 
                  end if 

                  call polprp (dim)

               end if 

            end if 

            x = x + dt1

            if (x.gt.x1) exit 

         end do  

         call finprp (dim,n5name,n6name,node) 

      else 
c                                   points from a data file:
         icoors = 1

         do 

            if (oned) then
   
               read (n7,*,iostat=ierr) xx(icoors)
               if (ierr.ne.0) exit 

               yy(icoors) = vmn(2)

            else

               read (n7,*,iostat=ierr) xx(icoors),yy(icoors)
               if (ierr.ne.0) exit 

            end if 

            if (xx(icoors).lt.vmn(1).or.xx(icoors).gt.vmx(1)) cycle 
            if (yy(icoors).lt.vmn(2).or.yy(icoors).gt.vmx(2)) cycle 

            icoors = icoors + 1

            if (icoors.gt.5*l5) call error (69,xx(1),icoors+100,
     *                                     'L5 (MODE4)')

         end do 
               
         icoors = icoors - 1

         close (n7)

         if (icoors.eq.0) then 
            write (*,*) 'file contains no points within data bounds'
            return
         end if 

         write (*,1310) icoors
         read (*,*) inc
c                                 select properties            
         call chsprp
c                                 write plot file header
         call tabhed (n5,r,r,k,dim,n5name,n6name)

         do i = 1, icoors, inc

            var(1) = xx(i)
            var(2) = yy(i)

            call polprp (dim)

         end do  

         call finprp (dim,n5name,n6name,node) 

      end if 

1280  format (/,'Enter the name of the file containing the path',
     *          ' ordinates:',/)
1290  format (/,'Path will be described by:',/,
     *          '   1 - a file containing a polynomial function',/,
     *          '   2 - a file containing a list of x-y points',/,
     *          'Enter 1 or 2:'/) 
1300  format (/,'Enter the file name:',/)
1310  format (/,'File contains ',i5,' points',/,
     *          'every nth plot will be plotted, enter n:',/)

      end

      subroutine avgcmp (isol,jdsol)
c-------------------------------------------------------------------
c makes the average composition and properties of from isol compositions
c of the phase indexed by the array jdsol
c-------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h' 

      integer jdsol(k5), index, isol, i, j

      double precision x(k5), ntot

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------
c                                 average immiscible compositions
c                                 get mole fractions
      ntot = 0d0 
      index = jdsol(1)

      do i = 1, isol
         ntot = ntot + props(16,jdsol(i))
      end do 

      do i = 1, isol
         x(i) = props(16,jdsol(i))/ntot
      end do 
c                                 set composition
      do j = 1, icomp

         pcomp(j,index) = x(1) * pcomp(j,index)

         do i = 2, isol
            pcomp(j,index) = pcomp(j,index) + x(i)*pcomp(j,jdsol(i))
         end do 
      end do 
c                                 set physical properties assuming molar
c                                 weighting (this is wrong for volumetric
c                                 properties!!!). 
      do j = 1, 17

         if (j.eq.16) cycle

         props(j,index) = x(1) * props(j,index)

         do i = 2, isol 
            props(j,index) = props(j,index) + x(i)*props(j,jdsol(i))
         end do 

      end do

      props(16,index) = ntot

      end

      subroutine allmod
c----------------------------------------------------------------
c computes modes of all stable phases, i.e., over the entire range
c of physical conditions. 
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, id, jk, ind(i11), jnd(i11), nsol, knd(i11), ksol

      logical stble(i11), quit

      double precision mode(3), smode(i11), dinc 

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer idsol,nrep,nph
      common/ cst38/idsol(k5,k3),nrep(k5,k3),nph(k3)

      save ind, knd, stble, ksol, nsol
c----------------------------------------------------------------------
      do i = 1, iprop
         prop(i) = nopt(7)
         smode(i) = 0d0
      end do 

      id = 0

      do i = 1, nph(ias)

         jk = 0

         do j = 1, istab
c                               jnd points to the column before 
c                               the j'th solution
            jnd(j) = jk

            if (idstab(j).eq.idsol(i,ias)) then

               do k = 1, nrep(i,ias)

                  id = id + 1

                  call gtmode (mode,id)
                  prop(jk+k) = mode(iopt(3)+1)
                  smode(j) = smode(j) + prop(jk+k)

               end do

            end if
c                                mode column pointer
            jk = jk + nstab(j)
 
         end do

      end do

      if (lopt(2)) then

         if (first) then

            knd(1) = 0

            do i = 1, istab
               ind(i) = i
               if (i.eq.1) cycle 
               knd(i) = knd(i-1) + nstab(i-1)
            end do 

            call rankem (smode,ind,istab,istab)

            nsol = nph(ias)

            first = .false.

         else

            do i = 1, istab

               if (.not.stble(i).and.smode(i).gt.0d0) then
c                                 check if it's already in ind:
                 quit = .false.

                 do j = 1, nsol
                    if (ind(j).eq.i) then
                       quit = .true.
                       exit
                    end if
                 end do 

                 if (quit) then 
                    stble(i) = .true.
                    cycle 
                 end if 

                 if (nph(ias).eq.ksol.and.lopt(45)) then
c                                 probably univariant, find the phase that
c                                 was stable, but isn't anymore
                    do j = 1, nsol

                       if (stble(ind(j)).and.
     *                     smode(ind(j)).gt.0d0) cycle
c                                 shift all higher indices up
                       do k = nsol, j, -1
                          ind(k+1) = ind(k)
                       end do

                       ind(j) = i
                       nsol = nsol + 1

                       exit 

                    end do

                 else
c                                 probably low variance
                    nsol = nsol + 1
                    ind(nsol) = i

                  end if

                  stble(i) = .true.

               end if 

            end do
c                                 set modes of phase that have dissappeared
c                                 to real zeros (not nopt(7)).
            do i = 1, nsol

               if (stble(ind(i)).and.smode(ind(i)).gt.0) cycle

               do j = 1, nstab(i)
                  prop(jnd(ind(i))+j) = 0d0
               end do 

            end do

         end if

         do i = 1, istab
            if (smode(i).eq.0d0) then 
               stble(i) = .false.
            else 
               stble(i) = .true.
            end if
         end do

         ksol = nph(ias)
c                                 convert to cumulative modes if
c                                 requested
         dinc = 0d0

         do j = nsol, 1, -1

            i = ind(j)
c                                 in case of multiple occurences run
c                                 over all the occurrences of ind(j):
            do k = 1, nstab(i)

               if (isnan(prop(knd(i)+k))) exit

               dinc = prop(knd(i)+k) + dinc

               prop(knd(i)+k) = dinc

            end do 

         end do

      end if
  
      end 

      subroutine gtmode (mode,id) 
c----------------------------------------------------------------
c function to extract vol/wt/mol mode (%) of phase id 
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision mode(3)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)
c----------------------------------------------------------------
      if (aflu.and.lflu.or.(.not.aflu).or.psys1(1).eq.0d0) then
c                     total mode:
c                     volume fraction
         mode(1) = props(1,id)*props(16,id)/psys(1)*1d2
c                     weight fraction 
         mode(2) = props(16,id)*props(17,id)/psys(17)*1d2
c                     mol fraction
         mode(3) = props(16,id)/psys(16)*1d2

      else 
c                     solid only mode:
         if (fluid(id)) then 
            mode(1) = 0d0
            mode(2) = 0d0
            mode(3) = 0d0
         else 
c                     volume fraction
            mode(1) = props(1,id)*props(16,id)/psys1(1)*1d2
c                     wt fraction
            mode(2) = props(16,id)*props(17,id)/psys1(17)*1d2
c                     mol fraction
            mode(3) = props(16,id)/psys1(16)*1d2
         end if 

      end if 

      end 

      subroutine outmod (dim,n6name,node)
c----------------------------------------------------------------
c reformat output from "all_modes" requests according to the 
c options set in perplex_option.dat and dimension dim. 
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 ynm, n6name*100, rec*1

      logical node, bdnum

      integer i, j, k, ipt, dim, ier

      double precision ymx,ymn,xl(i11),yl(i11),x(3),dx,dy(i11),xmx,xmn,
     *                 mzero

      character vnm*8
      common/ cxt18a /vnm(l3) 

      double precision vip
      common/ cst28 /vip(l2,k2)

      character*162 title
      common/ csta8 /title(4)

      integer ivar,ind
      common/ cst83 /ivar,ind

      integer inv
      character dname*14, titl1*162
      common/ cst76 /inv(i11),dname(i11),titl1

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      external bdnum
c----------------------------------------------------------------------
      if (dim.eq.1) then 

         rewind (n5)
c                                 skip header lines
         do i = 1, 9
            read (n5,'(a)') rec
         end do 

         open (n6,file=n6name)
c                                 read file to get limits, this 
c                                 must also be done for x if 
c                                 ind ne 1.
         ymn = 1d30
         ymx = -1d30
         xmn = ymn
         xmx = ymx
         ynm = 'y axis'
         ipt = 0 

         do i = 1, iprop
            dy(i) = -1d30
            yl(i) = -1d30
         end do 
c                                 read the data to get the range  
         do           

            read (n5,*,iostat=ier) (x(i),i=1,ivar), (prop(i),i=1,iprop)

            ipt = ipt + 1
            if (ipt.gt.k2) call error (999,ymn,ipt,'OUTMOD')

            if (node) then 
               vip(4,ipt) = ipt
            else 
               vip(4,ipt) = x(ind)
            end if 

            if (ier.ne.0) exit 

            do i = 1, iprop
c                                 don't count bad_number value modes
               if (bdnum(prop(i))) cycle

               if (x(ind).gt.xmx) xmx = x(ind)
               if (x(ind).lt.xmn) xmn = x(ind)
               if (prop(i).gt.ymx) ymx = prop(i)
               if (prop(i).lt.ymn) ymn = prop(i)
c                                 get a label position
               if (lopt(2)) then 
c                                 cumulative mode
                  if (i.eq.1) then 

                     if (prop(i).gt.2d0*yl(i)) then
                        yl(i) = prop(i)/2d0
                        xl(i) = x(ind)
                     end if 

                  else 
c                                 check for NaN badnumber value
                     if (bdnum(prop(i-1))) then
                        mzero = 0d0
                     else
                        mzero = prop(i-1)
                     end if 

                     if (prop(i)-mzero.gt.dy(i)) then
                        dy(i) = prop(i) - mzero
                        yl(i) = (prop(i) + mzero)/2d0 
                        xl(i) = x(ind)
                     end if
 
                  end if 

               else

                  if (prop(i).gt.yl(i)) then 
                     yl(i) = prop(i)
                     xl(i) = x(ind)
                  end if 

               end if 

            end do 

         end do

         dx = (xmx-xmn)/1d1
c                                 header of psvdraw file
         write (n6,1000) title,xmx,xmn,ymx,ymn,vnm(ind),ynm
c                                 make each column into a curve use of
c                                 the vip array might cause problems 
c                                 if someone uses chemical potential
c                                 variables and calculates multiple
c                                 property sets, as vip(4,5) will be 
c                                 reset and are not reread. 
         do i = 1, iprop

            rewind (n5)
c                                 skip header lines
            do j = 1, 9
               read (n5,'(a)') rec
            end do 

            do j = 1, ipt       
                   
               read (n5,*,iostat=ier) (x(k),k=1,ivar),
     *                                (prop(k),k=1,iprop)
               vip(5,j) = prop(i)

            end do
c                                 output the curve
            write (n6,1010) ipt*2,i,dname(i)
            write (n6,*) (vip(4,j),vip(5,j),j=1,ipt)

         end do
c                                 write trailers for psipts, pscurv
         write (n6,1010) 1,1,'trailer'
         write (n6,*) '0 0'
c                                 write label coordinates and text
         do i = 1, iprop
            if (lopt(2)) then 
               write (n6,*) xl(i)-dx,yl(i)
            else
               write (n6,*) xl(i),yl(i)
            end if 
            write (n6,'(a)') dname(i)
         end do 

      end if

      close (n6)

1000  format ('1',/,'0 0 0',/,'0 0 0 0 0 0',4(/,a162),/,'2 1 2 0 0',/,
     *        '0 0 0 0. 0. 0. 0. 0.',/,4(g14.7,1x),/,a,/,a)
1010  format (i5,' 1 ',i3,' 1 1 1 1 1 1 ',/,'0. ',a)

      end 

      subroutine allprp (dim)
c----------------------------------------------------------------
c output spreadsheet format properties (lop = 36 or 38).
c   icx = 0   - output just system properties
c   icx = 999 - system + phase properties
c   else      - phase icx properties
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, id, dim, dummy

      double precision mode(3), fwt, cprp(i11)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer ivar,ind
      common/ cst83 /ivar,ind

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      character pname*14
      common/ cxt21a /pname(k5)

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      double precision atwt
      common/ cst45 /atwt(k0)

      save cprp
c----------------------------------------------------------------------
      if (kop(1).eq.38) then 
c                                 custom property choices
         if (kcx(1).eq.999.or.kcx(1).eq.0) then 
c                                 get system props
            do i = 1, iprop
               call getprp (prop(i),kop(i+1),0,0,.true.)
            end do

            tname = 'system        '

            call outprp (dim)

         end if 

         if (kcx(1).ne.0) then 
c                                 properties of phases
            do j = 1, ntot
             
               if (kcx(1).eq.999) then
c                                 phemgp output, the phase
c                                 will be identified by name
                  id = j

                  do i = 1, iprop
                     call getprp (prop(i),kop(i+1),id,dummy,.true.)
                  end do

               else
c                                 kcx is true phase id, locate
c                                 the phase in the assemblage:
                  dummy = 0

                  do k = 1, ntot
                     if (kkp(k).eq.kcx(1)) then
                        dummy = dummy + 1
                        id = k
                     end if
                  end do

                  if (dummy.eq.0) then 
c                                 the phase is not present assign badnum
                     do i = 1, iprop
                        prop(i) = nopt(7)
                     end do

                  else
c                                 the solution is present
                     do i = 1, iprop
                        call getprp (prop(i),kop(i+1),id,dummy,.true.)
                     end do
c                                 but as more than one phase, write warning
                     if (dummy.gt.1) call warn (41,prop(1),dummy,
     *                                                        pname(id))

                  end if

               end if

               tname = pname(id)

               call outprp (dim)

               if (kcx(1).ne.999) exit 

            end do       

         end if  

      else
c                                 "all property" option, lop.eq.36
         if (first) then 
c                                 initialize cumulative props
            do i = 1, iprop
               cprp(i) = 0d0
            end do 

            first = .false.

         end if

         if (kcx(1).eq.999.or.kcx(1).eq.0) then 
c                                 only system properties requested:
            if (aflu.and.lflu.or.(.not.aflu)) then 
c                                 normal properties
               do i = 1, i8
                  prop(i) = psys(i)
               end do 
c                                 bulk composition and excess charge
               fwt = psys(17)

               do i = i8+4, i8+3+icomp
c                                 absolute molar bulk composition
                  prop(i) = fbulk(i-i8-3)

               end do 

            else
c                                 exclude fluid:
c                                 normal properties
               do i = 1, i8
                  prop(i) = psys1(i)
               end do 
c                                 bulk composition 
               fwt = psys1(17)

               do i = i8+4, i8+3+icomp
                  prop(i) = fbulk1(i-i8-3)
               end do 

            end if
c                                 chemical potentials
            do i = i8+icomp+4, i8+2*icomp+3
               prop(i) = mu(i-i8-3-icomp)
            end do 

c                                 dummy modes
            do i = i8+1, i8+3
               prop(i) = 1d2
            end do 
c                                 convert bulk composition if necessary
            if (lopt(23)) then 

               if (oxchg) then 
c                                 net nominal oxidation state 
                  prop(iprop) = 0d0

                  do i = i8+4, i8+3+icomp

                     prop(iprop) = prop(iprop) + cox(i-i8-3) * prop(i)

                  end do

               else 

                  prop(iprop) = nopt(7) 

               end if

               do i = i8+4, i8+3+icomp
c                                 convert to mass:
                  if (lopt(41)) then
c                                 absolute mass
                     prop(i) = prop(i)*atwt(i-i8-3)
                  else 
c                                 relative (%)
                     prop(i) = prop(i)*atwt(i-i8-3)/fwt*1d2
                  end if

               end do

            end if

            tname = 'system        '

            call outprp (dim)

         end if 

         if (kcx(1).ne.0) then
c                                 phase properties requested,
c                                 if kcx(1) = 999 then all phases
c                                 otherwise kcx(1) is the phase 
c                                 pointer
            do j = 1, ntot
c                                 properties of all phases
               if (kcx(1).eq.999) then 
                  id = j 
               else
c                                 all properties of a specific phase
c                                 find the phase index
                  call soltst (id,kcx(1))

                  if (id.eq.0) then 
c                                 the phase is not present
                     do i = 1, iprop
                        prop(i) = nopt(7)
                     end do 
                      
                     tname = 'not stable'

                     call outprp (dim) 

                     goto 99 

                  end if 

               end if 
c                                 normal properties
               do i = 1, i8
                  prop(i) = props(i,id)
               end do 
c                                 compute modes
               call gtmode (mode,id)

               do i = 1, 3
                  prop(i8+i) = mode(i)
               end do 
c                                 bulk phase compositions are wt% if iopt(2) = 1 
c                                 and molar if iopt(2) = 0.
               do i = i8+4, i8+3+icomp
                  prop(i) = pcomp(i-i8-3,id)
               end do 

               if (oxchg.and.iopt(2).eq.0) then 
c                                 net nominal oxidation state 
                  prop(iprop) = 0d0

                  do i = i8+4, i8+3+icomp

                     prop(iprop) = prop(iprop) + cox(i-i8-3) * prop(i)

                  end do

                  prop(iprop) = prop(iprop) * prop(16)

               else

                  prop(iprop) = nopt(7)

               end if
c                                 convert to absolute mol/mass if requested
               if (lopt(41)) then 

                  do i = i8+4, i8+3+icomp
c                                 prop(16) is the molar amount of the phase
                     prop(i) = prop(i) * prop(16)
c                                 mass ammount, undo percent
                     if (iopt(2).eq.1) prop(i) = prop(i)/2d2

                  end do

               end if 
c                                 chemical potentials
               do i = i8+icomp+4, i8+2*icomp+3
                  prop(i) = mu(i-i8-3-icomp)
               end do 

               tname = pname(id)
c                                 convert if absolute amounts requested 

c                                 convert compositions to cumulative, this only
c                                 makes sense for absolute amounts of a fractionated
c                                 phase
               if (lopt(42).and.kcx(1).ne.999.and.dim.eq.1) then

                  do i = i8+4, i8+3+icomp
                     prop(i) = cprp(i) + prop(i)
                     cprp(i) = prop(i)
                  end do

                  prop(iprop) = cprp(iprop) + prop(iprop)
                  cprp(iprop) = prop(iprop)

               end if 

               call outprp (dim) 

               if (kcx(1).ne.999) exit 

            end do

         end if

      end if

99    end 

      subroutine getind
c----------------------------------------------------------------
c get the plotting variable index (ind) for 1-d property plots
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier

      character vnm*8
      common/ cxt18a /vnm(l3) 

      integer ivar,ind
      common/ cst83 /ivar,ind
c----------------------------------------------------------------------
c                                 choose plotting variable
      write (*,1000) vnm(1)

      do 
         write (*,1030) (ier,vnm(ier),ier=1,ivar)
         read (*,'(bn,i80)',iostat=ier) ind
         if (ier.ne.0) cycle
         if (ind.ne.2.and.ind.ne.3) ind = 1
         exit 
      end do  

1000  format (/,'The independent variable for this calculation is: ',a)
1030  format (/,'Choose the independent variable for data plots:',/,
     *       4x,i1,' - ',a,' [default]',6(/,4x,i1,' - ',a))
 
      end 

      subroutine chsprp 
c----------------------------------------------------------------
c chsprp asks the user to choose properties to be extracted, it
c then creates a list of names for the properties (dname).

c   lop/kop  - flag indicating the property chosen
c   icx/kcx  - if lop = 6, the component chosen
c   icx/kcx  - if lop > 6, the identity of the solution chosen,
c              icx = -1 if a solution is not chosen.

c   lflu/kfl - .true. include fluids for bulk props.

c 1                 Specific enthalpy (J/m3)',
c 2                 Density (kg/m3)',
c 3                'Specific Heat capacity (J/K/m3)',
c 4                'Expansivity (1/K, for volume)',
c 5                'Compressibility (1/bar, for volume)',
c 6                'Weight percent of a component',
c 7                'Mode (Vol %) of a compound or solution',
c 8                'Composition of a solution'
c 9                 Grueneisen thermal ratio',
c 10               'Adiabatic bulk modulus (bar)',
c 11               'Adiabatic shear modulus (bar)'
c 12                Sound velocity (km/s)
c 13               'P-wave velocity (km/s)',
c 14                S-wave velocity (km/s)',
c 15                Vp/Vs
c 16               'Specific Entropy (J/K/m3)'
c 17               'Entropy (J/K/kg)'
c 18               'Enthalpy (J/kg)'
c 19               'Heat Capacity (J/K/kg)'
c 20               'Specific mass (kg/m3) of a phase'
c 21               'Poisson's Ratio'
c 22               'Molar Volume (J/bar)'
c 23                Dependendent potentials (J/mol)
c 24                Assemblage index
c 25                Modes of all phases (wt or vol%)
c 26                Sound velocity temperature derivative (km/s/K)
c 27                P-wave velocity temperature derivative (km/s/K)
c 28                S-wave velocity temperature derivative (km/s/K)
c 29                Adiabatic bulk modulus temperature derivative (bar/K)
c 30                Shear modulus temperature derivative (bar/K)
c 31                Sound velocity pressure derivative (km/s/bar)
c 32                P-wave velocity pressure derivative (km/s/bar)
c 33                S-wave velocity pressure derivative (km/s/bar)
c 34                Adiabatic bulk modulus pressure derivative (unitless)
c 35                Shear modulus pressure derivative (unitless)
c 36                All properties of a phase or the system
c 37                Absolute amounts
c 38                Multiple property grid for system and phases
c 39                heat capacity ratio (cp/cv)
c 40                Lagged or Back-calculated aqueous solute chemistry
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical phluid, readyn

      integer i, j, icx, kprop, ier, lop, komp, mprop

      parameter (kprop=40)

      character propty(kprop)*60, pname*10

      external readyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer ivar,ind
      common/ cst83 /ivar,ind

      character cname*5
      common/ csta4  /cname(k5)

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      integer icps, jcx, jcx1, kds
      logical stol, savg, spec
      double precision rcps, a0
      common/ comps /rcps(2*k5,m13),a0(m13,2),icps(2*k5,m13),jcx(m13),
     *               jcx1(m13),kds(m13),stol(i11),savg(i11),spec(m13)

      save propty

      data propty/'Specific Enthalpy (J/m3)',
     *            'Density (kg/m3)',
     *            'Specific heat capacity (J/K/m3)',
     *            'Expansivity (1/K, for volume)',
     *            'Compressibility (1/bar, for volume)',
     *            'Composition (Mol, Mass, or Wt%) of the system',
     *            'Mode (Vol, Mol, or Wt proportion) of a phase',
     *            'Composition (Mol, Mass, or Wt%) of a solution phase',
     *            'Grueneisen thermal ratio',
     *            'Adiabatic bulk modulus (bar)',
     *            'Adiabatic shear modulus (bar)',
     *            'Sound velocity (km/s)',
     *            'P-wave velocity (Vp, km/s)',
     *            'S-wave velocity (Vs, km/s)',
     *            'Vp/Vs',
     *            'Specific entropy (J/K/m3)',
     *            'Entropy (J/K/kg)',
     *            'Enthalpy (J/kg)',
     *            'Heat Capacity (J/K/kg)',
     *            'Specific mass of a phase (kg/m3-system)',
     *            'Poisson ratio','Molar Volume (J/bar)',
     *            'Dependent potentials (J/mol, bar, K)',
     *            'Assemblage Index',
     *            'Modes of all phases',
     *            'Sound velocity T derivative (km/s/K)',
     *            'P-wave velocity T derivative (km/s/K)',
     *            'S-wave velocity T derivative (km/s/K)',
     *            'Adiabatic bulk modulus T derivative (bar/K)',
     *            'Shear modulus T derivative (bar/K)',
     *            'Sound velocity P derivative (km/s/bar)',
     *            'P-wave velocity P derivative (km/s/bar)',
     *            'S-wave velocity P derivative (km/s/bar)',
     *            'Adiabatic bulk modulus P derivative (unitless)',
     *            'Shear modulus P derivative (unitless)',
     *            'All phase &/or system properties',
     *            'Absolute amount (Vol, Mol, or Wt) of a phase',
     *            'Multiple property output',
     *            'Heat capacity ratio (Cp/Cv)',
     *            'Lagged or back-calculated aqueous solute chemistry'/
c----------------------------------------------------------------------
      do i = 1, istab

         if (stol(i)) then 
c                                 doing a second run, with an 
c                                 existing solvus criterion, ask
c                                 whether to change.
            write (*,1030)
            if (.not.readyn()) stol(i) = .false.

         end if

      end do 
c                                 property counter
      iprop = 0
c                                 phase composition counter
      komp = 0
c                                 choose property
      do 

         icx = 0
         lflu = .false.

         if (iprop.eq.0) then 
            write (*,1050)
            do i = 1, kprop
               write (*,1060) i,propty(i)
            end do 
         else 
            write (*,1051)
         end if

         do 

            call rdnumb (prop(1),0d0,lop,999,.false.)
            if (lop.ne.999) exit
            write (*,'(a)') 'Select a property or enter 0 to finish...'

         end do 
  
         if (lop.lt.0.or.lop.gt.kprop) then 

            write (*,1020)
            cycle 

         else if (lop.eq.0) then 

            exit 

         end if

         if (lop.eq.7.or.lop.eq.20.or.lop.eq.37) then 
c                                 modes:
c                                 get phase name
             call rnam1 (icx,pname,2)

             phluid = .false.

             if (icx.gt.0) then
                if (fp(icx)) phluid = .true.
             else if (icx.lt.0) then 
                if (ifp(-icx).ne.0) phluid = .true.
             end if

             if (phluid) then

                lflu = .true.

             else if (gflu.and.lop.ne.37) then 
c                                 the phase isn't fluid, ask if fluid should be 
c                                 included in relative modes:
                write (*,1120) 

                if (readyn()) lflu = .true.

             end if 
c                                 write blurb about units
             if (lop.eq.7) then 
                write (*,1080)
             else if (lop.eq.37) then 
                write (*,1090)
             end if

         else if (lop.eq.40) then 

            kfl(1) = .false.

            if (iprop.gt.0) then 
c                                eject if other props already chosen:
               call warn (54,nopt(1),icx,'CHSPRP')
               cycle

            else if (.not.gflu) then 
c                                 eject if no fluid phase
               call warn (99,0d0,0,'no stable fluid phase was '//
     *                             'identified in this calculation')
               cycle

            else if (lopt(32).and.lopt(25)) then 
c                                  ask which result is to be output
               write (*,1160)

               if (readyn()) kfl(1) = .true.

            else if (.not.lopt(25)) then 
c                                 eject if no aqueous species
               call warn (99,0d0,0,'missing data for back-'//
     *                             'calculation of solute speciation')
               cycle

            else                  
c                                 back-calculated is the only option
               kfl(1) = .true. 

            end if 
c                                 identify the solvent
            icx = idaq
            kop(1) = lop
            kcx(1) = icx
             
         else if (lop.eq.25) then 
c                                eject if other props already chosen:
            if (iprop.gt.0) then 
               call warn (54,nopt(1),icx,'CHSPRP')
               cycle
            end if 
c                                 all modes
             write (*,1070)

             if (readyn()) then
c                                 warn about fancy_cumulative_modes
                if (lopt(45)) then 
                   write (*,1170)
                else 
                   write (*,1180)
                end if 

                lopt(2) = .true.

             else

                lopt(2) = .false.

             end if 
c                                 ask if fluid should be included:
             if (gflu) then 

                write (*,1120) 

                if (readyn()) lflu = .true.

             end if 
c                                 double loop necessary because solution 
c                                 i may occur as j coexisting phases
             do i = 1, istab
                do j = 1, nstab(i)

                   iprop = iprop + 1
                   call getnam (dname(iprop),idstab(i))

                end do 
             end do 
 
         else if (lop.eq.6.or.lop.eq.23) then
c                                 warn if no potentials
            if (jpot.eq.1.and.lop.eq.23) then

               call warn (31,nopt(1),iopt(1),'CHSPRP')
               cycle

            end if 
c                                 get component to be contoured
            do

                write (*,1000)

               if (lop.eq.23) then 

                  write (*,1010) (i, cname(i), i = 1, icp)

               else 

                  write (*,1010) (i, cname(i), i = 1, icomp)

               end if 

               read (*,*,iostat=ier) icx

               if (ier.eq.0) exit

               call rerr

            end do 
c                                 ask if fluids included
            if (gflu.and.lop.eq.6) then 

               write (*,1120) 

               if (readyn()) lflu = .true. 

            end if 

         else if (lop.eq.8) then
c                                 get solution identity
            do 

               call rnam1 (icx,pname,0)
               if (icx.gt.0) exit  
               write (*,1140)

            end do 
c                                 get user defined composition:
            komp = komp + 1

            if (komp.gt.m13) then 
               call warn (27,nopt(1),m13,'CHSPRP')
               cycle 
            end if 

            call mkcomp (komp,icx)

         else if (lop.eq.36.or.lop.eq.38) then   
c                                eject if other props already chosen:
            if (iprop.gt.0) then 
               call warn (54,nopt(1),icx,'CHSPRP')
               cycle
            end if 
c                                 multi-prop options, get case i: 
c                                 1 - system, 2 - phase
c                                 3 - system + phases
            write (*,1130)
            call rdnumb (nopt(1),0d0,i,1,.false.)
c                                 convert kop to the internal 
c                                 icx flag => icx = 0 system prop,
c                                 icx = 999 all props, else phase index
            if (i.eq.3) then 

               icx = 999 

            else if (i.eq.2) then 
c                                 get phase index
               call rnam1 (icx,pname,2)

            end if 

            if (gflu) then

               write (*,1120) 
               if (readyn()) lflu = .true.

            end if

            if (lop.eq.36) then 

               iprop = i8 + 3 + icomp + icp + 1
               if (iprop.gt.i11) call error (1,0d0,iprop,'I11')

            else 
c                                 lop.eq.38
               kop(1) = 38
c                                 custom list
               do

                  write (*,1150)
                  read (*,*,iostat=ier) i

                  if (ier.ne.0.or.i.gt.kprop-1.or.i.lt.0) then
                     write (*,1020)
                     cycle
                  else if (i.eq.8.or.i.eq.6.or.i.eq.23.or.i.eq.24.or.
     *                     i.eq.25.or.i.eq.36.or.i.eq.38) then 
                     write (*,1100) 
                     cycle
                  else if (i.eq.0) then 
                     exit
                  end if 
c                                 save property choice
                  iprop = iprop + 1
                  if (iprop.gt.i11) call error (1,0d0,iprop,'I11')
                  kop(iprop+1) = i                     

               end do

            end if 

         else if (lop.ne.6.and.lop.ne.8.and.lop.ne.24) then
               
            if (icx.eq.0) then   
c                                 ask if bulk or phase property
               write (*,1110)

               if (.not.readyn()) then 
c                                 it's a bulk property, ask if fluid
c                                 should be included:
                  if (gflu) then 
                     write (*,1120)
                     if (readyn()) lflu = .true. 
                  end if 

               else 
c                                 get phase name
                  do 
                     call rnam1 (icx,pname,2)
                     if (icx.lt.1.and.lop.eq.8) then
                        write (*,1140)
                        cycle 
                     end if 
                     exit 
                  end do  

               end if

            end if 

         end if 
c                                 make dependent variable names:
         if (lop.eq.25.or.lop.eq.36.or.lop.eq.38.or.lop.eq.40) then
c                                 multi prop options, only allowed as 
c                                 single choices.
            kop(1) = lop
            kcx(1) = icx
            kfl(1) = lflu
c                                 assign property names
c                                 lop = 25 -> all mode names are assigned above
            if (lop.eq.36) then 

c                                 "all" prop option
               mprop = i8 + 3
c                                 basic props
               do i = 1, mprop
                  call gtname (lop,icx,i,komp,pname)
               end do 
c                                 for all prop option make the
c                                 bulk composition and chemical 
c                                 potential variable names
               do i = mprop + 1, mprop + icomp
c                                 bulk compositions, lop = 6
                  call gtname (6,i-mprop,i,komp,pname)
               end do 

               do i = mprop + icomp + 1, i8 + 3 + icomp + icp
c                                 chemical potentials, lop = 23
                  call gtname (23,i-mprop-icomp,i,komp,pname)
               end do 

               dname(iprop) = 'nom_ox'

            else if (lop.eq.38) then 
c                                 "custom" prop option, kop
c                                 pointer is offset because 
c                                 kop(1) = 38
               if (icx.eq.999) pname = 'phase     '

               do i = 1, iprop
                  call gtname (kop(i+1),icx,i,komp,pname)
               end do 

            else if (lop.eq.40) then 

               call aqname

            end if 

            exit 

         else 
c                                 save the local choice options in the 
c                                 global arrays
            iprop = iprop + 1
            kop(iprop) = lop
            kcx(iprop) = icx
            kfl(iprop) = lflu
            k2c(iprop) = komp
c                                 make the name of the property and save
c                                 it in array dname
            call gtname (lop,icx,iprop,komp,pname)

         end if 

      end do 

1000  format (/,'Enter a component:')
1010  format (2x,i2,' - ',a5)
1020  format (/,'Invalid input, try again...',/)
1030  format (/,'Retain the compositional criteria you defined ',
     *          'earlier (y/n)?',/,'Answer yes only if you intend ',
     *          'to extract properties for the same phase.',/)
1050  format (/,'Select a property:')
1051  format (/,'Select an additional property or enter 0 to finish:')
1060  format (3x,i2,' - ',a60)
1070  format (/,'Output cumulative modes (y/n)?',/
     *         ,'(see www.perplex.ethz.ch/perplex_options.html'
     *         ,'#cumulative_modes)')
1080  format (//,'Fractions are Wt, Vol, or Mol depending on the '
     *         ,'perplex_option.dat proportions keyword.',//)
1090  format (//,'Amounts are kg, m3, or Mol per unit quantity system '
     *         ,'as specified by the',/
     *         ,'perplex_option.dat proportions keyword.',/)
1100  format (/,'Property not allowed for this option, try again...',/)
1110  format (/,'Calculate individual phase properties (y/n)?')
1120  format (/,'Include fluid in computation of aggregate ', 
     *          '(or modal) properties (y/n)?')
1130  format (/,'In this mode you may tabulate:',
     *     /,4x,'1 - properties of the system',
     *     /,4x,'2 - properties of a phase',   
     *     /,4x,'3 - properties of the system and its phases',/
     *         ,'Output for option 1 & 2 can be plotted with '
     *         ,'PSPLOT, PYWERAMI or MatLab.',/,'Output for '
     *         ,'option 3 can only be plotted with PHEMGP.',//
     *         ,'Select an option [default = 1]:')
1140  format (/,'Hey cowboy, that warn''t no solution, try again.',/)
1150  format (/,'Specify a property to be computed from the ',
     *          'list above [0 to end]')
1160  format (/,'This calculation was done with lagged speciation, but',
     *        ' back-calculated',/,'speciation calculations are also',
     *        ' enabled. Output the back-calculated',/,'results (y/n)?')
1170  format (/,'**warning ver210** the fancy_cumulative_modes option '
     *       ,'(perplex_option.dat) is T.',/,'This option may lead to '
     *       ,'inconsistent curves for coarsely sampled profiles.',/)
1180  format (/,'**warning ver210** the fancy_cumulative_modes option '
     *       ,'(perplex_option.dat) is F.',/,'Set this option to obtain'
     *       ,' more easily read plots of modal abundance.',/)
  
      end

      subroutine gtname (lop,icx,jprop,komp,pname)
c----------------------------------------------------------------
c makes the name of property iprop and saves it in dname(iprop)
c called only by chsprp, variable names as in chsprp.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical warned

      integer icx, jprop, lop, komp, l2p(39)

      character prname(45)*14, pname*10, temp*20

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      character cname*5
      common/ csta4  /cname(k5)

      save warned
      data warned/.false./

      save prname,l2p

      data prname/'V,J/bar/mol   ','H,J/mol       ','Gruneisen_T   ',
c                                 4-6
     *      'Ks,bar        ','Gs,bar        ','v0,km/s       ',
c                                 7-9
     *      'vp,km/s       ','vs,km/s       ','vp/vs         ',
c                                 10-12
     *      'rho,kg/m3     ','G,J/mol       ','cp,J/K/mol    ',
c                                 13-15
     *      'alpha,1/K     ','beta,1/bar    ','S,J/K/mol     ',
c                                 16-18
     *      'n,mol         ','N,g           ','Ks_{T},bar/K  ',
c                                 19-21
     *      'Gs_{T},bar/K  ','Ks_{P}        ','Gs_P          ',
c                                 22-24
     *      'v0_{T}        ','vp_{T}        ','vs_{T}        ',
c                                 25-27
     *      'v0_{P}        ','vp_P          ','vs_{P}        ',
c                                 28
     *      'cp/cv         ',
c                                 29-31
     *      'vol,%         ','wt,%          ','mol,%         ',
c                                 32-34
     *      'h,J/m3        ','cp,J/K/m3     ','blk_comp      ',
c                                 35-37
     *      'mode          ','composition   ','s,J/K/m3      ',
c                                 38-40
     *      's,J/K/kg      ','h,J/kg        ','cp,J/K/kg     ',
c                                 41-43
     *      'specific_mass ','poisson_ratio ','chemical_pot  ',
c                                 44-45
     *      'assemblage_i  ','extent        '/
c                                 l2p points from lop to prname, 
c                                 1-10
      data l2p/35,10,33,13,14,34,38,39, 3, 4,
c                                 11-20            
     *         5 , 6, 7, 8, 9,37,38,39,40,41,
c                                 21-30 
     *         42, 1,43,44, 0,22,23,24,18,19,
c                                 31-39
     *         25,26,27,20,21, 0, 45, 0, 28/
c----------------------------------------------------------------------
      if ((lop.eq.6.or.lop.eq.36).and..not.warned) then
c                                this warning could be shifted to where
c                                the phase/system option is chosen.
         if (lopt(23)) then

            if (iopt(2).eq.0) write (*,1010) 'mole','mass','mole'

         else

            if (iopt(2).eq.1) write (*,1010) 'mass','mole','mass'

         end if

         warned = .true.

      end if
c                                 make property name
      if (lop.eq.6) then
c                                 wt% or mol component icx
         if (lopt(23)) then

            if (lopt(41)) then
c                                 absolute mass, units hardwired here, but
c                                 are actually determined by component formula 
c                                 mass specified in the thermodynamic data file.
               temp = cname(icx)//',g       '

            else
c                                 mass fraction (%)
               temp = cname(icx)//',wt%     '

            end if

         else

            if (lopt(41)) then

               temp = cname(icx)//',mol,abs '

            else

               temp = cname(icx)//',mol,pfu '

            end if 

         end if

      else if (lop.eq.7) then
c                                 mode of a phase
         if (iopt(3).eq.0) then 
c                                 vol%
            temp = pname//',vo%'

         else if (iopt(3).eq.1) then 
c                                 wt%
            temp = pname//',wt%'

         else  
c                                 mol%
            temp = pname//',mo%'

         end if 

      else if (lop.eq.8) then 
c                                phase composition
        write (temp,'(a,i2,a)') 'C['//pname,komp,']'

      else if (lop.eq.23) then 
c                                chemical potential of a component
        temp = 'mu['//cname(icx)//'],J/mol'

      else if (lop.eq.36) then 
c                                allprp option, lop points directly
c                                to prmame
         temp = prname(jprop)

      else if (lop.eq.37) then
c                                extent of a phase
         if (iopt(3).eq.0) then 
c                                volume
            temp = pname//',m3 '

         else if (iopt(3).eq.1) then 
c                                 mass
            temp = pname//',kg '

         else  
c                                 mol
            temp = pname//',mol'

         end if   

      else if (lop.eq.38) then 
c                                 custom prop (phemgp)
         temp = prname(l2p(jprop))

      else 

         temp = prname(l2p(lop))
          
      end if    

      call unblnk(temp)

      read (temp,'(a14)') dname(jprop)

      if (lop.eq.8) write (*,1000) dname(jprop)

1000  format (/,'This composition will be designated: ',a,/)
1010  format (/,'**warning ver069** composition_phase is ',a,' but '
     *       ,'composition_system is ',a,'. As ',/, 
     *        'only one unit can be output for property choice 36 '
     *        'system units will be',/,'indicated; the ',
     *        'true units for phase compositions are ',a,/)

      end

      subroutine finprp (dim,n5name,n6name,node)
c----------------------------------------------------------------
c wrap up property output, writes blurb on processing.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer dim, i 

      logical node

      character*100 n5name,n6name

      integer inv
      character dname*14, titl1*162
      common/ cst76 /inv(i11),dname(i11),titl1

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname
c----------------------------------------------------------------------
c                                 write data ranges
      write (*,1090) nopt(7)
      write (*,'(5x,200(a,1x))') (dname(i),i=1,iprop)
      write (*,'(a3,2x,200(g14.6E3,1x))') 'min',(prmn(i),i=1,iprop)
      write (*,'(a3,2x,200(g14.6E3,1x))') 'max',(prmx(i),i=1,iprop)      

      if (kop(1).eq.25) then 
c                                 create plt format file for "all
c                                 modes" output option
         call outmod (dim,n6name,node)

         if (dim.eq.1) then 
            write (*,1000) n6name, n5name
            write (*,1030) 
            write (*,1010) dim,'tab'
            write (*,1020) 
         else 
            write (*,1040) dim,'tab',n5name
            write (*,1010) dim,'tab'
            write (*,1060)
         end if

      else if (kcx(1).eq.999) then 

         write (*,1040) dim,'phm',n5name
         write (*,1010) dim,'phm'
         if (dim.eq.1) then 
            write (*,1070)
         else 
            write (*,1080)
         end if

      else 

         write (*,1040) dim,'tab',n5name
         write (*,1010) dim,'tab'
        
         if (dim.eq.1) then 
            write (*,1020)
         else 
            write (*,1060) 
         end if 

      end if 

      close (n5)

1000  format (/,'Output has been written to two files:',//,
     *       5x,'plt format is in file: ',a,/,
     *       5x,'1d tab format is in file: ',a)
1010  format (/,i1,'d ',a,' format files can be processed with:',/)
1020  format (5x,'PSTABLE - a Perple_X plotting program',
     *      /,5x,'PERPLE_X_PLOT - a Matlab plotting script',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on tab format refer to:',/,5x,
     *       'perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format',
     *       '.txt',/)
1030  format (/,'plt format files can be plotted with:',//,
     *       5x,'PSVDRAW')
1040  format (/,'Output has been written to the ',i1,
     *          'd ',a,' format file: ',a)
1060  format (5x,'PSTABLE - a Perple_X plotting program',
     *      /,5x,'PERPLE_X_PLOT - a MATLAB plotting script',
     *      /,5x,'PYWERAMI - github.com/ondrolexa/pywerami',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on tab format refer to:',
     *      /,5x,'perplex.ethz.ch/perplex/faq/Perple_X_tab_file_format',
     *       '.txt',/)
1070  format (5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on phm format refer to:',
     *      /,5x,'perplex.ethz.ch/perplex/faq/Perple_X_phm_file_format',
     *           '.txt',/)
1080  format (5x,'PHEMGP - perplex.ethz.ch/phemgp',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on phm format refer to:',
     *      /,5x,'perplex.ethz.ch/perplex/faq/Perple_X_phm_file_format',
     *           '.txt',/)
1090  format (/,'Data ranges excluding values equal to bad_number ',
     *       '(',g10.3,') specified in perplex_option.dat:',/)
      end  

      subroutine aqname 
c----------------------------------------------------------------
c makes a list of property names for and saves then in dname(iprop)
c called only by chsprp for back-calculated and lagged aqueous speciation
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      character spec(9)*14

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      character cname*5
      common/ csta4  /cname(k5)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      save spec 
      data spec/'pH-pH_0','pH','error_pH','permittivity','I,m',
     *          'ref_chg','tot_solute_m','tot_molality','solvent_kgfw'/
c----------------------------------------------------------------------
      if (icomp+ns+aqct+9.gt.i11)
     *                           call error (1,rt,icomp+ns+aqct+9,'i11')
c                                 bulk composition, wt% or mol 
      do i = 1, icomp 

         if (iopt(2).eq.1) then

            if (lopt(41)) then
c                                 absolute mass, units hardwired here, but
c                                 are actually determined by component formula 
c                                 mass specified in the thermodynamic data file.
               dname(i) = cname(i)//',g       '

            else
c                                 mass fraction (%)
               dname(i) = cname(i)//',wt%     '

            end if

         else

            if (lopt(41)) then

               dname(i) = cname(i)//',mol,abs '

            else 

               dname(i) = cname(i)//',mol,pfu '

            end if

         end if

         call unblnk(dname(i))  

      end do

      iprop = icomp 
c                                 solvent composition, mol or molal
      do i = 1, ns 

         iprop = iprop + 1

         if (lopt(26)) then
c                                 mole fraction
            dname(iprop) = 'y_{'//names(jnd(i))//'}'
         else 
c                                 molality
            dname(iprop) = 'm_{'//names(jnd(i))//'}'
         end if  

         call unblnk(dname(iprop))   

       end do    
c                                 solute composition, mol or molal
      do i = 1, aqct 

         iprop = iprop + 1

         if (lopt(27)) then
c                                 mole fraction
            dname(iprop) = 'm_{'//aqnam(i)//'}'
         else 
c                                 molality
            dname(iprop) = 'y_{'//aqnam(i)//'}'
         end if  

         call unblnk(dname(iprop))   

      end do 
c                                 special variables
      do i = 1, 9 
         iprop = iprop + 1
         dname(iprop) = spec(i)
      end do 

      end

      subroutine lagprp (jd)
c-----------------------------------------------------------------------
c load interpolated lagged aqueous speciation results into the prop array 
c for WERAMI 

c   jd - is the pointer to the phase in the local assemblage if only 
c        one phase is to be loaded
c javg - counter for the jdsol pointer array if multiple phases are 
c        to be averaged. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, k, l, ld, jd

      double precision ntot, x, ximp

      integer javg,jdsol
      common/ cxt5 /javg,jdsol(k5)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)
c----------------------------------------------------------------------
      if (jd.eq.0) then 
c                                 no solvent phase stable
         do i = 1, iprop
            prop(i) = nopt(7)
         end do 
 
      else if (javg.le.1) then 
c                                 only one phase:
c                                 bulk composition
         if (lopt(41)) then
c                                 absolute composition
            x = props(16,jd)
         else
c                                 relative composition
            x = 1d0
         end if 

         do i = 1, icomp
            prop(i) = pcomp(i,jd) * x
         end do

         k = icomp

         do i = 1, ns 
c                                 solvent speciation
            k = k + 1

            if (lopt(26)) then
c                                 mole fraction
               prop(k) = caq(jd,i)
            else 
c                                 molality
               prop(k) = caq(jd,i)*caq(jd,na2)
            end if  

         end do

         do i = sn1, nsa 
c                                 solute speciation
            k = k + 1

            if (lopt(27)) then
c                                 molality
               prop(k) = caq(jd,i)
            else 
c                                 mole fraction
               prop(k) = caq(jd,i)/caq(jd,na2)
            end if  

         end do 
c                                  other properties:
c                                 pH - pH0
         prop(k+1) = caq(jd,na3+3)
c                                 pH
         prop(k+2) = caq(jd,na3+2)
c                                 err_log10(K_w)
         prop(k+3) = caq(jd,na3+1)
c                                 dielectric constant
         prop(k+4) = caq(jd,nat)
c                                 ionic strength
         prop(k+5) = caq(jd,na1)
c                                 ref_chg
         prop(k+6) = caq(jd,nat-1)
c                                 total solute molality
         prop(k+7) = caq(jd,nat-2)
c                                 total molality
         prop(k+8) = caq(jd,na2)
c                                 solvent molar mass, kg/mol-solvent-species
         prop(k+9) = caq(jd,na3)

      else 
c                                 averaging multiple phases:

c                                 bulk composition has already been 
c                                 averaged by avgcmp, load into prop:
c                                 bulk composition
         if (lopt(41)) then
c                                 absolute composition
            x = props(16,jd)
         else
c                                 relative composition
            x = 1d0
         end if 

         do i = 1, icomp
            prop(i) = pcomp(i,jd)*x
         end do
c                                 initialize remaining prop values
         do i = icomp+1, icomp+nat
            prop(i) = 0d0
         end do
c                                 get molar proportions:
         ntot = 0d0
         ximp = 0d0

         do i = 1, javg
            ntot = ntot + props(16,jdsol(i))
         end do 

         do l = 1, javg

            ld = jdsol(l)
            x = props(16,ld)/ntot
            if (caq(ld,na1).ne.0d0) ximp = ximp + x

            k = icp

            do i = 1, ns 
c                                 solvent speciation
               k = k + 1

               if (lopt(26)) then
c                                 mole fraction
                  prop(k) = prop(k) + x * caq(ld,i)
               else 
c                                 molality
                  prop(k) = prop(k) + x * caq(ld,i)*caq(ld,na2)
               end if

            end do

            do i = sn1, nsa
c                                 solute speciation
               k = k + 1

               if (lopt(27)) then
c                                 molality
                  prop(k) = prop(k) + x * caq(ld,i)
               else 
c                                 mole fraction
                  prop(k) = prop(k) + x * caq(ld,i)/caq(ld,na2)
               end if

            end do 
c                                 other properties:
c                                 pH - pH0
            prop(k+1) = prop(k+1) + x * caq(ld,na3+3)
c                                 pH
            prop(k+2) = prop(k+2) + x * caq(ld,na3+2)
c                                 err_log10(K_w)
            prop(k+3) = prop(k+3) + x * caq(ld,na3+1)
c                                 dielectric constant
            prop(k+4) = prop(k+4) + x * caq(ld,nat)
c                                 ionic strength
            prop(k+5) = prop(k+5) + x * caq(ld,na1)
c                                 ref_chg
            prop(k+6) = prop(k+6) + x * caq(jd,nat-1)
c                                 solute molality
            prop(k+7) = prop(k+7) + x * caq(ld,nat-2)

         end do

         if (ximp.lt.1d0.and.ximp.gt.0d0) then
c                                renomalize err_log_kw, pH, Delta_pH, epsilon
            prop(k+1) = prop(k+1)/ximp
            prop(k+2) = prop(k+2)/ximp
            prop(k+3) = prop(k+3)/ximp
            prop(k+4) = prop(k+4)/ximp

         end if

      end if

      end
