c----------------------------------------------------------------------
c rewrite 2010 is a program to rewrite *ver.dat thermodynami data files
c from before May 2010 to the current format.

c to run this code you must temporarily modify perplex_parameters such that:

c k5 = max number of components in the data base (<=k0)
c m7 = number of parameters in a transition (reduce to value in the data 
c      base to be translated and modify initialization of tstrg in block data)

c additionally, use/suppress hsc conversion if desired. 

c----------------------------------------------------------------------

      implicit none
 
      include 'perplex_parameters.h'

      integer i,j 

      character text(5)*1

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ic
      common/ cst42 /ic(k0)

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 10

      do i = 1, k5
         ic(i) = i
      end do 
c                               assign data files
      call sopen 
c                               Read THERMODYNAMIC DATA file (N2):
c                               read the data base header
c old format top
c      call otopn2 (4)
c new format top
      call topn2 (4)

      icomp = icmpn
c                               lenght of component names
      do i = 1, icmpn
         read (cmpnt(i),'(5a1)') text
         cl(i) = 5
         do j = 1, 5
            if (text(j).eq.' ') then 
               cl(i) = j-1
               exit
            end if 
         end do 
      end do 


c                               read and echo data cards with
c                               component conversion
      call getph 
      
      end
  
      subroutine getph 
c------------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character name*8 

      logical eof
 
      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      double precision cp
      common/ cst12 /cp(k5,k1)

      do 
c                                 need to eliminate hsc conversion.
         call ogtphi (name,eof)

         if (eof) exit

         call outdat (n8,1,1)

      end do 
 
      end



      subroutine ogtphi (name,eof)
c-----------------------------------------------------------------------
c old getphi, modified for routine outdat
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i, j, ier

      double precision ct

      logical eof

      character*8 name, oname, record*1
 
      double precision emodu
      common/ cst318 /emodu(k15)

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      double precision therlm,therdi
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      integer iamir
      common/ cst53 /iamir
 
      double precision delta
      common/ cst325 /delta(11)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer eos
      common/ cst303 /eos(k10)

      character*80 com
      common/delet/com 

      save oname
 
      data oname/' '/
c----------------------------------------------------------------------
      eof = .false.
      jlam = 0 
      lmda(1) = 0
      ltyp(1) = 0
      lct(1)  = 0 
      idis(1) = 0 

      do 
 
         read (n2,'(a)',iostat=ier) record

         if (ier.lt.0) then
 
            eof = .true.
            return

         else if (ier.gt.0) then
                 
            exit

         end if 
c                          check for comments, i.e., data 
c                          with a blank 1st character 
         if (record.eq.' ') then 

            backspace (n2)

            read (n2,'(a)') com
            if (com.ne.' ') write (n8,1111) com
1111     format (5x,'| ',a)

            cycle

         end if 

         backspace (n2)

1010     format (a,i2,i2,i2,i2,a)
         read (n2,1010) names(1),ieos,ikind,ilam,idiso,com          
         read (n2,*,iostat=ier) (cp(i,1), i= 1, icmpn), 
     *                          (thermo(i,1), i= 1, 18)

         if (ier.ne.0) exit 
         if (ilam.ne.0) then
c                                 determine number of transitions from
c                                 flag ilam:
            jlam = ilam

            if (ilam.eq.20) then 
c                                 hp2010
               jlam = 1
               ilam = 13
            else if (ilam.eq.10) then 
               jlam = ilam - 9
            else if (ilam.gt.6) then
               jlam = ilam - 6
            else if (ilam.gt.3) then
               jlam = ilam - 3
            end if 
 
            lct(1)  = jlam 
            ltyp(1) = ilam/3 + 1
            lmda(1) = 1
 
            do i= 1, jlam
               read (n2,*,iostat=ier) (tm(j,i), j = 1, m7 - 2)
               if (ier.ne.0) exit 
            end do 

         end if
 
         if (idiso.ne.0) then 
            idis(1) = 1
            read (n2,*,iostat=ier) (therdi(i,1),i=1,m8)
            if (ier.ne.0) exit
         end if  

         if (ikind.ne.0) then 
            read (n2,*,iostat=ier) (emod(i,1),i=1,3)
            if (ier.ne.0) exit
         end if 
c                                 read uncertainties for MC calculations
         if (iamir.eq.999) then 
            read (n2,*,iostat=ier) delta
            if (ier.ne.0) exit
         end if 
  
         oname = name
 
         exit 

      end do 
      
      if (ier.ne.0) then 

         if (oname.ne.' ') then
            call error (23,ct,i,oname)
         else
            call error (23,ct,i,'  none  ')
         end if
      
      end if 


      if (ieos.eq.8.or.ieos.eq.9.or.ieos.eq.14) then 
         eos(1) = ieos
      else if (ieos.lt.100) then 
         eos(1) = 0 
      else 
         eos(1) = ieos
      end if 

c EoS types

c   1    - normal polynomials on V, alpha, Cp (Helgeson et al 78/Berman)
c   2    - normal polynomials on alpha, Cp; Murnaghan on pressure (b8 = K' >= 3) (Holland & Powell 98)
c   3    - normal polynomial on Cp; Ghiorso polynomial on alpha (b6=0), Birch-Murnaghan on pressure (b8 = -K' <= -3)
c   4    - normal polynomials on alpha, Cp; Birch-Murnaghan on pressure (b8 = -K' <= -3) (Saxena & Fei?)
c   5    - Debye Mie-Gruneisen (Stixrude & Bukowinski JGR '93) (V0<0, S0>0)
c   6    - Debye Mie-Gruneisen (Stixrude & Lithgow-Bertelloni GJI '05) (V0<0, S0<0) 
c   7    - normal polynomials on alpha, Cp; exponential polynomial on pressure (Haas et al. 86; Gottshalk 97)

c Special EoS (internal) types

c   601  - Komabayashi & Fei JGR 2010 
c   602  - Komabayashi & Fei JGR 2010
c   603  - Komabayashi & Fei JGR 2010
c   604  - Komabayashi & Fei JGR 2010

c mock-Lambda transition types

c   1    - UBC style (Berman JPet 1988)
c   2    - Standard (USGS 1452, JANAF)
c   3    - Helgeson et al. AJS 1978
c   4    - Holland & Powell JMet 1998

      if (eos(1).ne.8.and.eos(1).ne.9.and.eos(1).ne.14) then 
         if (thermo(3,1).lt.0d0) then 
c                              negative "volume" signals one of the 
c                              stixrude EoS's, for these EoS's "s" is
c                              +/-n - number of atoms pfu.
            if (thermo(2,1).gt.0d0) then

               eos(1) = 5
c                              stixrude & bukowinski JGR '93
c            n = s

            else

                eos(1) = 6
c                              stixrude & lithgow-bertelloni GJI '05
c           n = -s

            end if 

         else if (thermo(18,1).eq.0d0) then

            eos(1) = 1
c                              normal polynomial vdp term:
         else

            if (ieos.lt.100) then
               if (thermo(16,1).eq.0d0) then
                    eos(1) = 3
               else if (thermo(18,1).ge.3d0) then
                    eos(1) = 2
c                              convert from k = k0 + b7*T
c                              back to k = ktr + b7*(T-Tr)
                    thermo(16,1) = thermo(16,1) + 298.15*thermo(17,1)
               else if (thermo(18,1).le.3d0) then
                    eos(1) = 4
               else 
                    eos(1) = 7
               end if
            end if
         end if 
      end if 

      end


      subroutine otopn2 (jopt)
c----------------------------------------------------------------------
c topn2 reads the header of the thermodynamic data file, if jopt = 1
c then data base choice is known (idbase), else if 1 asks console for a 
c choice else if > 3 echos data except components and atwts to n8
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character*5 tag*11, n2name*100, string*140

      integer icod, ig(3), jopt, i

      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision p,t,xco2,u1,u2,tr,pr,r,ps,v(l2)
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      integer io3,io4,io9
      common/ cst41 /io3,io4,io9

      double precision ctrans
      integer ictr,itrans
      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
  
      common/ cst41a /n2name

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      double precision atwt
      common/ cst45 /atwt(k0)

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec
c-----------------------------------------------------------------------
      rewind n2
c                               read the number of data bases
c                               obsolete since always 1:
      read (n2,*,err=90) i
c                               read the extrinsic variable names:
      read (n2,1000,err=90) (vname(i), i = 1, 3)
      read (n2,1010,err=90) (ig(i), i = 1, 3)

      if (jopt.lt.4) then 
         if ((ifug.ge.10.and.ifug.le.12).or.
     *       ifug.eq.15.or.ifug.eq.17.or.ifug.eq.18) then 
            vname(3) = ' X(O) ' 
         else if (ifug.eq.25) then 
            vname(3) = 'Y(CO2)*'
         else if (ifug.eq.13) then 
            vname(3) = 'X(H2)'
         end if 
      end if  
c                               read delt the finite difference
c                               increments for v1, v2, and v3. and dtol
c                               utol and ptol, the critical tolerances
c                               (in energy units) for determination of
c                               the stability of divariant, univariant
c                               and invariant equilibria or reactions.
      read (n2,*,err=90) delt, dtol
c                               read obsolete data base code and title
c                               and the data base reference state
c                               conditions consistent with v1 and v2.
      read (n2,*,err=90) icod,pr,tr
      read (n2,'(a)',err=90) dname   
c                               read number of data base comps 
      read (n2,*,err=90) icmpn
c                               read component names.
      read (n2,1180,err=90) (cmpnt(i), i = 1, icmpn)
c                               component atomic wts. 
      read (n2,*,err=90) (atwt(i), i = 1, icmpn)
c                               read pointers for water and co2 cmpnt.  
c                               if data file doesn't contain water
c                               or co2 dummy values
c                               must be included in the file (ne. 0).
      read (n2,*,end=90) idspe(1), idspe(2)

      v(1) = pr
      v(2) = tr
      v(3) = 0d0
      v(4) = 0d0
      v(5) = 0d0
      vname(4) = 'mu(C1)'
      vname(5) = 'mu(C2)'


      if (jopt.gt.3) then 
         write (n8,1020) 
1020  format (/,' | comments are indicated by the | character.',/,
     *     ' | check for warnings at the end of the header section.',/) 
c                                 echo formatted header data for ctransf/actcor:
         write (n8,'(a,a,/)') dname,' |<= data base title'
         write (n8,'(a,a)') 'begin_standard_variables |<= name (<9 ',
     *                      'characters), reference value, tolerance'
         do i = 1, l2 
            write (n8,'(a8,1x,f7.2,3x,g6.1E1)') vname(i),v(i),delt(i)
         end do 
         write (n8,'(a,/)') 'end_standard_variables'
         write (n8,'(a,g6.1E1,a,/)') 'tolerance  ',dtol,
     *         '  |<= DTOL for unconstrained minimization, energy units'
         write (n8,'(a,a)') 'begin_components |<= name (<5 characters),'
     *                      ,' molar weight (g)'
         write (n8,'(a5,1x,f9.4)') (cmpnt(i),atwt(i), i = 1, icmpn)
         write (n8,'(a,/)') 'end_components'
c         if (lopt(7)) then 
            write (n8,'(a)') 'begin_special_components'
            write (n8,'(a)') 'H2O'
            write (n8,'(a)') 'CO2'
            write (n8,'(a,/)') 'end_special_components'
c         end if 
      end if 
c                                 read and echo unformatted comments and make data                            
      do 

         read (n2,'(a)',end=90) string
         read (string,'(a)') tag
         if (jopt.gt.3) write (n8,'(a)') string

         if (tag.eq.'begin_makes') then
 
            call rmakes (jopt)

            cycle 

         else if (tag.ne.'end'.and.tag.ne.'END') then     

            cycle

         else 

            exit       

         end if

      end do  

      goto 99

90    call error (21,r,i,n2name)

1000  format (3(a8,18x))
1010  format (i2)
1030  format (6(g12.6,1x))
1180  format (6(a5,1X)/6(a5,1X)/6(a5,1X))
9010  format (10(i2,1x))
9020  format (8(g9.2,1x))
9030  format (i2,1x,8(g12.6,1x))

99    t = tr
      p = pr

      end 
