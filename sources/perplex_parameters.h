

      integer h4,h5,h6,h8,h9,h0
      integer i6,i7,i8,i9,i10,i11
      integer j3,j4,j5,j6,j9
      integer k0,k1,k2,k3,k4,k5,k7,k8,k9,k10,k13,k14,k15
      integer k16,k17,k18,k19,k20,k21,k22,k23,k24,kd2,k25
      integer l2,l3,l5,l6,l7,l8,l9,l10,lchar
      integer m0,m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m14,m15
      integer m16,m17,m18,m19,m20,m21,m22,m23,m24,m25
      integer msp,mst,mdim,ms1
      integer n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,nsp,nx,ny
      integer memory,k31,k32
!                                 n0  - starting LUN-1 for fractionation files, these files may 
!                                       have LUNs from n0+1 up to n0+k23
!                                 n1  - problem definition file.
!                                 n2  - thermodynamic data file.
!                                 n3  - print output file
!                                 n4  - plot output file 1
!                                 n5  - plot output file 2 (bplot)
!                                 n6  - reaction list file
!                                 n7  - unused
!                                 n8  - locally opened and closed files
!                                 n9  - solution model file 
!                                 n10 - autorefine file 1
!                                 n11 - autorefine file 2
!                                 n12 - failed result file (fld)
!                                 n13 - aq error file (pts)
      parameter (n10=7,n11=8,n1=11,n2=12,n3=13,n4=14,n5=15,n6=16,n7=17)
      parameter (n8=18,n9=19,n12=20,n13=21,n0=30)
!                                 msp - max number of species on a solution identisite
!                                 mst - max number of distinct identisites per solution
!                                 mdim - hard constraint on max number of dimensions
!                                        for a solution model composition space.
      parameter (mst=4,mdim=8,msp=mdim+6,ms1=msp-1)
!                                 h4  - max-number of polytopes in a composite composition space
!                                 h5  - max number of saturated components
!                                 h6  - max number of saturated composants in any subcomposition
!                                 h8  - max number of excluded phases
!                                 h9  - max number of solutions
!                                 h0 - h9 + 1, added to eliminate temporary solution model arrays 
      parameter (h4=5,h5=5,h6=500,h8=250,h9=30,h0=h9+1)
!                                 i6  - maximum number of independent chemical potentials (or 
!                                       fugacity/activities).
!                                 i7  - number of system props used in werami
!                                 i8  - number of properties saved for each phase and
!                                       for the bulk composite, less than i8 properties
!                                       may be read from the bulk property file. 
!                                       if i8 is changed then adjust:
!                                          pt2prp - frendly.f
!                                          l2p    - werami.f
!                                          prname - werami.f
!                                 i9  - max number of solutions in solution model file
!                                 i10 - max number of option values in perplex_option.dat
!                                 i11 - max number of dependent properties in a tab format file
      parameter (i6=2,i7=20,i8=28,i9=200,i10=70,i11=150)
!                                 j3  - max number of ordered species
!                                 j4  - max number of species in the definition of a dependent species
!                                 j5  - max number of stoichiometric limits on an ordered species
!                                 j6  - max number of terms in a stoichiometric limit on an ordered species
!                                 j9  - max number of divariant assemblages
      parameter (j3=4,j4=8,j5=8,j6=12,j9=160000)
!                                 k0  - max number of database components
!                                 k1  - max number of compounds
!                                 k2  - max number of invariant and univariant compound 
!                                       equilibria and max number of divariant compound
!                                       assemblages for constrained bulk composition
!                                       calculations
!                                 k3  - max number of distinct phase (as opposed to 
!                                       pseudocompound assemblages)
!                                 k4  - number of parameters in a Perple_X EoS
!                                 k5  - max total number of active components
!                                 k9  - maximum number of true compounds with lambda transitions
!                                 k10 - max number of true compounds
!                                 k13 - pseudocompound parameter array dimension
!                                 k14 - number of parameters in the data base EoS
!                                 k15 - number of elastic moduli parameters in the emod array
!                                 k16 - max number of make definitions
!                                 k17 - max number of entities in a make definition
!                                 k18 - max number of simplicial coordinates saved for static compositions.
!                                 k19 - max number of refinement points for adaptive optimization.
!                                 k20 - max number of "z" coordinates saved for pseudocompounds
!                                       generated by adaptive refinement, max value = k21*m4, 
!                                       but usually much smaller (see k18).
!                                 k21 - max number of dynamic compositions for adaptive refinement.
!                                 k23 - max number of phases to be fractionated.
!                                 k24 - max number of static simplicial coordinate sets (jcox).
!                                 k25 - max number of dynamic simplicial coordinate sets (jcoz).
c----------------------------------------------------------------------
!                                 The array dimensions for static compositions:
!                                    k1, k18, k24
!                                 bear the same relation to eachother as the 
!                                 the array dimensions for dynamic compositions:
!                                    k21, k20, k25.
!                                 Thus if k18 and k24 are expressed as multiples of k1
!                                    k18 = k1 * k31
!                                    k24 = k1 * k32.
!                                 It is to be expected that the required values of k20 and k25 will be
!                                    k20 = k21 * k31
!                                    k25 = k21 * k32.
!                                 An additional array (simp) is used for temporary storage during both dynamic
!                                 and static subdivision, its dimension k13 should be proportional to, but less than,
!                                 max(k1,k21) (i.e., it is expected that k33 >= 1)
!                                    k13 = imax(k1,k21)/k33.
!                                 The memory allocated for compositions is then 
!                                    memory = (1 + k31 + k32)*(k1 + k21) + imax(k1,k21)/k33
!                                 or if k1 is taken as the independent variable 
!                                    k21 = ((memory (-k31 - k32 - 1)*k1) - imax(k1,k21)/ k33) / (1 + k31 + k32).
!                                 Experience indicates that for most problems k1>k21, and that k21 approaches
!                                 k1 for complex (10-d) solutions with the default values of resolution_factor (2)
!                                 reach_increment (0). For larger resolution_factor or reach_increment the required
!                                 value for k21 may excede k1, whereas for problems with simple (low-dimension) 
!                                 solutions k21 may be << k1. The compromise adopted here is to equate k21 and k1
!                                 in which case 
!                                    k1 = memory / (2 * k31 + 2 * k32 + 2 + 1 / k33)
!                                 The value of memory is dependent on the other parameters set here as well as 
!                                 system/compiler limitations, which typically limit image size to 2 Gb. For the
!                                 present parameter choices, memory was found by trial and error (i.e., by varying 
!                                 increasing the value of memory until the compiler or system complained about 
!                                 image size) to be 0.78 Gb.

!                                 The program kays.f writes out the parameters k1,k21,k20,k18,k24,k25,k13 computed
!                                 by the above logic. 

!                                 uncomment this line to specify k1 independently of k21 (assumes k1>k21)
!     parameter (memory=78000000,k31=4,k32=10,k33=1,k1=4000000)
!                                 the next two lines make k1 = k21, comment these if the line above is uncommented.
!     parameter (memory=78000000,k31=4,k32=10,k33=1)
!     parameter (k1 = memory/(2*k31+2*k32+2+1/k33))
!                                 set k13 = k21/k33 if k21 > k1
!     parameter (k13=k1/k33)
!     parameter (k21=((memory-(k31+k32+1)*k1)-k13)/(1+k31+k32))
!     parameter (k18=k1*k31, k20=k21*k31, k24=k1*k32, k25=k21*k32)
!                                 laggit version static: 
!                                    k18 = k1 * k31
!                                    k24 = k1 * k32
!                                 laggit version dynamic: 
!                                    k20 = k18 (static and dynamic generate comparable simplicial coordinates)
!                                    k25 = k21 * k32
!                                  then 
!                                    memory = k1 + k21 + k18 + k20 + k24 + k25 + k1;
!                                  and solving for k21 
!                                    k21 = (memory - k1*(2*k31 + k32 + 2))/(k32+1)
! DEBUG691
      parameter(memory=70000000,k31=2,k32=10,k1=3000000)
!     parameter(memory=42000000,k31=2,k32=10,k1=1800000)
!                                  static
      parameter(k18=k1*k31,k24=k1*k32,k13=k1)
!                                  dynamic
      parameter(k21=(memory-(2*k31+k32+2)*k1)/(1+k32))
      parameter(k20=k18,k25=k21*k32)
c----------------------------------------------------------------------
      parameter (k0=25,k2=100000,k3=2000,k4=32,k5=14)
      parameter (k7=k5+1,k8=k5+2,k17=8,k19=3*k5)
      parameter (k9=50,k10=500,k14=18,k15=6,k16=150)
      parameter (k22=mdim*mst*h4*k19,k23=25)
!                                 l2 - max number of independent potential variables
!                                 l3 - max number of variables for gridded min and graphics (l2+2)
!                                 l5 - max number of coordinates along a univariant curve
!                                 l6 - OBSOLETE 691+ (replaced by LP_max_it)
!                                      max number of iterations in lp optimization, in 
!                                      theory this may be up to ca 5*(k1+k5), generally
!                                      convergence occurs in less than 100 iterations.
!                                 l7 - max number of grid points along an axis for
!                                      constrained minimization on a 2-d grid. this
!                                      array is not essential.  
!                                 l8 - max number of levels for multilevel grids    
!                                 l9 - max number of aqueous solute species in minimization programs.         
!                                l10 - max number of parameters stored in caq for each phase.
!                                nsp - max number of species in fluid speciation routines 

      parameter (l2=5,l3=l2+2,l5=1000,l6=500,l7=2048,l8=10,l9=150,
     *           nsp=17,l10=nsp+l9+4)
!                                 m0 - max number of terms for a species site fraction?
!                                 m1 - max number of terms in excess function
!                                 m2 - max order of term in excess function
!                                 m3 - number of parameters for excess function coefficients
!                                 m4 - max number of endmembers in a solution model,
!                                      should include dependent endmembers? this cannot
!                                      exceed msp**mst+1, but will generally be MUCH
!                                      smaller.
!                                 m6 - max number of "lambda" transitions per true compound
!                                 m7 - number of "lambda" parameters per transition
!                                 m8 - maximum number of parameters to describe T dependent ordering.
!                                 m9 - maximum number of true compounds with T-dependent ordering.
!                                m10 - maximum number of mixing sites in Sconf model
!                                m11 - maximum number of species + 1 in Sconf model
!                                m12 - maximum order of terms in Sconf model
!                                m13 - 
!                                m14 - maximum number of independent endmembers
!                                m15 - maximum number of dependent endmembers
!                                m16 - max number of parameters in a redkich-kistler L
!                                m17 - max order of redlich-kistler expansion
!                                m18 - max number of pairwise terms in a redlich-kistler expansion
!                                m19 - max number of independent endmember fractions m14 - 1
!                                m20 - max number of site fractions, m10*(m11-1)
!                                m21 - m20 + m19, NLP constraints
!                                m22 - 3*m19 + m20, NLP integer workspace
!                                m23 - 2*m19**2 + 20*m19 + 11*m20, NLP real workspace
!                                m24 - max number of dynamic compositions to be saved for use as 
!                                      static compositions during auto-refine stage 
!                                m25 - max number of dynamic compositional coordinates to be saved
!                                      for use as static compositions during auto-refine stage.
      parameter (m0=12, m1=80, m2=8, m3=3, m4=96, m6=6, m7=15, m8=9,
     *           m9=10,
     *           m10=6, m11=k5, m12=4, m14=14, m15=85, m16=6, m17=5,
     *           m18=6, m19=m14, m20=m10*(m11-2)+1, m21=m20+m19,
     *           m22=3*m19+m20, m23=2*m19**2+20*m19+11*m20, 
     *           m24 = 10*60*60*k5, m25 = m24*m14)
!                                 nx - number of x-grid nodes in a contour data grid
!                                 ny - number of y-grid modes in a contour data grid
!                                 NOTE: pstable requires that parameter L5 > max(NX,NY), i.e., if
!                                 NX or NY is increased, then it may be necessary to increase L5.
      parameter (nx=1000,ny=1000)
!                                 lchar - maximum length of character strings
      parameter (lchar=400)
!                                 frac2d parameters:
      integer maxbox,lay,mpol,mord
      parameter (maxbox=1760,lay=6,mpol=7,mord=mpol-1) 


! NOTE: increasing parameter K5 requires changes to the following
! format statements:

!                     1150 in routine input1
!                     4000 in main in program build


!                                 the following are dependent parameters (also k7,k8)
      parameter (kd2=k8*35)
c------------------------------------------------------------------------------
c                                 commons with globally consistent, non-conflicting, variable names:

      integer spx, icox, jcox
      double precision xco
      common/ cxt10 /xco(k18),spx(h4,mst),icox(k1),jcox(k24)

      double precision zco
      integer icoz, jcoz, jkp, zcoct
      common/ cxt13 /zco(k20),icoz(k21),jkp(k21),jcoz(k25), zcoct

      double precision simp
      common/ cxt86 /simp(k13)

      integer scoct, snp, sco, icoct, npol
      common/ junk0 /scoct, snp(mst), sco(k13), icoct, npol(h4)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)
c                                 temporary subdivision limits:
      double precision pxmn, pxmx, pxnc
      common/ cxt108 /pxmn(h4,mst,msp),pxmx(h4,mst,msp),pxnc(h4,mst,msp)

      double precision pwt
      common/ cxt44 /pwt(h4)
c                                 interim storage array
      integer lcoor, lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)

      integer jpoint, jiinc
      common/ cxt60 /jpoint,jiinc

      integer jbulk, kbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk,kbulk

      double precision ctot
      common/ cst3  /ctot(k1)
c                                 precision stuff used in lpnag 

c                                 outprt is a universal flag
c                                 that suppresses print output in exploratory
c                                 stage of two-stage auto-refine calculation.
      logical outprt
      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p,outprt

c                                 -------------------------------
c                                 local solution model variables:
c                                 -------------------------------

      integer isimp, ipvert, ivert, pimd
      common/ cst688 /isimp(h4),ipvert(h4),ivert(h4,mst),
     *                pimd(h4,mst,msp)
c                                 solution model counter, temp logical flag
      logical ltemp1
      integer isoct
      common/ cst79 /isoct, ltemp1
c                                 -------------------------------
c                                 global solution model variables:
c                                 -------------------------------
c lstot(i)   - number of independent disordered endmembers of solution i (kstot)
c mstot(i)   - number of vertices for the composition space of solution i (istot)
c nstot(i)   - number of independent endmembers of solution i
c ndep(i)    - number of dependent endmembers of solution i = mstot(i) - lstot(i)
c nord(i)    - number of ordered endmembers of solution i = nstot(i) - lstot(i)
c poly(i)    - number of polytopes for solution i, i = h0 local value for input.
c msite(i)   - number of sites on which mixing takes place
c zsp(i, msite)    - number of indepenedntly variable species on msite.
c zsp1(i, msite)   - total number of species on msite(zsp1 = zsp for temkin models)
c zmult(i, tsite)  - effective site multiplicity*R used for configurational entropy calculation
c tzmult(i, tsite) - true site mutliplicity
c isimp(poly(h0))  - number of simplices in each polytope
c ivert(poly(h0),isimp(poly(h0)) - number of vertices in each simplex
c ipvert(poly(h0)) - number of vertices in each polytope
c istot - total number of vertices
c jmsol(m4,1:isimp(poly(h0)) - pointer from the endmember m4 to its polytope vertex
c jpoint - index of the last compound endmember in the icp x jphct optimization matrix
c jiinc - iphct - jphct, increment between the icp x iphct data matrix index and the optimization matrix
c icox(phct) - locates static compositional coordinates for composition phct->gcind
c jcox(gcind + 1:nsimp)  - locates the compositional coordinates for each simplex in xco
c icoz(phct) - locates dynamic compositional coordinates for composition phct->gcind
c jcoz(gcind + 1:nsimp)  - locates the compositional coordinates for each simplex in zco
c stind(ii) - starting index - 1 of the simplex indices of polytope ii in the sco array
c npoly(ii) - number of compositions in polytope ii

c jmsol(i, j) - species on the simplex j of endmember / vertex i
c kdsol(i) - identifier / status flag of endmember / vertex i, 0 if missing, > 0 if independent
c                 endmember, -1 or 0 if ordered(non - vertex), -2 if dependent ? , -3 to be killed.
c dedpol(ii) - true if polytope ii has no valid endmembers
c pvptr(ii, 1:2) - beginning and ending indexes of polytope ii

      integer tstot,lstot,mstot,nstot,ndep,nord
      common/ cxt25 /tstot(h9),lstot(h9),mstot(h9),nstot(h9),ndep(h9),
     *               nord(h9)
c                                 -------------------------------
c                                 model type
      logical lorder, lexces, llaar, lrecip, specil, simple, deriv
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9),
     *               specil(h9),simple(h9),deriv(h9)

      logical stable, limit, noder, lorch, boundd
      integer badinv
      double precision xlo, xhi
      common / cxt11 / xlo(m4, mst, h4, h9), xhi(m4, mst, h4, h9),
     *       badinv(h9,2), stable(h9), limit(h9), noder(h9), lorch(h9),
     *       boundd(h9)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h4,h9)
c                                polytope composition variable names are in poname(id,ii,j,k)
c                                polytope names are in poname(id,ipoly+1,1,1..ipoly)
      character poname*10
      common/ cxt47 /poname(h0,h4,mst,msp)
c                                site species names are in znames(id,1:nsite,1:nsp), id = h0 on input
c                                site names are in znames(id,1:nsite, 0)
      character znames*3
      common/ cxt48 /znames(h0,m10,0:m11)

      double precision dcoef, scoef
      common/ cxt1r /dcoef(0:m0,m11,m10,h9),scoef(m4,h9)

      integer msite, zsp
      double precision zcoef, zmult
      common/ cxt1n /zcoef(0:m0,m11,m10,h0),zmult(h0,m10),
     *               msite(h0),zsp(h0,m10)

      character zuffix*60
      integer zsp1
      logical zform
      double precision tzmult
      common/ cxt1m /tzmult(h0,m10),zsp1(h0,m10),zform(h0),zuffix(h0)

      logical quack
      integer solc, isolc
      common/ cxt1 /solc(k5),isolc,quack(k21)

      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)

      integer ksmod, kmsol, knsp
      common/ cxt0  /ksmod(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer istg, ispg, imdg, poly, pvert, pop1, nsum
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg,
     *                 xmnh, xmxh
      common/ cxt6r /
     *      xmng(h9,h4,mst,msp),xmxg(h9,h4,mst,msp),xncg(h9,h4,mst,msp),
     *      xmno(h9,h4,mst,msp),xmxo(h9,h4,mst,msp),reachg(h9),
     *      xmnh(h9,h4,mst,msp),xmxh(h9,h4,mst,msp)

      common/ cxt6i /istg(h9,h4),ispg(h9,h4,mst),pop1(h9),nsum(h9),
     *               imdg(ms1,mst,h4,h9),poly(h0),pvert(h9,h4,2)

      logical bdx, nrf
      integer ldsol
      common/ cxt36 /ldsol(m4,h9),bdx(h9),nrf(h9)
c                                 -------------------------------
c                                  variables set from perplex_option.dat
      integer iopt
      logical lopt
      character valu*3
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10),valu(i10)
c                                 -------------------------------
c                                 local solution model variables:
      logical stck, norf, lowrch, badx
      integer xtyp
      double precision reach
      common/ cxt61 /reach,xtyp,stck,norf,lowrch,badx

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      logical dedpol
      integer jmsol,kdsol,pvptr
      common/ cst142 /jmsol(m4,mst),kdsol(m4),pvptr(h4,2),dedpol(h4)

      integer iemod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(h9),iemod(k10)

      double precision times, btime, etime
      common/ time /times(30),btime(30),etime(30)

      double precision thermo
      common/ cst1 /thermo(k4,k10)

      double precision uf
      integer ifug, iff, idss
      common/ cst10 /uf(2),iff(2),idss(h5),ifug

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      character names*8
      common / cst8 /names(k1)

      integer did, dct
      double precision dgee
      common/ dean /dgee(k10),did(k10),dct

      logical restrt, dead
      integer ophct
      common/ lop28 /ophct,restrt,dead

      double precision pa3
      common/ cstpa3 /pa3(k19,m14)

      logical equimo
      double precision deph,dydy,dnu
      common/ cxt3r /deph(3,j3,h9),dydy(m4,j3,h9),dnu(j3,h9),equimo(h9)

      double precision dzdp, ds0dp, dgex, dcdp, gend
c                                 derivatives with respect to the p' th endmember fraction:
c                                 derivative of site fraction (endmember, species, site, model)
c                                 derivative of the mechanical negentropy (endmember,model)
c                                 derivative of the excess function (endmember,factor,term,model)
c                                 derivative of the bulk composition (component,endmember,model)
      common/ cdzdp /dzdp(m11,m10,m14,h9),
     *               ds0dp(m14,h9), dgex(m14,m2,m1,h9), 
     *               dcdp(k5,m14,h9), gend(m14)

      double precision apc, endt, endc
      common/ cstp2c /apc(h9,k5,m14), endt(h9,m14), endc(h9,m14,k5)

      integer tpct, tcct, itxp, dkp, stpct
      double precision txco
      common/ csts2d /txco(m25), tpct, tcct, itxp(m24), dkp(m24), stpct

      logical rkwak, kwak0
      integer rids, rkds, rnpt
      double precision rcp, rsum, rsmo
      common/ cxt12a /rcp(k5),rsum,rsmo,rids,rkds,rnpt,rkwak,kwak0
c                                 flag for near stable static compositions 
      logical ststbl
      common/ cststb /ststbl(k1)
c                                 local solution model logical variables
      logical depmod,laar,ordmod,recip,modres,unbd
      common/ cst160 /depmod,laar,ordmod,recip,modres,unbd
c                                 compenent degeneracy flags/counters
      logical dispro
      integer idegen, idg, jdegen, jdg
      common/ cst315 /idegen, idg(k5), jdegen, jdg(k5), dispro(k0)
c
      integer jndq, jdqf, iq
      double precision dqfg, dq
      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)
c                                 lagged chemical potentials:
      double precision xmu
      common/ cxt43 /xmu(k8)
c                                 plot_option_file options
      logical plopt
      integer piopt
      common/ cst213 /piopt(5),plopt(5)
c                                 plot program internal options
      logical spline, half, tenth, lgrid, fill, label
      integer ifont, bbox
      double precision xfac, cscale, nscale, ascale, rlabel, width, 
     *                 tcont
      common/ ops /xfac,cscale,nscale,ascale,rlabel,width,tcont,bbox(4),
     *             ifont,spline,half,tenth,lgrid,fill,label
c                                 solptr(i) points to the position of 
c                                 solution model i in the input list
      integer solptr
      common/ cst212 /solptr(h9)
c                                 global assemblage pointers
      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      character*5 zname
      common/ cst29a /zname

      integer iwt
      common/ cst209 /iwt
c                                 previous LP result amounts for
c                                 dynamic optimization starting points
      integer lcpt, lspt, ldv, lsdv, lsst
      double precision lsamt, lamt
      common/ cst120 / lsamt(k19), lamt(k19), 
     *                 lcpt, lspt, ldv(k19), lsdv(k19), lsst(k19)
c                                 nlpsol blocks
      double precision bigbnd, bigdx, bndlow, bndupp, tolact, tolfea, 
     *                 tolrnk
      common/ ngg019 /bigbnd, bigdx, bndlow, bndupp, tolact, tolfea, 
     *                tolrnk

      double precision rhomax, rhonrm, rhodmp
      common/ ngg017 /rhomax, rhonrm, rhodmp

      double precision hfwd, hctl
      common/ cxt009 /hfwd(m14), hctl(m14)

      integer itmxnp
      common/ ngg020 /itmxnp

      logical outrpc, maxs
      common/ ngg015 /outrpc, maxs

      logical fdset, cntrl, numric, fdincs
      common/ cstfds /fdset, cntrl, numric, fdincs

      integer count, rcount, lcount
      common/ cstcnt /count, rcount(5), lcount(5)
