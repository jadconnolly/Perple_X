        !This is intended to be a module to define subroutines to call to Meemum
        !Might inculde vertex in this in the future, to be determined
        !This is designed with intent to implement it into Julia
        !Note that all array values (Including strings) passed into a subroutine 
        !must be the same size as the fortran array

        module perplexwrap
        

        implicit none
        include 'perplex_parameters.h'

        contains

            subroutine initMeemum(fileName,nameLength, componentNames,sysCompo,componentMass)
            !This subroutine initializes Meemum
            !It is included so that if a user desires iterative modification of system P,T, or composition
                !PerpleX will not need to read the datafile over and over again
            !This will take a filename for the .dat build file and return the system composition stored in the file
            !Note that componentNames will just be returned as a single very large string
            
                integer, intent(in) :: nameLength
                character (len=nameLength), intent(in) :: fileName
                character (len=5), dimension(k5), intent(out) :: componentNames
                double precision, dimension(3,k5), intent(out) :: sysCompo
                double precision, dimension(k0), intent(out) :: componentMass
               
            !----------------------------------------------------------------------
            !                                 these common blocks are necessary to 
            !                                 communicate between the MAIN program
            !                                 and perplex:

            !                                 jbulk -> the number of chemical components
            !                                 cblk -> the molar bulk composition
                integer jbulk
                double precision cblk
                common/ cst300 /cblk(k5),jbulk
            !                                 iam -> a variable indicating the Perple_X program
                integer iam
                common/ cst4 /iam
            !                                 perplex option values:
                integer iopt
                logical lopt
                character valu*3
                double precision nopt
                common/ opts /nopt(i10),iopt(i10),lopt(i10),valu(i10)
            !                                 ntot - number of phases stable
            !                                 np - number of solution phases
            !                                 ncpd - number of compounds
            !                                 kkp(ntot) - pointer to cpd or solution model
            !                                 cp3(1:jbulk,1:ntot) - molar phase compositions
                integer kkp,np,ncpd,ntot
                double precision cp3,amt
                common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
            !                                 phase and system properties loaded by 
            !                                 subroutine getloc, the call to getloc
            !                                 is unnecessary if only the molar phase
            !                                 proportions and compositions are of interest.
            !                                 phase (prop(i,1:ntot)) and system (psys(i))
            !                                 are for index i (this list may not be exhaustive):
            !                                  1  - molar volume
            !                                  2  - molar enthalpy
            !                                  3  - gruneisen thermal parm
            !                                  4  - K_S
            !                                  5  - Mu_S
            !                                  6  - v_phi
            !                                  7  - v_p
            !                                  8  - v_s
            !                                  9  - v_p/v_s
            !                                  10 - rho
            !                                  11 - G
            !                                  12 - cp
            !                                  13 - alpha
            !                                  14 - beta
            !                                  15 - S
            !                                  16 - molar amount
            !                                  17 - molar weight
            !                                  18 - KS_T
            !                                  19 - MuS_T
            !                                  20 - KS_P
            !                                  21 - MuS_P
            !                                  22 - vphi_T
            !                                  23 - vp_T
            !                                  24 - vs_T
            !                                  25 - vphi_P
            !                                  26 - vs_P
            !                                  27 - vp_P
            !                                  28 - heat capacity ratio (cp/cv)
                double precision props,psys,psys1,pgeo,pgeo1
                common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)
            !                                 atwt -> g-formula wts of the chemical components
                
                double precision atwt
                common/ cst45 /atwt(k0)
            !                                 v -> the standard perplex potential variables
                double precision v,tr,pr,r,ps
                common/ cst5  /v(l2),tr,pr,r,ps
            !                                 ipot -> the number of potential variables in use
            !                                 jv -> the indices of the potential variables
                integer ipot,jv,iv
                common / cst24 /ipot,jv(l2),iv(l2)
            !                                 vname -> the name of the potential variables
                character*8 vname,xname
                common/ csta2  /xname(k5),vname(l2)
            !                                 cname -> the names of the components
                character*5 cname
                common/ csta4 /cname(k5)
           
            !                                 phase compositions
                double precision pcomp
                common/ cst324 /pcomp(k0,k5)
            !                                 phase names
                character pname*14
                common/ cxt21a /pname(k5)
            !                                 iwt => 1 mass bulk comp, 0 molar bulk comp
                integer iwt
                common/ cst209 /iwt

                character*100 prject,tfname
                common/ cst228 /prject,tfname

                integer icont
                double precision dblk,cx
                common/ cst314 /dblk(3,k5),cx(2),icont
                
                

                
                
            !----------------------------------------------------------------------
                iam = 2
                
                !do i = 1, len(fileName)
               !     dummy = fileName(0,i)
               ! end do
                getInput = .false.
                sWarn =.false.
                
                
                prject = fileName
                !print*, prject
            !                  initialization, read files etc. 
                call iniprp
                
               
                componentNames = cname
                
                !print*, componentNames
                sysCompo = dblk
                componentMass = atwt
                meemumInit = .true.

            end subroutine initMeemum

            subroutine minimizePoint(componentNames, sysCompo, P, T,suppressWarn,cPotentials,phaseNames,phaseProps,phaseComps,sysProps)
                !Calls the main meemum minimization subroutine using sysCompo as the system composition
                !This assumes that meemum was initialized properly and ideally sysCompo should be fed in the same order as 
                !was returned from initMeemum but will check just in case
                !Assumes input composition is in moles
                !The output is a large character array that contains the normal terminal output of meemum to be parsed by the calling program
                !This isnt elegant but its less of a headache than figuring out where all the important variables are stored

                character (len=5), dimension(k5), intent(in) :: componentNames
                double precision, dimension(k5), intent(in) :: sysCompo
                double precision, intent(in) :: P, T !Pressure and temperature in bars and K
                
                logical(kind=1),intent(in) :: suppressWarn !MUST BE KIND=1 TO WORK WITH JULIA
                double precision, dimension(k8), intent(out) :: cPotentials
                character(len=14), dimension(k5),intent(out)::phaseNames
                double precision, dimension(i8,k5), intent(out)::phaseProps
                double precision, dimension(k0,k5), intent(out)::phaseComps
                double precision, dimension(i8), intent(out)::sysProps
                logical bad

                integer i, j

                character cprop*18
 !----------------------------------------------------------------------
            !                                 these common blocks are necessary to 
            !                                 communicate between the MAIN program
            !                                 and perplex:

            !                                 jbulk -> the number of chemical components
            !                                 cblk -> the molar bulk composition
                integer jbulk
                double precision cblk
                character amount*6, yes*1

                double precision num
                common/ cst300 /cblk(k5),jbulk
            !                                 iam -> a variable indicating the Perple_X program
                integer iam
                common/ cst4 /iam
            !                                 perplex option values:
                integer iopt
                logical lopt
                character valu*3
                double precision nopt
                common/ opts /nopt(i10),iopt(i10),lopt(i10),valu(i10)
            !                                 ntot - number of phases stable
            !                                 np - number of solution phases
            !                                 ncpd - number of compounds
            !                                 kkp(ntot) - pointer to cpd or solution model
            !                                 cp3(1:jbulk,1:ntot) - molar phase compositions
                integer kkp,np,ncpd,ntot
                double precision cp3,amt
                common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
            !                                 phase and system properties loaded by 
            !                                 subroutine getloc, the call to getloc
            !                                 is unnecessary if only the molar phase
            !                                 proportions and compositions are of interest.
            !                                 phase (prop(i,1:ntot)) and system (psys(i))
            !                                 are for index i (this list may not be exhaustive):
            !                                  1  - molar volume
            !                                  2  - molar enthalpy
            !                                  3  - gruneisen thermal parm
            !                                  4  - K_S
            !                                  5  - Mu_S
            !                                  6  - v_phi
            !                                  7  - v_p
            !                                  8  - v_s
            !                                  9  - v_p/v_s
            !                                  10 - rho
            !                                  11 - G
            !                                  12 - cp
            !                                  13 - alpha
            !                                  14 - beta
            !                                  15 - S
            !                                  16 - molar amount
            !                                  17 - molar weight
            !                                  18 - KS_T
            !                                  19 - MuS_T
            !                                  20 - KS_P
            !                                  21 - MuS_P
            !                                  22 - vphi_T
            !                                  23 - vp_T
            !                                  24 - vs_T
            !                                  25 - vphi_P
            !                                  26 - vs_P
            !                                  27 - vp_P
            !                                  28 - heat capacity ratio (cp/cv)
                double precision props,psys,psys1,pgeo,pgeo1
                common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)
            !                                 atwt -> g-formula wts of the chemical components
                double precision atwt
                common/ cst45 /atwt(k0)
            !                                 v -> the standard perplex potential variables
                double precision v,tr,pr,r,ps
                common/ cst5  /v(l2),tr,pr,r,ps
            !                                 ipot -> the number of potential variables in use
            !                                 jv -> the indices of the potential variables
                integer ipot,jv,iv
                common / cst24 /ipot,jv(l2),iv(l2)
            !                                 vname -> the name of the potential variables
                character*8 vname,xname
                common/ csta2  /xname(k5),vname(l2)
            !                                 cname -> the names of the components
                character*5 cname
                common/ csta4 /cname(k5)
            !                                chemical potentials
                logical mus
                double precision mu
                common/ cst330 /mu(k8),mus
            !                                 phase compositions
                double precision pcomp
                common/ cst324 /pcomp(k0,k5)
            !                                 phase names
                character pname*14
                common/ cxt21a /pname(k5)
            !                                 iwt => 1 mass bulk comp, 0 molar bulk comp
                integer iwt
                common/ cst209 /iwt

                character*100 prject,tfname
                common/ cst228 /prject,tfname

                integer icont
                double precision dblk,cx
                common/ cst314 /dblk(3,k5),cx(2),icont
                
                integer io3,io4,io9
                common / cst41 /io3,io4,io9

                double precision goodc, badc
                common/ cst20 /goodc(3),badc(3)
            !----------------------------------------------------------------------
                
                ! print*,suppressWarn
               
                sWarn = suppressWarn
                if (meemumInit) then
                    iam = 2
                    v(1) = P
                    v(2) = T
                    
                   
                    !Ensure that the names line up properly
                    do i = 1, size(componentNames)

                        do j = 1, size(componentNames)
                    
                            if (componentNames(i) .eq. cname(j)) then
                                cblk(j) = sysCompo(i)
                            end if
                        end do
                    end do
                    call meemum (bad)

                    if (bad) write (*,*) 'optimization failed'

                    if (.not.bad) then
                        !                                print summary to LUN 6
                        !Output pname, props, psys, pcomp, mu, psys1, fbulk, fbulk1
                        !Aqueous: spnams, ysp
                        ! call calpr0(6)
                        !Set the output variables to their corresponding values
                        cPotentials = mu
                        phaseNames = pname
                        phaseProps = props
                        phaseComps = pcomp
                        sysProps = psys

                        !reset values in case of another system call
                        mu = mu*0.0
                        props = props*0.0
                        psys = psys*0.0
                        pcomp = pcomp*0.0
                        do i = 1,k5
                            pname(i) = "              "
                        end do

            
                    end if 
            
                    if (goodc(1)+badc(1).gt.0d0) then
            
                        num = badc(1)/(badc(1)+goodc(1))*1d2
                        if (num.gt.1d-1) call warn (53,num,i,'MEEMUM')
                   
                    end if 

                else
                    print*,"Please initialize Meemum"
                end if


            end subroutine minimizePoint

            subroutine pseudosection(fileName, nameLength,grid,gridToAssem,assemToPhase,purePhases,solPhases,numPhases,xMin,xMax,yMin,yMax,xInc,yInc,xVar,yVar)
                !This is a sort of clone of the pssect program, it will read the output given from vertex
                !And provide the information needed to plot a pseudosection
                !This includes variable ranges, variable names, and the pseudosection grid (In all its complexities)

                implicit none

                
                
                integer, intent(in) :: nameLength
                character (len=nameLength), intent(in) :: fileName
                
                integer, dimension(l7,l7), intent(out) :: grid !Corresponds to igrd
                integer, dimension(k2), intent(out) :: gridToAssem !Corresponds to iap
                integer, dimension(k5,k3), intent(out) :: assemToPhase !Corresponds to idasls
                character (len=8), dimension(k1), intent(out) :: purePhases !Corresponds to names
                character (len=6), dimension(h9), intent(out) :: solPhases !Corresponds to aname
                integer,dimension(k3),intent(out) :: numPhases
                double precision, intent(out) :: xMin,xMax,yMin,yMax,xInc,yInc
                character (len = 8), intent(out) :: xVar, yVar

                integer :: i
                
                integer jop0
          
                logical first, err
          
                character yes*1
           
                integer iop0 
                common / basic /iop0
          
                integer iam
                common/ cst4 /iam

                character*100 prject,tfname
                common/ cst228 /prject,tfname

                integer igrd
                common/ cst311 /igrd(l7,l7)

                integer iap,ibulk
                common/ cst74  /iap(k2),ibulk

                character fname*10, aname*6, lname*22
                common/ csta7 /fname(h9),aname(h9),lname(h9)

                integer jvar
                double precision var,dvr,vmn,vmx
                common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

                character vnm*8
                common/ cxt18a /vnm(l3)

                integer jlow,jlev,loopx,loopy,jinc
                common/ cst312 /jlow,jlev,loopx,loopy,jinc 
                !idasls iavar names in header
          !----------------------------------------------------------------------- 
          !                                 iam is a flag indicating the Perple_X program
                iam = 7
          !                                 version info
                ! call vrsion (6)
          
                iop0 = 0 
          
                outprt = .false.
                first  = .false.
               
                prject = fileName
          !                                 read input from unit n1 (terminal/disk).
          !                                 subroutine input1 also initializes:
          !                                 equilibrium counters; units n2 and n4;
          !                                 and the limits for numerical results.
                call input1 (first,err)
                
          !                                 don't allow users to do anything
          !                                 other than gridded min
                if (icopt.lt.5) then 
                   call error (4,0d0,icopt,'PSVDRAW')
                else if (icopt.eq.12) then
                   call error (72,0d0,icopt,'0-d infiltration results can only  plotted in tab file format')
                end if 
          !                                 read thermodynamic data on unit n2:
                call input2 (first)
                
          !                                 read autorefine lists
                call setau1
                
          !                                 read data for solution phases on n9:
                call input9 (first)
               
          
                call setau2
                
          !                                 read plot option file, set
          !                                 default transformation
                call rdopt 
                
          !                                 read the plot/blk files
                call interm (.false.,err)
                
          !                                 organize variables 
                call getvar
                
          !                                 initialize the grid parameters
                call setvar
                
          !                                 open output file 
                ! call psopen
          !                                 ask for options
                ! write (*,1000) 
                ! read (*,'(a)') yes
          
                ! if (yes.eq.'y'.or.yes.eq.'Y') iop0 = 1
          
          !                                 get user options and read
          !                                 rest of plot file, draw data
        !         call psdplt (jop0)
        !   !                                 george's data plotting routine
        !   !                                 option "plot_extra_data" = T
        !         if (plopt(4)) call psdat
          
                ! call psclos
           
                close (n5)
                close (n4)
                grid = igrd
                gridToAssem = iap
                assemToPhase = idasls
                purePhases = names
                solPhases = aname
                numPhases = iavar(3,1:k3)
                xVar = vnm(1)
                yVar = vnm(2)
                xMin = vmn(1)
                yMin = vmn(2)
                xMax = vmx(1)
                yMax = vmx(2)
                xInc = dvr(1)*jInc
                yInc = dvr(2)*jInc
               
                !Reset values in case of another function call
                !Only matters for arrays of indeterminate number of useful values
                igrd = igrd*0
                iap = iap*0
                idasls = idasls*0
                iavar = iavar*0
                do i = 1, k1
                    names(i) = "        "
                end do

                do i =1,h9
                    aname(i) = "      "
                end do
        

            end subroutine pseudosection

        end module perplexwrap