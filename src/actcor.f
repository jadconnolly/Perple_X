      program actcor 
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *  actcor.may.1989     *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------
c a fortran program for making fixed activity corrections to the
c thermodynamic data file for vertex.  actcor creates a new data file 
c with the corrected data on unit n8
c-----------------------------------------------------------------------
c files (see vertex program documentation for additional information): 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i
                                                         
      character*8 blank8, name, y*1, test

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ic
      common/ cst42 /ic(k0)

      integer iam
      common/ cst4 /iam

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      logical eof
           
      data blank8/' '/ 
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 9
c                                 version info
      call vrsion (6)

      write (*,1300) 
c                                 open files
      call sopen 
c                                 read and echo file header
      call topn2 (4)
c                                 disable HSC conversion
      hscon = .false. 
c                                 mock pointers
      do i = 1, icmpn
         ic(i) = i
      end do
c                                 the no blurb
      write (*,1010) 
c                                 allow user to enter names:
      write (*,1030) 
      read (*,1000) y

      if (y.ne.'y'.and.y.ne.'Y') then 
c                             get the name:
100      write (*,1020) 
         read (*,1000) test

         if (test.eq.blank8) goto 999

         rewind n2
         call eohead (n2)

         do 
 
            call getphi (name,.false.,eof)
 
            if (eof) then 
               write (*,1050) test
               goto 100
            end if 

            if (name.eq.test) then
               call gotcha (name)
               goto 100
            end if 

         end do 

      else 
c                             read and modify individual entries  
         do 
 
            call getphi (name,.false.,eof)
 
            if (eof) exit

            write (*,1040) name
            read (*,1000) y
            if (y.eq.'y'.or.y.eq.'Y') call gotcha (name)   

         end do 
                     
      end if

1000  format (a) 
1010  format ('This program will create a new thermodynamic data',/,
     *        'file with (optionally) activity corrected entries.',/,
     *        'You must specify all phases that are to be included',/,
     *        'in the new data file (actcor.dat).',//)
1020  format ('Enter a phase to be included [<9 characters, blank to ',
     *        'finish]:')
1030  format ('Prompt for phases (y/n)?')
1040  format ('Include (y/n): ',a)
1050  format ('No such phase as: ',a)
1300  format (/,'NO is the default answer to all prompts',/)           

999   end 

      subroutine gotcha (name)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
                                                         
      character*8 blank8,name,y*1

      integer i

      double precision xmole,xmix,act

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer eos
      common/ cst303 /eos(k10)

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      integer ictr, itrans
      double precision ctrans
      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
                                                            
      data blank8/'        '/
c----------------------------------------------------------------------
      write (*,1060) name
      read (*,1000) y  
                       
      if (y.eq.'y'.or.y.eq.'Y') then 
c 
         write (*,1070) name
         read (*,1000) blank8
         write (*,1080) name
         write (*,2000) (cmpnt(i),i=1,icmpn)
         write (*,1040) (comp(i),i=1,icmpn) 
         write (*,1090)

         read (*,1000) y    

         if (y.eq.'y'.or.y.eq.'Y') then 
            write (*,1100) name,blank8

            read (*,*) xmole
               write (*,1110) name
               read (*,*) xmix

               act = xmole**xmix   

         else   
               write (*,1120) name
               read (*,*) act
         end if 

         write (*,1130) name,blank8,act 
         thermo(1,k10) = thermo(1,k10) + t * 8.314413 * dlog(act)
         thermo(2,k10) = thermo(2,k10) - 8.314413 * dlog(act)
         name = blank8

      end if 

      names(k10) = name

      eos(k10) = ieos

      lct(k10) = ilam

      ltyp(k10) = jlam

      idis(k10) = idiso

      call outdat (n8,k10,0)

1000  format (a)
1040  format (13(f5.2,1x))       
1060  format ('make an activity correction for ',a,' (y/n)?')
1070  format ('enter a unique name for the activity corrected version',
     *        ' of ',a,'(<9 characters):')
1080  format ('the stoichiometry of ',a,' is:')
1090  format (/,'ideal activity model (y/n)?')
1100  format ('enter mole fraction (x) of ',a,' in ',a,':')
1110  format ('activity of ',a,' will be computed as x**n',/,
     *        'enter number of mixing sites (n):')
1120  format ('enter activity of ',a,':')
1130  format (/,'activity of ',a,' in ',a,' is: ',g12.6)
2000  format (/,1x,13(a,1x),/,1x,13(a,1x))

      end

