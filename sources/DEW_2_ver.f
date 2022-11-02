c----------------------------------------------------------------------
c DEW_2_ver converts aqueous species data from the DEW spreadsheet 
c to perple_x format/units
c----------------------------------------------------------------------

      implicit none
 
      include 'perplex_parameters.h'

      integer i, ibeg, iend, siz, ier, iscan, iscnlt, jbeg(3), jend(3),
     *        nel, iel, ist, j, lun, d2v(13)

      logical elchk, good, reject, bad, hsc

      double precision nums(13),sel(50),stoich,elst(50),ost(50),ox(50), 
     *                 otot

      character rec*(lchar), name8*8, elem(50)*2, elname*2, number*8,
     *          oxname(50)*5, text(14)*1, two*2

      external iscan, iscnlt, elchk

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
c     iam = 11
c                                 the DEW data, space delimited
      open (n9, file='DEW_input.txt', status= 'old')
c                                 a case sensitive list of acceptable elements
c                                 species that contain elements not in this list
c                                 will be rejected.
      open (11,file='DEW_elements_input.txt',iostat=ier,status='old')
      nel = 0

      read (11,*,iostat=ier) reject 
c                                 use HSC convention for g0, otherwise use gf.
      read (11,*) hsc

      do
c                                 oxygen must be last in list
         read (11,'(a,1x,a5,1x,f2.0,1x,f2.0)',iostat=ier) elname, 
     *                               oxname(50), nums(1), nums(2)
         if (ier.ne.0.or.elname.eq.'en') exit
         nel = nel + 1
         elem(nel) = elname
         oxname(nel) = oxname(50)
         elst(nel) = nums(1)
         ost(nel) = nums(2)
      end do 

      ier = 0 

      if (nel.eq.0) then 
         write (*,*) 'gork'
         stop
      end if 

      lun = 10
      open (lun,file='DEW_2_ver_output.dat')
c                                 DEW data consists of a name, formula, thermodata, comment
c                                 pointers from dew spread sheet columns 
c                                 Gf, Hf, S0, v0, cp0,    a1, a2, a3, a4, c1, c2, omega, z comment
c                                  1   2   3   4    5      6   7   8   9  10  11    12  13   
c                                 to perplex EoS 16 variables (cp & v0 are unused):
c                                'G0 ','S0 ','V0 ','Cp0','w ','q ','a1 ','a2 ','a3 ','a4 ','c1 ','c2 ','HOH'
      do i = 1, 13
         d2v(i) = -1 
      end do 
c                                 g 
      d2v(1) = 1
c                                 s 
      d2v(2) = 3
c                                 omega
      d2v(5) = 12
c                                 z
      d2v(6) = 13
c                                 a1
      d2v(7) = 6
c                                 a2
      d2v(8) = 7
c                                 a3
      d2v(9) = 8
c                                 a4
      d2v(10) = 9
c                                 c1
      d2v(11) = 10
c                                 c2
      d2v(12) = 11

      do 

         bad = .false. 

         call redcd3 (n9,ier)
         if (ier.ne.0) exit
c                                find end of first name
         ibeg = 1 
         iend = iscan (ibeg,com,' ') - 1

         jbeg(1) = 1
         jend(1) = iend
c                                find the second name
         ibeg = iscnlt (iend+1,com,' ')
         iend = iscan (ibeg,com,' ') - 1 
c                                decode the name
         ist = ibeg
         sel = 0d0

         do
c                                assume first char is beginning of an element name
            if (chars(ist+1).le.'Z'.and.chars(ist+1).ge.'A'.or.
     *         chars(ist+1).eq.'('.or.
     *         chars(ist+1).eq.')'.or.
     *         chars(ist+1).eq.' '.or.
     *         chars(ist+1).eq.'+'.or.
     *         chars(ist+1).eq.'-'.or.
     *         chars(ist+1).le.'9'.and.chars(ist+1).ge.'1') then 
c                                one character element
               write(elname,'(2a)') chars(ist)
               ist = ist + 1
            else if (chars(ist+1).le.'z'.and.chars(ist+1).ge.'a') then
c                                 two character name
               write(elname,'(2a)') chars(ist:ist+1)
               ist = ist + 2
            else 
               write (*,*) 'fugga: ',chars(ibeg:iend)
               stop
            end if 
c                                 check the element
            good = elchk(elname,elem,nel,iel)
            if (.not.good) bad = .true.
c                                 skip an unused repeat group
            if (chars(ist).eq.')') ist = ist + 1 
c                                 find its stoichiometry
            if (chars(ist).le.'9'.and.chars(ist).ge.'1') then

               read (chars(ist),*) stoich

               if (chars(ist+1).le.'9'.and.chars(ist+1).ge.'0') then

                  if (chars(ist+1).eq.'0') then 
                     write (*,*) chars(1:jend(1)),' wudda been wrong689'
                  end if 
c                                 two digit numeral, no provision for larger numbers...
                  write (two,'(a,a)') chars(ist:ist+1)
                  read (two,*) stoich
                  ist = ist + 1

               end if

               ist = ist + 1
               if (chars(ist).eq.')') ist = ist + 1 

            else if (chars(ist).le.'Z'.and.chars(ist).ge.'A'.or.
     *               chars(ist).eq.'('.or.ist.eq.iend) then

               stoich = 1d0

            end if 
c                                 check if already at the end
            if (chars(ist).eq.'(') then

               if (iend-ist.gt.7) then
                  i = ist + 8
               else 
                  i = iend - 1
               end if 
          
               write(number,'(8a)') chars(ist+1:i)
               read (number,'(i8)',iostat=ier) j

               if (ier.eq.0.or.ist+1.eq.i) then 

                 ist = iend 

               else

                  write (*,*) 'found a repeat group ', 
     *                        'these seem to be useless ',
     *                        'will ignore, CHECK THIS FORMULA!'
                  write (*,*) 'mugga: ',chars(ibeg:iend)
                  ist = ist + 1

               end if 

            else if (chars(ist).eq.'+'.or.chars(ist).eq.'-') then

               if (chars(ist-1).le.'Z'.and.chars(ist-1).ge.'A'.or.
     *             chars(ist-1).le.'z'.and.chars(ist-1).ge.'a') then
c                                 the previous character was an element
c                                 and there is no stoichiometric coeff
                  stoich = 1d0
               end if

               ist = iend 

            end if 

            if (good) sel(iel) = sel(iel) + stoich
c                                 end of a repeat group
            if (chars(ist).eq.')')  ist = ist + 1

            if (ist.eq.iend.and.
     *          chars(ist).le.'Z'.and.chars(ist).ge.'A') then
               cycle

            end if

            if (ist.ge.iend) exit

         end do 

         if (reject.and.bad) cycle


         jbeg(2) = ibeg
         jend(2) = iend
c                                find a comment, if there is one
         ibeg = iscnlt (iend+1,com,'@') 
         jbeg(3) = ibeg
         jend(3) = com
c                                write the numbers to record
         write (rec,'(400a)') chars(jend(2)+1:com)
c                                read the numbers
         read (rec,*,iostat=ier) nums

         if (ier.ne.0) then
            write (*,*) 'probably missing enthalpy for:',
     *                  chars(jbeg(1):jend(1))
            read (rec,*,iostat=ier) nums(1),(nums(i),i=3,13)
            nums(2) = 0d0
            ier = 0 
         end if 

         if (chars(jbeg(1)).eq.'"') then 
            jbeg(1) = jbeg(1) + 1
            jend(1) = jend(1) - 1
         end if

         if (jbeg(3).eq.jend(3)) jend(3) = jbeg(3) - 1

         if (jend(1)-jbeg(1).gt.7) then
            write (*,*) 'truncated name from ',
     *                  chars(jbeg(1):jend(1)),' to ',
     *                  chars(jbeg(1):jbeg(1)+7)
            jend(1) = jbeg(1) + 7
         end if 

         write (name8,'(8a)') chars(jbeg(1):jend(1))

         write (lun,'(/,a,a,f9.0,400a)') name8,' EoS = 16 | Gf = ',
     *                                  nums(1)*4.184,
c                                 formula
     *                                chars(jbeg(2)-1:jend(2)+1),
c                                 comment
     *                                chars(jbeg(3)-1:jend(3))
c                                 get oxide stoichiometry
         ox = 0d0
         otot = 0d0 

         if (.not.bad) then 

            do i = 1, nel
               if (elst(i).eq.0d0) cycle
               ox(i) = sel(i)/elst(i)
               if (i.lt.nel) otot = otot + ox(i)*ost(i)
            end do 

            ox(nel) = ox(nel) - otot/2d0

         else 

            ox(nel) = -1d0

         end if 

         ibeg = 1
 
         do i = 1, nel 

            if (ox(i).ne.0) then 
c                                 load text name
               iend = ibeg + 4

               read (oxname(i),'(5a)') chars(ibeg:iend)
c                                 left parenthesis
               chars(iend + 1) = '('
c                                 get number
               call numtxt (ox(i),text,siz)
c                                 load number into chars
               ibeg = iend + 2
               iend = ibeg + siz - 1

               do j = ibeg, iend
                  chars(j) = text(j-ibeg+1)
               end do 

               chars(j) = ')'
 
               ibeg = iend + 2

            end if 

         end do 

         iend = iend + 1

         call ftext (1,iend)
c                                 write the formula
         write (lun,'(400a)') chars(1:iend)
c                                 write the DEW spreadsheet options as a comment:
c         write (lun,'(3(a,i1))') 
c     *         '   | DEW correlation options => omega: ',int(nums(8)),
c     *                               ' HKF parameters: ',int(nums(9)),
c     *                               ' a1 Correlation: ',int(nums(10))
c                                 =====================================
c                                 thermo data conversions:
         do i = 1, 13

            if (i.eq.4) then
c                                 unused parameters v0
               nums(i) = nums(i)/1d1
            else if (i.eq.13) then
c                                 charge
            else 
               nums(i) = 4.184d0*nums(i)
            end if

         end do 
c                                 scale omega, a1, a2, a4, c2
         nums(12) = nums(12)*1d5
         nums(6)  = nums(6)/1d1
         nums(7)  = nums(7)*1d2
         nums(9)  = nums(9)*1d4
         nums(11) = nums(11)*1d4

         ibeg = 1
 
         do i = 1, 6

            j = d2v(i)

            if (d2v(i).lt.0) cycle

            if (j.eq.1) then 
               if (HSC) nums(1) = nums(2) - 298.15*nums(3)
            end if 

            call outthr (nums(j),e16st(i),3,ibeg)

         end do

c                                 write G,S,V,cp,w,q
         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)
c                                 write a1-a4,c1-c2, HOH
         ibeg = 1
 
         do i = 7, 12
            j = d2v(i)
            call outthr (nums(j),e16st(i),3,ibeg)
         end do    

         if (name8.eq.'H+') then 
            call outthr (1d0,e16st(13),3,ibeg)
         else if (name8.eq.'OH-') then 
            call outthr (2d0,e16st(13),3,ibeg)
         end if 

         if (ibeg.gt.1) write (lun,'(400a)') chars(1:ibeg)

         write (lun,'(a,/)') 'end'

      end do

      end 

      logical function elchk (elname,elem,nel,iel)
c----------------------------------------------------------------------
c----------------------------------------------------------------------   
      integer nel, iel, i

      character elem(50)*2, elname*2

      elchk = .false.
 
      do i = 1, nel
         if (elname.eq.elem(i)) then 
            iel = i
            elchk = .true. 
            exit
         end if 
      end do

      end 

      subroutine redcd3 (nloc,ier)
c----------------------------------------------------------------------
c readcd - read 240 column card image from unit 9, strip out unwanted
c characters. ier = 1 no card found.
c----------------------------------------------------------------------    
      implicit none

      include 'perplex_parameters.h'

      integer i, ier, iscan, ict, iscnlt, ibeg, nloc

      external iscan, iscnlt

      character card*(lchar)

      integer length,com
      character chars*1
      common/ cst51 /length,com,chars(lchar)
c----------------------------------------------------------------------

      ier = 0 

      ibeg = 0
  
      com = 0 

      card = ' '

      do while (ibeg.ge.com) 

         read (nloc,'(a)',end=90) card

         if (card.ne.' ') then 

            read (card,'(400a)') chars
c                                 find end of data marker '|'
            com = iscan (1,240,'|') - 1
c                                 '|' in first column
            if (com.eq.0) cycle
c                                 find a non blank character
            ibeg = iscnlt (1,com,' ')

         end if 

      end do

      ict = 1

      do i = 2, com
c                                 strip out '+' and '*' chars
c                                 eliminate blanks after '/' and '-'
c                                 and double blanks
         if (chars(ict).eq.' '.and.chars(i  ).eq.' ') cycle
         ict = ict + 1
         chars(ict) = chars(i)

      end do 

      com = ict

      goto 99

90    ier = 3

99    end