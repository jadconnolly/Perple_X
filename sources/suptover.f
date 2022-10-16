      open (3,file='sup.junk')
      open (10,file='supver.dat')

      call pread1
      call pread2

      end
 
      subroutine pread2

      implicit double precision (a-h,o-z)

      common /eldata/ phcmp(120), weight(120), psyco(120), nccc

      character lname*20, name*8, dummy(200)*1, comment(80)*1,
     *          record*200,rname*20

c      read in phase (mineral or gas) properties
c      if a blank line, read next line
c      if first non-blank character is a *, go to section decode area
c      decide if it is in compositional space of reaction

c      read name, formula, abbreviation (skip blank lines)
c                check line for a '(' - if not present, keep reading
c                lines till one is found (this is the way that we know
c                that we are at the start of a new set of phase data.
c                this is also used to skip phases not in the right
c                compositonal space)
  
10    read( 3,1001,end=99) dummy
      call ignore (dummy,200,' ', 1, 1,200, ifound)
      if (ifound .eq. 0) go to 10
c                                   first card contains a long name
c                                   and unreadable formula
      call decode (lname,dummy,name,200,-2,ier)
c                                   second card contains a short name
c                                   and readable formula.
      read (3,1001,end=99) dummy

      call decode (rname,dummy,name,200,-1,ier)

      if (ier.eq.1) then
         write (*,*) 'rejected ',lname
         goto 10
      end if 

      ist = 0

      do i = 3, 199
         if (dummy(i).ne.' '.and.dummy(i-1).eq.' ') ist = 1
         if (ist.ne.0) then
            comment(ist) = dummy(i)
            ist = ist + 1
            if (dummy(i).eq.' ') goto 20
         end if
      end do

20    read (3,1001,end=99) dummy

      do i = 2, 80
         if (dummy(i-1).eq.' '.and.dummy(i).eq.' ') goto 30
         comment(ist) = dummy(i)
         if (comment(ist).ne.' ') ilast = ist
         ist = ist + 1
30    end do 

      read (3,*) g,h,s,v

      read (3,1001,end=99) dummy

      call howman (dummy,inum)
      write (record,'(200a1)') (dummy(i),i=1,200)

      if (inum.lt.4) then
         itran = 0
         read (record,*) a,b,c
         read (3,*) tmax
      else
         itran = 1
         read (record,*) a,b,c,tt1
         read (3,1001,end=99) dummy
         write (record,'(200a1)') (dummy(i),i=1,200)
         call howman (dummy,inum)
         if (inum.lt.4) then
            read (record,*) a1,b1,c1
            read (3,*) tmax
         else
            itran = 2
            read (record,*) a1,b1,c1,tt2
            read (3,1001,end=99) dummy
            write (record,'(200a1)') (dummy(i),i=1,200)
            call howman (dummy,inum)
            if (inum.lt.4) then 
               read (record,*) a2,b2,c2
               read (3,*) tmax
            else 
               itran = 3
               write (*,*) 'ugga bugga'
               stop 
            end if
         end if
      end if

      jlam = 0

      if (itran.ne.0) jlam = itran + 3

      write (10,2000) rname,1,0,jlam,0,lname,
     *                (comment(i),i=1,ilast)
c                                   convert to oxide stoichiometry
      phcmp(1) = phcmp(1)/2.
      phcmp(3) = phcmp(3)/2.
      phcmp(5) = phcmp(5)/2.
      phcmp(11) = phcmp(11)/2.
      phcmp(10) = phcmp(10) - phcmp(1) - phcmp(2) - 3. * phcmp(3)
     *             - 2. * phcmp(4) - phcmp(5) - phcmp(6) 
     *             - 2. * phcmp(7) - phcmp(8) - phcmp(9)
     *             - phcmp(11) - 2. * phcmp(12)
      phcmp(10) = phcmp(10) / 2.
c                                   a coeff conv:
      atov = 4.184d0
c                                   b coeff conv:
      btov = 4.184d-3
c                                   c coeff conv:
      ctov = 4.184d5

      write (10,2001) (phcmp(j), j = 1, nccc)
      write (10,2002) g*atov,s*atov,v/10.*atov,a*atov,
     *                b*btov,c*ctov,
     *                0.,0.,0.,0.,0.,0.,0.,0.,0.,0.

      if (itran.eq.1) then
         write (10,2002) tt1,0.,0.,a1*atov,b1*btov,
     *                   c1*ctov,0.,0.,0.
      else if (itran.eq.2) then
         write (10,2002) tt1,0.,0.,a1*atov,b1*btov,
     *                   c1*ctov,0.,0.,0.
         write (10,2002) tt2,0.,0.,a2*atov,b2*btov,
     *                   c2*ctov,0.,0.,0.
      else if (itran.eq.3) then
         write (10,2002) tt1,0.,0.,a1*atov,b1*btov,
     *                   c1*ctov,0.,0.,0.
         write (10,2002) tt2,0.,0.,a2*atov,b2*btov,
     *                   c2*ctov,0.,0.,0.
         write (10,2002) tt3,0.,0.,a3*atov,b3*btov,
     *                   c3*ctov,0.,0.,0.
      end if 

      goto 10

1001  format (200a1)
2000  format (a8,4(i2),1x,a20,1x,80a1)
2002  format (5(g13.7,1x))
2001  format (12(f5.2,1x))

99    end 

      subroutine howman (dummy,inum)

      character*1 dummy(200)

      inum = 0

      do i = 1, 199
          if (dummy(i).ne.' '.and.dummy(i+1).eq.' ') inum = inum + 1
      end do

      end 

      subroutine decode (name,ldummy,nabv,isize,iswt,ier)

      implicit double precision (a-h,o-z)

      character*20 name
      character*8 nabv, ncomp
      character*1 ldummy(isize), ltemp(200)

      common /namele/ ncomp(120)
      common /eldata/ phcmp(120), weight(120), psyco(120), nccc
 
c     the elemental symbols are passed in ncomp (a4 format)
c     the number of elements is nctz

c     ldummy optionally containing one or more of the following:
c                (a)     a mineral name,     3 blanks (or more),
c                (b)     a chemical formula, 3 blanks (or more),
c                (c)     an abreviation,     3 blanks (or more),
c                (d)     miss. information
c                all in a1 format.   (d) is always ignored.
c                dimensioned isize
c
c     iswt contains a number .ge. -2 .and. .le. 4
c               if iswt = -2, ldummy contains  (a)
c               "    "     -1,   "       "     (a), (c)
c               "    "      0,   "       "     (a), (b)
c               "    "      1,   "       "     (a), (b), (c)
c               "    "      2,   "       "     (b), (c)
c               "    "      3,   "       "     (b)
c               "    "      4,   "       "     (c)
c
c     the mineral name is returned in name (in a20 format)
c     the formula is returned in phcmp, in order specified by symbols
c         in ncomp.  dimensioned ndim
c     the abreviation is returned in nabv (in a8 format)
 
c     three blanks (   ) are used as deliminators between fields
 
c     single blanks are ignored
 
         iout = 6
 
         nabv = '        '
         name = '                    '

         do 5 i = 1, 120
5             phcmp(i) = 0.0d0

         istart = 1
         ifound = 1
         ifin   = isize
 
c    check if name is present
 
         if (iswt .ge. 2) go to 100
 
c    find first non blank character, and then 3 blanks -- name is
c                  between these two strings
 
         call ignore(ldummy,isize,' ',1,istart,ifin,ifound)
         if (ifound.eq.0) then
                write(iout,1001) (ldummy(i),i=1,isize)
                write(iout,1002)
                stop
             end if
         istart = ifound
         call findc(ldummy,isize,'   ',3,istart,ifin,ifound)
         if (ifound.eq.0) then
                write(iout,1001) (ldummy(i),i=1,isize)
                write(iout,1003)
                stop
             end if
         inum = ifound - istart
         if (inum .gt. 20) then
            write(iout,1008) (ldummy(i), i=1,isize)
            inum = 20
          endif
         call movec(inum,isize,200,istart,ldummy,ltemp)

         do 10 i = 1, inum
10               name(i:i) = ltemp(i)
 
c    test if formula comes next
 
100   if (iswt .eq. -2) return
      if (iswt .eq. 4) go to 200
 
c    find start of formula, then 3 blanks.  in between is formula,
c    decode with elsort
 
         istart = ifound
         call ignore(ldummy,isize,' ',1,istart,ifin,ifound)
         if (ifound.eq.0) then
            write(iout,1001) (ldummy(i),i=1,isize)
            write(iout,1004)
            stop
         end if
         istart = ifound
         call findc(ldummy,isize,'   ',3,istart,ifin,ifound)
         if (ifound.eq.0) then
                write(iout,1001) (ldummy(i),i=1,isize)
                write(iout,1005)
                stop
             end if
         inum = ifound - istart
         call movec(inum,isize,200,istart,ldummy,ltemp)
         if ( inum .lt. 200) then
               do 150 i = (inum+1), 200
                    ltemp(i) = ' '
  150         continue
           else if ( inum .gt. 200) then
               write(iout,1009) (ldummy(i), i=1,isize)
               stop
           end if

      call elsort (ltemp, 200, ier)

      if (ier.eq.1.or.iswt.lt.0) goto 99

c    test if abreviation comes next

200   if (iswt .eq.  0) return
      if (iswt .eq.  3) return

c    find start of abreviation, then 3 blanks.  abv is between

         istart = ifound
         call ignore(ldummy,isize,' ',1,istart,ifin,ifound)
         if (ifound.eq.0) then
                write(iout,1001) (ldummy(i),i=1,isize)
                write(iout,1006)
                stop
             end if
         istart = ifound
         call findc(ldummy,isize,'   ',3,istart,ifin,ifound)
         if (ifound.eq.0) then
                write(iout,1001) (ldummy(i),i=1,isize)
                write(iout,1007)
                stop
             end if
         inum = ifound - istart
         if (inum .gt. 8) then
             inum = 8
             write(iout,1010) (ldummy(i),i=1,isize)
           endif
         call movec(inum,isize,200,istart,ldummy,ltemp)
 
         do 250 i = 1, inum
250         nabv(i:i) = ltemp(i)

1001  format(' error in thermodynamic data bank.  line in error is',/,
     .             ' ',200a1)
1002  format(' in particular, the entire line is blank thus there',
     .         ' is no phase/species name.')
1003  format(' in particular, the name is too long or 3 blanks',
     .         ' were not found as a terminator.')
1004  format(' in particular, no formula was found.')
1005  format(' in particular, the formula is too long or 3 blanks',
     .         ' were not found as a terminator.')
1006  format(' in particular, no abreviation was found.')
1007  format(' in particular, the abreviation is too long or 3',
     .         ' blanks were not found as a terminator.')
1008  format(' truncation of phase name to 20 characters ',
     .         ' line in question is: ',/,
     .         200a1)
1009  format(' formula truncated - please correct data base ',/,
     .         ' and rerun.  line in question is: ',/,
     .        200a1)
1010  format(' abreviation truncated.  line in question is:',/,
     .        200a1)
99    end

      subroutine pread1

      implicit double precision (a-h,o-z)
 
      character*8 ncomp
      character*1 dummy(200)
 
      common /eldata/ phcmp(120), weight(120), psyco(120), nccc
      common /namele/ ncomp(120)
 
        read (3,*) nccc
        i = 0
        do 30 j = 1, (1 + ((nccc - 1)/7))
            read (3,1002) dummy
            istart = 1
            kk     = nccc - (7 * (j - 1))
            if (kk .gt. 7) kk = 7
            do 20 k = 1, kk
                call ignore(dummy,200, ' ', 1, istart,200, ifound)
                if (ifound .eq. 0) go to 20
                istart = ifound
                call findc (dummy,200, ' ', 1, istart,200, ifound)
                if (ifound .eq. 0) go to 20
                i = i + 1
                do 10 l = 1, 8
                    if (l .le. (ifound - istart) ) then
                        ncomp(i)(l:l) = dummy((istart + l - 1))
                    else
                        ncomp(i)(l:l) = ' '
                    endif
   10          continue
                istart = ifound
   20      continue
   30  continue
 
1002  format(200a1)
      end
 
 

      subroutine ignore (l1,l1dim,l2,l2dim,istart,ifin,ifound)

      implicit double precision (a-h,o-z)
c
c       this routine searchs for the the first occurance (starting at
c       istart and ending at ifin) of a l2dim long character string
c       in the array l1 (dimensioned l1dim) that does not match the
c       l2dim long character string in the array l2.
c       if a string that does not match is found, then ifound is set
c       to the first position of that string (in l1).
c       if the string is a perfect match, ifound is set to 0.

      character*1 l1(l1dim), l2(l2dim)

      iend = ifin + 1 - l2dim
      do 20 i = istart, iend, l2dim
            do 10 j = 1, l2dim
               ii = i + j - 1
               if( l1(ii) .ne. l2(j)) then
               ifound = i
               goto 99
            end if
10       continue
20    continue
      ifound = 0
99    end

      subroutine findc (l1,l1dim,l2,l2dim,istart,ifin,ifound)
c
c       this routine searchs the character array l1 (starting
c       at location istart, ending at location ifin) for the string
c       with length l2dim which is stored in the character array l2.
c       if found, ifound is set the starting position of the string,
c       if not, ifound is set to zero.

       character*1 l1(l1dim), l2(l2dim)

       ifound = 0
       iend = ifin - l2dim + 1

         do 20 i = istart, iend
            do 10 j = 1, l2dim
               ii = i + j - 1
10             if (l1(ii) .ne. l2(j)) go to 20
               ifound = i
               goto 99
20    continue
99    end

      subroutine movec (num,n1,n2,istart,l1,l2)

c       this routine transfers num characters from the character
c       array l1 (dimensioned n1) starting at the istart'th character
c       to the character array l2 (dimensioned n2).
c       this routine does not check for array overflow, so be careful.

      character*1 l1(n1), l2(n2)

      do 10 i = 1, num
         j = istart + i - 1
10       l2(i) = l1(j)
 
      end

      subroutine elsort (ldummy,len,ier)

      implicit double precision (a-h,o-z)

      character*20 buffer
      character*8 name1, ncomp
      character*1 ldummy(len), lname1(8), lblank
 
      common /eldata/ phcmp(120), weight(120), psyco(120), nccc
      common /namele/ ncomp(120)

c    this routine decodes a formula which has the general form
c                         ca(3)al(2)si(3)o(12)
c    see notes on format below

c    input:  ldummy - contains the formula (see notes on format below)
c            len    - dimension of ldummy
c            elname - elemental symbols
c            nchem  - dimension of comp and elname
c            numel  - number of elemental symbols in elname
c            iswt   - if iswt .ge. 0 and an unknown element if found
c                                    (symbol not in elname), adds symbol
c                                    to elname and adds 1 to numel
c                   - if iswt .lt. 0 and an unknown element if found
c                                    (symbol not in elname), an error
c                                    message is printed and program
c                                    stops
c            iout   - unit to which error messages are written

c    output: comp   - contains the coefficients for each of the elements
c                   - zero if an element is not present
c                   - ordered by the order of the atomic symbols in
c                     elname

c    notes on format of the formula:

c       elemental symbols must be separated by the sequence: left
c                         deliminator, coefficient, right deliminator
c       the left parenthesis '(' is the left deliminator of a
c                         coefficient
c       the right parenthesis ')' is the right deliminator of a
c                         coefficient
c       deliminators must be paired
c       the coefficients can be either integers or real numbers
c       elemental symbols/coefficients can be repeated as many times as
c                 need be.  ca(1)ca(3)ca(5.0)..... is prefectly legal

c       repeat groups can also be used.  the deliminators are [ and ],
c              respectively left and right, and must be paired.  they
c              must be followed by ( coeff ) sequence.

      data lblank/' '/

      do 1 i = 1, 120
1        phcmp(i) = 0.d0

      ier = 0

         ifin   = 0
         irepl  = 0
         irepr  = 0
         irepn  = 0
         irepn2 = 0
         repeat = 0.0
c
   10  istart = ifin + 1
         ifin = len
         if (istart .eq. irepl) istart = istart + 1
         if (istart .eq. irepr) istart = irepn  + 1
         if (istart .ge. len) return
c
c    look for [  if present, repeat group exists.
c       only look if a repeat group is not currently active
c
         if (istart .ge. irepl .and. istart .ge. irepr) then
             call findc(ldummy,len,'[',1,istart,ifin,ifound)
             if (ifound .le. 0) then
                irepl  = 0
                irepr  = len + 1
                irepn  = len + 1
                repeat = 1.0
             else
                irepl = ifound
                call findc(ldummy,len,']',1,irepl,ifin,ifound)
                if (ifound .gt. 0) then
                   irepr = ifound
                   irepn = ifound + 1
                else
                   write(6,1400) (ldummy(i),i=1,len)
                   stop
                endif
c
c              now  decode repeat coefficient
c
                if (ldummy(irepn) .ne. '(' ) then
                   write(6,1410) (ldummy(i),i=1,len)
                   stop
                endif
                call findc(ldummy,len,')',1,irepn,ifin,ifound)
                if (ifound.eq.0) then
                     write(6,1100) (ldummy(i),i=1,len)
                     stop
                end if
c
                irepn2 = ifound - 1
                irepn  = irepn + 1
                call findc(ldummy,len,'.',1,irepn,irepn2,ifound)
                if (ifound.ne.0) then
                   write(buffer,'(20a1)') (ldummy(i),i=irepn, irepn2)
                   read(buffer,'(f12.3)',err=100) repeat
                else
                   inum = 6 - (irepn2 - irepn)
                   if(inum.ne.0) write(buffer,'(20a1)')
     .                     (lblank,j=1,inum),(ldummy(i),i=irepn,irepn2)
                   if(inum.eq.0) write(buffer,'(20a1)')
     .                     (ldummy(i),i=irepn,irepn2)
                   read(buffer,'(i7)',err=100) i
                   repeat = float(i)
                endif
                irepn = irepn2 + 1
             endif
         endif
         ifin = len
c
c
c    find a "(".  value between the 0'th position or the last ")" and
c      this "(" is the symbol of an element.
c
         call findc(ldummy,len,'(',1,istart,ifin,ifound)
         if (ifound.eq.0) return
c
         inum = ifound - istart
         name1 = '        '
         call movec(inum,len,8,istart,ldummy,lname1)
         do 15 i = 1, inum
               name1(i:i) = lname1(i)
   15   continue
c
c    if it is a new element, store name and obtain index.  if it
c      has alread been input, get index.
c
      do 20 i = 1, nccc
         index = i
20       if(ncomp(i).eq.name1) go to 30
c                                     name not found
      ier = 1
      goto 99

c    find the next ")".  characters between last "(" and this ")"
c      correspond to the number of moles the the element just read
c      in.  this element has location index in the composition array
c      comp.  store number of moles here.

   30  istart = ifound + 1
         call findc(ldummy,len,')',1,istart,ifin,ifound)
         if (ifound.eq.0) then
                 write(6,1100) (ldummy(i),i=1,len)
                 stop
            end if
c
         ifin = ifound - 1
         call findc(ldummy,len,'.',1,istart,ifin,ifound)
         if (ifound.ne.0) then
              write(buffer,'(20a1)') (ldummy(i),i=istart,ifin)
              read(buffer,'(f12.3)',err=100) x
         else
              inum = 6 - (ifin - istart)
              if(inum.ne.0) write(buffer,'(20a1)') (lblank,j=1,inum),
     .                                       (ldummy(i),i=istart,ifin)
              if(inum.eq.0) write(buffer,'(20a1)') (ldummy(i),i=istart
     .                                                          ,ifin)
              read(buffer,'(i7)',err=100) i
              x = float(i)
         end if
         if (istart .ge. irepl .and. istart .le. irepr)
     .                        x  = x * repeat
         phcmp(index) = phcmp(index) + x
c
         ifin = ifin + 1
         go to 10
c
  100   write(6,1300) (ldummy(i),i=1,len)
         stop
c
 1100   format(' one of the coefficients in a formula does not',
     .         ' have a closing',/,
     .        ' bracket.  line in question is:',/,
     .        ' ',200a1)
 1200   format(' one of the elemental symbols could not be recognized.',
     .         /,'  symbol in question is: ',a8,/,
     .         ' and is contained in the line: ',/,' ',200a1)
 1300   format(' one of the coefficients',
     .         ' in a formula has non-numeric symbols ',
     .       /,' embeded in it.  line in question is: ',/,
     .         ' ',200a1)
 1400   format(' a repeat group does',
     .         ' not have a closing ] ',/,' line in questions is:',/,
     .         ' ',200a1)
 1410   format(' coefficient does',
     .         ' not immediately follow a repeat group.',/,
     .         ' line in question is: ',/,' ',200a1)
99    end
