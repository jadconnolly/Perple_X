c   program to read the original hp2010 thermodynamic datafile (DS6 thermocalc 
c   files) to old perplex format, the result must then be converted to the new
c   perplex format with the program rewrite_2010.f (which must be linked to 
c   tlib.f with a modified perplex_parameters.h file, see header comments in 
c   rewrite_2010.f). hp2010tover.f has no dependencies. JADC, 10/2022 

      program trans

      open (9, file='tc-ds636.txt', status= 'old')
      open (10,file='junk.dat')

      call rdin
      end

      subroutine rdin
      
      implicit none

c   this subroutine reads part of the H&P data


      character     text(132)*1, name*8, cnum*80, gnum*80, twod*80 
      character     snum*80, vnum*80, lname*40
      logical aq
      double precision rnum, comp(19), rgib, g, reas, s, reav,sfe,
     *                 tr,b1,b5,b6,b7,b8,dsf,atoms,l4,l5,l6,
     *                 catoms(19),patoms,lam,kp,kpp,dkdt
      double precision v, a, b, c, e, newa, k298, l1, l2, l3,szr,smn
      double precision ctoj,sk,sna,sca,sc,sti,sal,ssi,smg,so,sh,scl,sni
      integer       index, igst, igen, isst, isen, ivst, iven
      integer       ihptov(19),i,inen,nst,j,inst,icomp,ilam,itype
c     hp ind      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 
      data ihptov/4,7,3,9,2,8,6,1,5,13,14,15,12,19,10,11,18,16,17/            !    2016
      data catoms/3.,2.,5.,3.,3.,2.,3.,2.,2.,2.,3.,2.,2.,3.,3.,2.,5.,2.,
     *            0./

c component indices in hp vs ver:

c                    hp'10?       ver      hp'16 (622)
c na2o               4            1        8
c mgo                7            2        5
c al2o3              3            3        3 
c sio2               9            4        1
c k2o                2            5        9
c cao                8            6        7
c tio2               6            7        2
c mno                1            8        6
c feo                5            9        4
c nio                15           10       15
c zro2               16           11       16
c cl2                13           12       13
c o2                 10           13       10
c h2o                11           14       11
c co2                12           15       12
c cuo                18           16       18
c cr2o3              19           17       19 
c s2                 17           18       17
c electron           14           19       14 

1000  format(132a1)
2000  format(a8,4(i2),' H= ',g14.7)
2002  format(6(g14.7,1x))
2001  format(20(f5.2,1x))
c                          entropies of the elements at 298, 1bar
          icomp = 18

          sk = 64.68
          sna = 51.30
          sca = 41.63
          sc  = 5.74
          sti = 30.63
          sal = 28.35
          ssi = 18.81
          smg = 32.68
          smn = 32.01
          sfe = 27.28 
c                           values from Robie
          so = 205.20 /2.
          sh = 130.70 /2.
          ctoj = 4.184
          scl = ctoj*53.288/2.
          sni = ctoj*7.14
          szr = ctoj*9.32
       
5       read(9,1000,end=900) text

c   inst and inen specify the positions of initial and terminal 
c   character of the phase name in the array called text
c   the phase name is then written to the variable name

      inst = 0
      inen = 0

      do i = 1, 132
            if (text(i) .ne. ' ') then
                  inst = i
                  goto 10
            end if
      end do

10      do i = inst, 132
            if (text(i) .eq. ' ') then
                  inen = i-1
                  goto 20
            end if
      end do
c                          blank card, assume end of
c                          file.
      goto 900 

20      write(name,1000) (text(i), i = inst, inen)
c                          ds636 adds long name, locate first number
      do i = 21, 132
            if (text(i) .ge. '0'.and.text(i).le.'9') then

               if (text(i).ne.'0') then 
                  write (*,*) 'hm'
               end if
c                          the first integer is some code, jump it
               inst = i + 3
               exit 
            end if
      end do
c                          ds636 adds long name starting in column 8, 
c                          locate the next blank
      do i = 8, 132
            if (text(i) .eq. ' ') then

c                          the first integer is some code, jump it
               inen = i - 1
               exit 
            end if
      end do

      write(lname,1000) (text(i), i = 8, inen)
c                          now find the next number
      do i = inen + 1, 132
            if (text(i) .ge. '0'.and.text(i).le.'9') then

               if (text(i).ne.'0') then 
                  write (*,*) 'hm'
               end if
c                          the first integer is some code, jump it
               inst = i + 3
               exit 
            end if
      end do
       

c      reading indices and real volues of components from the array called text
c      indices are read and then written into the variable index as integers.
c      Then the value for the component corresponding to the actual index 
c      is read and written ("internal write" statement) to the varaible cnum.
c      cnum is written to rnum which is then assigned to the corresponding 
c      element of the array compon.

      atoms = 0
      comp(:) = 0d0

      do i = 1, 19
            nst = inst + (i-1)*12
c                                            end of composition signalled by zero


            if (text(nst).eq.'0'.and.text(nst-1).eq.' ') goto 50
           write(twod,1000) text(nst-1), text(nst)
           read(twod,*) index
            write(cnum,1000) (text(j), j = nst+1,nst+9)
            read(cnum,*) rnum
            comp(ihptov(index)) = rnum
            atoms = atoms + rnum
      end do
c                                            ds6.36: the first card is lam, g, s, v
c                                            and can be read in the standard fortran way
        
50       read (9,*) lam, g, s, v
c                                            end of first line in h&p data file

c                                            as the 2nd and 3rd line of the data file have a fixed format, 
c                                            they are read in a standard fortran way
         read(9,*) a, b, c, e
         read(9,*) newa, k298, kp

         if (k298.ne.0) then 
            kpp = -kp/k298
         else 
c                                            a gas
            kpp = 0
         end if

         dkdt = 0d0
         itype = 8
         aq = .false.

       ieos = 

       ilam = lam * 10

      if(lam.gt.0.and.ilam.lt.2) then
c                                            ds6.36 transition data is 
         backspace(9)

         if (lam.eq.1) then 

c            text = ' '
c            read (9,1000) text
c            read (text,*) newa, k298, kp, l1, l2, l3
c                             distinguish type by landau entropy 
            read(9,*) newa, k298, kp, l1, l2, l3
            write (*,*) 'landau ',name

         else if (lam.eq.2) then 

            read(9,*) newa, k298, kp, l1, l2, l3,l4,l5,l6
            write (*,*) 'bragg ',name
c                             correction for -fac
            if (l6.lt.0d0) l6 = (-l6*l5 + 1d0)/(l5 + 1d0)

         else 

            write (*,*) 'unknown transition code lam = ',lam,name
            call errpau

         end if 


      else if (lam.eq.-1) then 
c                             aqueous species
            aq = .true.
            backspace (9)
            read(9,*) newa, k298, kp, kpp, lam, c
            lam = 0


      else if (lam.lt.0.or.name.eq.'lcL') then 
c                             hp98 style eos
         backspace(9)
         read(9,*) newa, k298, kp, kpp, dkdt
         write (*,*) 'old style endmember ',name
         itype = 9
      else if (lam.gt.2) then 
         write (*,*) 'ugga bugga ',name
      end if

c                                            ds6.36 adds a blank line, skip it:
         read(9,*)

       g = g * 1000.
       s = s * 1000.
       a = a * 1000.
       b = b * 1000.
       c = c * 1000.
       e = e * 1000.
       tr =298.15d0

       b1 = newa
       b8 = kp
       b6 = k298*1d3
       b7 = kpp/1d3

       if (aq.or.v.eq.0d0) then 
c                                    v .eq. 0 should catch gas species
          b5 = 0

       else if (dkdt.ne.0d0) then 
c                                    dkdt
          b5 = dkdt*1000.

       else 
c                                    einstein T
          b5 = 10636./(s/atoms+6.44)

       end if 
c                                     type flag
c                                      0 - old
c                                      1 - new
c                                      2 - ideal gas
      if (aq) then
         itype = 15 
      else if (b8.eq.0d0) then
         itype = 1
      end if 

       write (10,2000) name,itype,0,ilam,0,g
c                                    convert hp H to G
       dsf = s        - comp(1) * sna - comp(2) * smg
     *                - comp(3) * sal - comp(4) * ssi
     *                - comp(5) * sk  - comp(6) * sca
     *                - comp(7) * sti - comp(8) * smn
     *                - comp(9) * sfe - comp(10) * sni
     *                - comp(11) * szr - comp(12) * scl
     *                - comp(13) * so - comp(14) * sh
     *                - comp(15) * sc 
c                                   use apparent gibbs energy
       g = g - 298.15 * s
c                                   convert to oxide stoichiometry
       comp(1) = comp(1)/2.
       comp(3) = comp(3)/2.
       comp(5) = comp(5)/2.
       comp(14) = comp(14)/2.
       comp(13) = comp(13) - comp(1) - comp(2) - 3. * comp(3)
     *             - 2. * comp(4) - comp(5) - comp(6) 
     *             - 2. * comp(7) - comp(8) - comp(9)
     *             - comp(10) - 2. * comp(11)
     *             - comp(14) - 2. * comp(15)
     *             - comp(16) - 3./2. * comp(17) 

       comp(13) = comp(13) / 2.
       comp(12) = comp(12) / 2.
c        cr2o3
       comp(17) = comp(17)/2.
c        cuo
       comp(16) = comp(16)
c        s2 
       comp(18) = comp(18)/2.

       patoms = 0d0

       do i = 1, 19
          patoms = patoms + comp(i)*catoms(i)
       end do 

       if (patoms.ne.atoms) then 
          write (*,*) 'correct stoich of:',name
c          write (10,*) 'correct stoich of:',name
       end if 

       write (10,2001) (comp(j),j=1,icomp)

       if (aq) then 
       write (10,2002) g,s,v,c,b,-comp(19),
     *                 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
       else 
       write (10,2002) g,s,v,a,b,c,0.,e,0.,0.,b1,0.,0.,0.,b5,b6,b7,b8
       end if 
      if (lam.eq.1d0) then
         write (10,2002) l1,l2*1d3,l3,0.,0.,0.,0.,0.,0.,0.
      else if (lam.eq.2) then 
         write (10,2002) l1*1d3,l2,l3*1d3,l4,l5,l6,0.,0.,0.,0.
      end if 

       do i = 1, 19
          comp(i) = 0.
       end do 

       goto 5

900     end