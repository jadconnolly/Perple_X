c----------------------------------------------------------------------
c STX_2_ver converts Stixrude's species data files into a single perplex
c thermodynamic data file.

c the species data file names are listed in species_names.txt

c the output is written to STXVER.dat
c----------------------------------------------------------------------

      implicit none
 
      include 'perplex_parameters.h'

      integer i, ier

      character spname*8

      character card*(lchar), key*22, val*3,
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40

      double precision parm(44)

      open (n1,file='species_names.txt',status='old')
      open (n3,file='STXVER.dat',status='unknown')

      do 

         read (n1,*,end=99) spname
         open (n2,file=spname,status='old')

         call redcd1 (n2,ier,key,val,nval1,nval2,nval3,strg,strg1)

         write (n3,'(a8,a,2(a,1x),/)') spname,' EoS = 6 | ', key, strg

         do i = 1, 43
            read (n2,*) parm(i)
         end do

        close (n2)

        if (parm(9).ne.0d0) call errdbg ('non-zero k"')
        if (parm(32).ne.0d0) call errdbg ('vinet')
c                                 F0, -n, V0
         write (n3,1000) parm(5)*1d3,-int(parm(1)),parm(6)/1d1
c                                 K0 (7), K' (8), Debeye T (10), 
c                                 gamma0(26), q0 (27), eta_S0 (37)
         write (n3,1010) parm(7)*1d4,parm(8),parm(10),parm(26),
     *                  parm(27),parm(37)
c                                 mu0 (35), mu' (36)
         write (n3,1020) parm(35)*1d4, parm(36)

         if (parm(41).ne.1d0) then 
            write (n3,*) 'van laar size = ',parm(41)
         end if

         if (parm(38).ne.0) then 
c                                 landau low, Tc (38), Sd (39), Vd (40)
            write (n3,1030) parm(38), parm(39), parm(40)/1d1
         end if

         write (n3,'(a,//)') 'end'

      end do

1000  format ('G0 = ',f11.2,' S0 = ',i3,' V0 = ',f8.4)
1010  format ('c1 = ',f9.1,' c2 = ',f8.5,' c3 = ',f10.5,' c4 = ',f7.5
     *       ,' c5 = ',f7.5,' c6 = ',f7.5)
1020  format ('m0 = ',f9.1,' m1 = ',f9.5)
1030  format ('transition = 1 type = 7  t1 = ',f11.5,' t2 = ',f10.5,
     *        ' t3 = ',f8.4)

99    close (n1)
      close (n3)

      end



