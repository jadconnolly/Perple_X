
      program unsplt
c----------------------------------------------------------------------
c       lun restrictions
c       n4 - local plot
c       n5 - local bulk
c       n6 - global plot
c       n7 - global bulk
c       n8 - test
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, fake, match, nomtch, reord, used(k5), err,
     *        ok

      integer h, i, j, k, l, m, n, o, p, gasct, ggrd(l7,l7), gi,
     *        gsoct, loc2gs(h9), loc2ga(k3), gdasls(k5,k3), gavar(3,k3),
     *        gap(k2), gbulk, gico(k2), gjco(k2), gicoor(k1), inct, gj,
     *        ioct, ids, jxco, kxco, gjxco, gkxco, gas, kbad,kb(1200,2),
     *        loc2gb(k2), ist, jst, gstg(h9), gspg(h9,mst), gcoor(h9)

      character gprjct*100, dir*100, gname(h9)*10

      double precision gbg(k5,k2), gmus(k5,k2), gxcoor(k18),
     *                 xt(k5,mst*msp), bt(k5)

      logical oned
      common/ cst82 /oned

      integer igrd
      common/ cst311 /igrd(l7,l7)

      integer jlow,jlev,loopx,loopy,jinc1
      common/ cst312 /jlow,jlev,loopx,loopy,jinc1

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer ivar,ind
      common/ cst83 /ivar,ind

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer iam
      common/ cst4 /iam

      logical fileio, flsh, anneal, short
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal,short

      character*100 prject,tfname
      common/ cst228 /prject,tfname
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer jx, jy, lev, xn, yn
      common/ cst58 /jx, jy, lev, xn, yn

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      double precision bg
      common/ cxt19 /bg(k5,k2)

      integer icog,jcog
      common/ cxt17 /icog(k2),jcog(k2)

      double precision amu
      common/ cst48 /amu(k8,k2)
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 13
c                                 version info
      call vrsion (6)
c                                 initialize some flags
      first = .true.
      outprt = .false.
      fake   = .false.
c                                 read input from unit n1 (terminal/disk).

c                                 input1 also initializes:
c                                 equilibrium counters; units n2 n4 and n6;
c                                 and the limits for numerical results.
      call input1 (first,err)
c                                 set ivar flag, this indicates the number
c                                 of possible independent plotting variables, jvar
c                                 indicates the number of thermodynamic variables
      ivar = 2
c                                 don't allow users to do anything
c                                 other than gridded min
      if (icopt.ne.5) call error (4,1d0,icopt,'UNSPLT')
c                                 read thermodynamic data on unit n2:
      call input2 (fake)
c                                 read autorefine lists
      call setau1
c                                 read data for solution phases on n9:
      call input9 (fake)

      call setau2
c                                 initialize, set global lists and project name
      gsoct = isoct

      do i = 1, gsoct
         gname(i) = fname(i)
         gstg(i) = istg(i,1)
         gcoor(i) = ncoor(i)
         do j = 1, istg(i,1)
            gspg(i,j) = gspg(i,j)
         end do 
      end do 

      gprjct = prject
      gasct = 0 
      gbulk = 0 
      gjxco = 0

      do i = 1, l7
         do j = 1, l7
            ggrd(i,j) = 0
         end do 
      end do

      write (*,'(/)')
c                                 at this point we have a complete list
c                                 of possible phases and solutions.

c                                 inialize segment counter:
      k = 0
c                                 and bad segment counter:
      kbad = 0
c                                 switch iam so input routines
c                                 know the program is doing segments
      iam = 14
c                                 loop through the segments:
c----------------------------------------------------------------------
c                                 MARK: if segements numbering is reordered
c                                 to increase columnwise switch the order
c                                 of the next two statements:
c----------------------------------------------------------------------
      do j = 1, jx
         do i = 1, jy
c                                 initialize err just in case
            err = .false.
c                                 segment index
            k = k + 1
c                                 segment folder
            write (dir,'(i4,a)') k
            call mertxt (dir,dir,'/',0)
            call mertxt (dir,'segment',dir,0)
            call mertxt (dir,'/',dir,0)
            call mertxt (dir,gprjct,dir,0)
c                                 local project
            call mertxt (prject,dir,gprjct,0)
c                                 for final version comment from here:
c----------------------------------------------------------------------
c                                 MARK: temporary fix overwrite above names
c                                 assuming the perplex directory (to 
c                                 which unsplt will write the output
c                                 contains a folder mark_IV in with the
c                                 programs and data files:
c            write (dir,'(i4,a)') k
c            call mertxt (dir,dir,'/',0)
c            call mertxt (dir,'yanliu8/working/segment',dir,0)
c            call mertxt (dir,'working/segment',dir,0)
c            call mertxt (prject,dir,gprjct,0)
c                                                            to here:
c----------------------------------------------------------------------

            write (*,'(a,a)') 'processing: ',dir
c                                 read input, this is probably not all
c                                 necessary, in particular the thermo data 
            call input1 (first,err)
c                                 if err the files haven't been written
c                                 by the splitter
            if (err) goto 10
c                                 read thermodynamic data on unit n2:
            call input2 (fake)
c                                 read autorefine lists
            call setau1
c                                 read data for solution phases on n9:
            call input9 (fake)
c                                 does what? 
            call setau2
c                                 make pointers from the local solution model 
c                                 indices to the global indices:
            do l = 1, isoct
               do m = 1, gsoct
                  if (fname(l).ne.gname(m)) cycle 
                  loc2gs(l) = m 
                  exit 
               end do
            end do  
c                                 read local plt/blk files => loopx, loopy, jinc1
            call interm (.false.,err)

10          if (err.or.iasct.eq.0.or.ibulk.eq.0) then
               kbad = kbad + 1
               kb(kbad,1) = i
               kb(kbad,2) = j
               write (*,1000) dir 
               cycle 
            end if

c                                 cycle through assemblages
            do l = 1, iasct
c                                 replace local solution pointer with global
               do m = 1, iavar(1,l)
                  idasls(m,l) = loc2gs(idasls(m,l))
               end do 

               match = .false.

               if (k.ne.1) then
c                                 check if local assemblage matches
c                                 a global assemblage
                  do m = 1, gasct

                     if (gavar(1,m).ne.iavar(1,l).or.
     *                   gavar(2,m).ne.iavar(2,l)) cycle 

                     nomtch = .false.

                     do n = 1, iavar(3,l)
c                                  two criterion, if ok every phase
c                                  was found, if nomatch every phase
c                                  was found with the same multiplicity
                        ok = .false. 

                        do o = 1, iavar(3,l) 

                           if (idasls(n,l).eq.gdasls(o,m)) then

                              ok = .true. 
c                                 check that the phase occurs the same 
c                                 number of times in each assemblage:
                              inct = 0 
                              ioct = 0

                              do p = 1, iavar(1,l)
                                 if (idasls(n,l).eq.idasls(p,l)) 
     *                                                  inct = inct + 1
                                 if (idasls(n,l).eq.gdasls(p,m)) 
     *                                                  ioct = ioct + 1
                              end do 

                              if (ioct.ne.inct) then
                                 nomtch = .true.
                                 exit
                              end if 
                           end if
                             
                        end do
                     
                        if (.not.ok) nomtch = .true.
                        if (nomtch) exit

                     end do 

                     if (nomtch) cycle
                     match = .true.
                     loc2ga(l) = m
                     ias = l
                     exit

                  end do 
 
                  if (match) then 
c                                 check if the phases are in the same order
                     reord = .false.

                     do m = 1, iavar(3,l)
                        if (idasls(m,l).eq.gdasls(m,loc2ga(l))) cycle
c                                 have to cycle through all ibulk results and
c                                 reorder the phases as in sorter.
                        reord = .true.
                        exit 
                     end do 
                      
                     if (reord) then 

                        do m = 1, ibulk

                           if (iap(m).ne.ias) cycle 
c                                 make arrays of the local amount bt, composition xt,
c                                 and number of compositions coords for each phase
                           do n = 1, iavar(3,ias)
                              bt(n) = bg(n,ias)
                              if (n.gt.iavar(1,ias)) cycle
c                                  load solution compositions into xt
                              jxco = icox(m)

                              do o = 1, iavar(1,ias)

                                 used(o) = .false.
                                 ids = idasls(o,ias)

                                 kxco = jxco + gcoor(ids)
                                 jxco = jxco + 1

                                 do h = jxco, kxco
                                    xt(o,h-jxco+1) = xco(h)    
                                 end do 

                                 jxco = kxco

                              end do 
                           end do  

                           gas = loc2ga(ias)
                           jxco = icox(m)

                           do n = 1, gavar(3,gas)
                              do o = 1, gavar(3,gas)
                                 if (used(o).or.
     *                               idasls(o,ias).ne.gdasls(n,gas)) 
     *                                                             cycle
                                 gbg(n,gas) = bt(o)
c                                 THIS NEEDS TO BE CHECKED, it was iavar(o,ias). 1/9/2019
                                 if (o.gt.iavar(1,ias)) cycle
                                 
                                 do p = 1, gcoor(idasls(o,ias))
                                    xco(jxco+p) = xt(o,p)
                                 end do 

                                 jxco = jxco + gcoor(idasls(o,ias))
                                 used = .true.

                              end do 
                           end do 

                        end do 

                     end if 

                  end if 

               end if
                
               if (k.eq.1.or..not.match) then 
c                                 add the assemblage to the global list        
                  gasct = gasct + 1
                  if (gasct.gt.k3) call error (184,0d0,k3,'UNSPLT')
                  loc2ga(l) = gasct
                  gavar(1,gasct) = iavar(1,l)
                  gavar(2,gasct) = iavar(2,l)
                  gavar(3,gasct) = iavar(3,l)
                  do m = 1, iavar(3,l)
                     gdasls(m,gasct) = idasls(m,l)
                  end do 
               end if 

            end do  
c                                 the assemblage data is sorted, now
c                                 relocate the ibulk analyses:
            do l = 1, ibulk
c                                 compute the global coordinate
               gi = icog(l) + (i-1)*loopx - (i - 1)
               gj = jcog(l) + (j-1)*loopy - (j - 1)
c                                 global assemblage
               ias = loc2ga(iap(l))

               if (j.gt.1.and.jcog(l).eq.1 .or.
     *             i.gt.1.and.icog(l).eq.1) then
c                                 check if the local analyses replicates
c                                 an analysis already in the global array
                  match = .false.

                  do m = 1, gasct

                     if (gi.ne.gico(m).or.gj.ne.gjco(m)) cycle
c                                 coordinates match, make sure the assemblages
c                                 are the same 
                     if (gap(m).ne.ias) exit 

                     loc2gb(l) = m
                     match = .true.
                     exit

                  end do 
 
                  if (match) cycle

               end if 

               gbulk = gbulk + 1
c                                 assemblage pointer
               gap(gbulk) = ias
c                                 global coordinates
               gico(gbulk) = gi
               gjco(gbulk) = gj
c                                 bulk pointer
               loc2gb(l) = gbulk
c                                 molar amounts
               do m = 1, gavar(3,ias)
                  gbg(m,gbulk) = bg(m,l)
               end do 
c                                 solution compositions
               gicoor(gbulk) = gjxco
               jxco = icox(l)

               do m = 1, gavar(1,ias)

                  ids = gdasls(m,ias)

                  gjxco = gjxco + 1
                  gkxco = gjxco + gcoor(ids) - 1
                  jxco = jxco + 1
                  kxco = jxco + gcoor(ids) - 1

                  if (gkxco.gt.k18) call error (61,0d0,k18,'UNSPLT')

                  o = 0 

                  do n = gjxco, gkxco
                     gxcoor(n) = xco(jxco+o)
                     o = o + 1
                  end do 

                  gjxco = gkxco
                  jxco = kxco

               end do 

c                                 read mu's if available
               if (jpot.ne.1) then
c                                 if error on read most probably its
c                                 because of NaN's for the chemical 
c                                 potentials
                  do m = 1, jbulk
                     gmus(m,gbulk) = amu(m,l)
                  end do 
  
               end if 

            end do
c                                 map igrd to ggrd
            if (i.eq.1) then
               ist = 1
            else 
               ist = 2
            end if 

            if (j.eq.1) then 
               jst = 1
            else 
               jst = 2
            end if 

            do l = ist, loopx
               do m = jst, loopy
             
                  n = l + (i-1)*loopx - (i-1)
                  o = m + (j-1)*loopy - (j-1)

                  ggrd(n,o) = loc2gb(igrd(l,m))

               end do
            end do     

         end do 
      end do

      if (kbad.eq.k) then 

         call error (72,bt(1),1,dir)

      else if (kbad.ne.0) then 
c                                 fill in bad segments with the null
c                                 assemblage k3
         do k = 1, kbad

            i = kb(k,1)
            J = kb(k,2)

            if (i.eq.1) then
               ist = 1
            else 
               ist = 2
            end if 

            if (j.eq.1) then 
               jst = 1
            else 
               jst = 2
            end if 

            do l = ist, loopx
               do m = jst, loopy
             
                  n = l + (i-1)*loopx - (i-1)
                  o = m + (j-1)*loopy - (j-1)

                  ggrd(n,o) = k2

               end do
            end do   

         end do 

      end if 
c                                 write the global array results to
c                                 local directory
      call mertxt (tfname,gprjct,'.plt',0)
      open (n4, file = tfname)
      call mertxt (tfname,gprjct,'.blk',0)
      open (n5, file = tfname)

      if (jx*jy.ne.1) then 
         loopx = jx*loopx - (jx-1)
         loopy = jy*loopy - (jy-1)

         if (loopx.gt.l7.or.loopy.gt.l7) call error (2,bt(1),j,
     *'grid coordinates, decrease final_grid_resolution or increase L7')

      end if 

      iasct = gasct

      do i = 1, iasct
                 
         iavar(1,i) = gavar(1,i)
         iavar(2,i) = gavar(2,i)
         iavar(3,i) = gavar(3,i)
         do j = 1, iavar(3,i)
            idasls(j,i) = gdasls(j,i)
         end do 

      end do 

      do i = 1, loopx
         do j = 1, loopy
            igrd(i,j) = ggrd(i,j)
         end do 
      end do 

      call outgrd (loopx,loopy,jinc1,n4,0)
c                                 and the blk file:
      do i = 1, gbulk
       
         write (n5,'(3(i8,1x))') gico(i),gjco(i),gap(i)
         write (n5,1010) (gbg(j,i),j=1,gavar(3,gap(i)))
         jxco = gicoor(i)
         do j = 1, gavar(1,gap(i))
            kxco = jxco + gcoor(gdasls(j,gap(i)))         
            jxco = jxco + 1
            write (n5,1010) (gxcoor(k), k = jxco, kxco)
            jxco = kxco
         end do 

         if (jpot.ne.1) write (n5,1010) (amu(j,i), j = 1, jbulk)

      end do 

      close (n4)
      close (n5)

1000  format (/,'**warning ver061** error processing: ',a,/,
     *       'most probably VERTEX has not finished the segment or ',
     *       'terminated prematurely',/,'because of an error. In the ',
     *       'latter case, correct the error and recalculate the ',/,
     *       'segment.',//,'UNSPLT will depict the failed segment ',
     *       'with the null assemblage (red).',/)
1010  format (20(g16.8,1x))

      end 
