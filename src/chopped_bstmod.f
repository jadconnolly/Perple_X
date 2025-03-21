         if (bstout.lt.2) then 
c                                 bstout = 0 print all converged results
            improv = i.eq.ibest.or.i.eq.jbest.or.bstout.eq.0
c                                 print loop
            lu = 6
            consol = .false.

            do
c                                 output some stats
               if (mcgrid.and.lu.eq.6) then 
c                                 george's minimal console output
                  write (*,'(a,i7,1h/,i7,1x,a,2(1x,1pg12.6),a,$)')
     *                'Try ',i,random(3),'best:',bstobj,bstbay,char(13)

               else if (.not.invxpt.and.improv) then
c                                thermobarometry
                  write (lu,'(/80(''-''))')
                  write (lu,2000) i
                  write (lu,2070) (vname(j),x(j), j = 1, ipot)
c                                write likelihood results
                  write (lu,2010) objf
               
                  if (i.eq.ibest) then
                     write (lu,2020)
                  else
                     write (lu,2030) bstobj, ibest
                  end if
c                                 write Bayeseian results
                  write (lu,2040) bay

                  if (i.eq.jbest) then
                     write (lu,2050)
                  else
                     write (lu,2060) bstbay, jbest
                  end if
c                                  unnecessary obj call for print output
                  fprint = .true.
                  call mcobj2 (x,objf,bad)
                  fprint = .false.
   
               else if (invxpt.and.
     *                 (.not.mcgrid.or.(mcgrid.and.improv))) then
c                                   george only outputs improved results         
                  if (iquad.gt.0.and.jcount.gt.0.and..not.mcgrid) then
                     write (lu,1125) i, igood, icount, jcount, objf, 
     *                              bstobj, ibest
                  else
                     write (lu,1120) i, igood, icount, objf, bstobj, 
     *               ibest
                  end if

                  write (lu,1130) ssp, bay, bstbay, jbest
                  write (lu,1080) sx(1:n)
                  write (lu,1085) x0(1:n)
                  write (lu,1030) x(1:n)

                  if (lu.eq.n6.and.invxpt) then 

                  write (lu,'(/,a,i7,a,/)') 'Scores for try = ',i,
     *                                      ' follow:'

                  do id = 1, mxpt

                     if (xptpt(id,5).eq.1d0) then
                        write (lu,1010)
     *                   id, xptnam(id),' score = ',scores(id)
                     else
                        write (lu,1010)
     *                   id, xptnam(id),' score = ',scores(id),
     *                   xptpt(id,5)
                     end if

                  end do

               end if

               write (lu,'(/80(''-''))')

            end if

            if (lu.eq.n6.or..not.n6out) exit

            lu = n6
            consol = .true.

         end do

         if (.not.mcgrid) then
c                               new random starting point
            do j = 1, n
               call random_number (sx(j))
            end do

         else
c                               new grid search point
            icount = random(2)
            random(2) = icount + 1
            
            do j = 1, n
               sx(j) = dble(mod(icount,pnum(j)))/max(1,pnum(j)-1)
               icount = icount / pnum(j)
            end do

         end if