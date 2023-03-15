      implicit none

      double precision phi, v, p, t

      do 

         write (*,*) 'enter p(bar), t(k): '
         read (*,*) p,t
         call crkpur (p,t,phi,v)
         write (*,'(/,a,1x,g12.6)') ' volume (J/bar) = ',v
         write (*,'(a,1x,g12.6,/)') ' fugacity coefficient = ',phi

      end do

      end 

      subroutine crkpur (p,t,phi,v)
c-----------------------------------------------------------------------

c subroutine to calculate the log fugacities and volume of single
c species fluids using the carnahan-starling-RK hybrid EoS

c input:

c        p,t     - input p(bars), t(K)

c output 

c        phi    - fugacity coefficient 
c        v      - volume 

c-----------------------------------------------------------------------
      implicit none

      double precision bv, v4b, rt, rt3, prt, f1, f2, f3, df1, df2, 
     *            df3, vi(2), fi(2), f, df, at2, fac1, v0, r, v, b,
     *            vdfmin, vdfmax, dfmax, dfmin, phi, nopt, p, t, a

      integer j, itic, ir(2), ict, jtic 

      logical fog, bad, case1

      integer iroots
c----------------------------------------------------------------------
      bad = .false.
      fog = .true.
c                                 gas constant
      r = 8.314 
c                                 argon a & b parameters (file hsrk_argon_pulmer.mws)
      a = 167809.2801d0
      b = .6201365306d0
c                                 relative convergence tolerance
      nopt = 1d-5

      rt = r*t 
      rt3 = rt*dsqrt(t)
      prt = p/rt

         at2 = a/dsqrt(t)
         v0 = rt/p
         fac1 = 1.5d0
         dfmin = 1d99
         vdfmin = 1d99
         dfmax = 0d0
         vdfmax = 0d0
         case1 = .true.
c                                 two cases: 1) 3 real roots, if this occurs
c                                 convergence to the false root can be 
c                                 avoided by expanding the distance of the
c                                 guesses from the false root (case 1 below).
c                                 2) real roots, in this case the scheme to 
c                                 solve case 1 may not yield any results, in 
c                                 this case the correct root should lie between
c                                 v_dfmax(1) and v_dfmin(2) of the 3 root search
c                         
10       ict = 0

         do j = 1, 2
c                                 iterate for a low and high root
            if (j.eq.1.and.case1) then
               vi(j) = fac1 * b
            else if (case1) then
               vi(j) = v0
            else 
c                                 hail mary
               vi(j) = (vdfmax + vdfmin)/2d0
            end if 

            ir(j) = 0
            itic = 0
            jtic = 0
c                                iteration loop
            do 

               bv = b/vi(j)
               v4b = vi(j) + 4d0*b

               f1 = 1d0 + bv*(1d0 + bv*(1d0 - bv))
               df1 = (3d0*bv**3 - 2d0*bv**2 - bv)/vi(j)
               f2 = (1d0 - bv)**3 
               df2 = 3d0*(1d0 - bv)**2 * bv/vi(j)
               f3 = vi(j)*v4b
               df3 = vi(j) + v4b

               f = rt/vi(j)*f1/f2 - at2/f3  - p

               df = (df1 - f1*(1d0/vi(j) + df2/f2))*rt/f2/vi(j) + 
     *               at2 / f3 ** 2 * df3

               if (case1) then 
                  if (j.eq.1.and.df.gt.dfmax) then 
                     dfmax = df
                     vdfmax = vi(j)
                  else if (j.eq.2) then 
                     if (vi(j).gt.vdfmax.and.df.lt.dfmin) then 
                        dfmin = df
                        vdfmin = vi(j)
                     end if 
                  end if 
               end if 

               if (vi(j)-f/df.lt.b) then 
                  vi(j) = 1.1d0*b
               else if (j.eq.2.and.dfmax.gt.0d0.and.
     *                  vi(j)-f/df.lt.vdfmax) then 
                  vi(j) = (vdfmax + vdfmin)/2d0
               else 
                  vi(j) = vi(j) - f/df
               end if 

               itic = itic + 1

               if (dabs(f/df/vi(j)).lt.nopt/100d0) then
c                                 converged
                  if (df.gt.0d0) then 
c                                 case 1 =>
c                                 got the middle root
                     if (j.eq.2) then 
c                                 try a new guess, assuming the root 
c                                 is < videal
                        v0 = vi(2) + 0.9*(v0-vi(2))
                        vi(2) = v0
                        jtic = jtic + 1

                     else if (j.eq.1) then 

                        fac1 = 1d0 + (fac1-1d0)*0.9
                        vi(1) = fac1 * b 
                        jtic = jtic + 1  

                     end if 
 
                     if (jtic.lt.10) then

                        cycle 

                     else  
c                                 ten tries and you're toast.
                        write (*,*) 'crkpur got wrong root on ',j,p,t
                        write (*,*) 'b,a,vi',b,a,vi
                        exit

                     end if 

                  end if 

                  ir(j) = 1
                  ict = ict + 1
                  v4b = vi(j) + 4d0*b
c                                 compute (log) fugacity
                  fi(j) = b*(3d0*b*(b - 3d0*vi(j))
     *                  + 8d0*vi(j)**2)/(vi(j) - b)**3  
     *                  + dlog(rt/vi(j))
     *                  + (dlog(vi(j)/v4b)/b/4d0 - 1d0/v4b)*a/rt3

                  exit

               else if (itic.gt.100) then 

                  exit

               end if 

            end do
                  
         end do

         j = 2
 
         if (ict.eq.2) then       
c                                 found high low roots, choose the stable one:
            if (dabs(vi(1)-vi(2))/vi(1).lt.2d0*nopt) then
               iroots = 1
            else 
               iroots = 3
               if (fi(1).lt.fi(2))  j = 1
            end if 

         else if (ict.eq.1.and.ir(1).eq.1) then

             j = 1

         end if    

         if (ict.ge.1) then

            v = vi(j)
            phi = dexp(fi(j))/p

         else if (case1) then 
c                                 failed
            if (fog) then
               write (*,*) 'crkpur got no root on ',itic,jtic
               write (*,*) 'b = ',b,', a = ',a,', p = ',p, 
     *                     ', t = ',t
               case1 = .false.
               goto 10
            end if 

         else 
c                                 failed
            write (*,*) 'failed crkpur case 2'
            write (*,*) p,t,vi,iroots
            v = 1d0/prt
            phi = 1d0

         end if 

      end
