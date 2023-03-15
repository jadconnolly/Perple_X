      subroutine mrkpur (ins, isp)
c-----------------------------------------------------------------------
c subroutine to calculate the log fugacities and volume of single
c species fluids using the hard sphere MRK EoS. 

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp     - the number of species to be calculated.
c        p,t     - input from cst5, p(bars), t(K)

c output (to common cstcoh):

c        g(i)    - fugacity coefficient of ith species
c        v(i)    - volume of the ith species

c species indices:

c         7 = O2
c        12 = O
c        13 = SiO
c        14 = SiO2
c        15 = Si  
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision bv, v4b, rt, rt3, prt, f1, f2, f3, f4, df1, df2, 
     *                 df3, dv, vi(2), fi(2), vmb

      integer i, j, itic, ir(2), ins(*), isp, k, ict
 
      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision fg
      common/ cst11 /fg(2) 

      double precision x,g,v
      common/ cstcoh /x(nsp),g(nsp),v(nsp)

      double precision vol
      common/ cst26 /vol

      double precision a, b
      common/ rkab /a(nsp),b(nsp)

      double precision pv, pvv
      integer iroots
      logical switch, rkmin, min
      common/ rkdivs /pv,pvv,iroots,switch,rkmin,min

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      rt = r*t 
      rt3 = rt*dsqrt(t)
      prt = p/rt

      call crkprm (ins,isp)

      do k = 1, isp

         i = ins(k)

         do j = 1, 2
c                                 iterate for a low and high root
            if (j.eq.1) then
               vi(j) = 2d0*b(i)
            else 
               vi(j) = 1d0/prt
            end if 

            ir(j) = 0
            itic = 0
            ict = 0
c                                iteration loop
            do 

               bv = b(i)/vi(j)
               v4b = vi(j) + 4d0*b(i)
               f1 = 1d0 + bv*(1d0 + bv*(1d0 + bv))
               df1 = (3d0*bv**3 + 2d0*bv**2 + bv)/vi(j)
               f2 = (1d0 - bv)**3 
               df2 = 3d0*(1d0 - bv)**2*bv/vi(j)
               f3 = 1d0/vi(j)/v4b
               df3 = -f3 * (1d0/vi(j) + 1d0/v4b)

               dv = (f1/vi(j)/f2 - a(i)/f3/rt3 - prt) / 
     *              (df1 - f1*(1d0/vi(j) + df2/f2)/f2/vi(j) 
     *                      + a(i)/f3**2/rt3*df3)

               vi(j) = vi(j) - dv
               itic = itic + 1

               if (dv/v(i).lt.nopt(5)) then
c                                 converged
                  ir(j) = 1
                  ict = ict + 1
                  vmb = vi(j) - b(i)
                  f4 = b(i)/vmb
c                                 compute fugacity
                  fi(j) = (6d0 + 8d0 *f4 + 4d0 * f4**2)*f4  
     *                  + dlog(vi(j)*rt/vmb**2)
     *                  + (dlog(vi(j)/v4b)/b(i)/4d0 - 1d0/v4b)*a(i)/rt3
                  exit

               else if (itic.gt.1000.or.vi(j).lt.0d0) then 

                  exit

               end if 

            end do
                  
         end do

         j = 2
 
         if (ict.eq.2) then       
c                                 found high low roots, choose the stable one:
            if ((vi(1)-vi(2))/vi(1).lt.2d0*nopt(5)
     *                    .or.fi(1).lt.fi(2))  j = 1

         else if (ict.eq.1.and.ir(1).eq.1) then

             j = 1

         end if    

         if (ict.ge.1) then 
            v(i) = vi(j)
            g(i) = dexp(fi(j))/p
         else 
c                                 failed
            write (*,*) 'failed'
            write (*,*) p,t,i,vi
            v(i) = 1d0/prt
            g(i) = 1d0

         end if 
c                                 next species:
      end do 

      end

      subroutine crkprm (ins, isp)
c-----------------------------------------------------------------------
c subroutine to return standard crk a and b terms. 

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp     - the number of species to be calculated.
c        t       - input from cst5, p(bars), t(K)

c species indices:

c         7 = O2
c        12 = O
c        13 = SiO
c        14 = SiO2
c        15 = Si  
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision brk(nsp), ark(nsp)

      integer ins(*), isp, i, k

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision a, b
      common/ rkab /a(nsp),b(nsp)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision fac, fac1
      common/junk/fac,fac1
c----------------------------------------------------------------------
c      if (first) then 
c         write (*,*) 'enter fac'
c         read (*,*) fac
c         if (fac.eq.0d0) fac=16.77
c         first = .false.
c      end if 

      do k = 1, isp

         i = ins(k)

         if (i.eq.14) then 

            if (iopt(28).eq.2) then 
c                                  shornikov 
               stop

            else if (iopt(28).eq.0.or.iopt(28).eq.1.) then         
c                                  HSC DP fit, with a0 correction
               b(14) = 0.1404308450D1
               a(14) = -0.749255053043097258D9 
     *                 + dsqrt(T) * T * 0.513703560938930059D4 
     *                 + dlog(T) * 0.139385980759038210D9 
     *                 - T * 0.268011409248898912D6 
     *                 - T**2 * 0.324444596182990921D2
c delta component :
c        *            +  fac*(t-1999.) + fac1*(t-1999.)**2

c-------------------------------------------------------------
            else if (iopt(28).eq.3.or.iopt(28).eq.4) then    
c                                 HSC decaying CP, june 2015.
               stop

            end if 

         else if (i.eq.15) then 
c                                 MRK dispersion term for Si 
            stop

         else 

            stop
 
         end if

      end do 

      
      if (iopt(28).ne.1.and.iopt(28).ne.4) then 
      
         a(13) = ark(14)/16.722d0
         b(13) = brk(14)/4.831d0

       else 
c                                    sio max a parms:
         b(13) = 16.06179418
         a(13) = 2205440030.

       end if 


      end 

      subroutine crkmix (ins, isp, iavg)
c-----------------------------------------------------------------------
c subroutine to calculate the log fugacities and volume of mixed
c species fluids using the  hard spehere MRK EoS. 

c input:

c        ins(i) -  pointers indicating the species are to be calculated.
c        isp    - the number of species to be calculated.
c        p,t    - input from cst5, p(bars), t(K)
c        iavg   - a flag indicating whether the a cross term is to 
c                 be computed as the geometric mean (iavg = 1), the
c                 arithmetic mean (iavg = 2), or the harmonic mean

c output (to common cstcoh):

c        g(i)    - fugacity coefficient of ith species
c        v(i)    - volume of the ith species

c species indices:

c         1 = H2O
c         2 = CO2
c         3 = CO
c         4 = CH4 
c         5 = H2
c         6 = H2S
c         7 = O2
c         8 = SO2
c         9 = COS
c        xx = C2H6 a = 90d6, b = 20.
c        10 = N2
c        11 = NH3
c        12 = O
c        13 = SiO
c        14 = SiO2
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical max, bad
 
      integer i, j, itic, ir(2), ins(*), isp, k, ict, iavg
 
      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision fg
      common/ cst11 /fg(2) 

      double precision x,g,v
      common/ cstcoh /x(nsp),g(nsp),v(nsp)

      double precision vol
      common/ cst26 /vol

      double precision a, b
      common/ rkab /a(nsp),b(nsp)

      double precision vrt
      integer irt
      logical sroot, nospe
      common/ rkroot /vrt,irt,sroot,nospe

      double precision pv, pvv
      integer iroots
      logical switch, rkmin, min
      common/ rkdivs /pv,pvv,iroots,switch,rkmin,min

      integer idspec
      double precision spec
      common/ tspec /spec(nsp,k5),idspec
 
      save max
c---------------------------------------------------------------------- 
      bad = .false.

      if (nospe) then 
c                                 load known composition to override
c                                 speciation routine calculations
         do k = 1, isp
            x(ins(k)) = spec(k,idspec) 
         end do 

      end if 

      call rkparm (ins,isp)

      bx = 0d0
      aij = 0d0

      do k = 1, isp

         i = ins(k)

         aik(i) = 0d0
         bx = bx + b(i)*x(i)

      end do 

      do k = 1, isp

         i = ins(k)

         do l = 1, isp

            j = ins(l)
 
            if (iavg.eq.1) then 
c                                 geometric mean mixing rule
               ax = dsqrt(a(i)*a(j))
            else if (iavg.eq.2) then 
c                                 arithmetic mean mixing rule
               ax = (a(i) + a(j))/2d0
            else 
c                                 harmonic mean mixing rule
               ax = 2d0/(1d0/a(i) + 1d0/a(j))
            end if 

            aij = aij + x(i)*x(j)*ax
            aik(i) = aik(i) + x(j)*ax

         end do 

      end do 
c                                 solve for high/low mixture molar volume
      rt = r*t 
      rt3 = rt*dsqrt(t)
      prt = p/rt
      at2 = aij/dsqrt(t)
      vvb = dlog(vol/(v

      ict = 0

      do j = 1, 2
c                                 iterate for a low and high root
         if (j.eq.1) then
            vi(j) = 2d0*b(i)
         else 
            vi(j) = 1d0/prt
         end if 

         ir(j) = 0
         itic = 0
c                                iteration loop
         do 

            bv = bx/vi(j)
            v4b = vi(j) + 4d0*b(i)
            f1 = 1d0 + bv*(1d0 + bv*(1d0 + bv))
            df1 = -(3d0*bv**3 + 2d0*bv**2 + bv)/vi(j)
            f2 = (1d0 - bv)**3 
            df2 = 3d0*(1d0 - bv)**2 * bv/vi(j)
            f3 = vi(j)*v4b
            df3 = vi(j) + v4b

            fdf = (rt/vi(j)*f1/f2 - at2/f3  - p) / 
     *            ((df1 - f1*(1d0/vi(j) + df2/f2)) *rt/f2/vi(j) + 
     *               at2 / f3 ** 2 * df3)

            vi(j) = vi(j) - dv
            itic = itic + 1

            if (dabs(dv/vi(j)).lt.nopt(5)) then
c                                 converged
               ir(j) = 1
               ict = ict + 1
               vmb = vi(j) - bx
               f4 = bx/vmb
c                                 compute fugacity
               fi(j) = (6d0 + 8d0 *f4 + 4d0 * f4**2)*f4  
     *                 + dlog(vi(j)*rt/vmb**2)
     *                 + (dlog(vi(j)/v4b)/bx/4d0 - 1d0/v4b)*aij/rt3
               exit

            else if (itic.gt.1000.or.vi(j).lt.0d0) then 

               exit

            end if 

         end do
                  
      end do

      j = 2
      rkmin = .false.
 
      if (ict.eq.2) then       
c                                 found high low roots, choose the stable one:
         if (dabs(vi(1)-vi(2))/vi(1).lt.2d0*nopt(5)
     *                    .or.fi(1).lt.fi(2))  j = 1

      else if (ict.eq.1.and.ir(1).eq.1) then

         j = 1

      end if    

      if (ict.ge.1) then 

         vol = vi(j)
         if (j.eq.1) rkmin = .true.

      else 
c                                 failed
         write (*,*) 'failed'
         write (*,*) p,t,i,vi
         vol = 1d0/prt

         bad = .true. 

      end if 
c                                 compute fugacities:
      y = bx/vol
      b4 = 4d0*bx
      v4b = vol + b4
      vvb = dlog(vol/v4b)
      ym1 = 1d0 - y
      c0 = (4d0-3d0*y)*y/ym1**2 - dlog(p*vol/rt)
      c1 = ((4d0-2d0*y)*y/ym1**3 - aij/rt32*(vvb/b4+1d0/v4b))/b
      c2 = 0.5d0/rt32/bx*vvb
 
      do i = 1, isp

         k = ins(i)
          
         if (x(k).gt.0d0) then
c                                 ln(g(k)) 
            g(k) = c0 + c1*b(k) + c2*aik(k)
c                                 f(k) is returned as the log of the fugacity
            f(k) = g(k) + dlog(p*x(k))
c                                 g(k) is returned as the fugacity coefficient
            g(k) = dexp(g(k))

         else 

            g(l) = 1d0
            f(l) = dlog(1d4*p)

         end if 

      end do 
  
      end