c                                       fluid flow in 2-d reactive comp med.
c                                       with fluid layer beneath reaction front.

C                                       changes from steve2 
c                                       1) move crack to middle
c                                       2) allow reaction to create layer 
c                                          before "cracking"
c                                       3) elastic component to compaction for
c                                          elastic shock. 
c                                       4) stiff upper crust
c                                       5) added noise and differential yielding.
c
c                                            jadc
      implicit double precision (a-h,k,m,o-z)  

      integer hcrk, icrust
c
      parameter (ix=601,iy=321,imx=601)

      character*6 fname*3, p1, p2, cnum*3, title*50
	 
      double precision stim(50), nh2o, dp(ix,iy), nenth, noise,
     *                 p01(iy), x0(ix), rxn(ix,iy), trxn(ix,iy)
c
      integer ny,nx,nym1,nxm1
      double precision pn,dx,dy
      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1

      common/ data2 /po(ix,iy),c0,c1,c2,c3,c4,betat,betap,rho,mu,dtdts
     *      / data3 /ty(iy),y0(iy),p0(iy),dpz(iy),dpf(iy)
     *      / data4 /phio(ix,iy),phin(ix,iy),tvol(ix,iy),dvol(ix,iy)
     *      / data5 /qh(ix,iy),qv(ix,iy)
      common/ blok1 /nh2o,vsol,fac,enth,vtot,entr
     *      / blok2 /vand,rng,vm,dv
     *      / blok3 /erxn,t0,rnrxn,v1,v2,area,icrust
     *      / comp1 /ac,q,rn,rni,phim
     *      / data6 /k0,por(ix,iy),s(ix,iy),k(ix,iy),nphi,hcrk(ix,iy)
     *      / rates /edot(ix,iy),r(ix,iy)

      open (20,file='steve5.dat')
c                                      read input:
      read (20,1030) fname
c
      read (20,1020) title
c                                      number of x and y nodes
      read (20,*) nx, ny
c                                      max depth, max width, crack width
c                                      top and bottom of layer (m)
      read (20,*) ymax, ymin, xmax, xcrk, ycrk
c                                      thickness of initially reacted layer (m)
      read (20,*) ylay
c                                      ycrust - thickness of non-compacting
c                                      upper layer and it's permeability.
      read (20,*) ycrust, kcrust 
c                                      permeability, crack permeability
      read (20,*) kcrk
c                                      rx and fluid density, porosity, 
c                                      compressibility, viscosity
      read (20,*) rhor, rho, beta, mu
c                                      pore compressibility
      read (20,*) betap
c                                      permeability-porosity parms
c                                      background perm.
c                                      initial phi (fraction)
c                                      porosity exponent
      read (20,*) kb, phi0, nphi
c                                      minimum phi
      read (20,*) phim
c                                      read flow law:
c                                      rn stress exp
c                                      A pre-exp factor (MPa**-(rn))
c                                      Q activation energy (J)
c                                      phim, lowest phi allowed
      read (20,*) rn,a,q
c                                      geotherm K/m, front velocity m/my
      read (20,*) dtdy, ssv
c                                      wt fraction water realeased
c                                      stoich of rlm
c                                      enthalpy, volume, entr, solid
c                                      volume changes per mole H2O
c                                      assume mole rxn yields 1 mole H2O
      read (20,*) wtp
      read (20,*) vand
c                                      use mks units:
      read (20,*) enth,vtot,entr,vsol
c                                      parameters for eq 3.4
c                                      mol/m2, J/mol, K, none.
      read (20,*) c0,erxn,t0,rnrxn
c                                      initial rlm vol (m3)
      read (20,*) v0
c                                      final grain d-dimension (m)
      read (20,*) df
c                                      molar volume of rlm (m3)
      read (20,*) vm
c                                      time interval for plots, and "instantaneous" rates
      read (20,*) dtav,dtinst

      read (20,*) t_end

      close (20)
c                                      rate plot file:
      write (p2,2001) fname,'.t'
      raw = ''
      open (200,file=p2)
c                                      constants (mks)
      drho = rhor - rho
      g = 9.8d0 
c                                      initialization:
      nym1 = ny - 1
      nxm1 = nx - 1
      iplot = 0 
      ixy = 1
      tictoc = 0 
      fcum = 0d0 
      dy = (ymax-ymin)/(nym1)
      dx = xmax/nxm1
      area = dx*dy
      xmin = 0d0
      odtav = dtav
      tplot = -1d6
      tinst = 0d0
      r = 0d0
c                                      time to make crack is time (y)
c                                      to react layer of thickness ylay
      tcrky = ylay/ssv*1e6
      dtav = tcrky/10d0
c                                      get permeability constant.
      k0 = kb/phi0**nphi
c                                      conversion from years to seconds:
      ys = 365.25 * 24. * 60. * 60.
c                                      get densification paramters:

c                                      convert to Pa, etc:
c                                      e = a*dp**n*exp(-q/r/tk)*1d6**(-n)
      a = a /1d6**(rn)
      q = -q/8.314d0
      eta = 1./(3d0*a)
      zeta = 0.75*eta/phi0
      cl = dsqrt(kb/mu*zeta)
      tau = (eta/0.02)/(cl*drho*g)/ys

      write (*,*) ' compaction time, length:',tau,cl
      write (*,*) ' it will take ', tcrky,tcrky/tau,' y until faulting'
      write (*,*) ' the first 20 plots at ',dtav,' y intervals'
c                                      now: e = a*dp**rn*exp(q/tk+fh2o*fm)
c                                      the compaction constants:
c                                      stage 2 (0.02>phi):
c                                      tubes.
c                               phid = ac*(1-phi)*phi*(dp/(1-phi**rni))**rn
      ac  = a*2.**(rn+1d0)/rn**rn
      rni = 1d0/rn
c                                      compaction rate is now:
c                                      ac*exp(q/tk)*phi*(dp/(1d0-phi**rni))**(rn)
c
c                                      now reaction constants:

c                                      max moles of h2o produced
c                                      2.8d3 is rock density kg/m3
      nh2o = wtp * rhor / 18.015d-3

      erxn = -erxn/8.314
      t0 = erxn/t0
c                                      number of grains
      rng = (v0 + vand*vm*nh2o)/2d0/df**3
c                                      use average grain size:
      a1 = rng*10d0*((2d0*df**3/4d0)**(1d0/3d0))**2
c                                      factor to convert r dot
c                                      mole h2o/ cm2 RLM -s to
c                                      m3 RLM/m2 and s.
      fac = a1*vm*vand*c0 
c                                      max change in volume of RLM m3/m3 rx
      dv = nh2o*vand*vm
c                                      convert m3 RLM to mole rxn:
      v1 = -vsol/vand/vm
c                                      convert m3 volume to m3 h2o
      v2 = -18.015d-3/vsol/rho
c                                      negative enthalpy
      nenth = -enth
c                                      general junk:
      c3 = (dy/dx)**2
      betat = beta + betap
      c1 = 2d0*dx**2*mu*betat
      c4 = 2d0*dy**2*mu*betat
      c5 = rho/4d0/mu
c                                      time increment based on stabilty
c                                      criteria kdt/(cvdxdy) < 0.25:  
      da = dy
      c6 = mu*betat*da*da/4d0
      if (dx.lt.dy) da = dx
      r1 = phi0/kb
      odts = c6*r1*1d2
      dts = odts
      dty = dts / ys 
c
      write (6,*) ' time increment in years =',dts/ys 
c                                      background values:
c                                      first horizontally homogeneous
c                                      properties:
      do j = 1, ny
         y0(j)  = ymin + dy*(j-1)
         tk  = 273.15 + dtdy*y0(j)
c                                      p0 is pieziometric head:
         p0(j)  = drho*g*y0(j) + 1e5

         if (y0(j).lt.ycrust) then 
            icrust = j
            p0(j) = 0d0
         end if 
c                                      add a small number to prevent 
c                                      divide by zero. 
         p01(j) = drho*g*y0(j) + 1e5
c                                      to get deformation head, 
c                                      dp = pz - dpz ~(pf-ps)
         dpz(j) = drho*g*y0(j)
c                                      to get pf, pf = rho*g*z + pz
         dpf(j) = 1d5 + rho*g*y0(j)
c                                      check if reacted:
         call grxn (tk,p0(j)+dpf(j),dg)

         if (dg.gt.0d0) then
c                                      material can react:
            tvol(1,j) = 0d0
            tyb = tk
            iyb = j

         else 
c                                      material is not reacted
            tvol(1,j) = dv

         end if

      end do
c                                      crack bottom:
      
c                                      shift the geotherm, so the 
c                                      lowest reactive node is at the
c                                      reaction temperature:
      call equit (teq,p0(iyb)+dpf(iyb))
c                                      background flux:
      qb = rho*kb/mu*drho*g
      qvb = qb/rho
      qr = qb 
      vb = qvb/phi0
c                                      if (ssv.ne.0) compute heating
c                                      rate:
      dtdts = 0d0
      ssvs = ssv/ys/1d6
      if (ssv.ne.0d0) then 
c
	   dtdt = -ssv*(rhor*vtot/entr - dtdy)
         dtdty = dtdt/1d6
         dtdts = dtdts/1d6
c                                      reaction flux:
         qrxn = wtp*rhor*ssvs
         qvrxn = qrxn/rho
         qr = qrxn
      end if       

      ty = 273.15 + dtdy*y0 + teq - tyb

      k = kb
      s = 0d0
      rxn = 0d0
      qx = 0d0
      qy = 0d0
      trxn = 0d0
c                                      1 no noise, 0 100% noise
      noise = 1 
c                                      initialization:
      do i = 1, nx
         x0(i) = (i-1)*dx
         do j = 1, ny
            po(i,j) = p0(j)
            pn(i,j) = p0(j)
            tvol(i,j) = tvol(1,j)
         end do
      end do  

      do j = 1, ny
         do i = 1, nx
            phio(i,j) = phi0*(noise+2d0*(1d0-noise)*rand(0))
            k(i,j) = k0*phio(i,j)**nphi
         end do 
      end do 

      phin = phio
      por = phio
c                                      time in years
      time = 0.0

      ddpmax = 0d0
      trvol = 0d0
      ttrvol = 0d0         
c                                       output conditions  
c                                       normalized porosity      
                open (200,file='./output/'//fname//'/f'//fname//'0  ')
                write (200,3000) ((phin(i,j),j=1,ny),i=1,nx)
                close (200) 
c                                       normalized pressure
c                                       1 = lithostatic
c                                       0 = hydrostatic
                open (200,file='./output/'//fname//'/p'//fname//'0  ')
                write (200,3000) ((pn(i,j)/p01(j),j=1,ny),i=1,nx) 
                close (200)
c                                       normalized water fertility
                open (200,file='./output/'//fname//'/h'//fname//'0  ')
                write (200,3000) ((r(i,j),j=1,ny),i=1,nx) 
                close (200)



c                                      plot parameters
c                                      output xmesh
      write (*,*) './output/'//fname//'/'//fname//'_x'
        
      open (200,file='./output/'//fname//'/'//fname//'_x')
      write (200,2000) (x0(i),i=1,nx)
      close (200)

      open (200,file='./output/'//fname//'/'//fname//'_y')
      write (200,2000) (y0(i),i=1,ny)
      close (200)

      open (200,file='./output/'//fname//'/'//fname//'_dx')
      write (200,2000) (dx,i=1,nx)
      close (200)

      open (200,file='./output/'//fname//'/'//fname//'_dy')
      write (200,2000) (dy,i=1,ny)
      close (200)

      open (200,file=fname)
      write (200,1000) nx,ny,0d0,0d0,0d0,0d0,odtav,
     *              0d0,0d0,0d0,0d0
      close (200)   
c                                      
      open (300,file='./output/'//fname//'/'//fname//'_fl_prod_inst')
      open (301,file='./output/'//fname//'/'//fname//'_fl_prod_cum')                            
c                                      solve for pressures (explicit adi):
2     if (ixy.eq.1) then
c                                      explicit in x direction
         call adixy (dts,rvol)
         ixy = 0
      else 
c                                      explicit in y direction
         call adiyx (dts,rvol)
         ixy = 1
      end if 
c                                      "instantaneous" fluid production 
      trvol = trvol + area*rvol*dts*v2
      tinst = tinst + dts/ys
c 
      dpmax = 0.0
      dtt = 1d99

      phimax = 0d0
      irxn = 0
      do i = 1, nx 
         do j = 1, ny
            if (phin(i,j).gt.phimax) phimax = phin(i,j)
            if (dvol(i,j).ne.0) irxn = j 
            if (c6*phin(i,j)/k(i,j).lt.dtt) dtt = c6*phin(i,j)/k(i,j)
            dp(i,j) = (pn(i,j) - p0(j))*1d-5
            ddp = (pn(i,j)-po(i,j))*1d-5

            if (dabs(ddp).gt.ddpmax) ddpmax = dabs(ddp)

            if (dabs(dp(i,j)).gt.dpmax) then
               im = i
               jm = j
               dpmax = dabs(dp(i,j))
            end if 
         end do 
      end do 


      tvol = tvol + dvol 
      rxn = rxn + dvol
      phio = phin
      po = pn
c                                 new plotting stuff

c                                 this is the instantaneous rate output
      if (tinst.gt.dtinst) then 

         fcum = fcum + trvol

         write (300,*) time, trvol/tinst
         write (301,*) time, fcum

         tinst = 0d0
         trvol = 0d0 

      end if 




c                                 this is the pressure/porosity plot section
      if (tplot.le.0d0) then 

                iplot = iplot + 1
                tplot = dtav 
                if (time.gt.t_end) goto 99 

                if (iplot.lt.10) then
                   write (cnum,'(i1,3x)') iplot
                else if (iplot.lt.100) then 
                   write (cnum,'(i2,2x)') iplot
                else if (iplot.lt.1000) then 
                   write (cnum,'(i3,1x)') iplot
                else 
                   write (*,*) ' too many plots'
                   goto 99
                end if


                do i = 1, nx
                   do j = 1, ny
                      if (j.ne.ny) then 
                         qv(i,j) = -k(i,j)/mu*(pn(i,j) - pn(i,j+1))/dy 
                      else 
                         qv(i,j) = -k(i,j)/mu*(pn(i,j-1) - pn(i,j))/dy 
                      end if 
                      if (i.ne.nx) then 
                         qh(i,j) = -k(i,j)/mu*(pn(i+1,j) - pn(i,j))/dx 
                      else 
                         qh(i,j) = -k(i,j)/mu*(pn(i,j) - pn(i-1,j))/dx 
                      end if 
                   end do 
                end do 
c                                       output conditions  
c                                       normalized porosity      
                open (200,file='./output/'//fname//'/f'//fname//cnum)
                write (200,3000) ((phin(i,j),j=1,ny),i=1,nx)
                close (200) 
c                                       normalized pressure
c                                       1 = lithostatic
c                                       0 = hydrostatic
                open (200,file='./output/'//fname//'/p'//fname//cnum)
                write (200,3000) ((pn(i,j)/p01(j),j=1,ny),i=1,nx) 
                close (200)
c                                       reaction rate
                open (200,file='./output/'//fname//'/r'//fname//cnum)
                write (200,3000) ((r(i,j),j=1,ny),i=1,nx) 
                close (200)
c                                       normalized water fertility
                open (200,file='./output/'//fname//'/h'//fname//cnum)
                write (200,3000) ((tvol(i,j)/dv,j=1,ny),i=1,nx) 
                close (200)
c                                       vertical flux
                open (200,file='./output/'//fname//'/qv'//fname//cnum)
                write (200,3000) ((qv(i,j),j=1,ny),i=1,nx) 
                close (200)
c                                       horizontal flux
                open (200,file='./output/'//fname//'/qh'//fname//cnum)
                write (200,3000) ((qh(i,j),j=1,ny),i=1,nx) 
                close (200)

      end if 
            
      tplot = tplot - dts/ys

      ty = ty + dtdty*dts/ys

      tictoc = tictoc + 1
      if (tictoc.gt.100) then 
         write (*,*) 'time, dts/ys',time, dts/ys, tplot, dtav, tcrky,
     *               trvol, iplot,time/tau
         tictoc = 0
      end if 

      if (time.ge.tcrky) then 
c                                        make the crack 
         imid = nxm1/2
         icrk = xcrk/dx 
         jcrk = ycrk/dy + 1

         do i = imid-icrk, imid+icrk
            do j = 1, jcrk
c                                        make the cracks k-phi = to the 
c                                        max value in the layer.
               if (phio(i,j).gt.phicrk) phicrk = phio(i,j)
c                                        set hcrk = 1 to make the crack
c                                        non reacting and non compacting.
               hcrk(i,j) = 1
            end do
         end do  

c       kcrk = k0*phicrk**nphi

         do i = imid-icrk, imid+icrk
            do j = 1, jcrk

               phio(i,j) = phicrk
               phin(i,j) = phicrk
               por(i,j) = phicrk
               k(i,j) = kcrk
               r(i,j) = 0d0
               dvol(i,j) = 0d0
               tvol(i,j) = dv

            end do
         end do  

         if (phi0/kcrk.lt.r1) r1 = phi0/kcrk
c                                      crack properties:
         odts = c6*r1*1d3
c
         dtav = odtav
         tplot = 0d0
         time = 0d0
         write (*,*) ' from now on the plot interval is (y):', dtav
         tcrky = 1d99

      end if 

      time = time + dts/ys
c
      dts = dts * 1.5
      if (dts.gt.odts) dts = odts

      goto 2

1000  format (g14.7)
1010  format (20(1x,f8.3)) 
1020  format (a132)
1030  format (a3,i1)
1040  format (a3,i2)
1050  format (8(a132,/),5(g12.6,1x),/,2(i4,1x),/,5(g12.6,1x))
1060  format (a200)
1070  format (20(1x,f8.1))
2000  format (g14.7)
2001  format (a3,a2)
3000  format (g21.14e3)
                                            
99    close (300)
      close (301)

      end

      subroutine adixy (dt,trvol) 
c--------------------------------------------------------------------------
c adi explicit in x direction and implicit in y, with isopotential 
c boundaries on y and no flow boundaries on x.
c--------------------------------------------------------------------------
      parameter (ix=601,iy=321,imx=601)

      implicit double precision (a-h,k,m,o-z)

      integer hcrk, icrust

      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1 
     *      / data2 /po(ix,iy),c0,c1,c2,c3,c4,betat,betap,rho,mu,dtdts
     *      / data3 /ty(iy),y0(iy),p0(iy),dpz(iy),dpf(iy)
     *      / data4 /phio(ix,iy),phin(ix,iy),tvol(ix,iy),dvol(ix,iy)
     *      / data6 /k0,por(ix,iy),s(ix,iy),k(ix,iy),nphi,hcrk(ix,iy)
     *      / rates /edot(ix,iy),r(ix,iy)
     *      / blok3 /erxn,t0,rnrxn,v1,v2,area,icrust

      double precision a(imx),c(imx),d(imx)

      odp = 99d99

1     dpmax = 0d0

      trvol = 0d0

      c0 = c4/dt
      c2 = c1/dt
      dtemp = dtdts*dt/2d0

      do i = 1, nx
         do j = 2, nym1
c                                      update properties:
c                                      get rates:
            if (hcrk(i,j).eq.0) then 
               p = (po(i,j) + pn(i,j))/2d0
               por(i,j) = (phio(i,j) + phin(i,j))/2d0
               tk = ty(j) + dtemp
c                                      reaction rate
               call rxnrat (tk,p+dpf(j),dt,tvol(i,j),
     *                   dvol(i,j),rvol)
c                                      deformation rate
c                                      explicit, correct when
c                                      iterative:
               if (j.gt.icrust) then 
                  call strrat (tk,p-dpz(j),por(i,j),dpor,edot(i,j),dt)
               else
                  edot(i,j) = 0d0
               end if
c                                      time rate of change in porosity
               phit = edot(i,j) + rvol
c                                      elastic component
               dphi = phit*dt + por(i,j)*betap*(pn(i,j)-po(i,j)) 

               phin(i,j) = phio(i,j) + dphi

               por(i,j) = phio(i,j) + dphi/2d0 

     	         k(i,j) = k0*por(i,j)**nphi

               r(i,j) = rvol/dt

               s(i,j) = (rvol*v2 - dphi/dt)/betat/por(i,j)

               trvol = trvol + rvol
            end if 
         end do
c                                      use a first order accurate
c                                      symmetry boundary on porosity
c                                      this could be second order
c        k(i,1)  = k(i,2)
c        phin(i,1) = phin(i,2)
c        phin(i,ny) = phin(i,nym1)
c        k(i,ny) = k(i,nym1)
      end do
c                             do columns first:
c                             constant potential boundaries,
c                             assume value constant along boundary rows
      d(1) = po(1,1)
      c(1) = 0d0

      do i = 1, nx

         ip = i + 1
         im = i - 1

         do j = 2, nym1
c                             compute a, b, c, d
            jp = j + 1
            jm = j - 1

            if (i.eq.1) then 
               im = ip
            else if (i.eq.nx) then
               ip = nxm1
            end if 

            dk = c0*por(i,j) + k(i,jp) + k(i,jm)
            a(j) = -k(i,jm)/dk
            c(j) = -k(i,jp)/dk
            d(j) = -(po(i,j)*(k(ip,j)+k(im,j)-c2*por(i,j))
     *             -k(im,j)*po(im,j)-k(ip,j)*po(ip,j)
     *             -s(i,j)*c1*por(i,j))*c3/dk

         end do
c                             solve for pn(i,2..ny-1)
         call tridiy (a,c,d,dpmax,i)
      end do

      if (dpmax.gt.1d3) then

        if (dpmax.gt.odp) then
           odp = 99.99d99   
           pn = po
           phin = phio       
           dt = dt / 2d0
           goto 1

        else
 
           odp = dpmax / 2d0

        end if

        goto 1

      end if 

      end

      subroutine adiyx (dt,trvol) 
c--------------------------------------------------------------------------
c adi explicit in y direction and implicit in x, with isopotential 
c boundaries on x and no flow boundaries on y.
c--------------------------------------------------------------------------
      parameter (ix=601,iy=321,imx=601)

      implicit double precision (a-h,k,m,o-z)    

      integer hcrk, icrust

      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1 
     *      / data2 /po(ix,iy),c0,c1,c2,c3,c4,betat,betap,rho,mu,dtdts
     *      / data3 /ty(iy),y0(iy),p0(iy),dpz(iy),dpf(iy)
     *      / data4 /phio(ix,iy),phin(ix,iy),tvol(ix,iy),dvol(ix,iy)
     *      / data6 /k0,por(ix,iy),s(ix,iy),k(ix,iy),nphi,hcrk(ix,iy)
     *      / rates /edot(ix,iy),r(ix,iy)
     *      / blok3 /erxn,t0,rnrxn,v1,v2,area,icrust

      double precision a(imx),c(imx),d(imx)

      odp = 99d99

1     dpmax = 0d0

      trvol = 0d0

      c0 = c4/dt
      c2 = c1/dt
      dtemp = dtdts*dt/2d0


      do i = 1, nx
         do j = 2, nym1
c                                      update properties:
c                                      get rates:
            if (hcrk(i,j).eq.0) then 
               tk = ty(j) + dtemp
               p = (po(i,j) + pn(i,j))/2d0
               por(i,j) = (phio(i,j) + phin(i,j))/2d0
c                                      reaction rate
               call rxnrat (tk,p+dpf(j),dt,tvol(i,j),
     *                   dvol(i,j),rvol)
c                                      deformation rate
c                                      explicit, correct when
c                                      iterative:
               if (j.gt.icrust) then 
                  call strrat (tk,p-dpz(j),por(i,j),dpor,edot(i,j),dt)
               else
                  edot(i,j) = 0d0
               end if 

               phit = edot(i,j) + rvol
c                                      elastic component
               dphi = phit*dt + por(i,j)*betap*(pn(i,j)-po(i,j)) 

               phin(i,j) = phio(i,j) + dphi

               por(i,j) = phio(i,j) + dphi/2d0 

     	         k(i,j) = k0*por(i,j)**nphi

               r(i,j) = rvol/dt

               s(i,j) = (rvol*v2 - dphi/dt)/betat/por(i,j)

               trvol = trvol + rvol
            end if 
         end do
c                                      use a first order accurate
c                                      symmetry boundary on porosity
c                                      this could be second order
c        k(i,1)  = k(i,2)
c        phin(i,1) = phin(i,2)
c        phin(i,ny) = phin(i,nym1)
c        k(i,ny) = k(i,nym1)
      end do
c                             do columns first:
c                             constant potential boundaries,
c                             assume value constant along boundary rows

      do j = 2, nym1

         jp = j + 1
         jm = j - 1

         do i = 2, nxm1
c                             compute a, b, c, d
            ip = i + 1
            im = i - 1

            dk = c2*por(i,j) + k(ip,j) + k(im,j)
            a(i) = -k(im,j)/dk
            c(i) = -k(ip,j)/dk
            d(i) = -(po(i,j)*(k(i,jp)+k(i,jm)-c0*por(i,j))
     *             -k(i,jm)*po(i,jm)-k(i,jp)*po(i,jp)
     *             -s(i,j)*c4*por(i,j))/c3/dk 

         end do
c                             set boundary coefficients:
         dk   = c2*por(1,j) + 2d0*k(2,j)
         c(1) = -2d0*k(2,j)/dk
         d(1) = -(po(1,j)*(k(1,jp)+k(1,jm)-c0*por(1,j))
     *            -k(1,jm)*po(1,jm)-k(1,jp)*po(1,jp)
     *            -s(1,j)*c4*por(1,j))/c3/dk

         dk   = c2*por(nx,j) + 2d0*k(nxm1,j)
         a(nx) = -2d0*k(nxm1,j)/dk
         d(nx) = -(po(nx,j)*(k(nx,jp)+k(nx,jm)-c0*por(nx,j))
     *            -k(nx,jm)*po(nx,jm)-k(nx,jp)*po(nx,jp)
     *            -s(nx,j)*c4*por(nx,j))/c3/dk

c                             solve for pn(i,2..ny-1)
         call tridix (a,c,d,dpmax,j)
      end do

      if (dpmax.gt.1d3) then

        if (dpmax.gt.odp) then
           odp = 99.99d99   
           pn = po
           phin = phio       
           dt = dt / 2d0
           goto 1

        else
 
           odp = dpmax / 2d0

        end if

        goto 1

      end if 

      end

      subroutine tridiy (a,c,d,dpmax,i)
c---------------------------------------------------------------
c tridi for constant potential boundaries, only runs over
c internal nodes. values of b(1),c(1),d(1) are never changed. 
c assumes values of p(1)=d(1) and p(ny) have been set.
c---------------------------------------------------------------
 
      parameter (ix=601,iy=321,imx=601)
     
      implicit double precision (a-h,k,m,o-z) 

      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1 

      double precision a(imx),b(imx),c(imx),d(imx)

      save b

      data b(1)/1d0/

      do j = 2, nym1
         jm = j - 1
         b(j) = 1d0 - a(j)*c(jm)/b(jm)
         d(j) = d(j) - a(j)*d(jm)/b(jm)
      end do

      do j = nym1, 2, -1
         pnew = (d(j) - c(j)*pn(i,j+1))/b(j)
         dp = dabs(pnew - pn(i,j))
         if (dp.gt.dpmax) dpmax = dp 
         pn(i,j) = pnew
      end do

      end 

      subroutine tridix (a,c,d,dpmax,j)
c---------------------------------------------------------------
c tridi for constant flow boundaries, only runs over
c boundary nodes.
c---------------------------------------------------------------
      parameter (ix=601,iy=321,imx=601)

      implicit double precision (a-h,k,m,o-z) 

      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1 

      double precision a(imx),b(imx),c(imx),d(imx)

      save b

      data b(1)/1d0/

      do i = 2, nx
         im = i - 1
         b(i) = 1d0 - a(i)*c(im)/b(im)
         d(i) = d(i) - a(i)*d(im)/b(im)
      end do

      pn(nx,j) = d(nx)/b(nx)

      do i = nxm1, 1, -1
         pnew = (d(i) - c(i)*pn(i+1,j))/b(i)
         dp = dabs(pnew - pn(i,j))
         if (dp.gt.dpmax) then 
            dpmax = dp 
            if (dpmax.gt.1d5) then
	         im = i
            end if 
         end if 
         pn(i,j) = pnew
      end do

      end 

      subroutine grxn (tk,p,dg)
c-----------------------------------------------------------------
      implicit double precision (a-h,k,n-z)  

      common/ blok1 /nh2o,vsol,fac,enth,vtot,entr

      dg = enth + p*vtot - entr*tk

      end

      subroutine equit (tk,p)
c-----------------------------------------------------------------
      implicit double precision (a-h,k,n-z)  

      common/ blok1 /nh2o,vsol,fac,enth,vtot,entr

      tk = (enth + p*vtot)/entr 

      end 

      subroutine rxnrat (tk,p,dts,tvol,dvol,rvol)
c------------------------------------------------------------------
c INPUT: 
c            tk   - temperature, kelvin
c            p    - pressure, pascal
c            dts  - time step, seconds
c            tvol - total amount of product that has
c                   been produced, m3/m3, cannot exceed dv.
c            denth - the random increment on enth.
c                   
c OUTPUT:
c            dvol - amount of RLM produced in time step
c                   m3/m3
c            rvol - porosity production m3/m3/s
c            rh2o - fluid production kg/m3/s
c            renth - enthalpy production kj/m3/s

c after each time step tvol should be incremented by dvol.
c-----------------------------------------------------------------
      implicit double precision (a-h,k,n-z) 

      integer icrust 

      common/ blok1 /nh2o,vsol,fac,enth,vtot,entr
     *      / blok2 /vand,rng,vm,dv
     *      / blok3 /erxn,t0,rnrxn,v1,v2,area,icrust
c
      rvol  = 0d0
      dvol  = 0d0
c                                     reactants gone, no prograde possible.
      if (tvol.ge.dv) goto 99
c                                     compute negative dg
      ndg = tk*entr - enth - p*vtot 
c                                     dg < 0, no dehydration possible,
c                                     this would have to be changed to
c                                     allow hydration.
      if (ndg.lt.1d0) goto 99
c                                     specific reaction rate
c                                     m3 RLM/s * dt
      dvol = fac*dexp(erxn/tk - t0)*(ndg)**rnrxn*dts
c                                     compare to maximum amt possible
      if (dvol.gt.dv - tvol) dvol = dv - tvol
c                                     volume m3 produced
      rvol = dvol * v1 / dts

99    end

      subroutine strrat (tk,dp,por,dpor,edot,dts)

      implicit double precision (a-h,k,o-z)  

      common / comp1 /ac,q,rn,rni,phim

      dpor = 0d0

      edot = 0d0

      ddp = dabs(dp)

      if (ddp.lt.1d-2) goto 99
c                                stage 2 compaction:
      edot = dp/ddp*ac*dexp(q/tk)*por*
     *       (1d0-por)*(ddp/(1d0-por**rni))**rn

      if (edot.gt.0d0) then 
         edot = 100.*edot
      else 
         edot = 0d0
      end if 

      dpor = edot * dts 

      if (por.le.phim) then
         if (dpor.lt.0d0) then
            dpor = 0.
            edot = 0.
         end if 
      end if 

      if (por+dpor.le.phim) then
         dpor = (por-phim)*0.5d0
         edot = -dpor/dts
      end if

99    end
