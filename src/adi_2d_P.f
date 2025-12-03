      subroutine adixy (dt,trvol) 
c--------------------------------------------------------------------------
c adi explicit in x direction and implicit in y, with isopotential 
c boundaries on y and no flow boundaries on x.
c--------------------------------------------------------------------------
      parameter (ix=601,iy=321,imx=601)

      implicit double precision (a-h,k,m,o-z)

      logical abort

      integer hcrk, icrust

      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1 
     *      / data2 /po(ix,iy),c0,c1,c2,c3,c4,betat,betap,rho,mu,dtdts
     *      / data3 /ty(iy),y0(iy),dpz(iy),dpf(iy)
     *      / data4 /phio(ix,iy),phin(ix,iy),tvol(ix,iy),dvol(ix,iy)
     *      / data6 /k0,por(ix,iy),s(ix,iy),k(ix,iy),nphi,hcrk(ix,iy)
     *      / rates /edot(ix,iy)
     *      / blok3 /erxn,t0,rnrxn,v1,v2,area,icrust

      double precision t,qy,qx,h,dx2,dy2,ktc,cfc,renth
      common/ data7 /t(ix,iy),qy(ix,iy),qx(ix,iy),h(ix,iy),dx2,dy2,
     *               ktc,cfc,renth

      double precision a(imx),c(imx),d(imx)

      odp = 99d99
      abort = .false.

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
c                                      reaction rate
               call rxnrat (t(i,j),p+dpf(j),dt,tvol(i,j),
     *                   dvol(i,j),rvol)
c                                      deformation rate
c                                      explicit, correct when
c                                      iterative:
               call strrat (t(i,j),p-dpz(j),por(i,j),dpor,edot(i,j),dt)
c                                      time rate of change in porosity
c                                      multiply rvol by v2 to emulate instantaneous 
c                                      cracking
               phit = edot(i,j) + rvol*v2
c                                      elastic component
               dphi = phit*dt + por(i,j)*betap*(pn(i,j)-po(i,j)) 

               phin(i,j) = phio(i,j) + dphi

               por(i,j) = phio(i,j) + dphi/2d0 

     	         k(i,j) = k0*por(i,j)**nphi

               h(i,j) = rvol*renth

               s(i,j) = (rvol*v2 - dphi/dt)/betat/por(i,j)

               trvol = trvol + rvol
            end if 
         end do
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
         call tridiy (a,c,d,dpmax,i,abort)

         if (abort) exit

      end do

      if (abort.or.dpmax.gt.1d3) then

        if (abort.or.dpmax.gt.odp) then
           odp = 99.99d99   
           pn = po
           phin = phio       
           dt = dt / 2d0
           abort = .false.

        else
 
           odp = dpmax / 2d0

        end if

        goto 1

      end if 
c                                 compute fluxes for heat flow
      call fluxes 

      end

      subroutine adiyx (dt,trvol) 
c--------------------------------------------------------------------------
c adi explicit in y direction and implicit in x, with isopotential 
c boundaries on x and no flow boundaries on y.
c--------------------------------------------------------------------------
      parameter (ix=601,iy=321,imx=601)

      implicit double precision (a-h,k,m,o-z)    

      integer hcrk, icrust

      logical abort

      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1 
     *      / data2 /po(ix,iy),c0,c1,c2,c3,c4,betat,betap,rho,mu,dtdts
     *      / data3 /ty(iy),y0(iy),dpz(iy),dpf(iy)
     *      / data4 /phio(ix,iy),phin(ix,iy),tvol(ix,iy),dvol(ix,iy)
     *      / data6 /k0,por(ix,iy),s(ix,iy),k(ix,iy),nphi,hcrk(ix,iy)
     *      / rates /edot(ix,iy)
     *      / blok3 /erxn,t0,rnrxn,v1,v2,area,icrust


      double precision t,qy,qx,h,dx2,dy2,ktc,cfc,renth
      common/ data7 /t(ix,iy),qy(ix,iy),qx(ix,iy),h(ix,iy),dx2,dy2,
     *               ktc,cfc,renth

      double precision a(imx),c(imx),d(imx)

      odp = 99d99
      abort = .false. 

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
c                                      reaction rate
               call rxnrat (t(i,j),p+dpf(j),dt,tvol(i,j),
     *                                dvol(i,j),rvol)
c                                      deformation rate
c                                      explicit, correct when
c                                      iterative:
               call strrat (t(i,j),p-dpz(j),por(i,j),dpor,edot(i,j),dt)
c                                      time rate of change in porosity
c                                      multiply rvol by v2 to emulate instantaneous 
c                                      cracking
               phit = edot(i,j) + rvol*v2
c                                      elastic component
               dphi = phit*dt + por(i,j)*betap*(pn(i,j)-po(i,j)) 

               phin(i,j) = phio(i,j) + dphi

               por(i,j) = phio(i,j) + dphi/2d0 

     	         k(i,j) = k0*por(i,j)**nphi

               h(i,j) = rvol/dt*renth

               s(i,j) = (rvol*v2 - dphi/dt)/betat/por(i,j)

               trvol = trvol + rvol

            end if 
         end do

      end do
c                             do columns first:
c                             constant potential boundaries,
c                             assume value constant along boundary rows

      do j = 2, nym1

         jp = j + 1
         jm = j - 1

         do i = 1, nx
c                             compute a, b, c, d
            ip = i + 1
            im = i - 1

            if (i.eq.nx) then 
               ip = im
            else if (i.eq.1) then 
               im = ip
            end if 

            dk = c2*por(i,j) + k(ip,j) + k(im,j)
            a(i) = -k(im,j)/dk
            c(i) = -k(ip,j)/dk
            d(i) = -(po(i,j)*(k(i,jp)+k(i,jm)-c0*por(i,j))
     *             -k(i,jm)*po(i,jm)-k(i,jp)*po(i,jp)
     *             -s(i,j)*c4*por(i,j))/c3/dk 

         end do

         a(nx) = a(nx) + c(nx)
         c(1)  = c(1)  + a(1) 
c                             solve for pn(i,2..ny-1)
         call tridix (a,c,d,dpmax,j,abort)

         if (abort) exit 

      end do

      if (abort.or.dpmax.gt.1d3) then

        if (abort.or.dpmax.gt.odp) then
           odp = 99.99d99   
           pn = po
           phin = phio       
           dt = dt / 2d0
           abort = .false.

        else
 
           odp = dpmax / 2d0

        end if

        goto 1

      end if 
c                                 compute fluxes for heat flow
      call fluxes 

      end
