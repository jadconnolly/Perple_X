      subroutine t2d_ady (dt) 
c--------------------------------------------------------------------------
c adi 2d solver for temperature, implicit in the y-direction
c mirror symmetry x boundaries, constant y boundaries
c upwind differencing
c--------------------------------------------------------------------------
      implicit none  

      integer i,ip,im,j,jp,jm,ix,iy,imx

      parameter (ix=601,iy=321,imx=601)

      double precision idt,dt,t0(ix,iy)

      integer ny,nx,nym1,nxm1
      double precision pn,dx,dy
      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1 

      double precision t,qy,qx,h,dx2,dy2,ktc,cfc,renth
      common/ data7 /t(ix,iy),qy(ix,iy),qx(ix,iy),h(ix,iy),dx2,dy2,
     *               ktc,cfc,renth

      double precision a(imx),b(imx),c(imx),d(imx),temp(imx)

      idt = 1d0/dt

      t0 = t

      do i = 1, nx

         ip = i + 1
         im = i - 1

         if (i.eq.1) then 
            im = ip
         else if (i.eq.nx) then 
            ip = im
         end if 

         do j = 2, nym1 

            jm = j - 1
            jp = j + 1

            if (qy(i,j).le.0d0.and.qx(i,j).le.0d0) then 

               a(j) = ktc / dy2
               b(j) = -2d0 * ktc / dy2 + cfc * qy(i,j) / dy - idt
               c(j) = ktc / dy2 - cfc * qy(i,jp) / dy
               d(j) = (qx(ip,j) * t0(ip,j) - qx(i,j) * t0(i,j)) / dx * 
     #cfc - h(i,j) + (-t0(ip,j) +2d0* t0(i,j) - t0(im,j)) * 
     #ktc / dx2 - t0(i,j) / dt

             else if (qy(i,j).le.0d0.and.qx(i,j).ge.0d0) then

               a(j) = ktc / dy2
               b(j)  = -2d0 * ktc / dy2 + cfc * qy(i,j) / dy - idt
               c(j) = ktc / dy2 - cfc * qy(i,jp) / dy
               d(j) = (qx(i,j) * t0(i,j) - qx(im,j) * t0(im,j)) / dx * 
     #cfc - h(i,j) + (-t0(ip,j) +2d0* t0(i,j) - t0(im,j)) * ktc /
     # dx2 - t0(i,j) / dt

             else if (qy(i,j).ge.0d0.and.qx(i,j).ge.0d0) then

               a(j) = ktc / dy2 + cfc * qy(i,jm) / dy
               b(j)  = -2d0 * ktc / dy2 - cfc * qy(i,j) / dy - idt
               c(j) = ktc / dy2
               d(j) = (qx(i,j) * t0(i,j) - qx(im,j) * t0(im,j)) / dx * 
     #cfc - h(i,j) + (-t0(ip,j) +2d0* t0(i,j) - t0(im,j)) * ktc /
     # dx2 - t0(i,j) / dt

             else if (qy(i,j).ge.0d0.and.qx(i,j).le.0d0) then

               a(j) = ktc / dy2 + cfc * qy(i,jm) / dy
               b(j)  = -2d0 * ktc / dy2 - cfc * qy(i,j) / dy - idt
               c(j) = ktc / dy2
               d(j) = (qx(ip,j) * t0(ip,j) - qx(i,j) * t0(i,j)) / dx * 
     #cfc - h(i,j) + (-t0(ip,j) +2d0* t0(i,j) - t0(im,j)) * 
     #ktc / dx2 - t0(i,j) / dt

             end if
         end do 
c                                 boundary conditions
         b(1)      = 1d0
         c(1)      = 0d0
         d(1)      = t0(i,1)
         a(ny)     = 0d0 
         b(ny)     = 1d0
         d(ny)     = t0(i,ny)
c                                 tridi solver 
         call tridi (temp,a,b,c,d,2,ny,nym1)       

         t(i,1:ny) = temp(1:ny)

      end do 

      end 

      subroutine t2d_adx (dt) 
c--------------------------------------------------------------------------
c adi 2d solver for temperature, implicit in the x-direction
c mirror symmetry x boundaries, constant y boundaries
c upwind differencing
c--------------------------------------------------------------------------
      implicit none  

      integer i,ip,im,j,jp,jm,ix,iy,imx

      parameter (ix=601,iy=321,imx=601)

      double precision idt,dt,t0(ix,iy),cx,cy

      integer ny,nx,nym1,nxm1
      double precision pn,dx,dy
      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1

      double precision t,qy,qx,h,dx2,dy2,ktc,cfc,renth
      common/ data7 /t(ix,iy),qy(ix,iy),qx(ix,iy),h(ix,iy),dx2,dy2,
     *               ktc,cfc,renth

      double precision a(imx),b(imx),c(imx),d(imx),temp(imx)

      idt = 1d0/dt

      t0 = t

      do j = 2, nym1

         jm = j - 1
         jp = j + 1

         do i = 1, nx 

            ip = i + 1
            im = i - 1

            if (i.eq.1) then 
               im = i+1
            else if (i.eq.nx) then 
               ip = i-1
            end if 

            cx = dabs(qx(i,j)*cfc*dt/dx)
            cy = dabs(qy(i,j)*cfc*dt/dy)
            if (cx.gt.1d0.or.cy.gt.1d0) then 
               write (*,*) 'oink oink ',qy(1,185),qy(2,185),
     *                     qy(3,185)
            end if 

            if (qy(i,j).le.0d0.and.qx(i,j).le.0d0) then 

               a(i) = ktc / dx2
               b(i) = -2d0 * ktc / dx2 + cfc * qx(i,j) / dx - idt
               c(i) = ktc / dx2 - cfc * qx(ip,j) / dx
               d(i) = (qy(i,jp) * t0(i,jp) - qy(i,j) * t0(i,j)) / dy * 
     #cfc - h(i,j) + (-t0(i,jp) +2d0* t0(i,j) - t0(i,jm)) *
     # ktc / dy2 - t0(i,j) / dt

             else if (qy(i,j).le.0d0.and.qx(i,j).ge.0d0) then

               a(i) = ktc / dx2 + cfc * qx(im,j) / dx
               b(i) = -2d0 * ktc / dx2 - cfc * qx(i,j) / dx - idt
               c(i) = ktc / dx2
               d(i) = (qy(i,jp) * t0(i,jp) - qy(i,j) * t0(i,j)) / dy *
     #cfc - h(i,j) + (-t0(i,jp) +2d0* t0(i,j) - t0(i,jm)) * 
     #ktc / dy2 - t0(i,j) / dt

             else if (qy(i,j).ge.0d0.and.qx(i,j).ge.0d0) then

               a(i) = ktc / dx2 + cfc * qx(im,j) / dx
               b(i) = -2d0 * ktc / dx2 - cfc * qx(i,j) / dx - idt
               c(i) = ktc / dx2
               d(i) = (qy(i,j) * t0(i,j) - qy(i,jm) * t0(i,jm)) / dy * 
     #cfc - h(i,j) + (-t0(i,jp) +2d0* t0(i,j) - t0(i,jm)) * ktc /
     # dy2 - t0(i,j) / dt

             else if (qy(i,j).ge.0d0.and.qx(i,j).le.0d0) then

               a(i) = ktc / dx2
               b(i) = -2d0 * ktc / dx2 + cfc * qx(i,j) / dx - idt
               c(i) = ktc / dx2 - cfc * qx(ip,j) / dx
               d(i) = (qy(i,j) * t0(i,j) - qy(i,jm) * t0(i,jm)) / dy * 
     #cfc - h(i,j) + (-t0(i,jp) +2d0* t0(i,j) - t0(i,jm)) * ktc /
     # dy2 - t0(i,j) / dt

             end if 

         end do 
c                                 boundary equations
         c(1) = c(1) + a(1)
         a(nx) = a(nx) + c(nx)
c                                  solve the system:
         call tridi (temp,a,b,c,d,1,nx,nx)

         t(1:nx,j) = temp(1:nx)

      end do 

      end 


      subroutine tridi (t,a,b,c,d,imin,imax,imax1)
c----------------------------------------------------------------------
c implicit tridi solver
c constant p=0 in y ix = 1, iy = ny, imin = 2, imax = ny, imax1 = ny - 1
c symmetry in x     ix = nx, iy = 1, imin = 1, imax = nx, imax1 = nx
c----------------------------------------------------------------------
      implicit none 
 
      integer imin, imax, imax1, j, jm, imx

      parameter (imx=601)

      double precision a(imx),b(imx),c(imx),d(imx),t(imx)

      do j = 2, imax1
         jm = j - 1
         b(j) = b(j) - a(j)*c(jm)/b(jm)
         d(j) = d(j) - a(j)*d(jm)/b(jm)
      end do

      t(1) = d(1)/b(1)
      t(imax) = d(imax)/b(imax)

      do j = imax-1, imin, -1
         t(j) = (d(j) - c(j)*t(j+1))/b(j)
      end do 

      end 

      subroutine fluxes

      implicit none 

      integer i,j,ix,iy,imx

      parameter (ix=601,iy=321,imx=601)

      integer ny,nx,nym1,nxm1
      double precision pn,dx,dy
      common/ data1 /pn(ix,iy),dx,dy,ny,nx,nym1,nxm1

      double precision t,qy,qx,h,dx2,dy2,ktc,cfc,renth
      common/ data7 /t(ix,iy),qy(ix,iy),qx(ix,iy),h(ix,iy),dx2,dy2,
     *               ktc,cfc,renth

      integer nphi, hcrk
      double precision k0,por,s,k
      common/ data6 /k0,por(ix,iy),s(ix,iy),k(ix,iy),nphi,hcrk(ix,iy)

      double precision po,c0,c1,c2,c3,c4,betat,betap,rho,mu,dtdts
      common/ data2 /po(ix,iy),c0,c1,c2,c3,c4,betat,betap,rho,mu,dtdts


      do i = 1, nx
         do j = 1, ny
            if (j.ne.ny.and.j.ne.1) then 
                qy(i,j) = -k(i,j)/mu*(pn(i,j+1) - pn(i,j-1))/2d0/dy 
            else if (j.eq.ny) then 
                qy(i,j) = -k(i,j)/mu*(pn(i,j) - pn(i,j-1))/dy
            else
                qy(i,j) = -k(i,j)/mu*(pn(i,j+1) - pn(i,j))/dy 
            end if 
            if (i.ne.nx.and.i.ne.1) then 
               qx(i,j) = -k(i,j)/mu*(pn(i+1,j) - pn(i-1,j))/2d0/dx 
            else 
               qx(i,j) = 0d0 
            end if 
         end do 
      end do 

      end 