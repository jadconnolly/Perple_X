      subroutine pscontor (cmin,ncon,dcon)
C----------------------------------------------------------------
C suncontour is a program to contour in x-y space.   
C----------------------------------------------------------------
C INPUT:
C
C NX, NY - number of x and y values, respectively.
C DX, DY - spacing of points in x and y directions, 
c          respectively.
C Z(I,J),I=1,NX,J=1,NY - value of the potential at each x-y 
C          point.
C----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer nseg,npts,npcs,mcon,i,np,j,ipts,istart,k,
     *         ncon,ipiece,iswit,iout

      double precision rline,cmin,dcon,thick,c

      parameter (nseg=100000,npts=250000,npcs=100000,mcon=50)

      character kontor*80, yes*1
    
      double precision cont(mcon),clinex(npts),
     *     cliney(npts),linex(npts),liney(npts),cline(2,npts),
     *     segs(4,nseg)

      integer ipieces(2,npcs),npiece(mcon),ifirst(mcon),
     *        next(nseg),ilast(mcon)

      integer ix,iy,mvar
      double precision z,zt 
      common/ dim   /z(nx,ny),ix,iy,mvar
      common/ dim1  /zt(nx,ny)

      double precision zmin,zmax
      common/ stuff /zmax,zmin

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      double precision xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
      common/ wsize /xmin,xmax,ymin,ymax,dcx,dcy,xlen,ylen
c----------------------------------------------------------------            
c                                  contor interval
      j = 0
      do i = 1, ncon
         c = cmin + (i-1) * dcon
         if (c.lt.zmin.or.c.gt.zmax) cycle
         j = j + 1
         cont(j) = cmin + (i-1) * dcon
      end do 

      if (j.eq.0) then 
         write (*,'(a)') 
     *   'no data within your contour limits, press enter to quit'
         read (*,'(a)') yes
         stop
      end if 

      ncon = j
      cmin = cont(1)

      do i = 1, ix 
         do j = 1, iy 
              zt(j,i)= z(i,j)
         end do
      end do 
c                                 write title info
      call pssctr (ifont,nscale,nscale,0d0)
      call pstext (xmin-2d0*dcx,ymax+15.5d0*dcy,title,162)
      write (kontor,1000) dcon, cmin, cont(ncon)  
      call pstext (xmin-2d0*dcx,ymax+12d0*dcy,kontor,0)
      write (kontor,1010) zmin, zmax
      call pstext (xmin-2d0*dcx,ymax+8.5d0*dcy,kontor,0)
      write (kontor,'(a)') 
     *     'Min/Max contours => thick solid/dotted curves'
      call pstext (xmin-2d0*dcx,ymax+5d0*dcy,kontor,0)

      call contra (xmin,xmax,ymin,
     *             ymax,ncon,cont,clinex,cliney,cline,segs,
     *             npts,nseg,npcs,ipieces,npiece,
     *             ifirst,next,ilast)

      write (*,1020) 
      read (*,'(a)') yes
      iout = 0
      if (yes.eq.'y'.or.yes.eq.'Y') iout = 1

      if (iout.eq.1) open (69,file='contor.dat')

      ipiece = 1
      iswit = 0
      do k = 1, ncon

         iswit = iabs (iswit-1)

         if (iswit.eq.0) then
            rline = 7d0
         else
            rline = 1d0
         end if

         if (k.eq.1.or.k.eq.ncon) then 

            thick = 2d0

            if (k.eq.1) then 
               rline = 1d0
            else
               rline = 9d0
            end if 

         else 
            thick = 0d0
         end if 

         if (iout.eq.1) write (69,*) 'contor: ',k

         np = npiece(k)
         if (np.gt.0) then 
            do i = 1, np
               ipts = ipieces(2,ipiece)
               if (iout.eq.1) write (69,*) 'segment: ',i
                
               if (ipts.ne.0) then 
                  do j = 1, ipts
                     istart = ipieces(1,ipiece)
                     linex(j) = clinex(istart+j-1)
                     liney(j) = cliney(istart+j-1)
                     if (iout.eq.1) write (69,*) liney(j),linex(j)
                  end do

                  call psbspl (linex,liney,ipts,rline,thick,0)
               end if 
               ipiece = ipiece + 1
            end do 
         end if
      end do 

1000  format ('contour interval: ',g10.4,'; range: ',g10.4,' => ',g10.4)
1010  format ('variable range: ',g10.4,' => ',g10.4)
1020  format ('Echo contour data to file contor.dat (Y/N)?')

      end     

      subroutine contra (xmin,xmax,ymin,ymax,ncon,cont,
     *                   clinex,cliney,cline,segs,npts,
     *                   nseg,npcs,ipieces,npiece,ifirst,next,
     *                   ilast)
c----------------------------------------------------------------------
c          cp      value to be contoured in a square array
c          ny     no of values in y direction
c          nx     no of values in x direction
c          xmin,xmax min and max values in x
c          ymin,ymax min and max value sin y
c          ncon    no of contopur levels.
c          cont    array containing contour levels.
c          clinex,cliney x,y pairs of points on contour segments.
c          ipieces array of pointers for contour segments.
c                  1,n start of coordinates in clinex cliney for
c                  current segment.
c                  2,n no of x,y pairs in current segment.
c          npiece  no of segments in each contour.
c
c       author: martin casey, 13-07-88

c       modified to use 4 data triangles, jadc 10-89
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ncon,npts,nseg,npcs

      logical iscon
  
      double precision c(3),x(3),y(3),cont(ncon),cline(2,npts),dx,dy,
     *                 clinex(npts),cliney(npts),segs(4,nseg),ymin,ymax,
     *                 xmin,xmax,dy2,dx2,xf,xn,yn,tol,yf,x1,x2,y2,y1,
     *                 xpc,ypc

      integer ipieces(2,npcs),npiece(ncon),ifirst(ncon),next(nseg),icf,
     *        ilast(ncon),i,j,k,nextg,numcon,iconlt,inow,ipiece,iseg,ic,
     *        lines,isegn,isego,icl,icf1

      integer ix,iy,mvar
      double precision z,zt 
      common/ dim   /z(nx,ny),ix,iy,mvar
      common/ dim1  /zt(nx,ny)


      dx = (xmax - xmin) / float(ix -1)
      dy = (ymax - ymin) / float(iy -1)

      do k = 1,ncon
         ifirst(k) = 0
         ilast(k) = 0
      end do 

      nextg = 0 

      do k = 1,nseg
         next(k) = 0
      end do

      dy2 = dy/2d0
      dx2 = dx/2d0
    
      do i = 1, iy-1
         do j = 1, ix-1
c                                            set the common vertex
          c(3) = (zt(i,j)+zt(i+1,j)+zt(i,j+1)+zt(i+1,j+1))/4.0
          x(3) = xmin + (float(j)-0.5)*dx
          y(3) = ymin + (float(i)-0.5)*dy
c                                            first triangle, (i,j) (i,j+1)
          c(1) = zt(i,j)
          c(2) = zt(i,j+1)
          x(1) = x(3)-dx2
          x(2) = x(3)+dx2
          y(1) = y(3)-dy2
          y(2) = y(1)

          call cfind (c,cont,ncon,iconlt,numcon,iscon)
          if (numcon .gt. 0) call cdraw (c,x,y,cont,ncon,iconlt,numcon,
     *                            nseg,segs,ifirst,next,ilast,nextg)     
c                                            second triangle, (i+1,j) (i+1,j+1)
          c(1) = zt(i+1,j)
          c(2) = zt(i+1,j+1)
          y(1) = y(3)+dy2
          y(2) = y(1)

          call cfind (c,cont,ncon,iconlt,numcon,iscon)
          if (numcon .gt. 0) call cdraw (c,x,y,cont,ncon,iconlt,numcon,
     *                            nseg,segs,ifirst,next,ilast,nextg)

c                                            third triangle, (i,j),(i+1,j)
          c(2) = c(1)
          c(1) = zt(i,j)
          x(2) = x(1)
          y(1) = y(3)-dy2 

          call cfind (c,cont,ncon,iconlt,numcon,iscon)
          if (numcon .gt. 0) call cdraw (c,x,y,cont,ncon,iconlt,numcon,
     *                            nseg,segs,ifirst,next,ilast,nextg)
c                                            fourth triangle, (i,j+1),(i+1,j+1)
          c(1) = zt(i,j+1)
          c(2) = zt(i+1,j+1)
          x(1) = x(3)+dx2
          x(2) = x(1)

          call cfind (c,cont,ncon,iconlt,numcon,iscon)
          if (numcon .gt. 0) call cdraw (c,x,y,cont,ncon,iconlt,numcon,
     *                            nseg,segs,ifirst,next,ilast,nextg)

         end do 
      end do 

      tol = 1d-10
      inow = 0
      ipiece = 1
      do 220 k = 1,ncon
      npiece(k) = 0
        if (ilast(k) .eq. 0) go to 220
        iseg = ifirst(k)
  145   continue
        ic = 2
        cline(1,1) = segs(1,iseg)
        cline(2,1) = segs(2,iseg)
        cline(1,2) = segs(3,iseg)
        cline(2,2) = segs(4,iseg)
        xf = segs(1,iseg)
        yf = segs(2,iseg)
        xn = segs(3,iseg)
        yn = segs(4,iseg)
        lines = 0
  150   continue
        isegn = next(iseg)
        lines = lines + 1
        if (iseg .eq. ilast(k)) go to 170
        do 160 i = 1, 1000000
c          if (isegn .eq. ilast(k)) go to 170
          x1 = segs(1,isegn)
          y1 = segs(2,isegn)
          x2 = segs(3,isegn)
          y2 = segs(4,isegn)
          if (((xn - x1)**2 + (yn - y1)**2) .lt. tol) then
            ic = ic + 1
            cline(1,ic) = x2
            cline(2,ic) = y2
            if (isegn .eq. ilast(k)) then
              ilast(k) = iseg
            else
              next(iseg) = next(isegn)
            endif
            iseg = ifirst(k)
            isegn = next(iseg)
            xn = x2
            yn = y2
          else if (((xn - x2)**2 + (yn - y2)**2) .lt. tol) then
              ic = ic + 1
              cline(1,ic) = x1
              cline(2,ic) = y1
              if (isegn .eq. ilast(k)) then
                ilast(k) = iseg
              else
                next(iseg) = next(isegn)
              endif
              iseg = ifirst(k)
              isegn = next(iseg)
              xn = x1
              yn = y1
            else
              iseg = isegn
              isegn = next(iseg)
          endif
          if (iseg .eq. ilast(k)) go to 170
  160   continue
  170   continue
        if (lines .eq. 2) go to 180
        if (((xn - xf)**2 + (yn - yf)**2) .lt. tol) go to 180
          icf = ic
          xn = xf
          yn = yf
          iseg = ifirst(k)
          go to 150
  180   continue
        iseg = ifirst(k)
        isego = iseg
        ifirst(k) = next(iseg)
        iseg = ifirst(k)
      npiece(k) = npiece(k) + 1

      if (ipiece.gt.npcs) call error (1,xn,npcs,'NPCS, CONTRA')

      ipieces(1,ipiece) = inow + 1

        if (lines .eq. 1) then
          do i = 1,ic
            xpc = cline(1,i)
            ypc = cline(2,i)
            inow = inow + 1
 
            if (inow.gt.npts) call error (1,xn,npts,'NPTS, CONTRA')
 
            clinex(inow) = xpc
            cliney(inow) = ypc
          end do 
          else
            icf1 = icf + 1
            icl = ic +1
            do i = icf1,ic
              inow = inow + 1

            if (inow.gt.npts) call error (1,xn,npts,'NPTS, CONTRA')

              icl = icl - 1
              xpc = cline(1,icl)
              ypc = cline(2,icl)
              clinex(inow) = xpc
              cliney(inow) = ypc
            end do 
            do i = 1,icf
              inow = inow + 1

            if (inow.gt.npts) call error (1,xn,npts,'NPTS, CONTRA')
              
              xpc = cline(1,i)
              ypc = cline(2,i)
              clinex(inow) = xpc
              cliney(inow) = ypc
             end do 
          end if
        ipieces(2,ipiece) = ic 
          ipiece = ipiece + 1
          if(isego .ne. ilast(k)) go to 145
  220 continue

      end

      subroutine cfind(c,cont,ncon,iconlt,numcon,iscon)
c
c      purpose
c
c       find if a contour crosses a data triangle
c
c      author:
c
c       martin casey 13-07-88
c
      implicit none 

      logical iscon

      integer i,k,ncon,iconlt,numcon,iconlo,iconhi

      double precision c(3),cont(ncon),cmin,cmax

      iscon = .false.
      numcon = 0
      cmax = -1d10
      cmin = 1d10
      do k = 1,3
         if(cmax .lt. c(k) ) cmax = c(k)
         if(cmin .gt. c(k) ) cmin = c(k)
      end do 
      if (cmax .eq. cmin ) return
      iconlo = 0

      do i = 1,ncon
         iconlo = iconlo + 1
         if(cont(iconlo) .ge. cmin) go to 120
      end do 

  120 if(cont(ncon) .lt. cmin) return

      iconhi = ncon

      do i = 1,ncon
        if(cont(iconhi) .lt. cmax) go to 140
        iconhi = iconhi - 1 
      end do

  140 if (cont(1) .gt. cmax) return
      if (iconhi .lt. iconlo) return
      iscon = .true.
      numcon = iconhi - iconlo + 1
      iconlt = iconlo

      end

      subroutine cdraw(c,x,y,cont,ncon,iconlt,numcon,nseg,
     *                 segs,ifirst,next,ilast,nextg)
c----------------------------------------------------------------
c      purpose:
c
c       draw contours across a data triangle
c
c      author:
c
c       martin casey 13-07-88
c-----------------------------------------------------------------
      implicit none

      integer ncon,iconlt,numcon,nseg,nextg,k,kmax,kmin,kx,kmid,ix

      double precision x(3),y(3),c(3),cont(ncon),segs(4,nseg),cmin,cmax,
     *                 delc,prop,x1,y1,x2,y2

      integer ifirst(ncon),next(nseg),ilast(ncon)

      cmax = -1d10
      cmin = 1d10

      kmax = 3
      kmin = 1

      do k = 1,3
        if(cmax .lt. c(k)) then
          cmax = c(k)
          kmax = k
        endif
        if(cmin .gt. c(k)) then
          cmin = c(k)
          kmin = k 
        endif
      end do

      do k = 1,3
        kmid = k
        if (k.ne.kmin.and.k.ne.kmax) go to 120
      end do

  120 k = iconlt - 1

      do kx = 1,numcon
        k = k + 1
        delc = cont(k) -cmin
        prop = delc /(cmax - cmin)
        x1 = x(kmin) + (x(kmax) - x(kmin)) * prop
        y1 = y(kmin) + (y(kmax) - y(kmin)) * prop

        if (c(kmid) .gt. cont(k)) then
          delc = cont(k) - cmin
          prop = delc / (c(kmid) - cmin)
          x2 = x(kmin) + (x(kmid) - x(kmin)) * prop
          y2 = y(kmin) + (y(kmid) - y(kmin)) * prop
        else
          delc = cont(k) - cmax
          prop = delc /(cmax - c(kmid))
          x2 = x(kmax) + (x(kmax) - x(kmid)) * prop
          y2 = y(kmax) + (y(kmax) - y(kmid)) * prop
        end if

        nextg = nextg + 1

        if (nextg.gt.nseg) call error (1,x2,nseg,'NSEG, CDRAW')

        if (ilast(k) .eq. 0 ) then
          ifirst (k) = nextg
          ilast(k) = nextg
        else
          ix = ilast(k)
          next(ix) = nextg
          ilast(k) = nextg
        endif

        segs(1,nextg) = x1
        segs(2,nextg) = y1
        segs(3,nextg) = x2
        segs(4,nextg) = y2
      end do 

      end