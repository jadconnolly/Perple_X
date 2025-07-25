      integer nvars, ivar, ipts, i
      double precision v(20), ptamin(20) 
      character*60 filename
      write (*,*) 'filename'
      read (*,*) filename

      open (10,file=filename,status='old')
      open (11,file='filtered.pts',status='unknown')

      ipts = 0
      ptamin = 1d99
      nvars = 6

      do 
         read (10,*,end= 90) ivar,(v(i),i=1,nvars)

         if (v(2).ge.573.and.v(2).le.973.and.
     *       v(1).ge.1e4.and.v(1).le.3e4) then
            ipts = ipts + 1
            write (11,1000) ivar, v(1:nvars)
            if (v(nvars).lt.ptamin(nvars).and.ivar.ne.1) then
               ptamin(1:nvars) = v(1:nvars)
            end if
         end if

      end do

90    write (11,1000) 1, (ptamin(i),i=1,nvars)
      write (*,*) ipts
      write (*,1000) 1, (ptamin(i),i=1,nvars)

1000  format (i1,20(2x,g14.6))

      close (10)
      close (11)

      end

c VK30eclm_grt_cpx_central.pts
c VK30eclm_grt_cpx_perturbed.pts