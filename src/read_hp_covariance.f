c   program to: 

c 1) read the covariance matrix and phase names from  ds6*.txt TC format file
c 2) names and diagonal elements are then ouput to ds6*.dia
c 3) the elements can then be inserted into a perplex data file with the
c    read_666_data

c because ds636 added blank line separators the ds6*.txt file must by modified
c by adding "covariance" in front of the covariance matrix. 

      program trans

      implicit none

      include 'perplex_parameters.h'

      write (*,'(a)') 'enter file name root'
      read (*,'(a)') prject

      call mertxt (n1name,prject,'.txt',0)
      open (9, file=n1name, status= 'old')
      call mertxt (n2name,prject,'.dat',0)
      open (10,file=n2name)

      write (*,'(a,a)') 'writing output to :',n2name

      call rdin

      end

      subroutine rdin

      implicit none

      include 'perplex_parameters.h'

c   this subroutine reads part of the H&P data


      character     text(132)*1, name*8, cnum*80, gnum*80, twod*80 
      character     snum*80, vnum*80
      logical aq

      double precision rnum, comp(19), rgib, g, reas, s, reav,sfe,
     *                 tr,b1,b5,b6,b7,b8,dsf,atoms,ll4,ll5,ll6,
     *                 w1(300*300),
     *                 catoms(19),patoms,lam,kp,kpp,dkdt, zbig(300,300)
      double precision v, a, b, c, e, newa, k298, ll1, ll2, ll3,szr,smn

      double precision ctoj,sk,sna,sca,sc,sti,sal,ssi,smg,so,sh,scl,sni

      integer       index, igst, igen, isst, isen, ivst, iven, nphase

      integer       ihptov(19),i,inen,nst,j,inst,ilam,itype,ict,ier

c     hp ind      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 
      data ihptov/4,7,3,9,2,8,6,1,5,13,14,15,12,19,10,11,18,16,17/                2016
      data catoms/3.,2.,5.,3.,3.,2.,3.,2.,2.,2.,3.,2.,2.,3.,3.,2.,5.,2.,
     *            0./


1000  format(132a1)
1010  format(a,1x,g12.6,1x,i4)
2000  format(a8,4(i2),' H= ',g14.7)
2002  format(6(g14.7,1x))
2001  format(20(f5.2,1x))

      nphase = 0

c                                 read names and count 'til 'covaria'
      do 

         read (9,'(a7)') name


         if (name.eq.'covaria') exit

c                                 and 4 data records
         read (9,*)
         read (9,*)
         read (9,*)
         read (9,*)

         nphase = nphase + 1
         names(nphase) = name

      end do 
c                                 write the ouput to *.cov
         call mertxt (n2name,prject,'.cov',0)
         open (19,file=n2name)
c                          for some reason the covariance matrix 
c                          elements are preceded by a number, so read
c                          the extra number and offset the counter
         read (9,*,iostat=ier) w1(1:((nphase**2-nphase)/2+nphase+1))

         write (19,'(i4,1x,a)') nphase,'| number of entities'
         write (19,'(a)') names(1:nphase)
         write (19,*) w1(1:((nphase**2-nphase)/2+nphase))
         close (19)

        ict = 1

         call mertxt (n2name,prject,'.dia',0)
         open (19,file=n2name)

         do i = 1, nphase
            do j = i, nphase
               ict = ict + 1
               zbig(i,j) = w1(ict)
               if (i.eq.j) then 
                  write (*,*) i,names(i),zbig(i,i)
                  write (19,1010) names(i),zbig(i,i)*1d3,i
               end if
            end do
         end do

         close (19)

         stop

900     end