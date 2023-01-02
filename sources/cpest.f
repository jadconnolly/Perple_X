 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ier,idead

      logical nodata, bulk

      character amount*6, yes*1,char*1,fname*30

      integer itri(4),jtri(4),ijpt

      double precision wt(3), num,sl2,sl1,sg2,sg1,scrit,rhol1,
     * rhol2,rhog1,rhog2,rcrit,p2,p1,mt,msl,msg,mrhol,mrhog,bsl,bsg,
     * brhol,brhog,tcrit,pcrit,pavg,t01,t02

      double precision beta, rl, rv, lt, rhocl, rhocv, b, lpc, lp3,
     *                 f, df, tc, r1, r2, r3, t1, t2, t3, lp1, lp2

      integer iwt
      common/ cst209 /iwt

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k5),ctotal,jdv(k19),npt,fulrnk

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)

      integer iam
      common/ cst4 /iam

      double precision tlv, dt, rho, tst, tlv1, tlv2, rho1, rho2, cp1, 
     *                 cp2, x1, dp, rhoc, x2, pv1, pv2, pvv1, pvv2,
     *                 spec1(5,2),n(2),x(2),molwt(2),specwt(5),nat,
     *                 prps(7,2),xb(2),no,nsi,tot,p,t,lnk1,lnk2,lnk3,
     *                 lnk4,lnk5,ravg,savg

      data specwt/60.084, 44.085, 15.999, 31.998, 28.086/

      integer imax,imin,ir1,ir2,tic,j,lprops(7),lun
      logical go, quit
c                                 N, H, S, Cp, Cp/Cv, rho, vphi
      data lprops/17,2,15,12,28,10,7/
      integer idspec
      double precision spec
      common/ tspec /spec(nsp,k5),idspec

      common/ rcrt /rhoc

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision vp,vvp
      integer rooti
      common/ srkdiv /vp(3),vvp(3),rooti(3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------- 

      t1 = 0d0
      t2 = 0d0
      t3 = 0d0

      lun = 112
      write (*,*) 'fname'
      read (*,'(a)') fname

      open (lun,file=fname)

      read (lun,'(a)') char
      read (lun,*) char
      read (lun,*) i
      read (lun,'(a)') char
      read (lun,*) r
      read (lun,*) r
      read (lun,*) i
      read (lun,*) r
      read (lun,'(a)') char
c      write (lun,'(9a)') 'log(p) p(bar) T(K) dtL(K) dtG(K) xL xG ',
c     *    'nL(mol) nG(mol) NL(g-mol) NG(g-mol) NA(g-at) DH(J/kg) ',
c     *    'SL SG CpL CpG Cp/CvL Cp/CvG rhoL rhoG v_phiL v_phiG ',
c     *    'yL_SiO2 yL_SiO yL_O yL_O2 yL_Si ', 
c     *    'yG_SiO2 yG_SiO yG_O yG_O2 yG_Si ',
c     *    'pvL pvvL irL pvG pvvG irG ln(K4) ln(C4) ln(K5) ln(C5) ravg ',
c     *    'savg ssG rsG'

         rhol1 = 1d0
         rhog1 = 1d0 
         sl1 = 1d0
         sg1 = 1d0 
         t01 = 1d0 
         p1 = -5


         do 

         read (lun,*,end=99) 
     *  v(1),p,tlv,tlv1,tlv2,
     *  x, n, molwt, prps(1,1), 
     *  prps(2,2), 
     *  ((prps(j,i),i=1,2),j=3,7),
     *  ((spec1(j,i),j=1,5),i=1,2),
     *   pv1,pvv1,ir1,pv2,pvv2,ir2,
     *   lnk4,lnk4,lnk5,lnk5,ravg,savg,
     *   props(15,1),props(10,1)

         rhol2 = prps(6,1)
         rhog2 = prps(6,2)
         sl2 = prps(3,1)
         sg2 = prps(3,2)
         t02 = tlv
         p2 = v(1)

         if (v(1).eq.0d0) then 
c                         boiling point props
         write (*,1001) tlv, rhol2, sl2, rhog2, sg2
1001  format ('tb(i) = ',g14.7,/,'rlb(i) = ',g14.7,/,'slb(i) = ',g14.7,/
     *        'rlg(i) = ',g14.7,/,'slg(i) = ',g14.7,/)

         end if 

         if (p2.eq.p1) cycle 

         bsg = (-sg1*p2+p1*sg2)/(-p2+p1)
         bsl = (-sl1*p2+p1*sl2)/(-p2+p1)
         brhog = (-rhog1*p2+p1*rhog2)/(-p2+p1)
         brhol = (-rhol1*p2+p1*rhol2)/(-p2+p1)
                    
         msg = (-sg2+sg1)/(-p2+p1)
         msl = (-sl2+sl1)/(-p2+p1)
         mrhog = (-rhog2+rhog1)/(-p2+p1)
         mrhol = (-rhol2+rhol1)/(-p2+p1)
         mt = (t02-t01)/(p2-p1)

         ps = -(-bsg+bsl)/(-msg+msl)
         pr = -(-brhog+brhol)/(-mrhog+mrhol)

         scrit    = msl*ps + bsl
         rcrit    = mrhol*pr + brhol
         pavg  = (ps + pr)/2d0
         tcrit = t02 + mt*(pr-p2)
            write (*,*)  '============================================'
      write (*,'(a,12(g12.6,1x))') 
     *            ' scrit rhocrit ',scrit,rcrit,tcrit,pr,ps,pavg

         write (*,1000) pr, tcrit, 1d1**pr, rcrit, scrit

1000  format ('lgp(i) = ',g14.7,/,'t(i) = ',g14.7,/,'p(i) = ',g14.7,/
     *        'rho(i) = ',g14.7,/,'s(i) = ',g14.7,/)

c-----------------------------------------------------------------------------------
c                                critical exponent stuff
         rho2 = rhog2
         rho1 = rhol2

         r3 = dlog(rho1-rho2)
         t3 = tlv
         lp3 = v(1)

         if (t2.gt.0d0.and.t1.gt.0d0.and.t3.ne.t1.and.t2.ne.t1.and.
     *       t2.ne.t3) then

            dt = dabs (t3 - t1)
            tc = t3 + dt


            do 

               f = ((r3 - r2) / (dlog(tc - t3) - 
     *                           dlog(tc - t2)) - (r2 - r1) / 
     *                          (dlog(tc - t2) - dlog(tc - t1))) / 
     *             (dlog(tc - t3) / 2D0 - dlog(tc - t1)/ 2d0) 

               df = ((r3 - r2) * (t2 - t3) / (-tc + t3) / (-tc + t2) / 
     *            (log(tc -t3) - dlog(tc - t2)) ** 2 - (r2 - r1) * 
     *             (t1 - t2) / (-tc + t2) / (-tc + t1) / 
     *            (dlog(tc - t2) - dlog(tc - t1)) ** 2) / (dlog(tc - t3)
     *             / 0.2D1 - dlog(tc - t1) / 0.2D1) + ((r3 - r2) / 
     *            (dlog(tc - t3) - dlog(tc - t2)) + (r1 - r2) / 
     *            (dlog(tc - t2) - dlog(tc - t1))) * (t1 - t3) 
     *            / (-tc + t3) / (-tc + t1) / (dlog(tc - t3) / 0.2D1 
     *            - dlog(tc - t1) / 0.2D1) ** 2 / 0.2D1 

               tc = tc - f/df

               if (isnan(tc)) exit

               if (tc.lt.0d0.or.dabs(f/df).lt.1d-3) exit

            end do 

            if (.not.isnan(tc)) then 
            beta = (-r1 + r3) / (log(tc - t3) - log(tc - t1))
            b = (dlog(tc - t3) * r1 - r3 * dlog(tc - t1)) / 
     *          (dlog(tc - t3) - dlog(tc - t1)) 

            b = exp(b)/2d0

            rhocv = b * (tc - t3) ** beta + rho2
            rhocl = -b * (tc - t3) ** beta + rho1
c                                 crit pressure
      lpc = lp2 + ((lp3 - lp2) / (t3 - t2) - (lp2 - lp1) / (t2 - t1)) /(
     #t3 / 0.2D1 - t1 / 0.2D1) * (-t2 ** 2 + tc ** 2) / 0.2D1 - ((lp3 - 
     #lp2) / (t3 - t2) - (lp2 - lp1) / (t2 - t1)) / (t3 / 0.2D1 - t1 / 0
     #.2D1) * t2 * (tc - t2) + (lp3 - lp1) / (t3 - t1) * (tc - t2)


            write (*,*)  '---------------------------------------'
            write (*,*)  ' '
            write (*,*) 'at p, t, rhol/g:', p, t3, rho2, rho1
            write (*,*) 'scaled t:',(tc-t3)/tc
            write (*,*) 'beta tc',beta, tc 
            write (*,*) 'crit s',(sl2+sg2)/2d0, msl*lpc + bsl, 
     *                   msg*lpc + bsg
            write (*,*) 'crit rho',rhocv, mrhol*lpc + brhol
            write (*,*) 'crit p',1d1**lpc
c write an lv tab entry:

        x(2) = x(1)
        n(1) = (n(1) + n(2))/2d0
        n(2) = n(1)  
        molwt(1) = (molwt(1) + molwt(2))/2d0
        molwt(2) = molwt(1)    
        do j = 3, 7
           prps(j,1) = (prps(j,1) + prps(j,2))/2d0
           prps(j,2) = prps(j,1)
        end do 
        do j = 1, 2
        prps(3,j) = (sl2+sg2)/2d0
        prps(6,j) = rhocv
        end do
        do j = 1, 5
           spec1(j,1) = (spec1(j,1) + spec1(j,2))/2d0
           spec1(j,2) = spec1(j,1)
        end do              

            write (*,*) 'crit x',x(1)
            write (*,*) 'crit n, 1/n',n(1),1d0/n(1)
            write (*,*) 'crit N',molwt(1)
            write (*,*) 'crit z',1d1**lpc*molwt(1)/rhocv/8.314/tc*1d2


        write (*,'(/,a,/)') 'mock entry II'

        write (*,'(120(g14.7,1x))') 
     *  lpc,1d1**lpc,tc,0d0,tlv2,
     *  x, n, molwt, prps(1,1), 
     *  0d0, 
     *  ((prps(j,i),i=1,2),j=3,7),
     *  ((spec1(j,i),j=1,5),i=1,2),
     *   0d0,0d0,0,0d0,0d0,0,
     *   0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0

            write (*,*)  '---------------------------------------'
            end if 

         end if

         t1 = t2
         lp1 = lp2
         r1 = r2

         t2 = t3
         lp2 = lp3
         r2 = r3 
c------------------------------------------------------------------------------
         rhol1 = rhol2
         rhog1 = rhog2 
         sl1 = sl2 
         sg1 = sg2 
         t01 = t02 
         p1 = p2

            write (*,*)  '============================================'

         end do 
c                                     critical point prop
99      x(2) = x(1)
        n(1) = (n(1) + n(2))/2d0
        n(2) = n(1)  
        molwt(1) = (molwt(1) + molwt(2))/2d0
        molwt(2) = molwt(1)    
        do j = 3, 7
           prps(j,1) = (prps(j,1) + prps(j,2))/2d0
           prps(j,2) = prps(j,1)
        end do 
        do j = 1, 2
        prps(3,j) = scrit
        prps(6,j) = rcrit  
        end do
        do j = 1, 5
           spec1(j,1) = (spec1(j,1) + spec1(j,2))/2d0
           spec1(j,2) = spec1(j,1)
        end do              

         write (*,1000) pr, tcrit, 1d1**pr, rcrit, scrit

        write (*,'(120(g14.7,1x))') 
     *  pr,1d1**pr,tcrit,0d0,tlv2,
     *  x, n, molwt, prps(1,1), 
     *  0d0, 
     *  ((prps(j,i),i=1,2),j=3,7),
     *  ((spec1(j,i),j=1,5),i=1,2),
     *   0d0,0d0,0,0d0,0d0,0,
     *   0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0


      end 
    
   