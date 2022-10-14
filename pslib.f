 
c PSLIB - a library of subroutines to generate postscript and idraw
c graphics primitives.
 
c Please do not distribute any part of this code.
 
c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy and
c Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c----------------------------------------------------------------
      subroutine pselip (xor,yor,dx,dy,rline,wide,ifill,ifg,ibg)

      implicit none

c pselip - subroutine to generate ellipse primitives.

      double precision rline,xor,yor,dy,dx,wide

      integer ifill,ixor,iyor,ifg,ibg

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1030)
 
      call psolin (rline,wide)
c     call psocfg (ifg,ibg)
      call psoclr
      call psofil (ifill)
      call psotrn
 
      call psscpt (xor,yor,ixor,iyor)

      write (nps,1020) ixor,iyor,int(dx * xscale),int(dy * yscale)
 
1020  format ('%I',/,4(i7,1x),' Elli',/,'End',/)
1030  format (/,'Begin %I Elli')
      end
c----------------------------------------------------------------
      subroutine psbspl (x,y,npts,rline,width,ifill)

      implicit none
 
c psbspl - subroutine to generate open bsplines.

      integer ifill,npts
 
      double precision width,x(npts),y(npts),rline

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1030)
 
      call psolin (rline,width)
      call psoclr
      call psofil (ifill)
      call psotrn
      call psopts (x,y,npts)
 
      write (nps,1020) npts
 
1020  format (i5,' BSpl',/,'End',/)
1030  format (/,'Begin %I BSpl')
      end
c----------------------------------------------------------------
      subroutine psrpgn (x1,y1,rx,ry,npts,rline,wide,ifill)
 
c psrpgn - subroutine to generate closed polygons, rel. coordinates.

      implicit none

      include 'perplex_parameters.h'

      integer ifill,npts,i,jpts
 
      double precision x(l5),y(l5),rx(npts),ry(npts),
     *                 x1,y1,rline,wide
 
      jpts = npts + 1
 
      if (jpts.gt.l5) call error (1,rline,l5,'L5 (PSRPGN)')
 
      x(1) = x1
      y(1) = y1
 
      do i = 2, jpts
         x(i) = x(i-1) + rx(i-1)
         y(i) = y(i-1) + ry(i-1)
      end do 
 
 
      call pspygn (x,y,jpts,rline,wide,ifill)
 
      end
c----------------------------------------------------------------
      subroutine pspygr (x,y,npts,rline,width,rfill)
 
c pspygr - subroutine to generate closed polygons, abs. coordinates.
c pspygr - with gray scale fill.

      implicit none

      integer npts
 
      double precision x(npts),y(npts),rline,width,rfill

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1030)
 
      call psolin (rline,width)
      call psoclr
      call psrfil (rfill)
      call psotrn
      call psopts (x,y,npts)
 
      write (nps,1020) npts
 
1020  format (i5,' Poly',/,'End',/)
1030  format (/,'Begin %I Poly')
      end
c----------------------------------------------------------------
      subroutine pspygn (x,y,npts,rline,width,ifill)
 
c pspygn - subroutine to generate closed polygons, abs. coordinates.

      implicit none

      integer npts,ifill
 
      double precision x(npts),y(npts),width,rline

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1030)
 
      call psolin (rline,width)
      call psoclr
      call psofil (ifill)
      call psotrn
      call psopts (x,y,npts)
 
      write (nps,1020) npts
 
1020  format (i5,' Poly',/,'End',/)
1030  format (/,'Begin %I Poly')
      end
c----------------------------------------------------------------
      subroutine pspyln (x,y,npts,rline,width,ifill)
 
c pspyln - subroutine to generate open polylines.

      implicit none

      integer npts,ifill
 
      double precision x(npts),y(npts),width,rline

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1030)
 
      call psolin (rline,width)
      call psoclr
      call psofil (ifill)
      call psotrn
      call psopts (x,y,npts)
 
      write (nps,1020) npts
 
1020  format (i5,' MLine',/,'End',/)
1030  format (/,'Begin %I MLine')
      end
c----------------------------------------------------------------
      subroutine pscspl (x,y,npts,rline,width,ifill)
 
c pscspl - subroutine to generate closed bsplines.

      implicit none

      integer npts,ifill
 
      double precision x(npts),y(npts),width,rline

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1030)
 
      call psolin (rline,width)
      call psoclr
      call psofil (ifill)
      call psotrn
      call psopts (x,y,npts)
 
      write (nps,1020) npts
 
1020  format (i5,' CBSpl',/,'End',/)
1030  format (/,'Begin %I CBSpl')
      end
c----------------------------------------------------------------
      subroutine psopts (x,y,npts)
 
c psopts - subroutine to output points.

      implicit none

      integer npts,i
 
      double precision x(npts),y(npts)

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,*) '%I ',npts
      write (nps,1010) ( int((x(i) - xmn) * xscale),
     *                   int((y(i) - ymn) * yscale),
     *                   i = 1, npts)
 
1010  format (10(i7,1x))
 
      end
c----------------------------------------------------------------
      subroutine psscpt (x,y,ix,iy)
 
c psscpt - subroutine to scale and fix real coordinates.

      implicit none

      integer ix,iy
      
      double precision x,y

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      ix = int((x - xmn) * xscale)
      iy = int((y - ymn) * yscale)
 
      end
c----------------------------------------------------------------
      subroutine psstrn (xs,ys,xt,yt,theta)
 
c psstrn - subroutine to set transformation matrix.

      implicit none

      double precision rad,theta,c0,s0,xs,ys,xt,yt

      double precision a,b,c,d,xxt,yyt
      common/ trans /a,b,c,d,xxt,yyt
c                                 RST transformation, theta is
c                                 the angle (degrees)
C                                 counterclockwise

      rad = theta*0.01745329251994d0

      c0 = dcos (rad)
      s0 = dsin (rad)
c                                 set lower limit for rotation
c                                 to 0.02 degree to avoid overflow
c                                 errors in adobe and coreldraw
      if (dabs(c0).lt.3.5d-4) c0 = 0d0
      if (dabs(s0).lt.3.5d-4) s0 = 0d0
 
      a = xs * c0
      b = ys * s0
      c = -xs * s0
      d = ys * c0
 
      xxt = xt
      yyt = yt
 
      end
c----------------------------------------------------------------
      subroutine psssc1 (ymin,ymax,xmin)
 
c psssc1 - subroutine to set scaling by y axis length

      implicit none

      double precision ymin,ymax,xmin

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      ymn = ymin
      xmn = xmin
      yscale = 3d3/(ymax-ymin)
      xscale = yscale
 
      end
c----------------------------------------------------------------
      subroutine psssc2 (xmin,xmax,ymin,ymax)

      implicit none
 
c psssc2 - subroutine to set scaling of axes independently

      double precision xmin,xmax,ymin,ymax

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      ymn = ymin
      xmn = xmin
      yscale = 3d3/(ymax-ymin)
      xscale = 3d3/(xmax-xmin)
 
      end
c----------------------------------------------------------------
      subroutine pssscm (xfac, yfac)
 
c pssscm - subroutine to modify scaling

      implicit none

      double precision xfac,yfac

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      yscale = yscale * yfac
      xscale = xscale * xfac
 
      end
c----------------------------------------------------------------
      subroutine psotrn
 
c psotrn - subroutine to output transformation matrix.

      implicit none

      double precision a,b,c,d,xt,yt
      common/ trans /a,b,c,d,xt,yt

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1000) a,b,c,d,xt,yt
 
1000  format ('%I t',/,'[',6(g9.3,1x),'] concat')
 
      end
c----------------------------------------------------------------
      subroutine psoclr
 
c psoclr - subroutine to output color choice.
 
      implicit none

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1000)
 
1000  format ('%I cfg Black',/,'0 0 0 SetCFg',/,'%I cbg White',/,
     *        '1 1 1 SetCBg')
 
      end
c----------------------------------------------------------------
      subroutine psored
 
c psoclr - subroutine to output red color.
 
      implicit none

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1000)
 
1000  format ('%I cfg Red',/,'1 0 0 SetCFg',/,'%I cbg Red',/,
     *        '1 0 0 SetCBg')
 
      end

      subroutine psocfg (ifg,ibg)
 
c psoclr - subroutine to output color fore/back-ground.
 
      implicit none

      integer nps, j, ifg, ibg

      real col(0:12,3)

      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps

      save col

      data col/
c                                 0 - black 
     *          0., 0., 0.,
c                                 1 - white 
     *          1., 1., 1.,
c                                 2 - red 
     *          1., 0., 0.,
c                                 3 - green
     *          0., 1., 0.,
c                                 4 - blue
     *          0., 0., 1.,
c                                 5 - purple
     *          0., 1., 1.,
c                                 6 - yellow
     *          1., 1., 0.,
c                                 7 - brown
     *          1., 0., 1.,
c                                 8 - orange
     *          1., .5, 0.,
c                                 9 - dark blue
     *          0., 0., .5,
c                                 10 - dark red
     *          .5, .0, .0,
c                                 11 - dark green
     *          0., .5, 0.,
c                                 12 - dark yellow
     *          0.5, 0.5, 0./
 
      write (nps,1000) (col(ifg,j),j=1,3),(col(ibg,j),j=1,3)
 
1000  format ('%I cfg Red',/,3(F3.1,1x),'SetCFg',/,'%I cbg Red',/,
     *        3(F3.1,1x),' SetCBg')
 
      end
c----------------------------------------------------------------
      subroutine psolin (rline,width)
 
c psolin - subroutine to ouput line header.
 
c width - width of line in points (i think).
c rline  - real line style indicator:

c      0 - no line
c      1 - solid line
c      2 - short dashed line
c      3 - uneven dashed line
c      4 - very sparse dotted line
c      5 - long even dashed line
c      6 - dashed line
c      7 - very short dashed line
c      8 - sparse dotted
c      9 - heavily dotted
c     10 - ultra sparse

c  for rline gt 9.0 relies on the idraw to generate a
c  line pattern from the value of rline which when converted
c  to hexadecimal gives the bit pattern for the line.

      implicit none

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps

      character plines(10)*28
 
      integer ilines(10),istyl,int

      double precision rline,width
 
      save ilines, plines
 
      data ilines/ 65535, 61680, 65080, 32896, 65520,
     *             65280, 52428, 34952, 43690, 32768/
 
      data plines/' 0 0 [] 0                 ',
     *            ' 0 0 [4 4 4 4] 17         ',
     *            ' 0 0 [7 3 3 3] 17         ',
     *            ' 0 0 [1 7 1 7] 17         ',
     *            ' 0 0 [12 4] 17            ',
     *            ' 0 0 [8 8] 17             ',
     *            ' 0 0 [2 2 2 2 2 2 2 2] 17 ',
     *            ' 0 0 [1 3 1 3 1 3 1 3] 17 ',
     *            ' 0 0 [1] 0                ',
     *            ' 0 0 [1 15] 17            '/
 
      istyl = int(rline)
 
      if (istyl.eq.0) then
 
         write (nps,1000)
 
      else if (istyl.gt.0.and.istyl.lt.11) then
 
         write (nps,1020) ilines(istyl),width,plines(istyl)
 
      else
 
          write (nps,1010) istyl, width
 
      end if
 
1000  format ('none SetB %I b n')
1010  format ('%I b ',i5,/,f5.2,' 0 0 [] 0 SetB')
1020  format ('%I b ',i5,/,f5.2,a28,'SetB')
 
      end
c----------------------------------------------------------------
      subroutine psrfil (rfill)
 
c psrfil - subroutine to output gray scale fill 

      implicit none

      double precision rfill

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps

      write (nps,1010) rfill
 
1010  format ('%I p',/,f6.4,' SetP')
 
      end

c----------------------------------------------------------------
      subroutine psofil (ifill)
 
c psofil - subroutine to output fill type

      implicit none

      integer ifill 

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      character*30 fill(15)
 
      save fill
c
c                               first seven fills vary in intensity
c                               from white to black
      data fill/'1','0.95','0.85','0.75','0.5','0.25','0',
c                               next 8 fills are patterned
     *          '< 11 22 44 88 11 22 44 88 > -1',
     *          '< 88 44 22 11 88 44 22 11 > -1',
     *          '< ff 00 00 00 ff 00 00 00 > -1',
     *          '< 88 88 88 88 88 88 88 88 > -1',
     *          '< ff 88 88 88 ff 88 88 88 > -1',
     *          '< 88 55 22 55 88 55 22 55 > -1',
     *          '< cc cc 33 33 cc cc 33 33 > -1',
     *          '< 77 bb ee dd 77 bb ee dd > -1'/
 
      if (ifill.eq.0) then
          write (nps,1000) 
      else if (ifill.le.15) then
         write (nps,1010) fill(ifill)
      else
         write (*,*) 'invalid fill choice'
         stop
      end if
 
1000  format ('none SetP %I p n')
1010  format ('%I p',/,a30,' SetP')
 
      end
 
c----------------------------------------------------------------
      subroutine psclos
 
c psclos subroutine to output epilogue and close LUN nps.

      implicit none

      character epips(4)*10

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps

      data epips/'End %I eop','showpage','%%Trailer','end'/

c                                     write epilog 
      write (nps,'(a)') epips
 
      close (nps)
 
      end
c----------------------------------------------------------------------
      subroutine psopen 
 
c psopen - subroutine to open LUN nps, write prologue

      implicit none

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer nps
      double precision xscale, yscale, xmn, ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
c----------------------------------------------------------------------
      nps = 50

      call mertxt (tfname,prject,'.ps',0)
 
      open (nps,file=tfname,status='unknown')
c                                 output eps prolog
      call psprol (nps)       
c                                 write output file name
      write (*,1010) tfname
 
1010  format (/,'PostScript will be written to file: ',a)
 
      end

c----------------------------------------------------------------
      subroutine psprol (nps)
 
c psprol - subroutine to write (EPS) postscript prolog.

      implicit none

      include 'perplex_parameters.h'

      integer i,nps

      character props(191)*63

      character font*40
      common/ myfont /font

      data (props(i),i=1,95)/
     *'%!PS-Adobe-2.0 EPSF-1.2',
     *'%%Pages: 1','%%EndComments',
     *'50 dict begin','/arrowHeight 8 def','/arrowWidth 4 def',
     *'/none null def','/numGraphicParameters 17 def',
     *'/stringLimit 65535 def',
     *'/Begin{save numGraphicParameters dict begin}def',
     *'/End{end restore}def',
     *'/SetB{dup type/nulltype eq{pop',
     *'false/brushRightArrow idef false/brushLeftArrow idef',
     *'true/brushNone idef}{/brushDashOffset idef',
     *'/brushDashArray idef 0 ne/brushRightArrow idef',
     *'0 ne/brushLeftArrow idef',
     *'/brushWidth idef false/brushNone idef}ifelse}def',
     *'/SetCFg{/fgblue idef/fggreen idef/fgred idef}def',
     *'/SetCBg{/bgblue idef/bggreen idef/bgred idef}def',
     *'/SetF{/printSize idef/printFont idef}def',
     *'/SetP{dup type/nulltype eq{pop true/patternNone idef}',
     *'{/patternGrayLevel idef',
     *'patternGrayLevel -1 eq{/patternString idef}if',
     *'false/patternNone idef}ifelse}def',
     *'/BSpl{0 begin storexyn newpath',
     *'n 1 gt{0 0 0 0 0 0 1 1 true subspline',
     *'n 2 gt{0 0 0 0 1 1 2 2 false subspline',
     *'1 1 n 3 sub{/i exch def',
     *'i 1 sub dup i dup i 1 add dup i 2 add dup false subspline}for',
     *'n 3 sub dup n 2 sub dup n 1 sub dup 2 copy false subspline}if',
     *'n 2 sub dup n 1 sub dup 2 copy 2 copy false subspline',
     *'patternNone not brushLeftArrow not brushRightArrow not and ',
     *'and{ifill}if',
     *'brushNone not{istroke}if 0 0 1 1 leftarrow',
     *'n 2 sub dup n 1 sub dup rightarrow}if',
     *'end}dup 0 4 dict put def','/Circ{newpath 0 360 arc ',
     *'patternNone not{ifill}if brushNone not{istroke}if}def',
     *'/CBSpl{0 begin dup 2 gt{storexyn newpath',
     *'n 1 sub dup 0 0 1 1 2 2 true subspline',
     *'1 1 n 3 sub{/i exch def',
     *'i 1 sub dup i dup i 1 add dup i 2 add dup false subspline}for',
     *'n 3 sub dup n 2 sub dup n 1 sub dup 0 0 false subspline',
     *'n 2 sub dup n 1 sub dup 0 0 1 1 false subspline',
     *'patternNone not{ifill}if',
     *'brushNone not{istroke}if}{Poly}ifelse',
     *'end}dup 0 4 dict put def',
     *'/Elli{0 begin newpath 4 2 roll translate scale',
     *'0 0 1 0 360 arc patternNone not{ifill}if',
     *'brushNone not{istroke}if end}dup 0 1 dict put def',
     *'/Line{0 begin 2 storexyn newpath',
     *'x 0 get y 0 get moveto x 1 get y 1 get lineto',
     *'brushNone not{istroke}if',
     *'0 0 1 1 leftarrow 0 0 1 1 rightarrow end}dup 0 4 dict put def',
     *'/MLine{0 begin storexyn newpath',
     *'n 1 gt{x 0 get y 0 get moveto 1 1 n 1 sub{/i exch def',
     *'x i get y i get lineto}for',
     *'patternNone not brushLeftArrow not brushRightArrow not and ',
     *'and{ifill}if','brushNone not{istroke}if',
     *'0 0 1 1 leftarrow n 2 sub dup n 1 sub dup rightarrow}if',
     *'end}dup 0 4 dict put def',
     *'/Poly{3 1 roll newpath moveto -1 add',
     *'{lineto}repeat closepath patternNone not{ifill}if',
     *'brushNone not{istroke}if}def',
     *'/Rect{0 begin/t exch def/r exch def/b exch def/l exch def',
     *'newpath l b moveto l t lineto r t lineto r b lineto closepath',
     *'patternNone not{ifill}if',
     *'brushNone not{istroke}if end}dup 0 4 dict put def',
     *'/Text{ishow}def',
     *'/idef{dup where{pop pop pop}{exch def}ifelse}def',
     *'/ifill{0 begin gsave',
     *'patternGrayLevel -1 ne{fgred bgred fgred sub patternGrayLevel ',
     *'mul add',
     *'fggreen bggreen fggreen sub patternGrayLevel mul add',
     *'fgblue bgblue fgblue sub patternGrayLevel mul add setrgbcolor',
     *'eofill}{eoclip originalCTM setmatrix',
     *'pathbbox/t exch def/r exch def/b exch def/l exch def',
     *'/w r l sub ceiling cvi def','/h t b sub ceiling cvi def',
     *'/imageByteWidth w 8 div ceiling cvi def',
     *'/imageHeight h def',
     *'bgred bggreen bgblue setrgbcolor eofill',
     *'fgred fggreen fgblue setrgbcolor',
     *'w 0 gt h 0 gt and{l b translate w h scale',
     *'w h true [w 0 0 h neg 0 h]{patternproc}imagemask}if}ifelse',
     *'grestore end}dup 0 8 dict put def',
     *'/istroke{gsave brushDashOffset -1 eq{[] 0 setdash',
     *'1 setgray}{brushDashArray brushDashOffset setdash',
     *'fgred fggreen fgblue setrgbcolor}ifelse',
     *'brushWidth setlinewidth originalCTM setmatrix',
     *'stroke grestore}def',
     *'/ishow{0 begin gsave fgred fggreen fgblue setrgbcolor',
     *'/fontDict printFont findfont printSize scalefont dup ',
     *'setfont def'/
      data (props(i),i=96,189)/
     *'/descender fontDict begin 0 [FontBBox] 1 get FontMatrix end',
     *'transform exch pop def',
     *'/vertoffset 0 descender sub printSize sub printFont/Courier ne',
     *'printFont/Courier-Bold ne and{1 add}if def{0 vertoffset moveto ',
     *'show /vertoffset vertoffset printSize sub def}forall',
     *'grestore end}dup 0 3 dict put def',
     *'/patternproc{0 begin',
     *'/patternByteLength patternString length def',
     *'/patternHeight patternByteLength 8 mul sqrt cvi def',
     *'/patternWidth patternHeight def',
     *'/patternByteWidth patternWidth 8 idiv def',
     *'/imageByteMaxLength imageByteWidth imageHeight mul',
     *'stringLimit patternByteWidth sub min def',
     *'/imageMaxHeight imageByteMaxLength imageByteWidth idiv ',
     *'patternHeight idiv',
     *'patternHeight mul patternHeight max def',
     *'/imageHeight imageHeight imageMaxHeight sub store',
     *'/imageString imageByteWidth imageMaxHeight mul ',
     *'patternByteWidth add string def',
     *'0 1 imageMaxHeight 1 sub{/y exch def',
     *'/patternRow y patternByteWidth mul patternByteLength mod def',
     *'/patternRowString patternString patternRow patternByteWidth ',
     *'getinterval def',
     *'/imageRow y imageByteWidth mul def',
     *'0 patternByteWidth imageByteWidth 1 sub{/x exch def',
     *'imageString imageRow x add patternRowString putinterval}for}',
     *'for imageString end}dup 0 12 dict put def',
     *'/min{dup 3 2 roll dup 4 3 roll lt{exch}if pop}def',
     *'/max{dup 3 2 roll dup 4 3 roll gt{exch}if pop}def',
     *'/arrowhead{0 begin transform originalCTM itransform',
     *'/taily exch def',
     *'/tailx exch def transform originalCTM itransform',
     *'/tipy exch def/tipx exch def',
     *'/dy tipy taily sub def/dx tipx tailx sub def',
     *'/angle dx 0 ne dy 0 ne or{dy dx atan}{90}ifelse def',
     *'gsave originalCTM setmatrix tipx tipy translate angle rotate',
     *'newpath 0 0 moveto arrowHeight neg arrowWidth 2 div lineto',
     *'arrowHeight neg arrowWidth 2 div neg lineto closepath',
     *'patternNone not{originalCTM setmatrix',
     *'/padtip arrowHeight 2 exp 0.25 arrowWidth 2 exp mul add sqrt ',
     *'brushWidth mul arrowWidth div def',
     *'/padtail brushWidth 2 div def',
     *'tipx tipy translate angle rotate padtip 0 translate',
     *'arrowHeight padtip add padtail add arrowHeight div dup scale',
     *'arrowheadpath ifill}if',
     *'brushNone not{originalCTM setmatrix',
     *'tipx tipy translate angle rotate arrowheadpath istroke}if',
     *'grestore end}dup 0 9 dict put def',
     *'/arrowheadpath{newpath 0 0 moveto',
     *'arrowHeight neg arrowWidth 2 div lineto',
     *'arrowHeight neg arrowWidth 2 div neg lineto','closepath}def',
     *'/leftarrow{0 begin y exch get/taily exch def',
     *'x exch get/tailx exch def y exch get/tipy exch def',
     *'x exch get/tipx exch def',
     *'brushLeftArrow{tipx tipy tailx taily arrowhead}if',
     *'end}dup 0 4 dict put def',
     *'/rightarrow{0 begin y exch get/tipy exch def',
     *'x exch get/tipx exch def y exch get/taily exch def',
     *'x exch get/tailx exch def',
     *'brushRightArrow{tipx tipy tailx taily arrowhead}if',
     *'end}dup 0 4 dict put def',
     *'/midpoint{0 begin/y1 exch def/x1 exch def/y0 exch def/x0 exch ',
     *'def x0 x1 add 2 div y0 y1 add 2 div end}dup 0 4 dict put def',
     *'/thirdpoint{0 begin/y1 exch def/x1 exch def/y0 exch def',
     *'/x0 exch def x0 2 mul x1 add 3 div y0 2 mul y1 add 3 div',
     *'end}dup 0 4 dict put def','/subspline{0 begin',
     *'/movetoNeeded exch def',
     *'y exch get/y3 exch def x exch get/x3 exch def',
     *'y exch get/y2 exch def x exch get/x2 exch def',
     *'y exch get/y1 exch def x exch get/x1 exch def',
     *'y exch get/y0 exch def x exch get/x0 exch def',
     *'x1 y1 x2 y2 thirdpoint/p1y exch def/p1x exch def',
     *'x2 y2 x1 y1 thirdpoint/p2y exch def/p2x exch def',
     *'x1 y1 x0 y0 thirdpoint p1x p1y midpoint',
     *'/p0y exch def/p0x exch def x2 y2 x3 y3 thirdpoint',
     *'p2x p2y midpoint/p3y exch def/p3x exch def',
     *'movetoNeeded{p0x p0y moveto}if p1x p1y p2x p2y p3x p3y curveto',
     *'end}dup 0 17 dict put def',
     *'/storexyn{/n exch def/y n array def/x n array def n 1 ',
     *'sub -1 0{/i exch def y i 3 2 roll put x i 3 2 roll put}for}def',
     *'%%EndProlog','%I Idraw 7 Grid 8','%%Page: 1 1','Begin',
     *'%I b u','%I cfg u','%I cbg u','%I f u','%I p u','%I t',
     *'[ .8 0 0 .8 0 0 ] concat',
     *'/originalCTM matrix currentmatrix def'/

      write (nps,'(a)') (props(i),i=1,2)
      write (nps,1000) font
      write (nps,1010) bbox
      write (nps,'(a)') (props(i),i=3,189)

1000  format ('%%IncludeFont: ',a)
1010  format ('%%BoundingBox: ',4(i4,1x))
      end 

c----------------------------------------------------------------
      subroutine pschct (nchar,ict,record)
 
c pschct - subroutine to count to first non-blank character.

      implicit none

      integer nchar,ict,i
 
      character*1 record(*)
 
      nchar = 1
 
      do i = 2, ict
         if ((record(i-1).eq.' ').and.(record(i).eq.' ')) exit
         nchar = nchar + 1
      end do 
 
      end
c----------------------------------------------------------------
      subroutine psrect (x1,x2,y1,y2,rline,width,ifill)
 
c psrect - subroutine to output a rectangle, with integer fill

      implicit none

      double precision x1,x2,y1,y2,rline,x(4),y(4),width

      integer ifill
 
      x(1) = x1
      x(2) = x1
      x(3) = x2
      x(4) = x2
      y(1) = y1
      y(2) = y2
      y(3) = y2
      y(4) = y1
 
      call pspygn (x,y,4,rline,width,ifill)
 
      end
c----------------------------------------------------------------
      subroutine psrecr (x1,x2,y1,y2,rline,width,rfill)
 
c psrecr - subroutine to output a rectangle, with real fill

      implicit none

      double precision x1,x2,y1,y2,rline,x(4),y(4),rfill,width
 
      x(1) = x1
      x(2) = x1
      x(3) = x2
      x(4) = x2
      y(1) = y1
      y(2) = y2
      y(3) = y2
      y(4) = y1
 
      call pspygr (x,y,4,rline,width,rfill)
 
      end
c----------------------------------------------------------------
      subroutine psrecb (x1,x2,y1,y2,rline,width)
 
c psrecb - subroutine to output a red rectangle for bad results in PSSECT

      implicit none

      double precision x1,x2,y1,y2,rline,x(4),y(4),width

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      x(1) = x1
      x(2) = x1
      x(3) = x2
      x(4) = x2
      y(1) = y1
      y(2) = y2
      y(3) = y2
      y(4) = y1

      write (nps,1030)
 
      call psolin (rline,width)
      call psored
      call psofil (1)
      call psotrn
      call psopts (x,y,4)
 
      write (nps,1020) 4
 
1020  format (i5,' Poly',/,'End',/)
1030  format (/,'Begin %I Poly')

      end
c------------------------------------------------------------------
      subroutine psline (x1,y1,x2,y2,rline,width)
 
c psline - subroutine to output a line (absolute).
 
      implicit none

      double precision x1,y1,x2,y2,rline,width

      integer nps
      double precision xscale, yscale, xmn, ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      write (nps,1000)
 
      call psolin (rline,width)
      call psoclr
 
      write (nps,1010)
 
      call psotrn
 
      write (nps,1020) int ((x1-xmn) * xscale), int ((y1-ymn) * yscale),
     *                 int ((x2-xmn) * xscale), int ((y2-ymn) * yscale)
 
1000  format ('Begin %I Line')
1010  format ('%I p',/,'0 SetP')
1020  format ('%I',/,4(I6,1x),' Line',/,'End',/)
 
      end
c------------------------------------------------------------------
      subroutine psmove (x1,y1)
 
c subroutine to make an absolute move.
 
      implicit none
      double precision x1,y1

      double precision x,y
      common/ posit /x,y

      x = x1
      y = y1
      end
c------------------------------------------------------------------
      subroutine psrmov (dx,dy)
 
c subroutine to make a relative move.
 
      implicit none

      double precision dx,dy

      double precision x,y
      common/ posit /x,y

      x = x + dx
      y = y + dy
      end
c------------------------------------------------------------------
      subroutine psrlin (dx,dy,rline,width)
 
c psrlin - subroutine to draw a relative line.
 
      implicit none
      
      double precision dx,dy,rline,width

      double precision x,y
      common/ posit /x,y
 
      call psline (x,y,x+dx,y+dy,rline,width)
 
      x = x + dx
      y = y + dy

      end
c------------------------------------------------------------------
      subroutine pswtod (x1,y1,x2,y2)
 
c pswtod - subroutine to do something i can't remember

      implicit none

      double precision x1,y1,x2,y2,x0,y0

      double precision a,b,c,d,xt,yt
      common/ trans /a,b,c,d,xt,yt

      integer nps
      double precision xscale, yscale, xmn, ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      x0 = ((x1-xmn) * xscale)
      y0 = ((y1-ymn) * yscale)
 
      x2 = x0 * a + y0 * c + xt
      y2 = x0 * b + y0 * d + yt
 
      end
c-------------------------------------------------------------------
      subroutine pssctr (kfont,xs,ys,theta)
 
c pssctr - subroutine to set character rotation-scale transformation
c          to be concatenated with the transformation set by psstrn

      implicit none

      include 'perplex_parameters.h'

      double precision x,y,xs,ys,theta,c0,s0,rad

      integer kfont

      integer jfont
      double precision a,b,c,d
      common / chars /a,b,c,d,jfont

      jfont = kfont

      rad = theta*0.01745329251994d0
 
      c0 = dcos (rad)
      s0 = dsin (rad)

      if (dabs(c0).lt.3.5d-4) c0 = 0d0
      if (dabs(s0).lt.3.5d-4) s0 = 0d0

      x = xs*cscale
      y = ys*cscale

      a = x * c0
      b = y * s0
      c = -x * s0
      d = y * c0
 
      end
c-------------------------------------------------------------------
      subroutine pstext (x,y,text,jchar)
 
c pstext - subroutine to output text strings.
 
c     x, y - data coordinates of the label
c     text - character string to be output
c     jchar - length of character string, 0 if unknown.

      implicit none

      double precision x,y,x0,y0,xtr,ytr

      integer jchar,nchar,i,ict
 
      character*(*) text
 
      character*1 itsy(400),ifonts(13)*33,bitsy(400)
 
      character font*40
      common/ myfont /font
      
      double precision a,b,c,d,xt,yt
      common/ trans /a,b,c,d,xt,yt

      integer ifont
      double precision ac,bc,cc,dc
      common/ chars /ac,bc,cc,dc,ifont

      integer nps
      double precision xscale,yscale,xmn,ymn
      common/ scales /xscale,yscale,xmn,ymn,nps
 
      save ifonts
 
c      data fonts/'/Times-Italic 14 SetF',
c     *           '/Times-Bold 14 SetF',
c     *           '/Times-Roman 14 SetF',
c     *           '/Times-Roman 12 SetF',
c     *           '/Helvetica-Oblique 14 SetF',
c     *           '/Helvetica-Bold 14 SetF',
c     *           '/Helvetica 14 SetF',
c     *           '/Helvetica 12 SetF',
c     *           '/Courier-Bold 12 SetF',
c     *           '/Courier 10 SetF',
c     *           '/Courier 8 SetF',
c     *           '/Symbol 12 SetF',
c     *           '/Symbol 14 SetF'/
 
      data ifonts/'%I f *-times-medium-i-*-140-*',
     *            '%I f *-times-bold-r-*-140-*',
     *            '%I f *-times-medium-r-*-140-*',
     *            '%I f *-times-medium-r-*-120-*',
     *            '%I f *-helvetica-medium-o-*-140-*',
     *            '%I f *-helvetica-bold-r-*-140-*',
     *            '%I f *-helvetica-medium-r-*-140-*',
     *            '%I f *-helvetica-medium-r-*-120-*',
     *            '%I f *-courier-bold-r-*-120-*',
     *            '%I f *-courier-medium-r-*-100-*',
     *            '%I f *-courier-medium-r-*-80-*',
     *            '%I f Symbol-12',
     *            '%I f Symbol-14'/
 
      if (jchar.eq.0) then
         nchar = len(text)
      else
         nchar = jchar
      end if
 
      if (nchar.gt.398) nchar = 398
 
      read (text,1020) (bitsy(i),i=2,nchar+1)
c                                 scan for '(' or ')'
      ict = 1
      do i = 2, nchar + 1
         if ((bitsy(i).eq.')').or.(bitsy(i).eq.'('))then
            ict = ict + 1
c                                 the double backslash
c                                 here is necessary for SUNS
c                                 because the first is read
c                                 as an escape character. On
c                                 convex this cause a compile
c                                 time warning, and on Mac and
c                                 IBM compilers there is no
c                                 warning.
            itsy(ict) = '\\'
            ict = ict + 1
            itsy(ict) = bitsy(i)
         else
            ict = ict + 1
            itsy(ict) = bitsy(i)
         end if
      end do 

      if (ict.gt.399) ict = 399
 
      itsy(1) = '('
      itsy(ict+1) = ')'
 
      x0 = (x-xmn) * xscale
      y0 = (y-ymn) * yscale
 
      xtr = x0 * a + y0 * c + xt
      ytr = x0 * b + y0 * d + yt
 
      write (nps,1010) ifonts(ifont),font,
     *                 ac,bc,cc,dc,xtr,ytr
 
      write (nps,1020) (itsy(i),i=1,ict+1)
      write (nps,1030)
 
1010  format ('Begin %I Text',/,
     *        '%I cfg Black',/,'0 0 0 SetCFg',/,a,/,
     *        '/',a,' 14 SetF',/,
     *        '%I t',/,'[',6(g9.3,1x),'] concat',/,
     *        '%I',/,'[')
1020  format (400a)
1030  format ('] Text',/,'End',/)
      end
c-------------------------------------------------------------------
      subroutine psublk (text,jchar)
 
c psublk - subroutine to remove double and leading blanks from text
 
c     text - character string 
c     jchar - length of unblanked character string, 0 if unknown.

      implicit none

      integer i,jchar,ict,ist
 
      character text*(*)
 
      character*1 itsy(255), bitsy(255) 
    
      if (jchar.eq.0) jchar = len(text)
 
      if (jchar.gt.255) jchar = 255
 
      read (text,1000) (bitsy(i), i = 1, jchar)
c                                 find first non-blank
      ist = 0 
      do i = 1, jchar
         if (bitsy(i).eq.' ') cycle
         ist = i 
         exit
      end do 

      if (ist.gt.0) then 
c                                 scan for double blanks:
        ict = 1
        itsy(1) = bitsy(i)

        do i = ist+1, jchar
           if (bitsy(i-1).ne.' '.or.bitsy(i).ne.' ') then
               ict = ict + 1
               itsy(ict) = bitsy(i)
            end if 
         end do 

         jchar = ict

         write (text,1000) (itsy(i), i = 1, jchar)

      else 

         text = ' '

      end if 
 
1000  format (400a)
      end

