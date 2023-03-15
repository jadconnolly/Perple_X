c-----------------------------------------------------------------------
      program GEOPT

      implicit double precision (a-g,o-y),integer (h-m,z)

      character*8 uname,y*1,min(90),color,x,dbase*40

c-----------------------------------------------------------------------
      write (*,1040) 
      read (*,1000) uname 
      write (*,1050) uname
      read (*,1020) y
      if (y.ne.'y'.and.y.ne.'Y') then
         write (*,1060) uname
         goto 99
      end if 
      write (*,*) 'enter the name of the thermodynamic data base: '
      read (*,1030) dbase
      write (*,*) uname,' enter an x, if you also want to estimate ',
     *            'XCO2: '
      read (*,1020) x
      write (*,*) uname,' how many minerals in your rock?'
      read (*,*) ict
      do 10 i = 1, ict 
         write (*,*) 'Enter name of mineral ',i
10       read (*,1000) min(i)
      write (*,*) uname,' what color is your rock? '
      read (*,1000) color
      write (*,*) ' Enter an estimate for the',
     *            ' Fe(2+)/(Fe(2+)+Fe(3+)) ratio?'
      read (*,*) fefe
      write (*,*) uname,' enter a rough guess of the P(b) and T(C): '
      read (*,*) p,t
      write (*,*) 'what are the uncertainties for your P and T ',
     *             'estimates?'
      read (*,*) dp,dt
      if (x.eq.'x') then 
         write (*,*) uname,' enter a rough guess of the XCO2: '
         read (*,*) XCO2
         write (*,*) 'what is the uncertainty for your XCO2 estimate?'
         read (*,*) dx
         write (*,2450)
         read (*,*) ifug
      end if 

      dx = -dx/3.7
      dt = dt/1.5
      dp = -dp/2.1
      fo2 = -15.3 - fefe*9.

      write (*,*) 'ok ',uname,' i am gonna do my best, but it may ',
     *            'take some time '
      itic = 0
      jtic = 0 

      do while (jtic.lt.100) 
         itic = itic + 1
         if (itic.eq.25000000) then 
            write (*,*)
            write (*,*) 'working....'
            write (*,*) 
            itic = 0
            jtic = jtic + 1
         end if
      end do 

      write (*,*) uname,' a miracle has occurred! your rock works'
      write (*,1070) color   
      do  30 i = 1, ict 
         write (*,1080) min(i) 
30    continue 
      write (*,1090) p+dp,t+dt,fo2
      if (x.eq.'x') write (*,1100) xco2+dx
      write (*,1110) dbase
    
1000  format (a8)
1010  format (4(1x,i1))
1020  format (a1)
1030  format (a40)
1040  format (/,' Hi! i am a UNIVERSAL GEOTHERMOBAROMETER ',
     *         ' program.',/, ' what is your name? ')
1050  format (/,' Hi ',a8,'!, i like you and i hope we have fun!',/,
     *          ' can you say "therm-o-dy-nam-ics" (y/n)? ')
1060  format (/,' Weeell ', a8,' why dont you go read gibbs and try me',
     *          ' again later.')
1070  format (/,' Your ',a8,' colored rock containing the minerals:',/)
1080  format (12x,a8)
1090  format (/,' is estimated to have equilbrated at:',/,
     *        /,'         P(b)    = ',g13.6,/
     *         ,'         T(C)    = ',g13.6,/
     *         ,'         log fO2 = ',g13.6)
1100  format (  '         XCO2    = ',g13.6,/) 
1110  format (/,' based on calculations with the ',a40,/,' data base ',
     *        'and Markov Chain theory.',/)
2450  format (/,' Enter the number identifying the equation of state', 
     *        ' for the fluid phase',/
     *        '  (1) modified redlich kwong (mrk)',/,
     *        '  (2) hard sphere mrk (hsmrk)',/,
     *        '  (3) hybrid mrk-hsmrk',/,
     *        '  (4) saxena and fei 1987 virial expansion',/,
     *        '  (5) bottinga and richet redlich kwong 1982',/,
     *        '  (6) holland and powell 1990',/,
     *        '  (7) hybrid haar-hsmrk (trkmrk)',/)
99    pause

      end 
