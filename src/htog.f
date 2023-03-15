
          implicit real (a-z)
          tr = 298.15 
          k = 64.68
          na = 51.30
          ca = 41.63
          c  = 5.74
          cr = 25.543
          al = 28.35
          si = 18.81
          mg = 32.68
          fe = 27.28 
          o2 = 205.15
          h2 = 130.68
          c = 5.74 
 10       write (*,*)
     *     'enter nal, nsi, nca, nmg, nfe, nk, nna, no2, nh2, nc, ncr'
          read (*,*) nal, nsi, nca, nmg, nfe, nk, nna, no2,nh2, nc, ncr
          write (*,*) 'enter h (or g) and s for the phase:'
          read (*,*) h,s
          dsf = s - 
     *    (nal*al + nsi*si + nmg*mg + nfe*fe + nk*k + nca*ca +
     *     nna*na + no2*o2 + nh2*h2 + nc * c + ncr * cr)
          gf = h - tr * dsf
          hf = h + tr * dsf
          write (*,*) 'gf (if hf entered) =',gf
          write (*,*) 'hf (if gf entered) =',hf
          goto 10  
          end 
