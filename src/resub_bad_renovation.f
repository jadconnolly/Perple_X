      subroutine resub (iter)
c----------------------------------------------------------------------
c subroutine to generate new pseudocompounds around each refinenent 
c point during adaptive optimization. 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical swap, bad, badsol

      integer i, ids, lds, id, kd, iter, idif, ifail, help

      double precision gg, gsol1

      external gsol1, badsol

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ikp
      common/ cst61 /ikp(k1)

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      character tname*10
      logical refine, lresub
      common/ cxt26 /refine,lresub,tname

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer igood(h9), ibad(h9)
      data igood, ibad/h9*0,h9*0/
      save igood, ibad
c----------------------------------------------------------------------
c                                 reset refinement point flags
      do i = 1, jpoint
         hkp(i) = 0
      end do

c     if (.not.lopt(59)) then
c                                 if ~keep_all_rpcs, then reset counters
c        jphct = jpoint
c        zcoct = 0
c     end if
c                                 loop on previous stable phases
c                                 refine as necessay:
      lds = 0
      lsdv(1:npt) = 0

      do kd = 1, npt

         if (iter.eq.1) then 
c                                 static array pointer is
            id = jdv(kd) + istct - 1
c                                 solution model pointer is
            ids = ikp(id)
c                                 refine if a solution
            if (ids.eq.0) cycle
c                                 reject if bad endmember
            if (badsol(ids)) cycle
c                                 get the refinement point composition
            if (id.gt.ipoint) then

               call setxyp (ids,id,bad)
c                                 save the composition for autorefine
               ststbl(id) = .true.

            else

               if (nrf(ids)) cycle
               call endpa (kd,id,ids)

            end if

            rkds = kd

         else
c                                 use pointer array lkp this uses 
c                                 negative values to index static
c                                 compounds and positive values to
c                                 point to solution models
            id = lkp(kd)
            rkds = mkp(kd)
c                                 save the location so that the 
c                                 amount can be initialized
            lsdv(kd) = jdv(kd)

            if (id.lt.0) then

               ids = ikp(-id)

               if (ids.eq.0) cycle
c                                 reject if bad endmember
               if (badsol(ids)) cycle

               if (nrf(ids)) cycle

               rkds = id
c                                 endmember refinement point:
               call endpa (kd,-id,ids)

            else

               ids = id
c                                 reject if bad endmember
               if (badsol(ids)) cycle
c                                 solution refinement point:
               call getpa (ids,kd)

            end if

         end if
c                                 set local model and refinement point pointers
         rids = ids
c                                 set solution model parameters for
c                                 gsol1, don't call if the previous
c                                 refinement point was the same solution.
         if (ids.ne.lds) then
            call ingsol (ids)
            if (deriv(ids)) call ingend (ids)
         end if

         if (.not.lopt(32).or.ksmod(ids).ne.39) then 
c                                 not lagged speciation, or not the electrolytic phase
            if (lopt(61)) call begtim (15)
c                                  normal solution
            call minfrc

            if (lopt(61)) call endtim (15,.false.,'minfrc')

         else
c                                 lagged speciation,call gsol1 get the initial
c                                 state (solute-free or pure solvent)
            gg = gsol1 (ids,.false.)
c                                 set kwak0 to record the state of the
c                                 refinement point
            kwak0 = rkwak

            if (nstot(ids).gt.1) then 
c                                 multispecies solvent
               if (lopt(61)) call begtim (15)
c                                 optimize solvent composition
               call minfrc

               if (lopt(61)) call endtim (15,.false.,'minfrc')

            else if (.not.rkwak) then

               idif = 0
c                                 save with 0-threshold
               call savkwk (gg,0d0,swap,idif)

            end if

         end if

         lds = ids

      end do

c     write (*,*) 'end of resub'

      end
