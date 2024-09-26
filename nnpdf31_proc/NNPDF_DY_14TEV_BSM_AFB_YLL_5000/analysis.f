cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer nwgt
      character*(*) weights_info(*)

      call HwU_inithist(nwgt,weights_info)
      call HwU_book(1,'lmlp yrap',25, 0d0, 1d0)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(dummy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      double precision dummy
      call HwU_write_file
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'nexternal.inc'
      include 'cuts.inc'
      integer istatus(nexternal)
      integer iPDG(nexternal)
      integer ibody
      integer i
      integer j
      double precision p(0:4,nexternal)
      double precision wgts(*)
      double precision ppl(0:3), pplb(0:3), ppv(0:3)
      logical is_nextph_iso(nexternal),is_nextph_iso_reco(nexternal)
      double precision xyll, getabsy
      external getabsy

      double precision p_reco(0:4,nexternal)
      integer iPDG_reco(nexternal)



      is_nextph_iso(:) = .false.
      call recombine_momenta(rphreco, etaphreco, lepphreco, quarkphreco,
     $                       p, iPDG, is_nextph_iso, p_reco, iPDG_reco,
     $                       is_nextph_iso_reco)

      do j = nincoming+1, nexternal
        if (iPDG_reco(j).eq.13) ppl(0:3)=p_reco(0:3,j)
        if (iPDG_reco(j).eq.-13) pplb(0:3)=p_reco(0:3,j)
      enddo
      do i=0,3
        ppv(i)=ppl(i)+pplb(i)
      enddo

      xyll=getabsy(ppv(0),ppv(3))

      !if (xyll > 1.0d0) then
      !  write (*,*) "error: event outside bin", xyll
      !  stop 1
      !endif

      call HwU_fill(1,xyll,wgts)
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function getabsy(en,pl)
      implicit none
      real*8 getabsy,en,pl,tmp
c
      tmp=pl/en
      if (abs(tmp) < 1d0) then
        tmp=abs(atanh(tmp))
      else
        print *, 'Attempt to compute atanh(x) with x > 1'
        stop
      endif
      getabsy=tmp
      return
      end
