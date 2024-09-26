cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer nwgt
      character*(*) weights_info(*)

      call HwU_inithist(nwgt,weights_info)
      call HwU_book(1,'lmlp yrap',50, -1d0, 1d0)

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
      double precision xcos, getcostheta
      external getcostheta

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

      xcos=getcostheta(ppl(0),ppl(1),ppl(2),ppl(3),pplb(0),pplb(1),
     &                 pplb(2),pplb(3))

      call HwU_fill(1,xcos,wgts)
      end

      function getcostheta(enl,ptxl,ptyl,pzl,enal,ptxal,ptyal,pzal)
      implicit none
      real*8 getcostheta,enl,ptxl,ptyl,pzl,enal,ptxal,ptyal,pzal
      real*8 pzll,mll,p1p,p1m,p2p,p2m,pt2ll
c     implementation of first formula on page 6 of https://arxiv.org/abs/1710.05167
      p1p = enl + pzl
      p1m = enl - pzl
      p2p = enal + pzal
      p2m = enal - pzal
      pzll = pzl + pzal
      pt2ll = (ptxl + ptxal) * (ptxl + ptxal) +
     &        (ptyl + ptyal) * (ptyl + ptyal)
      mll = sqrt((enl + enal) * (enl + enal) - (pt2ll + pzll * pzll))
      getcostheta = (p1p*p2m-p1m*p2p)/sqrt(mll*mll+pt2ll)/mll*
     &              sign(1d0,pzll)

      if (abs(getcostheta) > 1d0) then
        print *, 'Calculating cos with cos > 1'
        stop
      endif
      end
