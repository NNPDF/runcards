cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer nwgt
      character*(*) weights_info(*)

      call HwU_inithist(nwgt,weights_info)
      call HwU_book(1,'total rate      ', 1,0.5d0,1.5d0)

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
      integer istatus(nexternal)
      integer iPDG(nexternal)
      integer ibody
      integer i
      double precision p(0:4,nexternal)
      double precision wgts(*)
      double precision pttx(0:3), mtt, yt
      double precision dot
      double precision getabsy
      external getabsy
      double precision var

      do i=0,3
         pttx(i)=p(i,3)+p(i,4)
      enddo

      mtt = dsqrt(dot(pttx, pttx))
      yt = abs(getabsy(p(0,3), p(3,3)))

      var=1d0
      if(mtt.gt.340d0.and.mtt.lt.1500d0)then
         if(yt.lt.3.5d0)then
            call HwU_fill(1,var,wgts)
         endif
      endif

 999  return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function getabsy(en,pl)
      implicit none
      real*8 getabsy,en,pl,tmp
      tmp=pl/en
      if(abs(tmp).lt.1d0)then
        tmp=abs(atanh(tmp))
      else
        write(*,*)'Attempt to compute atanh(x) with x > 1'
        stop
      endif
      getabsy=tmp
      return
      end
