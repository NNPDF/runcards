#!/usr/bin/env python3

import os.path
import sys

cuts_code = {
    'mmllmax': '''c     cut for mmllmax (SFOS lepton pairs)
      do i=1,nexternal-1
        if (is_a_lm(i) .or. is_a_lp(i)) then
          do j=i+1,nexternal
            if (ipdg_reco(i) .eq. -ipdg_reco(j)) then
              if (invm2_04(p_reco(0,i),p_reco(0,j),1d0) .gt. {}**2) then
                passcuts_user=.false.
                return
              endif
            endif
          enddo
        endif
      enddo

''',
    'ptl1min': '''c     cut for ptl1min (leading lepton)
      j = 0
      do i=1,nexternal
        if (is_a_lm(i) .or. is_a_lp(i)) then
          if (j.eq.0 .or. pt_04(p_reco(0,i)).ge.pt_04(p_reco(0,j))) then
            j = i
          endif
        endif
      enddo
      if (pt_04(p_reco(0,j)) .lt. {}) then
        passcuts_user=.false.
        return
      endif

''',
    'yll': '''c     cut on the rapidity of the two leading leptons
      j = 0 ! leading lepton index
      mm = 0 ! subleading lepton index
      do i=1,nexternal
        if (is_a_lm(i) .or. is_a_lp(i)) then
          tmpvar = pt_04(p_reco(0,i))
          if (j.eq.0) then
            j = i
          else if (mm.eq.0) then
            if (tmpvar.ge.pt_04(p_reco(0,j))) then
              mm = j
              j = i
            else if (tmpvar.ge.pt_04(p_reco(0,mm))) then
              mm = i
            endif
          else if (tmpvar.ge.pt_04(p_reco(0,j))) then
            mm = j
            j = i
          else if (tmpvar.ge.pt_04(p_reco(0,mm))) then
            mm = i
          endif
        endif
      enddo
      if (abs(atanh((p_reco(3,j)+p_reco(3,mm))
     &    /(p_reco(0,j)+p_reco(0,mm)))) .gt. {}) then
        passcuts_user=.false.
        return
      endif

''',
    'ptzmin': '''c     cut on the pt of SFOS lepton pairs
      do i=1,nexternal-1
        if (is_a_lm(i) .or. is_a_lp(i)) then
          do j=i+1,nexternal
            if (ipdg_reco(i) .eq. -ipdg_reco(j)) then
              if (((p_reco(1,i)+p_reco(1,j))**2+
     &            (p_reco(2,i)+p_reco(2,j))**2) .lt. {}**2) then
                passcuts_user=.false.
                return
              endif
            endif
          enddo
        endif
      enddo

''',
    'ptzmax': '''c     cut on the pt of SFOS lepton pairs
      do i=1,nexternal-1
        if (is_a_lm(i) .or. is_a_lp(i)) then
          do j=i+1,nexternal
            if (ipdg_reco(i) .eq. -ipdg_reco(j)) then
              if (((p_reco(1,i)+p_reco(1,j))**2+
     &            (p_reco(2,i)+p_reco(2,j))**2) .gt. {}**2) then
                passcuts_user=.false.
                return
              endif
            endif
          enddo
        endif
      enddo

''',
    'yz': '''c     cut on the rapidity of SFOS lepton pairs
      do i=1,nexternal-1
        if (is_a_lm(i) .or. is_a_lp(i)) then
          do j=i+1,nexternal
            if (ipdg_reco(i) .eq. -ipdg_reco(j)) then
              if (abs(atanh((p_reco(3,i)+p_reco(3,j))
     &            /(p_reco(0,i)+p_reco(0,j)))) .gt. {}) then
                passcuts_user=.false.
                return
              endif
            endif
          enddo
        endif
      enddo

''',
}

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Error: wrong number of arguments: {}'.format(sys.argv))
        exit(1)

    # the name of the file we want to patch
    filename = sys.argv[1]

    if not os.path.exists(filename):
        print('Error: cut file `{}` does not exist'.format(filename))
        exit(2)

    for i in zip(sys.argv[2::3], sys.argv[3::3], sys.argv[4::3]):
        name = i[0]

        # check if the cut is recognised
        if not name in cuts_code:
            print('Error: unrecognised cut: {}'.format(name))
            exit(3)

        if i[1] != '=':
            print('Error: wrong argument format')
            exit(4)

    with open(filename, 'r') as fd:
        contents = fd.readlines()

    insertion_marker = 'USER-DEFINED CUTS'
    marker_pos = -1

    for lineno, value in enumerate(contents):
        if insertion_marker in value:
            marker_pos = lineno
            break

    if marker_pos == -1:
        print('Error: could not find insertion marker `{}` in cut file `{}`'
            .format(insertion_marker, filename))
        exit(5)

    # skip some lines with comments
    marker_pos = marker_pos + 4
    # insert and empty line
    contents.insert(marker_pos - 1, '\n')

    for i in zip(reversed(sys.argv[2::3]), reversed(sys.argv[4::3])):
        name = i[0]
        value = i[1]
        code = cuts_code[name].format(float(value))
        contents.insert(marker_pos, code)

    with open(filename, 'w') as fd:
        fd.writelines(contents)
