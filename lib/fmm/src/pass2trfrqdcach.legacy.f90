#include "fmm.h"

      subroutine pass2trfrqdcach(nmultipoles,nsqmultipoles,jaddress,&
      jacc,romega1,iomega1,romega2,iomega2,romega1a,iomega1a,romega2a,&
      iomega2a,romega1b,iomega1b,romega2b,iomega2b,romega1c,iomega1c,&
      romega2c,iomega2c,mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11,&
      mu12,mu13,mu14,mu15,mu16,cmphi,smphi,cmphipi,smphipi,csmphi,&
      csmphipi,sg,fr,d2,d3f,scr1,scr2)
!
      use fmmkinds
!
      implicit none
!
      real(kind=fmm_real) romega1(*),iomega1(*),romega2(*),iomega2(*),&
      romega1a(*),iomega1a(*),romega2a(*),iomega2a(*),romega1b(*),&
      iomega1b(*),romega2b(*),iomega2b(*),romega1c(*),iomega1c(*),&
      romega2c(*),iomega2c(*),mu1(*),mu2(*),mu3(*),mu4(*),mu5(*),&
      mu6(*),mu7(*),mu8(*),mu9(*),mu10(*),mu11(*),mu12(*),mu13(*),&
      mu14(*),mu15(*),mu16(*),cmphi(0:*),smphi(0:*),cmphipi(0:*),&
      smphipi(0:*),csmphi(*),csmphipi(*),sg(0:*),fr(0:*),d2(*),d3f(*),&
      scr1(*),scr2(*),a,aa,aaa,aaaa,b,bb,bbb,bbbb,c,cc,ccc,cccc,d,dd,&
      ddd,dddd,gl,g,glm
!
      integer(kind=fmm_integer) nmultipoles,nsqmultipoles,jaddress(*),i,&
      j,k,l,m,mmmm,mmm,n,nn,mm,mmmmm,mmmmmm,nnn
!
      logical(kind=fmm_logical) jacc(*)
!
      real(kind=fmm_real) zero
      parameter(zero=0.e0_fmm_real)
      real(kind=fmm_real) one
      parameter(one=1.e0_fmm_real)
!
!      call f_hpm_start_i8_(3,__LINE__,__FILE__,'rot')
!

!      call alignx(16,d2(1))
!      call alignx(16,d3f(1))
!      call alignx(16,scr1(1))
!      call alignx(16,scr2(1))
!
!     rotate about z
!
      scr1(1) = romega1(1)
      scr1(2) = romega1a(1)
      scr1(3) = romega1b(1)
      scr1(4) = romega1c(1)
!
      scr1(5) = iomega1(1)
      scr1(6) = iomega1a(1)
      scr1(7) = iomega1b(1)
      scr1(8) = iomega1c(1)
!
      scr1(9) = romega2(1)
      scr1(10) = romega2a(1)
      scr1(11) = romega2b(1)
      scr1(12) = romega2c(1)
!
      scr1(13) = iomega2(1)
      scr1(14) = iomega2a(1)
      scr1(15) = iomega2b(1)
      scr1(16) = iomega2c(1)
!
      i = 1
      j = 1
!
      do 1 l = 1,nmultipoles
         i = i+1
         j = j+16
         k = l+l
!
         scr1(j) = romega1(i)
         scr1(j+1) = romega1a(i)
         scr1(j+2) = romega1b(i)
         scr1(j+3) = romega1c(i)
!
         scr1(j+4) = iomega1(i)
         scr1(j+5) = iomega1a(i)
         scr1(j+6) = iomega1b(i)
         scr1(j+7) = iomega1c(i)
!
         scr1(j+8) = romega2(i)
         scr1(j+9) = romega2a(i)
         scr1(j+10) = romega2b(i)
         scr1(j+11) = romega2c(i)
!
         scr1(j+12) = iomega2(i)
         scr1(j+13) = iomega2a(i)
         scr1(j+14) = iomega2b(i)
         scr1(j+15) = iomega2c(i)
!
         do 2 m = 1,k,2
            i = i+1
            j = j+16
            scr1(j)=csmphi(m)*romega1(i)-csmphi(m+1)*iomega1(i)
            scr1(j+1)=csmphi(m)*romega1a(i)-csmphi(m+1)*iomega1a(i)
            scr1(j+2)=csmphi(m)*romega1b(i)-csmphi(m+1)*iomega1b(i)
            scr1(j+3)=csmphi(m)*romega1c(i)-csmphi(m+1)*iomega1c(i)
            scr1(j+4)=csmphi(m)*iomega1(i)+csmphi(m+1)*romega1(i)
            scr1(j+5)=csmphi(m)*iomega1a(i)+csmphi(m+1)*romega1a(i)
            scr1(j+6)=csmphi(m)*iomega1b(i)+csmphi(m+1)*romega1b(i)
            scr1(j+7)=csmphi(m)*iomega1c(i)+csmphi(m+1)*romega1c(i)
            scr1(j+8)=csmphipi(m)*romega2(i)-csmphipi(m+1)*iomega2(i)
            scr1(j+9)=csmphipi(m)*romega2a(i)-csmphipi(m+1)*iomega2a(i)
            scr1(j+10)=csmphipi(m)*romega2b(i)-csmphipi(m+1)*iomega2b(i)
            scr1(j+11)=csmphipi(m)*romega2c(i)-csmphipi(m+1)*iomega2c(i)
            scr1(j+12)=csmphipi(m)*iomega2(i)+csmphipi(m+1)*romega2(i)
            scr1(j+13)=csmphipi(m)*iomega2a(i)+csmphipi(m+1)*romega2a(i)
            scr1(j+14)=csmphipi(m)*iomega2b(i)+csmphipi(m+1)*romega2b(i)
            scr1(j+15)=csmphipi(m)*iomega2c(i)+csmphipi(m+1)*romega2c(i)
 2       continue
 1    continue
!
!     rotate about y
!
      mmmm = nmultipoles+1
!
      mmm = 0
!
      scr2(1) = scr1(1)
      scr2(2) = scr1(2)
      scr2(3) = scr1(3)
      scr2(4) = scr1(4)
!
      scr2(5) = scr1(5)
      scr2(6) = scr1(6)
      scr2(7) = scr1(7)
      scr2(8) = scr1(8)
!
      scr2(9) = scr1(9)
      scr2(10) = scr1(10)
      scr2(11) = scr1(11)
      scr2(12) = scr1(12)
!
      scr2(13) = scr1(13)
      scr2(14) = scr1(14)
      scr2(15) = scr1(15)
      scr2(16) = scr1(16)
!
      i = 1
      j = 1
      k = 1
!
      gl = one
!
      do 3 l = 1,nmultipoles
         gl = -gl
         i = i+1
         j = j+16
         k = k+16*l
         n = k
!
         mmm = mmm+1
         a = d2(mmm)*scr1(n)
         aa = d2(mmm)*scr1(n+1)
         aaa = d2(mmm)*scr1(n+2)
         aaaa = d2(mmm)*scr1(n+3)
!
         c = gl*d2(mmm)*scr1(n+8)
         cc = gl*d2(mmm)*scr1(n+9)
         ccc = gl*d2(mmm)*scr1(n+10)
         cccc = gl*d2(mmm)*scr1(n+11)
!
#ifdef FMM_ALLOCALIGNED
         mmm = mmm+1
#endif
!
         nn = n+16
         n = n+16*l
!
         g = gl
!
         do 4 mm = nn,n,16
            g = -g
            mmm = mmm+1
            a = a+d2(mmm)*scr1(mm)
            aa = aa+d2(mmm)*scr1(mm+1)
            aaa = aaa+d2(mmm)*scr1(mm+2)
            aaaa = aaaa+d2(mmm)*scr1(mm+3)
            c = c+g*d2(mmm)*scr1(mm+8)
            cc = cc+g*d2(mmm)*scr1(mm+9)
            ccc = ccc+g*d2(mmm)*scr1(mm+10)
            cccc = cccc+g*d2(mmm)*scr1(mm+11)
#ifdef FMM_ALLOCALIGNED
            mmm = mmm+1
#endif
 4       continue
!
         scr2(j) = a
         scr2(j+1) = aa
         scr2(j+2) = aaa
         scr2(j+3) = aaaa
!
         scr2(j+4) = zero
         scr2(j+5) = zero
         scr2(j+6) = zero
         scr2(j+7) = zero
!
         scr2(j+8) = c
         scr2(j+9) = cc
         scr2(j+10) = ccc
         scr2(j+11) = cccc
!
         scr2(j+12) = zero
         scr2(j+13) = zero
         scr2(j+14) = zero
         scr2(j+15) = zero
!
         mmmmm = mmmm
         mmmmmm = i
!
         glm = gl
!
         do 5 m = 1,l
            glm = -glm
            n = k
!
            a = d2(mmm+1)*scr1(n)
            aa = d2(mmm+1)*scr1(n+1)
            aaa = d2(mmm+1)*scr1(n+2)
            aaaa = d2(mmm+1)*scr1(n+3)
!
            b = d2(mmm+2)*scr1(n+4)
            bb = d2(mmm+2)*scr1(n+5)
            bbb = d2(mmm+2)*scr1(n+6)
            bbbb = d2(mmm+2)*scr1(n+7)
!
            c = glm*d2(mmm+1)*scr1(n+8)
            cc = glm*d2(mmm+1)*scr1(n+9)
            ccc = glm*d2(mmm+1)*scr1(n+10)
            cccc = glm*d2(mmm+1)*scr1(n+11)
!
            d = -glm*d2(mmm+2)*scr1(n+12)
            dd = -glm*d2(mmm+2)*scr1(n+13)
            ddd = -glm*d2(mmm+2)*scr1(n+14)
            dddd = -glm*d2(mmm+2)*scr1(n+15)
!
            mmm = mmm+2
            nn = n+16
            n = n+16*l
!
            g = glm
!
            do 6 nnn = nn,n,16
               g = -g
               a = a+d2(mmm+1)*scr1(nnn)
               aa = aa+d2(mmm+1)*scr1(nnn+1)
               aaa = aaa+d2(mmm+1)*scr1(nnn+2)
               aaaa = aaaa+d2(mmm+1)*scr1(nnn+3)
               b = b+d2(mmm+2)*scr1(nnn+4)
               bb = bb+d2(mmm+2)*scr1(nnn+5)
               bbb = bbb+d2(mmm+2)*scr1(nnn+6)
               bbbb = bbbb+d2(mmm+2)*scr1(nnn+7)
               c = c+g*d2(mmm+1)*scr1(nnn+8)
               cc = cc+g*d2(mmm+1)*scr1(nnn+9)
               ccc = ccc+g*d2(mmm+1)*scr1(nnn+10)
               cccc = cccc+g*d2(mmm+1)*scr1(nnn+11)
               d = d-g*d2(mmm+2)*scr1(nnn+12)
               dd = dd-g*d2(mmm+2)*scr1(nnn+13)
               ddd = ddd-g*d2(mmm+2)*scr1(nnn+14)
               dddd = dddd-g*d2(mmm+2)*scr1(nnn+15)
               mmm = mmm+2
 6          continue
!
            mmmmmm = mmmmmm+mmmmm
            mmmmm = mmmmm-1
            mm = mmmmmm-m
!
            nn = 16*mm-15
!
            scr2(nn) = a
            scr2(nn+1) = aa
            scr2(nn+2) = aaa
            scr2(nn+3) = aaaa
!
            scr2(nn+4) = b
            scr2(nn+5) = bb
            scr2(nn+6) = bbb
            scr2(nn+7) = bbbb
!
            scr2(nn+8) = c
            scr2(nn+9) = cc
            scr2(nn+10) = ccc
            scr2(nn+11) = cccc
!
            scr2(nn+12) = d
            scr2(nn+13) = dd
            scr2(nn+14) = ddd
            scr2(nn+15) = dddd
 5       continue
 3    continue
!
!     perform shift
!
      i = -15
!
      a = zero
      aa = zero
      aaa = zero
      aaaa = zero
!
      c = zero
      cc = zero
      ccc = zero
      cccc = zero
!
      do 7 j = 0,nmultipoles
        i = i+16
        a = a+fr(j)*scr2(i)
        aa = aa+fr(j)*scr2(i+1)
        aaa = aaa+fr(j)*scr2(i+2)
        aaaa = aaaa+fr(j)*scr2(i+3)
        c = c+fr(j)*scr2(i+8)
        cc = cc+fr(j)*scr2(i+9)
        ccc = ccc+fr(j)*scr2(i+10)
        cccc = cccc+fr(j)*scr2(i+11)
 7    continue
!
      scr1(1) = c
      scr1(2) = cc
      scr1(3) = ccc
      scr1(4) = cccc
!
      scr1(5) = zero
      scr1(6) = zero
      scr1(7) = zero
      scr1(8) = zero
!
      scr1(9) = a
      scr1(10) = aa
      scr1(11) = aaa
      scr1(12) = aaaa
!
      scr1(13) = zero
      scr1(14) = zero
      scr1(15) = zero
      scr1(16) = zero
!
      i = 1
!
      do 8 l = 1,nmultipoles
        i = i+16*l
        j = -15
        k = nmultipoles+l
!
        a = zero
        aa = zero
        aaa = zero
        aaaa = zero
!
        c = zero
        cc = zero
        ccc = zero
        cccc = zero
!
        do 9 m = l,k
          j = j+16
          a = a+fr(m)*scr2(j)
          aa = aa+fr(m)*scr2(j+1)
          aaa = aaa+fr(m)*scr2(j+2)
          aaaa = aaaa+fr(m)*scr2(j+3)
          c = c+fr(m)*scr2(j+8)
          cc = cc+fr(m)*scr2(j+9)
          ccc = ccc+fr(m)*scr2(j+10)
          cccc = cccc+fr(m)*scr2(j+11)
 9      continue
!
        g = sg(l)
!
        scr1(i) = g*c
        scr1(i+1) = g*cc
        scr1(i+2) = g*ccc
        scr1(i+3) = g*cccc
!
        scr1(i+4) = zero
        scr1(i+5) = zero
        scr1(i+6) = zero
        scr1(i+7) = zero
!
        scr1(i+8) = g*a
        scr1(i+9) = g*aa
        scr1(i+10) = g*aaa
        scr1(i+11) = g*aaaa
!
        scr1(i+12) = zero
        scr1(i+13) = zero
        scr1(i+14) = zero
        scr1(i+15) = zero
 8    continue
!
      mm = 16*nmultipoles
!
      i = 1
      n = mm+1
!
      do 10 m = 1,nmultipoles
         i = i+16*m
         j = i
!
         do 11 l = m,nmultipoles
            j = j+16*l
            nn = n
            k = m+l
            mmm = nmultipoles+l
!
            a = zero
            aa = zero
            aaa = zero
            aaaa = zero
!
            b = zero
            bb = zero
            bbb = zero
            bbbb = zero
!
            c = zero
            cc = zero
            ccc = zero
            cccc = zero
!
            d = zero
            dd = zero
            ddd = zero
            dddd = zero
!
            do 12 mmmm = k,mmm
               nn = nn+16
               a = a+fr(mmmm)*scr2(nn)
               aa = aa+fr(mmmm)*scr2(nn+1)
               aaa = aaa+fr(mmmm)*scr2(nn+2)
               aaaa = aaaa+fr(mmmm)*scr2(nn+3)
               b = b-fr(mmmm)*scr2(nn+4)
               bb = bb-fr(mmmm)*scr2(nn+5)
               bbb = bbb-fr(mmmm)*scr2(nn+6)
               bbbb = bbbb-fr(mmmm)*scr2(nn+7)
               c = c+fr(mmmm)*scr2(nn+8)
               cc = cc+fr(mmmm)*scr2(nn+9)
               ccc = ccc+fr(mmmm)*scr2(nn+10)
               cccc = cccc+fr(mmmm)*scr2(nn+11)
               d = d-fr(mmmm)*scr2(nn+12)
               dd = dd-fr(mmmm)*scr2(nn+13)
               ddd = ddd-fr(mmmm)*scr2(nn+14)
               dddd = dddd-fr(mmmm)*scr2(nn+15)
 12         continue
!
            g = sg(k)
!
            scr1(j) = g*c
            scr1(j+1) = g*cc
            scr1(j+2) = g*ccc
            scr1(j+3) = g*cccc
!
            scr1(j+4) = g*d
            scr1(j+5) = g*dd
            scr1(j+6) = g*ddd
            scr1(j+7) = g*dddd
!
            scr1(j+8) = g*a
            scr1(j+9) = g*aa
            scr1(j+10) = g*aaa
            scr1(j+11) = g*aaaa
!
            scr1(j+12) = g*b
            scr1(j+13) = g*bb
            scr1(j+14) = g*bbb
            scr1(j+15) = g*bbbb
 11      continue
         n = n+mm
         mm = mm-16
 10   continue
!
!     rotate back expansion
!
      mmm = 0
!
      scr2(jaddress(1)+1) = scr1(1)
      scr2(jaddress(2)+1) = scr1(2)
      scr2(jaddress(3)+1) = scr1(3)
      scr2(jaddress(4)+1) = scr1(4)
!
      scr2(jaddress(5)+1) = scr1(5)
      scr2(jaddress(6)+1) = scr1(6)
      scr2(jaddress(7)+1) = scr1(7)
      scr2(jaddress(8)+1) = scr1(8)
!
      scr2(jaddress(9)+1) = scr1(9)
      scr2(jaddress(10)+1) = scr1(10)
      scr2(jaddress(11)+1) = scr1(11)
      scr2(jaddress(12)+1) = scr1(12)
!
      scr2(jaddress(13)+1) = scr1(13)
      scr2(jaddress(14)+1) = scr1(14)
      scr2(jaddress(15)+1) = scr1(15)
      scr2(jaddress(16)+1) = scr1(16)
!
      i = 1
      j = 1
!
      gl = one
!
      do 13 l = 1,nmultipoles
         gl = -gl
         i = i+1
         j = j+16*l
         n = j
         mmmm = l+l
!
         mmm = mmm+1
         a = d3f(mmm)*scr1(n)
         aa = d3f(mmm)*scr1(n+1)
         aaa = d3f(mmm)*scr1(n+2)
         aaaa = d3f(mmm)*scr1(n+3)
!
         c = gl*d3f(mmm)*scr1(n+8)
         cc = gl*d3f(mmm)*scr1(n+9)
         ccc = gl*d3f(mmm)*scr1(n+10)
         cccc = gl*d3f(mmm)*scr1(n+11)
!
#ifdef FMM_ALLOCALIGNED
         mmm = mmm+1
#endif
!
         nn = n+16
         n = n+16*l
!
         g = gl
!
         do 14 k = nn,n,16
            g = -g
            mmm = mmm+1
            a = a+d3f(mmm)*scr1(k)
            aa = aa+d3f(mmm)*scr1(k+1)
            aaa = aaa+d3f(mmm)*scr1(k+2)
            aaaa = aaaa+d3f(mmm)*scr1(k+3)
            c = c+g*d3f(mmm)*scr1(k+8)
            cc = cc+g*d3f(mmm)*scr1(k+9)
            ccc = ccc+g*d3f(mmm)*scr1(k+10)
            cccc = cccc+g*d3f(mmm)*scr1(k+11)
#ifdef FMM_ALLOCALIGNED
            mmm = mmm+1
#endif
 14      continue
!
         scr2(jaddress(1)+i) = a
         scr2(jaddress(2)+i) = aa
         scr2(jaddress(3)+i) = aaa
         scr2(jaddress(4)+i) = aaaa
!
         scr2(jaddress(5)+i) = zero
         scr2(jaddress(6)+i) = zero
         scr2(jaddress(7)+i) = zero
         scr2(jaddress(8)+i) = zero
!
         scr2(jaddress(9)+i) = c
         scr2(jaddress(10)+i) = cc
         scr2(jaddress(11)+i) = ccc
         scr2(jaddress(12)+i) = cccc
!
         scr2(jaddress(13)+i) = zero
         scr2(jaddress(14)+i) = zero
         scr2(jaddress(15)+i) = zero
         scr2(jaddress(16)+i) = zero
!
         glm = gl
!
         do 15 m = 1,mmmm,2
            glm = -glm
            i = i+1
            n = j
!
            a = d3f(mmm+1)*scr1(n)
            aa = d3f(mmm+1)*scr1(n+1)
            aaa = d3f(mmm+1)*scr1(n+2)
            aaaa = d3f(mmm+1)*scr1(n+3)
!
            b = d3f(mmm+2)*scr1(n+4)
            bb = d3f(mmm+2)*scr1(n+5)
            bbb = d3f(mmm+2)*scr1(n+6)
            bbbb = d3f(mmm+2)*scr1(n+7)
!
            c = glm*d3f(mmm+1)*scr1(n+8)
            cc = glm*d3f(mmm+1)*scr1(n+9)
            ccc = glm*d3f(mmm+1)*scr1(n+10)
            cccc = glm*d3f(mmm+1)*scr1(n+11)
!
            d = -glm*d3f(mmm+2)*scr1(n+12)
            dd = -glm*d3f(mmm+2)*scr1(n+13)
            ddd = -glm*d3f(mmm+2)*scr1(n+14)
            dddd = -glm*d3f(mmm+2)*scr1(n+15)
!
            mmm = mmm+2
            nn = n+16
            n = n+16*l
!
            g = glm
!
            do 16 k = nn,n,16
               g = -g
               a = a+d3f(mmm+1)*scr1(k)
               aa = aa+d3f(mmm+1)*scr1(k+1)
               aaa = aaa+d3f(mmm+1)*scr1(k+2)
               aaaa = aaaa+d3f(mmm+1)*scr1(k+3)
               b = b+d3f(mmm+2)*scr1(k+4)
               bb = bb+d3f(mmm+2)*scr1(k+5)
               bbb = bbb+d3f(mmm+2)*scr1(k+6)
               bbbb = bbbb+d3f(mmm+2)*scr1(k+7)
               c = c+g*d3f(mmm+1)*scr1(k+8)
               cc = cc+g*d3f(mmm+1)*scr1(k+9)
               ccc = ccc+g*d3f(mmm+1)*scr1(k+10)
               cccc = cccc+g*d3f(mmm+1)*scr1(k+11)
               d = d-g*d3f(mmm+2)*scr1(k+12)
               dd = dd-g*d3f(mmm+2)*scr1(k+13)
               ddd = ddd-g*d3f(mmm+2)*scr1(k+14)
               dddd = dddd-g*d3f(mmm+2)*scr1(k+15)
               mmm = mmm+2
 16         continue
!
            scr2(jaddress(1)+i) = csmphipi(m)*a-csmphipi(m+1)*b
            scr2(jaddress(2)+i) = csmphipi(m)*aa-csmphipi(m+1)*bb
            scr2(jaddress(3)+i) = csmphipi(m)*aaa-csmphipi(m+1)*bbb
            scr2(jaddress(4)+i) = csmphipi(m)*aaaa-csmphipi(m+1)*bbbb
            scr2(jaddress(5)+i) = csmphipi(m)*b+csmphipi(m+1)*a
            scr2(jaddress(6)+i) = csmphipi(m)*bb+csmphipi(m+1)*aa
            scr2(jaddress(7)+i) = csmphipi(m)*bbb+csmphipi(m+1)*aaa
            scr2(jaddress(8)+i) = csmphipi(m)*bbbb+csmphipi(m+1)*aaaa
            scr2(jaddress(9)+i) = csmphi(m)*c-csmphi(m+1)*d
            scr2(jaddress(10)+i) = csmphi(m)*cc-csmphi(m+1)*dd
            scr2(jaddress(11)+i) = csmphi(m)*ccc-csmphi(m+1)*ddd
            scr2(jaddress(12)+i) = csmphi(m)*cccc-csmphi(m+1)*dddd
            scr2(jaddress(13)+i) = csmphi(m)*d+csmphi(m+1)*c
            scr2(jaddress(14)+i) = csmphi(m)*dd+csmphi(m+1)*cc
            scr2(jaddress(15)+i) = csmphi(m)*ddd+csmphi(m+1)*ccc
            scr2(jaddress(16)+i) = csmphi(m)*dddd+csmphi(m+1)*cccc
 15      continue
 13   continue
!
      i = 0
!
      if(jacc(1)) then
         do 17 j = 1,nsqmultipoles
            i = i+1
            mu1(j) = mu1(j)+scr2(i)
 17      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 18 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 18      continue
      endif
!
      if(jacc(2)) then
         do 19 j = 1,nsqmultipoles
            i = i+1
            mu2(j) = mu2(j)+scr2(i)
 19      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 20 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 20      continue
      endif
!
      if(jacc(3)) then
         do 21 j = 1,nsqmultipoles
            i = i+1
            mu3(j) = mu3(j)+scr2(i)
 21      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 22 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 22      continue
      endif
!
      if(jacc(4)) then
         do 23 j = 1,nsqmultipoles
            i = i+1
            mu4(j) = mu4(j)+scr2(i)
 23      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 24 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 24      continue
      endif
!
      if(jacc(5)) then
         do 25 j = 1,nsqmultipoles
            i = i+1
            mu5(j) = mu5(j)+scr2(i)
 25      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 26 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 26      continue
      endif
!
      if(jacc(6)) then
         do 27 j = 1,nsqmultipoles
            i = i+1
            mu6(j) = mu6(j)+scr2(i)
 27      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 28 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 28      continue
      endif
!
      if(jacc(7)) then
         do 29 j = 1,nsqmultipoles
            i = i+1
            mu7(j) = mu7(j)+scr2(i)
 29      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 30 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 30      continue
      endif
!
      if(jacc(8)) then
         do 31 j = 1,nsqmultipoles
            i = i+1
            mu8(j) = mu8(j)+scr2(i)
 31      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 32 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 32      continue
      endif
!
      if(jacc(9)) then
         do 33 j = 1,nsqmultipoles
            i = i+1
            mu9(j) = mu9(j)+scr2(i)
 33      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 34 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 34      continue
      endif
!
      if(jacc(10)) then
         do 35 j = 1,nsqmultipoles
            i = i+1
            mu10(j) = mu10(j)+scr2(i)
 35      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 36 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 36      continue
      endif
!
      if(jacc(11)) then
         do 37 j = 1,nsqmultipoles
            i = i+1
            mu11(j) = mu11(j)+scr2(i)
 37      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 38 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 38      continue
      endif
!
      if(jacc(12)) then
         do 39 j = 1,nsqmultipoles
            i = i+1
            mu12(j) = mu12(j)+scr2(i)
 39      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 40 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 40      continue
      endif
!
      if(jacc(13)) then
         do 41 j = 1,nsqmultipoles
            i = i+1
            mu13(j) = mu13(j)+scr2(i)
 41      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 42 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 42      continue
      endif
!
      if(jacc(14)) then
         do 43 j = 1,nsqmultipoles
            i = i+1
            mu14(j) = mu14(j)+scr2(i)
 43      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 44 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 44      continue
      endif
!
      if(jacc(15)) then
         do 45 j = 1,nsqmultipoles
            i = i+1
            mu15(j) = mu15(j)+scr2(i)
 45      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 46 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 46      continue
      endif
!
      if(jacc(16)) then
         do 47 j = 1,nsqmultipoles
            i = i+1
            mu16(j) = mu16(j)+scr2(i)
 47      continue
      else
         j = i+nsqmultipoles
         k = j+nsqmultipoles
         j = j+1
         do 48 l = j,k
            i = i+1
            scr2(l) = scr2(l)+scr2(i)
 48      continue
      endif
!      call f_hpm_stop_i8_(3,__LINE__,__FILE__)
      return
      end subroutine pass2trfrqdcach
