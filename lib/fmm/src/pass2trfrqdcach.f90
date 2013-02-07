!
! mpixlf2003_r -WF,-DHAVE_FCONFIG_H -I. -I..   -I. -qintsize=8  -qfixed -O3 -qhot=nosimd -qmaxmem=-1 -qarch=450d -qtune=450 -c -o pass2trfrqdcach.o -c
! -qsuffix=cpp=f -qlist -qsource -qreport -qattr pass2trfrqdcach.f
! -qhot=nosimd bringt 3-4% gegenüber ohne -qhot
! -qhot, -O4, -O5 immer langsamer
#include "fmm.h"
!
      subroutine pass2trfrqdcach(nmultipoles,nsqmultipoles,jaddress,&
      jacc,romega1,iomega1,romega2,iomega2,romega1a,iomega1a,romega2a,&
      iomega2a,romega1b,iomega1b,romega2b,iomega2b,romega1c,iomega1c,&
      romega2c,iomega2c,mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11,&
      mu12,mu13,mu14,mu15,mu16,cmphi,smphi,cmphipi,smphipi,csmphi,&
      csmphipi,sg,fr,d2,d3f,scr1,scr2)
!
!     use fmmkinds
!
      implicit none
!
      real(kind=8) romega1(*),iomega1(*),romega2(*),iomega2(*),&
      romega1a(*),iomega1a(*),romega2a(*),iomega2a(*),romega1b(*),&
      iomega1b(*),romega2b(*),iomega2b(*),romega1c(*),iomega1c(*),&
      romega2c(*),iomega2c(*),mu1(*),mu2(*),mu3(*),mu4(*),mu5(*),&
      mu6(*),mu7(*),mu8(*),mu9(*),mu10(*),mu11(*),mu12(*),mu13(*),&
      mu14(*),mu15(*),mu16(*),cmphi(0:*),smphi(0:*),cmphipi(0:*),&
      smphipi(0:*),csmphi(*),csmphipi(*),sg(0:*),fr(0:*),d2(*),d3f(*),&
      scr1(*),scr2(*),gl,g,glm
!
      complex(kind=8) :: reg00,reg01,reg02,reg03,reg04,reg05,reg06,reg07
      complex(kind=8) :: reg08,reg09,reg10,reg11,reg12,reg13,reg14,reg15
      complex(kind=8) :: reg16,reg17     
      
      integer(kind=8) nmultipoles,nsqmultipoles,jaddress(*)
      integer(kind=4) i,j,l,m,mmmm,mmm,k,n,nn,mm,mmmmm,mmmmmm,nnn
!
      logical(kind=8) jacc(*)
!
      real(kind=8) zero
      parameter(zero=0.e0_8)
      real(kind=8) one
      parameter(one=1.e0_8)
      
      real(kind=8), dimension(16) :: a
 
!
!      call f_hpm_start_i8_(3,__LINE__,__FILE__,'rot')
!
      call alignx(16,d2(1))
      call alignx(16,d3f(1))
      call alignx(16,scr1(1))
      call alignx(16,scr2(1))
!
!     rotate about z
!
      scr1( 1) = romega1(1)
      scr1( 2) = romega1a(1)
      scr1( 3) = romega1b(1)
      scr1( 4) = romega1c(1)
!
      scr1( 5) = iomega1(1)
      scr1( 6) = iomega1a(1)
      scr1( 7) = iomega1b(1)
      scr1( 8) = iomega1c(1)
!
      scr1( 9) = romega2(1)
      scr1(10) = romega2a(1)
      scr1(11) = romega2b(1)
      scr1(12) = romega2c(1)
!
      scr1(13) = iomega2(1)
      scr1(14) = iomega2a(1)
      scr1(15) = iomega2b(1)
      scr1(16) = iomega2c(1)

      i = 1
      j = 1
!
      do 1 l = 1,nmultipoles
         i = i+1
         j = j+16
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
!         do 2 m = 1,l
!-ik
         do 2 m = 1,2*l,2
!-ik
            i = i+1
            j = j+16
!-ik
#ifndef QWERTZ
            call alignx(16,csmphi(m))
            call alignx(16,csmphipi(m))

            reg16 = loadfp(csmphi(m))
            reg17 = loadfp(csmphipi(m))

!            reg16 = cmplx(cmphi(m),smphi(m),8)
!            reg17 = cmplx(cmphipi(m),smphipi(m),8)

            reg00 = cmplx(romega1( i),romega1a(i),8)
            reg01 = cmplx(romega1b(i),romega1c(i),8)
            reg02 = cmplx(iomega1( i),iomega1a(i),8)
            reg03 = cmplx(iomega1b(i),iomega1c(i),8)

            reg04 = cmplx(romega2( i),romega2a(i),8)
            reg05 = cmplx(romega2b(i),romega2c(i),8)
            reg06 = cmplx(iomega2( i),iomega2a(i),8)
            reg07 = cmplx(iomega2b(i),iomega2c(i),8)

            reg08 = fxpmul(reg00,real(reg16))
            reg08 = fxcpnmsub(reg08,reg02,imag(reg16))
            reg09 = fxpmul(reg01,real(reg16))
            reg09 = fxcpnmsub(reg09,reg03,imag(reg16))

            reg10 = fxpmul(reg02,real(reg16))
            reg10 = fxcpmadd(reg10,reg00,imag(reg16))
            reg11 = fxpmul(reg03,real(reg16))
            reg11 = fxcpmadd(reg11,reg01,imag(reg16))

            reg12 = fxpmul(reg04,real(reg17))
            reg12 = fxcpnmsub(reg12,reg06,imag(reg17))
            reg13 = fxpmul(reg05,real(reg17))
            reg13 = fxcpnmsub(reg13,reg07,imag(reg17))

            reg14 = fxpmul(reg06,real(reg17))
            reg14 = fxcpmadd(reg14,reg04,imag(reg17))
            reg15 = fxpmul(reg07,real(reg17))
            reg15 = fxcpmadd(reg15,reg05,imag(reg17))

        
!            reg08 = fxsmul(reg02,-imag(reg16))
!            reg08 = fxcpmadd(reg08,reg00,real(reg16))
!            reg09 = fxsmul(reg03,-imag(reg16))
!            reg09 = fxcpmadd(reg09,reg01,real(reg16))

!            reg10 = fxsmul(reg00,imag(reg16))
!            reg10 = fxcpmadd(reg10,reg02,real(reg16))
!            reg11 = fxsmul(reg01,imag(reg16))
!            reg11 = fxcpmadd(reg11,reg03,real(reg16))

!            reg12 = fxsmul(reg06,-imag(reg17))
!            reg12 = fxcpmadd(reg12,reg04,real(reg17))
!            reg13 = fxsmul(reg07,-imag(reg17))
!            reg13 = fxcpmadd(reg13,reg05,real(reg17))

!            reg14 = fxsmul(reg04,imag(reg17))
!            reg14 = fxcpmadd(reg14,reg06,real(reg17))
!            reg15 = fxsmul(reg05,imag(reg17))
!            reg15 = fxcpmadd(reg15,reg07,real(reg17))
	    
!            reg08 = fxpmul(reg00,real(reg16))
!            reg08 = fxcpmadd(reg08,reg02,-imag(reg16))
!            reg09 = fxpmul(reg01,real(reg16))
!            reg09 = fxcpmadd(reg09,reg03,-imag(reg16))

!            reg10 = fxpmul(reg02,real(reg16))
!            reg10 = fxcpmadd(reg10,reg00,imag(reg16))
!            reg11 = fxpmul(reg03,real(reg16))
!            reg11 = fxcpmadd(reg11,reg01,imag(reg16))

!            reg12 = fxpmul(reg04,real(reg17))
!            reg12 = fxcpmadd(reg12,reg06,-imag(reg17))
!            reg13 = fxpmul(reg05,real(reg17))
!            reg13 = fxcpmadd(reg13,reg07,-imag(reg17))

!            reg14 = fxpmul(reg06,real(reg17))
!            reg14 = fxcpmadd(reg14,reg04,imag(reg17))
!            reg15 = fxpmul(reg07,real(reg17))
!            reg15 = fxcpmadd(reg15,reg05,imag(reg17))
 
!            scr1(j   ) = real(reg08)
!            scr1(j+ 1) = imag(reg08)	    
!            scr1(j+ 2) = real(reg09)
!            scr1(j+ 3) = imag(reg09)	    
!            scr1(j+ 4) = real(reg10)
!            scr1(j+ 5) = imag(reg10)	    
!            scr1(j+ 6) = real(reg11)
!            scr1(j+ 7) = imag(reg11)	    
!            scr1(j+ 8) = real(reg12)
!            scr1(j+ 9) = imag(reg12)	    
!            scr1(j+10) = real(reg13)
!            scr1(j+11) = imag(reg13)	    
!            scr1(j+12) = real(reg14)
!            scr1(j+13) = imag(reg14)	    
!            scr1(j+14) = real(reg15)
!            scr1(j+15) = imag(reg15)	    

            call storefp(scr1(j   ),reg08)
            call storefp(scr1(j+ 2),reg09)
            call storefp(scr1(j+ 4),reg10)
            call storefp(scr1(j+ 6),reg11)
            call storefp(scr1(j+ 8),reg12)
            call storefp(scr1(j+10),reg13)
            call storefp(scr1(j+12),reg14)
            call storefp(scr1(j+14),reg15)
!-ik
#else

            scr1(j   ) = cmphi(m)*romega1( i)-smphi(m)*iomega1(i)
            scr1(j+ 1) = cmphi(m)*romega1a(i)-smphi(m)*iomega1a(i)
            scr1(j+ 2) = cmphi(m)*romega1b(i)-smphi(m)*iomega1b(i)
            scr1(j+ 3) = cmphi(m)*romega1c(i)-smphi(m)*iomega1c(i)
            scr1(j+ 4) = cmphi(m)*iomega1( i)+smphi(m)*romega1(i)
            scr1(j+ 5) = cmphi(m)*iomega1a(i)+smphi(m)*romega1a(i)
            scr1(j+ 6) = cmphi(m)*iomega1b(i)+smphi(m)*romega1b(i)
            scr1(j+ 7) = cmphi(m)*iomega1c(i)+smphi(m)*romega1c(i)

            scr1(j+ 8) = cmphipi(m)*romega2( i)-smphipi(m)*iomega2(i)
            scr1(j+ 9) = cmphipi(m)*romega2a(i)-smphipi(m)*iomega2a(i)
            scr1(j+10) = cmphipi(m)*romega2b(i)-smphipi(m)*iomega2b(i)
            scr1(j+11) = cmphipi(m)*romega2c(i)-smphipi(m)*iomega2c(i)
            scr1(j+12) = cmphipi(m)*iomega2( i)+smphipi(m)*romega2(i)
            scr1(j+13) = cmphipi(m)*iomega2a(i)+smphipi(m)*romega2a(i)
            scr1(j+14) = cmphipi(m)*iomega2b(i)+smphipi(m)*romega2b(i)
            scr1(j+15) = cmphipi(m)*iomega2c(i)+smphipi(m)*romega2c(i)
#endif
 2       continue
 1    continue
!
!     rotate about y
!
      call rotatey(nmultipoles,scr1,scr2,d2)
!
      call performshift(nmultipoles,scr1,scr2,fr,sg)
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
!
         mmm = mmm+1
	 
         call alignx(16,scr1(n))

         reg00 = loadfp(scr1(n))
         reg01 = loadfp(scr1(n+2))
         reg04 = loadfp(scr1(n+8))
         reg05 = loadfp(scr1(n+10))
	 
         reg08 =  fxpmul(reg00,d3f(mmm))
         reg09 =  fxpmul(reg01,d3f(mmm))	 

         reg04 =  fxpmul(reg04,gl)
         reg12 =  fxpmul(reg04,d3f(mmm))	 
         reg05 =  fxpmul(reg05,gl)
         reg13 =  fxpmul(reg05,d3f(mmm))	 

!         a(1) = d3f(mmm)*scr1(n)
!         a(2) = d3f(mmm)*scr1(n+1)
!         a(3) = d3f(mmm)*scr1(n+2)
!         a(4) = d3f(mmm)*scr1(n+3)
!
!         a(9) = gl*d3f(mmm)*scr1(n+8)
!         a(10) = gl*d3f(mmm)*scr1(n+9)
!         a(11) = gl*d3f(mmm)*scr1(n+10)
!         a(12) = gl*d3f(mmm)*scr1(n+11)
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
            call alignx(16,scr1(k))

            reg00 = loadfp(scr1(k))
            reg01 = loadfp(scr1(k+2))
            reg04 = loadfp(scr1(k+8))
            reg05 = loadfp(scr1(k+10))
	    
            reg08 = fxcpmadd(reg08,reg00,d3f(mmm))
            reg09 = fxcpmadd(reg09,reg01,d3f(mmm))
            
            reg04 = fxpmul(reg04,g)
            reg12 = fxcpmadd(reg12,reg04,d3f(mmm))	    	    	    

            reg05 = fxpmul(reg05,g)
            reg13 = fxcpmadd(reg13,reg05,d3f(mmm))	    	    	    
	    
!            a(1) = a(1)+d3f(mmm)*scr1(k)
!            a(2) = a(2)+d3f(mmm)*scr1(k+1)
!            a(3) = a(3)+d3f(mmm)*scr1(k+2)
!            a(4) = a(4)+d3f(mmm)*scr1(k+3)
!            a(9) = a(9)+g*d3f(mmm)*scr1(k+8)
!            a(10) = a(10)+g*d3f(mmm)*scr1(k+9)
!            a(11) = a(11)+g*d3f(mmm)*scr1(k+10)
!            a(12) = a(12)+g*d3f(mmm)*scr1(k+11)
#ifdef FMM_ALLOCALIGNED
            mmm = mmm+1
#endif
 14      continue
 
#ifndef QWERTZ
!-ik
         scr2(jaddress(1)+i) = real(reg08)
         scr2(jaddress(2)+i) = imag(reg08)
         scr2(jaddress(3)+i) = real(reg09)
         scr2(jaddress(4)+i) = imag(reg09)
!
         scr2(jaddress(5)+i) = zero
         scr2(jaddress(6)+i) = zero
         scr2(jaddress(7)+i) = zero
         scr2(jaddress(8)+i) = zero
!
         scr2(jaddress(9)+i) = real(reg12)
         scr2(jaddress(10)+i) = imag(reg12)
         scr2(jaddress(11)+i) = real(reg13)
         scr2(jaddress(12)+i) = imag(reg13)
!
         scr2(jaddress(13)+i) = zero
         scr2(jaddress(14)+i) = zero
         scr2(jaddress(15)+i) = zero
         scr2(jaddress(16)+i) = zero

!-ik
#else	 	 
!
         call storefp(a(1),reg08)
         call storefp(a(3),reg09)

         call storefp(a(9),reg12)
         call storefp(a(11),reg13) 

         scr2(jaddress(1)+i) = a(1)
         scr2(jaddress(2)+i) = a(2)
         scr2(jaddress(3)+i) = a(3)
         scr2(jaddress(4)+i) = a(4)
!
         scr2(jaddress(5)+i) = zero
         scr2(jaddress(6)+i) = zero
         scr2(jaddress(7)+i) = zero
         scr2(jaddress(8)+i) = zero
!
         scr2(jaddress(9)+i) = a(9)
         scr2(jaddress(10)+i) = a(10)
         scr2(jaddress(11)+i) = a(11)
         scr2(jaddress(12)+i) = a(12)
!
         scr2(jaddress(13)+i) = zero
         scr2(jaddress(14)+i) = zero
         scr2(jaddress(15)+i) = zero
         scr2(jaddress(16)+i) = zero
#endif
!
         glm = gl
!
!         do 15 m = 1,l
         do 15 m = 1,2*l,2
            glm = -glm
            i = i+1
            n = j
!
            call alignx(16,scr1(n))

            reg00 = loadfp(scr1(n))
            reg01 = loadfp(scr1(n+2))
            reg02 = loadfp(scr1(n+4))
            reg03 = loadfp(scr1(n+6))
            reg04 = loadfp(scr1(n+8))
            reg05 = loadfp(scr1(n+10))
            reg06 = loadfp(scr1(n+12))
            reg07 = loadfp(scr1(n+14))
	    
            reg08 = fxpmul(reg00,d3f(mmm+1))	    
            reg09 = fxpmul(reg01,d3f(mmm+1))
            reg10 = fxpmul(reg02,d3f(mmm+2))
            reg11 = fxpmul(reg03,d3f(mmm+2))	    

            reg04 = fxpmul(reg04,glm)
            reg12 = fxpmul(reg04,d3f(mmm+1))	    
            reg05 = fxpmul(reg05,glm)
            reg13 = fxpmul(reg05,d3f(mmm+1))	    
            reg06 = fxpmul(reg06,-glm)
            reg14 = fxpmul(reg06,d3f(mmm+2))	    
            reg07 = fxpmul(reg07,-glm)
            reg15 = fxpmul(reg07,d3f(mmm+2))	    


!            a(1) = d3f(mmm+1)*scr1(n)
!            a(2) = d3f(mmm+1)*scr1(n+1)
!            a(3) = d3f(mmm+1)*scr1(n+2)
!            a(4) = d3f(mmm+1)*scr1(n+3)
!
!            a(5) = d3f(mmm+2)*scr1(n+4)
!            a(6) = d3f(mmm+2)*scr1(n+5)
!            a(7) = d3f(mmm+2)*scr1(n+6)
!            a(8) = d3f(mmm+2)*scr1(n+7)
!
!            a(9) = glm*d3f(mmm+1)*scr1(n+8)
!            a(10) = glm*d3f(mmm+1)*scr1(n+9)
!            a(11) = glm*d3f(mmm+1)*scr1(n+10)
!            a(12) = glm*d3f(mmm+1)*scr1(n+11)
!
!            a(13) = -glm*d3f(mmm+2)*scr1(n+12)
!            a(14) = -glm*d3f(mmm+2)*scr1(n+13)
!            a(15) = -glm*d3f(mmm+2)*scr1(n+14)
!            a(16) = -glm*d3f(mmm+2)*scr1(n+15)
!
            mmm = mmm+2
            nn = n+16
            n = n+16*l
!
            g = glm
!
            do 16 k = nn,n,16
               g = -g
               call alignx(16,scr1(k))
	       
               reg00 = loadfp(scr1(k))
               reg01 = loadfp(scr1(k+2))
               reg02 = loadfp(scr1(k+4))
               reg03 = loadfp(scr1(k+6))
               reg04 = loadfp(scr1(k+8))
               reg05 = loadfp(scr1(k+10))
               reg06 = loadfp(scr1(k+12))
               reg07 = loadfp(scr1(k+14))
	       
               reg08 = fxcpmadd(reg08,reg00,d3f(mmm+1))
               reg09 = fxcpmadd(reg09,reg01,d3f(mmm+1))
               reg10 = fxcpmadd(reg10,reg02,d3f(mmm+2))
               reg11 = fxcpmadd(reg11,reg03,d3f(mmm+2))	       	       	       
               
               reg04 = fxpmul(reg04,g)
               reg12 = fxcpmadd(reg12,reg04,d3f(mmm+1))	       	       	       
               reg05 = fxpmul(reg05,g)
               reg13 = fxcpmadd(reg13,reg05,d3f(mmm+1))	       	       	       
               reg06 = fxpmul(reg06,-g)
               reg14 = fxcpmadd(reg14,reg06,d3f(mmm+2))	       	       	       
               reg07 = fxpmul(reg07,-g)
               reg15 = fxcpmadd(reg15,reg07,d3f(mmm+2))	       	       	       
	       
!               a( 1) = a( 1)+d3f(mmm+1)*scr1(k)
!               a( 2) = a( 2)+d3f(mmm+1)*scr1(k+1)
!               a( 3) = a( 3)+d3f(mmm+1)*scr1(k+2)
!               a( 4) = a( 4)+d3f(mmm+1)*scr1(k+3)
!               a( 5) = a( 5)+d3f(mmm+2)*scr1(k+4)
!               a( 6) = a( 6)+d3f(mmm+2)*scr1(k+5)
!               a( 7) = a( 7)+d3f(mmm+2)*scr1(k+6)
!               a( 8) = a( 8)+d3f(mmm+2)*scr1(k+7)
!               a( 9) = a( 9)+g*d3f(mmm+1)*scr1(k+8)
!               a(10) = a(10)+g*d3f(mmm+1)*scr1(k+9)
!               a(11) = a(11)+g*d3f(mmm+1)*scr1(k+10)
!               a(12) = a(12)+g*d3f(mmm+1)*scr1(k+11)
!               a(13) = a(13)-g*d3f(mmm+2)*scr1(k+12)
!               a(14) = a(14)-g*d3f(mmm+2)*scr1(k+13)
!               a(15) = a(15)-g*d3f(mmm+2)*scr1(k+14)
!               a(16) = a(16)-g*d3f(mmm+2)*scr1(k+15)	
	    
               mmm = mmm+2
 16         continue
!
!-ik
#ifndef QWERTZ
            call alignx(16,csmphi(m))
            call alignx(16,csmphipi(m))

            reg16 = loadfp(csmphi(m))
            reg17 = loadfp(csmphipi(m)) 

!            reg16 = cmplx(cmphi(m),smphi(m),8)
!            reg17 = cmplx(cmphipi(m),smphipi(m),8)

            reg00 = fxpmul(reg08,real(reg17))
            reg00 = fxcpnmsub(reg00,reg10,imag(reg17))
            reg01 = fxpmul(reg09,real(reg17))
            reg01 = fxcpnmsub(reg01,reg11,imag(reg17))
            reg02 = fxpmul(reg10,real(reg17))
            reg02 = fxcpmadd(reg02,reg08,imag(reg17))
            reg03 = fxpmul(reg11,real(reg17))
            reg03 = fxcpmadd(reg03,reg09,imag(reg17))

            reg04 = fxpmul(reg12,real(reg16))
            reg04 = fxcpnmsub(reg04,reg14,imag(reg16))
            reg05 = fxpmul(reg13,real(reg16))
            reg05 = fxcpnmsub(reg05,reg15,imag(reg16))
            reg06 = fxpmul(reg14,real(reg16))
            reg06 = fxcpmadd(reg06,reg12,imag(reg16))
            reg07 = fxpmul(reg15,real(reg16))
            reg07 = fxcpmadd(reg07,reg13,imag(reg16))

            scr2(jaddress( 1)+i) = real(reg00)
            scr2(jaddress( 2)+i) = imag(reg00)
            scr2(jaddress( 3)+i) = real(reg01)
            scr2(jaddress( 4)+i) = imag(reg01)
            scr2(jaddress( 5)+i) = real(reg02)
            scr2(jaddress( 6)+i) = imag(reg02)
            scr2(jaddress( 7)+i) = real(reg03)
            scr2(jaddress( 8)+i) = imag(reg03)

            scr2(jaddress( 9)+i) = real(reg04)
            scr2(jaddress(10)+i) = imag(reg04)
            scr2(jaddress(11)+i) = real(reg05)
            scr2(jaddress(12)+i) = imag(reg05)
            scr2(jaddress(13)+i) = real(reg06)
            scr2(jaddress(14)+i) = imag(reg06)
            scr2(jaddress(15)+i) = real(reg07)
            scr2(jaddress(16)+i) = imag(reg07)
!-ik
#else
            call storefp(a( 1),reg08)
            call storefp(a( 3),reg09)
            call storefp(a( 5),reg10)
            call storefp(a( 7),reg11)
            call storefp(a( 9),reg12)
            call storefp(a(11),reg13)
            call storefp(a(13),reg14)
            call storefp(a(15),reg15)
    	    	    	    	    	    	    
            scr2(jaddress( 1)+i) = cmphipi(m)*a(1)-smphipi(m)*a(5)
            scr2(jaddress( 2)+i) = cmphipi(m)*a(2)-smphipi(m)*a(6)
            scr2(jaddress( 3)+i) = cmphipi(m)*a(3)-smphipi(m)*a(7)
            scr2(jaddress( 4)+i) = cmphipi(m)*a(4)-smphipi(m)*a(8)
            scr2(jaddress( 5)+i) = cmphipi(m)*a(5)+smphipi(m)*a(1)
            scr2(jaddress( 6)+i) = cmphipi(m)*a(6)+smphipi(m)*a(2)
            scr2(jaddress( 7)+i) = cmphipi(m)*a(7)+smphipi(m)*a(3)
            scr2(jaddress( 8)+i) = cmphipi(m)*a(8)+smphipi(m)*a(4)
            scr2(jaddress( 9)+i) = cmphi(m)*a( 9)-smphi(m)*a(13)
            scr2(jaddress(10)+i) = cmphi(m)*a(10)-smphi(m)*a(14)
            scr2(jaddress(11)+i) = cmphi(m)*a(11)-smphi(m)*a(15)
            scr2(jaddress(12)+i) = cmphi(m)*a(12)-smphi(m)*a(16)
            scr2(jaddress(13)+i) = cmphi(m)*a(13)+smphi(m)*a( 9)
            scr2(jaddress(14)+i) = cmphi(m)*a(14)+smphi(m)*a(10)
            scr2(jaddress(15)+i) = cmphi(m)*a(15)+smphi(m)*a(11)
            scr2(jaddress(16)+i) = cmphi(m)*a(16)+smphi(m)*a(12)
#endif
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

      subroutine rotatey(nmultipoles,scr1,scr2,d2)
      implicit none
      
      real(kind=8) :: scr1(*),scr2(*),d2(*)
      integer(kind=8) :: nmultipoles
      real(kind=8), parameter  :: one = 1.0d0
      real(kind=8), parameter  :: zero = 0.0d0
      real(kind=8) :: g,gl,glm
      real(kind=8), dimension(16) :: a
      integer(kind=4) :: mmmmmm,mmmmm,mmmm,mmm,mm,m
      integer(kind=4) :: i,j,k,l,n,nn,nnn
      
      complex(kind=8) :: reg00,reg01,reg02,reg03,reg04,reg05,reg06,reg07
      complex(kind=8) :: reg08,reg09,reg10,reg11,reg12,reg13,reg14,reg15
      complex(kind=8) :: reg16
      mmmm = nmultipoles+1
!
      mmm = 0
!
#ifdef QWERTZ
      call alignx(16,scr1(1))
      call alignx(16,scr2(1))
      
      reg00 = loadfp(scr1(1))
      reg01 = loadfp(scr1(3))
      reg02 = loadfp(scr1(5))
      reg03 = loadfp(scr1(7))
      reg04 = loadfp(scr1(9))
      reg05 = loadfp(scr1(11))
      reg06 = loadfp(scr1(13))
      reg07 = loadfp(scr1(15))
      
      call storefp(scr2(1),reg00)
      call storefp(scr2(3),reg01)
      call storefp(scr2(5),reg02)
      call storefp(scr2(7),reg03)
      call storefp(scr2(9),reg04)
      call storefp(scr2(11),reg05)
      call storefp(scr2(13),reg06)
      call storefp(scr2(15),reg07)	      
#else
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
#endif
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
         call alignx(16,scr1(n))
         call alignx(16,d2(mmm))	 
	 
         reg00 = loadfp(scr1(n))
         reg01 = loadfp(scr1(n+2))
	 
         reg04 = loadfp(scr1(n+8))
         reg05 = loadfp(scr1(n+10))
! gl und d2(*) gleich in D matrizen abspeichern	 
         reg08 = fxpmul(reg00,d2(mmm))
         reg09 = fxpmul(reg01,d2(mmm))
	 
         reg12 = fxpmul(reg04,d2(mmm))
         reg12 = fxpmul(reg12,gl)
         reg13 = fxpmul(reg05,d2(mmm))
         reg13 = fxpmul(reg13,gl) 

!         a(1) = d2(mmm)*scr1(n)
!         a(2) = d2(mmm)*scr1(n+1)
!         a(3) = d2(mmm)*scr1(n+2)
!         a(4) = d2(mmm)*scr1(n+3)
!
!         a(9) = gl*d2(mmm)*scr1(n+8)
!         a(10) = gl*d2(mmm)*scr1(n+9)
!         a(11) = gl*d2(mmm)*scr1(n+10)
!         a(12) = gl*d2(mmm)*scr1(n+11)
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

            call alignx(16,scr1(mm))
            call alignx(16,d2(mmm))
	   
            reg00 = loadfp(scr1(mm))
            reg01 = loadfp(scr1(mm+2))
            reg04 = loadfp(scr1(mm+8))
            reg05 = loadfp(scr1(mm+10))
	   
            reg08 = fxcpmadd(reg08,reg00,d2(mmm))
            reg09 = fxcpmadd(reg09,reg01,d2(mmm))
            reg04 = fxpmul(reg04,g)
            reg12 = fxcpmadd(reg12,reg04,d2(mmm))
            reg05 = fxpmul(reg05,g)
            reg13 = fxcpmadd(reg13,reg05,d2(mmm))
	   
!            a(1) = a(1)+d2(mmm)*scr1(mm)
!            a(2) = a(2)+d2(mmm)*scr1(mm+1)
!            a(3) = a(3)+d2(mmm)*scr1(mm+2)
!            a(4) = a(4)+d2(mmm)*scr1(mm+3)
!            a(9) = a(9)+g*d2(mmm)*scr1(mm+8)
!            a(10) = a(10)+g*d2(mmm)*scr1(mm+9)
!            a(11) = a(11)+g*d2(mmm)*scr1(mm+10)
!            a(12) = a(12)+g*d2(mmm)*scr1(mm+11)


	 	    
#ifdef FMM_ALLOCALIGNED
            mmm = mmm+1
#endif
 4       continue
!    
         call alignx(16,scr2(j))
	 
         call storefp(scr2(j),reg08)
         call storefp(scr2(j+2),reg09)
!         call storefp(scr2(j+4),(zero,zero))
!         call storefp(scr2(j+6),(zero,zero))

         call storefp(scr2(j+8),reg12)
         call storefp(scr2(j+10),reg13)	 
!         call storefp(scr2(j+12),(zero,zero))
!         call storefp(scr2(j+14),(zero,zero))
	 	 
!         scr2(j) = a(1)
!         scr2(j+1) = a(2)
!         scr2(j+2) = a(3)
!         scr2(j+3) = a(4)
!c
         scr2(j+4) = zero
         scr2(j+5) = zero
         scr2(j+6) = zero
         scr2(j+7) = zero
!c
!         scr2(j+8) = a(9)
!         scr2(j+9) = a(10)
!         scr2(j+10) = a(11)
!         scr2(j+11) = a(12)
!c
         scr2(j+12) = zero
         scr2(j+13) = zero
         scr2(j+14) = zero
         scr2(j+15) = zero
!c
         mmmmm = mmmm
         mmmmmm = i
!
         glm = gl
!
         do 5 m = 1,l
            glm = -glm
            n = k
!
            call alignx(16,scr1(n))
            call alignx(16,d2(mmm+1))

            reg00 = loadfp(scr1(n))
            reg01 = loadfp(scr1(n+2))
            reg02 = loadfp(scr1(n+4))
            reg03 = loadfp(scr1(n+6))
            reg04 = loadfp(scr1(n+8))
            reg05 = loadfp(scr1(n+10))
            reg06 = loadfp(scr1(n+12))
            reg07 = loadfp(scr1(n+14))

            reg16 = loadfp(d2(mmm+1))
	    
            reg08 = fxpmul(reg00,real(reg16))
            reg09 = fxpmul(reg01,real(reg16))
	    
            reg10 = fxpmul(reg02,imag(reg16))
            reg11 = fxpmul(reg03,imag(reg16))	    	    	    

            reg04 = fxpmul(reg04,glm)
            reg12 = fxpmul(reg04,real(reg16))
            reg05 = fxpmul(reg05,glm)
            reg13 = fxpmul(reg05,real(reg16))
	    
            reg06 = fxpmul(reg06,-glm)
!           reg06 = fpneg(reg06)
            reg14 = fxpmul(reg06,imag(reg16))
            reg07 = fxpmul(reg07,-glm)
!           reg07 = fpneg(reg07)
            reg15 = fxpmul(reg07,imag(reg16))	    

	    
!            a(1) = d2(mmm+1)*scr1(n)
!            a(2) = d2(mmm+1)*scr1(n+1)
!            a(3) = d2(mmm+1)*scr1(n+2)
!            a(4) = d2(mmm+1)*scr1(n+3)
!
!            a(5) = d2(mmm+2)*scr1(n+4)
!            a(6) = d2(mmm+2)*scr1(n+5)
!            a(7) = d2(mmm+2)*scr1(n+6)
!            a(8) = d2(mmm+2)*scr1(n+7)
!
!            a(9) = glm*d2(mmm+1)*scr1(n+8)
!            a(10) = glm*d2(mmm+1)*scr1(n+9)
!            a(11) = glm*d2(mmm+1)*scr1(n+10)
!            a(12) = glm*d2(mmm+1)*scr1(n+11)
!
!            a(13) = -glm*d2(mmm+2)*scr1(n+12)
!            a(14) = -glm*d2(mmm+2)*scr1(n+13)
!            a(15) = -glm*d2(mmm+2)*scr1(n+14)
!            a(16) = -glm*d2(mmm+2)*scr1(n+15)
!
            mmm = mmm+2
            nn = n+16
            n = n+16*l
!
            g = glm
!
            do 6 nnn = nn,n,16
               g = -g

               call alignx(16,scr1(nnn))
               call alignx(16,d2(mmm+1))
      	       call alignx(16,a(1))
	       
               reg00 = loadfp(scr1(nnn))
               reg01 = loadfp(scr1(nnn+2))
               reg02 = loadfp(scr1(nnn+4))
               reg03 = loadfp(scr1(nnn+6))
               reg04 = loadfp(scr1(nnn+8))
               reg05 = loadfp(scr1(nnn+10))
               reg06 = loadfp(scr1(nnn+12))
               reg07 = loadfp(scr1(nnn+14))

               reg16 = loadfp(d2(mmm+1))
	    	       
               reg08 = fxcpmadd(reg08, reg00, real(reg16))
               reg09 = fxcpmadd(reg09, reg01, real(reg16))
               reg10 = fxcpmadd(reg10, reg02, imag(reg16))
               reg11 = fxcpmadd(reg11, reg03, imag(reg16))
	       
               reg04 = fxpmul(reg04,g)
               reg12 = fxcpmadd(reg12,reg04,real(reg16))
               reg05 = fxpmul(reg05,g)
               reg13 = fxcpmadd(reg13,reg05,real(reg16))

               reg06 = fxpmul(reg06,g)
               reg14 = fxcpnmsub(reg14,reg06,imag(reg16))
               reg07 = fxpmul(reg07,g)
               reg15 = fxcpnmsub(reg15,reg07,imag(reg16))    	       	       	       	       	       	       

!               a(1) = a(1)+d2(mmm+1)*scr1(nnn)
!               a(2) = a(2)+d2(mmm+1)*scr1(nnn+1)
!               a(3) = a(3)+d2(mmm+1)*scr1(nnn+2)
!               a(4) = a(4)+d2(mmm+1)*scr1(nnn+3)
!               a(5) = a(5)+d2(mmm+2)*scr1(nnn+4)
!               a(6) = a(6)+d2(mmm+2)*scr1(nnn+5)
!               a(7) = a(7)+d2(mmm+2)*scr1(nnn+6)
!               a(8) = a(8)+d2(mmm+2)*scr1(nnn+7)
	       
!               a(9) = a(9)+g*d2(mmm+1)*scr1(nnn+8)
!               a(10) = a(10)+g*d2(mmm+1)*scr1(nnn+9)
!               a(11) = a(11)+g*d2(mmm+1)*scr1(nnn+10)
!               a(12) = a(12)+g*d2(mmm+1)*scr1(nnn+11)
!               a(13) = a(13)-g*d2(mmm+2)*scr1(nnn+12)
!               a(14) = a(14)-g*d2(mmm+2)*scr1(nnn+13)
!               a(15) = a(15)-g*d2(mmm+2)*scr1(nnn+14)
!               a(16) = a(16)-g*d2(mmm+2)*scr1(nnn+15)
               mmm = mmm+2
 6          continue
!
            mmmmmm = mmmmmm+mmmmm
            mmmmm = mmmmm-1
            mm = mmmmmm-m
!
            nn = 16*mm-15
!
            call alignx(16,scr2(nn))
	    
            call storefp(scr2(nn),reg08)
            call storefp(scr2(nn+2),reg09)
            call storefp(scr2(nn+4),reg10)
            call storefp(scr2(nn+6),reg11)
            call storefp(scr2(nn+8),reg12)
            call storefp(scr2(nn+10),reg13)
            call storefp(scr2(nn+12),reg14)
            call storefp(scr2(nn+14),reg15)	
!	       
!            scr2(nn) = a(1)
!            scr2(nn+1) = a(2)
!            scr2(nn+2) = a(3)
!            scr2(nn+3) = a(4)
!
!            scr2(nn+4) = a(5)
!            scr2(nn+5) = a(6)
!            scr2(nn+6) = a(7)
!            scr2(nn+7) = a(8)
!
!            scr2(nn+8) = a(9)
!            scr2(nn+9) = a(10)
!            scr2(nn+10) = a(11)
!            scr2(nn+11) = a(12)
!
!            scr2(nn+12) = a(13)
!            scr2(nn+13) = a(14)
!            scr2(nn+14) = a(15)
!            scr2(nn+15) = a(16)
 5       continue
 3    continue
      end subroutine rotatey

      subroutine performshift(nmultipoles,scr1,scr2,fr,sg)

      implicit none
      
      real(kind=8) :: scr1(*),scr2(*),fr(0:*),sg(0:*)
      integer(kind=8) :: nmultipoles
      real(kind=8), parameter  :: one = 1.0d0
      real(kind=8), parameter  :: zero = 0.0d0
      real(kind=8) :: g,gl,glm
      real(kind=8) :: a,aa,aaa,aaaa
      real(kind=8) :: b,bb,bbb,bbbb
      real(kind=8) :: c,cc,ccc,cccc      
      real(kind=8) :: d,dd,ddd,dddd      
      integer(kind=4) :: mmmmmm,mmmmm,mmmm,mmm,mm,m
      integer(kind=4) :: i,j,k,l,n,nn,nnn
      
      complex(kind=8) :: reg00,reg01,reg02,reg03,reg04,reg05,reg06,reg07
      complex(kind=8) :: reg08,reg09,reg10,reg11,reg12,reg13,reg14,reg15
      complex(kind=8) :: reg16      
      
!     perform shift
!
      i = -15
!
	       
      reg08 = (zero,zero)
      reg09 = (zero,zero)
      reg12 = (zero,zero)
      reg13 = (zero,zero)
      
!      a = zero
!      aa = zero
!      aaa = zero
!      aaaa = zero
!
!      c = zero
!      cc = zero
!      ccc = zero
!      cccc = zero
!
      do 7 j = 0,nmultipoles
        i = i+16
        call alignx(16,scr2(i))
	
        reg00 = loadfp(scr2(i))
        reg01 = loadfp(scr2(i+2))
        reg04 = loadfp(scr2(i+8))
        reg05 = loadfp(scr2(i+10))	       	       

        reg08 = fxcpmadd(reg08,reg00,fr(j))
        reg09 = fxcpmadd(reg09,reg01,fr(j))
        reg12 = fxcpmadd(reg12,reg04,fr(j))
        reg13 = fxcpmadd(reg13,reg05,fr(j))
		
!        a = a+fr(j)*scr2(i)
!        aa = aa+fr(j)*scr2(i+1)
!        aaa = aaa+fr(j)*scr2(i+2)
!        aaaa = aaaa+fr(j)*scr2(i+3)
!        c = c+fr(j)*scr2(i+8)
!        cc = cc+fr(j)*scr2(i+9)
!        ccc = ccc+fr(j)*scr2(i+10)
!        cccc = cccc+fr(j)*scr2(i+11)
 7    continue
!
      call alignx(16,scr1(1)) 
      call storefp(scr1( 1),reg12)
      call storefp(scr1( 3),reg13) 
      call storefp(scr1( 5),(zero,zero)) 
      call storefp(scr1( 7),(zero,zero))             
      call storefp(scr1( 9),reg08)
      call storefp(scr1(11),reg09)      
      call storefp(scr1(13),(zero,zero)) 
      call storefp(scr1(15),(zero,zero))             
      
!      scr1(1) = c
!      scr1(2) = cc
!      scr1(3) = ccc
!      scr1(4) = cccc
!
!      scr1(5) = zero
!      scr1(6) = zero
!      scr1(7) = zero
!      scr1(8) = zero
!
!      scr1(9) = a
!      scr1(10) = aa
!      scr1(11) = aaa
!      scr1(12) = aaaa
!
!      scr1(13) = zero
!      scr1(14) = zero
!      scr1(15) = zero
!      scr1(16) = zero
!
      i = 1
!
      do 8 l = 1,nmultipoles
        i = i+16*l
        j = -15
        k = nmultipoles+l
!

        reg08 = (zero,zero)
        reg09 = (zero,zero)
        reg12 = (zero,zero)
        reg13 = (zero,zero)
      
!        a = zero
!        aa = zero
!        aaa = zero
!        aaaa = zero
!c
!        c = zero
!        cc = zero
!        ccc = zero
!        cccc = zero
!c
        do 9 m = l,k
          j = j+16

          call alignx(16,scr2(j))

          reg00 = loadfp(scr2(j))
          reg01 = loadfp(scr2(j+2))
          reg04 = loadfp(scr2(j+8))
          reg05 = loadfp(scr2(j+10))

          reg08 = fxcpmadd(reg08,reg00,fr(m))
          reg09 = fxcpmadd(reg09,reg01,fr(m))
          reg12 = fxcpmadd(reg12,reg04,fr(m))
          reg13 = fxcpmadd(reg13,reg05,fr(m))	  
		  
!          a = a+fr(m)*scr2(j)
!          aa = aa+fr(m)*scr2(j+1)
!          aaa = aaa+fr(m)*scr2(j+2)
!          aaaa = aaaa+fr(m)*scr2(j+3)
!          c = c+fr(m)*scr2(j+8)
!          cc = cc+fr(m)*scr2(j+9)
!          ccc = ccc+fr(m)*scr2(j+10)
!          cccc = cccc+fr(m)*scr2(j+11)
 9      continue
!
        g = sg(l)
!
        call alignx(16,scr1(i))
	
        reg12 = fxpmul(reg12,g)
        call storefp(scr1(i),reg12)
        
        reg13 = fxpmul(reg13,g)
        call storefp(scr1(i+2),reg13) 
	
        call storefp(scr1(i+4),(zero,zero)) 
        call storefp(scr1(i+6),(zero,zero))             
        
        reg08 = fxpmul(reg08,g)
        call storefp(scr1(i+8),reg08)
	
        reg09 = fxpmul(reg09,g)
        call storefp(scr1(i+10),reg09)      
        
        call storefp(scr1(i+12),(zero,zero)) 
        call storefp(scr1(i+14),(zero,zero)) 	
	
!        scr1(i) = g*c
!        scr1(i+1) = g*cc
!        scr1(i+2) = g*ccc
!        scr1(i+3) = g*cccc
!
!        scr1(i+4) = zero
!        scr1(i+5) = zero
!        scr1(i+6) = zero
!        scr1(i+7) = zero
!
!        scr1(i+8) = g*a
!        scr1(i+9) = g*aa
!        scr1(i+10) = g*aaa
!        scr1(i+11) = g*aaaa
!
!        scr1(i+12) = zero
!        scr1(i+13) = zero
!        scr1(i+14) = zero
!        scr1(i+15) = zero
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
            reg08 = (zero,zero)
            reg09 = (zero,zero)
            reg10 = (zero,zero)
            reg11 = (zero,zero)
            reg12 = (zero,zero)
            reg13 = (zero,zero)
            reg14 = (zero,zero)
            reg15 = (zero,zero)
      
!            a = zero
!            aa = zero
!            aaa = zero
!            aaaa = zero
!c
!            b = zero
!            bb = zero
!            bbb = zero
!            bbbb = zero
!c
!            c = zero
!            cc = zero
!            ccc = zero
!            cccc = zero
!c
!            d = zero
!            dd = zero
!            ddd = zero
!            dddd = zero
!c
            do 12 mmmm = k,mmm
               nn = nn+16

               call alignx(16,scr2(nn))
	       
               reg00 = loadfp(scr2(nn))
               reg01 = loadfp(scr2(nn+2))
               reg02 = loadfp(scr2(nn+4))	       
               reg03 = loadfp(scr2(nn+6))
               reg04 = loadfp(scr2(nn+8))
               reg05 = loadfp(scr2(nn+10))	       	       
               reg06 = loadfp(scr2(nn+12))	       	       
               reg07 = loadfp(scr2(nn+14))   
	       	       	       
               reg08 = fxcpmadd(reg08,reg00,fr(mmmm))
               reg09 = fxcpmadd(reg09,reg01,fr(mmmm))
               reg10 = fxcpnmsub(reg10,reg02,fr(mmmm))	       
               reg11 = fxcpnmsub(reg11,reg03,fr(mmmm))	       	       	       

               reg12 = fxcpmadd(reg12,reg04,fr(mmmm))
               reg13 = fxcpmadd(reg13,reg05,fr(mmmm))
               reg14 = fxcpnmsub(reg14,reg06,fr(mmmm))	       
               reg15 = fxcpnmsub(reg15,reg07,fr(mmmm))	       	       	       
	       
!               a = a+fr(mmmm)*scr2(nn)
!               aa = aa+fr(mmmm)*scr2(nn+1)
!               aaa = aaa+fr(mmmm)*scr2(nn+2)
!               aaaa = aaaa+fr(mmmm)*scr2(nn+3)
!               b = b-fr(mmmm)*scr2(nn+4)
!               bb = bb-fr(mmmm)*scr2(nn+5)
!               bbb = bbb-fr(mmmm)*scr2(nn+6)
!               bbbb = bbbb-fr(mmmm)*scr2(nn+7)
!               c = c+fr(mmmm)*scr2(nn+8)
!               cc = cc+fr(mmmm)*scr2(nn+9)
!               ccc = ccc+fr(mmmm)*scr2(nn+10)
!               cccc = cccc+fr(mmmm)*scr2(nn+11)
!               d = d-fr(mmmm)*scr2(nn+12)
!               dd = dd-fr(mmmm)*scr2(nn+13)
!               ddd = ddd-fr(mmmm)*scr2(nn+14)
!               dddd = dddd-fr(mmmm)*scr2(nn+15)
 12         continue
!
            g = sg(k)
!
            call alignx(16,scr1(j))
	     
            reg12 = fxpmul(reg12,g)
            call storefp(scr1(j),reg12)
        
            reg13 = fxpmul(reg13,g)
            call storefp(scr1(j+2),reg13) 

            reg14 = fxpmul(reg14,g) 	
            call storefp(scr1(j+4),reg14) 

            reg15 = fxpmul(reg15,g)        	   
            call storefp(scr1(j+6),reg15)             
        
            reg08 = fxpmul(reg08,g)
            call storefp(scr1(j+8),reg08)
	
            reg09 = fxpmul(reg09,g)
            call storefp(scr1(j+10),reg09)      

            reg10 = fxpmul(reg10,g)        
            call storefp(scr1(j+12),reg10) 

            reg11 = fxpmul(reg11,g)        	    
            call storefp(scr1(j+14),reg11)
		     
!            scr1(j) = g*c
!            scr1(j+1) = g*cc
!            scr1(j+2) = g*ccc
!            scr1(j+3) = g*cccc
!
!            scr1(j+4) = g*d
!            scr1(j+5) = g*dd
!            scr1(j+6) = g*ddd
!            scr1(j+7) = g*dddd
!
!            scr1(j+8) = g*a
!            scr1(j+9) = g*aa
!            scr1(j+10) = g*aaa
!            scr1(j+11) = g*aaaa
!
!            scr1(j+12) = g*b
!            scr1(j+13) = g*bb
!            scr1(j+14) = g*bbb
!            scr1(j+15) = g*bbbb
 11      continue
         n = n+mm
         mm = mm-16
 10   continue
      end subroutine performshift
