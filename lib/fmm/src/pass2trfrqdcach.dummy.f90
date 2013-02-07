#include "fmm.h"

subroutine rotateback(nmultipoles,scr1,scr2,d3f,jaddress,cmphi,smphi,cmphipi,smphipi)

implicit none

real(kind=8) :: scr1(*),scr2(*),d3f(*)
real(kind=8) :: cmphi(0:*),smphi(0:*)
real(kind=8) :: smphipi(0:*),cmphipi(0:*)

integer(kind=8) :: nmultipoles,jaddress(*)

real(kind=8), parameter  :: one = 1.0d0
real(kind=8), parameter  :: zero = 0.0d0
real(kind=8) :: g,gl,glm
real(kind=8), dimension(16) :: a
integer(kind=4) :: mmmmmm,mmmmm,mmmm,mmm,mm,m
integer(kind=4) :: i,j,k,l,n,nn,nnn

complex(kind=8) :: reg00,reg01,reg02,reg03,reg04,reg05,reg06,reg07
complex(kind=8) :: reg08,reg09,reg10,reg11,reg12,reg13,reg14,reg15
complex(kind=8) :: reg16 


      i = 1
      j = 1

      gl = one

      do 13 l = 1,nmultipoles
         gl = -gl
         i = i+1
         j = j+16*l
         n = j

         mmm = mmm+1
     
         call alignx(16,scr1(n))

         reg00 = loadfp(scr1(n))
         reg01 = loadfp(scr1(n+2))
         reg04 = loadfp(scr1(n+8))
         reg05 = loadfp(scr1(n+10))
     
         reg08 = fxpmul(reg00,d3f(mmm))
         reg09 = fxpmul(reg01,d3f(mmm))     

         reg04 = fxpmul(reg04,gl)
         reg12 = fxpmul(reg04,d3f(mmm))     
         reg05 = fxpmul(reg05,gl)
         reg13 = fxpmul(reg05,d3f(mmm))     

#ifdef FMM_ALLOCALIGNED
         mmm = mmm+1
#endif

         nn = n+16
         n = n+16*l

         g = gl

 
         call storefp(a(1),reg08) 
         call storefp(a(3),reg09)

         call storefp(a(9),reg12)
         call storefp(a(11),reg13)          

         scr2(jaddress(1)+i) = a(1)
         scr2(jaddress(2)+i) = a(2)
         scr2(jaddress(3)+i) = a(3)
         scr2(jaddress(4)+i) = a(4)

         scr2(jaddress(5)+i) = zero
         scr2(jaddress(6)+i) = zero
         scr2(jaddress(7)+i) = zero
         scr2(jaddress(8)+i) = zero

         scr2(jaddress(9)+i) = a(9)
         scr2(jaddress(10)+i) = a(10)
         scr2(jaddress(11)+i) = a(11)
         scr2(jaddress(12)+i) = a(12)

         scr2(jaddress(13)+i) = zero
         scr2(jaddress(14)+i) = zero
         scr2(jaddress(15)+i) = zero
         scr2(jaddress(16)+i) = zero

         glm = gl

         do 15 m = 1,l
            glm = -glm
            i = i+1
            n = j

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

            mmm = mmm+2
            nn = n+16
            n = n+16*l

            g = glm


            call storefp(a( 1),reg08)
            call storefp(a( 3),reg09)
            call storefp(a( 5),reg10)
            call storefp(a( 7),reg11)
            call storefp(a( 9),reg12)
            call storefp(a(11),reg13)
            call storefp(a(13),reg14)
            call storefp(a(15),reg15)

                                                    
            scr2(jaddress(1)+i) = cmphipi(m)*a(1)-smphipi(m)*a(5)
            scr2(jaddress(2)+i) = cmphipi(m)*a(2)-smphipi(m)*a(6)
            scr2(jaddress(3)+i) = cmphipi(m)*a(3)-smphipi(m)*a(7)
            scr2(jaddress(4)+i) = cmphipi(m)*a(4)-smphipi(m)*a(8)
            scr2(jaddress(5)+i) = cmphipi(m)*a(5)+smphipi(m)*a(1)
            scr2(jaddress(6)+i) = cmphipi(m)*a(6)+smphipi(m)*a(2)
            scr2(jaddress(7)+i) = cmphipi(m)*a(7)+smphipi(m)*a(3)
            scr2(jaddress(8)+i) = cmphipi(m)*a(8)+smphipi(m)*a(4)
            scr2(jaddress(9)+i) = cmphi(m)*a(9)-smphi(m)*a(13)
            scr2(jaddress(10)+i) = cmphi(m)*a(10)-smphi(m)*a(14)
            scr2(jaddress(11)+i) = cmphi(m)*a(11)-smphi(m)*a(15)
            scr2(jaddress(12)+i) = cmphi(m)*a(12)-smphi(m)*a(16)
            scr2(jaddress(13)+i) = cmphi(m)*a(13)+smphi(m)*a(9)
            scr2(jaddress(14)+i) = cmphi(m)*a(14)+smphi(m)*a(10)
            scr2(jaddress(15)+i) = cmphi(m)*a(15)+smphi(m)*a(11)
            scr2(jaddress(16)+i) = cmphi(m)*a(16)+smphi(m)*a(12)
 15      continue
 13   continue

end subroutine rotateback
