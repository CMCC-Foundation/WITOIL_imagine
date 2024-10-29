c
c     sets codes for de-crypting map and bath files	  
c
      subroutine setcodes(regn1)
	  implicit real*8(a-h,o-z)
      character regn1*4
	  common /encr/ aloncorr,ia,ib,ic

      if(regn1.eq.'mltc'.or.regn1.eq.'mlts') then
        ia = 9331220          
        ib = 19380209
        ic = 6791
        aloncorr = 0.d0

      else if(regn1.eq.'medf') then
        ia = 6340856            
        ib = 13849450
        ic = 1267
        aloncorr = 8.13579d0
      
      else if(regn1.eq.'sici') then
        ia = 9047211            
        ib = 22851092
        ic = 5003
        aloncorr = 0.d0
      
      else if(regn1.eq.'syri') then
        ia = 7400116            
        ib = 16387024
        ic = 4680
        aloncorr = 0.d0
      
      else if(regn1.eq.'emed') then
        ia = 32180772            
        ib = 25770147
        ic = 2055
        aloncorr = 0.d0
      
      else if(regn1.eq.'cyba'.or.regn1.eq.'cyse') then
        ia = 25410893            
        ib = 3176880
        ic = 3855
        aloncorr = 0.d0
      
c      else if(regn1.eq.'adri'.or.regn1.eq.'adno'.or.
c     &         regn1.eq.'anza') then
c        ia = 6682015            
c        ib = 15882209
c        ic = 5601
c        aloncorr = 0.d0
      
      endif

      return
      end

       
