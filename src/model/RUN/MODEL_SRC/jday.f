c get a day in the format yyyymmdd as first argument
c and return a date in the same format plus or minus
c a number of days given as second argument.
      
      program jjday
      integer nday,date,dd,mm,yyyy,jday
      character*8 arg1,arg2

      indx = iargc( )
      if(indx.ne.2)then
        print*,' Usage: jday date increment'
        print*,'        date      must be in YYYYMMDD format'
        print*,'        increment must be a + or - integer'
        stop
      endif
      
      call getarg(1,arg1)
      call getarg(2,arg2)
      
      read(arg1,'(I8)') date
      read(arg2,'(I8)') nday
      
      yyyy = date/10000
      mm   = mod(date,10000)/100
      dd   = mod(date,100)

      jday = julday(mm,dd,yyyy)
     
      jday = jday + nday
      call caldat(jday,mm,dd,yyyy)
      
      date = yyyy*10000 + mm*100 + dd
      write(*,'(I8)') date
    
      end
      
      
      
      function julday(mm,id,iyyy)
      integer julday,id,iyyy,mm,IGREG
      parameter (IGREG=15+31*(10+12*1582))
      integer ja,jm,jy
      jy=iyyy
c      if (jy.eq.0) pause 'julday: there is no year zero'
      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
          jm=mm+1
      else
          jy=jy-1
	  jm=mm+13
      endif
      julday=365*jy+int(0.25d0*jy+2000.d0)+int(30.6001d0*jm)+id+1718995
      if (id+31*(mm+12*iyyy).ge.IGREG) then
          ja=int(0.01d0*jy)
	  julday=julday+2-ja+int(0.25d0*ja)
      endif
      return
      end
      
      subroutine caldat(julian,mm,id,iyyy)
      integer id,iyyy,julian,mm,IGREG
      parameter (IGREG=2299161)
      integer ja,jalpha,jb,jc,jd,je
      if (julian.ge.IGREG) then
          jalpha=int(((julian-1867216)-0.25d0)/36524.25d0)
          ja=julian+1+jalpha-int(0.25d0*jalpha)
      else if (julian.lt.0) then
          ja=julian+36525*(1-julian/36525)
      else
          ja=julian
      endif
      jb=ja+1524
      jc=int(6680.0d0+((jb-2439870)-122.1d0)/365.25d0)
      jd=365*jc+int(0.25d0*jc)
      je=int((jb-jd)/30.6001d0)
      id=jb-jd-int(30.6001d0*je)
      mm=je-1
      if (mm.gt.12) mm=mm-12
      iyyy=jc-4715
      if (mm.gt.2) iyyy=iyyy-1
      if (iyyy.le.0) iyyy=iyyy-1
      if (julian.lt.0) iyyy=iyyy-100*(1-julian/36525)
      return
      end
