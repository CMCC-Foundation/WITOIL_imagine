      character mode*2,day*2,month*2,year*4,hour*4
      character duration*4,name*3
      real*8 :: dist,splon,splat,degrad,dep
      real alon1,alon2,alat1,alat2
      real lat_degree,lat_minute,lon_degree,lon_minute
      real length,step
      integer numfiles
      integer day_int,count
      open(100,file='medslik5.inp',status='old')
      read(100,*) mode
      read(100,*) restart
      read(100,*) day,month,year
      read(100,*) hour
      read(100,*) duration
      read(100,*) lat_degree,lat_minute
      read(100,*) lon_degree,lon_minute
      read(100,*) name
      read(100,*) length
      read(100,*) step
      read(100,*) icurrents
C      close(100)
      dist=6.5d0*length
      dist=dist/60
      splon=lon_degree+lon_minute/60.d0
      splat=lat_degree+lat_minute/60.d0
      pi=4.*datan(1.d0)
      degrad=180.d0/pi
      dep=dist/dcos(splat/degrad)
C      print *, dist,splon,splat,dep,degrad
      alon1=splon-dep
      alon2=splon+dep
      alat1=splat-dist
      alat2=splat+dist
      numfiles=length/24+2
C      prdate=year(3:4)//month//day//'24'
      read(day,*) day_int
      open(101,file='medslik1.tmp')
      write(101,*) alon1,'   ',alon2,'        ','Min & Max Longitudes'
      write(101,*) alat1,'   ',alat2,'        ','Min & Max Latitudes'
           
      end
