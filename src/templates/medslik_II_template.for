cc-----------------------------------------------------------------------------------
c  MEDSLIK-II_2.02
c  oil spill fate and transport model 
c-----------------------------------------------------------------------------------
c  medslik_II.for
c  This program simulates the tranport and weathering of an oil spill.
c-----------------------------------------------------------------------------------
c  Copyright (C) <2012>
c  This program was originally written
c  by Robin Lardner and George Zodiatis.
c  Subsequent additions and modifications
c  have been made by Michela De Dominicis. 
c----------------------------------------------------------------------------------
c  The development of the MEDSLIK-II model is supported by a formal agreement
c  Memorandum of Agreement for the Operation and Continued Development of MEDSLIK-II
c  signed by the following institutions:
c  INGV - Istituto Nazionale di Geofisica e Vulcanologia
c  OC-UCY - Oceanography Center at the University of Cyprus
c  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per 
c  lo Studio dell’Ambiente Marino Costiero
c  CMCC - Centro Euro-Mediterraneo sui Cambiamenti Climatici
c 
c  This program is free software: you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation, either version 3 of the License, or
c  any later version.
c
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c  You should have received a copy of the GNU General Public License
c  along with this program.  If not, see <http://www.gnu.org/licenses/>.
c-----------------------------------------------------------------------------------

      implicit real*8(a-h,o-z)
	character dummy*80

c----------------------------------------------------------------------------
c	DEFITION OF BASIC NETCDF PARAMETERS
C----------------------------------------------------------------------------
      include 'netcdf.inc'
C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      parameter (FILE_NAME = 'spill_properties.nc')
      integer ncid
C     Setting up dimensions --- n_parcel and time
      integer NDIMS
      parameter (NDIMS = 2)
      character*(*) PARCEL_NAME, TIME_NAME
      parameter (PARCEL_NAME = 'parcel_id')
      parameter (TIME_NAME = 'time')
C     Error handling.
      integer nc_status	

c     Set the common block
      common /ncfile/ ncid,nc_status

C     Create the file. 
      nc_status = nf_create(FILE_NAME, NF_NETCDF4, ncid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)


C----------------------------------------------------------------------
C	NETCDF FILE HAS BEEN CREATED
C-----------------------------------------------------------------------
c----------------------------------------------
      open(1,file='config2.txt')
	do k=1,43
	  read(1,'(a80)') dummy
	enddo
c     if(.not.eof(1)) read(1,'(a80)') dummy
      read(1,'(a80)') dummy
        if(dummy(1:6).eq.'Comput') then
	  read(1,*) nstph
        read(1,*) npl
        delt = 1.d0 / dfloat(nstph) 
      else
	  delt = 0.5d0
        npl = 2000
      endif
      close(1)     
c      write(6,*) delt,npl
           
      call main(delt,npl)
      
      end    
c----------------------------------------------------------------------
c     begin the main program
c----------------------------------------------------------------------
      
      subroutine main(delt, npl)
      implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX,
     &          ntm=2000,npc=100000,nss=350000,msp=1200)
c     &          ntm=2000,npc=100000,nss=350000,npl=4000,msp=1200)

c----------------------------------------------------------------------
c     dimension declarations for environment
c----------------------------------------------------------------------
      dimension 
     &      h(mm,nm),itype(mm,nm),uadv(mm,nm),vadv(mm,nm)
c
      dimension 
     &      wndtim(0:ntm),wndvel(0:ntm),wnddir(0:ntm),
     &      wcx(0:11,50),wcy(0:11,50),wcu(0:11,50),wcv(0:11,50),
     &      winx(mm,nm),winy(mm,nm),wdrftx(mm,nm),wdrfty(mm,nm),
     &      stoku(mm,nm),stokv(mm,nm),wavenum(mm,nm),
     &      wvel_vec(720),wvel_mean(720),wdir_vec(720),wdir_mean(720),
     &      freq_sp(700),fang(700),erre(700),sp_exp1(700), spectra(700),
     &      wave_num(700),stoke_sp(700),hwave_d(700),stoke_d(700),
     &      sp_exp2(700)
c
      dimension 
     &      nwc(0:10),uu(24),vv(24),ss(3),
     &      curtim(0:ntm),curvel(0:ntm),curdir(0:ntm),
     &      crntim(5),crnu(5),crnv(5),
     &      ifcstfn(720),fcsttim(720),iwfcstfn(30),wfcsttim(30,24)
c     
      dimension 
     &      bmtim(20),bmlat1(20),bmlon1(20),bmlat2(20),bmlon2(20),
     &      bmx1(20),bmy1(20),bmx2(20),bmy2(20),ibmeff(20)
c----------------------------------------------------------------------
c     dimension declarations for oil slick
c----------------------------------------------------------------------
      dimension 
     &      blat_filt(10000),blon_filt(10000),
     &      c1p(npc),c2p(npc),px(npc),py(npc),pz(npc),poll(npl,npl),
     &      is(npc),ib(npc),itmp(npc),bd1(nss,4),ibd(nss),
     &      seg(nss,4),alngt(nss),segx(nss),segy(nss),ns0(nss),
     &      seep(nss),prel(nss),sfrac(nss),vcst(nss),
     &      den(msp),vis(msp),visem(msp),tre(msp),c1ms(msp),c2ms(msp),
     &      atn(msp),atk(msp),ato(msp),ttk(msp),tto(msp),
     &      vtn(msp),vtk(msp),vto(msp),xcl(msp),xcs(msp),
     &      vtne(msp),vtke(msp),vte(msp),vtnd(msp),vtkd(msp),vtd(msp),
     &      ftk(msp),ftn(msp),fw(msp),pcte(msp),pctd(msp)
c----------------------------------------------------------------------
c     common declarations
c----------------------------------------------------------------------
      common /blk1/ wndtim,wndvel,wnddir,nwind
      common /blk2/ wcx,wcy,wcu,wcv,nwc
      common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad
      common /blk4/ along1,alatg1,along2,alatg2,dlong,dlatg
      common /blk5/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2
      common /encr/ aloncorr,iaencr,ibencr,icencr
c
      common /fcst/ sst(mm,nm),wx(mm,nm),wy(mm,nm),usrf(mm,nm),
     &      vsrf(mm,nm),u10(mm,nm),v10(mm,nm),u30(mm,nm),v30(mm,nm),
     &      u120(mm,nm),v120(mm,nm),dsst(mm,nm),dwx(mm,nm),dwy(mm,nm),
     &      dusrf(mm,nm),dvsrf(mm,nm),du10(mm,nm),dv10(mm,nm),
     &	  du30(mm,nm),dv30(mm,nm),du120(mm,nm),dv120(mm,nm)
c----------------------------------------------------------------------
c     common declarations for oil & spill characteristics
c----------------------------------------------------------------------
	common /spill/  idd,imm,iyr,ispill,tstart,tcomp,x0,y0
      common /evap/ ce1,ce,vappr,fmaxe
      common /disp/ cd1,cd3,cd4,cd5,vl,vs1,um,stk,stn,fmaxd
      common /emul/ cm1,cm2,cm3,visemx
      common /sprd/ cs1,cs2,cs3
      common /phys/ deno,denk,cdt,viso,visk,cvt,tvk0,denw,tk,tk0,wvel
c---------------------------------------------------------------------
c     common declarations of the netcdf file
c---------------------------------------------------------------------
      common /ncfile/ ncid,nc_status
c----------------------------------------------------------------------
c     character declarations
c---------------------------------------------------------------------- 
      character empty*80, a0(3)*4, a1*1, a2*2, a3*3, a4*4, 
     &          ay*2, am*2, ad*2, ah*2, d24*3, d06*3, vlunit*4, pref*3,
     &          nore*2, nora*2,ora*2,ore*2,indate(30)*8,dummy*80 
      character fn(3)*11,fnsim(30)*17, fnadm(30)*22, fncym(30)*19
      character seas1*1, seas2*1, fcstcurdir*14,list*11 
      character fcstfn(720)*16, wfcstfn(30)*14, date*10
      integer nt,ix,count
      
      logical ex
c----------------------------------------------------------------------
c	DECLARING NETCDF VARIABLE AND DIMENSION IDS
c----------------------------------------------------------------------
C     We recommend that each variable carry a "units" attribute.
c     Oil density and number of parcels was also added for 
c     further calculations
      character*(*) UNITS,OIL_DENSITY,PARCEL_VOL
      character*(*) SP_LON,SP_LAT
      parameter (UNITS = 'units')
      parameter (OIL_DENSITY = 'oil_density')
      parameter (PARCEL_VOL = 'volume_of_parcel')
      parameter (SP_LON = 'initial_position_x')
      parameter (SP_LAT = 'initial_position_y')
      character*(*) TIME_UNITS, EVOL_UNITS, LAT_UNITS, LON_UNITS
      character*(*) NVOL_UNITS, PRTS_UNITS, WC_UNITS,TFXD_UNITS
      character*(*) VEM_UNITS,VIS_UNITS,DEN_UNITS,VOLR_UNITS  
      character*(*) SEG_UNITS
      parameter (SEG_UNITS = 'segment_id')
      parameter (TIME_UNITS = 'hours')
      parameter (EVOL_UNITS = 'm3')
      parameter (NVOL_UNITS = 'm3')             ! check that 
      parameter (PRTS_UNITS = '0,1,2,3,4,-nsg,9')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')
      parameter (WC_UNITS = 'percentage')
      parameter (TFXD_UNITS = 'tonnes')
      parameter (VEM_UNITS = 'Pa.s')
      parameter (VIS_UNITS = 'Pa.s')
      parameter (DEN_UNITS = 'kg/m3')
      parameter (VOLR_UNITS = 'percentage')
 
C     The start and count arrays will tell the netCDF library where to
C     write our data.	
      integer NPARCEL,NTIME
      integer NDIMS
      parameter (NDIMS = 2)
      integer prcl_dimid, time_dimid,nc_time
      integer ncpointer(NDIMS),ncounter(NDIMS),dimids(NDIMS)
      REAL x_coordinate,y_coordinate
c----------------------------------------------------------------------
c     transformations betweenthe medslik grid coords and lat/lon
c----------------------------------------------------------------------
      xgrid(alon)=(alon-along1)/dlong+1.d0
      ygrid(alat)=(alat-alatg1)/dlatg+1.d0
      glon(x)=along1+(x-1.d0)*dlong
      glat(y)=alatg1+(y-1.d0)*dlatg

c----------------------------------------------------------------------
c----------------------------------------------------------------------
	data hmin /1.d0/,    hmax /2500.d0/
      data a0 /'.cst','.srf','.dsp'/

      pi=4.*datan(1.d0)
      degrad=180.d0/pi
      open(90,file='medslik.log')
c      open(91,file='smag.log')
c-----------------------------------------------------------------------
c     list of input & output files and unit numbers
c-----------------------------------------------------------------------
c      open(39,file='medslik5.par',status='old')
c      open(40,file='medslik5.inp',status='old')
c      open(41,file='medslik.wnd',status='old')
c      open(43,file='medslik.crr',status='old')
c      open(44,file='medslik.nuc',status='old')
c      open(46,file='medslik.crn',status='old')
c      open(48,file='medslik.bms',status='old')
c      open(30,file='flag.tmp',status='old'
c      open(50,file='/DTM_INP/dtm.bath',status='old')
c      open(51,file='/DTM_INP/dtm.map',status='old')
c      open(52,file='/DTM_INP/dtmcst1.d',dispose='delete')
c      open(71,file=fcstcurdir//fcstfile)
c      open(72,file=fcstwinddirr//fcstfile)
c
c      open(90,file='medslik.log')
c      open(99,file='output/medslik.fte')
c
c      open(38,file='initial.txt') satellite data file
c---------------------------------------------------------------------------------
c      model options (M. De Dominicis)

c----------------------------------------------------------------------
c     model parameters        ! default values
c     istoke= 0:no stokes drift calculation; 1: stoke drift calculation using Jonswap spectrum (M. De Dominicis)
c     alpha, beta = empirical parameters for wind drift of slick
c     ibetared = 1 if beta reduces by 50% at wind speed = halfspeed
c     beta=beta0*(1.-ibetared*wvel/(wvel+halfspeed))
c     iwindred = 1 if fraction of forecast wind is subtracted in drift formula
c     wredfrac = fraction of fcst wind to be subtracted when fcst currents are used 
c     ismag = 1 if Smagorinski scheme is used for horiz diffusivity
c     horizk = horizontal diffusivity
c     vertk1,2 = vertical diffusivities above & below thermocline
c     thermocl = depth of thermocline
c     ntot = no of parcels used to model diffusion and dispersion
c     fcstdepth1,2,3 = depths at which currents are printed in forecast files
c----------------------------------------------------------------------
      open(39,file='config2.txt',status='old')
      read(39,*) empty 
      read(39,*) istoke        ! 01
      print *, 'istoke',istoke
      read(39,*) alpha        ! 0.031
	read(39,*) beta0        ! 13.
	read(39,*) ibetared     ! 00
	read(39,*) halfspeed    ! 0.0
	read(39,*) iwindred     ! 00
	read(39,*) wredfrac     ! 0.0
	read(39,*) ismag        ! 00
	read(39,*) horizk       ! 2.0
	read(39,*) vertk1       ! 0.01
	read(39,*) vertk2       ! 0.0001
	read(39,*) thermocl     ! 30.0
	read(39,*) ntot         ! 10000
	read(39,*) fcstdep1,fcstdep2,fcstdep3     ! 10.,30.,120.
	read(39,*) idepth       ! 02
	if(ntot.gt.100000) then
	  write(6,*) 'Total number of parcels cannot exceed 100,000'
	  write(6,*) 'Return to Parameters Form and change this number'
c	  pause
	  stop
	endif
c----------------------------------------------------------------------
c     evaporation constants from Mackay et al
c     ce=coeff accounts for drop in vapour pressure with evaporation (ce=10-20)
c          ce1=akew*(wvel*3.6)**gamma_evap = evaporative exposure to wind
c     visk = coeff for change of viscosity with evaporation
c----------------------------------------------------------------------
      read(39,*) empty
      read(39,*) ce     ! 12.0
      read(39,*) akew   ! 0.000033
      read(39,*) gamma_evap  ! 0.78
      read(39,*) visk   ! 4.
c----------------------------------------------------------------------
c     emulsion constants from Mackay et al
c     cm1 controls the effect of water fraction on mousse viscosity
c     cm2 controls the rate of water absorption
c     cm3 controls maximum water fraction in the mousse (decreased for heavy oils)
c     icm3 = 1 if max water fraction increases/decreases for light/heavy oils
c     visemx = maximum mousse viscosity after which emulsification stops
c----------------------------------------------------------------------
      read(39,*) empty
      read(39,*) cm1     ! 0.65
      read(39,*) cm2     ! 1.6e-06
      read(39,*) cm3     ! 1.333
      read(39,*) icm3    ! 1
      visemx=100000.d0
c----------------------------------------------------------------------
c     dispersion constants from Mackay et al
c     cd1=downward diffusion velocity of small droplets (m/s)
c     cd3=controls the rate of dispersion of all droplets by waves
c     cd4=controls the fraction of droplets below the critical size
c     cd5=controls the dispersion from the thin slick (sheen)
c     vs1=rising velocity of small droplets (m/s)
c     vl=rising velocity of large droplets (m/s)
c     um=controls depth of well-mixed surface layer (m)
c     st=interfacial tension between oil and water
c     fmaxd=max dispersive fraction
c----------------------------------------------------------------------
      read(39,*) empty
      read(39,*) cd1     ! 0.001
      read(39,*) cd3     ! 0.8e-05
      read(39,*) cd4     ! 50.0
      read(39,*) cd5     ! 2000.0
      read(39,*) vl      ! 0.08
      read(39,*) vs1     ! 0.0003
      read(39,*) um      ! 0.5
      read(39,*) st      ! 24.0
      read(39,*) fmaxd   ! 1.0
      stk=st
      stn=st
c----------------------------------------------------------------------
c     spreading constants from Mackay et al
c     sprdmx = max time for spreading of spill (hours)
c     for new parameter files read the following:
c     seepmx = limiting concentration on coast (bbls/km)
c     apicoeff = coeff of reduction of coastalretention rate for heavy oils
c----------------------------------------------------------------------
      read(39,*) empty
      read(39,*) cs1     ! 1.0
      read(39,*) cs2     ! 150.0
      read(39,*) cs3     ! 0.0015
      sprdmx=24.d0

	data seepmx /5000.d0/
C       if(.not.eof(39)) then
        read(39,*) empty
        read(39,*) seepmx     ! 5000.0
        read(39,*) apicoeff   ! 0.0
C	else
C	  seepmx = 1.d10
C	  apicoeff = 0.d0
C	endif
	close(39)
c----------------------------------------------------------------------
c     delt = time increment for spill re-computations (hrs)
c	stph, nstph = no of steps per hour
c----------------------------------------------------------------------
c      delt=0.5d0
        deltmn=delt*60.d0
	deltsc=delt*3600.d0
	stph=1.d0/delt
	nstph=stph+0.001d0

c----------------------------------------------------------------------
c     mode = 0/1/2 for trajectory/oil spill/pollutant transport simulation
c     numspills = number of spills (mode = 1 or 2 only)
c     irestart = 1 if a previous run is to be restarted
c     ihrestart = number of hours of this previous run
c----------------------------------------------------------------------



      icencr = 1 ! initialization done in 2.01; bug in previous versions?

      open(40,file='medslik5.inp',status='old')


      read(40,*) mode,numspills
      if(mode.ne.1) then
	  write(6,*) '   The file medslik5.inp has incorrect mode.'
	  write(6,*) 'Modify this file using the input interface if you'
	  write(6,*) '       wish to make a spill simulation,'
        stop
      endif 
      if(numspills.ne.1) then
	  write(6,*) 'File medslik5.inp is for multiple spill simulation'
	  write(6,*) 'Modify this file using the input interface if you'
	  write(6,*) '       wish to simulate a single spill,'
        stop
      endif 

      read(40,*) irestart,ihrestart
	if(irestart.eq.0) ihrestart=0
c----------------------------------------------------------------------
c     idd/imm/iyr = date of spill
c     istart = time of day at start of spill
c     tstart = nearest hour to start of spill
c     tspill (ispill) = duration of spill (hours)
c     splrte = spill rate (bbls per hour)
c     tcomp (icomp) = duration of computation from start of spill (hrs)
c     lat,alat, lon,alon = geographical location of spill (deg,min)
c     pref = 3-letter prefix for labelling output files
c	iprs = interval for output (hrs)
c	icrn = 1/0 if observation of spill posn is (is not) to be applied
c----------------------------------------------------------------------
	  
      read(40,'(i2,1x,i2,1x,i4)') idd,imm,iyr

      read(40,'(i4)') istart
      read(40,'(i4)') ispill
      read(40,*) lat,alat
      read(40,*) lon,alon


      read(40,'(a3)') pref
      read(40,'(i4)') icomp
      read(40,'(i3)') iprs

      read(40,'(i2)') icrn

      tstart=dfloat(istart/100)+dfloat(istart-100*(istart/100))/60.d0
      tspill=dfloat(ispill)
      tcomp=dfloat(icomp)
      splat=dfloat(lat) + alat/60.d0
      splon=dfloat(lon) + alon/60.d0

c----------------------------------------------------------------------
c     vlunit = unit of volume (tons, cu.m, bbls, gals)
c     splrte = spill rate (no units per hour) or total volume for ispill=0 
c     api = API number of the oil
c     deno = oil density  (g/cm**3)
c     den2 = density of residual part (g/cm**3)
c     respc = residual percentage
c     viso = initial oil viscosity
c     tem0 = reference temperature for viscosity
c     vappr = oil vapour pressure from file
c     max water content reduced/increased for heavy/light oils
c----------------------------------------------------------------------
      read(40,'(f9.2)') splrte
      read(40,'(a4)') vlunit
      print *, 'splrte', splrte
c---------------------------------------------------------------------------------
c M. De Dominicis
c      isat= 0:point source; 1: areal source of spill (from manual slick contour
c      or from satellite data)
c      iage= age of the slick (in hours, 0, 24, 48)
c---------------------------------------------------------------------------------
      read(40,*) iage       
      read(40,*) isat
      print *, 'isat', isat, 'iage', iage
      read(40,'(a80)') empty
      read(40,'(f5.2)') api
      read(40,'(f5.3)') deno
      read(40,'(f5.3)') den2
      read(40,'(f5.2)') respc
      read(40,'(f5.1)') viso
      read(40,'(f5.1)') tem0
      read(40,'(f5.1)') vappr

      
      tvk0=tem0+273.d0
      if(icm3.eq.1.and.cm3.eq.1.333d0) then
        cm3=(10.d0/9.d0)*(2.d0-1.d0/(1.d0+4.d0**(1.d0-api/17.d0)))
	end if
      
	po=vappr
c----------------------------------------------------------------------
c	isst = index that specifies source of SST data
c	tc = SST in deg C - read only if isst = 8, otherwise zero
c	ibooms = 1/0 if booms are (are not) deployed
c     al5 = output pixel size (m) for counting up slick parcels
c----------------------------------------------------------------------
      read(40,'(i2)') isst
      read(40,'(f4.1)') tc
      read(40,'(i2)') ibooms
      read(40,'(f5.1)') al5

      gridkm=al5/1000.d0
      area=gridkm*gridkm
c----------------------------------------------------------------------
c	fcstfn() = name of forecast file of currents
c	ifcstfn() = 1/0 if the file is (is not) available
c	wfcstfn() = name of forecast file of winds
c	iwfcstfn() = 1/0 if the wind file is (is not) available
c----------------------------------------------------------------------
	  iprod=1
	  read(40,'(i4)') nfcst  
c----------------------------------------------------------------------
c       Hourly forecast files (M. De Dominicis)
c----------------------------------------------------------------------
        open(41,file='medslik.tmp')
        do n=1,2
        read(41,*) dummy
        enddo
        read(41,*) numfiles
        do n=1,numfiles+1
          read(41,'(a8)') indate(n)
        enddo
        close(41)

	  do i=1,nfcst
	 
	        do k=1,3
		read(40,'(a11,i2)') fn(k),ifcstfn(i)
		enddo

c 	date list for hourly currents files 
c       starting from 01:00

                list=fn(1)
		do nt=1,24
        	 if(nt.lt.10) then
        	 write(nore,'(i2)') nt
		 nora='0'//nore(2:2)
		 else 
		 write(nora,'(i2)') nt
		 end if 
  
                if(nt.le.12) then
                 count=i
                write(ore,'(i2)') nt+12
                ora=ore(1:2)
                endif
                if(nt.gt.12.and.nt.lt.22) then
                 count=i+1
                write(ore,'(i2)') nt-12
                ora='0'//ore(2:2)
                endif
                if(nt.ge.22.and.nt.le.24) then
                 count=i+1
                write(ore,'(i2)') nt-12
                ora=ore(1:2)
                endif

	        fcstfn((i-1)*24+nt)='merc'//list(1:6)//nora(1:2)//'.mrc'

	  		ifcstfn((i-1)*24+nt)=ifcstfn(i)	
	        enddo
	  
		iprod=iprod*ifcstfn(i)
		write(6,'(a16,i2)') fcstfn(i),ifcstfn(i)

	  end do

	  if(iprod.eq.0) then
	    write(6,*) '*************************************************'
	    write(6,*) 'WARNING: The following forecast files are missing'
	    do i=1,nfcst
	  if(ifcstfn(i).eq.0)then
	      write(6,*) fn(1),', ',fn(2),', ',fn(3)
	    write(6,*) fcstfn(i)
	      end if
	    end do
	    write(6,*) 'The simulation will extrapolate from the last'
	    write(6,*) 'available file and so will lose reliability.'
	    write(6,*) '*************************************************'
	  end if
      
        if(ismag.eq.1) then
	    write(6,*) '*************************************************'
	    write(6,*) 'Smagorinsky Diffusivity Model can ONLY be used'
	    write(6,*) '    when FCORECAST Water Currents are used.'
	    write(6,*) 'Return to the Startup Screen and enter a value' 
	    write(6,*) '         for Horizontal Diffusivity' 
	    write(6,*) '*************************************************'
	    stop
        end if   
      
	
	  iprod=1
	  read(40,'(i3)') nwfcst
	  do i=1,nwfcst
		read(40,*) wfcstfn(i),iwfcstfn(i)

		iprod=iprod*iwfcstfn(i)
	
	  end do

	  if(iwfcstfn(1).eq.0.or.(nwfcst.gt.1.and.iwfcstfn(2).eq.0)) then
	    write(6,*) '*************************************************'
	    write(6,*) 'The first two wind forecast files MUST be present'
	    write(6,*) 'The required files are: '
	    write(6,*) '               ',wfcstfn(1),', ',wfcstfn(2)
	    write(6,*) '*************************************************'
c	    pause
	    stop
	  end if

	  if(iprod.eq.0) then
	    write(6,*) '*************************************************'
	    write(6,*) 'WARNING: The following forecast files are missing'
	    do i=1,nwfcst
	      if(iwfcstfn(i).eq.0) write(6,*) wfcstfn(i)
	    end do
	    write(6,*) 'The simulation will extrapolate from the last'
	    write(6,*) 'available file and so will lose reliability.'
	    write(6,*) '*************************************************'
	  end if
      
      close(40)
      
c----------------------------------------------------------------------
c     assign directory for forecast current data
c     first check if the region is an added region
c----------------------------------------------------------------------


	  a4 = fcstfn(1)(1:4)
      fcstcurdir='TEMP//////OCE/'			! str

c----------------------------------------------------------------------
c     Construct basic grid parameters and bathymetry.
c     delx, dely = bathymetry grid spacing in metres
c         (set grid spacing delx, at mid-latitude of region)
c     restrict region: assuming
c        slick travels at less than 2.5 kts from spill site.
c     rax, ray = ratios of grid spacing to output pixel size
c----------------------------------------------------------------------
    
        open(50,file='../DTM_INP/dtm.bath',status='old')


      read(50,'(a80)') empty 
      read(50,*) along1,along2,alatg1,alatg2
      read(50,*) mmax,nmax          

      dlong = (along2 - along1) / dfloat(mmax-1)
      dlatg = (alatg2 - alatg1) / dfloat(nmax-1)
      avlat = (alatg1 + alatg2) / 2.d0


c - original medslik setup
c          dist = 2.5d0 * tcomp                ! distance in NM
c - adapted by augusto neves
          dist = 6.5d0 * tcomp                ! distance in NM
c - Expected advection of the oil spill (2 NM/h) where tcomp is the time of simulation
          dist = dist / 60.d0                 ! convert to deg of latitude
          dep = dist / dcos(splat/degrad)    ! same dist in deg longitude

c - In case dep is larger than the lagrangean grid (mm is the # of cells in the meridional direction:
          if(dep.gt.dlong * (mm-4) / 2.d0) then
            dep = dlong * (mm-4) / 2.d0
            dist = dep * dcos(splat/degrad)
          endif   
c - Longitude limit = lon0Spill + expected advection

          alon1 = splon - dep
          alon2 = splon + dep

          if(alon1.lt.along1) alon1 = along1
          if(alon2.gt.along2) alon2 = along2
          m1 = int( xgrid(alon1) )
	    m2 = int( xgrid(alon2) + 1.d0 )
          mmax1 = m2 - m1 + 1
	    alon1 = glon(dfloat(m1))
        

          alat1 = splat - dist
          alat2 = splat + dist

          if(alat1.lt.alatg1) alat1 = alatg1
          if(alat2.gt.alatg2) alat2 = alatg2
          n1 = int( ygrid(alat1) )
	    n2 = int( ygrid(alat2) + 1.d0 )
	  
          nmax1 = n2 - n1 + 1
          alat1 = glat(dfloat(n1)) 

c Loading bathymetry
          do n=nmax,1,-1
          do m=1,mmax

          read(50,'(i4)') ih1     
          if(n.ge.n1.and.n.le.n2.and.m.ge.m1.and.m.le.m2) then
            mp = m - m1 + 1 
            np = n - n1 + 1

          if(ih1.ne.9999) then

          itype(mp,np) = ih1

     	  endif

          endif 

	      enddo
          enddo
c


          mmax = mmax1
          nmax = nmax1
  
          along1 = alon1      
 
          along2 = alon1 + (mmax-1) * dlong      
          alatg1 = alat1      
          alatg2 = alat1 + (nmax-1) * dlatg      

 
        

      close(50)


      avlat = (alatg1 + alatg2) / 2.d0
  

      dely = dlatg*60.d0*1849.d0
      delx = dlong*60.d0*1849.d0*dcos(avlat/degrad)

      rax = delx / al5
      ray = dely / al5

      hmax = 0.
      do n=nmax,1,-1
        do m=1,mmax
          if(itype(m,n).eq.9999) itype(m,n) = 0
          h(m,n) = float( itype(m,n) )
          if(itype(m,n).gt.0) itype(m,n) = 1
          if(h(m,n).lt.hmin.and.itype(m,n).ne.0) h(m,n) = hmin
          if(h(m,n).gt.hmax.and.itype(m,n).ne.0) hmax = h(m,n)
        enddo
      enddo



c----------------------------------------------------------------------
c     vertd1,2 = vertical  diffusion distance during time delt
c     horizd = horizontal diffusion distance during time delt
c     check stability of horizontal diffusion simulation (ismag=0)
c     ismag = 1 if Smagorinski scheme used for horizd (forecast data only)
c----------------------------------------------------------------------
      vertd1=dsqrt(6.d0*vertk1*delt*3600.d0)
      vertd2=dsqrt(6.d0*vertk2*delt*3600.d0)

      if(ismag.eq.0) then
        horizd=dsqrt(6.d0*horizk*delt*3600.d0)
        dtmx1=delx**2/(6.d0*horizk*3600.d0)
        dtmx2=dely**2/(6.d0*horizk*3600.d0)
        if(delt.gt.dtmx1.or.delt.gt.dtmx2) then
          write(6,*)'delt must not exceed ',dtmx1,' or ',dtmx2,' hrs'
c	    pause
          stop
        end if
      end if
c----------------------------------------------------------------------
c     x0,y0 = location of spill in coords of the bathymetry grid
c     Rough check that spill is inside the given water body
c     write headings to screen and log file
c----------------------------------------------------------------------
      x0=xgrid(splon)
      y0=ygrid(splat)
	
        m0=int(x0)
	n0=int(y0)
	isum=itype(m0,n0)+itype(m0+1,n0)+itype(m0,n0+1)+itype(m0+1,n0+1)

	if(m0.lt.1.or.m0.gt.mmax-1.or.n0.lt.1.or.n0.gt.nmax-1.or.
     &   isum.eq.0) then
	  write(6,*) '************************************************'
	  write(6,*) ' Spill location is outside the given water body'
	  write(6,*) ' You should continue only if using restart file'
	  write(6,*) '************************************************'
c	  pause
      else if(isum.le.2) then
	  write(6,*) '************************************************'
	  write(6,*) 'WARNING:  Spill location is very close to a coast'
	  write(6,*) '   and  may be outside the given water body.'
	  write(6,*) 'Check its latitude and longitude and if necessary'
	  write(6,*) '      move it further into the water body'
	  write(6,*) '************************************************'
c	  pause
	endif


      write(6,*) 'The spill simulated has the following parameters:'

      write(6,600) idd,imm,iyr,istart
      write(6,601) splat,splon,x0,y0,tspill,icomp 
      if(ispill.eq.0) write(6,603) splrte,vlunit
      if(ismag.eq.0) write(6,604) horizd
      

  600 format(' Date:           ',i2,'/',i2,'/',i4,'  hour: ',i4.4)    
  601 format(' Location:  lat  ',f8.4,',    lon ',f8.4/
     &       ' Grid coords: x  ',f8.4,',      y ',f8.4/
     &       ' Spill duration: ',f7.0,' hours'/
     &       ' Length of run: ',i4,' hours')
  602 format(' Spill rate:    ',f9.2,' ',a4,'/hr')    
  603 format(' Total volume:  ',f9.2,' ',a4)    
  604 format(' Horiz diffusion distance:  ',f12.5)    

       

  770 format(/'Spill positions will be corrected from observations')

      write(90,'(/'' Welcome to the MEDSLIK run module''/)')
      write(90,*) 'The spill simulated has the following parameters:'
      write(90,600) idd,imm,iyr,istart
      write(90,601) splat,splon,x0,y0,tspill,icomp 
      if(ispill.gt.0) write(90,602) splrte,vlunit
      if(ispill.eq.0) write(90,603) splrte,vlunit
      if(ismag.eq.0) write(90,604) horizd
      

	if(icrn.eq.1) then
        write(6,770)
	write(90,770)
	endif

c--------------------------------------------------------------------
c       Initial particle positions from satellite data or contour
c       data (M. De Dominicis & R. Lardner, 2009)
c--------------------------------------------------------------------
 
       if(isat.eq.1) then
	call seedmedslik(ix)
        call readsat(ix,ntot,px,py)
       write(6,*) 'READING SLICK CONTOUR DATA'
	do i=1,ntot
	  px(i) = xgrid(px(i))
	  py(i) = ygrid(py(i))
        enddo

	nst0 = iage / delt   !no of timesteps for aging the oil properties
        endif
c----------------------------------------------------------------------
c     For case when region is whole mediterranean, print grid details 

c     for the selected subregion +/- 'dist' from spill.
c     subroutine setcst selectes coastal segments within the subregion 
c----------------------------------------------------------------------


        write(90,*) ' '
        write(90,*) 'Grid details in selected region'
        write(90,750) mmax,nmax
        write(90,751) along1,along2,alatg1,alatg2
        write(90,*) ' '
        write(90,*) 'Land/water mask:'
        do n=nmax,1,-1
          write(90,'(400i1)') (itype(m,n),m=1,mmax)
        enddo
        write(90,*) ' '
        write(90,*) 'Bathymetry in selected region'
        do n=nmax,1,-1
          write(90,'(400i5)') (int(h(m,n)),m=1,mmax)
        enddo

  750 format(' Dimensions: ',2i5)
  751 format(' Longitude limits:',2f9.5/' Latitude limits: ',2f9.5)

c----------------------------------------------------------------------
c     convert spill volume to barrels
c     rbm3 = cu m per barrel, gpm3 = gals per cu m.
c	ix=random 6-digit seed for random number generator
c----------------------------------------------------------------------
      rbm3=0.158987d0
      gpm3=219.970d0
      if(vlunit.eq.'bbls') volfac=1.d0
      if(vlunit.eq.'cu.m') volfac=1.d0/rbm3
      if(vlunit.eq.'tons') volfac=1.d0/(rbm3*deno)
      if(vlunit.eq.'gals') volfac=1.d0/(rbm3*gpm3)
      splrte=splrte*volfac
c
	call seedmedslik(ix)

c----------------------------------------------------------------------
c     totbbl = total no of barrels spilt
c     nmini = number of incremental mini-spills one mini-spill per delt)
c     vmspl=volume per mini-spill (in cu m)
c
c     ntot = total no of lagrangian parcels spilt (initialy from par-file)
c     bpp = bbls per lagrangian parcel
c     nppms = no of parcels per mini-spill
c     sscale=scale factor for spill printout (= 1 if totbbl < 500000)
c----------------------------------------------------------------------
      if(ispill.gt.0) totbbl=splrte*tspill
      if(ispill.eq.0) totbbl=splrte
c      
      nmini=int(tspill / delt + 0.0001d0) + 1
      if(tspill.gt.(nmini-1)*delt) nmini=nmini+1
      vmspl=totbbl*rbm3/nmini
c
      nppms=int( dfloat(ntot)/dfloat(nmini) + 0.01d0 )
	ntot=nppms*nmini
      bpp=totbbl/ntot
      bblpms=nppms*bpp
c      
      sscale=1.0d0
      i=int(1.d0+dlog10(totbbl/500000.d0))
      if(i.ge.1) sscale=10**i

	write(90,*) 'Total bbls released = ',totbbl
	write(90,*) 'Total no of parcels = ',ntot
	write(90,*) 'Bbls per parcel     = ',bpp
	write(90,*) 'No of minispills    = ',nmini
	write(90,*) 'Parcels per m_spill = ',nppms
c----------------------------------------------------------------------
c     c1i = fraction of evaporative part
c     c2i = fraction of residual part,   c1i + c2i = 1
c     den1 = density of evaporative part (kg/m**3)
c     den2 = density of residual part (kg/m**3)
c     deno = oil density (kg/m**3):    deno = c1i*den1 + c2i*den2
c     tsat = saturation temp for evaporative component
c     fmaxe=max evaporative fraction
c
c     denk = coeff for change of density with evaporation
c     cdt, cvt = coeffs for change of density & viscosity with temperature
c     denw0 = sea water density at temp tk0 (kg/m3)
c     tk0 = reference temperature (deg Kelvin)
c----------------------------------------------------------------------
      deno=deno*1000.d0
      den2=den2*1000.d0
      c2i=respc/100.d0
      c1i=1.d0-c2i
      den1=(deno-c2i*den2)/c1i
      tsat=400.d0-3.d0*api
      fmaxe=c1i
      fmaxd=c2i

      denk=(den2-den1)/deno   !0.18
	cdt=0.008d0
	cvt=5000.d0
      denw0=1026.d0
      denw=1026.d0
      row=(denw0-deno)/deno
	tk0=293.d0
c----------------------------------------------------------------------
c     maxst = no of steps of computation of length delt each
c     ned = no of calls to evap/disp subroutine per delt
c     tedsec = time interval in secs between evap/disp calls
c     tiprs = time interval (hours) for printing output
c     nprs = step number at which output printing begins
c     iprs = no of steps between output printing
c----------------------------------------------------------------------
	maxst = tcomp*stph+0.001d0
      if(irestart.eq.1) maxst = maxst - ihrestart*stph+0.001d0
c
      ned=30
      tedsec=deltsc/ned
c
      tiprs=dfloat(iprs)
      nprs=tiprs*stph+0.001d0
      iprs=tiprs*stph+0.001d0
c----------------------------------------------------------------------
c	initialise wind velocity and wind forecast velocity
c----------------------------------------------------------------------
	do m=1,mm
	do n=1,nm
	  wx(m,n)=0.d0
	  wy(m,n)=0.d0
	  winx(m,n)=0.d0
	  winy(m,n)=0.d0
	end do
	end do               

c----------------------------------------------------------------------
c     if booms deployed, read boom data
c     bmtim(k) = hour of kth boom deployment relative to 0 hrs on nday 0
c     bmx1(k),bmy1(k),bmx2(k),bmy2(k) = 5 km grid coords of boom ends 
c     ibmeff(k) = efficiency (%) of kth boom 
c----------------------------------------------------------------------


c       (M. De Dominicis, booms deployment has not been used and tested in the MEDSLIK-II_1.0)

	nbooms=0


      if(ibooms.eq.1) then
        open(48,file='medslik.bms',status='old')
        read(48,'(a80)') empty
        read(48,'(i2)') nbooms
        read(48,'(a80)') empty
        do k=1,nbooms
          read(48,'(i2,1x,i2,1x,i4,2x,i2,3x,i5,4(3x,i2,2x,f5.2),2x,i3)') 
     &            id,im,iy,ihb,lenbm,latdg1,alatm1,londg1,alonm1,
     &            latdg2,alatm2,londg2,alonm2,ibmeff(k)
          
          nday = jdiff(idd,imm,iyr,id,im,iy)
          bmtim(k)=nday*24.+ihb
          bmlat1(k)=dfloat(latdg1)+alatm1/60.d0
          bmlon1(k)=dfloat(londg1)+alonm1/60.d0
          bmlat2(k)=dfloat(latdg2)+alatm2/60.d0
          bmlon2(k)=dfloat(londg2)+alonm2/60.d0

          bmx1(k)=xgrid(bmlon1(k))
          bmy1(k)=ygrid(bmlat1(k))
          bmx2(k)=xgrid(bmlon2(k))
          bmy2(k)=ygrid(bmlat2(k))
        enddo    
        close(48)
      end if

c----------------------------------------------------------------------
c     If icrn = 1, read spill correction data and construct correction 
c     velocities in m/s to shift centre of slick from the previous
c     computed position (xold,yold) to the observed position (xnew,ynew)
c     For 1st correction, t0avg = mean time of oil release
c----------------------------------------------------------------------
      if(icrn.eq.1) then
        open(46,file='medslik.crn',status='old')
        read(46,'(a80)') empty
        read(46,'(i2)') ncrns
        do k=1,ncrns
          read(46,'(a80)') empty 
          read(46,'(i4)') kkk
          crntim(k)=dfloat(kkk)
          read(46,'(i2)') nvert 
          do j=1,nvert+1
            read(46,'(a80)') empty
          end do
          
          read(46,'(f6.3,3x,f6.3)') alat,alon
          xnew=xgrid(alon)
          ynew=ygrid(alat)
                   
          read(46,'(f6.3,3x,f6.3)') alat,alon
          xold=xgrid(alon)
          yold=ygrid(alat)
                    

          deltax=(xnew-xold)*delx
          deltay=(ynew-yold)*dely
          if(k.eq.1) then
	      t0avg = tspill/2.d0
            if(crntim(1).lt.tspill) t0avg = crntim(1)/2.d0
            deltat = (crntim(1) - t0avg)*3600.d0
          endif
          if(k.ge.2) deltat = (crntim(k) - crntim(k-1))*3600.d0
          if(deltat.ne.0) then
            crnu(k)=deltax/deltat
            crnv(k)=deltay/deltat
            if(k.gt.1) then
              crnu(k)=crnu(k)+crnu(k-1)
              crnv(k)=crnv(k)+crnv(k-1)
            end if
          else
            write(6,*) 'Two spill corrections for the same times'
            write(6,*) 'Times are: ',crntim(k)
c	      pause
            stop
          end if
          read(46,'(a80)') empty 
        end do
        close(46)
      end if  

c----------------------------------------------------------------------
c     ttn=thickness of thin slick (m) (= 10 microns)
c     ttki=initial thickness of thick slick (m) (= 2 cm)
c     afac=initial area ratio of thin to thick slick areas (afac=4-8)
c     compute initial areas, volumes & radii of thick and thin slicks
c----------------------------------------------------------------------
      ttn=0.00001d0
      ttki=0.02d0
      afac=4.0d0
	atk0=vmspl/(ttki+afac*ttn)
	atn0=atk0*afac
	vtk0=atk0*ttki
	vtn0=atn0*ttn
	rtk0=dsqrt(atk0/pi)
	rtn0=dsqrt((atn0+atk0)/pi)
	write(90,*) ' '
	write(90,*) 'initial mini-slick radii = ',rtk0,rtn0
	write(90,*) ' '
c----------------------------------------------------------------------
c     initial assignment of mini-spill properties
c     tre = time from release
c     ttk/atk/vtk = thickness/area/volume of thick slick
c     ttn/atn/vtn = thickness/area/volume of thick slick
c     tt0/at0/vt0 = average thickness/total area/total volume of 2 slicks
c     den = average density of oil in minispill i
c     vis/visem = viscosity of oil/mousse in minispill i
c     vtke/vtne/vte = evaporated volumes in thick/thin/total slick
c     vtkd/vtnd/vtd = dispersed volumes in thick/thin/total slick
c     ftk/ftn = evaporated fractions
c     fw = water fraction in mousse
c     xcl/xcs = volumes of large/small oil droplets dispersed below slick
c     pcte/pctd = percentages evaporated/dispersed
c----------------------------------------------------------------------
      do i=1,nmini
        tre(i)=0.d0
        ttk(i)=ttki
        atk(i)=atk0
        atn(i)=atn0
        ato(i)=atk(i)+atn(i)
        vtk(i)=vtk0
        vtn(i)=vtn0
        vto(i)=vtk(i)+vtn(i)
        tto(i)=vto(i)/ato(i)
        den(i)=deno
        vis(i)=viso
        visem(i)=viso

        vtke(i) =0.d0
        vtne(i) =0.d0
        vte(i)  =0.d0
        vtkd(i) =0.d0
        vtnd(i) =0.d0
        vtd(i)  =0.d0
        ftk(i)  =0.d0
        ftn(i)  =0.d0
        fw(i)   =0.d0
        xcl(i)  =0.d0
        xcs(i)  =0.d0
        pcte(i) =0.d0
        pctd(i) =0.d0
	
	write(90,*) 'Initial minispill properties:'
        write(90,*) tre(i)
        write(90,*) ttk(i),atk(i),vtk(i) 
        write(90,*) ttn   ,atn(i),vtn(i) 
        write(90,*) tto(i),ato(i),vto(i) 
        write(90,*) ftk(i),ftn(i),fw(i) 
        write(90,*) ttk(i),atk(i),vtk(i) 
        write(90,*) vtke(i),vtne(i),vte(i) 
        write(90,*) pcte(i),pctd(i) 

         enddo
	
c----------------------------------------------------------------------
c     initial assignment of status parameter is for each parcel
c        is=0 parcel not released,
c        is=1 in the spreading surface slick 
c        is=2 on surface but not spreading
c        is=3 dispersed into water column
c        is=-nsg stuck on shore segment number nsg
c        is=9 beyond boundary limit
c     ib(i) = k or -k indicates parcel is stuck on kth boom
c     c1p(i) = bbls of evaporative component left in ith parcel
c     c2p(i) = bbls of non-evaporative component left in ith parcel
c     px(i), py(i) = horizontal grid coordinates of ith parcel
c     pz(i) = vertical sigma-coordinate of ith parcel (0/1 on bottom/surface)
c----------------------------------------------------------------------
      if(isat.eq.1) then
      do i=1,ntot
        is(i)=1
        ib(i)=0
        c1p(i)=c1i*bpp
        c2p(i)=(1.d0-c1i)*bpp
c        px(i)=x0
c        py(i)=y0:
         pz(i)=1.d0
      enddo
      elseif(isat.eq.0) then 
      do i=1,ntot
        is(i)=0
        ib(i)=0
        c1p(i)=c1i*bpp
        c2p(i)=(1.d0-c1i)*bpp
        px(i)=x0
        py(i)=y0
        pz(i)=1.d0
      enddo
      endif
c----------------------------------------------------------------------
c    open file for fate parameters and write headings
c----------------------------------------------------------------------
      if(irestart.eq.0) then
        open(99,file='../OUT/medslik.fte')
        write(99,600) idd,imm,iyr,istart
        write(99,601) splat,splon,x0,y0,tspill,icomp 
        if(ispill.gt.0) write(99,602) splrte,vlunit
        if(ispill.eq.0) write(99,603) splrte,vlunit
        write(99,680)
                
  680   format('   time   vol spilt  %evap    %srf   %srftot   %disp',
     &       '  %cstfxd   %csttot   visem1   visem2  visoil1  visoil2',
     &       '  denoil1  denoil2   wfrac1  wfrac2  volratio')
c    &       '   wfrac1  wfrac2  volratio denw')

      elseif(irestart.eq.1) then
        open(99,file='../OUT/medslik.fte',access='APPEND')
	endif
c----------------------------------------------------------------------
c     npcl = number of parcels released so far
c     nspill = number of minispills which have occured
c     nseg = number of coastal segments in vicinity of spill
c     totucrn,totvcrn = total observational corrections to u,v up to current time 
c     utot, vtot, nuv used to compute an average slick velocity for diagnostics
c----------------------------------------------------------------------
      npcl=0
      nspill=0
      nseg=0
      icur=0
      totucrn=0.d0
      totvcrn=0.d0
c
      utot=0.d0
      vtot=0.d0
      nuv=0 
c----------------------------------------------------------------------
c     For hot restart read restart file
c----------------------------------------------------------------------
      if(irestart.eq.1) then
        if(iyr.ge.2000) write(ay,'(i2.2)') iyr-2000
        if(iyr.lt.2000) write(ay,'(i2.2)') iyr-1900
	  write(am,'(i2.2)') imm
	  write(ad,'(i2.2)') idd
	  write(ah,'(i2.2)') istart/100
	  write(a3,'(i3.3)') ihrestart
        open(98,file=ay//am//ad//ah//'_'//a3//'.rso',form='unformatted')

	  read(98) jtime,npcl,nspill,pcevp,pcsrf,pcdsp,pccst,deno,viso

	  if(jtime.ne.ihrestart) then
          write(6,*) 'Restart file time ',jtime,' does not agree ',
     &          'with restart time in input file ',ihrestart
c	    pause
          stop
        endif 

        do i=1,npcl
          read(98) is(i),ib(i),c1p(i),c2p(i),alon,alat,pz(i),seep(i) 
          px(i) = xgrid(alon)
          py(i) = ygrid(alat)
        enddo

        do i=1,nspill
          read(98) den(i),vis(i),visem(i),tre(i),c1ms(i),c2ms(i),
     &      atn(i),atk(i),ato(i),ttk(i),tto(i),
     &      vtn(i),vtk(i),vto(i),xcl(i),xcs(i),
     &      vtne(i),vtke(i),vte(i),vtnd(i),vtkd(i),vtd(i),
     &      ftk(i),ftn(i),fw(i),pcte(i),pctd(i) 
	  enddo
        

	  nsgr=1
    1   continue
          read(98,end=2) ns0(nsgr),xx1,yy1,vcst(nsgr)
          segx(nsgr) = (xx1 - along1) / dlong + 1.d0
          segy(nsgr) = (yy1 - alatg1) / dlatg + 1.d0 
	    nsgr = nsgr + 1
          go to 1 
    2   continue
	  close(98)
        nsgr = nsgr - 1 
        
        write(90,*) 'Impacted coastal segments '
        write(90,*) 'Num impacted segments = ',nsgr
	  write(90,*) '    No Old-num Numpcls   Seep'
        do ns=1,nsgr
          num=0
          do i=1,npcl
            if(is(i).eq.-ns0(ns)) num=num+1
          enddo
          write(90,'(3i7,f13.6)') ns,ns0(ns),num,vcst(ns)
        enddo     
        
        totseep=0.
        do ns=1,nsgr
            totseep =totseep + vcst(ns)
        enddo
        write(90,*) 'total beached = ',totseep
        write(90,*) ' '

        pccst1=100.*(totseep)/(npcl*bpp)
	  write(90,*) 'percentages beached = ',pccst,pccst1
           
	endif

c----------------------------------------------------------------------
c	ADDING NETCDF DIMENSION AND VARIABLES
C----------------------------------------------------------------------
      NPARCEL=ntot
      NTIME=icomp/tiprs

c     Dimensions
      nc_status = nf_def_dim(ncid,'parcel_id', NPARCEL, prcl_dimid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_dim(ncid, 'time', NTIME, time_dimid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      dimids(1) = prcl_dimid
      dimids(2) = time_dimid

C     Variables
      nc_status = nf_def_var(ncid, 'latitude', 5, NDIMS, dimids, 
     +     lat_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'longitude', 5, NDIMS, dimids, 
     +     lon_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'evaporative_volume', 5,  
     +     NDIMS, dimids,evol_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'non_evaporative_volume', 5,  
     +     NDIMS, dimids,nvol_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'water_fraction', 5,  
     +     NDIMS, dimids,wc_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'particle_status', 5, 
     +     NDIMS, dimids,prtt_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid,'time', 5, 1,time_dimid, 
     +     time_dimid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'segment_id', 5, 
     +     NDIMS, dimids,seg_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

c	SPILL VARIABLES
      nc_status = nf_def_var(ncid, 'total_fixed_oil',5,
     +     1, dimids(2),tfxd_varid)
      nc_status = nf_def_var(ncid, 'viscosity_emulsion_1',5,
     +     1, dimids(2),vemi_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'viscosity_emulsion_2',5,
     +     1, dimids(2),vemf_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'viscosity_oil_1',5,1,
     +        dimids(2),visi_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'viscosity_oil_2',5,1,
     +        dimids(2),visf_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'density_emulsion_1',5,1,
     +        dimids(2),deni_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'density_emulsion_2',5,1,
     +        dimids(2),denf_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'water_fraction_1',5,1,
     +        dimids(2),fwi_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'water_fraction_2',5,1,
     +        dimids(2),fwf_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_def_var(ncid, 'volume_ratio',5,1,
     +        dimids(2),volr_varid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

C     Assign units attributes to the netCDF variables.
      nc_status = nf_put_att_text(ncid, lat_varid, UNITS, 
     +	   len(LAT_UNITS), LAT_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, lon_varid, UNITS,  
     +     len(LON_UNITS),LON_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, evol_varid, UNITS, 
     +     len(EVOL_UNITS),EVOL_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, nvol_varid, UNITS,  
     +     len(NVOL_UNITS),NVOL_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, wc_varid, UNITS,
     +     len(WC_UNITS), WC_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, prtt_varid, UNITS,
     +     len(PRTS_UNITS), PRTS_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, time_dimid, UNITS,
     +     len(TIME_UNITS), TIME_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

      nc_status = nf_put_att_text(ncid, tfxd_varid, UNITS,
     +     len(TFXD_UNITS), TFXD_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, vemi_varid, UNITS,
     +     len(VEM_UNITS), VEM_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, vemf_varid, UNITS,
     +     len(VEM_UNITS), VEM_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, visi_varid, UNITS,
     +     len(VIS_UNITS), VIS_UNITS)
      nc_status = nf_put_att_text(ncid, visf_varid, UNITS,
     +     len(VIS_UNITS), VIS_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, deni_varid, UNITS,
     +     len(DEN_UNITS), DEN_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, denf_varid, UNITS,
     +     len(DEN_UNITS), DEN_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, fwi_varid, UNITS,
     +     len(WC_UNITS), WC_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, fwf_varid, UNITS,
     +     len(WC_UNITS), WC_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_text(ncid, volr_varid, UNITS,
     +     len(VOLR_UNITS), VOLR_UNITS)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)


c     Other attributes - oil density and volume of each parcel
      nc_status = nf_put_att_double(ncid, evol_varid, OIL_DENSITY, 
     +     5,1,deno)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_double(ncid, evol_varid, PARCEL_VOL,  
     +     5,1,bpp)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

      nc_status = nf_put_att_double(ncid, nvol_varid, OIL_DENSITY, 
     +     5,1,deno)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_double(ncid, nvol_varid, PARCEL_VOL,  
     +     5,1,bpp)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

      nc_status = nf_put_att_double(ncid, nvol_varid, SP_LON,  
     +     5,1,splon)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_double(ncid, nvol_varid, SP_LAT,  
     +     5,1,splat)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_double(ncid, evol_varid, SP_LON,  
     +     5,1,splon)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status = nf_put_att_double(ncid, evol_varid, SP_LAT,  
     +     5,1,splat)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

C     Compressing variables
      nc_status=nf_def_var_deflate(ncid,lat_varid,0,1,4)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status=nf_def_var_deflate(ncid,lon_varid,0,1,4)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status=nf_def_var_deflate(ncid,evol_varid,0,1,4)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status=nf_def_var_deflate(ncid,nvol_varid,0,1,4)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status=nf_def_var_deflate(ncid,wc_varid,0,1,4)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
      nc_status=nf_def_var_deflate(ncid,prtt_varid,0,1,4)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)
C     End define mode.
      nc_status = nf_enddef(ncid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

cc----------------------------------------------------------------------

c     begin the main computational loop;  nst=time step
c     timehr = time in hours from start of spill
c     tsec = time in seconds
c----------------------------------------------------------------------
      write(6,161)
      if(irestart.eq.1) write(6,162) ihrestart
       
  161 format(//' Program initialization is now complete. '/
     &       ' Spill simulation will now commence.'/)
  162 format(' Run will restart from ',i4,' hours after spill'//)
c
      do 94 nst=1,maxst
       
        timehr=nst*delt
      
        if(irestart.eq.1) timehr = timehr + ihrestart
	  tsec=timehr*3600.d0
c----------------------------------------------------------------------
c     select area that contains all the parcels: (mzb,nzb)-(mzf,nzf) for 
c         selection of coastal segments
c	(xavg,yavg) approximates the grid coords of current centre of slick
c     select area (mzb2,nzb2)-(mzf2,nzf2) for interpolation of forecast
c         data, extending distance = 'dist' nm  from spill --
c ****    assumes spill does not travel faster than 1.5 nm/hr on average
c----------------------------------------------------------------------
        if(npcl.eq.0.and.isat.eq.0) then
	    xmin=x0
          xmax=x0
          ymin=y0
          ymax=y0
	  else

c for areal source of spill (added by M. De Dominicis) 

	  xmin=10000.d0
          xmax=-10000.d0
          ymin=10000.d0
          ymax=-10000.d0
          do 34 i=1,npcl
            if(px(i).lt.xmin) xmin=px(i)
            if(px(i).gt.xmax) xmax=px(i)
            if(py(i).lt.ymin) ymin=py(i)
            if(py(i).gt.ymax) ymax=py(i)
   34     continue
        end if
	  
	  mzb=int(xmin)-2
        mzf=int(xmax+1)+2
        nzb=int(ymin)-2
        nzf=int(ymax+1)+2
        if(mzb.lt.1) mzb=1
        if(nzb.lt.1) nzb=1
        if(mzf.gt.mmax) mzf=mmax
        if(nzf.gt.nmax) nzf=nmax
c
c	  xavg=(xmax+xmin)/2.
c	  yavg=(ymax+ymin)/2.
	  xavg=dfloat(mzb+mzf)/2.d0
	  yavg=dfloat(nzb+nzf)/2.d0

	  dist = 1.5 * (tcomp - timehr)
	  mex = int( dist * 1849. / delx )
        nex = int( dist * 1849. / dely )

	  mzb2 = int(xmin) - mex
        mzf2 = int(xmax+1) + mex
        nzb2 = int(ymin) - nex
        nzf2 = int(ymax+1) + nex
        if(mzb2.lt.1) mzb2=1
        if(nzb2.lt.1) nzb2=1
        if(mzf2.gt.mmax) mzf2=mmax
        if(nzf2.gt.nmax) nzf2=nmax

c----------------------------------------------------------------------
c     compute water currents
c     if ismag = 1, compute horiz diffusivity from Smagorinsky model
c         (done later for each minispill)
c----------------------------------------------------------------------
          call fcstcur_1hr(nst,timehr+tstart,delt,nfcst,fcstfn,ifcstfn,
     &                   fcsttim,fcstcurdir)



      m1 = xavg + 0.5d0
      n1 = yavg + 0.5d0
	u7 = usrf(m1,n1)
	v7 = vsrf(m1,n1)
c	write(90,*) m1,n1
c	write(90,*) timehr,u7,v7
c----------------------------------------------------------------------
c     select water current field to be used for slick transport; 
c----------------------------------------------------------------------
        do m=1,mmax
        do n=1,nmax
          if(idepth.eq.0) then
            uadv(m,n)=usrf(m,n)
            vadv(m,n)=vsrf(m,n)
          else if(idepth.eq.1) then
            uadv(m,n)=u10(m,n)
            vadv(m,n)=v10(m,n)
          else if(idepth.eq.2) then
            uadv(m,n)=u30(m,n)
            vadv(m,n)=v30(m,n)
          else if(idepth.eq.3) then
            uadv(m,n)=u120(m,n)
            vadv(m,n)=v120(m,n)
          end if   
        end do
        end do 


c----------------------------------------------------------------------
c     compute wind speed and direction
c	(winx,winy) = wind velocity components (m/s)
c	(wx,wy) = wind velocity that was used to compute forecast 
c         currents. This is already fixed in subroutine fcstcur
c     set evaporation rate parameter that depends on wind speed.
c----------------------------------------------------------------------

	    call ski_ecmwf_wind(xavg,yavg,nst,timehr+tstart,delt,
     &             nwfcst,wfcstfn,iwfcstfn,wfcsttim,winx,winy,wvel,wdir)




        ce1=akew*(wvel*3.6d0)**gamma_evap


c      m1 = x0 + 0.5
c      n1 = y0 + 0.5
c      write(90,*) m1,n1
c      write(90,*) timehr,winx(m1,n1),winy(m1,n1),wvel,wdir

c----------------------------------------------------------------------
c     compute the wind-induced drift velocity for the oil parcels.
c     slick speed = alpha * wind speed at angle beta to right of wind.
c     if ibetared = 1, beta reduces with wind velocity, else const. 
c	subtract a fraction wredfrac of the wind already incorporated in 
c     the forecast current in case when iwindred = 1.
c----------------------------------------------------------------------
	  if(wvel.eq.0.d0) then
          beta=beta0
        else
          beta=beta0*(1.d0-ibetared*wvel/(wvel+halfspeed))
        end if 
	  csbeta=dcos(beta/degrad)
	  snbeta=dsin(beta/degrad)

	  do m=1,mm
	  do n=1,nm
	    winu=winx(m,n) - iwindred*wredfrac*wx(m,n)
	    winv=winy(m,n) - iwindred*wredfrac*wy(m,n)
	    wdrftx(m,n)=alpha*( winu*csbeta+winv*snbeta)
	    wdrfty(m,n)=alpha*(-winu*snbeta+winv*csbeta)
	  end do
	  end do


	 
c----------------------------------------------------------------------
c         Stokes drift velocity calculation
c         using JONSWAP spectrum parameterization 
c         (M. De Dominicis)
c---------------------------------------------------------------------- 
        
c        istoke= 1 !(0: no stoke drift, 1: instantaneous wind)


c--------------------------------------------------------------------
        if(istoke.eq.0) then
        write(6,*) 'ISTOKE= 0 NO STOKE DRIFT CALCULATION'
        elseif(istoke.eq.1) then
        wdirstoke=wdir
        write(6,*) 'STOKE DRIFT CALCULATION'
        endif
	
	
	if(istoke.eq.1.or.istoke.eq.2) then
        call calcfetch(xavg,yavg,wdirstoke,fetch)
	xavg_lon=glon(xavg)
	yavg_lat=glat(yavg)

c       fetch=20000
        grav=9.8
        pi=4.*datan(1.d0)
         
           if(istoke.eq.2) then
           wwdir=wdir_mean(nst)
           elseif(istoke.eq.1) then
	
	   m0=int(xavg+0.5d0)
	   n0=int(yavg+0.5d0)
	   wxx=winx(m0,n0)
	   wyy=winy(m0,n0)
	   
	   wvel=dsqrt(wxx*wxx+wyy*wyy)
	   wwdir=0.d0
	   if(wxx.eq.0.) then
	    wwdir=0.d0
	    if(wyy.gt.0) wwdir=180.d0
	    if(wyy.le.0) wwdir=0.d0
	   else
	    wwdir=datan(wyy/wxx) * degrad
	    if(wxx.lt.0.d0) wwdir=wwdir+180.d0
	    wwdir=270.d0-wwdir
	   endif
	  
           endif 
	 
	        
         do m=1,mm
	 do n=1,nm
          if(istoke.eq.2) then
          ww=wvel_mean(nst)
          elseif(istoke.eq.1) then
          ww=dsqrt((winx(m,n)*winx(m,n))+(winy(m,n)*winy(m,n))) !wind velocity module
          endif         
         
	 if(ww.eq.0)then
         stoke=0
         stoku(m,n)=0
         stokv(m,n)=0
         else
	  

	 gamma=3.3
	 alfa=0.076*((ww**2)/(fetch*grav))**(0.22)
	 fang_peak=22*((grav**2)/(ww*fetch))**(0.3333333333)
	 freq_sp(1)=0.001
	  do k=2,700
	  freq_sp(k)=freq_sp(k-1)+0.001
	  fang(k)=2*pi*freq_sp(k)
	    if (fang(k).ge.fang_peak) sigma=0.09
	    if (fang(k).lt.fang_peak) sigma=0.07
          erre(k)=exp(-((fang(k)-fang_peak)**2)
     &/(2*(sigma**2)*(fang_peak**2)))
          sp_exp2(k)=exp(-1.25*(fang_peak/fang(k))**4)
          sp_exp1(k)=alfa*grav**2*fang(k)**(-5)
	  spectra(k)=sp_exp1(k)*sp_exp2(k)*(gamma**erre(k))
	  wave_num(k)=fang(k)*fang(k)/grav
	  stoke_sp(k)=2*spectra(k)*fang(k)*wave_num(k)
	  enddo
	 hwave_tot=0
	 stoke_tot=0
          do k=1,699
	  hwave_d(k)=0.001*2*pi*(spectra(k)+spectra(k+1))/2
	  stoke_d(k)=0.001*2*pi*(stoke_sp(k)+stoke_sp(k+1))/2
	  stoke_tot=stoke_tot+stoke_d(k)
	  hwave_tot=hwave_tot+hwave_d(k)
	  enddo
	 hwave=4*sqrt(hwave_tot)
	 stoke=stoke_tot

	          
          

            wwangle = (270.d0-wwdir) / degrad
            cswwdir = dcos(wwangle)
	    snwwdir = dsin(wwangle)

          stoku(m,n)=stoke*cswwdir 
	  stokv(m,n)=stoke*snwwdir

         endif
       enddo
       enddo

       endif
       
c----------------------------------------------------------------------
c     compute sea surface temperature
c     tc, tk = sea surface temp in deg centigrade and kelvin
c     if isst = 0 tc is obtaines from climatology
c     if isst = 1 tc = sst(m,n) is obtained in subroutine fcstcur
c     if isst = 8 tc is already read from file medslik.inp
c
c     set evap and disp properties that depend on sst
c	denw = density of sea water at temp tk
c     viso = oil viscosity at temp tk; tvk0 = ref temp for initial viscosity
c     re-initialize minispill oil properties that depend on temp
c     write initial data into fate output file for time 0
c----------------------------------------------------------------------
	  m1=int(xavg+0.5d0)
	  n1=int(yavg+0.5d0)
	  tc=sst(m1,n1)

        
	  tk=tc+273.d0
        if(nst.eq.1.and.irestart.eq.0) then
          fac = dexp(cvt*((1.d0/tk)-(1.d0/tvk0)))
          viso = viso * fac
c          write(90,*) viso,fac,cvt,tk,tvk0

            do i=1,nmini
              den(i)=deno
              vis(i)=viso
              visem(i)=viso
	      enddo

	      timout=0.
	      nbblr=0.
	      if(ispill.eq.0) nbblr = totbbl
	      pcevap   =0.d0
	      pcsrf    =100.d0
	      pcsrtot  =100.d0
	      pcdsp    =0.d0
	      pccstfxd =0.d0
	      pccsttot =0.d0
	      wfr      =0.d0
	      volratio =1.d0

            write(99,'(f9.2,f9.0,6f9.4,6f9.2,2f8.3,f10.5)') timehr,
     &        nbblr/volfac,pcevp,pcsrf,pcsrftot,pcdsp,pccstfxd,pccsttot,
     &        visem(1),viso,viso,viso,deno,deno,wfr,wfr,volratio
	  endif

c	  denw=denw0*(1.-cdt*(tk-tk0))
c	  m0=x0+0.5
c	  n0=y0+0.5
c	  write(90,*) tc,wvel,wdir,winx(m0,n0),winy(m0,n0)
c	  if((nst/8)*8.eq.nst) write(90,978) timehr,tk,denw
c  978   format('time = ',f6.1,' hrs    water temp = ',f5.1,' deg K'/
c     &         '           water dens = ',f6.1)
c----------------------------------------------------------------------
c     compute corrections to oil drift velocity from slick observation
c----------------------------------------------------------------------
        if(icrn.eq.1) then
          if(timehr.le.crntim(1)) k1=1
          do k=2,ncrns
            if(timehr.le.crntim(k).and.timehr.gt.crntim(k-1)) k1=k
          end do
          if(timehr.gt.crntim(ncrns)) k1=ncrns
          totucrn=crnu(k1)
          totvcrn=crnv(k1)
        end if

c----------------------------------------------------------------------
c     subroutine coast reads coast segments and coastal types and 
c     selects those segments that lie within the spill region
c     
c     alngt(i) = length of coastal segment i in metres
c     prel(i) = prob of oil being released in interval delt
c     sfrac(i) = frac of oil seeping into segment i in delt
c
c     For restart locate beached parcels on their new coastal segments
c----------------------------------------------------------------------

        call coast(delt,seg,sfrac,prel,nseg,alngt,api,apicoeff)
        
        if(irestart.eq.1.and.nst.eq.1) then
          write(90,*) 'Reassign beached parcels to ',nseg,' new segmnts'
          do ns=1,nss
            seep(ns) = 0.d0
          enddo
          do i=1,npcl
            itmp(i) = 0
          enddo
          
          do n=1,nsgr
            dmin=100000.d0
            do ns=1,nseg
              dx = segx(n) - seg(ns,1)
              dy = segy(n) - seg(ns,2)
              d = dx*dx + dy*dy
              if(d.lt.dmin) then
                dmin = d
                nsmin = ns
              endif
            enddo
	      seep(nsmin) = vcst(n)
c	write(90,'(i5,4d18.8)')n,segx(n),segy(n),seg(nsmin,1),seg(nsmin,2)
c 
            num=0
            do i=1,npcl
              if(is(i).eq.-ns0(n)) then
                itmp(i) = -nsmin
                num=num+1
              endif   
	      enddo
c	      write(90,*) 'Num parcels: ',num,vcst(n)
          enddo
           
          do i=1,npcl
            if(is(i).lt.0) is(i) = itmp(i)
	    enddo

          write(90,*)
          write(90,*) 'New1y impacted coastal segments '
          write(90,*) 'Segment Numpcls Seep'

          totseep = 0.d0
          do ns=1,nss
            if(seep(ns).gt.0.d0) then
              totseep =totseep + seep(ns)
              num=0
              do i=1,npcl
                if(is(i).eq.-ns) num=num+1
              enddo
              write(90,'(2i7,f13.6)') ns,num,seep(ns)
            endif
          enddo     
          write(90,*) 'total beached 1 = ',totseep
          write(90,*) ' '
        endif
c----------------------------------------------------------------------
c     release new mini-spill
c     initially locate new parcels at spill site
c     initialize totals
c----------------------------------------------------------------------
        npcl0=npcl
        if (nspill .lt. nmini) then
          nspill=nspill+1
          npcl=npcl+nppms
          if(isat.eq.1) then
          do 54 i=npcl0+1,npcl
            is(i)=1
c            px(i)=x0
c            py(i)=y0
   54     continue
          elseif(isat.eq.0)then
          do 99 i=npcl0+1,npcl
            is(i)=1
            px(i)=x0
            py(i)=y0
   99     continue
         endif  
       end if
        en1ps=0.d0
        en2ps=0.d0
        en3ps=0.d0
        en4ps=0.d0
        en5ps=0.d0
	  volsrf=0.d0
	  tvolsrf=0.d0
c----------------------------------------------------------------------
c     start mini spill loop
c     mini-spill ns contains parcels l1 through l2
c     tre(ns): time in hours since release of mini spill ns
c     nsps: no of parcels still spreading from mini spill ns
c     xcm, ycm: centre of spreading part of slick mini spill ns
c     xcmall, ycmall: centre of whole mini spill ns
c----------------------------------------------------------------------
	  if(nst.eq.1.and.irestart.eq.0) then
          write(90,*) 
	    write(90,*) 'Release of parcels:'
	  endif
        if(timehr.le.tspill+delt) write(90,*) 'Time hrs = ',timehr

        do 70 ns=1,nspill
          l1=(ns-1)*nppms+1
          l2=ns*nppms
          tre(ns)=tre(ns)+delt
          nsps=0
          c1ms(ns)=0.d0
          c2ms(ns)=0.d0
          xcm=0.d0
          ycm=0.d0
          nspsall = 0
          xcmall = 0.d0
          ycmall = 0.d0
          if(timehr.le.tspill+delt) then 
	      write(90,*) '    Minispill no. ',ns
	      write(90,*) '    Parcels ',l1,' to ',l2
	    endif
c----------------------------------------------------------------------
c     stop minispill spreading after time sprdmx
c     compute centres of spreading part and whole mini-slick
c     compute Smagorinsky diffusivity at centre of each mini-slick
c----------------------------------------------------------------------
          do 58 i=l1,l2
            if (tre(ns).ge.sprdmx .and. is(i).eq.1) is(i)=2
            if (is(i) .eq. 1) then
              xcm=xcm+px(i)
              ycm=ycm+py(i)
              nsps=nsps+1
            end if

            xcmall = xcmall + px(i)
            ycmall = ycmall + py(i)
            nspsall = nspsall + 1
   58     continue

          if (nsps .ge. 1) then
            xcm=xcm/nsps
            ycm=ycm/nsps
          else
            xcm = x0
            ycm = y0
          end if

          if (nspsall .ge. 1) then
            xcmall = xcmall / nspsall 
            ycmall = ycmall / nspsall 
          else
            xcmall = xavg
            ycmall = yavg
          end if
c
        if(ismag.eq.1) then
	    call smag(xcmall, ycmall, usrf, vsrf, horizk)

          horizd=sqrt(6.d0*horizk*delt*3600.d0)
c	    if(mod(nst,4).eq.0) 
c     &        write(91,*) 'Smagorinsky Horiz Diffus: ',nst/2,ns,horizk
        end if   
c----------------------------------------------------------------------
c     subroutine ed computes evaporating and dispersing fractions
c	and rate of mousse formation; called ned times per step
c     evaporation, dispersion and emulsification by mckay
c     pctd,probd = percent dispersed and probability of disperion
c     pcte,frace = percent and fraction evaporated
c	fw(ns) = fraction of water in minispill ns
c
c     rratio = ratio of radii of thick slick before and after time step
c----------------------------------------------------------------------
          pctd0=pctd(ns)
	    atkns0 = atk(ns)
          rtkns0 = dsqrt(atk(ns) / pi) 

          do k=1,ned

            call ed(tedsec,vmspl,den(ns),vis(ns),visem(ns),ttk(ns),ttn,
     &         tto(ns),atk(ns),atn(ns),ato(ns),vtk(ns),vtn(ns),vto(ns),
     &         ftk(ns),ftn(ns),ft,fw(ns),
     &         vtke(ns),vtne(ns),vte(ns),vtkd(ns),vtnd(ns),vtd(ns),
     &         xcl(ns),xcs(ns),pcte(ns),pctd(ns))
          enddo
          
c       Aged spill properties R. Lardner and M. De Dominicis (2009)	
            
          if(isat.eq.1.and.nst.eq.nst0) then 
	    write(90,*) 'Aged minispill properties:'
            write(90,*) tre(ns)
            write(90,*) ttk(ns),atk(ns),vtk(ns) 
            write(90,*) ttn   ,atn(ns),vtn(ns) 
            write(90,*) tto(ns),ato(ns),vto(ns) 
            write(90,*) ftk(ns),ftn(ns),fw(ns) 
            write(90,*) ttk(ns),atk(ns),vtk(ns) 
            write(90,*) vtke(ns),vtne(ns),vte(ns) 
            write(90,*) pcte(ns),pctd(ns) 
          endif

	  
          probd=(pctd(ns)-pctd0)/(100.d0-pctd0)
          frace=pcte(ns)/100.d0
          rtkns = dsqrt(atk(ns) / pi) 
	    rratio = rtkns / rtkns0
c          write(90,*) 'Slick radii: ',rtkns0,rtkns
c----------------------------------------------------------------------
c     displace and transform the lagrangian parcels
c     first evaporation and dispersion check
c----------------------------------------------------------------------
          do 66 i=l1,l2
            c1p(i)=(c1i-frace)*bpp
	     if(isat.eq.1.and.nst.le.nst0) go to 65
            if(is(i).eq.9) go to 66
            rrr=randmedslik(ix)
            if((is(i).eq.1.or.is(i).eq.2).and.rrr.lt.probd) is(i)=3           
c----------------------------------------------------------------------
c     compute spreading displacement (thick slick contains most of the oil
c       since mechanical spreading is stopped after sprdmx = 24 hrs)
c     new parcels randomly & uniformly distributed within initial thick
c         slick radius rtk0. (Thin slick volume is initially v small)
c----------------------------------------------------------------------
            xds=0.d0
            yds=0.d0
            if (is(i).eq.1) then
              if(i.le.npcl0) then
                xds = (rratio - 1.d0) * (px(i) -xcm) 
                yds = (rratio - 1.d0) * (py(i) -ycm) 
              else
                rnd=randmedslik(ix)
                radd=rtk0*dsqrt(rnd)
                phi=2.d0*pi*randmedslik(ix)
                xds=radd*dcos(phi)
                yds=radd*dsin(phi)
              end if
            end if    
c----------------------------------------------------------------------
c     compute diffusion displacement: (m,n) = nearest grid pt to parcel
c----------------------------------------------------------------------
            m=int(px(i)+0.5d0)
            n=int(py(i)+0.5d0)


            xdd=(2.d0*randmedslik(ix)-1.d0)*horizd
            ydd=(2.d0*randmedslik(ix)-1.d0)*horizd 

            zdd=0.d0

    
            if(is(i).eq.3) then
              call intrpl0(px(i),py(i),h,hint)
              dep=(1.d0-pz(i))*hint
              zdd=(2.d0*randmedslik(ix)-1.d0)*
     &vertd(dep,vertd1,vertd2,thermocl)
            end if  
c----------------------------------------------------------------------
c     compute advective displacements of surface & dispersed parcels
c     include wind drift & correction term from observation(s) of spill
c----------------------------------------------------------------------
c     Bilinear interpolation for wind & stoke drift velocity (added by 
c     M. De Dominicis) 


c---------------------------------------------------------------------
        
	    if(is(i).ne.3) then

              call intrpl(px(i),py(i),uadv,itype,ui)
              call intrpl(px(i),py(i),vadv,itype,vi)
              call intrpl(px(i),py(i),wdrftx,itype,wui)
              call intrpl(px(i),py(i),wdrfty,itype,wvi)
	      call intrpl(px(i),py(i),stoku,itype,sui)
              call intrpl(px(i),py(i),stokv,itype,svi)
c	      call intrpl(px(i),py(i),hwave_vec,itype,hi)


              if(istoke.eq.0) then
c             ui=ui+wdrftx(m,n)
c             vi=vi+wdrfty(m,n)
	      
	      ui=ui+wui
              vi=vi+wvi
              elseif(istoke.ne.0) then
c             ui=ui+stoku(m,n)+wdrftx(m,n)
c             vi=vi+stokv(m,n)+wdrfty(m,n)
              ui=ui+sui+wui
              vi=vi+svi+wvi

              endif
c******************************************************                
              if(is(i).eq.1.or.is(i).eq.2) then
                utot=utot+ui+totucrn
                vtot=vtot+vi+totvcrn
                nuv=nuv+1
	        endif
c******************************************************                
             
            else if(is(i).eq.3) then
              dep=(1.d0-pz(i))*hint

              if(dep.le.fcstdep1) then
	      
                ui=(usrf(m,n)*(fcstdep1-dep)+u10(m,n)*dep)/10.d0
                vi=(vsrf(m,n)*(fcstdep1-dep)+v10(m,n)*dep)/10.d0
		
              else if(dep.gt.fcstdep1.and.dep.le.fcstdep2) then
	      
                if(hint.ge.fcstdep2) then
                  denom=fcstdep2-fcstdep1
	            fac1=(dep-fcstdep1)
                  fac2=(fcstdep2-dep) 
                  ui=(u30(m,n)*fac1+u10(m,n)*fac2)/denom
                  vi=(v30(m,n)*fac1+v10(m,n)*fac2)/denom
			  else
                  denom=hint-fcstdep1
	            fac1=(dep-fcstdep1)
                  fac2=(hint-dep) 
                  ui=(u30(m,n)*fac1+u10(m,n)*fac2)/denom
                  vi=(v30(m,n)*fac1+v10(m,n)*fac2)/denom
                end if
		
	

              else if(dep.gt.fcstdep2.and.dep.le.fcstdep3 ) then
	      
                if(hint.ge.fcstdep3) then
                  denom=fcstdep3-fcstdep2
	            fac1=(dep-fcstdep2)
                  fac2=(fcstdep3-dep) 
                  ui=(u120(m,n)*fac1+u30(m,n)*fac2)/denom
                  vi=(v120(m,n)*fac1+v30(m,n)*fac2)/denom

	        else if(hint.lt.fcstdep2) then  !----- original

                  ui=u30(m,n)
                  vi=v30(m,n)

                else
                  denom=hint-fcstdep2
	            fac1=(dep-fcstdep2)
                  fac2=(hint-dep) 
                  ui=(u120(m,n)*fac1+u30(m,n)*fac2)/denom
                  vi=(v120(m,n)*fac1+v30(m,n)*fac2)/denom

                end if

              else if(hint.lt.fcstdep3) then !--- original
    
                ui=u120(m,n)
                vi=v120(m,n)

		
              end if

              end if

            xdc=(ui+totucrn)*delt*3600.d0
            ydc=(vi+totvcrn)*delt*3600.d0
c----------------------------------------------------------------------
c     displace parcels on surface and dispersed parcels
c     (m21,n21) = old grid coordinates on bathy grid (nearest grid point)
c	(ppx,ppy) = new grid coordinates of parcel i
c----------------------------------------------------------------------
c            if(iter.eq.1) then
            m21=int(px(i)+0.5d0)
            n21=int(py(i)+0.5d0)
             
            xdispl = (xds+xdc+xdd)
	    ydispl = (yds+ydc+ydd)
            ppx = px(i) + xdispl/delx
            ppy = py(i) + ydispl/dely

c----------------------------------------------------------------------
c     check if beached parcels are released
c     sfrac = fraction of beached oil seeping into sand per step at 'ns'
c     prel = probability of beached oil being washed off at 'ns'
c     seep(nsg) = volume of oil at the impacted coastal segment 'nsg'
c     seepage slows when it approaches carrying capacity seepmx (bbls/km)
c     (xdispl,ydispl) must be to left of (dxseg,dyseg)
c----------------------------------------------------------------------
            if(is(i).lt.0.and.ib(i).eq.0) then
              nsg=-is(i)
              seepge=c2p(i)*sfrac(nsg)

	        fac = seep(nsg) / ( seepmx * alngt(nsg) )
	if(fac.gt.10.d0) write(90,*) nst,i,fac
	
              reducn = dexp( - fac )
              seepge = seepge * reducn

              c2p(i)=c2p(i)-seepge
              seep(nsg)=seep(nsg)+seepge
              probr=prel(nsg)

	      
              dxseg = seg(nsg,3) - seg(nsg,1)
	        dyseg = seg(nsg,4) - seg(nsg,2)
	        cross = dxseg * (ppy - py(i)) - dyseg * (ppx - px(i))

              if((randmedslik(ix).lt.probr).and. cross.gt.0.d0) then
                pz(i)=1.d0
                is(i)=2
              end if
            end if
c----------------------------------------------------------------------
c     check if parcels stuck on booms are released
c----------------------------------------------------------------------
c            if((is(i).eq.1.or.is(i).eq.2).and.ib(i).ne.0) then
c              k=iabs(ib(i))
c              dxbm = bmx2(k) - bmx1(k)            
c              dybm = bmy2(k) - bmy1(k)
c              cross = dxbm * (ppy - py(i)) - dybm * (ppx - px(i))            
c              if(cross*ib(i).gt.0) ib(i)=0
c            end if
c----------------------------------------------------------------------
c     check if surface parcel hits one of the booms 
c     ib(i) = +k indicates parcel is on right side of boom k
c     ib(i) = -k indicates parcel is on left side of boom k
c----------------------------------------------------------------------
            if((is(i).eq.1.or.is(i).eq.2).and.ib(i).eq.0) then
              xi1 = px(i)            
              eta1 = py(i)
              xi2 = ppx
              eta2 = ppy
              dmin = 100.d0
              nbmin=0
              ibmin=0
              do 64 k=1,nbooms
                if(timehr+tstart.lt.bmtim(k)) go to 64
                xx1 = bmx1(k)
                yy1 = bmy1(k)
                xx2 = bmx2(k)
                yy2 = bmy2(k)
                ddel  = (xx2-xx1)*(eta2-eta1)-(yy2-yy1)*(xi2-xi1)
                ddel1 = (eta2-eta1)*(xi1-xx1)-(xi2-xi1)*(eta1-yy1)
                ddel2 = (yy2-yy1)*(xi1-xx1)-(xx2-xx1)*(eta1-yy1)

                if(ddel.eq.0.d0) go to 64
                if( 100*randmedslik(ix).gt.ibmeff(k) ) go to 64
                 
                alam=ddel1/ddel
                alamp=ddel2/ddel
                if(alam.ge.0.d0.and.alam.le.1.d0.and.
     &                       alamp.ge.0.d0.and.alamp.le.1.d0) then
                   xx=xx1+alam*(xx2-xx1)
                   yy=yy1+alam*(yy2-yy1)
                   dd1=dabs(xx-xi1)+dabs(yy-eta1)
                   if(dd1.lt.dmin) then
                     dmin=dd1
                     pppx=xx
                     pppy=yy
                     nbmin=k
                     if(ddel.gt.0.d0) ibmin=1
                     if(ddel.lt.0.d0) ibmin=-1
                   end if  
                end if
   64         continue
c              
              if(nbmin.ne.0) then
                px(i)=pppx
                py(i)=pppy
                ib(i)=ibmin*nbmin
              end if
            end if
c----------------------------------------------------------------------
c     check if parcel hits any coastal segment
c     if parcel 'i' hits coastal segment 'ns', set is(i) = -ns
c----------------------------------------------------------------------
            if(is(i).gt.0.and.ib(i).eq.0) then
              xi1 = px(i)            
              eta1 = py(i)
              xi2 = ppx
              eta2 = ppy
              dmin = 100.d0
              nsmn=0
c            
              do 62 nsg=1,nseg

                xx1 = seg(nsg,1)
                yy1 = seg(nsg,2)
                xx2 = seg(nsg,3)
                yy2 = seg(nsg,4)

                ddel  = (xx2-xx1)*(eta2-eta1)-(yy2-yy1)*(xi2-xi1)
                ddel1 = (eta2-eta1)*(xi1-xx1)-(xi2-xi1)*(eta1-yy1)
                ddel2 = (yy2-yy1)*(xi1-xx1)-(xx2-xx1)*(eta1-yy1)
                if(ddel.eq.0.d0) go to 62
                 
                alam=ddel1/ddel

                alamp=ddel2/ddel
                if(alam.ge.0.d0.and.alam.lt.1.d0.and.
     &                       alamp.ge.0.d0.and.alamp.lt.1.d0) then
                   xx=xx1+alam*(xx2-xx1)
                   yy=yy1+alam*(yy2-yy1)
c                   dd1=dabs(xx-xi1)+dabs(yy-eta1)
                   dd1=dsqrt( (xx-xi1)*(xx-xi1) + (yy-eta1)*(yy-eta1) )
                   if(dd1.lt.dmin) then
                     dmin=dd1
                     pppx=xx
                     pppy=yy
                     nsmn=nsg
                   end if  
                end if
   62         continue
              
			if(nsmn.ne.0) then
                px(i)=pppx
                py(i)=pppy  
                is(i)=-nsmn
              else
                px(i)=ppx
                py(i)=ppy
              end if
            end if   

c----------------------------------------------------------------------
c     vertical diffusion of dispersed parcels
c	ppz = temporary new sigma coordinate of parcel (1 at surface)
c	if ppz > 1 reflect displacement from surface
c     if ppz*caph < 0.2 (< 20 cm from bottom) parcel is sedimented
c----------------------------------------------------------------------
            if(is(i).eq.3) then
              caph=h(m21,n21)
              ppz=1.0d0
             
              if(caph.gt.0.2d0) ppz=pz(i)+zdd/caph
c (michela)   if(caph.gt.0.2d0) ppz=pz(i)+zdd/caph+zdispl/caph
c              if(ppz.lt.0.d0) ppz=abs(ppz)
c              if(ppz.gt.1.d0) ppz=abs(ppz-2.d0*int((ppz+1.d0)/2.d0))
              if(ppz.gt.1.) ppz=2. - ppz
              if(ppz*caph.lt.0.2) then
                ppz = 0.
                is(i) = 4
              endif   
              pz(i)=ppz
              px(i)=ppx
              py(i)=ppy
            end if

c----------------------------------------------------------------------
c     count fractions of light and heavy components left in minispill ns
c----------------------------------------------------------------------
   65       continue
	    c1ms(ns)=c1ms(ns)+c1p(i)
            c2ms(ns)=c2ms(ns)+c2p(i)
            
   66     continue
          c1ms(ns)=c1ms(ns)/bblpms
          c2ms(ns)=c2ms(ns)/bblpms
c----------------------------------------------------------------------
c     compute totals for output:
c         en1ps = vol of oil still spreading 
c         en2ps = vol of oil on surface but no longer spreading 
c         en3ps = vol of oil dispersed 
c         en4ps = vol of oil on the coast but not permanently there
c         volsrf = vol of oil-water mousse = (vol of oil)/(1-fw) 
c         tvolsrf = total vol of oil-water mousse incl releasable oil on coast
c     Compute spill centre and shift to nearby al5-grid point
c----------------------------------------------------------------------
          do i=l1,l2
            if(is(i).eq.1) en1ps=en1ps+c1p(i)+c2p(i)
            if(is(i).eq.2) en2ps=en2ps+c1p(i)+c2p(i)
            if(is(i).eq.3) en3ps=en3ps+c1p(i)+c2p(i)
            if(is(i).lt.0) en4ps=en4ps+c1p(i)+c2p(i)
            if(is(i).eq.1.or.is(i).eq.2) 
     &                volsrf=volsrf+(c1p(i)+c2p(i))/(1.d0-fw(ns))
            if(is(i).ne.0.and.is(i).le.2) 
     &                tvolsrf=tvolsrf+(c1p(i)+c2p(i))/(1.d0-fw(ns))
          enddo
c		volsrf  =volsrf+(en1ps+en2ps)/(1.d0-fw(ns))
c		tvolsrf =tvolsrf+(en1ps+en2ps+en4ps)/(1.d0-fw(ns))
   70   continue


        volratio=tvolsrf/(npcl*bpp)
        sumx=0.d0
	  sumy=0.d0
	  do i=1,npcl
          sumx=sumx+px(i)
          sumy=sumy+py(i)
        end do

	  xavg=sumx/npcl
	  yavg=sumy/npcl

        mavg = int( (xavg-1.d0) * rax )
        xavg = 1.d0 + mavg / rax
        navg = int( (yavg-1.d0) * ray )
        yavg = 1.d0 + navg /ray

c----------------------------------------------------------------------
c     print spill output. first summary statistics of fate percentages
c     ttseep = vol of oil permanently on the coast -> pccstfxd = %
c     pcsrftot & tvolsrf include oil on coast but not permanently there
c     pcsrf &  volsrf exclude such oil
c     pccsttot includes all oil on coast both permanent and temporary
c     pccstfxd incudes only oil permanently attached to coast
c----------------------------------------------------------------------
        ttseep=0.d0
        do 72 nsg=1,nseg
          ttseep=ttseep+seep(nsg)
 
   72   continue
c

        c1tot=0.d0
        c2tot=0.d0
        do 74 ns=1,nspill
          c1tot=c1tot+c1ms(ns)
          c2tot=c2tot+c2ms(ns)
   74   continue
        c1tot=c1tot/dfloat(nspill)
        c2tot=c2tot/dfloat(nspill)
c
        pcevp    =100.d0*(c1i-c1tot)
        pcsrftot =100.d0*(en1ps+en2ps+en4ps)/(npcl*bpp)
        pcsrf    =100.d0*(en1ps+en2ps)/(npcl*bpp)
        pcdsp    =100.d0*en3ps/(npcl*bpp)
        pccstfxd =100.d0*(ttseep)/(npcl*bpp)
        pccsttot =100.d0*(en4ps+ttseep)/(npcl*bpp)
        nbblr=int(npcl*bpp+0.5d0)


c----------------------------------------------------------------------
c	ADDING VALUES TO NC FILE
C----------------------------------------------------------------------	
C     These settings tell netcdf to write one timestep of data. (The
C     setting of start(2) inside the loop below tells netCDF which
C     timestep to write.)

        if(nst.eq.nprs) then

        ncounter(1) = 1
        ncounter(2) = 1

	ncpointer(2)=nst/nstph        

	do i=1,ntot
         ncpointer(1)=i
         x_coordinate=glon(dble(px(i)))
         y_coordinate=glat(dble(py(i)))
         nc_time=INT(nst/nstph)

c	PARTICLES
         nc_status = nf_put_vara(ncid, lat_varid, ncpointer,
     +        ncounter,y_coordinate)
         if (nc_status .ne. nf_noerr) call handle_err(nc_status)
         nc_status = nf_put_vara(ncid, lon_varid, ncpointer,
     +        ncounter,x_coordinate)
         if (nc_status .ne. nf_noerr) call handle_err(nc_status)
         nc_status = nf_put_vara_double(ncid, evol_varid, ncpointer,
     +        ncounter,c1p(i))
         if (nc_status .ne. nf_noerr) call handle_err(nc_status)
         nc_status = nf_put_vara_double(ncid, nvol_varid, ncpointer,
     +        ncounter,c2p(i))
         if (nc_status .ne. nf_noerr) call handle_err(nc_status)
         nc_status = nf_put_vara_double(ncid, wc_varid, ncpointer,
     +        ncounter,fw(ns))
         if (nc_status .ne. nf_noerr) call handle_err(nc_status)
c	BEACHED VOLUME         
    	 if(is(i).ge.0) then
         nc_status = nf_put_vara_double(ncid, prtt_varid, ncpointer,
     +        ncounter,dble(is(i)))
	 else if(is(i).lt.0) then
         nc_status = nf_put_vara_double(ncid, prtt_varid, ncpointer,
     +        ncounter,-seep(-is(i)))
	 endif
         if (nc_status .ne. nf_noerr) call handle_err(nc_status)
C	TO WHICH BEACH SEGMENT
    	 if(is(i).ge.0) then
         nc_status = nf_put_vara_double(ncid, seg_varid, ncpointer,
     +        ncounter,0)
	 else if(is(i).lt.0) then
         nc_status = nf_put_vara_double(ncid, seg_varid, ncpointer,
     +        ncounter,dble(is(i)))
	 endif
         if (nc_status .ne. nf_noerr) call handle_err(nc_status)
	end do
c	TIME
        nc_status = nf_put_vara_int(ncid, time_dimid,ncpointer(2),
     +        ncounter(2),nc_time)
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)

c	SPILL PROPERTIES
        nc_status = nf_put_vara_double(ncid,tfxd_varid,
     +        ncpointer(2),ncounter(2),pccstfxd)
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, vemi_varid,ncpointer(2),
     +        ncounter(2),visem(1))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, vemf_varid,ncpointer(2),
     +        ncounter(2),visem(nmini))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)	
        nc_status = nf_put_vara_double(ncid, visi_varid,ncpointer(2),
     +        ncounter(2),vis(1))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, visf_varid,ncpointer(2),
     +        ncounter(2),vis(nmini))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, deni_varid,ncpointer(2),
     +        ncounter(2),den(1))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, denf_varid,ncpointer(2),
     +        ncounter(2),den(nmini))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, fwi_varid,ncpointer(2),
     +        ncounter(2),fw(1))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, fwf_varid,ncpointer(2),
     +        ncounter(2),fw(nmini))
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)
        nc_status = nf_put_vara_double(ncid, volr_varid,ncpointer(2),
     +        ncounter(2),volratio)
        if (nc_status .ne. nf_noerr) call handle_err(nc_status)

        endif

c---------------------------------------------------------------------
c----------------------------------------------------------------------
c    write fate parameters
c       time  volrel  %evap   %srf  %disp   %cst  visc1   visc2
c----------------------------------------------------------------------
        write(99,'(f9.2,f9.0,6f9.4,6f9.2,2f8.3,f10.5)') timehr,
     &     nbblr/volfac,pcevp,pcsrf,pcsrftot,pcdsp,pccstfxd,pccsttot,
     &     visem(1),visem(nmini),vis(1),vis(nmini),
     &     den(1),den(nmini),fw(1),fw(nmini),volratio
c----------------------------------------------------------------------
c     print spill distributions on coast, surface and dispersed in water
c     first open files and write headings
c----------------------------------------------------------------------
        if(nst.lt.nprs) go to 92

          jtime=timehr+0.001d0

	    write(a4,'(i4.4)') jtime
            fn=pref//a4//a0(1)

c            open(81,file='output/'//fn(1))
c            write(81,801) vlunit  

c            write(81,810) jtime,nbblr/volfac,vlunit,splat,splon,al5,
c     &          pcevp,pcsrf,pcdsp,pccsttot,volsrf/volfac,tvolsrf/volfac        

  801     format(a4,' of oil on coast')
  802     format(a4,' of oil on surface')
  803     format(a4,' of dispersed oil')
  810     format(i4,'     : Hours after start of spill'/
     &           f8.0,' : ',a4,' of oil released so far'/
     &           2f8.3,'  : Lat & Long of spill location'/
     &           f6.1,'    : Pixel size (m) for spill plotting'/'  ',
     &           '%evap   %srf   %disp    %cst  srf vol  Tot srf vol'/ 
     &           4f9.4,2e14.5)
c----------------------------------------------------------------------
c     count up and print barrels permanently stuck on coast
c	use the coastal segments to count & print parcel densities
c     vcst(nsg) = vol of oil per km at coastal segment 'nsg'
c----------------------------------------------------------------------
          kount=0        
          sum=0.d0
          do 76 nsg=1,nseg
            vcst(nsg)=0.d0
            if(seep(nsg).eq.0.d0) go to 76
            sum=sum+seep(nsg)
            dist=alngt(nsg)
            vcst(nsg)=seep(nsg)/dist
            kount=kount+1        
   76     continue
c
c          write(81,839) kount
c          write(81,840) vlunit

c          do 79 nsg=1,nseg
c            if(vcst(nsg).eq.0.) go to 79
c              blon1=glon(seg(nsg,1))
c              blat1=glat(seg(nsg,2))
c              blon2=glon(seg(nsg,3))
c              blat2=glat(seg(nsg,4))

c              write(81,842) blat1,blon1,blat2,blon2,vcst(nsg)/volfac
c   79     continue
c          write(81,843) vlunit,sum/volfac
          
          kbooms=0
          do 80 k=1,nbooms
            if(timehr+tstart.lt.bmtim(k)) go to 80
            kbooms=kbooms+1
   80     continue  
c          write(81,845) kbooms
          do 81 k=1,nbooms
            if(timehr+tstart.lt.bmtim(k)) go to 81
c            write(81,846) bmlat1(k),bmlon1(k),bmlat2(k),bmlon2(k),
c     &                    ibmeff(k)
   81     continue  
c          close(81)
  839     format(i6,'     : Number of data points')
  840     format('   Lat      Lon      Lat      Lon        ',a4,'/km')

  841     format('    Lat      Lon         ',a4,'/km2')
  852     format('    Lat      Lon             m3/km2')
  842     format(4f11.7,e20.5)
  843     format('Total ',a4,' = ',f12.3)
  844     format(2f9.4,e20.5)
  845     format(i2,'       : Number of booms')
  846     format(4f9.4,i5)
  847     format(2f9.4)
  850     format('   Lat      Lon        u        v   ',
     &            2f9.6,'   dlat,dlon')
  851     format(2f12.6,2f9.4)
c----------------------------------------------------------------------
c     j=2: count up and print barrels on surface
c     j=3: count up and print barrels dispersed
c     use the grid with pixel size al5 to count & print
c----------------------------------------------------------------------
          
		nprs=nprs+iprs
c
   92   continue

        jtime=nst/nstph
        if(jtime*nstph.eq.nst) then
	    if(irestart.eq.0) then
            write(6,*) 'Simulation has now completed ',jtime,' hours'
          else 
            write(6,*) 'Simulation has now completed ',jtime,' hours',
     &               ' after restarting'   
          end if 
        end if 
        
c----------------------------------------------------------------------
c     next step
c----------------------------------------------------------------------
   94 continue
      write(6,*) 'Simulation has completed successfully'
C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
      nc_status = nf_close(ncid)
      if (nc_status .ne. nf_noerr) call handle_err(nc_status)

c******************************************************                
      if(nuv.ne.0) then
        write(6,*) 'Average velocity components of surface slick (m/s):'
        write(6,*) 'East ',utot/nuv,'   North ',vtot/nuv
      end if
c******************************************************                
c----------------------------------------------------------------------
c     Write restart file
c----------------------------------------------------------------------
c	if(irestart.eq.0) then
      if(iyr.ge.2000) write(ay,'(i2.2)') iyr-2000
      if(iyr.lt.2000) write(ay,'(i2.2)') iyr-1900
	write(am,'(i2.2)') imm
	write(ad,'(i2.2)') idd
	write(ah,'(i2.2)') istart/100
	write(a3,'(i3.3)') jtime + ihrestart
      open(98,file=ay//am//ad//ah//'_'//a3//'.rso',form='unformatted')

	write(98) jtime+ihrestart,npcl,nspill,pcevp,pcsrf,pcdsp,pccst,
     &          deno,viso

      do i=1,npcl
          vcst(i)=0.d0
          if(is(i).lt.0) vcst(i) = seep(-is(i))
          alon = glon(px(i))
          alat = glat(py(i))
          write(98) is(i),ib(i),c1p(i),c2p(i),alon,alat,pz(i),vcst(i) 
      enddo

      do i=1,nspill
        write(98) den(i),vis(i),visem(i),tre(i),c1ms(i),c2ms(i),
     &    atn(i),atk(i),ato(i),ttk(i),tto(i),
     &    vtn(i),vtk(i),vto(i),xcl(i),xcs(i),
     &    vtne(i),vtke(i),vte(i),vtnd(i),vtkd(i),vtd(i),
     &    ftk(i),ftn(i),fw(i),pcte(i),pctd(i) 
	enddo

	do ns=1,nseg
	  if(seep(ns).gt.0.d0) then
          xx1 = glon(seg(ns,1)) 
          yy1 = glat(seg(ns,2)) 
          write(98) ns,xx1,yy1,seep(ns)
        endif   
	enddo

	close(98)

      end


c**********************************************************************
c	  LIST OF SUBROUTINES
c**********************************************************************
c     OIL FATE routines
c         subroutine setcst
c         subroutine coast
c         function hlseep
c         function hlwash
c         subroutine ed
c
c     DIFFUSIVITY routines
c         function vertd
c         subroutine smag
c
c     SST routines
c         subroutine readsstfcst
c
c     WIND routines

c         subroutine ski_ecmwf_wind
c         subroutine readwind
c
c     WATER CURRENT routines
c         subroutine fcstcur_1hr
c         subroutine readfcst
c         subroutine readfcst_1hr
c
c     UTILITY routines
c         subroutine interpol
c         subroutine intrpl
c         subroutine hsmoot
c         function julday
c         subroutine  date
c         subroutine merparam
c         subroutine ll2mer
c         subroutine mer2ll
c         function randmedslik
c         subroutine seedmedslik
c
c         
c         subroutine readsat  
c         subroutine calcfetch  
c
c**********************************************************************
c     OIL FATE routines
c**********************************************************************

c**********************************************************************
c     subroutine setcst constructs coastal segments from regional map 
c----------------------------------------------------------------------
	subroutine setcst()
	implicit real*8(a-h,o-z)
      parameter(ndim=1000, nss=350000)
c
      dimension bd1(nss,4),ibd(nss),
     &          cstlat(ndim),cstlon(ndim),icst(ndim)

      character empty*80
	logical ex

      common /blk4/ along1,alatg1,along2,alatg2,dlong,dlatg
	common /encr/ aloncorr,iaencr,ibencr,icencr


        open(51,file='../DTM_INP/dtm.map',status='old')

	
      open(52,file='../DTM_INP/dtmcst1.d')

      read(51,*) ncontours

      k=1
	do ni=1,ncontours
        read(51,*) isle

          read(51,*) x1,y1

        do i=2,isle

            read(51,*) x2,y2

          if( (x1.ge.along1.and.x1.le.along2.and.
     &         y1.ge.alatg1.and.y1.le.alatg2).or. 
     &        (x2.ge.along1.and.x2.le.along2.and.
     &         y2.ge.alatg1.and.y2.le.alatg2) ) then
            bd1(k,1) = x1 
            bd1(k,2) = y1 
            bd1(k,3) = x2 
            bd1(k,4) = y2
	
	      k = k + 1
          endif           
	    x1 = x2
	    y1 = y2
        enddo 
	enddo

	nsegt = k-1
	write(90,*) 'No of coastal segments = ',nseg
      close(51)      

      do ns=1,nsegt
	  ibd(ns) = 1
	enddo
c
c     read beach types
c
      ncstty = 0
      inquire(file='data/dtm.cst',EXIST=ex)
	if(ex) then
        write(6,*) 'Setting beach type data'
        write(90,*) 'Setting beach type data'
        open(80,file='../DTM_INP/dtm.cst')
        do k=1,11
          read(80,'(a80)') empty
        end do
        k = 0    
   95   continue
          read(80,*,end=96) clon,clat,ind
          k = k + 1
          cstlon(k) = clon
          cstlat(k) = clat
          icst(k) = ind
          go to 95
   96   continue 
        close(80)
        ncstty = k
c
c     set beach types for each segment
c
        do ns=1,nsegt
	    ibd(ns) = 1
          x2 = (bd1(ns,1) + bd1(ns,3)) / 2.d0
          y2 = (bd1(ns,2) + bd1(ns,4)) / 2.d0
          dmin = 1000.d0
          do k=1,ncstty
            dist = dabs(x2 - cstlon(k)) + dabs(y2 - cstlat(k))
            if(dist.lt.dmin) then
              dmin=dist
              kmin=k
            end if  
          end do 
c          write(90,'(i5,2f8.4,i5,3f8.4,i5)') 
c     &     ns,alon,alat,kmin,dmin,cstlon(kmin),cstlat(kmin),icst(kmin)
          if(dmin.lt.0.01d0) ibd(ns) = icst(kmin)
        end do
      endif
c
c      output results
c
      write(52,'(''no. of segments = '',i6)') nsegt
      do ns=1,nsegt
        write(52,100) (bd1(ns,j),j=1,4),ibd(ns)
      end do
      
c ---------------------------      
c old format -- Augusto Neves
c  100 format(4f11.7,i2) 
c new format -- Augusto Neves  
  100 format(4f12.7,i2) 
c ---------------------------   
      rewind(52)
      
      return
      end

c**********************************************************************
c     subroutine coast reads coastal segment data
c     selects those coastal segments that are close enough to spill
c----------------------------------------------------------------------
      subroutine coast(dt,seg,sfrac,prel,nseg,alngt,api,apicoeff)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX,
     &          ntm=2000,npc=100000,nss=350000,msp=1200)
c	
      save ind,iused,ibd,bd1,nsegt,dkm
c
      integer*2 iused(nss)
      dimension bd1(nss,4),ibd(nss),dkm(nss),itype(mm,nm),
     &          seg(nss,4),prel(nss),sfrac(nss),alngt(nss)

      common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad
	common /blk4/ along1,alatg1,along2,alatg2,dlong,dlatg
	common /blk5/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

      data ind /0/

      xgrid(alon)=(alon - along1) / dlong + 1.d0
      ygrid(alat)=(alat - alatg1) / dlatg + 1.d0
      glon(x)=along1+(x-1.d0)*dlong
      glat(y)=alatg1+(y-1.d0)*dlatg
 
c     coeff for reduction of coastal retention rate for heavy oils

      cseep = 1.d0
      if(api.lt.30.) cseep = 1. + (30.-api) * apicoeff
c
c   ....convert the boundary segments to grid coords....
c   ....ibd(n) = shore type of segment n (e.g. sand, rock, mangrove etc.)
c   ....dkm(n) = length of segment n in km.
c
      if(ind.eq.0) then

        call setcst()

        read(52,'(18x,i6)') nsegt
        do ns=1,nsegt
          read(52,100) (bd1(ns,j),j=1,4),ibd(ns)
          do j=1,3,2
            bd1(ns,j)   = xgrid(bd1(ns,j))
            bd1(ns,j+1) = ygrid(bd1(ns,j+1))
          enddo      

          dx = ( bd1(ns,3) - bd1(ns,1) ) * delx
          dy = ( bd1(ns,4) - bd1(ns,2) ) * dely
          dkm(ns) = dsqrt(dx*dx + dy*dy) / 1000.d0
        enddo
        close(52)
      end if
      ind=1

c --------------------------
c old format -- Augusto Neves
c  100 format(4f11.7,i2)   
c updated format -- Augusto Neves
  100 format(4f12.7,i2) 
c --------------------------

c
c   ....xb - xf, yb - yf: expected limits of spill....
c
      xb=mzb
      xf=mzf
      yb=nzb
      yf=nzf
c
c     select boundary segments that are close enough to spill site....
c     prel(j) = prob of release from coastal segment j in time step dt
c     sfrac(j) = fraction absorbed onto coast in time step dt
c     alngt(j) = length of segment j in km
c
      j=nseg
      do 20 ns=1,nsegt
        if(iused(ns).eq.1) go to 20
        x1=bd1(ns,1)
        y1=bd1(ns,2)

        x2=bd1(ns,3)
        y2=bd1(ns,4)
        if((x1.ge.xb.and.x1.le.xf.and.y1.gt.yb.and.y1.le.yf).or.
     &     (x2.ge.xb.and.x2.le.xf.and.y2.gt.yb.and.y2.le.yf)) then
          j=j+1
          if(j.gt.100000) then
            write(6,*) 'Too many coastal segments are impacted:' 
            write(6,*) '           reduce length of computation'
c	      pause
            stop
          end if

          iused(ns)=1  
          seg(j,1)=x1
          seg(j,2)=y1
          seg(j,3)=x2
          seg(j,4)=y2
c      write(90,'(4d25.20)') glon(seg(j,1)),glat(seg(j,2)),
c     &                      glon(seg(j,3)),glat(seg(j,4))
          ib=ibd(ns)
          tseep=hlseep(ib) * cseep
          twash=hlwash(ib)
	    if(twash.eq.0.) then
            prel(j)=1.0d0
            sfrac(j)=0.0d0
          else   
            prel(j)  = 1.0d0 - 0.5d0**(dt/twash)
            sfrac(j) = 1.0d0 - 0.5d0**(dt/tseep)
          end if
          alngt(j)=dkm(ns)
        end if
   20 continue
      nseg=j  
  199 format(i3,4f10.3,i6,2f10.5)     
c

      return
      end             

c----------------------------------------------------------------------

c   ....half-life for seepage into the beach; beach type = ib
c   ....coastal types ib:                   
c     1   sand beach                                                            
c     2   sand and gravel beach                                                 
c     3   cobble beach                                                          
c     4   rocky shore                                                           
c     5   seawall; concrete wharf etc                                                    
c     6   exposed headland                                                      
c     7   sheltered sand or gravel beach                                        
c     8   sheltered rocky shore                                                 
c     9   sheltered marsh or mud flats                                          
c----------------------------------------------------------------------
      function hlseep(ib)
	implicit real*8(a-h,o-z)
      if(ib.eq.1) hlseep= 24.d0 
      if(ib.eq.2) hlseep= 36.d0
      if(ib.eq.3) hlseep= 48.d0
      if(ib.eq.4) hlseep= 96.d0 
      if(ib.eq.5) hlseep= 96.d0
      if(ib.eq.6) hlseep= 96.d0 
      if(ib.eq.7) hlseep= 24.d0
      if(ib.eq.8) hlseep= 96.d0
      if(ib.eq.9) hlseep= 24.d0

      return
      end

c----------------------------------------------------------------------
c   ....half-life for washing off the beach; beach type = ib
c----------------------------------------------------------------------
      function hlwash(ib)
	implicit real*8(a-h,o-z)
      if(ib.eq.1) hlwash= 24.d0
      if(ib.eq.2) hlwash= 24.d0
      if(ib.eq.3) hlwash= 24.d0
      if(ib.eq.4) hlwash= 18.d0
      if(ib.eq.5) hlwash=  0.d0
      if(ib.eq.6) hlwash=  1.d0 
      if(ib.eq.7) hlwash=120.d0
      if(ib.eq.8) hlwash=120.d0
      if(ib.eq.9) hlwash=120.d0
      return
      end      

c**********************************************************************
c     subroutine ed computes evaporation probability (eprob) and
c           dispersion probability (dprob) and emulsification rate
c     Program is modified from Mackay, Paterson and Trudel, 'Mathematical 
c         Model of Oil Spill Behaviour', Dec 1980
c----------------------------------------------------------------------
      subroutine ed(dtim,vmspl,den,vis,visem,ttk,ttn,tto,
     &              atk,atn,ato,vtk,vtn,vto,ftk,ftn,ft,fw,
     &              vtke,vtne,vte,vtkd,vtnd,vtd,
     &              xcl,xcs,pcte,pctd)

	implicit real*8(a-h,o-z)
      common /evap/ ce1,ce,vappr,fmaxe
      common /disp/ cd1,cd3,cd4,cd5,vl,vs1,um,stk,stn,fmaxd
      common /emul/ cm1,cm2,cm3,visemx
      common /sprd/ cs1,cs2,cs3
      common /phys/ deno,denk,cdt,viso,visk,cvt,tvk0,denw,tk,tk0,wsms
c----------------------------------------------------------------------
c     evaporation; ftk/ftn = fraction of oil evaporated in thick/thin slicks
c----------------------------------------------------------------------
      dvtke=0.d0
      dftk=0.d0
c      if(pcte.ge.fmaxe*100.d0) then
c        pcte=fmaxe*100.d0
c        dvtne=0.0d0
c        go to 10
c      end if
c
      if(vtk.gt.0.d0) then
c        dextk=ce1*atk*vmol*dtim*(1.0-ftk)/(r*tk*vtk)
        dextk=ce1*atk*dtim*(1.d0-ftk)/(tk*vtk)
        vpoil=vappr*exp(-ce*ftk)
        dftk=dextk*vpoil
	  dftkmax=fmaxe-ftk
	  if(dftk.gt.dftkmax/5.d0) dftk = dftkmax / 5.d0
        dvtke=vtk*dftk/(1.d0-ftk)
      end if
      vtke=vtke+dvtke
c
      dvtne=0.d0
      if(ftn.le.fmaxe) dvtne=vtn*(fmaxe-ftn)/(1.d0-ftn)
      ftn=fmaxe
      vtne=vtne+dvtne
      vte=vtke+vtne
      pcte=100.d0*vte/vmspl
      if(pcte.ge.fmaxe*100.d0) pcte=fmaxe*100.d0
c----------------------------------------------------------------------
c     check for dispersion if it exceeds max amount
c	initial/old values
c----------------------------------------------------------------------
   10 continue
      if(pctd.ge.fmaxd*100.d0) go to 30
      dvtkd=0.d0
      dvtnd=0.d0
	xclo=xcl
	xcso=xcs
      if(ttk.le.0.d0) go to 20
c----------------------------------------------------------------------
c     dispersion in thick slick
c     f = volume dispersed per per second per unit vol of slick
c     fb = fraction of small droplets
c     rbl = total vol of large drops dispersed per second
c     rbs = total vol of small drops dispersed per second
c----------------------------------------------------------------------
      f=cd3*(wsms+1.0)**2
      if(ttk.le.0.d0) go to 20
      fb=1.d0/(1.d0+cd4*(visem/10.d0)**(0.5d0)*(ttk/0.001d0)**1.d0*
     &                                           (stk/24.d0))
	rb=f*ttk
	rbt=rb*atk
	rbl=rb*(1.d0-fb)
	rbs=rb*fb

c      write(90,*) ttk,atk,vtk
c      write(90,*) ttn,atn,vtn
c      write(90,*) f,fb
c      write(90,*) rbl,rbs
c----------------------------------------------------------------------
c     cl, xcl = concentration and amount of large drops dispersed
c     cs, xcs = concentration and amount of small drops dispersed
c     ct, xct = total concentration and amount dispersed  (ppm)
c----------------------------------------------------------------------

	cl=rbl/vl
      xcl=cl*um*atk
      cs=2.d0*rbs/(vs1+cd1)
      xcs=cs*um*atk
      ct=cl+cs
      xct=xcl+xcs
c      write(90,*) xcl,xcs,cl,cs
c----------------------------------------------------------------------
c     dxll = amount lost to lower layer per time step dtim
c     dxcl,dxcs = change in amount of large and small drops during dtim
c----------------------------------------------------------------------
	rd=0.5d0*(cd1-vs1)*cs
	rdt=rd*atk
      dxll=rdt*dtim
      dxcl=xcl-xclo
      dxcs=xcs-xcso
c      write(90,*) dxll,dxcl,dxcs
c----------------------------------------------------------------------
c     volume loss from spill: do not include change in large droplets
c	keep old concentrations
c----------------------------------------------------------------------
      dvtkd=dxll+dxcs !+dxcl
      vtkd=vtkd+dvtkd
      xclo=xcl
      xcso=xcs
c----------------------------------------------------------------------
c     dispersion in thin slick
c----------------------------------------------------------------------
   20 continue
      rtn=f*ttn/(1.d0+cd5*stn/24.d0)
	rtnt=rtn*atn
      dvtnd=rtnt*dtim
      vtnd=vtnd+dvtnd
c
c      ratd=rdt+rtnt
	vtd=vtkd+vtnd
      pctd=vtd*100.d0/vmspl
	dxm=dxll+dvtnd
c      write(90,*) dvtkd,vtkd,dvtnd,vtnd
c      write(90,*) vtd,vmspl,pctd
c      write(90,*) ' '
c----------------------------------------------------------------------
c     calculate mousse formation: fw = water fraction
c----------------------------------------------------------------------
   30 continue
      visr=exp(2.5d0*fw/(1.d0-cm1*fw))
      dfw=cm2*(wsms+1.d0)**2*(1.d0-cm3*fw)*dtim
      fw=fw+dfw
	if(fw.gt.1.d0/cm3) fw=1./cm3
      pcw=100.d0*fw
      visem=vis*visr
	if(visem.gt.visemx) visem=visemx
c----------------------------------------------------------------------
c     spreading: use Fay model for areas of both thick and thin slicks
c     calculate new volumes, areas, properties, etc.
c----------------------------------------------------------------------
   40 continue
      datns=0.
      dvtns=0.
      if(vtk.gt.0.) then
        datos=cs1*(ato**0.333d0)*exp(-cs3/(ttk+0.00001d0))*dtim
        dvtns=ttn*datos
        dvtks=-dvtns
        datks=dvtks/ttk + cs2*atk**0.333d0*ttk**1.333d0*dtim
        vtk=vtk-dvtke-dvtkd-dvtns
        atk=atk+datks
        ttk=vtk/atk
      end if
      vtn=vtn-dvtne-dvtnd+dvtns
c----------------------------------------------------------------------
c     transfer thick slick and droplet clouds to thin slick if ttk <= ttn
c----------------------------------------------------------------------
      if(ttk.le.ttn) then
        vtn=vtn+vtk+xcl+xcs
        vtk=0.d0
        atk=0.d0
        ttk=0.d0
        xcl=0.d0
        xcs=0.d0
      end if
      atn=vtn/ttn
c
      vto=vtn+vtk
      ato=atk+atn
      tto=vto/ato
c----------------------------------------------------------------------
c     calculate new compositions of slicks
c----------------------------------------------------------------------
      ftn=fmaxe
      ftk=ftk+dftk
	if(ftk.gt.fmaxe) ftk=fmaxe
      ft=(ftk*vtk+ftn*vtn)/(vtk+vtn)
      ftn=ftn-dvtns*(ftn-ftk)/vtn
c----------------------------------------------------------------------
c     Calculate new oil parameters. Note, denw & den do not contain 
c	temperature expansion effects so only the ratio den/denw is correct.
c	The ratio only is displayed on the output interface.
c     temperature effect on viscosity is included in main program
c----------------------------------------------------------------------
      if(ttk.gt.0.d0) then

c	  den=denw*fw+deno*(1.d0-fw)*(1.d0-cdt*(tk-tk0))*(1.d0+denk*ftk)
	  den=denw*fw+deno*(1.d0-fw)*(1.d0+denk*ftk)
c        vis=viso*exp(visk*ftk+cvt*((1.d0/tk)-(1.d0/tvk0)))
        fac = exp(visk * ftk)
        vis = viso * fac
c        write(90,*) ftk,fac,viso,vis,visr,visem

      end if
c
      return
      end




c**********************************************************************
c     DIFFUSIVITY routines
c**********************************************************************

c**********************************************************************
c     computes the mean diffusion step at depth 'dep'
c----------------------------------------------------------------------
      function vertd(dep,vertd1,vertd2,thermocl)
	implicit real*8(a-h,o-z)
c
      if(dep.lt.thermocl) then
        vertd=vertd1
      else
        vertd=vertd2
      end if

      return
      end             

c**********************************************************************
c     compute horizontal diffusity from Smagorinsky model
c----------------------------------------------------------------------
      subroutine smag(x0,y0,us,vs,horizk)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX)
      dimension itype(mm,nm),us(mm,nm),vs(mm,nm)
      common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad

      data c /0.1d0/

	m=int(x0+0.5d0)
	n=int(y0+0.5d0)

	if(itype(m,n).eq.0) then
	  ux=0.d0
	  vx=0.d0
	else if(itype(m+1,n).ne.0.and.itype(m-1,n).ne.0) then
	  ux=(us(m+1,n)-us(m-1,n))/(2.d0*delx)
	  vx=(vs(m+1,n)-vs(m-1,n))/(2.d0*delx)
      else if(itype(m+1,n).ne.0.and.itype(m-1,n).eq.0) then
	  ux=(us(m+1,n)-us(m,n))/(delx)
	  vx=(vs(m+1,n)-vs(m,n))/(delx)
      else if(itype(m+1,n).eq.0.and.itype(m-1,n).ne.0) then
	  ux=(us(m,n)-us(m-1,n))/(delx)
	  vx=(vs(m,n)-vs(m-1,n))/(delx)
      else 
	  ux=0.d0 
	  vx=0.d0
	end if

	if(itype(m,n).eq.0) then
	  uy=0.d0
	  uy=0.d0
	else if(itype(m,n+1).ne.0.and.itype(m,n-1).ne.0) then
	  uy=(us(m,n+1)-us(m,n-1))/(2.d0*dely)
	  vy=(vs(m,n+1)-vs(m,n-1))/(2.d0*dely)
      else if(itype(m,n+1).ne.0.and.itype(m,n-1).eq.0) then
	  uy=(us(m,n+1)-us(m,n))/(dely)
	  vy=(vs(m,n+1)-vs(m,n))/(dely)
      else if(itype(m,n+1).eq.0.and.itype(m,n-1).ne.0) then
	  uy=(us(m,n)-us(m,n-1))/(dely)
	  vy=(vs(m,n)-vs(m,n-1))/(dely)
      else 
	  uy=0.d0
	  vy=0.d0
	end if

      fac=ux*ux+0.5d0*(vx+uy)*(vx+uy)+vy*vy
	horizk=c*delx*dely*dsqrt(fac)

	if(horizk.lt.0.5d0) horizk=0.5d0
	if(horizk.gt.15.d0) horizk=15.d0

      return
	end


c**********************************************************************
c     SST routines
c**********************************************************************
c     read forecast SST data and interpolate to Medslik grid
c----------------------------------------------------------------------
      subroutine readsstfcst(fcstcurdir,fn,sst)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX)
c
      dimension sst(mm,nm)
	character fcstcurdir*14,fn*16,dummy*80

      open(71,file=fcstcurdir//fn)
      read(71,*) dummy
	read(71,*) dummy
	read(71,*) alon1,alon2,alat1,alat2,mmax,nmax
      read(71,*) ndata
	read(71,*) dummy

	dlon=(alon2-alon1)/dfloat(mmax-1)
	dlat=(alat2-alat1)/dfloat(nmax-1)

	do m=1,mm
	do n=1,nm
	  sst(m,n)=0.d0
	enddo
	enddo

	do i=1,ndata
	  read(71,*) alat,alon,st
	  m=int((alon-alon1)/dlon+1.01d0)
	  n=int((alat-alat1)/dlat+1.01d0)
	  sst(m,n)=st
	end do

	close(71)
c
c	 now interpolate data onto medslik grid
c	
	call interpol(sst,alon1,alat1,dlon,dlat)
c
	return
	end


c**********************************************************************
c     WIND routines
c**********************************************************************                               
c**********************************************************************
c     constructs wind velocity field in case of forecast wind
c     (winx,winy) = forecast wind velocity at grid point (m,n)
c     (dwinx,dwiny) = increment per time step
c     (modified by M. De Dominicis)
c----------------------------------------------------------------------
      subroutine ski_ecmwf_wind(xavg,yavg,nst,time,delt,nwfcst,
     &                  wfcstfn,iwfcstfn,wfcsttim,winx,winy,wvel,wdir)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX,
     &          ntm=2000,npc=100000,nss=350000,npl=2000,msp=1200)
      save isub,readdata,ifile,iwindrec,dwinx,dwiny
c
	dimension winx(mm,nm),winy(mm,nm),dwinx(mm,nm),dwiny(mm,nm)
      dimension iwfcstfn(30),wfcsttim(30,24),temp(24),itype(mm,nm)

	common /spill/  idd,imm,iyr,ispill,tstart,tcomp,x0,y0
	common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad !,ihcst
	common /blk4/ along1,alatg1,along2,alatg2,dlong,dlatg
	common /temp/ m0,n0

      character dirr*14,wfcstfn(30)*14,fn*14
	logical readdata

      glon(x)=along1+(x-1.d0)*dlong
      glat(y)=alatg1+(y-1.d0)*dlatg

	data isub /1/
	if(isub.eq.0) return         ! end of available data files

	m0=int(xavg+0.5d0)           ! centre of slick
	n0=int(yavg+0.5d0)
	slon=glon(dfloat(m0))
	slat=glat(dfloat(n0))


	nrecs=24 
         
      dtwind = 24.d0 / dfloat(nrecs)             
	timew = time - delt / 2.d0             ! timew = time at middle of step

	frac = dtwind / delt          !no of steps between successive wind records
      frac2 = frac + 0.5d0 


	dirr='TEMP//////MET/'							!str

c
c     On 1st step set times of forecast records: wfcsttim(i,k) = time of 
c         the kth record in ith file measured from 0 hrs on spill date
c	Then set the initial wind forecast file and record
c     fac = no of timesteps from 1st wind record to middle of 1st step
c
      if(nst.eq.1) then

	     
	  write(90,*) ' '
	  write(90,*) 'Wind forecast files and times of wind records:'
	  do i=1,nwfcst
          read(wfcstfn(i)(5:6),'(i2)') iy
          read(wfcstfn(i)(7:8),'(i2)') im
          read(wfcstfn(i)(9:10),'(i2)') id

          ndaym1 = jdiff(idd,imm,iyr,id,im,iy+2000)
		do k = 1,nrecs
            wfcsttim(i,k) = ndaym1 * 24.d0 + (k-1) * dtwind
	      temp(nrecs + 1 - k) = 24.d0 - wfcsttim(i,k)
		end do

		if(ihcst.eq.1) then
            do k = 1,nrecs
              wfcsttim(i,k) = temp(k)
		  end do
	    endif
         write(90,'(a14,2x,24f6.0)')wfcstfn(i),(wfcsttim(i,k),k=1,nrecs)
        end do

	  readdata=.true.
	  ifile=1
	  do k=1,nrecs
	    if(timew.ge.wfcsttim(ifile,k)) iwindrec=k
        end do
c	  fac=(tstart + delt/2.-wfcsttim(ifile,iwindrec))/delt
	  fac=(timew - wfcsttim(ifile,iwindrec))/delt

      end if
c
c	read a fresh pair of wind records
c
	if(timew.ge.wfcsttim(ifile,iwindrec)) then
        readdata=.true.
	  write(90,*) timew,wfcsttim(ifile,iwindrec)
	endif

	if(readdata) then

        readdata=.false.
	  fn=wfcstfn(ifile)
	  iw = iwindrec
	  if(ihcst.eq.1) iw = nrecs + 1 - iwindrec
	  write(6,*) 'Reading wind from file '//fn//', record ',iw
	  write(90,*) 'Reading wind from file '//fn//', record ',iw
	  
	  call readwind(dirr,fn,iw,nrecs,winx,winy)
	  

	  if(iwindrec.eq.nrecs) then
		iwindrec=1
		ifile=ifile+1
	  else
		iwindrec=iwindrec+1
	  end if

	  if(iwfcstfn(ifile).eq.0) then
		isub=0
		return
	  end if

        fn=wfcstfn(ifile)
	  iw = iwindrec
	  if(ihcst.eq.1) iw = nrecs + 1 - iwindrec
	  write(6,*) 'Reading wind from file '//fn//', record ',iw
	  write(90,*) 'Reading wind from file '//fn//', record ',iw	

	  call readwind(dirr,fn,iw,nrecs,dwinx,dwiny)

c
c	calculate the increments for each spill re-computation
c	on initial step compute the values for the time of the spill
c     After 1st step, winx, winy are at half timestep behind 'time'
c
	  do m=1,mmax
	  do n=1,nmax
c
	    if(nst.eq.1) then

		  dwinx(m,n) = (dwinx(m,n)-winx(m,n)) / frac
	      dwiny(m,n) = (dwiny(m,n)-winy(m,n)) / frac
	      winx(m,n) = winx(m,n) + dwinx(m,n) * fac
	      winy(m,n) = winy(m,n) + dwiny(m,n) * fac

		elseif(nst.gt.1) then

		  dwinx(m,n) = (dwinx(m,n)-winx(m,n)) / frac2
	      dwiny(m,n) = (dwiny(m,n)-winy(m,n)) / frac2

	    end if
c
	  end do
	  end do
c
      else
c
c	increment the values for each re-computation of the spill
c
	  do m=1,mmax
	  do n=1,nmax
	    winx(m,n)  = winx(m,n)  + dwinx(m,n)
	    winy(m,n)  = winy(m,n)  + dwiny(m,n)
	  end do
	  end do

	end if
c
c	compute wind speed and direction at centre of spill
c

	wxx=winx(m0,n0)
	wyy=winy(m0,n0)
	wvel=dsqrt(wxx*wxx+wyy*wyy)
	wdir=0.d0
	if(wxx.eq.0.d0) then
	  wdir=0.d0
	  if(wyy.gt.0.d0) wdir=180.d0
	  if(wyy.le.0.d0) wdir=0.d0
	else
	  wdir=datan(wyy/wxx) * degrad
	  if(wxx.lt.0.d0) wdir=wdir+180.d0
	  wdir=270.d0-wdir
	end if
c
	return
	end

c**********************************************************************
c     read wind forecast data and interpolate to Medslik grid
c     optional interpolation also over land points
c----------------------------------------------------------------------
      subroutine readwind(dirr,fn,iwindrec,nrecs,wfx,wfy)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX)
c
      dimension wfx(mm,nm),wfy(mm,nm),wfx2(mm,nm),wfy2(mm,nm),
     &          itype(mm,nm),wx1(24),wy1(24)
	character dirr*14,fn*14,dummy*80,acall*2
	
      common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad
	common /blk4/ along1,alatg1,along2,alatg2,dlong,dlatg
      common /temp/ m0,n0

	data icall /0/
	icall = icall + 1
      write(acall,'(i2.2)') icall

      open(72,file=dirr//fn)
      read(72,*) dummy
	read(72,*) dummy
	read(72,*) alon1,alon2,alat1,alat2,mmaxs,nmaxs
      read(72,*) ndat
	read(72,*) dummy

	dlon=(alon2-alon1)/dfloat(mmaxs-1)
	dlat=(alat2-alat1)/dfloat(nmaxs-1)
      
	do i=1,ndat

	  read(72,*) alat,alon,(wx1(k),wy1(k),k=1,nrecs)
	  m=int((alon-alon1)/dlon+1.01d0)
	  n=int((alat-alat1)/dlat+1.01d0)
	  
	  wfx(m,n)=wx1(iwindrec)
	  wfy(m,n)=wy1(iwindrec)
	  
	end do

	close(72)



      nwp = 0
      do 10 m=1,mm

	do 10 n=1,nm
        wfx2(m,n)=0.d0
        wfy2(m,n)=0.d0
c        if(itype(m,n).eq.0) go to 10          ! (*) 
        nwp = nwp + 1

	  x = along1 + dfloat(m-1) * dlong
	  y = alatg1 + dfloat(n-1) * dlatg

	  xdata = (x - alon1) / dlon + 1.d0
	  ydata = (y - alat1) / dlat + 1.d0
	  mdata=int(xdata)
	  ndata=int(ydata)
	  
	  if(mdata.lt.1.or.ndata.lt.1.or.mdata.gt.mmaxs.or.ndata.gt.nmaxs)
     &	   go to 10

        wfx2(m,n)=(wfx(mdata,ndata)*(dfloat(mdata+1)-xdata)+
     &                wfx(mdata+1,ndata)*(xdata-dfloat(mdata)))
     &                                       * (dfloat(ndata+1)-ydata)
     &           +(wfx(mdata,ndata+1)*(dfloat(mdata+1)-xdata)+
     &             wfx(mdata+1,ndata+1)*(xdata-dfloat(mdata))) 
     &                                       * (ydata-dfloat(ndata))
        wfy2(m,n)=(wfy(mdata,ndata)*(dfloat(mdata+1)-xdata)+
     &                wfy(mdata+1,ndata)*(xdata-dfloat(mdata)))
     &                                       * (dfloat(ndata+1)-ydata)

     &           +(wfy(mdata,ndata+1)*(dfloat(mdata+1)-xdata)+
     &             wfy(mdata+1,ndata+1)*(xdata-dfloat(mdata))) 
     &                                       * (ydata-dfloat(ndata))

   10 continue

             
	do m=1,mm
	do n=1,nm
	  wfx(m,n)=wfx2(m,n)
	  wfy(m,n)=wfy2(m,n)
      end do
	end do


    
	return
	end


c**********************************************************************
c     WATER CURRENT routines
c**********************************************************************                  

c**********************************************************************
c     constructs current field in case of HOURLY forecast data
c     time = current time in hours after 0 hrs on date of spill
c     (added by M. De Dominicis)
c----------------------------------------------------------------------
      subroutine fcstcur_1hr(nst,time,delt,nfcst,fcstfn,ifcstfn,fcsttim,
     &                   fcstcurdir)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX,
     &          ntm=2000,npc=100000,nss=350000,msp=1200)
      save ifile,isub
c
      dimension ifcstfn(720), fcsttim(720), itype(mm,nm)
      character fcstcurdir*14, fcstfn(720)*16, fn*16, a2*2
c     character fcstcurdir*9, fcstfn(720)*16, fn*16, a2*2			!TEMPG

	common /spill/  idd,imm,iyr,ispill,tstart,tcomp,x0,y0
      common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad
	common /blk5/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

      common /fcst/ sst(mm,nm),wx(mm,nm),wy(mm,nm),us(mm,nm),
     &      vs(mm,nm),u10(mm,nm),v10(mm,nm),u30(mm,nm),v30(mm,nm),
     &      u120(mm,nm),v120(mm,nm),dsst(mm,nm),dwx(mm,nm),dwy(mm,nm),     
     &      dus(mm,nm),dvs(mm,nm),du10(mm,nm),dv10(mm,nm),   
     &	     du30(mm,nm),dv30(mm,nm),du120(mm,nm),dv120(mm,nm)     
     
C      common /fcst/ sst(mm,nm),us(mm,nm),vs(mm,nm),
C     &      u10(mm,nm),v10(mm,nm),u30(mm,nm),v30(mm,nm),u120(mm,nm),
C     &      v120(mm,nm),dsst(mm,nm),dus(mm,nm),
C     &      dvs(mm,nm),du10(mm,nm),dv10(mm,nm),du30(mm,nm),dv30(mm,nm),
C     &	    du120(mm,nm),dv120(mm,nm)
c     &      ,wx(mm,nm),wy(mm,nm),dwx(mm,nm),dwy(mm,nm)
	common /temp/ m0,n0

	data ifile /1/,   isub /1/

	if(isub.eq.0) return         ! end of available data files
c
c     on 1st step set the times of the forecast files (hrs after 0 hrs 
c	 on spill date) then open the initial forecast file and read initial data
c     On 1st step, fac is used later to interpolate to mid-step
c     For hindcast all times are backwards from 2400 hrs on start date
c
      m0=x0+0.49
      n0=y0+0.49

      if(nst.eq.1) then

        write(90,*) ' '
        write(90,*) 'Frcst files and file times from 0 hrs on spill day'
      

    
	nfcst=nfcst*24
	do i=1,nfcst
	a2=fcstfn(i)(5:6)

        a2=fcstfn(i)(7:8)
        read(a2,'(i2)') im
        a2=fcstfn(i)(9:10)
        read(a2,'(i2)') id
        a2=fcstfn(i)(11:12)
        read(a2,'(i2)') ih
        
	if(dfloat(ih).eq.int(tstart).and.idd.eq.id) then
	ifile=i

	endif
	enddo

	
	do i=ifile,nfcst

          a2=fcstfn(i)(5:6)
          read(a2,'(i2)') iy
          a2=fcstfn(i)(7:8)
          read(a2,'(i2)') im
          a2=fcstfn(i)(9:10)
          read(a2,'(i2)') id
          a2=fcstfn(i)(11:12)
          read(a2,'(i2)') ih
          nday = jdiff(idd,imm,iyr,id,im,iy+2000)


	  fcsttim(i)=nday*24.d0 + dfloat(ih)-0.5 !average

	    write(90,*) i,'   ',fcstfn(i),fcsttim(i)

        enddo
        write(90,*) ' '

        fn=fcstfn(ifile)
	  write(6,*) 'Reading forecast currents from file ',fn
  	  write(90,*) 'Forecast current directory = ',fcstcurdir
	  write(90,*) 'Reading forecast currents from file ',fn
        write(90,*) ' '

     	  call readfcst_1hr(fcstcurdir,fn,sst,us,vs,u10,v10,
     &                                        u30,v30,u120,v120)

c        
	  fac=(time - 0.5d0*delt - fcsttim(ifile)) / delt
c
      end if  
c
c	open the next forecast file and read new data (on 1st step file #2)
c     time has already been advanced to the end of the current step
c
      if(time - 0.5d0*delt .ge. fcsttim(ifile)) then
	    
	  ifile=ifile+1

	  if(ifcstfn(ifile).eq.0) then      ! past the last available file
		isub=0                      
          write(6,*) 'WARNING: Not enough forecast data is available.' 
          write(6,*) 'Forecast data will be kept constant from now on' 
		go to 22
	  end if

        fn=fcstfn(ifile)
	  write(6,*) 'Reading forecast currents from file ',fn
	  write(90,*) 'Reading forecast currents from file ',fn
        write(90,*) ' '

             call readfcst_1hr(fcstcurdir,fn,dsst,dus,dvs,du10,dv10,
     &                                        du30,dv30,du120,dv120)

	dtfcst=fcsttim(ifile)-fcsttim(ifile-1)
        frac = dtfcst / delt          !no of steps between successive fcsts
        frac2 = frac + 0.5d0 

c       write(90,*) dus(m0,n0),dvs(m0,n0),dwx(m0,n0),dwy(m0,n0)
c
c	calculate the increments for each spill re-computation
c     remember after 1st step, sst, wx, etc are at mid-step
c       delt/2 behind the time of the newly-read dsst, dwx, etc
c	initial step compute values for mid-step =  spilltime + delt/2
c
	  
	  do m=1,mmax
	  do n=1,nmax

		if(nst.eq.1) then

		  dsst(m,n)  = (dsst(m,n)  -sst(m,n))  / frac
		  dus(m,n)   = (dus(m,n)   -us(m,n))   / frac
		  dvs(m,n)   = (dvs(m,n)   -vs(m,n))   / frac
		  du10(m,n)  = (du10(m,n)  -u10(m,n))  / frac
		  dv10(m,n)  = (dv10(m,n)  -v10(m,n))  / frac
		  du30(m,n)  = (du30(m,n)  -u30(m,n))  / frac
		  dv30(m,n)  = (dv30(m,n)  -v30(m,n))  / frac
		  du120(m,n) = (du120(m,n) -u120(m,n)) / frac
		  dv120(m,n) = (dv120(m,n) -v120(m,n)) / frac

		  sst(m,n)  = sst(m,n)  + dsst(m,n)  * fac
		  us(m,n)   = us(m,n)   + dus(m,n)   * fac
		  vs(m,n)   = vs(m,n)   + dvs(m,n)   * fac
		  u10(m,n)  = u10(m,n)  + du10(m,n)  * fac
		  v10(m,n)  = v10(m,n)  + dv10(m,n)  * fac
		  u30(m,n)  = u30(m,n)  + du30(m,n)  * fac
		  v30(m,n)  = v30(m,n)  + dv30(m,n)  * fac
		  u120(m,n) = u120(m,n) + du120(m,n) * fac
		  v120(m,n) = v120(m,n) + dv120(m,n) * fac

		elseif(nst.gt.1) then

		  dsst(m,n)  = (dsst(m,n)  -sst(m,n))  / frac2
		  dus(m,n)   = (dus(m,n)   -us(m,n))   / frac2
		  dvs(m,n)   = (dvs(m,n)   -vs(m,n))   / frac2
		  du10(m,n)  = (du10(m,n)  -u10(m,n))  / frac2
		  dv10(m,n)  = (dv10(m,n)  -v10(m,n))  / frac2
		  du30(m,n)  = (du30(m,n)  -u30(m,n))  / frac2
		  dv30(m,n)  = (dv30(m,n)  -v30(m,n))  / frac2
		  du120(m,n) = (du120(m,n) -u120(m,n)) / frac2
		  dv120(m,n) = (dv120(m,n) -v120(m,n)) / frac2

          end if

	  end do
	  end do
c       write(90,*) 'Increments:'
c       write(90,*) dus(m0,n0),dvs(m0,n0),dwx(m0,n0),dwy(m0,n0)
c       write(90,*) 'Updates:'

      end if  
c
c	increment the values for each re-computation of the spill except first
c
	if(nst.ne.1) then
	  do m=1,mmax
	  do n=1,nmax
          sst(m,n)  = sst(m,n)  + dsst(m,n)
c	    wx(m,n)   = wx(m,n)   + dwx(m,n)
c	    wy(m,n)   = wy(m,n)   + dwy(m,n)
	    us(m,n)   = us(m,n)   + dus(m,n)
	    vs(m,n)   = vs(m,n)   + dvs(m,n)
	    u10(m,n)  = u10(m,n)  + du10(m,n)
	    v10(m,n)  = v10(m,n)  + dv10(m,n)
	    u30(m,n)  = u30(m,n)  + du30(m,n)
	    v30(m,n)  = v30(m,n)  + dv30(m,n)
	    u120(m,n) = u120(m,n) + du120(m,n)
	    v120(m,n) = v120(m,n) + dv120(m,n)
	  end do
	  end do
	endif

   22 continue
c     write(90,*) us(m0,n0),vs(m0,n0),wx(m0,n0),wy(m0,n0)
c      
      return
      end
 
c**********************************************************************
c     read forecast hourly current data and interpolate to Medslik grid
c     (added by M. De Dominicis)
c----------------------------------------------------------------------
      subroutine readfcst_1hr(fcstcurdir,fn,sst,us,vs,u10,v10,
     &                     u30,v30,u120,v120)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX)
c
      dimension sst(mm,nm),
     &          us(mm,nm),vs(mm,nm),u10(mm,nm),v10(mm,nm),
     &          u30(mm,nm),v30(mm,nm),u120(mm,nm),v120(mm,nm)
	character fcstcurdir*14,fn*16,dummy*80
	common /temp/ m0,n0

      open(71,file=fcstcurdir//fn)
      read(71,*) dummy
	read(71,*) dummy
	read(71,*) alon1,alon2,alat1,alat2,mmax,nmax
      read(71,*) ndata
	read(71,*) dummy

	dlon=(alon2-alon1)/dfloat(mmax-1)
	dlat=(alat2-alat1)/dfloat(nmax-1)

	do m=1,mm
	do n=1,nm
	  sst(m,n)=0.d0
	  us(m,n)=0.d0
	  vs(m,n)=0.d0
	  u10(m,n)=0.d0
	  v10(m,n)=0.d0
	  u30(m,n)=0.d0
	  v30(m,n)=0.d0
	  u120(m,n)=0.d0
	  v120(m,n)=0.d0
	enddo
	enddo

      do i=1,ndata
	  read(71,*) alat,alon,st,u,v,u1,v1,u3,v3,u4,v4
	  m=int((alon-alon1)/dlon+1.1d0)
	  n=int((alat-alat1)/dlat+1.1d0)
        if(m.gt.mm.or.n.gt.nm.or.m.lt.1.or.n.lt.1) go to 1 
	  sst(m,n)=st
c	  wx(m,n)=wx1
c	  wy(m,n)=wy1
	  us(m,n)=u
	  vs(m,n)=v
	  u10(m,n)=u1
	  v10(m,n)=v1
	  u30(m,n)=u3
	  v30(m,n)=v3
	  u120(m,n)=u4
	  v120(m,n)=v4
    1   continue
 	end do

	close(71)
c
c	 now interpolate data onto medslik grid
c	
	call interpol(sst,alon1,alat1,dlon,dlat)
	call interpol(us,alon1,alat1,dlon,dlat)
	call interpol(vs,alon1,alat1,dlon,dlat)
	call interpol(u10,alon1,alat1,dlon,dlat)
	call interpol(v10,alon1,alat1,dlon,dlat)
	call interpol(u30,alon1,alat1,dlon,dlat)
	call interpol(v30,alon1,alat1,dlon,dlat)
	call interpol(u120,alon1,alat1,dlon,dlat)
	call interpol(v120,alon1,alat1,dlon,dlat)
c

	return
	end


c**********************************************************************
c     UTILITY routines
c**********************************************************************

c**********************************************************************
c     interpolate an array to the medslik grid - dble precision
c----------------------------------------------------------------------
	subroutine interpol(array,alon1,alat1,dlon,dlat)
	implicit real*8 (a-h,o-z)
      parameter(mm=NMAX,nm=NMAX)
c	
	dimension array(mm,nm),array2(mm,nm)
      dimension itype(mm,nm)

      common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad
	common /blk4/ along1,alatg1,along2,alatg2,dlong,dlatg
	common /blk5/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2
	common /temp/ m0,n0
c
c	first extrapolate to the whole grid (i.e. over the land points)

      knt0=0
    1 continue
        knt=0
	  do m=mzb2,mzf2
	  do n=nzb2,nzf2
		array2(m,n)=array(m,n)
          if(array2(m,n).ne.0.d0) go to 5
		knt=knt+1
	
		sum=0.d0
		no=0
		do i=m-1,m+1
		do j=n-1,n+1
		  if(array(i,j).ne.0.d0.and.i.ge.1.and.j.ge.1.and.
     &	                    i.le.mm.and.j.le.nm) then
		    sum=sum+array(i,j)
			no=no+1
		  end if
		end do
		end do
		if(no.gt.0) array2(m,n)=sum/dfloat(no)

    5   continue
        end do
	  end do
        
	  do m=mzb2,mzf2
	  do n=nzb2,nzf2
	    array(m,n)=array2(m,n)
        end do
	  end do

        knt0=knt0+1

	if(knt.gt.0.and.knt0.le.1) go to 1	
c
c	now interpolate from the data grid to the medslik grid 
c

	do m=mzb2,mzf2
	do n=nzb2,nzf2
        array2(m,n)=0.d0
c        if(itype(m,n).eq.0) go to 10
	  x = along1 + dfloat(m-1) * dlong
	  y = alatg1 + dfloat(n-1) * dlatg
	  xdata = (x - alon1) / dlon + 1.d0
	  ydata = (y - alat1) / dlat + 1.d0
	  mdata=int(xdata)
	  ndata=int(ydata)
	  if(mdata.lt.1.or.ndata.lt.1.or.mdata.gt.mm-1.or.ndata.gt.nm-1)
     &	   go to 10

        array2(m,n)=(array(mdata,ndata)*(dfloat(mdata+1)-xdata)+
     &                array(mdata+1,ndata)*(xdata-dfloat(mdata)))
     &                                       * (dfloat(ndata+1)-ydata)
     &           +(array(mdata,ndata+1)*(dfloat(mdata+1)-xdata)+
     &             array(mdata+1,ndata+1)*(xdata-dfloat(mdata))) 
     &                                       * (ydata-dfloat(ndata))


c        m0 = (35.8d0 - dble(along1)) / dble(dlong) + 1.02d0
c        n0 = (34.8d0 - dble(alatg1)) / dble(dlatg) + 1.02d0
c        if(m.eq.m0.and.n.eq.n0) then
c	    x0 = dble(along1) + dfloat(m0-1) * dble(dlong)
c	    y0 = dble(alatg1) + dfloat(n0-1) * dble(dlatg)
c	    write(90,*) 'Interpolation data'
c	    write(90,*) m,n,mdata,ndata
c	    write(90,*) x0,y0,xdata,ydata
c	    write(90,*) array(mdata,ndata),(mdata+1-xdata),(ndata+1-ydata)
c		write(90,*) array(mdata+1,ndata),(xdata-mdata),(ndata+1-ydata)
c	    write(90,*) array(mdata,ndata+1),(mdata+1-xdata),(ydata-ndata)
c		write(90,*) array(mdata+1,ndata+1),(xdata-mdata),(ydata-ndata)
c		write(90,*) array2(m,n)
c	  end if

   10 continue
      end do
	end do
        
	do m=mzb2,mzf2
	do n=nzb2,nzf2
        array(m,n) = array2(m,n)
      enddo
      enddo

      return
      end	         

c**********************************************************************
c     subroutine intrpl0 interpolates values of an array q
c     from grid values to a point (x,y)
c----------------------------------------------------------------------
      subroutine intrpl0(x,y,q,qint)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX)
c
      dimension q(mm,nm),itype(mm,nm)
      m=int(x)
      n=int(y)
      q1=q(m,n)
      q2=q(m+1,n)
      q3=q(m+1,n+1)
      q4=q(m,n+1)
      qint=(q1*(m+1-x)+q2*(x-m))*(n+1-y)+
     &	       (q4*(m+1-x)+q3*(x-m))*(y-n)
      return
      end

c**********************************************************************
c     subroutine intrpl interpolates values of an array q from grid 
c     values to a point (x,y) taking account of land/water mask
c----------------------------------------------------------------------
      subroutine intrpl(x,y,q,itype,qint)
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX,
     &          ntm=2000,npc=100000,nss=350000,msp=1200)
c     &          ntm=2000,npc=100000,nss=350000,npl=2000,msp=1200)
c
      dimension q(mm,nm),itype(mm,nm)
      m=int(x)
      n=int(y)
      it1=itype(m,n)
      it2=itype(m+1,n)
      it3=itype(m+1,n+1)
      it4=itype(m,n+1)
      q1=q(m,n)
      q2=q(m+1,n)
      q3=q(m+1,n+1)
      q4=q(m,n+1)


      if(it1.eq.0.and.it2.ne.0.and.it3.ne.0.and.it4.ne.0) then
        q1=q2+q4-q3
      else if(it1.ne.0.and.it2.eq.0.and.it3.ne.0.and.it4.ne.0) then
        q2=q1+q3-q4
      else if(it1.ne.0.and.it2.ne.0.and.it3.eq.0.and.it4.ne.0) then
        q3=q2+q4-q1
      else if(it1.ne.0.and.it2.ne.0.and.it3.ne.0.and.it4.eq.0) then
        q4=q1+q3-q2

      else if(it1.eq.0.and.it2.eq.0.and.it3.ne.0.and.it4.ne.0) then
        if(itype(m,n+2).ne.0) then
          q1=2.*q4-q(m,n+2)
        else
          q1=q4
        end if
        if(itype(m+1,n+2).ne.0) then  
          q2=2.*q3-q(m+1,n+2)
        else
          q2=q3
        end if
      
      else if(it1.ne.0.and.it2.eq.0.and.it3.eq.0.and.it4.ne.0) then
        if(itype(m-1,n).ne.0) then
          q2=2.*q1-q(m-1,n)
        else
          q2=q1
        end if  
        if(itype(m-1,n+1).ne.0) then

          q3=2.*q4-q(m-1,n+1)
        else
          q3=q4
        end if
      
      else if(it1.ne.0.and.it2.ne.0.and.it3.eq.0.and.it4.eq.0) then
        if(itype(m+1,n-1).ne.0) then
          q3=2.*q2-q(m+1,n-1)
        else
          q3=q2
        end if
        if(itype(m,n-1).ne.0) then
          q4=2.*q1-q(m,n-1)
        else
          q4=q1
        end if
      
      else if(it1.eq.0.and.it2.ne.0.and.it3.ne.0.and.it4.eq.0) then
        if(itype(m+2,n+1).ne.0) then
          q4=2.*q3-q(m+2,n+1)
        else
          q4=q3
        end if        
        if(itype(m+2,n).ne.0) then
          q1=2.*q2-q(m+2,n)
        else
          q1=q2
        end if        
      
      else if(it1.ne.0.and.it2.eq.0.and.it3.ne.0.and.it4.eq.0) then
        if(itype(m-1,n).ne.0) then
          q2=2.*q1-q(m-1,n)
        else
          q2=(q1+q3)/2.d0
        end if        
        if(itype(m+2,n+1).ne.0) then
          q4=2.*q3-q(m+2,n+1)
        else
          q4=(q1+q3)/2.d0
        end if
      
      else if(it1.eq.0.and.it2.ne.0.and.it3.eq.0.and.it4.ne.0) then
        if(itype(m+2,n).ne.0) then
          q1=2.*q2-q(m+2,n)
        else
          q1=(q2+q4)/2.d0
        end if
        if(itype(m-1,n+1).ne.0) then
          q3=2.*q4-q(m-1,n+1)
        else
          q3=(q2+q4)/2.d0
        end if
      
      else if(it1.ne.0.and.it2.eq.0.and.it3.eq.0.and.it4.eq.0) then
        if(itype(m-1,n).ne.0) then
          q2=2.*q1-q(m-1,n)
        else
          q2=q1
        end if
        if(itype(m,n-1).ne.0) then
          q4=2.*q1-q(m,n-1)
        else
          q4=q1
        end if
        q3=q2+q4-q1

      else if(it1.eq.0.and.it2.ne.0.and.it3.eq.0.and.it4.eq.0) then
        if(itype(m+2,n).ne.0) then
          q1=2.*q2-q(m+2,n)
        else
          q1=q2
        end if
        if(itype(m+1,n-1).ne.0) then
          q3=2.*q2-q(m+1,n-1)
        else
          q3=q2
        end if
        q4=q1+q3-q2

      else if(it1.eq.0.and.it2.eq.0.and.it3.ne.0.and.it4.eq.0) then
        if(itype(m+1,n+2).ne.0) then
          q2=2.*q3-q(m+1,n+2)
        else
          q2=q3
        end if
        if(itype(m+2,n+1).ne.0) then
          q4=2.*q3-q(m+2,n+1)
        else
          q4=q3
        end if
        q1=q2+q4-q3

      else if(it1.eq.0.and.it2.eq.0.and.it3.eq.0.and.it4.ne.0) then
        if(itype(m,n+2).ne.0) then
          q1=2.*q4-q(m,n+2)
        else
          q1=q4
        end if
        if(itype(m-1,n+1).ne.0) then
          q3=2.*q4-q(m-1,n+1)
        else
          q3=q4
        end if
        q2=q1+q3-q4
      end if
c
      qint=(q1*(m+1-x)+q2*(x-m))*(n+1-y)+
     &	       (q4*(m+1-x)+q3*(x-m))*(y-n)
      return
      end

c**********************************************************************
c     smooth the bathymetry - slmin = max bottom slope
c----------------------------------------------------------------------
      subroutine hsmoot(h) 
	implicit real*8(a-h,o-z)
      parameter(mm=NMAX,nm=NMAX,
     &          ntm=2000,npc=100000,nss=350000,msp=1200)
c     &          ntm=2000,npc=100000,nss=350000,npl=2000,msp=1200)
      dimension itype(mm,nm),h(mm,nm) 
      common /blk3/ mmax,nmax,delx,dely,itype,pi,degrad
      
      slmin=0.2d0
      
      do 10 k=1,10
c         sweep right
        do 3 n=1,nmax-1
          do 1 m=1,mmax-1
            if(itype(m,n).eq.0.or.itype(m+1,n).eq.0) go to 1
            sl=dabs(h(m+1,n)-h(m,n))/(h(m+1,n)+h(m,n))
            if(sl.lt.slmin) go to 1
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m+1,n))
            sn=-1.d0
            if(h(m+1,n).gt.h(m,n)) sn=1
            h(m+1,n)=h(m+1,n)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    1     continue
C           sweep left
          do 2 m=mmax-1,1,-1
            if(itype(m,n).eq.0.or.itype(m+1,n).eq.0) go to 2
            sl=dabs(h(m+1,n)-h(m,n))/(h(m+1,n)+h(m,n))
            if(sl.lt.slmin) go to 2
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m+1,n))
            sn=-1.d0
            if(h(m+1,n).gt.h(m,n)) sn=1
            h(m+1,n)=h(m+1,n)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    2     continue
    3   continue
   
c         sweep up
        do 6 m=1,mmax-1
          do 4 n=1,nmax-1
            if(itype(m,n).eq.0.or.itype(m,n+1).eq.0) go to 4
            sl=dabs(h(m,n+1)-h(m,n))/(h(m,n+1)+h(m,n))
            if(sl.lt.slmin) go to 4
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m,n+1))
            sn=-1.d0
            if(h(m,n+1).gt.h(m,n)) sn=1
            h(m,n+1)=h(m,n+1)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    4     continue
C           sweep down
          do 5 n=nmax-1,1,-1
            if(itype(m,n).eq.0.or.itype(m,n+1).eq.0) go to 5
            sl=dabs(h(m,n+1)-h(m,n))/(h(m,n+1)+h(m,n))
            if(sl.lt.slmin) go to 5
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m,n+1))
            sn=-1.d0
            if(h(m,n+1).gt.h(m,n)) sn=1
            h(m,n+1)=h(m,n+1)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    5     continue
    6   continue
   
   10 continue

      return
      end

c**********************************************************************
c  function julday computes the julian day for a given date
c----------------------------------------------------------------------
      function julday(idd,imm,iyr)
      dimension js(12),je(12)
c
      ly=0
      if( ( ((iyr/4)*4.eq.iyr) .and. ((iyr/100)*100.ne.iyr) ) .or.
     &    ( (iyr/400)*400.eq.iyr) ) ly=1
      js(1)=1
      je(1)=31
      js(2)=32
      je(2)=59+ly
      k=1
      do 5 i=3,12
        js(i)=je(i-1)+1
        je(i)=js(i)+29+k
        k=1-k
        if(i.eq.7) k=1
    5 continue
    
      julday=js(imm)+idd-1
    
      return
      end     

c**********************************************************************
c  subroutine date computes the date for a given julian day no and year
c----------------------------------------------------------------------
      subroutine  date(jd,iyr,id1,im1)
      dimension js(12),je(12)
c
    1 continue
      ly=0
      if( ( ((iyr/4)*4.eq.iyr) .and. ((iyr/100)*100.ne.iyr) ) .or. 
     &      ((iyr/400)*400.eq.iyr) ) ly=1
      js(1)=1
      je(1)=31
      js(2)=32
      je(2)=59+ly
      k=1
      do 5 i=3,12
        js(i)=je(i-1)+1
        je(i)=js(i)+29+k
        k=1-k
        if (i.eq.7) k=1
    5 continue
c
      iy1=iyr
      if(jd.gt.je(12)) then
        iyr=iyr+1
        jd=jd-je(12)
        go to 1
      else if(jd.lt.0) then
        iyr=iyr-1
        jd=jd+je(12)
	  if( ( ((iyr/4)*4.eq.iyr) .and. ((iyr/100)*100.ne.iyr) ) .or. 
     &      ((iyr/400)*400.eq.iyr) ) jd=jd+1
     	  go to 1
      end if
c          
      do 10 i=1,12
        if(jd.ge.js(i).and.jd.le.je(i)) im1=i
   10 continue
      id1=jd-js(im1)+1
c
      return
      end           


c**********************************************************************
c  function jdiff gives the no of days from id1/im1/iy1 to id2/im2/iy2 
c      the day-difference may be negative and the years differ by 1
c----------------------------------------------------------------------
      function  jdiff(id1,im1,iy1,id2,im2,iy2)

      if(iy2.eq.iy1) jdiff = julday(id2,im2,iy2) - julday(id1,im1,iy1)
	if(iy2.eq.iy1+1) jdiff = julday(id2,im2,iy2) + julday(31,12,iy1)
     &                                  - julday(id1,im1,iy1)
	if(iy2.eq.iy1-1) jdiff = julday(id2,im2,iy2) - julday(31,12,iy2)
     &                                  - julday(id1,im1,iy1)

	return
	end
c**********************************************************************
c    ....parameters for mercator metric coords....
c----------------------------------------------------------------------
      subroutine merparam(avlat)
c
      implicit real*8(a-h,o-z)
      common /data1/ pi,degrad1,fi,ec,ecc,demod,a
      
      pi=3.1415926536d0
      degrad1=pi/180.d0
      fi=avlat*degrad1
      ec=0.082094438d0
      ecc=ec/2.d0
      e1=ec*dsin(fi)
      demod=dsqrt(1.d0-e1*e1)/dcos(fi)
      a=6378137.d0

      return
      end

c**********************************************************************
c    ....converts latitude, longitude to mer coords....
c----------------------------------------------------------------------
      subroutine ll2mer(avlat,alat,alon,x,y,nst,maxst)
      implicit real*8 (a-h,o-z)
      save ind      
c
      common /data1/ pi,degrad1,fi,ec,ecc,demod,a
      
      data ind /0/
      if(ind.eq.0) call merparam(avlat)
      ind=1

      rlat=alat*degrad1
      rlon=alon*degrad1
      x=a*rlon/demod
      ec1=ec*dsin(rlat)
      vy=( (1.d0-ec1)/(1.d0+ec1) ) ** ecc
      vy=vy*dtan(pi/4.d0+rlat/2.d0)
      y=(a/demod)*dlog(vy)
      
      return
      end

c**********************************************************************
c    ....converts mercator metric coords to latitude, longitude....
c----------------------------------------------------------------------
      subroutine mer2ll(avlat,x,y,alat,alon)
      implicit real*8(a-h,o-z)
      save ind
c
      common /data1/ pi,degrad1,fi,ec,ecc,demod,a
      
      data ind /0/
      if(ind.eq.0) call merparam(avlat)
      ind=1

      rlon=x*demod/a
      rlat=0.68d0
    1 continue
        
        rlat1=rlat
        ec1=ec*dsin(rlat1)
        vy=( (1.d0-ec1)/(1.d0+ec1) ) ** ecc
        ve=dexp(y*demod/a) / vy
        rlat=2.d0*datan(ve) - pi/2.d0
        if(dabs((rlat-rlat1)/rlat).gt.0.0000001) go to 1
        
      alat=rlat/degrad1
      alon=rlon/degrad1
      return
      end

c**********************************************************************
c      function rand computes random numbers between 0.0 and 1.0
c----------------------------------------------------------------------
      function randmedslik(ix)
	implicit real*8(a-h,o-z)
    
      integer a,p,ix,b15,b16,xhi,xalo,leftlo,fhi,k

      data a/16807/,b15/32768/,b16/65536/,p/2147483647/
      xhi=ix/b16
      xalo=(ix-xhi*b16)*a
      leftlo=xalo/b16
      fhi=xhi*a+leftlo
      k=fhi/b15
      ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if(ix.lt.0) ix=ix+p
      randmedslik=dfloat(ix)*4.656612875d-10

      return
      end

c**********************************************************************
c      calculates an approximately random first seed for rand
c      (modified by M. De Dominicis)
c----------------------------------------------------------------------
        subroutine seedmedslik(ix)
        real*8 hms
c        character date*10
c        call date_and_time(time=date)
c        read(date,'(f10.3)') hms
	character date*8, time*10
	call date_and_time(date,time)
        call date_and_time(DATE=date)
        call date_and_time(TIME=time)
	read(time,'(f10.3)') hms
        ix=int(hms*100)
        iz=ix-(ix/1000000)*1000000
        iy=1000000
        ix=0
        do k=1,6
          ix0=iz
          iz=iz/10
          ix1=ix0-iz*10
          iy=iy/10
          ix=ix+ix1*iy
        end do
        write(90,*) 'Random seed = ',ix
        return
        end
      include 'codes.h'

c**********************************************************************
c      reads satellite contour
c      (R. Lardner, M. De Dominicis, 2009)
c      (A. Sepp, 2021)
c----------------------------------------------------------------------
	subroutine readsat(ix,ntot,px,py)
      implicit real*8 (a-h,o-z)
      parameter(npc=100000)

	dimension segx(17500,2), segy(17500,2), px(npc), py(npc)
      character datestamp*10, aslik*2, empty*80

	data datestamp /'0711020357'/, aslik /'01'/

      open(38,file='initial.txt')
      read(38,*) empty
c   load number of data points (segments)
      read(38,*) ndata
      read(38,*) empty

c   nsegs is the index for the segment of interest
	  nsegs = 1

c   create a counter
      c = 1

c   load fresh pair of coordinates

      read(38,*) ystart, xstart


c   assing those to variables segx and segx
	  segx(nsegs,1) = xstart
	  segy(nsegs,1) = ystart

c   create 'reference points' 
	  xlast = xstart
	  ylast = ystart

c   if you have less than 50 points in your polygons...
      if(ndata.lt.50) then


       do k=3,ndata
	  read(38,*) y,x
	   
	  if(dabs(x-xlast).lt.1.and.dabs(y-ylast).lt.1
     &.and.xpre.ne.xini.and.ypre.ne.yini) then
	    segx(nsegs,2) = x
	    segy(nsegs,2) = y
          nsegs = nsegs + 1
	    segx(nsegs,1) = x
	    segy(nsegs,1) = y
          xlast = x
          ylast = y
	  xpre=x
	  ypre=y
        else   
	    segx(nsegs,2) = xstart
	    segy(nsegs,2) = ystart
            xini=x
	    yini=y  
	    write(6,*) 'new seg at nsegs, ndata = ',nsegs,k
          nsegs = nsegs + 1
	    xstart = x
	    ystart = y
          segx(nsegs,1) = xstart
	    segy(nsegs,1) = ystart
	    xlast = xstart
	    ylast = ystart
        endif
	  if(k.eq.ndata) then
	    segx(nsegs,2) = xstart
	    segy(nsegs,2) = ystart

         endif  
      enddo   
    
      else

c   comecando a ler desde a parte que interessa      
	  do while (c.le.ndata-2)

c     read second fresh pair of coordinates

	  read(38,*) y,x
c     update counter
      c = c+1


c     if x,y are identical to xstart,ystart
c     we close the polygon
        if(dabs(x-xstart).eq.0.and.dabs(y-ystart).eq.0) then



c       se diz que o termino do segmento eh, na verdade, o ponto inicial do poligono
	    segx(nsegs,2) = xstart
	    segy(nsegs,2) = ystart

        if (c.gt.ndata-3) then 

        go to 106  
        endif
c	    write(6,*) 'new seg at nsegs, ndata = ',nsegs,k

c       now we update the segment ID
        nsegs = nsegs + 1

c       discard repeated values and update counter

        read(38,*) empty
        c = c+1

c       and give it an initial coordinate based on the following coordinates
        read(38,*) y,x
        c = c+1
        xstart=x
        ystart=y


c       which will also become the new reference points
        segx(nsegs,1) = xstart
	    segy(nsegs,1) = ystart

c       if x,y differ from xstart,ystart
c       it means we are still drawing the polygon
	    else 

c     se diz que que a segunda coluna do segmento corresponde ao x,y
	    segx(nsegs,2) = x
	    segy(nsegs,2) = y

c     se muda o indice do numero de segmento
        nsegs = nsegs + 1
c     que inicia novos x,y
	    segx(nsegs,1) = x
	    segy(nsegs,1) = y

        endif

        enddo
        endif

106	continue

c      open(1,file='test')
c      write(1,*) nsegs
c      do k=1,nsegs
c        write(1,'(4f8.4)') (segx(k,j),j=1,2),(segy(k,j),j=1,2)     
c      enddo
      
c       agora, criamos um box com base na very first linha dos 
c       segmentos (xini,yini)
        box_xmax = segx(1,1)
	   box_ymax = segy(1,1)
	    box_xmin = segx(1,1)
	    box_ymin = segy(1,1)

c   passamos por todos os segmentos identificados	 
	do i=1,nsegs
	do j=1,2
c     box_x,y max/min sao, na verdade, as coordenadas
c       minimas e maximas de x e y observadas em toda a mancha
	  if(box_xmax.lt.segx(i,j)) box_xmax = segx(i,j)
	  if(box_ymax.lt.segy(i,j)) box_ymax = segy(i,j)
	  if(box_xmin.gt.segx(i,j)) box_xmin = segx(i,j)
	  if(box_ymin.gt.segy(i,j)) box_ymin = segy(i,j)
	enddo   !j
      enddo   !i 
      
c	write(1,*) 'box='
c	write(1,*) box_xmin,box_xmax
c      write(1,*) box_ymin,box_ymax
c      close(1)

	npcl_count=0
107	continue

c   agora populamos com particulas o poligono
c   calculamos uma posicao x,y aleatoria com base na box determinada
        randx = randmedslik(ix) * (box_xmax-box_xmin) + box_xmin
        randy = randmedslik(ix) * (box_ymax-box_ymin) + box_ymin
        ind = 0
c   uma vez calculada a posicao, vemos se esta dentro do poligono:
c   para cada segmento
        do i=1,nsegs
c   vemos se a particula esta entre os limites de um determinado 
c   segmento
	    if( (randx.ge.segx(i,1).and.randx.lt.segx(i,2)).or.
     &        (randx.le.segx(i,1).and.randx.gt.segx(i,2)) ) then
c   se esta respeitando os limites
            y_num = (randx - segx(i,1)) * (segy(i,2) - segy(i,1)) +
     &              segy(i,1) * (segx(i,2) - segx(i,1))
            y_dem = segx(i,2) - segx(i,1)
            if(y_dem.ne.0.d0) then
              y_int = y_num / y_dem
	        if(y_int.gt.randy) ind = ind + 1
c	      else
c              ind = ind + 1
            endif 
          endif
        enddo

        if (mod(ind,2).eq.1) then
	    npcl_count=npcl_count+1
	    px(npcl_count) = randx
          py(npcl_count) = randy
	    if(npcl_count.ge.ntot) go to 108
        endif     
      go to 107    
      
  108 continue     
         
c      open(1,file='test.lst')
c      write(1,*) 'Parcel list'
c      write(1,*) ' '
c      write(1,*) ' '
c      write(1,'(i6)') npcl_count
c	write(1,*) '   Lon     Lat'
c	do k=1,npcl_count
c	  write(1,'(2f10.6)') px(k), py(k)
c      enddo

c      close(1)
c      stop

      return
	end


c**********************************************************************
c       fetch calculation
c       (R. Lardner, M. De Dominicis, 2009)
c----------------------------------------------------------------------
	subroutine calcfetch(xavg,yavg,wdirstoke,fetch)

      implicit real*8 (a-h,o-z)
      parameter(npts=400000, imx=700, jmx=300)
      dimension alon(npts), alat(npts), seg(npts,5)
c     &          ih(imx,jmx), ih1(3000), xlon(imx), ylat(jmx),
c     &          angledeg(8), fetch(imx,jmx,8)

      character dummy*80, a3*3
      
      
	logical ex
        glon(x)=along1+(x-1.d0)*dlong
        glat(y)=alatg1+(y-1.d0)*dlatg

        pi = 4.d0 * datan(1.d0)
	degrad = pi / 180.d0
	
	
	
	m0=int(xavg+0.5d0)           ! centre of slick
	n0=int(yavg+0.5d0)
	xavg_lon=glon(dfloat(m0))
	yavg_lat=glat(dfloat(n0))
c
c     read map points
c
      open(100,file='makefetch.log')

      open(1,file='../DTM_INP/dtm.map') ! 2.01 now there is a coastline of our region
      read(1,*) ncontours

      k = 1
	nseg = 0
	do ni=1,ncontours
	  read(1,*) isle
        do i=1,isle
          read(1,*) alon(i), alat(i) 
        enddo 
	  if(isle.lt.50) go to 1
c        nseg = nseg + isle - 1
c	  write(99,*) ni,isle,nseg

        seg(k,1) = alon(1) 
        seg(k,2) = alat(1) 
        do i=2,isle
          seg(k,3) = alon(i) 
          seg(k,4) = alat(i)

          k = k + 1
          seg(k,1) = alon(i) 
          seg(k,2) = alat(i) 
        enddo 
    1   continue
      enddo
	nseg = k - 1
c
      close(1)

          cs2 = dcos(yavg_lat * degrad) **2
          xi1  = xavg_lon
          eta1 = yavg_lat

          angledeg = wdirstoke
          angle = angledeg * degrad
          csangle = dcos(angle)
          snangle = dsin(angle)

          xi2 = xi1 + 50.d0 * snangle
          eta2 = eta1 + 50.d0 * csangle
          dmin = 100.d0
          nsmin = 0

          do 62 ns = 1,nseg
           
            xx1 = seg(ns,1)
            yy1 = seg(ns,2)
            xx2 = seg(ns,3)
            yy2 = seg(ns,4)

            ddel  = (xx2-xx1)*(eta2-eta1)-(yy2-yy1)*(xi2-xi1)
            ddel1 = (eta2-eta1)*(xi1-xx1)-(xi2-xi1)*(eta1-yy1)
            ddel2 = (yy2-yy1)*(xi1-xx1)-(xx2-xx1)*(eta1-yy1)
            if(ddel.eq.0.d0) go to 62
                 
            alam=ddel1/ddel
            alamp=ddel2/ddel
            if(alam.ge.0.d0.and.alam.le.1.d0.and.alamp.gt.0.d0) then
              xx=xx1+alam*(xx2-xx1)
              yy=yy1+alam*(yy2-yy1)
              dd1=dsqrt( cs2*(xx-xi1)*(xx-xi1) + (yy-eta1)*(yy-eta1) )
              if(dd1.lt.dmin) then
                dmin=dd1
                nsmn=ns
              end if  
            end if

   62     enddo    !ns
           
             fetch = dmin
             fetch = fetch*60*1852
             if(fetch.gt.20000) fetch=20000              


       return
       end
c-------------------------------------------------------
c ERROR HANDLING SUBROUTINE
c-------------------------------------------------------
      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode
      if(errcode.ne.0) then
	print *, 'Error: ', nf_strerror(errcode)
      	stop 2
      end if
      end
