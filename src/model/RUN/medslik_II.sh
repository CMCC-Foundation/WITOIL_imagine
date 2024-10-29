#!/bin/sh
#-----------------------------------------------------------------------------------
#  MEDSLIK-II_3.01
#  oil spill fate and transport model
#-----------------------------------------------------------------------------------
#  medslik_II.sh
#  This script coordinates the model run
#-----------------------------------------------------------------------------------
#  Copyright (C) <2012>
#  This program was originally written
#  by Robin Lardner and George Zodiatis.
#  Subsequent additions and modifications
#  have been made by Michela De Dominicis.
#----------------------------------------------------------------------------------
#  The development of the MEDSLIK-II model is supported by a formal agreement
#  Memorandum of Agreement for the Operation and Continued Development of MEDSLIK-II
#  signed by the following institutions:
#  INGV - Istituto Nazionale di Geofisica e Vulcanologia
#  OC-UCY - Oceanography Center at the University of Cyprus
#  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per
#  lo Studio dell’Ambiente Marino Costiero
#  CMCC - Centro Euro-Mediterraneo sui Cambiamenti Climatici
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# str re-structured
#-----------------------------------------------------------------------------------

rm *.rso
rm medslik5.inp
rm medslik.tmp
rm ../OUT/*.fte

rm flag*.tmp
rm tmp*.tmp

rm oil_file.txt
rm initial*.txt

HOME_MEDSLIK=.
F_DATA=${HOME_MEDSLIK}/METOCE_INP/ #str

source config1.txt

python read_oil_data.py $OIL "$OIL_TYPE"

###############################################################################
#0.                               OPTIONS
###############################################################################
if [ "$ContourSlick" == "YES" ] || [ "$SAT_DATA" == "YES" ]; then
isat=1; else
isat=0
fi

###############################################################################
#1.                               SPILL INFO
###############################################################################
step_output=001       # DO NOT CHANGE IT! in hours, 3 character example: 001
output_name=out       # DO NOT CHANGE IT! 3 letters, example:out


if [ $isat -eq 0 ] || [ "$ContourSlick" == "YES" ]
then
restart=0
#hrestart=0
day=$day                   # 2 character example: 07
month=$month               # 2 character example: 08
year=$year                 # 2 character example: 08
hour=$hour                 # 2 character example: 09
minutes=$minutes           # 2 character example: 05
duration=$duration         # in hours, 4 character example: 0024
lat_degree=$lat_degree     # degrees, 2 character example: 06
lat_minutes=$lat_minutes   # minutes, example: 06.62
lon_degree=$lon_degree     # degrees, 2 character example: 06
lon_minutes=$lon_minutes   # minutes, example: 06.62
spillrate=$spillrate       # in tons/hours, example: 000055.00
spillvolume=$spillvolume   # in tons example: 000055.00
fi
iage=$age


###############################################################################
#2.                  READ CONTOUR from the input file
###############################################################################
if [ "$ContourSlick" == "YES" ]
then
echo '20'$year'/'$month'/'$day' '$hour':'$minutes >> initial.txt
p=0
N=1
while [ $N -le $NSlick ]
 do

 var_name=$`echo '{#S'$N'lon[*]}'`
 apply_count="element_count=$var_name"
 eval $apply_count

  index=1
  while [ "$index" -le "$element_count" ]
  do
  var_name_lon=$`echo '{S'$N'lon[$index]}'`
  apply_name_lon="Slon[$index]=$var_name_lon"
  eval $apply_name_lon

  var_name_lat=$`echo '{S'$N'lat[$index]}'`
  apply_name_lat="Slat[$index]=$var_name_lat"
  eval $apply_name_lat

  echo ${Slat[$index]} ${Slon[$index]} >> initial0.txt
  index=`expr  $index + 1`
  p=`expr $p + 1`

  done
  N=`expr $N + 1`
  p=`expr $p + 1`
  echo ${Slat[1]} ${Slon[1]} >> initial0.txt
done
  echo $p'        Number of data points' >> initial.txt
  echo '  lat     lon' >> initial.txt
  cat initial0.txt>>initial.txt
fi
rm initial0.txt
################################################################################
#3.                   READ INPUT DATA FROM SATELLITE DATA
################################################################################
if [ "$SAT_DATA" == "YES" ]
then

python MODEL_SRC/ReadSatData_EMSA.py $namefileGML $N_OS


myFile="medslik_sat.inp"
count=0
for var in `cat $myFile` ; do
count=`expr $count + 1`
   if [ $count -eq 2 ]; then  day=$var; fi
   if [ $count -eq 4 ]; then  month=$var; fi
   if [ $count -eq 6 ]; then  year=$var; fi
   if [ $count -eq 8 ]; then  hour=$var; fi
   if [ $count -eq 10 ]; then  minutes=$var; fi
   if [ $count -eq 12 ]; then  lat_degree=$var; fi
   if [ $count -eq 14 ]; then  lat_minutes=$var; fi
   if [ $count -eq 16 ]; then  lon_degree=$var; fi
   if [ $count -eq 18 ]; then  lon_minutes=$var; fi
   if [ $count -eq 20 ]; then  spillrate=$var; fi
   if [ $count -eq 22 ]; then  duration=$var; fi
   if [ $count -eq 24 ]; then  output_name=$var; fi
   if [ $count -eq 26 ]; then  step_output=$var; fi

done
fi
################################################################################
#4. SAVE INPUT DATA
################################################################################

     if [ $iage -eq 24 ]
     then
        Day_24=20$year$month$day
        Day_0=`./jday $Day_24 -1`
        year=`echo $Day_0 |cut -c3-4`
        month=`echo $Day_0 |cut -c5-6`
        day=`echo $Day_0 |cut -c7-8`
        sim_length=`expr $sim_length + 24`
     fi
     if [ $iage -eq 48 ]
     then
        Day_24=20$year$month$day
        Day_0=`./jday $Day_24 -2`
        year=`echo $Day_0 |cut -c3-4`
        month=`echo $Day_0 |cut -c5-6`
        day=`echo $Day_0 |cut -c7-8`
        sim_length=`expr $sim_length + 48`
     fi


multiple=01 #number of sumperimposed spills

# NUMBER OF FILES NEEDED
declare -i extra_count=0
if [ $minutes -gt 0 ]
then
hour_=$((10#$hour + 1))
else
hour_=$((10#$hour))
fi

my_int=$((10#$sim_length/24*24))
rem_hours=$((10#$sim_length-10#$my_int))
late_cases=$((10#$rem_hours+10#$hour_))
num_days=$((10#$sim_length/24))



if [ $num_days -gt 0 -a $late_cases -ge 0 ]
then

let "extra_count=$extra_count + 1"
fi

if [ $num_days -eq 0 ]
then

let "extra_count=$extra_count + 1"
fi

if [ $late_cases -ge 24 ]
then

let "extra_count=$extra_count + 1"
fi

let "numfiles=$num_days + $extra_count"

# DONE


if [ $iage -eq 24 ]
then
numfiles=`expr $numfiles + 1`
fi


#INITIAL DATE
FCStart1=20$year$month$day
FcStart1=$FCStart1
FcStart1=`echo $FcStart1 |cut -c1-8`

cp MODEL_SRC/medslikYYYY.inp medslik0.inp
sed -e "s/giorno/$day/"\
    -e "s/mese/$month/"\
    -e "s/anno/$year/"\
    -e "s/ora/$hour/"\
    -e "s/minuti/$minutes/"\
    -e "s/durata/$duration/"\
    -e "s/lat_gradi/$lat_degree/"\
    -e "s/lat_primi/$lat_minutes/"\
    -e "s/lon_gradi/$lon_degree/"\
    -e "s/lon_primi/$lon_minutes/"\
    -e "s/nome/$output_name/"\
    -e "s/lunghezza/$sim_length/"\
    -e "s/step/$step_output/"\
    -e "s/eta/$iage/"\
    -e "s/sat/$isat/"\
    -e "s/portata/$spillrate/"\
    -e "s/multi/$multiple/"\
    -e "s/riinizio/$restart/"\
    -e "s/griglia/$grid_size/"\
    -e "s/numero_files/$numfiles/"\
    medslik0.inp>medslik1.inp
rm medslik0.inp

sed -n '1,15 p' medslik1.inp >medslik5.inp
cat oil_file.txt>>medslik5.inp
sed -n '24,28 p' medslik1.inp >medslik0.inp
cat medslik0.inp>>medslik5.inp
rm medslik[01].inp

#####################################################################
#5. PRE-PROCESSING OF CURRENTS & WIND FILES NEEDED FOR SIMULATION
#####################################################################
# some info on how your current files are named and organized
dir='OCE'
U='vozocrtx'
V='vomecrty'
T='votemper'
model=''
pre_name='MDK_ocean_'
tail_name=''


PREPROC=PREPROC
FD=$F_DATA/PREPROC/$dir

rm tmp*.tmp

n=0
loop=$numfiles

while [ $n != $loop ]; do
n=`expr $n + 1`
nn=`expr $n - 1`
Datafc=`./jday $FcStart1 +$nn`
DataFC=`echo $Datafc |cut -c3-8`
DataFc=$DataFC
echo 'CHECK CURRENTS' $DataFc
DataFc_out=${DataFc}


if [ -f $FD/$pre_name${DataFc}$tail_name'_U.nc' ]
then
   ln -s $FD/$pre_name${DataFc}$tail_name'_U.nc' $FD/${DataFc}'_U.nc'
   echo "${DataFc}_U.nc 1" >> tmp1.tmp
else
   echo "${DataFc}_U.nc 1" >> tmp1.tmp
   echo "For this run you need the file" $pre_name${DataFc}$tail_name'_U.nc'
fi

if [ -f $FD/$pre_name${DataFc}$tail_name'_V.nc' ]
then
   ln -s $FD/$pre_name${DataFc}$tail_name'_V.nc' $FD/${DataFc}'_V.nc'
   echo "${DataFc}_V.nc 1" >> tmp1.tmp
else
   echo "${DataFc}_V.nc 1" >> tmp1.tmp
   echo "For this run you need the file" $pre_name${DataFc}$tail_name'_V.nc'
fi

if [ -f $FD/$pre_name${DataFc}$tail_name'_T.nc' ]
then
   ln -s $FD/$pre_name${DataFc}$tail_name'_T.nc' $FD/${DataFc}'_T.nc'
   echo "${DataFc}_T.nc 1" >> tmp1.tmp
else
   echo "${DataFc}_T.nc 1" >> tmp1.tmp
   echo "For this run you need the file" $pre_name${DataFc}$tail_name'_T.nc'
fi

dir_wind='MET'
FD_wind_in=$F_DATA/$PREPROC/$dir_wind
FD_wind=$F_DATA/$PREPROC/$dir_wind

Datafc=$Datafc

DataFC=`echo $Datafc |cut -c3-8`
DataFc=$DataFC

Data_wind="20"${DataFc}
File_wind="erai"${DataFc}".eri"

echo 'CHECK WINDS' $Data_wind
if [ -f $FD_wind_in/${Data_wind}.nc ]
then
   echo "${File_wind} 1" >> tmp2.tmp
   echo $FD_wind_in/${Data_wind}.nc ": File exists"
else
   echo "${File_wind} 1" >> tmp2.tmp
   echo "For this run you need the file" ${Data_wind}.nc
fi

done


cat tmp1.tmp>>medslik5.inp
echo $numfiles >> medslik5.inp

cat tmp2.tmp>>medslik5.inp
cat medslik_multi.tmp>>medslik5.inp

#############################################################
#6. AREA SELECTION (MIN/MAX Longitudes & Latitudes)
#############################################################
cp MODEL_SRC/medslikYYYY.tmp medslik1.tmp

./lat_lon.exe
mv medslik1.tmp medslik.tmp


echo $numfiles >> medslik.tmp
numfiles_tmp=$numfiles
# done
n=0
while [ $n != $numfiles_tmp ]; do
n=`expr $n + 1`
nn=`expr $n - 1`
Datafc=`./jday $FcStart1 +$nn`
DataFC=`echo $Datafc |cut -c3-8`
DataFc=$DataFC
echo ${DataFc}24 >> medslik.tmp
done

echo $numfiles_tmp >> medslik.tmp


n=0
while [ $n != $numfiles_tmp ]; do
n=`expr $n + 1`
nn=`expr $n - 1`
Datafc=`./jday $FcStart1 +$nn`
DataFC=`echo $Datafc |cut -c3-8`
DataFc=$DataFC
echo ${Datafc} >> medslik.tmp
done

echo " 0" >> medslik.tmp


tr -d "\015" < medslik5.inp > medslik52.inp
cp medslik52.inp medslik5.inp
rm medslik52.inp

tr -d "\015" < medslik.tmp > medslik2.tmp
cp medslik2.tmp medslik.tmp
rm medslik2.tmp

#######################################################################
# CREATE OUTPUT DIRECTORY
#######################################################################
DIR_output='MDK_SIM_20'$year'_'$month'_'$day'_'$hour$minutes'_'$SIM_NAME
mkdir ../OUT/$DIR_output
#######################################################################
#7. EXTRACT CURRENTS, WIND AND SST DATA
#####################################################################
echo 'READING CURRENTS & WIND DATA'
# ./Extract_II.exe $F_DATA
######################################################################
#8. RUN
#####################################################################
echo 'PLEASE WAIT: SIMULATION IS RUNNING'
./medslik_II.exe
#####################################################################
#9. ARCHIVE OUTPUT FILES
#####################################################################
rm tmp*.tmp

cp medslik5.inp  ../OUT/$DIR_output
cp config2.txt  ../OUT/$DIR_output
cp medslik.tmp ../OUT/$DIR_output
cp initial.txt ../OUT/$DIR_output
cp config1.txt ../OUT/$DIR_output
mv ../OUT/*.tot ../OUT/$DIR_output
mv ../OUT/*.fte ../OUT/$DIR_output
mv ../OUT/*.cst ../OUT/$DIR_output
mv spill_properties.nc ../OUT/$DIR_output

mkdir ../OUT/$DIR_output/MET
mkdir ../OUT/$DIR_output/OCE

mv obs* ../OUT/$DIR_output

mv TEMP/MET/*.* ../OUT/$DIR_output/MET
mv TEMP/OCE/*.* ../OUT/$DIR_output/OCE

