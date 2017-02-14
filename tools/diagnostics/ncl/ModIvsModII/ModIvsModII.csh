echo ''
echo 'This script calls all the *ModIvsModII.ncl scripts'
echo 'and produces plots for all available/listed plot_types'
echo 'Note: All ncl scripts assumes that the input data on'
echo 'the listed directories are on a integer number times'
echo '12 nc-files for monthly model data for model I and II,'
echo 'and that no other files are present there.' 
echo ''
echo 'If the number of years is so large that some of the'
echo 'ncl scripts run out of memory, with the error message'
echo 'systemfunc: cannot create child process:[errno=12],'
echo 'then make climatological input for each mponth in,'
echo 'advance, e.g. by use of the ncea command...'
echo ''
echo '     (Created by Alf Kirkevåg, April 2014)'
echo ''

# ncl 'dataFile=addfile("./modelData.nc", "r")' plot_type=0 Emis_ModIvsModII.ncl
#************************* To be edited by the user ********************************************
# Plot type and plot output format:
plotf=png  # chosen output format for figures (ps, eps, pdf, png)
#
# Paths and names of input files from model version I and II (PD case): 
#pthI=/vol/fou/atmos2/alfk/hexagon/aerocomA2noresm_r128/tests/CTRL2000/
#pthI=/media/BackupAK/aerocomA2r128-tester/CTRL2000/
#fnmI=aerocomA2r128_2006.cam2.h0.0007-01.nc     # first file to be used
#fnmpI=aerocomA2r128_2006.cam2.h0.0007-         # all files fnmpI* are used
#pthI=/lustre/mnt/alfg/condTimeStep/R554SoaNucl/atm/hist/
#fnmI=R554SoaNucl.cam.h0.1981-01.nc                 # first file to be used
#fnmpI=R554SoaNucl.cam.h0.1981                        # all files fnmpI_PD* are used
#pthI=/media/BackupAK/NorESM2output/PDndgMG15MegVadSOA/
#fnmI=PDndgMG15MegVadSOA.cam.h0.1983-01.nc 
#fnmpI=PDndgMG15MegVadSOA.cam.h0.19
#pthI=/media/BackupAK/NorESM2output/PD_MG15MegVadSOA/clim3-12corr/
#fnmI=PD_MG15MegVadSOA.cam.h0.climyr3-12_01.nc
#fnmpI=PD_MG15MegVadSOA.cam.h0.clim
#pthI=/media/BackupAK/NorESM2output/PDaug16UVPSndg/yr18-19/
#fnmI=PDaug16UVPSndg.cam.h0.0018-01.nc
#fnmpI=PDaug16UVPSndg.cam.h0.00
#pthI=/media/BackupAK/NorESM2output/PD_ERA_2001-2015/
#fnmI=ERA_2001-2015.cam.h0.2012-01.nc
#fnmpI=ERA_2001-2015.cam.h0.20
pthI=/media/BackupAK/NorESM2output/PDnewBCfixedSOA/
fnmI=PDnewBCfixedSOA.cam.h0.2010-01.nc
fnmpI=PDnewBCfixedSOA.cam.h0.20
#
#pthII=/lustre/mnt/alfg/condTimeStep/R566SoaNucl/atm/hist/
#fnmII=R566SoaNucl.cam.h0.1981-01.nc                 # first file to be used
#fnmpII=R566SoaNucl.cam.h0.1981                        # all files fnmpII* are used
#pthII=/media/BackupAK/NorESM2output/PD_MG15MegVadSOA/
#fnmII=PD_MG15MegVadSOA.cam.h0.0010-01.nc
#fnmpII=PD_MG15MegVadSOA.cam.h0.00
#pthII=/media/BackupAK/NorESM2output/PD_MG15MegVadSOA/clim3-12corr/
#fnmII=PD_MG15MegVadSOA.cam.h0.climyr3-12_01.nc
#fnmpII=PD_MG15MegVadSOA.cam.h0.clim
pthII=/media/BackupAK/NorESM2output/PD-opticsINSITU_Vilje/
fnmII=PD-opticsINSITU.cam.h0.2010-01.nc
fnmpII=PD-opticsINSITU.cam.h0.20
#
# Paths and names of input files necessary for forcing plots (PI case)):
#pthI_PI=/vol/fou/atmos2/alfk/hexagon/aerocomA2noresm_r128/
#fnmI_PI=aerocomA2r128_1850.cam2.h0.0007-01.nc     # first file to be used
#fnmpI_PI=aerocomA2r128_1850.cam2.h0.0007-         # all files fnmpI_PI* are used
#pthI_PI=/media/BackupAK/NorESM2output/SOA_r610_PI/
#fnmI_PI=SOA_r610_PI.cam.h0.1984-12.nc
#fnmpI_PI=SOA_r610_PI.cam.h0.19
#pthI_PI=/media/BackupAK/NorESM2output/PIndgMG15MegVadSOA/
#fnmI_PI=PIndgMG15MegVadSOA.cam.h0.1983-01.nc
#fnmpI_PI=PIndgMG15MegVadSOA.cam.h0.19
#pthI_PI=/media/BackupAK/NorESM2output/PI_MG15MegVadSOA/clim3-12corr/
#fnmI_PI=PI_MG15MegVadSOA.cam.h0.climyr3-12_01.nc
#fnmpI_PI=PI_MG15MegVadSOA.cam.h0.clim
#pthI_PI=/media/BackupAK/NorESM2output/PIaug16UVPSndg/yr18-19/
#fnmI_PI=PIaug16UVPSndg.cam.h0.0018-01.nc
#fnmpI_PI=PIaug16UVPSndg.cam.h0.00
#pthI_PI=/media/BackupAK/NorESM2output/PI_ERA_2001-2015/
#fnmI_PI=PI_ERA_2001-2015.cam.h0.2012-01.nc
#fnmpI_PI=PI_ERA_2001-2015.cam.h0.20
pthI_PI=/media/BackupAK/NorESM2output/PInewBCfixedSOA/
fnmI_PI=PInewBCfixedSOA.cam.h0.2010-01.nc
fnmpI_PI=PInewBCfixedSOA.cam.h0.20

#pthII_PI=/media/BackupAK/NorESM2output/PI1850R516/
#fnmII_PI=PI1850R516.cam.h0.1984-12.nc                 # first file to be used
#fnmpII_PI=PI1850R516.cam.h0.19                        # all files fnmpII* are used
#pthII_PI=/media/BackupAK/NorESM2output/PI_MG15MegVadSOA/
#fnmII_PI=PI_MG15MegVadSOA.cam.h0.0010-01.nc
#fnmpII_PI=PI_MG15MegVadSOA.cam.h0.00
#pthII_PI=/media/BackupAK/NorESM2output/PI_MG15MegVadSOA/clim3-12corr/
#fnmII_PI=PI_MG15MegVadSOA.cam.h0.climyr3-12_01.nc
#fnmpII_PI=PI_MG15MegVadSOA.cam.h0.clim
pthII_PI=/media/BackupAK/NorESM2output/PI-opticsINSITU_Vilje/
fnmII_PI=PI-opticsINSITU.cam.h0.2010-01.nc
fnmpII_PI=PI-opticsINSITU.cam.h0.20
#
#ModelI=CAM4-Oslo  # gives CAM4-Oslo vs. new CAM5-Oslo comparison plots
ModelI=CAM5-Oslo  # gives CAM5-Oslo Revision N vs. CAM5-Oslo Revision M comparison plots
#**********************************************************************************************
#No changes by the user should be necessary below...


echo ''
echo 'Running Emis_ModIvsModII.ncl'
echo ''
for I in {0..5};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" Emis_ModIvsModII.ncl
done

echo ''
echo 'Running Cld2d_ModIvsModII.ncl'
echo ''
for I in {0..8};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" Cld2d_ModIvsModII.ncl 
done

echo ''
echo 'Running Load_ModIvsModII.ncl'
echo ''
for I in {0..19};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" Load_ModIvsModII.ncl
done

echo ''
echo 'Running Ext_ModIvsModII.ncl'
echo ''
for I in {1..20};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" Ext_ModIvsModII.ncl 
done

echo ''
echo 'Running AODratio_ModIvsModII.ncl'
echo ''
for I in {0..6};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" AODratio_ModIvsModII.ncl
done

echo ''
echo 'Running AOD_ModIvsModII.ncl'
echo ''
for I in {0..6};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" AOD_ModIvsModII.ncl
done

echo ''
echo 'Running ZonalAero_ModIvsModII.ncl'
echo ''
for I in {0..7};do
 ncl  plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalAero_ModIvsModII.ncl
done

echo ''
echo 'Running ZonalRHCl_ModIvsModII.ncl'
echo ''
for I in {0..6};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalRHCl_ModIvsModII.ncl
done

echo ''
echo 'Running Lifetimes_ModIvsModII.ncl'
echo ''
for I in {0..5};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" Lifetimes_ModIvsModII.ncl
done

echo ''
echo 'Running WetDepRat_ModIvsModII.ncl'
echo ''
for I in {0..5};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" WetDepRat_ModIvsModII.ncl
done

echo ''
echo 'Running EffDryRad_ModIvsModII.ncl'
echo ''
for I in {0..2};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" EffDryRad_ModIvsModII.ncl 
done

echo ''
echo 'Running ZonalModepar_ModIvsModII.ncl'
echo ''
for I in {1..9};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalModepar_ModIvsModII.ncl
done

echo ''
echo 'Running PM_ModIvsModII.ncl'
echo ''
for I in {1..5};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" PM_ModIvsModII.ncl 
done

echo ''
echo 'Running RadBudg_ModIvsModII.ncl'
echo ''
for I in {1..3};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" RadBudg_ModIvsModII.ncl 
done

echo ''
echo 'Running divPD-PI_ModIvsModII.ncl'
echo ''
for I in {0..7};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" divPD-PI_ModIvsModII.ncl 
done

echo ''
echo 'Running divPD-PI_Zonal_ModIvsModII.ncl'
echo ''
for I in {0..6};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" divPD-PI_Zonal_ModIvsModII.ncl 
done

echo ''
echo 'Running ERF_ModIvsModII.ncl'
echo ''
for I in {0..6};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" ERF_ModIvsModII.ncl 
done

echo ''
echo 'Running Mass-budget_ModIvsModII.ncl'
echo ''
for I in {0..5};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" Mass-budget_ModIvsModII.ncl
done

echo ''
echo 'All ncl script runs completed'
echo ''

echo "trim whitespace in images"
for i in `ls *.png`
do
  convert -trim $i $i
done

exit
