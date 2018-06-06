echo ' Note: on norstore cruncher, run "module load ncl" first!'
echo ''
echo 'This script calls all the *ModIvsModII.ncl scripts'
echo 'and produces plots for all available/listed plot_types'
echo 'Note: All ncl scripts assumes that the input data on'
echo 'the listed directories are on a integer number times'
echo '12 nc-files for monthly model data for model I and II,'
echo 'and that no other files on the same name form are present.' 
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
#
#pthI=/projects/NS2345K/noresm/cases/aerocomA2noresm_r128/tests_feb06-feb08/CTRL2000/yr7/
#fnmI=aerocomA2r128_2006.cam2.h0.0007-01.nc
#fnmpI=aerocomA2r128_2006.cam2.h0.0007
#
#pthI=/projects/NS2345K/noresm/cases/CCNCOMP06V_yr2011/
#fnmI=CCNCOMP06V.cam.h0.2011-01.nc
#fnmpI=CCNCOMP06V.cam.h0.2011-
#
#pthI=/projects/NS2345K/noresm/cases/53OSLO_PD_UNTUNED_yr2006-2010/
#fnmI=53OSLO_PD_UNTUNED.cam.h0.2006-2010_01.nc
#fnmpI=53OSLO_PD_UNTUNED.cam.h0.2006-2010_
#
#pthI=/projects/NS2345K/noresm/cases/53OSLO_PDcorrlifeNMR12/
#fnmI=53OSLO_PDcorrlifeNMR12.cam.h0.2006-2010_01.nc
#fnmpI=53OSLO_PDcorrlifeNMR12.cam.h0.2006-2010_
#
pthI=/projects/NS2345K/noresm/cases/F2000ERF/mnthclim_yr3-30/
fnmI=F2000ERF.cam.h0.0003-0030_01.nc
fnmpI=F2000ERF.cam.h0.0003-0030_
#
#pthI=/projects/NS2345K/noresm/cases/NF1998HygCfree/orig/
#fnmI=NF1998HygCfree.cam.h0.0001-05.nc
#fnmpI=NF1998HygCfree.cam.h0.000
#
#pthI=/projects/NS2345K/noresm/cases/F2000ERF/2-07to3-06/
#fnmI=F2000ERF.cam.h0.0002-07.nc
#fnmpI=F2000ERF.cam.h0.000
#
#pthI=/projects/NS2345K/noresm/cases/7310AMIP20002_mnthclim_yr3-30/
#fnmI=7310AMIP20002.cam.h0.0003-0030_01.nc
#fnmpI=7310AMIP20002.cam.h0.0003-0030_
#
#pthII=/projects/NS2345K/noresm/cases/F2000ERF/mnthclim_yr3-30/
#fnmII=F2000ERF.cam.h0.0003-0030_01.nc
#fnmpII=F2000ERF.cam.h0.0003-0030_
#
#pthII=/projects/NS2345K/noresm/cases/53OSLO_PDnewDust/
#fnmII=53OSLO_PDnewDust.cam.h0.2006-2010_01.nc
#fnmpII=53OSLO_PDnewDust.cam.h0.2006-2010_
#
#pthII=/projects/NS2345K/noresm/cases/53OSLO_PDnoDust/
#fnmII=53OSLO_PDnoDust.cam.h0.2006-2010_01.nc
#fnmpII=53OSLO_PDnoDust.cam.h0.2006-2010_
#
#pthII=/projects/NS2345K/noresm/cases/AlfKtests/53OsloTstNuc/tst3_7yrs/yr6-10/
#fnmII=53OsloTstNuc.cam.h0.2006-2010_01.nc
#fnmpII=53OsloTstNuc.cam.h0.2006-2010_
#
#pthII=/projects/NS2345K/noresm/cases/NF1998ERF/mnthclim_yr3-30/
#fnmII=NF1998ERF.cam.h0.0003-0030_01.nc
#fnmpII=NF1998ERF.cam.h0.0003-0030_
#
#pthII=/projects/NS2345K/noresm/cases/NF1998HygCfree/
#fnmII=NF1998HygCfree.cam.h0.0001-05.nc
#fnmpII=NF1998HygCfree.cam.h0.000
#
#pthII=/projects/NS2345K/noresm/cases/AlfKtests/CAM6-Oslo/NF2kNucl-all/mnth5-17b4corrsconc/
#fnmII=NF2kNucl.cam.h0.0001-05.nc
#fnmpII=NF2kNucl.cam.h0.000
#
pthII=/projects/NS2345K/noresm/cases/AlfKtests/CAM6-Oslo/NF2kNucl-all/
fnmII=NF2kNucl.cam.h0.0003-01.nc
fnmpII=NF2kNucl.cam.h0.000
#
#pthII=/projects/NS2345K/noresm/cases/53OSLO_PDcorrlifeNMR12/
#fnmII=53OSLO_PDcorrlifeNMR12.cam.h0.2006-2010_01.nc
#fnmpII=53OSLO_PDcorrlifeNMR12.cam.h0.2006-2010_
#
# Paths and names of input files necessary for forcing plots (PI case)):
#
#pthI_PI=/projects/NS2345K/noresm/cases/aerocomA2noresm_r128/tests_feb06-feb08/CTRL2000/PRE1850_yr7/
#fnmI_PI=aerocomA2r128_1850.cam2.h0.0007-01.nc
#fnmpI_PI=aerocomA2r128_1850.cam2.h0.0007
#
#pthI_PI=/projects/NS2345K/noresm/cases/CCNCOMP06V_yr2011/
#fnmI_PI=CCNCOMP06V.cam.h0.2011-01.nc
#fnmpI_PI=CCNCOMP06V.cam.h0.2011-
#
#pthI_PI=/projects/NS2345K/noresm/cases/53OSLO_PI_wPDoxi/
#fnmI_PI=53OSLO_PI_wPDoxi.cam.h0.2006-2010_01.nc
#fnmpI_PI=53OSLO_PI_wPDoxi.cam.h0.2006-2010_
#
#pthI_PI=/projects/NS2345K/noresm/cases/53OSLO_PIcorrlifeNMR12/
#fnmI_PI=53OSLO_PIcorrlifeNMR12.cam.h0.2006-2010_01.nc
#fnmpI_PI=53OSLO_PIcorrlifeNMR12.cam.h0.2006-2010_
#
#pthI_PI=/projects/NS2345K/noresm/cases/NF2000PIERF/mnthclim_yr3-30/
#fnmI_PI=NF2000PIERF.cam.h0.0003-0030_01.nc
#fnmpI_PI=NF2000PIERF.cam.h0.0003-0030_
#
#pthI_PI=/projects/NS2345K/noresm/cases/NF2000PIERF/2-07to3-06/
#fnmI_PI=NF2000PIERF.cam.h0.0002-07.nc
#fnmpI_PI=NF2000PIERF.cam.h0.000
#
pthI_PI=/projects/NS2345K/noresm/cases/NF2000PIERF/mnthclim_yr3-30/
fnmI_PI=NF2000PIERF.cam.h0.0003-0030_01.nc
fnmpI_PI=NF2000PIERF.cam.h0.0003-0030_
#
#pthI_PI=/projects/NS2345K/noresm/cases/NF1998HygCfree/orig/
#fnmI_PI=NF1998HygCfree.cam.h0.0001-05.nc
#fnmpI_PI=NF1998HygCfree.cam.h0.000
#
#pthII_PI=/projects/NS2345K/noresm/cases/AlfKtests/CAM6-Oslo/NF2kNuclPI-all/2-07to3-06/
#fnmII_PI=NF2kNuclPI.cam.h0.0002-07.nc
#fnmpII_PI=NF2kNuclPI.cam.h0.000
#
#pthII_PI=/projects/NS2345K/noresm/cases/AlfKtests/CAM6-Oslo/NF2kNucl-all/
#fnmII_PI=NF2kNucl.cam.h0.0001-05.nc
#fnmpII_PI=NF2kNucl.cam.h0.000
#
pthII_PI=/projects/NS2345K/noresm/cases/AlfKtests/CAM6-Oslo/NF2kNuclPI-all/
fnmII_PI=NF2kNuclPI.cam.h0.0003-01.nc
fnmpII_PI=NF2kNuclPI.cam.h0.000
#
#pthI_PI=/projects/NS2345K/noresm/cases/F2000ERF/
#fnmI_PI=F2000ERF.cam.h0.0001-05.nc
#fnmpI_PI=F2000ERF.cam.h0.000
#
#pthI_PI=/projects/NS2345K/noresm/cases/53OSLO_PInewDust/
#fnmI_PI=53OSLO_PInewDust.cam.h0.2006-2010_01.nc
#fnmpI_PI=53OSLO_PInewDust.cam.h0.2006-2010_
#
#pthII_PI=/projects/NS2345K/noresm/cases/53OSLO_PIcorrlifeNMR12/
#fnmII_PI=53OSLO_PIcorrlifeNMR12.cam.h0.2006-2010_01.nc
#fnmpII_PI=53OSLO_PIcorrlifeNMR12.cam.h0.2006-2010_
#
#pthII_PI=/projects/NS2345K/noresm/cases/53OSLO_PDnoDust/
#fnmII_PI=53OSLO_PDnoDust.cam.h0.2006-2010_01.nc
#fnmpII_PI=53OSLO_PDnoDust.cam.h0.2006-2010_
#
#pthII_PI=/projects/NS2345K/noresm/cases/53OSLO_PI_wPDoxi/
#fnmII_PI=53OSLO_PI_wPDoxi.cam.h0.2006-2010_01.nc
#fnmpII_PI=53OSLO_PI_wPDoxi.cam.h0.2006-2010_
#
#pthII_PI=/projects/NS2345K/noresm/cases/AMIP_PI_wPDoxi/atm/hist/mnthclim_yr3-30/
#fnmII_PI=AMIP_PI_wPDoxi.cam.h0.0003-0030_01.nc
#fnmpII_PI=AMIP_PI_wPDoxi.cam.h0.0003-0030_
#
#pthII_PI=/projects/NS2345K/noresm/cases/AlfKtests/PI53OsloTstNuc/tst3_7yrs/yr6-10/
#fnmII_PI=PI53OsloTstNuc.cam.h0.2006-2010_01.nc
#fnmpII_PI=PI53OsloTstNuc.cam.h0.2006-2010_
#
#pthII_PI=/projects/NS2345K/noresm/cases/F2000PIERFfram/
#fnmII_PI=F2000PIERF.cam.h0.0001-05.nc
#fnmpII_PI=F2000PIERF.cam.h0.000
#pthII_PI=/projects/NS2345K/noresm/cases/NF2000PIERF/mnthclim_yr3-30/
#fnmII_PI=NF2000PIERF.cam.h0.0003-0030_01.nc
#fnmpII_PI=NF2000PIERF.cam.h0.0003-0030_
#
#pthII_PI=/projects/NS2345K/noresm/cases/NF1998HygCfree/
#fnmII_PI=NF1998HygCfree.cam.h0.0001-05.nc
#fnmpII_PI=NF1998HygCfree.cam.h0.000
#
#pthII_PI=/projects/NS2345K/noresm/cases/AlfKtests/CAM6-Oslo/NF2kNucl-all/mnth5-17b4corrsconc/
#fnmII_PI=NF2kNucl.cam.h0.0001-05.nc
#fnmpII_PI=NF2kNucl.cam.h0.000
#
#pthII_PI=/projects/NS2345K/noresm/cases/53OSLO_PIinclLndUse/
#fnmII_PI=53OSLO_PIinclLndUse.cam.h0.2006-2010_01.nc
#fnmpII_PI=53OSLO_PIinclLndUse.cam.h0.2006-2010_
#pthII_PI=/projects/NS2345K/noresm/cases/B05NudgePI4_2013-2015/
#fnmII_PI=B05NudgePI4.cam.h0.2010-01.nc
#fnmpII_PI=B05NudgePI4.cam.h0.2010      
#
#ModelI=CAM4-Oslo  # gives CAM4-Oslo vs. new CAM5-Oslo comparison plots
ModelI=CAM5-Oslo  # gives CAM5/6-Oslo Revision N vs. CAM5/6-Oslo Revision M comparison plots
#
# Note: when using CAM5.5 or higher for ModI, include also GasdepI options, as for ModII below
GasdepI=Neu   # use diagnostics as in the Moz or in the Neu (>CAM5.3) gas deposition scheme for ModI
GasdepII=Neu  # use diagnostics as in the Moz or in the Neu (>CAM5.3) gas deposition scheme for ModII
# Note: when using 32 level version for ModI or ModII, include also the LevModI & LevModII options L32
# (this does not work if ModI is CAM4-Oslo), for use (only) in LevelCloudProp_ModIvsModII.ncl 
LevModI=L32 
LevModII=L32 
#**********************************************************************************************


# Run TMP version x instead of the above for ModII data before correction w.r.t. mmr vs. vmr in model output of individual tracers
# or the y version if the ModI data have this shortcoming, or z if both ModI and ModII do: 
 echo ''
 echo 'Running ZonalAero_ModIvsModII_TMP.ncl'
 echo ''
 for I in {1..7};do
 ncl  plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep#_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalAero_ModIvsModII.ncl
# ncl  plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep#_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalAero_ModIvsModII_TMPx.ncl
# ncl  plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalAero_ModIvsModII_TMPy.ncl
# ncl  plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalAero_ModIvsModII_TMPz.ncl
 done

#No changes by the user should be necessary below...

echo ''
echo 'Running Emis_ModIvsModII.ncl'
echo ''
for I in {-1..7};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" Emis_ModIvsModII.ncl
done

echo ''
echo 'Running diffTOAbalance.ncl'
echo ''
for I in {0..9};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" diffTOAbalance.ncl 
done

echo ''
echo 'Running LevelCloudProp_ModIvsModII.ncl'
echo ''
for I in {1..10};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" LevModI=\"$LevModI\" LevModII=\"$LevModII\" LevelCloudProp_ModIvsModII.ncl
done

echo ''
echo 'Running Cld2d_ModIvsModII.ncl'
echo ''
for I in {0..9};do
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
for I in {0..8};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" AOD_ModIvsModII.ncl
done

echo ''
echo 'Running ZonalRHCl_ModIvsModII.ncl'
echo ''
for I in {0..7};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalRHCl_ModIvsModII.ncl
done

echo ''
echo 'Running Lifetimes_ModIvsModII.ncl'
echo ''
for I in {0..5};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" GdepI=\"$GasdepI\" GdepII=\"$GasdepII\" Lifetimes_ModIvsModII.ncl
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
echo 'Running ZonalN_ModIvsModII.ncl'
echo ''
for I in {1..13};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" ZonalN_ModIvsModII.ncl
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
for I in {0..9};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" divPD-PI_ModIvsModII.ncl 
done

echo ''
echo 'Running divPD-PI_Zonal_ModIvsModII.ncl'
echo ''
for I in {0..6};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" divPD-PI_Zonal_ModIvsModII.ncl 
done

echo ''
echo 'Running ERFsurf_ModIvsModII.ncl'
echo ''
for I in {0..4};do
ncl  plot_type=$I format=\"$plotf\" filepathPD_I=\"$pthI\" filenamePD_I=\"$fnmI\" filepathPD_II=\"$pthII\" filenamePD_II=\"$fnmII\" filenamepPD_I=\"$fnmpI\" filenamepPD_II=\"$fnmpII\" filepathPI_I=\"$pthI_PI\" filenamePI_I=\"$fnmI_PI\" filepathPI_II=\"$pthII_PI\" filenamePI_II=\"$fnmII_PI\" filenamepPI_I=\"$fnmpI_PI\" filenamepPI_II=\"$fnmpII_PI\" ModI=\"$ModelI\" ERFsurf_ModIvsModII.ncl 
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
for I in {-1..5};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI\" filename_I=\"$fnmI\" filepath_II=\"$pthII\" filename_II=\"$fnmII\" filenamep_I=\"$fnmpI\" filenamep_II=\"$fnmpII\" ModI=\"$ModelI\" GdepI=\"$GasdepI\" GdepII=\"$GasdepII\" Mass-budget_ModIvsModII.ncl
done

echo ''
echo 'Running Mass-budget_ModIvsModII.ncl for PI emissions'
echo ''
for I in {-1..5};do
 ncl plot_type=$I format=\"$plotf\" filepath_I=\"$pthI_PI\" filename_I=\"$fnmI_PI\" filepath_II=\"$pthII_PI\" filename_II=\"$fnmII_PI\" filenamep_I=\"$fnmpI_PI\" filenamep_II=\"$fnmpII_PI\" ModI=\"$ModelI\" GdepI=\"$GasdepI\" GdepII=\"$GasdepII\" Mass-budget_ModIvsModII.ncl
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


