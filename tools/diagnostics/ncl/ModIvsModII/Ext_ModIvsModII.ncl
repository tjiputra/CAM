  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; This ncl script reads in aerosol loads and AOD from two versions of 
; NorESM/CAM-Oslo and makes global plots of the annually averaged 
; aerosol extinction coefficient, defined as AOD/load, including global 
; average as a number in the title line for each figure. In this updated
; version of the script the global average is be calculated both as area 
; averaged EXT and as area averaged AOD / area averaged aerosol column 
; burden (load), with its calculated value given in brackets.

; Model independent constants
 g=9.80665
 pi=3.1415926
 re=6378.39e3            ; earth radius in m
 coffa=pi*re^2./180.      
 area1=4.*pi*re^2          

; *************************************************************************
; **** To be edited by the user if the ncl script is run interactively ****
;
; Define plot type and plot output format
 if (.not. isvar("plot_type")) then      ; is plot_type on command line?
  plot_type = 1                      ; 1 => SO4 extinction
                                     ; 2 => BC  extinction
                                     ; 3 => POM extinction
                                     ; 4 => SS  extinction
                                     ; 5 => DU  extinction
 end if 
 if (.not. isvar("format")) then        ; is format on command line?
  format = "ps"
;  format = "eps"
;  format = "png"
;  format = "pdf"
 end if
;
; *************************************************************************
; No changes by the user should be necessary below...
; *************************************************************************
;old  all_files_I  = systemfunc ("ls /media/BackupAK/aerocomA2r128-tester/CTRL2000/aerocomA2r128_2006.cam2.h0.0007-*.nc")
      all_files_I  = systemfunc ("ls " + filepath_I + filenamep_I + "*")
      all_files_II = systemfunc ("ls " + filepath_II + filenamep_II + "*")
      f0_I      = addfile (filepath_I+filename_I, "r")
      f0_II     = addfile (filepath_II+filename_II, "r")
      f1_I      = addfiles (all_files_I, "r")   ; note the "s" of addfile
      f1_II     = addfiles (all_files_II, "r")   ; note the "s" of addfile

; Reading Gaussian weights and other required model variables
 gw0_I=doubletofloat(f0_I->gw)
 gw0_II=doubletofloat(f0_II->gw)

 lon_I=f0_I->lon
 dlon_I=360./dimsizes(lon_I)
 lon_II=f0_II->lon
 dlon_II=360./dimsizes(lon_II)
  
; Initialization (and obtain correct variable dimensions)
 tmp_I=f1_I[:]->PS
 tmp_II=f1_II[:]->PS         
 ext_I=tmp_I
 ext_II=tmp_II
 aod_I=tmp_I
 aod_II=tmp_II
 load_I=tmp_I
 load_II=tmp_II

 if (plot_type.eq.1) then
 var="EXT_SO4"                             ; name of input-variable and plot
 varname="SO~B~4~N~ extinction coefficient"      ; variable name used in text string:
 aod_I=(/f1_I[:]->D550_SO4/)                                 ; variable to be plotted from I
 aod_II=(/f1_II[:]->D550_SO4/)                               ; variable to be plotted from II 
 if(ModI.eq."CAM4-Oslo") then
   load_I=(/f1_I[:]->C_SO4/)*1.e3                              ; variable to be plotted from I
   ext_I=(/f1_I[:]->D550_SO4/)/(/(f1_I[:]->C_SO4)/)*1.e-3      ; variable to be plotted from I
 else
   load_I=(/(f1_I[:]->cb_SO4_A1)/)/3.06 + (/(f1_I[:]->cb_SO4_A2)/)/3.59 + (/(f1_I[:]->cb_SO4_AC)/)/3.06 + (/(f1_I[:]->cb_SO4_NA)/)/3.06 + (/(f1_I[:]->cb_SO4_PR)/)/3.06
   load_I=load_I*1.e3 
   ext_I=aod_I/load_I
 end if
 load_II=(/(f1_II[:]->cb_SO4_A1)/)/3.06 + (/(f1_II[:]->cb_SO4_A2)/)/3.59 + (/(f1_II[:]->cb_SO4_AC)/)/3.06 + (/(f1_II[:]->cb_SO4_NA)/)/3.06 + (/(f1_II[:]->cb_SO4_PR)/)/3.06
 load_II=load_II*1.e3 
 ext_II=aod_II/load_II
 else if (plot_type.eq.2) then
 var="EXT_BC"                              ; name of input-variable and plot
 varname="BC extinction coefficient"       ; variable name used in text string:
 aod_I=(/f1_I[:]->D550_BC/)                                  ; variable to be plotted from I
 aod_II=(/f1_II[:]->D550_BC/)                                ; variable to be plotted from II 
 if(ModI.eq."CAM4-Oslo") then
   ext_I=(/f1_I[:]->D550_BC/)/(/(f1_I[:]->C_BC)/)*1.e-3        ; variable to be plotted from I
   load_I=(/f1_I[:]->C_BC/)*1.e3                               ; variable to be plotted from I
 else
   ext_I=(/f1_I[:]->D550_BC/)/(/(f1_I[:]->cb_BC)/)*1.e-3     ; variable to be plotted from I 
   load_I=(/f1_I[:]->cb_BC/)*1.e3                             ; variable to be plotted from I 
 end if
 ext_II=(/f1_II[:]->D550_BC/)/(/(f1_II[:]->cb_BC)/)*1.e-3     ; variable to be plotted from II 
 load_II=(/f1_II[:]->cb_BC/)*1.e3                             ; variable to be plotted from II 
 else if (plot_type.eq.3) then
 var="EXT_POM"                             ; name of input-variable and plot
 varname="POM extinction coefficient"      ; variable name used in text string:
 aod_I=(/f1_I[:]->D550_POM/)                                 ; variable to be plotted from I
 aod_II=(/f1_II[:]->D550_POM/)                               ; variable to be plotted from II 
 if(ModI.eq."CAM4-Oslo") then
   ext_I=(/f1_I[:]->D550_POM/)/(/(f1_I[:]->C_POM)/)*1.e-3      ; variable to be plotted from I
   load_I=(/f1_I[:]->C_POM/)*1.e3                              ; variable to be plotted from I
 else
   ext_I=(/f1_I[:]->D550_POM/)/(/(f1_I[:]->cb_OM)/)*1.e-3   ; variable to be plotted from I 
   load_I=(/f1_I[:]->cb_OM/)*1.e3                            ; variable to be plotted from I 
 end if
 ext_II=(/f1_II[:]->D550_POM/)/(/(f1_II[:]->cb_OM)/)*1.e-3   ; variable to be plotted from II 
 load_II=(/f1_II[:]->cb_OM/)*1.e3                            ; variable to be plotted from II 
 else if (plot_type.eq.4) then
 var="EXT_SS"                              ; name of input-variable and plot
 varname="Sea-salt extinction coefficient" ; variable name used in text string:
 aod_I=(/f1_I[:]->D550_SS/)                                  ; variable to be plotted from I
 aod_II=(/f1_II[:]->D550_SS/)                                ; variable to be plotted from II 
 if(ModI.eq."CAM4-Oslo") then
   ext_I=(/f1_I[:]->D550_SS/)/(/(f1_I[:]->C_SS)/)*1.e-3        ; variable to be plotted from I
   load_I=(/f1_I[:]->C_SS/)*1.e3                               ; variable to be plotted from I
 else
   ext_I=(/f1_I[:]->D550_SS/)/(/(f1_I[:]->cb_SALT)/)*1.e-3     ; variable to be plotted from I 
   load_I=(/f1_I[:]->cb_SALT/)*1.e3                             ; variable to be plotted from I
 end if
 ext_II=(/f1_II[:]->D550_SS/)/(/(f1_II[:]->cb_SALT)/)*1.e-3     ; variable to be plotted from II 
 load_II=(/f1_II[:]->cb_SALT/)*1.e3                             ; variable to be plotted from II 
 else if (plot_type.eq.5) then
 var="EXT_DUST"                            ; name of input-variable and plot
 varname="Dust extinction coefficient"     ; variable name used in text string:
 aod_I=(/f1_I[:]->D550_DU/)                                  ; variable to be plotted from I
 aod_II=(/f1_II[:]->D550_DU/)                                ; variable to be plotted from II 
 if(ModI.eq."CAM4-Oslo") then
   ext_I=(/f1_I[:]->D550_DU/)/(/(f1_I[:]->C_DUST)/)*1.e-3      ; variable to be plotted from I
   load_I=(/f1_I[:]->C_DUST/)*1.e3                             ; variable to be plotted from I
 else
   ext_I=(/f1_I[:]->D550_DU/)/(/(f1_I[:]->cb_DUST)/)*1.e-3   ; variable to be plotted from I 
   load_I=(/f1_I[:]->cb_DUST/)*1.e3                           ; variable to be plotted from I 
 end if
 ext_II=(/f1_II[:]->D550_DU/)/(/(f1_II[:]->cb_DUST)/)*1.e-3   ; variable to be plotted from II 
 load_II=(/f1_II[:]->cb_DUST/)*1.e3                           ; variable to be plotted from II 
 end if
 end if
 end if
 end if
 end if

; Calculating area weighted extinctions

 ext_Ia=ext_I              ; initialization of global average variables
 ext_IIa=ext_II
 aod_Ia=ext_I 
 aod_IIa=ext_II
 load_Ia=ext_I
 load_IIa=ext_II

 xdims_I = dimsizes(gw0_I)
 ;print(xdims_I)
 ydims_I = dimsizes(ext_Ia)
 ;print(ydims_I)
 do i=0,dimsizes(gw0_I)-1
 ext_Ia(:,i,:)=ext_I(:,i,:)*coffa*dlon_I*gw0_I(i)
 aod_Ia(:,i,:)=aod_I(:,i,:)*coffa*dlon_I*gw0_I(i)
 load_Ia(:,i,:)=load_I(:,i,:)*coffa*dlon_I*gw0_I(i)
 end do
 extave_I=sum(dim_avg_n(ext_Ia,0))/area1
 aodave_I=sum(dim_avg_n(aod_Ia,0))/area1
 loadave_I=sum(dim_avg_n(load_Ia,0))/area1  

 xdims_II = dimsizes(gw0_II)
 ;print(xdims_I)
 ydims_II = dimsizes(ext_IIa)
 ;print(ydims_II)
 do i=0,dimsizes(gw0_II)-1
 ext_IIa(:,i,:)=ext_II(:,i,:)*coffa*dlon_II*gw0_II(i)
 aod_IIa(:,i,:)=aod_II(:,i,:)*coffa*dlon_II*gw0_II(i)
 load_IIa(:,i,:)=load_II(:,i,:)*coffa*dlon_II*gw0_II(i)
 end do
 extave_II=sum(dim_avg_n(ext_IIa,0))/area1
 aodave_II=sum(dim_avg_n(aod_IIa,0))/area1
 loadave_II=sum(dim_avg_n(load_IIa,0))/area1  

; Defining color scale
  digg=(/1,2,3,4,5,7,10,15,20,30/)

;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    Make the Plot
;
;;;;;;;;;;;;;;;;;;;;;;;;;

  wks = gsn_open_wks(format,var)

  gsn_define_colormap(wks,"amwg_blueyellowred")
;  gsn_define_colormap(wks,"BlueDarkRed18")
;  gsn_define_colormap(wks,"precip2_15lev")
;  gsn_define_colormap(wks,"gui_default")
;  gsn_define_colormap(wks,"hotres")
  plot=new(2,graphic)
  
  res                   = True                  ; plot mods desired
  res@gsnSpreadColors   = False                 ; use full colormap
  res@mpFillOn          = False                 
  res@cnFillOn          = True                  ; color fill 
  res@cnLinesOn         = False                 ; no contour lines
  res@cnLineLabelsOn    = False                 
  res@gsnFrame          = False                 ; Do not draw plot
  res@gsnDraw           = False                 ; Do not advance frame
  res@lbLabelBarOn      = False 
  res@tmXBOn =False
  res@tmXTOn =False
  res@tmYLOn =False
  res@tmYROn =False
  res@cnMissingValFillPattern = 0
  res@cnMissingValFillColor = 16
  res@tiMainFontHeightF = 0.03
  res@tiMainFontThicknessF = 2
  res@txFontHeightF       = 0.02
  res@cnFillMode        = "RasterFill"  ; Turn on raster fill
  res@tiMainFont = "helvetica"
  res@tmYRMode             = "Automatic"
  res@cnInfoLabelOn        = False
  res@cnLevelSelectionMode      = "ExplicitLevels"       ; manual levels
;  res@cnFillColors   = (/3,4,5,6,7,8,9,0,10,11,12,13,14,15,16/)   ; gir hvitt midt i ?
;  res@cnFillColors   = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/) 
  res@cnFillColors   = (/3,5,6,8,9,10,11,12,13,14,15,16/) 
;  res@cnLevels       =   sprintf("%4.1f",digg)             ; min level
  res@cnLevels       =   sprintf("%5.3f",digg)             ; min level

;  res@tiMainString   = "CAM4-Oslo"
if (var .eq. "EXT_SO4") then
  res@gsnRightString   = "avg = "+sprintf("%5.2f",extave_I)+" ("+sprintf("%5.2f",aodave_I/loadave_I)+") m~S~2~N~ (g S)~S~-1~N~"
else
  res@gsnRightString   = "avg = "+sprintf("%5.2f",extave_I)+" ("+sprintf("%4.2f",aodave_I/loadave_I)+") m~S~2~N~ g~S~-1~N~" 
end if
  res@gsnLeftString   = varname
  plot(0) = gsn_csm_contour_map_ce(wks,dim_avg_n_Wrap(ext_I,0),res)    ; create the plot

;  res@tiMainString   = "CAM5-Oslo"
if (var .eq. "EXT_SO4") then
;old version:  res@gsnRightString  = "avg = "+sprintf("%5.2f",extave_II)+" ("+sprintf("%5.2f",aodave_II/loadave_II)+") m~S~2~N~ (g S)~S~-1~N~"
;err  res@gsnRightString  = "avg = "+sprintf("%5.2f",extave_II)+" ("+sprintf("%5.2f",aodave_II/loadave_II)+") m~S~2~N~ (g sulfate)~S~-1~N~"
  res@gsnRightString  = "avg = "+sprintf("%5.2f",extave_II)+" ("+sprintf("%5.2f",aodave_II/loadave_II)+") m~S~2~N~ (g S)~S~-1~N~"
else
  res@gsnRightString  = "avg = "+sprintf("%5.2f",extave_II)+" ("+sprintf("%4.2f",aodave_II/loadave_II)+") m~S~2~N~ g~S~-1~N~" 
end if
  res@gsnLeftString   = varname
  plot(1) = gsn_csm_contour_map_ce(wks,dim_avg_n_Wrap(ext_II,0),res)   ; create the plot

  pres                            = True        ; panel plot mods desired
;  pres@gsnMaximize                = True        ; fill the page
;  pres@txString   =   var
  pres@txFontHeightF              =0.015
  pres@txFontThicknessF              =2
  pres@gsnPanelLabelBar = True
  pres@lbLabelFontHeightF = 0.01
;  pres@lbOrientation        ="Vertical"

  gsn_panel(wks,plot,(/1,2/),pres)              ; create panel plot

end
