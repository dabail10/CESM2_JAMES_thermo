load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "popRemap.ncl"
;**************************************************
; contributed.ncl must be loaded before popRemap.ncl
;**************************************************
begin

; This plots 2 scalar fields and their difference from 2
; cases of CCSM output.  Differences are of the form
; new - prev

wkdir = "./"
case_new = "b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1"
case_prev  = "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1"
path_prev = "/glade/p/cesm/pcwg/dbailey/archive/" + case_prev + "/atm/proc/tseries/month_1/"
path_new = "/glade/p/cesm/pcwg/dbailey/archive/" + case_new + "/atm/proc/tseries/month_1/"
dir = path_prev

;---------USER DEFINED QUANTITIES------------------
c_hemisphere = (/"NH", "SH"/)
min_lat   = (/55., -90./) 
max_lat   = (/90., -55./)

var_to_plot = (/"TREFHT"/)

nvar = dimsizes(var_to_plot)

ncntrs = 14
NH_Pcntrs = new((/nvar,ncntrs/), float)
SH_Pcntrs = new((/nvar,ncntrs/), float)

; TREFHT
NH_Pcntrs(0,:)=(/230.,235.,240.,245.,250.,255.,260.,265.,270.,275.,280.,285.,290.,295./)
SH_Pcntrs(0,:)=(/230.,235.,240.,245.,250.,255.,260.,265.,270.,275.,280.,285.,290.,295./)

;-------END USER DEFINED QUANTITIES----------------

 do k = 0,nvar-1                      ; Loop through variables to be plotted
  ;--------------------------------------------------------------------
  ; Construct file name
  ;--------------------------------------------------------------------
    new_file_name1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cam.h0.TREFHT.000101-005112.nc"
    prev_file_name1= "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cam.h0.TREFHT.000101-005012.nc"

  ;--------------------------------------------------------------------
  ; Read variables into memory
  ;--------------------------------------------------------------------
    a1 = addfile(path_prev + prev_file_name1,"r")     ; Read in both cases
    b1 = addfile(path_new  + new_file_name1 ,"r")

    var_new1 = b1->TREFHT
    var_prev1 = a1->TREFHT

    lat = b1->lat
    lon = b1->lon

    var_new_ann = month_to_annual(var_new1,1)
    var_new_jas = month_to_season(var_new1,"JAS")
    var_prev_ann = month_to_annual(var_prev1,1)
    var_prev_jas = month_to_season(var_prev1,"JAS")

    var_new = dim_avg_n(var_new_ann, 0)
    var_prev = dim_avg_n(var_prev_ann, 0)
    var_new_var = dim_variance_n(var_new_ann, 0)
    var_prev_var = dim_variance_n(var_prev_ann, 0)

    prob = ttest(var_new,var_new_var,50,var_prev,var_prev_var,50,False,False)
    significance = where(prob.lt.0.05,1.,0.)
    var_diff = significance*(var_new-var_prev)


    var_new!0 = "lat"
    var_new!1 = "lon"
    var_new&lat = lat
    var_new&lon = lon
    var_prev!0 = "lat"
    var_prev!1 = "lon"
    var_prev&lat = lat
    var_prev&lon = lon
    var_diff!0 = "lat"
    var_diff!1 = "lon"
    var_diff&lat = lat
    var_diff&lon = lon

;--------------------------------------------------------------------
; Loop through the hemispheres to create plots
;--------------------------------------------------------------------

    do ihem = 0,1         ; Loop through hemispheres

      wks = gsn_open_wks("pdf", wkdir + "diff_con_" + \
                        var_to_plot(k) + "_cice_"+c_hemisphere(ihem))

      plot = new(3,graphic)  

      cmap = read_colormap_file("MPL_RdBu")
      cmap = cmap(::-1,:)

      res                     = True
      res@gsnSpreadColors     = True              ; use full colormap
      res@cnFillPalette = cmap
      res@gsnSpreadColorStart  = 64 
      res@gsnSpreadColorEnd    = 127

      res@gsnDraw              = False            ; do not draw picture
      res@gsnFrame             = False            ; do not advance frame
      res@gsnPolar             = c_hemisphere(ihem)
      res@gsnTickMarksOn       = False            ; don't write long. labels

      res@cnFillOn             = True 
      res@cnLinesOn            = True
      res@cnLineLabelInterval  = 0                ; More labels on lines
      res@mpMinLatF            = min_lat(ihem)    ; specify min lat
      res@mpMaxLatF            = max_lat(ihem)    ; specify max lat
;     res@mpLandFillColor      = "grey"           ; choose color of continents.
;     res@mpFillDrawOrder      = "PostDraw"       ; Fill continents later to cover up lines
;     res@mpPerimDrawOrder     = "PostDraw"       ; Perimeter gets covered up by land
      res@lbLabelFontHeightF   = 1.3
      res@lbOrientation        = "Vertical"
  
      res@txFontHeightF        = 0.02

      res@cnLevelSelectionMode = "ExplicitLevels" ; set manual contour levels
; Set contour levels for NH
      if (c_hemisphere(ihem).eq."NH") then
        res@cnLevels   =  NH_Pcntrs(k,:)
      end if
; Set contour levels for SH
      if (c_hemisphere(ihem).eq."SH") then
        res@cnLevels   =  SH_Pcntrs(k,:)
      end if

      res@gsnLeftString       = " "
      res@gsnRightString      = " "
      res@gsnAddCyclic        = True

      res@txFontHeightF        = 0.025

      res@lbLabelBarOn         = False            ; Turns off label bar for first plot
      res@lbLabelStride        = 1                ; Label bar stride
      res@lbLabelFontHeightF   = 0.02             ; Label bar font height

      res@tiMainFontHeightF    = 0.03             ; size of title main string

      plot(0) = gsn_csm_contour_map_polar (wks,var_new(:,:), res)

      res@lbLabelBarOn         = True            ; Turn on label bar for second plot
      plot(1) = gsn_csm_contour_map_polar (wks,var_prev(:,:), res)

;---------------------------------------------------
; Define new scale and colormap for difference plots
;---------------------------------------------------

      delete(res@cnLevels)
      delete(res@cnLevelSelectionMode)
      delete(res@cnFillPalette)
      delete(cmap)

      cmap = read_colormap_file("MPL_RdBu")
      cmap = cmap(::-1,:)

      res                     = True
      res@gsnSpreadColors     = True              ; use full colormap
      res@cnLinesOn           = True
      res@cnFillPalette = cmap
      res@gsnSpreadColorStart  = 0
      res@gsnSpreadColorEnd    = 127

      stddev_got_ice = stddev(var_diff(:,:))

      minXY   = -4.*stddev_got_ice
      maxXY   =  4.*stddev_got_ice
      mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
      res@cnMinLevelValF       = mnmxint(0)
      res@cnMaxLevelValF       = mnmxint(1)
      res@cnLevelSpacingF      = mnmxint(2) * 0.5

      plot(2) = gsn_csm_contour_map_polar (wks,var_diff(:,:), res)

; Set panel parameters; draw panel without white space

      resP = True
;     resP@gsnPaperOrientation = "portrait"  
      resP@gsnPaperOrientation = "landscape"  
      resP@gsnMaximize         = True              ; maximize plot area
      resP@txFontHeightF       = 0.02              ; Font height of top label
      resP@txFontColor         = "black"           ; Font color of top label
;     resP@txString            = main_labels(jj)

;     resP@gsnPanelYWhiteSpacePercent = 5
      resP@gsnPanelRight = 0.90
      gsn_panel(wks,plot,(/2,2/),resP)

      delete(res@cnFillPalette)
      delete(cmap)
  end do    ; Loop through hemispheres

  delete(var_new)
  delete(var_prev)
  delete(var_diff)

end do   ; End loop through variables

end
