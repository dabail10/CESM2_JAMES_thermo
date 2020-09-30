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
case_new = "b.e21.B1850.f09_g17.CMIP6-piControl.001_snowmin"
case_prev  = "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snowmin"
path_prev = "/glade/p/cesm/pcwg/dbailey/archive/" + case_prev + "/ice/proc/tseries/month_1/"
path_new = "/glade/p/cesm/pcwg/dbailey/archive/" + case_new + "/ice/proc/tseries/month_1/"
dir = path_prev

;---------USER DEFINED QUANTITIES------------------
c_hemisphere = (/"NH", "SH"/)
min_lat   = (/55., -90./) 
max_lat   = (/90., -55./)

var_to_plot = (/"hi"/)

nvar = dimsizes(var_to_plot)

ncntrs = 10
NH_Pcntrs = new((/nvar,ncntrs/), float)
SH_Pcntrs = new((/nvar,ncntrs/), float)

; hi contour levels
NH_Pcntrs(0,:)=fspan(0.4,4.0,10)
SH_Pcntrs(0,:)=fspan(0.2,2,10)

;-------END USER DEFINED QUANTITIES----------------

 do k = 0,nvar-1                      ; Loop through variables to be plotted
  ;--------------------------------------------------------------------
  ; Construct file name
  ;--------------------------------------------------------------------
    new_file_name1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001_snowmin.cice.h.hi.000101-005112.nc"
    prev_file_name1= "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snowmin.cice.h.hi.000101-005012.nc"

  ;--------------------------------------------------------------------
  ; Read variables into memory
  ;--------------------------------------------------------------------
    a1 = addfile(path_prev + prev_file_name1,"r")     ; Read in both cases
    b1 = addfile(path_new  + new_file_name1 ,"r")
    grid = addfile("/glade/p/cesm/omwg/grids/gx1v7_grid.nc","r")

    var_new1 = b1->hi
    var_prev1 = a1->hi

    lat2d = grid->TLAT
    lon2d = grid->TLONG

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

    var_new@lat2d = lat2d
    var_new@lon2d = lon2d
    var_prev@lat2d = lat2d
    var_prev@lon2d = lon2d
    var_diff@lat2d = lat2d
    var_diff@lon2d = lon2d

;--------------------------------------------------------------------
; Loop through the hemispheres to create plots
;--------------------------------------------------------------------

    do ihem = 0,1         ; Loop through hemispheres

      wks = gsn_open_wks("pdf", wkdir + "diff_con_" + \
                        var_to_plot(k) + "_cice_"+c_hemisphere(ihem))

      plot = new(3,graphic)  

      cmap = read_colormap_file("cmocean_ice")

      res                      = True
      res@gsnSpreadColors     = True              ; use full colormap
      res@cnFillPalette = cmap
      res@gsnSpreadColorStart  = 0
      res@gsnSpreadColorEnd    = 255
      res@gsnDraw              = False            ; do not draw picture
      res@gsnFrame             = False            ; do not advance frame
      res@gsnPolar             = c_hemisphere(ihem)
      res@gsnTickMarksOn       = False            ; don't write long. labels

      res@cnFillOn             = True 
      res@cnLinesOn            = True            
      res@cnLineLabelInterval  = 0                ; More labels on lines
      res@cnLevelSelectionMode = "ExplicitLevels" ; set manual contour levels
      res@mpMinLatF            = min_lat(ihem)    ; specify min lat
      res@mpMaxLatF            = max_lat(ihem)    ; specify max lat
      res@mpLandFillColor      = "grey"           ; choose color of continents.
      res@mpFillDrawOrder      = "PostDraw"       ; Fill continents later to cover up lines
      res@mpPerimDrawOrder     = "PostDraw"       ; Perimeter gets covered up by land
      res@lbLabelFontHeightF   = 1.3
      res@lbOrientation        = "Vertical"
      res@trGridType           = "TriangularMesh"
  
      res@txFontHeightF        = 0.02
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
      delete(res@cnFillPalette)
      delete(res@cnLevelSelectionMode)
      delete(cmap)

      cmap = read_colormap_file("MPL_RdBu")
      cmap = cmap(::-1,:)
;     n_colors_one = 17
;     n_colors_two = 101

      stddev_got_ice = stddev(var_diff(:,:))

      minXY   = -3.*stddev_got_ice
      maxXY   =  3.*stddev_got_ice
      mnmxint = nice_mnmxintvl( minXY, maxXY, ncntrs, False)

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
      res@cnMinLevelValF       = mnmxint(0)
      res@cnMaxLevelValF       = mnmxint(1)
      res@cnLevelSpacingF      = mnmxint(2) * 0.5

;     color_start              = n_colors_one + 1
;     color_end                = n_colors_one+n_colors_two-1
      res@cnFillPalette = cmap
      res@gsnSpreadColors     = True              ; use full colormap
      res@gsnSpreadColorStart  = 0
      res@gsnSpreadColorEnd    = 127
      res@gsnAddCyclic        = True
      res@cnLinesOn           = True

      plot(2) = gsn_csm_contour_map_polar (wks,var_diff(:,:), res)
      res2 = True
      res2@gsnDraw              = False            ; do not draw picture
      res2@gsnFrame             = False            ; do not advance frame
      res2@cnLevelSelectionMode =  "ExplicitLevels"
      res2@cnInfoLabelOn    =  False
      res2@gsnLeftString    =  ""
      res2@cnLevels         =  (/0.99/)
      res2@tiXAxisString    = ""
      res2@tiYAxisString    = ""
      res2@cnLineThicknessF = 3.0      ; thicker contours
      res2@cnLineLabelsOn   = False     ; no line labels
      ploto = gsn_csm_contour(wks,significance, res2)
      overlay(plot(2),ploto)


; Set panel parameters; draw panel without white space

      resP = True
      resP@gsnFrame             = False            ; do not advance frame
;     resP@gsnPaperOrientation = "portrait"  
      resP@gsnPaperOrientation = "landscape"  
      resP@gsnMaximize         = True              ; maximize plot area
      resP@txFontHeightF       = 0.02              ; Font height of top label
      resP@txFontColor         = "black"           ; Font color of top label
;     resP@txString            = main_labels(jj)

;     resP@gsnPanelYWhiteSpacePercent = 5
      resP@gsnPanelRight = 0.90
      gsn_panel(wks,plot,(/2,2/),resP)

      frame(wks)

      delete(res@cnFillPalette)
      delete(cmap)
  end do    ; Loop through hemispheres

  delete(var_new)
  delete(var_prev)
  delete(var_diff)

end do   ; End loop through variables

end