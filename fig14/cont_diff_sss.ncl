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
path_prev = "/glade/p/cesm/pcwg/dbailey/archive/" + case_prev + "/ocn/proc/tseries/month_1/"
path_new = "/glade/p/cesm/pcwg/dbailey/archive/" + case_new + "/ocn/proc/tseries/month_1/"
dir = path_prev

;---------USER DEFINED QUANTITIES------------------
c_hemisphere = (/"NH", "SH"/)
min_lat   = (/55., -90./) 
max_lat   = (/90., -55./)

var_to_plot = (/"SALT"/)

nvar = dimsizes(var_to_plot)

ncntrs = 14
NH_Pcntrs = new((/nvar,ncntrs/), float)
SH_Pcntrs = new((/nvar,ncntrs/), float)

; SSS contour levels
NH_Pcntrs(0,:)=(/0,4,8,12,16,18,20,22,24,26,28,30,32,34/) ; aice
SH_Pcntrs(0,:)=(/0,4,8,12,16,18,20,22,24,26,28,30,32,34/) ; aice

;-------END USER DEFINED QUANTITIES----------------

 do k = 0,nvar-1                      ; Loop through variables to be plotted
  ;--------------------------------------------------------------------
  ; Construct file name
  ;--------------------------------------------------------------------
    new_file_name1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.pop.h."+var_to_plot(k)+".000101-005112.nc"
    prev_file_name1= "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.pop.h."+var_to_plot(k)+".000101-005012.nc"

  ;--------------------------------------------------------------------
  ; Read variables into memory
  ;--------------------------------------------------------------------
    a1 = addfile(path_prev + prev_file_name1,"r")     ; Read in both cases
    b1 = addfile(path_new  + new_file_name1 ,"r")

    var_new1 = b1->$var_to_plot(k)$
    var_prev1 = a1->$var_to_plot(k)$

    lat2d = b1->TLAT
    lon2d = b1->TLONG

    var_new_ann = month_to_annual(var_new1(:,0,:,:),1)
    var_new_jas = month_to_season(var_new1(:,0,:,:),"JAS")
    var_prev_ann = month_to_annual(var_prev1(:,0,:,:),1)
    var_prev_jas = month_to_season(var_prev1(:,0,:,:),"JAS")

    var_new = dim_avg_n(var_new_jas(29:49,:,:), 0)*1000.
    var_prev = dim_avg_n(var_prev_jas(29:49,:,:), 0)*1000.
    var_new_var = dim_variance_n(var_new_jas(29:49,:,:), 0)*1000.
    var_prev_var = dim_variance_n(var_prev_jas(29:49,:,:), 0)*1000.

    prob = ttest(var_new,var_new_var,21,var_prev,var_prev_var,21,False,False)
;   significance = 100.*(1.-prob)
    significance = where(prob .le. 0.05, 1., 0.)
    var_diff = significance*(var_new-var_prev)

    print(var_new_jas(:,342,310))
    print(var_prev_jas(:,342,310))

    print(var_new_var(342,310))
    print(var_prev_var(342,310))
    print(prob(342,310))
    print(significance(342,310))

    var_new@lat2d = lat2d
    var_new@lon2d = lon2d
    var_prev@lat2d = lat2d
    var_prev@lon2d = lon2d
    var_diff@lat2d = lat2d
    var_diff@lon2d = lon2d
;   significance@lat2d = lat2d
;   significance@lon2d = lon2d
;   prob@lat2d = lat2d
;   prob@lon2d = lon2d

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
      res@trGridType           = "TriangularMesh"
      res@mpMinLatF            = min_lat(ihem)    ; specify min lat
      res@mpMaxLatF            = max_lat(ihem)    ; specify max lat
      res@mpLandFillColor      = "grey"           ; choose color of continents.
      res@mpFillDrawOrder      = "PostDraw"       ; Fill continents later to cover up lines
      res@mpPerimDrawOrder     = "PostDraw"       ; Perimeter gets covered up by land
      res@lbLabelFontHeightF   = 1.3
      res@lbOrientation        = "Vertical"
  
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

      delete(res@cnLevelSelectionMode)
      delete(res@cnFillPalette)
      delete(cmap)

      cmap = read_colormap_file("MPL_RdBu")
      cmap = cmap(::-1,:)

      res@cnFillPalette = cmap
      res@gsnSpreadColors     = True              ; use full colormap
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

      res@gsnAddCyclic        = True
      res@cnLinesOn           = True

      plot(2) = gsn_csm_contour_map_polar (wks,var_diff(:,:), res)
; set resources for 95% statistical significance overlay
;     res_95                            = True               ; reset options
;     res_95@gsnLeftString              = ""
;     res_95@gsnRightString             = ""
;     res_95@gsnFrame                   = False             ; Do not draw plot
;     res_95@gsnDraw                    = False             ; Do not advance fr
;     res_95@cnInfoLabelOn              = False             ; Turn off info on contours at bottom
;     res_95@lbLabelBarOn               = False             ; Turn off label bar
;     res_95@cnFillOn                   = True              ; turn on filled contours
;     res_95@cnMonoFillColor            = True
;     res_95@cnMonoFillPattern          = False
;     res_95@cnLevelSelectionMode       = "ManualLevels"    ; set manual contour
;     res_95@cnMinLevelValF             = 99.0              ; set min contour level specified above
;     res_95@cnMaxLevelValF             = 100.0              ; set max contour level specified above
;     res_95@cnLevelSpacingF            = 0.25             ; set contour spacing specified above
;     res_95@cnFillColors               = (/"white","gray","gray","gray"/)
;     res_95@cnFillPatterns             = (/-1,2,2,2/)      ; stipling where significant and above threshold
;     ploto = gsn_csm_contour(wks,significance,res_95)

;     overlay(plot(2),ploto)

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

      delete(res@cnLevelSelectionMode)
      delete(res@cnFillPalette)
      delete(cmap)

  end do    ; Loop through hemispheres

  delete(var_new)
  delete(var_prev)
  delete(var_diff)

end do   ; End loop through variables

end
