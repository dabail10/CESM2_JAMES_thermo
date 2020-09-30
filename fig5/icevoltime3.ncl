
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

wks_type = "pdf"
wks_type@wkOrientation = "landscape"


;wks   = gsn_open_wks (wks_type,"ice_vol_sh")            ; open workstation
;yaxis = "SH Ice Volume 10~S~13 ~N~m~S~3" 
wks   = gsn_open_wks (wks_type,"ice_area_nh")            ; open workstation
yaxis = "NH Ice Area 10~S~12 ~N~m~S~2" 

f1    = addfile ("ice_vol_b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1_1-50.nc","r")                 ; note the "s" of addfile
f2    = addfile ("ice_vol_b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1_1-50.nc","r")                 ; note the "s" of addfile

ice_vol1  = f1->vhi_mo_nh
ice_area1 = f1->vai_mo_nh
ice_vol2  = f2->vhi_mo_nh
ice_area2 = f2->vai_mo_nh

ice_vol1_ann = dim_avg(reshape(ice_vol1,(/50,12/)))
ice_vol2_ann = dim_avg(reshape(ice_vol2,(/50,12/)))
ice_area1_ann = dim_avg(reshape(ice_area1,(/50,12/)))
ice_area2_ann = dim_avg(reshape(ice_area2,(/50,12/)))

ice_vol1_mean = dim_avg(ice_vol1_ann)
ice_vol2_mean = dim_avg(ice_vol2_ann)
ice_vol1_var = variance(ice_vol1_ann)
ice_vol2_var = variance(ice_vol2_ann)

ice_area1_mean = dim_avg(ice_area1_ann)
ice_area2_mean = dim_avg(ice_area2_ann)
ice_area1_var = variance(ice_area1_ann)
ice_area2_var = variance(ice_area2_ann)

tstat = ttest(ice_vol1_mean,ice_vol1_var,50,ice_vol2_mean,ice_vol2_var,50,True,True)
fstat = ftest(ice_vol1_var,50,ice_vol2_var,50,0)
tstat2 = ttest(ice_area1_mean,ice_area1_var,50,ice_area2_mean,ice_area2_var,50,True,True)
fstat2 = ftest(ice_area1_var,50,ice_area2_var,50,0)

print(ice_vol1_mean)
print(ice_vol2_mean)
print(ice_vol1_var)
print(ice_vol2_var)
print(tstat)
print(fstat)

print(ice_area1_mean)
print(ice_area2_mean)
print(ice_area1_var)
print(ice_area2_var)
print(tstat2)
print(fstat2)

res              = True                     ; plot mods desired

res@tmXMajorGrid = True
;res@tmYMajorGrid = True
res@tmXMajorGridLineDashPattern = 2
;res@tmYMajorGridLineDashPattern = 2

res@trXMaxF = 50.
res@trXMinF = 1.

res@xyLineThicknesses = (/3,3/)
res@xyLineColors = (/"blue","red"/)
res@xyDashPattern = (/0,0/) 

;res@tiMainString = title          ; add title
res@tiXAxisString = "Year"
res@tiYAxisString = yaxis
res@gsnMaximize  = True
res@vpHeightF= 0.3                    ; change aspect ratio of plot
res@vpWidthF = 0.8

;res@trYMaxF = 3.
;res@trYMinF = 0.

x = fspan(1,50,50)

y = (/ice_vol1_ann,ice_vol2_ann/)
;y = (/ice_area1_ann,ice_area2_ann/)

plot  = gsn_csm_xy(wks,x,y,res)

end
