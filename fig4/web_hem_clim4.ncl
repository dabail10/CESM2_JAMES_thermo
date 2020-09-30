
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "readobs.ncl"

begin


ncases1 = 2
case1a = "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snowmin"
case1b = "b.e21.B1850.f09_g17.CMIP6-piControl.001_snowmin"

f1a = addfile("ice_vol_"+case1a+"_1-50.nc","r")
f1b = addfile("ice_vol_"+case1b+"_1-50.nc","r")

ext_mo_nh_tmp = f1a->ext_mo_nh
ndim = dimsizes(ext_mo_nh_tmp)
ntimes = ndim(0)

chem = (/"NH","SH"/)

ext_mo_nh_in = new((/ntimes,ncases1/),float)
ext_mo_sh_in = new((/ntimes,ncases1/),float)

ext_mo_nh_in(:,0) = f1a->ext_mo_nh
ext_mo_sh_in(:,0) = f1a->ext_mo_sh
ext_mo_nh_in(:,1) = f1b->ext_mo_nh
ext_mo_sh_in(:,1) = f1b->ext_mo_sh

ext_mo_nh = new((/ntimes/12,12,ncases1/),float)
ext_mo_sh = new((/ntimes/12,12,ncases1/),float)

do n=0,ncases1-1
   ext_mo_nh(:,:,n) = onedtond(ext_mo_nh_in(:,n),(/ntimes/12,12/))
   ext_mo_sh(:,:,n) = onedtond(ext_mo_sh_in(:,n),(/ntimes/12,12/))
end do

ext_mo_nh!0 = "years"
ext_mo_nh!1 = "months"
ext_mo_sh!0 = "years"
ext_mo_sh!1 = "months"

ext_data = new((/2,12,ncases1/),float)
ext_data(0,:,0:ncases1-1) = dim_avg_n(ext_mo_nh,0)
ext_data(1,:,0:ncases1-1) = dim_avg_n(ext_mo_sh,0)

ext = new((/2,37,12/),float)
readobs(ext)

ext@long_name = "NSIDC sea ice extent"
ext@_FillValue = -9999

ext_clim = dim_avg_n(ext,1)

res = True
res@txFont   = "helvetica"
res@txFuncCode = "~"
res@tiMainFont          = "helvetica"
res@tiMainFuncCode          = "~"
res@tiXAxisFont          = "helvetica"
res@tiXAxisFuncCode          = "~"
res@tiYAxisFont          = "helvetica"
res@tiYAxisFuncCode          = "~"
res@tmXBLabelFont          = "helvetica"
res@tmYLLabelFont          = "helvetica"

res@trXMinF = 1.
res@trXMaxF = 12.

res@trYMinF = 0.

x = fspan(1,12,12)
x!0 = "month"
x@long_name = "Month"

do ihem=0,1

   y = new((/ncases1+1,12/),float)
   do n=0,ncases1-1
      y(n,:) = ext_data(ihem,:,n)
   end do
   y(ncases1,:) = ext_clim(ihem,:)

   res@xyExplicitLegendLabels = (/"BL99","MUSHY","NSIDC"/)
   res@xyLineColors = (/"blue","red","black"/)
   res@xyDashPatterns    = (/3,2,0/)
   res@xyLineThicknesses = (/3,3,3/)

   plotname = "clim_mo_"+chem(ihem)+"_diff"
   wks1 = gsn_open_wks("pdf",plotname)

   if (ihem .eq. 0) then
      y@long_name = "NH Sea Ice Extent 10~S~12~N~ m~S~2~N~"
   else
      y@long_name = "SH Sea Ice Extent 10~S~12~N~ m~S~2~N~"
   end if

;
; Add a legend.
;
   res@trYMaxF = 18.+int2flt(ihem)*5.
   res@pmLegendDisplayMode    = "Always"                 ; Turn on legend
   res@pmLegendWidthF         = 0.325                     ; Resize legend width
   res@pmLegendHeightF        = 0.15                     ; and height
;  res@pmLegendOrthogonalPosF = -1.15+int2flt(ihem)*0.65 ; Move legend up
   res@pmLegendOrthogonalPosF = -1.15 ; Move legend up
   res@pmLegendParallelPosF   = 0.6-int2flt(ihem)*0.3   ; and to the right.
   res@lgAutoManage = False
   res@lgBoxMinorExtentF = 0.1
   res@lgLabelFontHeightF = 0.015
   res@lgPerimOn = False

   plot1 = gsn_csm_xy(wks1,x,y,res)

   delete(wks1)
   delete(y)
   delete(res@xyExplicitLegendLabels)
   delete(res@xyLineColors)
   delete(res@xyLineThicknesses)
   delete(res@xyDashPatterns)

end do ; ihem

end
