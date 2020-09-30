
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   

begin

  a1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h.apeff_ai.000101-005012.nc","r")
  a2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h.apeff_ai.000101-005112.nc","r")
  b1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h.aice.000101-005012.nc","r")
  b2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h.aice.000101-005112.nc","r")

  apeffai1tmp = a1->apeff_ai
  apeffai2tmp = a2->apeff_ai
  aice1 = b1->aice
  aice2 = b2->aice

  tareatmp = a1->tarea
  tlat = a1->TLAT
  tlon = a1->TLON

  ndims = dimsizes(apeffai1tmp)
  ntimes = ndims(0)

  apeff1 = new((/ntimes/),"float")
  apeff2 = new((/ntimes/),"float")


  ice_mask = where(tlat.ge.70. .and. tlat.le.85. .and. tlon .ge. 180. .and. tlon .le. 230., 1., 0.) 
  
  do n=0,ntimes-1
     aice1tmp = where(aice1(n,:,:).ne.0.,aice1(n,:,:),aice1@_FillValue)
     aice2tmp = where(aice2(n,:,:).ne.0.,aice2(n,:,:),aice2@_FillValue)
     apeff1tmp = apeffai1tmp(n,:,:)/aice1tmp(:,:)
     apeff2tmp = apeffai2tmp(n,:,:)/aice2tmp(:,:)
     tarea = where(aice1(n,:,:).gt.0.15 .and. aice2(n,:,:).gt.0.15,tareatmp,0.)
     tot_tarea = dim_sum(dim_sum(tarea*ice_mask))
     apeff1(n) = dim_sum(dim_sum(apeff1tmp(:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
     apeff2(n) = dim_sum(dim_sum(apeff2tmp(:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
  end do

  apeff1a = reshape(apeff1,(/50,12/))
  apeff2a = reshape(apeff2,(/50,12/))

  apeff1clim = dim_avg_n(apeff1a,0)
  apeff2clim = dim_avg_n(apeff2a,0)
  apeff1std = dim_stddev_n(apeff1a,0)
  apeff2std = dim_stddev_n(apeff2a,0)

  wks2 = gsn_open_wks("pdf", "apeff_nh")    ; Open PS file "contour3a.ps"

  res = True
  res@gsnDraw = True
  res@gsnFrame = True
  res@gsnMaximize = True
  res@trXMinF = 1
  res@trXMaxF = 12
  res@trYMinF = -0.1
  res@trYMaxF = 1.1
  res@xyLineThicknesses  = (/3,3,3,3,3,3/)
  res@xyLineColors       = (/"blue","red","lightblue","lightblue","pink","pink"/)
  res@xyDashPatterns     = (/0,0,1,1,1,1/)
; res@tiMainString       = "BL99 (solid) MUSHY (dashed)"
  res@tiXAxisString = "Month"
  res@tiYAxisString = "Pond Fraction"


  x = fspan(1,12,12)
  apeff = (/apeff1clim,apeff2clim,apeff1clim-apeff1std,apeff1clim+apeff1std,apeff2clim-apeff2std,apeff2clim+apeff2std/)
  plot1 = gsn_csm_xy(wks2,x,apeff,res)

end

