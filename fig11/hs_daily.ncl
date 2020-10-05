
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   

begin

  a1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h1.hs_d.00010101-00501231.nc","r")
  a2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h1.hs_d.00010101-00511231.nc","r")
  b1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h1.aice_d.00010101-00501231.nc","r")
  b2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h1.aice_d.00010101-00511231.nc","r")

  hsai1tmp = a1->hs_d
  hsai2tmp = a2->hs_d
  aice1 = b1->aice_d
  aice2 = b2->aice_d

  tareatmp = a1->tarea
  tlat = a1->TLAT
  tlon = a1->TLON

  ndims = dimsizes(hsai1tmp)
  ntimes = ndims(0)

  hs1 = new((/ntimes/),"float")
  hs2 = new((/ntimes/),"float")


  ice_mask = where(tlat.ge.70. .and. tlat.le.85. .and. tlon .ge. 180. .and. tlon .le. 230., 1., 0.) 
; ice_mask = where(tlat.ge.0.,1.,0.)
  
  do n=0,ntimes-1
     aice1tmp = where(aice1(n,:,:).ne.0.,aice1(n,:,:),aice1@_FillValue)
     aice2tmp = where(aice2(n,:,:).ne.0.,aice2(n,:,:),aice2@_FillValue)
     hs1tmp = hsai1tmp(n,:,:)/aice1tmp(:,:)
     hs2tmp = hsai2tmp(n,:,:)/aice2tmp(:,:)
     tarea = where(aice1(n,:,:).gt.0.15 .and. aice2(n,:,:).gt.0.15,tareatmp,0.)
     tot_tarea = dim_sum(dim_sum(tarea*ice_mask))
     hs1(n) = dim_sum(dim_sum(hs1tmp(:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
     hs2(n) = dim_sum(dim_sum(hs2tmp(:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
  end do

  hs1a = reshape(hs1,(/50,365/))
  hs2a = reshape(hs2,(/50,365/))

  hs1clim = dim_avg_n(hs1a,0)
  hs2clim = dim_avg_n(hs2a,0)
  hs1std = dim_stddev_n(hs1a,0)
  hs2std = dim_stddev_n(hs2a,0)

  wks2 = gsn_open_wks("pdf", "hs_beaufort")    ; Open PS file "contour3a.ps"

  res = True
  res@gsnDraw = True
  res@gsnFrame = True
  res@gsnMaximize         = True              ; maximize plot area
  res@trXMinF = 1
  res@trXMaxF = 365
  res@trYMinF = 0.0
  res@trYMaxF = 0.3
  res@xyLineThicknesses  = (/3,3,3,3,3,3,3/)
  res@xyLineColors       = (/"blue","red","lightblue","lightblue","pink","pink","black"/)
  res@xyDashPatterns     = (/0,0,1,1,1,1,0/)
; res@tiMainString       = "BL99 (solid) MUSHY (dashed)"
  res@tiXAxisString = "Day"
  res@tiYAxisString = "Snow Thickness"


  x = fspan(1,365,365)
  hs = (/hs1clim,hs2clim,hs1clim-hs1std,hs1clim+hs1std,hs2clim-hs2std,hs2clim+hs2std,abs(hs2clim-hs1clim)/)
  plot1 = gsn_csm_xy(wks2,x,hs,res)

end

