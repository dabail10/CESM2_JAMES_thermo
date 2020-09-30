
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   

begin

  a1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h1.fswdn_d.00010101-00501231.nc","r")
  a2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h1.fswdn_d.00010101-00511231.nc","r")
  b1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h1.fswabs_d.00010101-00501231.nc","r")
  b2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h1.fswabs_d.00010101-00511231.nc","r")
  c1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h1.aice_d.00010101-00501231.nc","r")
  c2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/day_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h1.aice_d.00010101-00511231.nc","r")

  fswdn1 = a1->fswdn_d
  fswdn2 = a2->fswdn_d
  fswabs1 = b1->fswabs_d
  fswabs2 = b2->fswabs_d
  aice1 = c1->aice_d
  aice2 = c2->aice_d
  tareatmp = a1->tarea
  tlat = a1->TLAT
  tlon = a1->TLON

  ndims = dimsizes(fswdn1)
  ntimes = ndims(0)

  alb1 = new((/ntimes/),"float")
  alb2 = new((/ntimes/),"float")


  ice_mask = where(tlat.ge.70. .and. tlat.le.85. .and. tlon .ge. 180. .and. tlon .le. 230., 1., 0.)

  do n=0,ntimes-1
     fswdn1tmp = where(fswdn1(n,:,:).ne.0.,fswdn1(n,:,:),fswdn1@_FillValue)
     fswdn2tmp = where(fswdn2(n,:,:).ne.0.,fswdn2(n,:,:),fswdn2@_FillValue)
     alb1tmp = (fswdn1(n,:,:)-aice1(n,:,:)*fswabs1(n,:,:))/fswdn1tmp
     alb2tmp = (fswdn2(n,:,:)-aice2(n,:,:)*fswabs2(n,:,:))/fswdn2tmp
     tarea = where(aice1(n,:,:).ge.0.15 .and. aice2(n,:,:).ge.0.15,tareatmp,0.)
     tot_tarea = dim_sum(dim_sum(tarea*ice_mask))
     alb1(n) = dim_sum(dim_sum(alb1tmp(:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
     alb2(n) = dim_sum(dim_sum(alb2tmp(:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
  end do

  alb1a = reshape(alb1,(/50,365/))
  alb2a = reshape(alb2,(/50,365/))

  alb1clim = dim_avg_n(alb1a,0)
  alb2clim = dim_avg_n(alb2a,0)
  alb1std = dim_stddev_n(alb1a,0)
  alb2std = dim_stddev_n(alb2a,0)
  alb1clim = where(alb1clim .gt. 0., alb1clim, fswdn1@_FillValue)
  alb2clim = where(alb2clim .gt. 0., alb2clim, fswdn2@_FillValue)
  alb1clim = where(alb1clim .le. 1., alb1clim, fswdn1@_FillValue)
  alb2clim = where(alb2clim .le. 1., alb2clim, fswdn2@_FillValue)
  alb1std = where(alb1std .gt. 0., alb1std, fswdn1@_FillValue)
  alb2std = where(alb2std .gt. 0., alb2std, fswdn2@_FillValue)
  alb1std = where(alb1std .le. 1., alb1std, fswdn1@_FillValue)
  alb2std = where(alb2std .le. 1., alb2std, fswdn2@_FillValue)

  alb1clim@_FillValue = fswdn1@_FillValue
  alb2clim@_FillValue = fswdn1@_FillValue
  alb1std@_FillValue = fswdn1@_FillValue
  alb2std@_FillValue = fswdn1@_FillValue

  wks2 = gsn_open_wks("pdf", "alb_nh")    ; Open PS file "contour3a.ps"

  res = True
  res@gsnDraw = True
  res@gsnFrame = True
  res@gsnMaximize = True
  res@trXMinF = 1
  res@trXMaxF = 365
  res@trYMinF =  0.
  res@trYMaxF = 1.
  res@xyLineThicknesses  = (/3,3,3,3,3,3/)
  res@xyLineColors       = (/"blue","red","lightblue","lightblue","pink","pink","black"/)
  res@xyDashPatterns     = (/0,0,1,1,1,1,0/)
; res@tiMainString       = "BL99 (solid) MUSHY (dashed)"
  res@tiXAxisString = "Day"
  res@tiYAxisString = "Albedo"


  x = fspan(1,365,365)
  alb = (/alb1clim,alb2clim,alb1clim-alb1std,alb1clim+alb1std,alb2clim-alb2std,alb2clim+alb2std,abs(alb1clim-alb2clim)/)
  plot1 = gsn_csm_xy(wks2,x,alb,res)

end

