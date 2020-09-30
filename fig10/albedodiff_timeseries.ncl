
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   

begin

  a1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h.fswdn.000101-005012.nc","r")
  a2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h.fswdn.000101-005112.nc","r")
  b1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h.fswabs.000101-005012.nc","r")
  b2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h.fswabs.000101-005112.nc","r")
  c1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h.aice.000101-005012.nc","r")
  c2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h.aice.000101-005112.nc","r")
  d1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h.snowfrac.000101-005012.nc","r")
  d2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h.snowfrac.000101-005112.nc","r")
  e1 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99snow1.cice.h.apeff_ai.000101-005012.nc","r")
  e2 = addfile("/glade/p/cesm/pcwg/dbailey/archive/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1/ice/proc/tseries/month_1/b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1.cice.h.apeff_ai.000101-005112.nc","r")

  fswdn1 = a1->fswdn
  fswdn2 = a2->fswdn
  fswabs1 = b1->fswabs
  fswabs2 = b2->fswabs
  aice1 = c1->aice
  aice2 = c2->aice
  snowfrac1 = d1->snowfrac
  snowfrac2 = d2->snowfrac
  apeff_ai1 = e1->apeff_ai
  apeff_ai2 = e2->apeff_ai

  tareatmp = a1->tarea
  tlat = a1->TLAT
  tlon = a1->TLON

  ndims = dimsizes(fswdn1)
  ntimes = ndims(0)

  alb1 = new((/ntimes/),"float")
  alb2 = new((/ntimes/),"float")
  sf1 = new((/ntimes/),"float")
  sf2 = new((/ntimes/),"float")
  ap1 = new((/ntimes/),"float")
  ap2 = new((/ntimes/),"float")


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
     sf1(n) = dim_sum(dim_sum(snowfrac1(n,:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
     sf2(n) = dim_sum(dim_sum(snowfrac2(n,:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
     ap1(n) = dim_sum(dim_sum(apeff_ai1(n,:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
     ap2(n) = dim_sum(dim_sum(apeff_ai2(n,:,:)*tarea(:,:)*ice_mask(:,:)))/tot_tarea
  end do

  alb1a = reshape(alb1,(/50,12/))
  alb2a = reshape(alb2,(/50,12/))
  sf1a = reshape(sf1,(/50,12/))
  sf2a = reshape(sf2,(/50,12/))
  ap1a = reshape(ap1,(/50,12/))
  ap2a = reshape(ap2,(/50,12/))

  alb1clim = dim_avg_n(alb1a,0)
  alb2clim = dim_avg_n(alb2a,0)
  alb1clim(0:2) = 0.0
  alb1clim(9:11) = 0.0
  alb2clim(0:2) = 0.0
  alb2clim(9:11) = 0.0
  sf1clim = dim_avg_n(sf1a,0)
  sf2clim = dim_avg_n(sf2a,0)
  ap1clim = dim_avg_n(ap1a,0)
  ap2clim = dim_avg_n(ap2a,0)

  ones = ap1clim*0. + 1.
  ice1clim = ones-sf1clim-ap1clim
  ice2clim = ones-sf2clim-ap2clim

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

  wks2 = gsn_open_wks("pdf", "albdiff_nh")    ; Open PS file "contour3a.ps"

  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize = True
  res@trXMinF = 1
  res@trXMaxF = 12
  res@trYMinF = -0.2
  res@trYMaxF = 0.2
  res@xyLineThicknesses  = (/3,3,3,3,3,3,3/)
  res@xyLineColors       = (/"black","red","blue","green","magenta"/)
  res@xyDashPatterns     = (/0,0,0,0,0,0,0/)
; res@tiMainString       = "BL99 - MUSHY"
  res@tiXAxisString = "Month"
  res@tiYAxisString = "Albedo Difference"


  x = fspan(1,12,12)
  albsno = (/0.0,0.0,0.0,0.8,0.78,0.78,0.78,0.78,0.0,0.0,0.0,0.0/)
  albpnd = (/0.0,0.0,0.0,0.4,0.4,0.4,0.4,0.4,0.0,0.0,0.0,0.0/)
  albice = (/0.0,0.0,0.0,0.6,0.6,0.6,0.6,0.6,0.0,0.0,0.0,0.0/)

  albdiff = (sf1clim-sf2clim)*albsno + (ap1clim-ap2clim)*albpnd + (ice1clim-ice2clim)*albice
  print(alb1clim-alb2clim)
  print(albdiff)
; print((sf1clim-sf2clim)*albsno)
; print((ap1clim-ap2clim)*albpnd)
; print((ice1clim-ice2clim)*albice)

  alb = (/alb1clim-alb2clim,(sf1clim-sf2clim)*albsno,(ap1clim-ap2clim)*albpnd,(ice1clim-ice2clim)*albice,albdiff/)
  plot1 = gsn_csm_xy(wks2,x,alb,res)

; Add a legend.
;
  gres = True
  gres@YPosPercent = 30.    ; expressed as %, 0->100, sets position of top border of legend
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 60.    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  lineres@lgLineColors = res@xyLineColors ; line colors
  lineres@lgLineThicknesses = res@xyLineThicknesses                    ; line thicknesses
; lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

  textres = True
  textres@lgLabels = (/"total albedo","snow","pond","bare ice","sum"/); legend labels (required)
; textres@lgLabelColors = res@xyLineColors
  textres@lgLabelFontHeights = (/0.02,0.02,0.02,0.02,0.02,0.02,0.02/)

  plot2 = simple_legend(wks2,plot1,gres,lineres,textres)

  draw(plot2)
  frame(wks2)

end

