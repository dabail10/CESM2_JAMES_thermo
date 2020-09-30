
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "readobs.ncl"

begin

dir_asd = "/gpfs/fs1/p/cesm/community/ASD-HIGH-RES-CESM1/RCP8.5/"
;case_asd = "BRCP85C5CN_ne120_t12_pop62.c13b17.asdphys.001"
;files_asd = systemfunc("ls "+dir_asd+case_asd+".cice.h.aice_nh.*")
;case_asd = "BRCP85C5CN_ne120_t12_pop62.c13b17.asdphys.001a"
;files_asd = systemfunc("ls "+dir_asd+case_asd+".cice.h.aice.*.nc")
case_asd = "BRCP85C5CN_ne120_t12_pop62.c13b17.asdphys.001a"
files_asd = systemfunc("ls "+dir_asd+case_asd+".cice.h.aice.209*.nc")

;case_asd = "BRCP85C5CN_ne120_t12_pop62.c13b17.asdphys.001"
;files_asd = systemfunc("ls "+dir_asd+case_asd+".cice.h.aice_sh.*")
;case_asd = "BRCP85C5CN_ne120_t12_pop62.c13b17.asdphys.001a"
;files_asd = systemfunc("ls "+dir_asd+case_asd+".cice.h.aice.*.nc")
;case_asd = "BRCP85C5CN_ne120_t12_pop62.c13b17.asdphys.001a"
;files_asd = systemfunc("ls "+dir_asd+case_asd+".cice.h.aice.209*.nc")

f_asd = addfiles(files_asd,"r")

tarea_asd = f_asd[0]->tarea
tlat_asd = f_asd[0]->TLAT
mask_asd = where(tlat_asd.gt.0.,1.,0.)
;mask_asd = where(tlat_asd.lt.0.,1.,0.)
aice_asd = f_asd[:]->aice

ndims_asd = dimsizes(aice_asd)
ntimes_asd = ndims_asd(0)

extent_asd = new((/ntimes_asd/),"float")
do n=0,ntimes_asd-1
   extent_asd(n) = dim_sum(dim_sum(where(aice_asd(n,:,:)*mask_asd.ge.15.,tarea_asd,0.)))*1.0e-12
end do

;extent_asd_2D = onedtond(extent_asd,(/ntimes_asd/12,12/))
;extent_asd_clim = dim_avg_n(extent_asd_2D(4:11,:),0)
;extent_asd_2D = onedtond(extent_asd(6:ntimes_asd-6),(/(ntimes_asd-6)/12,12/))
;extent_asd_clim = dim_avg_n(extent_asd_2D(0:7,:),0)
extent_asd_2D = onedtond(extent_asd,(/ntimes_asd/12,12/))
extent_asd_clim = dim_avg_n(extent_asd_2D(0:7,:),0)

dir_slices = "/glade/p/cgd/oce/projects/nsfsubmeso/archive/"
case1_slices = "g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.027"
case2_slices = "g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.029"
case3_slices = "g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.028"

files1_slices = systemfunc("ls "+dir_slices+case1_slices+"/ice/hist/*.nc")
files2_slices = systemfunc("ls "+dir_slices+case2_slices+"/ice/hist/*.nc")
files3_slices = systemfunc("ls "+dir_slices+case3_slices+"/ice/hist/*.nc")

f1_slices = addfiles(files1_slices,"r")
f2_slices = addfiles(files2_slices,"r")
f3_slices = addfiles(files3_slices,"r")

tarea_slices = f1_slices[0]->tarea
tlat_slices = f1_slices[0]->TLAT
mask_slices = where(tlat_slices.gt.0.,1.,0.)
;mask_slices = where(tlat_slices.lt.0.,1.,0.)

aice1_slices = f1_slices[:]->aice
aice2_slices = f2_slices[:]->aice
aice3_slices = f3_slices[:]->aice

ndims_slices = dimsizes(aice1_slices)
ntimes_slices = ndims_slices(0)

extent1_slices = new((/ntimes_slices-12/),"float")
extent2_slices = new((/ntimes_slices-12/),"float")
extent3_slices = new((/ntimes_slices-12/),"float")

do n=11,ntimes_slices-2
   extent1_slices(n-11) = dim_sum(dim_sum(where(aice1_slices(n,:,:)*mask_slices.ge.0.15,tarea_slices,0.)))*1.0e-12
   extent2_slices(n-11) = dim_sum(dim_sum(where(aice2_slices(n,:,:)*mask_slices.ge.0.15,tarea_slices,0.)))*1.0e-12
   extent3_slices(n-11) = dim_sum(dim_sum(where(aice3_slices(n,:,:)*mask_slices.ge.0.15,tarea_slices,0.)))*1.0e-12
end do

extent1_slices_2D = onedtond(extent1_slices,(/(ntimes_slices-12)/12,12/))
extent1_slices_clim = dim_avg_n(extent1_slices_2D,0)
extent2_slices_2D = onedtond(extent2_slices,(/(ntimes_slices-12)/12,12/))
extent2_slices_clim = dim_avg_n(extent2_slices_2D,0)
extent3_slices_2D = onedtond(extent3_slices,(/(ntimes_slices-12)/12,12/))
extent3_slices_clim = dim_avg_n(extent3_slices_2D,0)


dir_lens = "/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/"
case_lens = "b.e11.BRCP85C5CNBDRD.f09_g16."
runs_lens = (/"001","002","003","004","005","006","007","008","009","010", \
              "011","012","013","014","015","016","017","018","019","020", \
              "021","022","023","024","025","026","027","028","029","030", \
              "031","032","033","034","035","101","102","103","104","105"/)

; Loop over LENS members

;extent_lens_clim = new((/40,12/),"float")
extent_lens_clim = new((/33,12/),"float")

;do n=0,39
do n=0,32

;files_lens = systemfunc("ls "+dir_lens+"aice/"+case_lens+runs_lens(n)+"*.aice_nh.2006*.nc")
files_lens = systemfunc("ls "+dir_lens+"aice/"+case_lens+runs_lens(n)+"*.aice_nh.2081*.nc")
;files_lens = systemfunc("ls "+dir_lens+"aice/"+case_lens+runs_lens(n)+"*.aice_sh.2006*.nc")
;files_lens = systemfunc("ls "+dir_lens+"aice/"+case_lens+runs_lens(n)+"*.aice_sh.2081*.nc")
f_lens = addfile(files_lens,"r")

tarea_lens = f_lens->tarea
tlat_lens = f_lens->TLAT
mask_lens = where(tlat_lens.gt.0.,1.,0.)
;mask_lens = where(tlat_lens.lt.0.,1.,0.)

aice_lens = f_lens->aice

ndims_lens = dimsizes(aice_lens)
ntimes_lens = ndims_lens(0)
;ntimes_lens = 660

extent_lens = new((/ntimes_lens/),"float")
do m=0,ntimes_lens-1
   extent_lens(m) = dim_sum(dim_sum(where(aice_lens(m,:,:)*mask_lens.ge.15.,tarea_lens,0.)))*1.0e-12
end do

extent_lens_2D = onedtond(extent_lens,(/ntimes_lens/12,12/))
;extent_lens_clim(n,:) = dim_avg_n(extent_lens_2D(4:11,:),0)
;extent_lens_clim(n,:) = dim_avg_n(extent_lens_2D(44:51,:),0)
extent_lens_clim(n,:) = dim_avg_n(extent_lens_2D(10:17,:),0)

delete(extent_lens_2D)
delete(extent_lens)
delete(ndims_lens)
delete(ntimes_lens)
delete(tlat_lens)
delete(tarea_lens)
delete(mask_lens)
delete(aice_lens)

end do

x = fspan(1,12,12)
x!0 = "month"
x@long_name = "Month"

ext = new((/2,37,12/),float)
readobs(ext)

ext@long_name = "NSIDC sea ice extent"
ext@_FillValue = -9999

ext_clim = dim_avg_n(ext(0,29:36,:),0)
;ext_clim = dim_avg_n(ext(1,29:36,:),0)

y3 = (/extent_asd_clim,ext_clim/)
y2 = extent3_slices_clim
y1 = extent_lens_clim

wks1 = gsn_open_wks("pdf","extent")

res = True
res@gsnDraw = False
res@gsnFrame = False
res@trYMaxF = 22.
res@trYMinF = 0.
res@trXMaxF = 12.
res@trXMinF = 1.
res@xyLineThicknesses = (/3,3/)
res@xyDashPattern = (/0,0/)

res1 = res
res2 = res

res2@xyLineColors = (/"red","black"/)
res1@xyLineColors = (/"blue"/)
res@xyLineColor = (/"gray80"/)

plot1 = gsn_csm_xy3(wks1,x,y1,y2,y3,res,res1,res2)

; Add a legend.
;
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 50.    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
; gres@XPosPercent = 10.    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  lineres@lgLineColors = (/"red","blue","gray80","black"/) ; line colors
  lineres@lgLineThicknesses = (/3,3,3,3/)                    ; line thicknesses
  lineres@lgDashIndexes = (/0,0,0,0/)                    ; line thicknesses
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

  textres = True
  textres@lgLabels = (/"ASD RCP8.5","Timeslice","LENS","NSIDC"/); legend labels (required)
  textres@lgLabelFontHeights = (/0.010,0.010,0.010/)

  plot2 = simple_legend(wks1,plot1,gres,lineres,textres)

  draw(plot2)
  frame(wks1)
end
