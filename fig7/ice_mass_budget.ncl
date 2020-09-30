
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   

begin

  mass1 = readAsciiTable("NCAR_CESM2f09g17_BL99_001_iceant.txt",13,"float",4)
  mass2 = readAsciiTable("NCAR_CESM2f09g17_MUSHY_001_iceant.txt",13,"float",4)

  mass1b = reshape(mass1,(/50,12,13/))
  mass2b = reshape(mass2,(/51,12,13/))

  mass1_avg = dim_avg_n(mass1b(25:49,:,4:12),0)
  mass2_avg = dim_avg_n(mass2b(25:49,:,4:12),0)
  mass1_var = dim_variance_n(mass1b(25:49,:,4:12),0)
  mass2_var = dim_variance_n(mass2b(25:49,:,4:12),0)

  printVarSummary(mass1_avg)
  printVarSummary(mass2_avg)

  wks2 = gsn_open_wks("pdf", "ice_mass_budget")    ; Open PS file "contour3a.ps"
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize = True
  res@trXMinF = 1
  res@trXMaxF = 12
  res@trYMinF = -3.0e+9
  res@trYMaxF = 2.0e+9
; res@trYMinF = -5.0e+8
; res@trYMaxF = 5.0e+8
  res@xyLineThicknesses  = (/3,3,3,3,3,3,3,3,3/)
  res@xyLineColors       = (/"cyan","red","blue","green","magenta","black"/)
  res@xyDashPatterns     = (/0,0,0,0,0,0,0,0,0/)
; res@tiMainString       = "BL99 (solid) MUSHY (dashed)"
; res@tiMainString       = "MUSHY - BL99"

  res2 = res
  res2@xyDashPatterns     = (/1,1,1,1,1,1,1,1,1/)
  res2@tmYRLabelsOn = False

  prob = ttest(mass1_avg,mass1_var,25,mass2_avg,mass2_var,25,False,False)

  masstotp1 = 0.
  masstotn1 = 0.
  masstot2 = 0.
  masstot3 = 0.
  masstot4 = 0.
  masstot5 = 0.
  masstot6 = 0.

  do n=0,11
     if (mass1_avg(n,8) .gt. 0.) then
        masstotp1 = masstotp1 + mass1_avg(n,8)
     end if
     if (mass1_avg(n,8) .lt. 0.) then
        masstotn1 = masstotn1 - mass1_avg(n,8)
     end if
     masstot2 = masstot2 + mass1_avg(n,1)
     masstot3 = masstot3 + mass1_avg(n,2)
     masstot4 = masstot4 + mass1_avg(n,3)
     masstot5 = masstot5 + mass1_avg(n,5)
     masstot6 = masstot6 + mass1_avg(n,0)
  end do

  print("BL99")
  print(masstotp1*1.0e-9)
  print(masstotn1*1.0e-9)
  print(masstot2*1.0e-9)
  print(masstot3*1.0e-9)
  print(masstot4*1.0e-9)
  print(masstot5*1.0e-9)
  print(masstot6*1.0e-9)

  masstotp1 = 0.
  masstotn1 = 0.
  masstot2 = 0.
  masstot3 = 0.
  masstot4 = 0.
  masstot5 = 0.
  masstot6 = 0.

  do n=0,11
     if (mass2_avg(n,8) .gt. 0.) then
        masstotp1 = masstotp1 + mass2_avg(n,8)
     end if
     if (mass2_avg(n,8) .lt. 0.) then
        masstotn1 = masstotn1 - mass2_avg(n,8)
     end if
     masstot2 = masstot2 + mass2_avg(n,1)
     masstot3 = masstot3 + mass2_avg(n,2)
     masstot4 = masstot4 + mass2_avg(n,3)
     masstot5 = masstot5 + mass2_avg(n,5)
     masstot6 = masstot6 + mass2_avg(n,0)
  end do

  print("MUSHY")
  print(masstotp1*1.0e-9)
  print(masstotn1*1.0e-9)
  print(masstot2*1.0e-9)
  print(masstot3*1.0e-9)
  print(masstot4*1.0e-9)
  print(masstot5*1.0e-9)
  print(masstot6*1.0e-9)

  x = fspan(1,12,12)
  y1 = (/mass1_avg(:,0),mass1_avg(:,1),mass1_avg(:,2),mass1_avg(:,3),mass1_avg(:,5),mass1_avg(:,8)/)
  y2 = (/mass2_avg(:,0),mass2_avg(:,1),mass2_avg(:,2),mass2_avg(:,3),mass2_avg(:,5),mass2_avg(:,8)/)
  y3 = (/prob(:,0),prob(:,1),prob(:,2),prob(:,3),prob(:,5),prob(:,8)/)
  ydiff = where(y3.le.0.05,y2-y1,0.)
; ydiff = y2-y1

  plot1 = gsn_csm_xy2(wks2,x,y1,y2,res,res2)
; plot1 = gsn_csm_xy(wks2,x,ydiff,res)

; Add a legend.
;
  gres = True
  gres@YPosPercent = 35.    ; expressed as %, 0->100, sets position of top border of legend
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5.    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  lineres@lgLineColors = res@xyLineColors ; line colors
  lineres@lgLineThicknesses = res@xyLineThicknesses                    ; line thicknesses
; lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

  textres = True
  textres@lgLabels = (/"congelation","frazil","top melt","bottom melt","snow ice","total"/); legend labels (required)
; textres@lgLabelColors = res@xyLineColors
  textres@lgLabelFontHeights = (/0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02/)

  plot2 = simple_legend(wks2,plot1,gres,lineres,textres)

  draw(plot2)
  frame(wks2)
end

