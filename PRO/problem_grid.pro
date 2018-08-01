;--------------------------------------------------------------------
;                         MAKE THE R-GRID
;--------------------------------------------------------------------
function make_rgrid,rin,rout,nr,rrefine=rrefine
if not keyword_set(rrefine) then begin
    r = rin*(rout/rin)^(findgen(nr)/(nr-1.d0))
endif else begin
    nextra = rrefine.nspanr*rrefine.nlevr*(2^rrefine.nstepr-1)
    if nextra gt nr-5 then begin
       print,"Sorry not enough nr for this refinement level"
       stop
    endif
    r = rin*(rout/rin)^(findgen(nr-nextra)/(nr-1.d0-nextra))
    for ilev=1,rrefine.nlevr do begin
        for ispan=rrefine.nspanr,1,-1 do begin
            rins=fltarr(2^rrefine.nstepr-1)
            fr=(r(ispan)/r(ispan-1))^(0.5^rrefine.nstepr)
            for i=1,2^rrefine.nstepr-1 do rins(i-1)=r(ispan-1)*fr^i
            r=[r(0:ispan-1),rins,r(ispan:*)]
        endfor
    endfor
 endelse
;
; Check the DeltaR/R
;
nnr  = n_elements(r)
if nnr ne nr then stop
drr  = (r[nr-1]-r[nr-2])/r[nr-2]
if drr gt 0.15 then begin
   print,'ERROR: Radial grid too coarse...'
   stop
endif
if drr lt 0.05 then begin
    print,'Radial grid too finely spaced... This will take too much'
    print,'computational time... Are you sure this is okay? (Type 1)'
    read,i
    if i ne 1 then stop
endif
return,r
end


;--------------------------------------------------------------------
;                       MAKE THE THETA-GRID
;
; Produce the theta grid for the RADMC simulation
;
; NOTE: Contrary to an earlier version (not used here) this routine
;       assumes that nt is the total numberr of theta points above the
;       equator. That is: n_elements(theta).eq.nt
; EXCEPTION: This is no longer true if zrefine is used. That is 
;       because zrefine is really necessary for a very fine midplane
;       layer of dust or so. You then don't want to lose other 
;       theta points.
;--------------------------------------------------------------------
function make_tgrid,hrgrid,nt,hrgmax=hrgmax,ntex=ntex,$
                    hrlg=hrlg,zrefine=zrefine
if keyword_set(ntex) then ntt=nt-ntex else ntt=nt
if ntt lt 4 then stop
if not keyword_set(hrlg) then begin
    thmax = !pi/2.d0
    thmin = !pi/2.d0 - hrgrid
    theta = (thmax-thmin)*findgen(ntt)/(ntt-1.d0+0.5d0)+thmin
    if keyword_set(hrgmax) then begin
        if not keyword_set(ntex) then begin
            print,'PROBLEM: If you define hrgmax, must also define nextra'
            stop
        endif 
        if hrgmax le hrgrid then begin
            print,'PROBLEM: hrgmax must be larger than hrgrid'
            stop
        endif
        thmax = !pi/2.d0 - hrgrid
        thmin = !pi/2.d0 - hrgmax
        thex=(thmax-thmin)*findgen(ntex)/(ntex*1.d0)+thmin
        theta = [thex,theta]
    endif
endif else begin
    B     = alog(hrgrid)
    A     = (B-alog(hrlg))/(nt-1.d0)
    theta = !pi/2.d0-exp(-A*findgen(nt)+B)
endelse
if keyword_set(zrefine) then begin
    ;hunt,theta[*],!pi/2-zrefine.zrref,iz
    theta = [theta[0:iz-1],!pi/2-rotate(dindgen(zrefine.nzref)+0.5,2)*$
             (!pi/2-theta[iz-1])/(1.d0*zrefine.nzref+1.)]
endif
return,theta
end


;------------------------------------------------------------------------
;                     REFINE THE RADIAL GRID
;------------------------------------------------------------------------
pro refine,r,irstart,nlevr=nlevr,nspanr=nspanr,nstepr=nstepr
if irstart gt 0 then begin
   rold=r[0:irstart-1]
   r=r[irstart:*]
endif
for ilev=1,nlevr do begin
   for ispan=nspanr,1,-1 do begin
      rins=fltarr(2^nstepr-1)
      fr=(r(ispan)/r(ispan-1))^(0.5^nstepr)
      for i=1,2^nstepr-1 do rins(i-1)=r(ispan-1)*fr^i
      r=[r(0:ispan-1),rins,r(ispan:*)]
   endfor
endfor
if irstart gt 0 then begin
   r=[rold,r]
endif
end
