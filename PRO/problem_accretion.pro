
@problem_natconst.pro
@problem_params.pro


;-----------------------------------------------------------------------------------
; If there is no accretion the stellarsource.dat will be deleted if
; existed at all
;-----------------------------------------------------------------------------------
acc_lumtot = 0
if acc_rate eq 0. and acc_lumtot eq 0 then begin
    res = findfile('stellarsource.dat', count=count)
    if count gt 0 then spawn, 'rm stellarsource.dat'

    goto, the_end
endif

;-----------------------------------------------------------------------------------
; Begin the real accretion stuff
;-----------------------------------------------------------------------------------

;-----------------------------------------------------------------------------------
; Read the grid and the density for grid refinement
;-----------------------------------------------------------------------------------
openr, 1, 'frequency.inp'
  readf, 1, nfreq
  nu = dblarr(nfreq)
  readf, 1, nu
close, 1

openr, 1, 'radius.inp'
  readf, 1, nr_old
  r_old = dblarr(nr_old)
  readf, 1, r_old
close, 1

openr, 1, 'theta.inp'
  readf, 1, nt_old, dum
  theta_old = dblarr(nt_old)
  readf, 1, theta_old
close, 1

openr,1,'dustdens.inp'
readf,1,nrspec,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
dens=dblarr(sizer,sizetbig,nrspec)
idum=0
for ispec=1,nrspec do begin
    for ir=1,sizer do begin
        for it=1,sizet do begin
            dum = 0.d0
            readf,1,dum
            dens(ir-1,it-1,ispec-1) = dum
            if(imirt eq 1) then begin
                itt = 2*sizet + 1 - it
                dens(ir-1,itt-1,ispec-1) = dum
            endif
        endfor
    endfor
endfor
close,1

;-----------------------------------------------------------------------------------
; Define the new radial grid
;-----------------------------------------------------------------------------------
rin_old = min(r_old)
r_add   = acc_rin * (rin_old/acc_rin)^(dindgen(acc_addnr)/acc_addnr)
r_new   = [r_add, r_old]
nr_new  = n_elements(r_new)
;-----------------------------------------------------------------------------------
; Define the new theta grid
;-----------------------------------------------------------------------------------
nt_new = nt_old+1
theta_new = [theta_old, !dpi/2d0*0.999d0]
;-----------------------------------------------------------------------------------
; Write out the new radius and theta grid and density distribution
;-----------------------------------------------------------------------------------

openw, 1, 'theta.inp'
printf, 1, fix(nt_new), 1
printf, 1, '  '
for it=0, nt_new-1 do $
   printf, 1, theta_new(it)
close, 1
sizet = nt_new

openw, 1, 'radius.inp'
printf, 1, fix(nr_new)
printf, 1, ' '
for ir=0, nr_new-1 do $
   printf, 1, r_new(ir)
close, 1

openw,1,'dustdens.inp'
printf,1,fix(nrspec),fix(nr_new),fix(sizet),fix(imirt)
for ispec=1,nrspec do begin
   for ir=1,nr_new do begin

      if r_new(ir-1) lt rin_old then begin         
         for it=1,sizet do begin
            printf, 1, 0d0
         endfor
      endif else begin
         for it=1,sizet do begin
            printf, 1, dens(ir-1-acc_addnr, it-1, ispec-1)
         endfor
      endelse        

    endfor
endfor
close,1

;-----------------------------------------------------------------------------------
; Calculate the temperature distribution according to Eq(5.73)
;  in L. Hartmann : Accretion processes in star formation
;-----------------------------------------------------------------------------------


if keyword_set(acc_lumtot) then begin
   acc_rate = acc_lumtot*2.*rstar/(gg*mstar)   
endif


tacc = ((3d0 * gg * mstar * acc_rate(0)) / (8d0 * !dpi * r_new^3. * ss) * $
       (1d0 - (rstar/r_new)^0.5d0))^0.25d0

lnu = dblarr(nfreq)

;-----------------------------------------------------------------------------------
; Write out the temperature distribution
;-----------------------------------------------------------------------------------

openw, 1, 'acc_temp.inp'
printf, 1, n_elements(r_new)

for ir=0, n_elements(r_new)-1 do begin
    if r_new(ir) le acc_rout then begin
        printf, 1, r_new(ir), tacc(ir)
    endif else begin
        printf, 1, r_new(ir), 0.0
    endelse
endfor
close, 1


;-----------------------------------------------------------------------------------
; We should specify the amount of enegy which comes out from a unit
; volume during a unit of time. Since in the following we will
; calculate the power output of a unit surface area, here we calculate
; the conversion factor to normalize the power output to unit volume.
;-----------------------------------------------------------------------------------

res = compute_vol()
r_i = res.r_i
nri = n_elements(r_i)
cfact = dblarr(nr_new)

for ix=0, nr_new-1 do begin   
   surf  = (r_i(ix+1)^2d0 - r_i(ix)^2d0)*!dpi
   cfact(ix) = res.vol(ix)/surf
endfor

;-------------------------------------------------------------------------
; Calculate the SED of the Shakura-Sunyaev disk by the analytical formulas
;  (Only for debugging purposes)
;-------------------------------------------------------------------------

lnu  = 0d0
clnu = 0d0
for ir=0, nr_new-1 do begin
   if r_new(ir) le acc_rout then begin
       temp  = sqrt(tacc(ir+1)*tacc(ir))
       da    = (r_i(ir+1)^2d0 - r_i(ir)^2d0)*!dpi
       bnu   = 2.d0*hh*nu^3d0/cc^2d0/(exp(hh*nu/kk/temp)-1d0) 
       lnu   = lnu + bnu*da 
       clnu  = clnu + bnu*da*cfact
   endif
endfor

;-------------------------------------------------------------------------
; Calculate the accretion luminosity in the midplane as it will be put in
;    the RADMC
;-------------------------------------------------------------------------

acc_lum = dblarr(nfreq, sizet, nr_new)
for ir=1, nr_new do begin
   if r_new(ir-1) lt acc_rout then begin
       bnu = 2.d0*hh*nu^3d0/cc^2d0/(exp(hh*nu/kk/tacc(ir-1))-1d0)
       for it=1, sizet do begin
           if it eq sizet then begin
               for inu=1, nfreq do begin
                   acc_lum(inu-1, it-1, ir-1) = bnu(inu-1) / cfact(ir-1) * 0.5
               endfor
           endif
       endfor   
   endif 
endfor

;-----------------------------------------------------------------------------------
; Write out the accretion luminosity 
;-----------------------------------------------------------------------------------

openw, 1, 'stellarsource.dat'
printf, 1, fix(nfreq), fix(nr_new), fix(sizet), 1
printf, 1, ' '

for inu=1, nfreq do begin
   for it=1, sizet do begin
         if it eq sizet then begin
            for ir=1, nr_new do begin
               printf, 1, acc_lum(inu-1, it-1, ir-1) 
            endfor
         endif else begin
            for ir=1, nr_new do begin
               printf, 1, 0d0 
            endfor
         endelse
   endfor
endfor
close,1
;-----------------------------------------------------------------------------------
; Calculate the luminosity of the SS disk
;-----------------------------------------------------------------------------------

sslum = 0d0
for i=0, n_elements(nu)-2 do begin
    sslum = sslum + 2.0 * (nu(i+1)-nu(i)) * (lnu(i+1)+lnu(i))*0.5D0 
endfor


;----------------------------------------------------------------------------------
; Add the radiation field of the hot spot on the stellar surface if we set it
;
; Assumptions:
;                - The luminosity of the spot is half of the total accretion
;                  luminosity (acc_type=0)
;                - The luminosity of the spot is (1-rstar/acc_rin) x Lacc
;                  (acc_type=1)
;                - The surface area of the spot is acc_ssr * 4. * pi * rstar^2
;----------------------------------------------------------------------------------
if add_spot eq 1 then begin
    
    sratio = acc_ssr

;-----------------------------------------------------------------------------------
; Calculate the spot temperature
;-----------------------------------------------------------------------------------
    acclum_full = gg * mstar * acc_rate / rstar

    spotsurf = 4.d0 * !dpi * rstar^2. * acc_ssr

    if acc_type eq 0 then spot_lumratio = 0.5
    if acc_type eq 1 then spot_lumratio = (1d0 - rstar/acc_rin)
;
; The temperature of the spot is calculated by simply summing the
; luminosities of the accretion and that of the star
;
    star_lum_cont = spotsurf * ss * tstar^4.
    spottemp = ((acclum_full*spot_lumratio + star_lum_cont) / spotsurf / ss)^0.25

    starlum  = 4.d0 * !dpi * rstar^2. * ss * tstar^4.
             
    readcol, 'starspectrum.inp', nu, starspec    
    nu  = double(nu)
    wav = 2.9979d14/nu
    
;-----------------------------------------------------------------------------------
; Calculate the radiation of the spot
;-----------------------------------------------------------------------------------
    
    bnu = 2.*hh*nu^3./cc^2./(exp(hh*nu/kk/spottemp(0))-1D0)
    spot_spec = bnu * spotsurf(0) / 4.d0 / pc^2. 
        
;-----------------------------------------------------------------------------------
; Add the spot radiation field to the stellar radiation field and write it out
;-----------------------------------------------------------------------------------
    starspec_new = starspec * (1.-acc_ssr) + spot_spec

    openw, 1, 'starspectrum.inp'
    printf, 1, n_elements(starspec_new)
    for i=0, n_elements(starspec_new)-1 do $
      printf, 1, nu(i), starspec_new(i)
    
    close, 1
    sflux = starspec*pc^2. * (1.-acc_ssr)
    spot_spec = spot_spec*pc^2.
endif else begin
    spot_spec = 0.0
    readcol, 'starspectrum.inp', nu, sflux
    sflux = sflux*pc^2
endelse
;-----------------------------------------------------------------------------------
; Plot the results
;-----------------------------------------------------------------------------------

window,2,retain=2, xs=800, ys=600
!p.multi=[2,0,2]
spec = read_spectrum(file='starspectrum.inp')
wav  = 2.9979d14/spec.freq(1:spec.nfr)
nu   = 2.9979d14/wav

tlnu = lnu + sflux + spot_spec

red = 1. / (140.*pc^2.)

plot_oo, wav, tlnu*nu*red, xr=[0.3, 1d2], xs=1, xtitle='!4k!X [!4l!Xm]', $
  ytitle='!4m!L!D!4m!X!N [erg/s/cm/cm]', title='SED at 140pc', charsize=1.3
oplot, wav, spot_spec*nu*red, col=11119999
oplot, wav, sflux*nu*red, col=118899
oplot, wav, lnu*nu*red, col=190

legend, ['Total', 'Spot', 'Star', 'Disk'], psym=0, col=[15777215, 11119999, 118899, 190],$
  charsize=1.3, /right_legend


save, file='acc_spec.sav', wav, lnu, sflux, spot_spec, tacc

plot_oo, r_new/au, tacc, xtitle='R [AU]', ytitle='T [K]', $
  title=' Temperature in the SS disk'

;-----------------------------------------------------------------------------------
; Print some parameters
;-----------------------------------------------------------------------------------
sim_spotlum = 0d0
for inu=0, n_elements(nu)-2 do $
  sim_spotlum = sim_spotlum + (nu(inu+1) - nu(inu)) * (spot_spec(inu+1) + spot_spec(inu)) * 0.5

sim_disklum = 0d0
for inu=0, n_elements(nu)-2 do $
  sim_disklum = sim_disklum + (nu(inu+1) - nu(inu)) * (lnu(inu+1) + lnu(inu)) * 0.5

sim_starlum = 0d0
for inu=0, n_elements(nu)-2 do $
  sim_starlum = sim_starlum + (nu(inu+1) - nu(inu)) * (sflux(inu+1) + sflux(inu)) * 0.5


print, '***********************************************'
print, ' Stellar Luminosity      : ', sim_starlum*4.*!dpi
print, ' Accretion Luminosity    : ', sim_disklum *!dpi*2.
if add_spot eq 1 then begin
    print, ' Spot         Luminosity : ', sim_spotlum * 4. * !dpi
    print, ' Spot/Stellar Surface    : ', spotsurf/(4.*!dpi*rstar^2.)
    print, ' Spot temperature        : ', spottemp, ' K'
endif
print, '***********************************************'
;-----------------------------------------------------------------------------------
; The End
;-----------------------------------------------------------------------------------
the_end:


