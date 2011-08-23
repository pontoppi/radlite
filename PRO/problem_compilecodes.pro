;==========================================================================
;    SCRIPT FOR COMPILATION OF RADICAL, RADMC, DIFFUSION AND VERTSTRUCT
;                  AND ALSO FOR THE DUST OPACITY MAKER
;                 (PARAMETERS ARE IN PROBLEM_PARAMS.PRO)
;==========================================================================
@problem_makeopac.pro
@problem_comp_files.pro
;
; Natural constants
;
AU     = 1.496d13
cc     = 2.9979d10
RS     = 6.96d10
MS     = 1.99d33
LS     = 3.8525d33     
pc     = 3.08572d18
ss     = 5.6703d-5  
kk     = 1.3807d-16   
mp     = 1.6726d-24  
GG     = 6.672d-8 
;
; Get hostname
;
spawn,'echo `hostname`',hostname
hostname=strcompress(hostname,/remove_all)
hostname=hostname[0]
;
; Read the parameters
;
@problem_params.pro
;
; Check if the executable directory exists
;
result = findfile(bindir,count=count)
if count eq 0 then begin
   spawn,'mkdir '+bindir
endif
;
; Find the number of frequency points
;
useopac,infile,pll,fresmd=fresmd,/onlynf,nf=nf,nspec=nspec
;
; If nr_compile is defined, and not smaller than nr itself, then
; use this one. Else use nr.
;
if n_elements(nr_compile) ne 0 then begin
   if nr_compile lt nr then begin
      print,'ERROR: nr_compile < nr ! Cannot compile this way!'
      stop
   endif
endif else begin
   nr_compile = nr
endelse
;
; If nt_compile is defined, and not smaller than nt, then
; use this one. Else use nt.
;
if n_elements(nt_compile) ne 0 then begin
;   if nt_compile lt nt+ntex then begin
   if nt_compile lt nt then begin
;      print,'ERROR: nt_compile < nt+ntex ! Cannot compile this way!'
      print,'ERROR: nt_compile < nt ! Cannot compile this way!'
      stop
   endif
endif else begin
;   nt_compile = nt+ntex
   nt_compile = nt
endelse
;
; If nf_compile is defined, and not smaller than nf itself, then
; use this one. Else use nf.
;
if n_elements(nf_compile) ne 0 then begin
   if nf_compile lt nf then begin
      print,'ERROR: nf_compile < nf ! Cannot compile this way!'
      stop
   endif
endif else begin
   nf_compile = nf
endelse
;;;
;;; If nsize_compile is defined, and not smaller than nsz itself, then
;;; use this one. Else use nf.
;;;
;;if n_elements(nsize_compile) ne 0 then begin
;;   if nsize_compile lt nsz then begin
;;      print,'ERROR: nsize_compile < nsz ! Cannot compile this way!'
;;      stop
;;   endif
;;endif else begin
;;   nsize_compile = nsz
;;endelse
nsize_compile = 1
;
; If nlevr_compile is defined, and not smaller than nlevr itself, then
; use this one. Else use nlevr.
;
if n_elements(nlevr_compile) ne 0 then begin
   if n_elements(nlevr) gt 0 then begin
      if nlevr_compile lt nlevr then begin
         print,'ERROR: nlevr_compile < nlevr ! Cannot compile this way!'
         stop
      endif
   endif
endif else begin
   ;;nlevr_compile = nlevr
   nlevr_compile = 10 ; Dummy value for RADICAL, which is disabled anyway
endelse
;
; If nspec_compile is defined, and not smaller than nspec itself, then
; use this one. Else use nspec.
;
if n_elements(nspec_compile) ne 0 then begin
   if nspec_compile lt nspec then begin
      print,'ERROR: nspec_compile < nspec ! Cannot compile this way!'
      stop
   endif
endif else begin
   nspec_compile = nspec
endelse
;
; Find current working directory
;
cd,current=current
;
;
;;
;;======================================================================
;;                        COMPILE RADICAL 
;;======================================================================
;;
;; Check if RADICAL source directory exists
;;
;result = findfile(src_rdcl,count=count)
;if(count eq 0) then begin
;   print,'ERROR: Directory for RADICAL source files not found.'
;   print,'   I tried: ',src_rdcl
;   cd,current
;   stop
;endif
;;
;; Now move to RADICAL directory
;;
;cd,src_rdcl
;spawn,'cd '+src_rdcl
;;
;; Message
;;
;print,'cd '+src_rdcl
;print,'===== Compiling RADICAL (ray-trace mode) ====='
;;
;; A Second RADICAL configure.h file, this time only for ray-tracing
;;
;write_conf_radical,nr_compile,nt_compile,nf_compile,nlevr_compile,$
;                   nspec_compile,/onlyray
;;
;; Remove radicalprog
;;
;spawn,'rm -f radicalprog'
;;
;; Now compile RADICAL
;;
;spawn,'make'
;;
;; Check if radicalprog was generated
;;
;result = findfile('radicalprog',count=count)
;if count eq 0 then begin
;   print,'ERROR: Somehow RADICAL (in ray-tracing mode) did not compile'
;   cd,current
;   stop
;endif
;result = findfile('radical',count=count)
;if count eq 0 then begin
;   print,'ERROR: Somehow RADICAL (in ray-tracing mode) did not compile'
;   cd,current
;   stop
;endif
;;
;; Leave these executables here, so that the current version of RADICAL
;; is the ray-tracing version! That saves a lot of memory. 
;;
;; Now remove all .o files, since for each new compilation we need to
;; recompile everything anyway.
;;
;spawn,'rm -f *.o'
;;
;
;======================================================================
;                        COMPILE RAYTRACE
;======================================================================
;
; Check if RAYTRACE source directory exists
;
result = findfile(src_rayt,count=count)
if(count eq 0) then begin
   print,'ERROR: Directory for RAYTRACE source files not found.'
   print,'   I tried: ',src_rayt
   cd,current
   stop
endif
;
; Now move to RAYTRACE directory
;
cd,src_rayt
spawn,'cd '+src_rayt
;
; Message
;
print,'cd '+src_rayt
print,'===== Compiling RAYTRACE ====='
;
; A configure.h file similar to RADICAL
;
write_conf_radical,nr_compile,nt_compile,nf_compile,nlevr_compile,$
                   nspec_compile,/onlyray
;
; Remove raytraceprog
;
spawn,'rm -f raytraceprog'
;
; Now compile RAYTRACE
;
spawn,'make'
;
; Check if raytraceprog was generated
;
result = findfile('raytraceprog',count=count)
if count eq 0 then begin
   print,'ERROR: Somehow RAYTRACE did not compile'
   cd,current
   stop
endif
result = findfile('raytrace',count=count)
if count eq 0 then begin
   print,'ERROR: Somehow RAYTRACE did not compile'
   cd,current
   stop
endif
;
; Copy this executable to the bin/ directory
; BUGFIX 04.04.07: also copy raytraceprog
;
spawn,'cp -f raytrace '+bindir
spawn,'cp -f raytraceprog '+bindir
;
; Now remove all .o files, since for each new compilation we need to
; recompile everything anyway.
;
spawn,'rm -f *.o'
;
;
;======================================================================
;                        COMPILE RADMC
;======================================================================
;
; Check if RADMC source directory exists
;
result = findfile(src_rdmc,count=count)
if(count eq 0) then begin
   print,'ERROR: Directory for RADMC source files not found.'
   print,'   I tried: ',src_rdmc
   cd,current
   stop
endif
;
; Now move to RADMC directory
;
cd,src_rdmc
spawn,'cd '+src_rdmc
;
; Message
;
print,'cd '+src_rdmc
print,'===== Compiling RADMC ====='
;
; Make the configure.h file (which is the same as the RADICAL one)
;
write_conf_radical,nr_compile,nt_compile,nf_compile,nlevr_compile,$
                   nspec_compile,ntempdb=3000,ndiff=400
;
; Remove radmc
;
spawn,'rm -f radmc'
;
; Now compile RADMC
;
spawn,'make'
;
; Check if radmc was generated
;
result = findfile('radmc',count=count)
if count eq 0 then begin
   print,'ERROR: Somehow RADMC did not compile'
   cd,current
   stop
endif
;
; Copy this executable to the bin/ directory
;
spawn,'cp -f radmc '+bindir
;
; Now remove all .o files, since for each new compilation we need to
; recompile everything anyway.
;
spawn,'rm -f *.o'
;
;
;======================================================================
;                 COMPILE THE DUST OPACITY MAKER CODE
;              (CODE TAKEN FROM MODUST BY ALEX DE KOTER)
;======================================================================
;
; Check if the dust opacity generator code directory exists
;
result = findfile(src_dust,count=count)
if(count eq 0) then begin
   print,'ERROR: Directory for dust code source files not found.'
   print,'   I tried: ',src_dust
   cd,current
   stop
endif
;
; Now move to dust code directory
;
cd,src_dust
spawn,'cd '+src_dust
;
; Message
;
print,'cd '+src_dust
print,'===== Compiling DUST OPACITY CODE ====='
;
; Make the configure.h file
;
if strlen(dustdata) lt 60 then begin
   dustdata1=dustdata
   dustdata2=''
endif else begin
   if strlen(dustdata) ge 120 then begin
      print,'SORRY: dust data path too long...'
      stop
   endif
   dustdata1=strmid(dustdata,0,60)
   dustdata2=strmid(dustdata,60,60)
endelse
write_conf_dust,nf_compile,nsize_compile,dustdata1,dustdata2
;
; Remove original dust code
;
spawn,'rm -f makedust'
;
; Now compile dust opacity code
;
spawn,'make'
;
; Check if dust code was generated
;
result = findfile('makedust',count=count)
if count eq 0 then begin
   print,'ERROR: Somehow the dust opacity code did not compile'
   cd,current
   stop
endif
;
; Copy this executable to the bin/ directory
;
spawn,'cp -f makedust '+bindir
;
; Now remove all .o files, since for each new compilation we need to
; recompile everything anyway.
;
spawn,'rm -f *.o'
;
;
;
;======================================================================
;                  COMPILE THE CHOPPING ROUTINE
;======================================================================
;
; Check if the chopping code directory exists
;
result = findfile(src_chop,count=count)
if(count eq 0) then begin
   print,'ERROR: Directory for chopping code source files not found.'
   print,'   I tried: ',src_chop
   cd,current
   stop
endif
;
; Now move to dust code directory
;
cd,src_chop
spawn,'cd '+src_chop
;
; Message
;
print,'cd '+src_chop
print,'===== Compiling CHOPPING CODE ====='
;
; Write configure.h file
;
write_conf_chop,nr_compile,nt_compile,nf_compile,nspec_compile
;
; Remove original chopping code
;
spawn,'rm -f chopdens'
;
; Now compile chopping code
;
spawn,'g77 chopdens.F dust.F -o chopdens'
;
; Check if chopping code was generated
;
result = findfile('chopdens',count=count)
if count eq 0 then begin
   print,'ERROR: Somehow the chopping code did not compile'
   cd,current
   stop
endif
;
; Copy this executable to the bin/ directory
;
spawn,'cp -f chopdens '+bindir
;
; Now remove all .o files, since for each new compilation we need to
; recompile everything anyway.
;
;spawn,'rm -f *.o'
;
;
;
;======================================================================
;                      COMPILE PAH CODE
;======================================================================
;
; Check if PAHCODE source directory exists
;
result = findfile(src_pahc,count=count)
if(count eq 0) then begin
   print,'ERROR: Directory for PAHCODE source files not found.'
   print,'   I tried: ',src_pahc
   cd,current
   stop
endif
;
; Now move to PAHCODE directory
;
cd,src_pahc
spawn,'cd '+src_pahc
;
; Message
;
print,'cd '+src_pahc
print,'===== Compiling PAHCODE ====='
;
; A configure.h file similar to RADICAL
;
write_conf_pah,nr_compile,nt_compile,nf_compile,nspec_compile
;
; Remove pahcode
;
spawn,'rm -f pahcode'
;
; Now compile PAHCODE
;
spawn,'make'
;
; Check if pahcode was generated
;
result = findfile('pahcode',count=count)
if count eq 0 then begin
   print,'ERROR: Somehow PAHCODE did not compile'
   cd,current
   stop
endif
;
; Copy this executable to the bin/ directory
;
spawn,'cp -f pahcode '+bindir
;
; Now remove all .o files, since for each new compilation we need to
; recompile everything anyway.
;
spawn,'rm -f *.o'
;
;
;
;======================================================================
;                             DONE
;======================================================================
;
; Back to 'current' dir
;
cd,current
spawn,'cd '+current
;
; Done...
;
print,' '
print,'  ====== ALL COMPILATIONS SUCCESFUL ======'
print,' '
end

