;-------------------------------------------------------------------------
;                WRITE THE BASIC PART OF THE RADLITE.INP FILE
;-------------------------------------------------------------------------
pro write_radlite_base
printf,1,'104             Input format version'
printf,1,'================================================================'
printf,1,'1               Maximum number iterations                                    '
printf,1,'0               Iteration method          (0=LI,1=ALI,2=ALI+Ng,-1=IS)        '
printf,1,'0               Flux conservation trick   (0.e0=disable)'
printf,1,'1.000e-08       Convergence tolerance                    '                    
printf,1,'0               Convergence crit type     (0=J_nu,1=j_nu) '                   
printf,1,'0               Initial guess type        (0=triv)         '                  
printf,1,'============================================================================'
printf,1,'31              Nr of mu-angle points'
printf,1,'16              Nr of phi-angle points'
printf,1,'2               Type of mu gridding       (0=regular,1=dr/r,2=dr/rspecial)'
printf,1,'0.1             Largest allowed error     '
printf,1,'1.0             Dmu wrt dR                (1.0=tangent ray method)'
printf,1,'3               Extra mus around mu=0'
printf,1,'================================================================'
printf,1,'0               Type of outer boundary    (0=empty,2=microwavebg)'
printf,1,'2               Type of inner boundary    (0=solid,1=empty?,2=central star,3=mirror)'
printf,1,'0               Type of equator boundary  (0=no,1=disk)                      '
printf,1,'================================================================'
printf,1,'57.3            Default inclination angle                        '            
printf,1,'100              Nr Phi-points circular CCD                        '           
printf,1,'1               The number of b_i of the image per r_i of the grid '          
printf,1,'-20             The number of extra b_i inside the inner radius     '         
printf,1,'============================================================================'
printf,1,'2               Save NLTE                                                    '
printf,1,'0               Save intensity at inu                                        '
printf,1,'-1              Save source at inu                                           '
printf,1,'-1              Save moments             (-1=save,2=append,3=onlylast) '
printf,1,'-1              Save dusttemp/etc        (-1=save,2=append,3=onlylast)'
printf,1,'0               Save ALI operator'
printf,1,'1               Save phys vars of medium'
printf,1,'1               Save flux conservation data'
printf,1,'============================================================================'
end


;-------------------------------------------------------------------------
;                 WRITE THE RADLITE.INP FILE - FOR LINES
;-------------------------------------------------------------------------
pro write_radlite_line,niter,noisrf=noisrf,cir_np=cir_np,b_per_r=b_per_r,b_extra=b_extra

IF NOT KEYWORD_SET(cir_np) THEN cir_np = 300
IF NOT KEYWORD_SET(b_per_r) THEN b_per_r = 1
IF NOT KEYWORD_SET(b_extra) THEN b_extra = 20

openw,lun,'radlite.inp',/get_lun
printf,lun,'104             Input format version'
printf,lun,'================================================================'
printf,lun,niter,'        Maximum number iterations                                    '
printf,lun,'2               Iteration method          (0=LI,1=ALI,2=ALI+Ng,-1=IS)        '
printf,lun,'0               Flux conservation trick   (0=disable)'
printf,lun,'1.000e-06       Convergence tolerance                    '                    
printf,lun,'0               Convergence crit type     (0=J_nu,1=j_nu) '                   
printf,lun,'1               Initial guess type        (0=triv)         '                  
printf,lun,'============================================================================'
printf,lun,'41              Nr of mu-angle points'
printf,lun,'16              Nr of phi-angle points'
printf,lun,'2               Type of mu gridding       (0=regular,1=dr/r,2=dr/rspecial)'
printf,lun,'0.1             Largest allowed error     '
printf,lun,'1.0             Dmu wrt dR                (1.0=tangent ray method)'
printf,lun,'3               Extra mus around mu=0'
printf,lun,'================================================================'
if keyword_set(noisrf) then begin
    printf,lun,'0               Type of outer boundary    (0=empty,2=microwavebg)'
endif else begin
    printf,lun,'3               Type of outer boundary    (0=empty,2=microwavebg)'
endelse
printf,lun,'2               Type of inner boundary    (0=solid,1=empty?,2=central star,3=mirror)'
printf,lun,'0               Type of equator boundary  (0=no,1=disk)                      '
printf,lun,'================================================================'
printf,lun,'57.3            Default inclination angle                        '            
printf,lun, STRTRIM(STRING(cir_np),2)  + '              Nr Phi-points circular CCD                        '           
printf,lun, STRTRIM(STRING(b_per_r),2) + '               The number of b_i of the image per r_i of the grid '          
printf,lun, STRTRIM(STRING(b_extra),2) + '             The number of extra b_i inside the inner radius     '         
printf,lun,'============================================================================'
printf,lun,'1               Save NLTE                                                    '
printf,lun,'0               Save intensity at inu                                        '
printf,lun,'-1              Save source at inu                                           '
printf,lun,'0               Save moments             (-1=save,2=append,3=onlylast) '
printf,lun,'-1              Save dusttemp/etc        (-1=save,2=append,3=onlylast)'
printf,lun,'0               Save ALI operator'
printf,lun,'1               Save phys vars of medium'
printf,lun,'0               Save flux conservation data'
printf,lun,'============================================================================'
printf,lun,'-551            Type of setup'
printf,lun,'0               Do dust?'
printf,lun,'1               Do lines?'
printf,lun,'1               Include dust in lines?'
printf,lun,'1               Include star pumping?'
printf,lun,'============================================================================'
close,lun
free_lun, lun
end



