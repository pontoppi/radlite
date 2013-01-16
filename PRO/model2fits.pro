PRO model2fits, filename=filename
  @natconst.pro

  IF ~KEYWORD_SET(filename) THEN filename='radlite_run.fits'
  LMAX = 100000L
  ;
  ;Determine number of line spectrum files
  spawn, 'ls linespectrum_moldata_*.dat', linefiles
  spawn, 'ls moldata_*.dat', molfiles
  nfiles = N_ELEMENTS(linefiles)
  
  test_read = read_line(linefiles[0])
  nfreq  = test_read.nfreq

  lines  = dblarr(LMAX,nfreq)
  velos  = dblarr(LMAX,nfreq)
  cfreqs = dblarr(LMAX)
  Eupper = dblarr(LMAX)
  Aud    = dblarr(LMAX)
  gupper = dblarr(LMAX)
  glower = dblarr(LMAX)
  lcount = 0L

  FOR ff=0,nfiles-1 DO BEGIN
     line_str = read_line(linefiles[ff])
     mol_str  = read_molecule_lambda(molfiles[ff])
     
     nl       = line_str.nlines
     nfreq    = line_str.nfreq
     cfreqs[lcount:lcount+nl-1]   = line_str.cfreqs
     Eupper[lcount:lcount+nl-1]   = mol_str.Eupper
     Aud[lcount:lcount+nl-1]      = mol_str.Aud
     gupper[lcount:lcount+nl-1]   = mol_str.g[mol_str.iup]
     glower[lcount:lcount+nl-1]   = mol_str.g[mol_str.idown]

     FOR i=0,nl-1 DO BEGIN
        velos[lcount+i,0:nfreq-1]  = line_str.velo[*,i]
        lines[lcount+i,0:nfreq-1]  = line_str.flux[*,i]
     ENDFOR
     lcount = lcount + nl
  ENDFOR
  
  ;
  ;remove duplicate lines if present
  cfreqs = cfreqs[0:lcount-1]
  Eupper = Eupper[0:lcount-1]
  Aud    = Aud[0:lcount-1]
  gupper = gupper[0:lcount-1]
  glower = glower[0:lcount-1]
  
  uniqsubs  = UNIQ(cfreqs, SORT(cfreqs))
  cfreqs    = cfreqs[uniqsubs]
  Eupper    = Eupper[uniqsubs]
  Aud       = Aud[uniqsubs]
  gupper    = gupper[uniqsubs]
  glower    = glower[uniqsubs]

  velos     = velos[uniqsubs,*]
  lines     = lines[uniqsubs,*]
  lcount    = N_ELEMENTS(uniqsubs)
    
  c_lines = fltarr(lcount)
  
  ;
  ;Calculate integrated line fluxes
  intensity = FLTARR(lcount)
  FOR i=0,lcount-1 DO BEGIN
     nvel = N_ELEMENTS(velos[i,*])
     cont = (lines[i,0]+lines[i,nvel-1])/2.
     freq = (REFORM(1d5*velos[i,*])/cc+1)*cfreqs[i]
     intensity[i] = INT_TABULATED(freq,lines[i,*]-cont[0],/SORT, /DOUBLE)
  ENDFOR

  MWRFITS, dummy, filename, /CREATE, /SILENT
  MWRFITS, {velo:velos,flux:lines,intensity:intensity,cfreqs:cfreqs,Eupper:Eupper,$
            Aud:Aud,gupper:gupper,glower:glower,species:mol_str.species,$
            Eupper_unit:'cm-1',Aud_unit:'s^-1',flux_unit:'erg/cm^2/s/Hz',velo_unit:'km/s'}, filename, /SILENT
  
stop
END
