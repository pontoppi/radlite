PRO calc_partition_function, isotop, Qstr=Qstr
  @natconst
  nt = 3000
  tmin = 1.
  tmax = 3000.
  Q = FLTARR(nt)
  lambdarange = [0.1,1d90]
  hitran_extract, isotop=isotop, lambdarange=lambdarange
  mol = read_molecule_lambda('moldata.dat')
  T = FINDGEN(nt)/(nt-1)*(tmax-tmin)+tmin
  FOR i=0,nt-1 DO BEGIN
     Q[i] = TOTAL(mol.g*exp(-mol.energy_in_k/T[i]))
  ENDFOR
  Qstr = {T:T, Q:Q}
  stop
END
