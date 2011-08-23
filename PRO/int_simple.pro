FUNCTION int_simple, X, F, cum_str=cum_str

  nx  = N_ELEMENTS(x)
  cum = DBLARR(nx)
  FOR i=1,nx-1 DO BEGIN
     cum[i] = cum[i-1] + 0.5d0*(X[i]-X[i-1])*(F[i]+F[i-1])
  ENDFOR
  RETURN, cum
END
