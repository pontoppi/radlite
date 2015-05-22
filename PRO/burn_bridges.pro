pro burn_bridges, bridges
	ncpus = n_elements(bridges)
	for cpu=0,ncpus-1 do $
		obj_destroy, bridges[ncpus-cpu-1] ;burn in opposite order as created or IDL hangs.
end
