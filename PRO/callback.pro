pro callback, status, error, bridge, ud
	COMPILE_OPT hidden
		
	npop     = bridge->getvar('npop')
	ini_npop = bridge->getvar('ini_npop')
	lte_npop = bridge->getvar('lte_npop')
	
	(*(ud.p_npop_all))[*,*, ud.i]     = npop
	(*(ud.p_npop_ini_all))[*,*, ud.i] = ini_npop
	(*(ud.p_npop_lte_all))[*,*, ud.i] = lte_npop
end
