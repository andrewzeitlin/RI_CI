// Program to parse lists of treatment dimensions and corresponding variables 
program define parse_tx , sclass 
	syntax namelist, [ filename(string) KEYvar(varname) ] 
	sreturn local txvars `namelist' 
	if "`filename'" ~= "" {
		sreturn local txfile `filename' 
		sreturn local keyvar `keyvar'
	}
end
