/*

Program to impose a treatment effect, either adding it (treating existing data like y0) or subtracting it (treating existing data like a mixture of y1,y0)

*/


program define impose_tx, rclass 
	syntax [if] [in] , ///
		dgp(string) ///
		Treatmenteffect(real) /// for now, assuming just one treatment variable
		y(varname) /// variable to hold name of newly generated outcome y 
		[SUBtract] //  for y --> y0, as opposed to y0 --> y 

	parse_dgp `"`dgp'"' 
	// di as err `"The outcome is `s(outcome)'"'
	// di as err `"The treatment is `s(rhs)' "'  // for now, assuming just one rhs variable.
	local outcome `s(outcome)'
	local rhs `s(rhs)'

	//  Count # of RHS variables in DGP. 
	//  For now, program accommodates only one-dimensional treatment.
	local K = wordcount("`rhs'")
	if `K' > 1 {
		di as err "Sorry: must specify only one treatment variable on RHS."
		exit 
	}

	// Return variable of interest.  
	if ("`subtract'" ~= "" ) qui replace `y' = `outcome' - `treatmenteffect'*`rhs'
	else qui replace `y' = `outcome' + `treatmenteffect'*`rhs'
end


//  Program to interpret Wilkinson-Rogers notation for the DGP 
//  See for reference http://www.jerous.org/att/2016/05/11/wilkinson-rogers/wilkinson2formula.html
program define parse_dgp , sclass 
	syntax anything 
	tokenize `0', parse("~")
	sreturn local outcome "`1'"
	sreturn local rhs "`3'"
end		
