/*
Program to provide point estimates and RI tests, for estimation procedure of the user's choosing

Allows specifying a set of parameters for each randomization dimension, or combination thereof  **

Borrows some ideas from Simon Heb's -ritest-, adapts to our specific problem.

___ToDo:___
(*) allow custom test statistics 
(*) allow non-zero nulls 

*/
program define ri_estimates, rclass 
	preserve 
	
	//  parse command, and pass only the "before" part to -syntax- for local parsing 
	_on_colon_parse `0'
	local estimator `"`s(after)'"'
	local 0 `"`s(before)'"'
	// di as err `"The syntax for the command to be executed is: `s(after)'"'

	//  Unpack meta options for ri 
	syntax [if] [in]  , ///
		Permutations(integer) /// permutations: number of permutations of (each) treatment vector
		/// Key(name) /// key variable for merging treatments onto main dataset
		t1(string) ///  first dimension of randomization.  Required.
		teststat(string) /// what to base p-values off of.
		[	///
			t2(string) ///  Second dimension of randomization. Optional. 
			POINTestimates ///  Estimate model and return test statistsics for the actual realization of the treatment
			VALues(string asis) ///  value of test statistic(s), if externally estimated. provide this as a list.
			PVALues /// 		Return p-value
			dgp(string asis) /// list of scalars corresponding to variables in t1(). For now, allowing additive treatment effects only.
			treatmenteffect(numlist max=1) /// treatment effects to be associated with variables on the RHS of DGP.  Passed to impose_tx(). Currently allows only one such (scalar) value. 
		]

	//  Checking input 
	if "`teststat'" ~= "b" & "`teststat'" ~= "t" {
		di as err "Must specify basis for p values as either b or t."
		exit
	}
	if ("`pointestimates'" == "" & "`values'" == "" ) & "`pvalues'" ~= "" {
		di as err "Option -pvalues- requires EITHER option -pointestimates- OR option -value- to be specified."
		exit
	}
	if (`"`dgp'"' ~= "" ) & `"`t2'"' ~= "" {
		di as err "Can only specify one of options dgp() and td()."
		exit 
	}

	tempfile actuals  //  tempfile for holding dataset with actual randomization in memory
	tempname RESULTS  //  tempname for matrix holding results (b,t,p, etc.)

	//  Unpacking treatment variables and corresponding assignments 
	parse_tx `t1'
	local t1vars `s(txvars)'
	local t1file `s(txfile)'
	local t1key `s(keyvar)'
	// display as err `"Assignments for variables `t1vars' can be found in file `t1file'"'

	if `"`t2'"' ~= "" {
		parse_tx `t2'
		local t2vars `s(txvars)' 
		local t2file `s(txfile)'
		local t2key `s(keyvar)' 
		// display as err `"Assignments for variables `t2vars' can be found in file `t2file'"'
	}

	//  What will be the actual executing syntax for the command?
	gettoken cmd spec : estimator // executing command 
	gettoken depvar spec : spec  // dependent variable
	local p = strpos("`spec'", ",")
	if `p' == 0 { // if there are no options specified for the estimating command.
		local p0 = length("`spec'") 
		local options "" 
	}
	else { // if options are specified.
		local p0 = `p' - 1 // if there are options
		local p1 = `p' + 1
		di as err `"Local spec is `spec'"'
		local options = substr("`spec'",`p1',.) 
	}
	local rhs = substr("`spec'", 1,`p0') 
	
	// Identify interactions in RHS 
	local rhs = trim(itrim("`rhs'"))
	local rhs = subinstr("`rhs'", " * ", "*", .) 
	local tx ""
	local controls "" 
	local interactions ""
	foreach w in `rhs' {
		// di "Word is: `w'"
		// check if word is an interaction. 
		local p=strpos("`w'", "*")
		//  if so, file it in the interactions macro
		if (`p' > 0) local interactions "`interactions' `w'"
		// otherwise, check if word is in either treatment list
		else {
			local p=strpos("`t1vars' `t2vars'", "`w'" )
			//  if so, assign it to the list of all treatment variables
			if (`p' > 0 ) local tx "`tx' `w'"
			// otherwise, assign this variable to list of controls.
			else local controls "`controls' `w'" 
		}
	} 
	// di as err "RHS: treatments are `tx'"
	// di as err "RHS: controls are `controls'"
	// di as err "RHS: interactions are `interactions'"

	//  Merge randomization(s) into the data, if required. 
	if "`t1file'" ~= "" {
		quietly {
			merge 1:1 `t1key' using `t1file', nogen assert(3) 
			if ( `"`t2file'"' ~= "" ) merge 1:1 `t2key' using `t2file', nogen assert(3) update replace // allowing for the possibility that t_0 is already in the data.
		}
	}

	//  restrict to sample of interest 
	if `"`if'"' ~= "" keep `if' 
	if `"`in'"' ~= "" keep `in' 

	//  Parse interactions
	if "`interactions'" ~= "" {
		tokenize `interactions' 
		local I = wordcount("`interactions'") 
		// local i = 1 
		// while "`1'" ~= "" {
		forvalues i= 1/`I' {
			// di "Now parsing this string: ``i''"
			parse_interactions ``i'' // `1'
			local first_`i' `s(first)'
			local second_`i' `s(second)'
			local interaction_`i' `first_`i''_x_`second_`i''
			// di "For the `i'th interaction, the first term is `first_`i'' and the second is `second_`i''."
			// di "The resulting variable will be called `interaction_`i''."
			local interaction_vars "`interaction_vars' `interaction_`i''"
			// mac shift 
			//  local I = `i' // number of interactions
			// local ++i
		}
	}
	// di as err "There are `I' interactions in the model."
	// di as err "The list of interaction variables is `interaction_vars'"
 	
	//  Containers 
	tempname B1 T1 // name for matrix to hold results on t1vars below 
	if `"`t2'"' == "" { // | "`interaction_vars'" == "" {
		local t1only "`t1vars' `interaction_vars'" 
	}
	else {
		local t1only "`t1vars'" // starting point. 
		// For each interaction variable, determine if it has elements from T1 only, from T2 only, or from both.
		if "`interaction_vars'" ~= "" {
			forvalues i = 1/`I' {
				local p11 : list posof "`first_`i''" in t1vars 
				local p21 : list posof "`second_`i''" in t1vars 
				if `p11' | `p21' {
					local p12 : list posof "`first_`i''" in t2vars 
					local p22 : list posof "`second_`i''" in t2vars 
					if ~`p12' & ~`p22' local t1only "`t1only' `interaction_`i''"
				} 
			}
		}
	}
	local K1 = wordcount("`t1only'") 
	mat `B1' = J(`permutations', `K1',.) 
	mat coln `B1' = `t1only'
	mat `T1' = J(`permutations', `K1',.) 
	mat coln `T1' = `t1only' 

	if `"`t2'"' ~= "" {
		tempname B2 T2 
		if ("`interaction_vars'" == "") local t2only "`t2vars'" 
		else {
			local t2only "`t2vars'" 
			forvalues i=1/`I' {
				local p12 : list posof "`first_`i''" in t2vars 
				local p22 : list posof "`second_`i''" in t2vars 
				if `p12' | `p22' {
					local p11 : list posof "`first_`i''" in t1vars 
					local p21 : list posof "`second_`i''" in t1vars 
					if ~`p11' & ~`p21' local t2only "`t2only' `interaction_`i''"
				} 
			}
		}
		local K2 = wordcount("`t2only'")
		mat `B2' = J(`permutations', `K2',.) 
		mat coln `B2' = `t2only' 
		mat `T2' = J(`permutations', `K2',.) 
		mat coln `T2' = `t2only' 
	}

	//  Gathering point estimates and test statistics 
	//  Pull in treatment variables from specification and construct interactions as needed.
	foreach v in `tx' {
		capture drop `v' 
		quietly gen `v' = `v'_0 
	}
	if "`interactions'" ~= "" {
		forvalues i=1/`I' {
			capture drop `interaction_`i''

			//  source for the first term 
			local p : list posof "`first_`i''" in tx 
			if `p' local x1 `first_`i''_0 
			else local x1 `first_`i''

			//  source for the second term
			local p : list posof "`second_`i''" in tx 
			if `p' local x2 `second_`i''_0 
			else local x2 `second_`i''

			quietly gen `interaction_`i'' = `x1' * `x2' 
		}
	}
	local K = wordcount("`tx' `interaction_vars'")
	
	//  Point estimates using actual realization of the treatment
	if "`pointestimates'" ~= "" {

		di as err "Point estimates, analytical SEs"

		//  Results matrix. Definitely collect b,t. Optionally collect p -- that part comes later.
		mat `RESULTS' = J(`K',2,.)
		mat rown `RESULTS' = `tx' `interaction_vars' 
		mat coln `RESULTS' = b t  
		`cmd' `depvar' `tx' `interaction_vars' `controls', `options' 
		ret scalar N = `e(N)' 
		foreach x in `tx' `interaction_vars' {
			mat `RESULTS'[rownumb(`RESULTS',"`x'"),1] = _b[`x']
			mat `RESULTS'[rownumb(`RESULTS',"`x'"),2] = _b[`x']/_se[`x'] 
		}
	}
	save `actuals' //  preserving dataset merged this way in memory to restore at start of each randomization type

	//	Inference:  first dimension of randomization ONLY 
	di as err "RI: Permutation of T1..."

	tempvar ystar //  this will hold the outcome net of any non-zero sharp null imposed.
	qui ge `ystar' = .
	
	// Extracting parameters for DGP here.
	tokenize `dgp' 
	local k = 1 // TODO: check to make sure l.c. `k' is not holding anything important.
	foreach v in `t1vars' {
		local `v' = `k' // assign effect of variable `k' to a local variable under its own name.
		local ++k 
	}


	forvalues p=1/`permutations' {

		//  Bring in the pth assignment. 
		foreach v in `t1vars' {
			quietly replace `v' = `v'_`p' 
		}
		// Loop over interactions. 
		// For any element in t1vars, replace it with the permutation-p assignment.
		// For any element in t2vars, replace it with the permutation 0 (actual) assignment.
		// For any control/observational attribute, use as is.
		if "`interactions'" ~= "" {
			forvalues i=1/`I' {
				local in_t1 : list posof "`first_`i''" in t1vars 
				if `in_t1' local x1 `first_`i''_`p' 
				else {
					local in_t2 : list posof "`first_`i''" in t2vars 
					if `in_t2' local x1 `first_`i''_0 
					else local x1 `first_`i''
				}

				local in_t1 : list posof "`second_`i''" in t1vars
				if `in_t1' local x2 `second_`i''_`p' // if this is part of what is being permuted
				else {
					local in_t2 : list posof "`second_`i'" in t2vars
					if `in_t2' local x2 `second_`i''_0 
					else local x2 `second_`i''
				}
				// now, construct the interaction variable 
				quietly replace `interaction_`i'' = `x1' * `x2' 
			}
		}

		//  Impose DGP for non-sharp zero nulls. 
		if `"`dgp'"' == "" qui replace `ystar' = `depvar' // zero null: no dgp specified.
		else {
			impose_tx, dgp(`dgp') treatmenteffect(`treatmenteffect') y(`ystar')			
		}
		qui replace `ystar' = `depvar'
		if `"`dgp'"' ~= "" {
			foreach k in `t1vars' {
				qui replace `ystar' = `ystar' + `k'*``k''
			}
		}

		//  Extract test statistic 
		quietly `cmd' `ystar' `tx' `interaction_vars' `controls', `options'
		foreach v of varlist `t1only' { 
			mat `B1'[`p',colnumb(`B1',"`v'")] = _b[`v']
			mat `T1'[`p',colnumb(`T1',"`v'")] = _b[`v']/_se[`v'] 
		}
	}
	di as err "... done."
	// mat li `B1' 


	//  [Optionally]:  Inference:  second dimension of randomization.
	if `"`t2'"' ~= "" {
		di as err "RI: Permutation of T2..."
		use `actuals' , clear 
		forvalues p=1/`permutations' {
			foreach v in `t2vars' {
				quietly replace `v' = `v'_`p' 
			}
			// Loop over interactions. 
			// For any element in t2vars, replace it with the permutation-p assignment.
			// For any element in t1vars, replace it with the permutation 0 (actual) assignment.
			// For any control/observational attribute, use as is.
			if "`interactions'" ~= "" {
				forvalues i=1/`I' {
					local in_t2 : list posof "`first_`i''" in t2vars 
					if `in_t2' local x1 `first_`i''_`p' 
					else {
						local in_t1 : list posof "`first_`i''" in t1vars 
						if `in_t1' local x1 `first_`i''_0 
						else local x1 `first_`i''
					}

					local in_t2 : list posof "`second_`i''" in t2vars
					if `in_t2' local x2 `second_`i''_`p' // if this is part of what is being permuted
					else {
						local in_t1 : list posof "`second_`i'" in t1vars
						if `in_t1' local x2 `second_`i''_0 
						else local x2 `second_`i''
					}
					// now, construct the interaction variable 
					quietly replace `interaction_`i'' = `x1' * `x2' 
				}
			}
			quietly `cmd' `depvar' `tx' `interaction_vars' `controls', `options'
			foreach v of varlist `t2only' { 
				mat `B2'[`p',colnumb(`B2',"`v'")] = _b[`v']
				mat `T2'[`p',colnumb(`T2',"`v'")] = _b[`v']/_se[`v'] 
			}
		}
		// mat li `B2' 
		di as err "... done."
	}

	//  [Optionally]:  Any test stats that require randomization in BOTH dimensions?
	local needboth = 0 // defining outside the below to condition on it for output later.
	if `"`t2'"' ~= "" & "`interactions'" ~= "" {
		// di as err "Now checking if there are any dual-treatment interactions: `interaction_vars'"
		local tx_interacted = "" 
		forvalues i = 1/`I' {
			local p1 : list posof "`first_`i''" in tx 
			local p2 : list posof "`second_`i''" in tx 
			if `p1' & `p2' {
				local needboth = 1
				local tx_interacted "`tx_interacted' `interaction_`i''"
			}
		}
		if `needboth' {
			di as err "RI: Simulatneous permutation of T1, T2"

			//  Matrix to store results 
			tempname B12 T12
			local k12 = wordcount("`tx_interacted'")
			mat `B12' = J(`permutations', `k12',.) 
			mat coln `B12' = `tx_interacted' 
			mat `T12' = J(`permutations', `k12',.) 
			mat coln `T12' = `tx_interacted' 
			//  mat li `B12'

			use `actuals', clear 
			forvalues p=1/`permutations' {
				foreach v in `t1vars' `t2vars' {
					quietly replace `v' = `v'_`p' 
				}
				forvalues i=1/`I' {
					local in_t1 : list posof "`first_`i''" in t1vars 
					local in_t2 : list posof "`first_`i''" in t2vars 
					if (`in_t1' | `in_t2') local x1 `first_`i''_`p' 
					else  local x1 `first_`i''

					local in_t1 : list posof "`second_`i''" in t1vars
					local in_t2 : list posof "`second_`i'" in t2vars
					if (`in_t1' | `in_t2') local x2 `second_`i''_`p' // if this is part of what is being permuted
					else local x2 `second_`i''
					
					// now, construct the interaction variable 
					quietly replace `interaction_`i'' = `x1' * `x2' 
				}
				quietly `cmd' `depvar' `tx' `interaction_vars' `controls', `options'
				foreach v of varlist `tx_interacted' { 
					mat `B12'[`p',colnumb(`B12',"`v'")] = _b[`v']
					mat `T12'[`p',colnumb(`T12',"`v'")] = _b[`v']/_se[`v'] 
				}
			}
			di as err "... done."
			// mat li `B12'
		}
	}

	//  Accumulate results 
	tempname B0 T0
	mat `B0' = `B1' 
	mat `T0' = `T1' 
	if (`"`t2'"' ~= "") {
		mat `B0' = [`B0',`B2']	
		mat `T0' = [`T0',`T2']	
	} 
	if `needboth' {
		mat `B0' = [`B0', `B12'] 
		mat `T0' = [`T0', `T12'] 
	}


	// calculate p-values and return, if requested
	if "`pvalues'" ~= "" {
		//  Container
		tempname Pvalues 
		mat `Pvalues' = J(`K',1,.)
		mat coln `Pvalues' = p
		mat rown `Pvalues' = `tx' `interaction_vars' 
		
		// Save RI results to Stata's data in memory
		drop _all 
		set obs 0 
		if "`teststat'" == "b" qui svmat `B0', names(col) 
		else qui svmat `T0', names(col) 

		//  Foreach variable about which RI is conducted, compute p-value.
		local k = 1 
		foreach x in  `tx' `interaction_vars' {
			// If values of test statistic were provided rather than estimated, these are the default.
			if ("`values'" ~= "") local teststat = word("`values'",`k')
			else {
				if "`teststat'" == "b" local teststat = `RESULTS'[rownumb(`RESULTS',"`x'"),1]
				else local teststat = `RESULTS'[rownumb(`RESULTS',"`x'"),2]
			}
			qui count if abs(`teststat') < abs(`x')
			qui count if `x' < `teststat' 
			local p_left = `r(N)' / `permutations' 
			qui count if `x' > `teststat' 
			local p_right = `r(N)' / `permutations' 

			mat `Pvalues'[rownumb(`Pvalues',"`x'"),1] = min(2*min(`p_left',`p_right'),1) // following Green Lab SOP. Works for cases where not centered on zero.  

			local ++k
		}
		// Add p-values to results matrix.
		// mat li `Pvalues'
		if ( "`pointestimates'" ~= "" ) mat `RESULTS' = `RESULTS',`Pvalues'
		else mat `RESULTS' = `Pvalues' 
	}

	// Return results 
	ret mat B0 = `B0'
	ret mat T0 = `T0' 
	// if ( "`pointestimates'" ~= "" ) ret mat RESULTS = `RESULTS',`Pvalues'
	// else ret mat RESULTS = `Pvalues'
	ret mat RESULTS = `RESULTS' 

	//  Restore original dataset
	restore 
end

// Program to parse lists of treatment dimensions and corresponding variables 
program define parse_tx , sclass 
	syntax namelist, [ filename(string) KEYvar(variable) ] 
	sreturn local txvars `namelist' 
	if "`txfile'" ~= "" {
		sreturn local txfile `filename' 
		sreturn local keyvar `keyvar'
	}
end


// Program to parse interaction pairs as given 
program define parse_interactions, sclass 
	syntax anything 
	tokenize `0' , parse("*")
	sreturn local first `1'
	sreturn local second `3' 	// n.b. tokenize keeps the parse character as part of the sequence, hence using `3' rather than `2' here
end
