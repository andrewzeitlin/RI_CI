/*

Program to calculate confidence intervals.

ri_ci ///
	[if] [in] 	
	, dgp() ///  	DGP
	ci0() /// 		initial guess 
	: regress y t 	/// model 
	, cl(schoolid) /// model options  

*/

program define ri_ci, rclass
	//  preserve 

	//  parse command, with elements passed to ri_estimates().
	_on_colon_parse `0'
	local estimator `"`s(after)'"' // estimator to pass to 
	local 0 `"`s(before)'"'

	//  unpack meta options for ri_ci	
	syntax , ///
		permutations(integer) ///
		t1(string asis) ///  only one dimension of randomization allowed at a time here.  Sub-options (filename() keyvar())
		teststat(string) /// what to base p-values off of.
		dgp(string asis) ///
		ci0(numlist sort min=2 max=2) ///  	For now, require manual specification of ci0
		[ ///
			NUMTRIals(integer 0) ///  	number of trials -- one possible stopping rule. 
			TOLerance(real 0) /// 	alternative stopping rule -- as a fraction of the width of the initial guess. Initial version is quite large.
			SIGnificance_level(real 0.05) /// significance level for CI.
		] 


	//  check input 
	if wordcount("`ci0'") ~= 2 {
		di as err "Confidence interval requires two scalar values"
		exit 
	}
	if real(word("`ci0'",1))== . | real(word("`ci0'",2)) == . {
		di as err "Both initial CI values must be real numbers."
		exit 
	}
	if `numtrials' == 0 & `tolerance' == 0 {
		di as err "User must specify either a number of trials or a value for the tolerance."
	}

	//  Declarations 


	//  Extract parameters of the search 
	//  Consider:  start at midpoint of search value?  Or conduct an internal estimate, check that the point estimate is on the interior of the interval, and then use that estimate as the midpoint?  
	//	For now, will take the midpoint and assume that the starting point for the CI is well chosen.
	tokenize `ci0' 
	local ci0_lb = `1' // min(real(word("`ci0'",1)),real(word("`ci0'",2)))
	local ci0_ub = `2' // max(real(word("`ci0'",1)),real(word("`ci0'",2)))
	di as err "The suggested starting interval runs from lb `ci0_lb' to ub `ci0_ub'"
	local m0 = (`ci0_lb' + `ci0_ub') / 2 

	//  Parse t1
	parse_tx `t1'
	local t1vars `s(txvars)'
	local t1file `s(txfile)'
	local t1key `s(keyvar)'
	// ToDo: save data here so that ri_estimates() can bring in t1file again?
	// Or modify ri_estimates() to make merging of t1file optional?

	if `"`t2'"' ~= "" {
		parse_tx `t2'
		local t2vars `s(txvars)' 
		local t2file `s(txfile)'
		local t2key `s(keyvar)' 
		// display as err `"Assignments for variables `t2vars' can be found in file `t2file'"'
	}

	//  ToDo:  check whether we need to merge in the treatments in order to estimate using the actual assignment--or is the actual assignment variable already present. 
	//  Require it to be present and then update ri_estimates() to assume this is the case?  
	if "`t1file'" ~= "" {
		quietly {
			merge 1:1 `t1key' using `t1file', nogen assert(3) 
			if ( `"`t2'"' ~= "" ) merge 1:1 `t2key' using `t2file', nogen assert(3) update replace // allowing for the possibility that t_0 is already in the data.
		}
	}

	//  restrict to sample of interest 
	if `"`if'"' ~= "" keep `if' 
	if `"`in'"' ~= "" keep `in' 

	//  Evaluate the model using the original data, to obtain the test statistic against which all RI will be conducted.
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

	//  The estimate
	`cmd' `depvar' `tx' `interaction_vars' `controls', `options' 

	//  Saving the test statistic of interest 
	foreach x in `tx' `interaction_vars' {
		if "`teststat'" == "b" local test_`x' = _b[`x']
		else if "`teststat'" == "t" local test_`x' = _b[`x']/_se[`x'] 
		else {
			di as err "Cannot currently accommodate test statistics other than regression coefficients or t-stats."
			exit 
		} 
		di as err "The test statistic for variable `x' is `test_`x''" 
	}

	****************************************************************************
	/*  identifying the confidence interval  */
	****************************************************************************

	tempvar y0 // this will hold the implied `control' outcome under hypothesized treatment effect tau0.  
	ge `y0' = .

	//  Let's start by confirming the upper side of the 95% CI.	  
		local trialno = 1 // counter for successfully executed trials.

		//  Evaluate p-value at upper bound.   
		local tau0 = `ci0_ub'
		evaluate_trial, dgp(`dgp') treatment(`ci0_ub') y0(`y0') estimator(`estimator') depvar(`depvar') permutations(`permutations') teststat(`teststat') values(`test_`t1vars'') t1vars(`t1vars')
		mat TRIALS_UB = r(THISTRIAL) // initialize list of trial outcomes for upper bound of 95% CI.

		//  Confirm initial value for upper bound of CI is big enough.
		if el(TRIALS_UB,`trialno',colnumb(TRIALS_UB,"pvalue")) >  `significance_level' / 2 {
			di as err "Initial value for upper bound of CI not big enough.  p-value `thispvalue' associated with treatment effect `tau0'."
			exit
		}

		//  Figure out next trial.
		local bottom `m0' 
		local top `ci0_ub'
		local middle = 0.5*(`m0' + `ci0_ub') 

		//  Loop with number of trials as the stopping rule.
		while `trialno' < `numtrials' {
			//  update counter
			local ++trialno 

			//  evaluate current candidate (middle) value
			evaluate_trial, dgp(`dgp') treatment(`middle') y0(`y0') estimator(`estimator') depvar(`depvar') permutations(`permutations') teststat(`teststat') values(`test_`t1vars'') t1vars(`t1vars')

			//  store results
			mat TRIALS_UB = TRIALS_UB \ r(THISTRIAL) 
			local thispvalue = el(TRIALS_UB,`trialno',colnumb(TRIALS_UB,"pvalue"))
			di "p-value for trial `trialno', treatment `middle', is `thispvalue'"

			//  update bottom, middle, top depending on outcome above.
			if (`thispvalue' > `significance_level' / 2 ) local bottom = `middle' //  move right 
			else local top = `middle' // move left.
			local middle = 0.5*(`bottom' + `top') // next trial value.
		}

		//  TODO:  loop with step size as stopping rule.

	//  Now identify lower bound of the 95% CI. 

	//  restore 

****************************************************************************
end
****************************************************************************


****************************************************************************
/*  Sub-functions  */
****************************************************************************

//  program to evaluate p-value for non-zero sharp null. Wrapper.
program define evaluate_trial , rclass
	
	syntax, dgp(string asis) t1vars(varname) treatment(real) y0(varname) estimator(string asis) depvar(varname) permutations(integer) teststat(string) values(real) 

	impose_tx, dgp( `dgp' ) treatment(`treatment') y(`y0') subtract 
	local dgp0 = subinstr(`"`dgp'"', word("`dgp'",1), "`y0'",1) // updating DGP for this to be based on the new dependent variable 
	local estimator0 = subinstr(`"`estimator'"', "`depvar'" , "`y0'", 1)
	ri_estimates, permutations(`permutations') t1( `t1vars' ) ///
		teststat(`teststat') pvalues values(`values')  ///
		dgp(`dgp0') treatmenteffect(`treatment')  ///
		:  `estimator0'

	local thispvalue = el(r(RESULTS),rownumb(r(RESULTS),"`t1vars'"),colnumb(r(RESULTS),"p"))
	di as err "The p-value at candidate treatment effect `treatment' on variable `t1vars' is `thispvalue'"

	//  Store results of first trial, and return to calling program.
	tempname THISTRIAL 
	mat `THISTRIAL' = J(1,3,.)
	mat coln `THISTRIAL' = tau0 permutations pvalue 		
	mat `THISTRIAL'[1,1] = `treatment'
	mat `THISTRIAL'[1,2] = `permutations'
	mat `THISTRIAL'[1,3] = `thispvalue'
	return matrix THISTRIAL = `THISTRIAL'

end

// Program to parse lists of treatment dimensions and corresponding variables 
program define parse_tx , sclass 
	syntax namelist, [ filename(string) KEYvar(varname) ] 
	sreturn local txvars `namelist' 
	if "`filename'" ~= "" {
		sreturn local txfile `filename' 
		sreturn local keyvar `keyvar'
	}
end
