//  Some preliminaries for stata
clear all
set seed 1234567890
set matsize 11000
set scheme s1mono

//  Add /ado/ subdirectory to top of Stata's search path.
//  Assume we are EITHER in the root of the github repo OR in the /demo/ folder to begin.
if regexm("`c(pwd)'","demo") local pwd = subinstr(subinstr("`c(pwd)'","\","/",.),"/demo", "", 1)
else local pwd "`c(pwd)'"
di as err "The current working directory is `pwd'"
adopath ++`pwd'/ado 

// Force reload of Stata programs
capture program drop ri_ci
capture program drop ri_estimates
capture program drop impose_tx 

qui set obs 500
ge i = _n // an identifier for observations, used later to merge potential treatments
global tau = 1
ge y0 = rnormal()
ge y1 = y0 + $tau

global R = 500 //  number of randomizations to consider.

ge t_0 = (runiform() >= 0.5)
ge y = y0 + t_0*(y1-y0)

/* Other randomizations will be saved in a temp file, to be merged in by ri_ci().  */
tempfile T0 
preserve 
keep i t_0
forvalues r=1/$R {
    ge t_`r' = (runiform() >= 0.5 )
}
qui save `T0'
restore 

capture program drop ri_estimates
ri_estimates, permutations($R) ///  Number of permutations to consider, for each candidate value.
    t1(t , filename(`T0') key(i)) ///  dataset containing alternative values of the treatment effect
    teststat(t) ///
    pointestimates ///  <- also report point estimates for the estimating model a
    pvalues ///  <-  calculate p-value for sharp null under consideration (default: tau = 0)
    : regress y t  //  This is the model from which the test statistic is derived.
mat li r(RESULTS)
global tstat = el(r(RESULTS),rownumb(r(RESULTS),"t"),colnumb(r(RESULTS),"t"))
global pval  = el(r(RESULTS),rownumb(r(RESULTS),"t"),colnumb(r(RESULTS),"p"))
mat T0 = r(T0)
di "The t statistic is $tstat, with p-value $pval"

// Visualizing result
preserve
qui drop _all 
qui svmat T0, names(tstat)
twoway (kdensity tstat) ///
    (scatteri 0.4 $tstat 0 $tstat, recast(line) lcolor(blue) lpattern(dash)) ///
    , /// xline($tstat, lcolor(blue))  ///
    xtitle("t-statistics") /// xscale(range($tstat)) /// 
    legend(order(1 "Distribution under sharp null" 2 "Estimate") cols(1) position(11) ring(0))
restore

capture program drop impose_tx
capture drop yminus
capture drop yplus

//  Example where treatment effect is being subtracted
qui ge yminus = .
impose_tx , dgp(y ~ t_0) treatment(1) y(yminus) subtract

//  Example where treatment effect is being added
qui ge yplus = .
impose_tx, dgp(y ~ t_0) treatment(1) y(yplus)

tw (sc yminus y if t_0 == 1, mcolor(red)) ///
    (sc yplus y if t_0 == 1, mcolor(green)) ///
    (sc yminus y if t_0 == 0) ///
    (function y=x , range(-4 4)) /// (line y y, sort) ///
    , legend(order(1 "treated; treatment effect added" 2 "treated; treatment effect subtracted" 3 "control") cols(1)) ///
    ytitle("new outcome") xtitle("actual outcome")

capture program drop ri_estimates
ri_estimates, permutations($R) t1(t , filename(`T0') key(i))  dgp(y ~ t) treatmenteffect(1) /// imposing the *true* DGP
    teststat(t) pvalues values($tstat) /// using previously estimated t-statistic to compute a p-value for this sharp null.
    : regress y t
mat li r(RESULTS) // show p-value 
mat T0alt = r(T0)  // capture distribution of test statistic under sharp null.

// Visualizing result
preserve
drop _all 
qui svmat T0alt, names(tstat)
twoway (kdensity tstat) ///
    (scatteri 0.4 $tstat 0 $tstat, recast(line) lcolor(blue) lpattern(dash)) ///
    , /// xline($tstat, lcolor(blue))  ///
    xtitle("t-statistics") /// xscale(range($tstat)) /// 
    legend(order(1 "Distribution under sharp null" 2 "Estimate") cols(1) position(11) ring(0))
restore

global numTrials = 10 // number of trials to run

//  demo the ri_ci command.
capture program drop ri_ci
ri_ci, numtrials($numTrials) permutations($R) ///
    t1(t, filename(`T0') key(i) ) ///
    teststat(t) ///
    dgp(y ~ t ) ///
    ci0( 5 -5) checkboundary ///
    : reg y t // estimation command delivering the test statistic.

mat TRIALS_UB = r(TRIALS_UB)
mat TRIALS_LB = r(TRIALS_LB)
ret li

di "Trial results in search for upper bound:"
mat li TRIALS_UB
di "Trial results in search for lower bound:"
mat li TRIALS_LB 

//  Visualize results
preserve
clear 
tempfile trials 
svmat TRIALS_UB, names(col)
save `trials' 

clear 
svmat TRIALS_LB, names(col)
append using `trials'

tw sc pvalue tau0 [w=permutations] ///
    , msymbol(oh) ///
    yline(0.025, lcolor(red)) ytitle("p-values") ///
    xline(1, lcolor(blue) lpattern(dash)) xtitle("treatment effect") ///

restore 

capture drop t_*
capture program drop ri_ci
ri_ci, numtrials($numTrials) permutations($R) ///
    t1(t, filename(`T0') key(i) ) ///
    teststat(t) ///
    dgp(y ~ t ) ///
    analytic_initial checkboundary /// ci0( 10 -10) /// noci ///
    /// pzero /// <- to report p-value for zero null
    : reg y t // estimation command delivering the test statistic.
ret list 
mat TRIALS_UB = r(TRIALS_UB)
mat TRIALS_LB = r(TRIALS_LB)

di "Trial results in search for upper bound:"
mat li TRIALS_UB
di "Trial results in search for lower bound:"
mat li TRIALS_LB 

preserve
clear 
tempfile trials 
svmat TRIALS_UB, names(col)
save `trials' 

clear 
svmat TRIALS_LB, names(col)
append using `trials'

tw sc pvalue tau0 [w=permutations] ///
    , msymbol(oh) ///
    yline(0.025, lcolor(red)) ytitle("p-values") ///
    xline(1, lcolor(blue) lpattern(dash)) xtitle("treatment effect") ///

restore 



