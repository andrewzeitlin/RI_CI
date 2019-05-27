

//  Some preliminaries for stata
clear all
set seed 12345
set matsize 11000
set scheme s1mono

//  Add /ado/ subdirectory to top of Stata's search path.
//  Assume we are EITHER in the root of the github repo OR in the /demo/ folder to begin.
if regexm("`c(pwd)'","demo") local pwd = subinstr(subinstr("`c(pwd)'","\","/",.),"/demo", "", 1)
else local pwd "`c(pwd)'"
di as err "The current working directory is `pwd'"
adopath ++`pwd'/ado 

// Force reload of Stata programs
capture program drop ri_estimates
capture program drop impose_tx 

set obs 500
global tau = 1
ge y0 = rnormal()
ge y1 = y0 + $tau

global R = 100 //  number of randomizations to consider.

ge t_0 = (runiform() >= 0.5)
ge y = y0 + t_0*(y1-y0)

ge i = _n // an identifier for observations, used later to merge potential treatments
tempfile T0 
preserve 
keep i t_0
forvalues r=1/$R {
    ge t_`r' = (runiform() >= 0.5 )
}
save `T0'
restore 

capture program drop ri_estimates
ri_estimates, permutations($R) t1(t , filename(`T0') key(i)) teststat(t) pointestimates pvalues : regress y t
mat li r(RESULTS)
global tstat = el(r(RESULTS),rownumb(r(RESULTS),"t"),colnumb(r(RESULTS),"t"))
global pval  = el(r(RESULTS),rownumb(r(RESULTS),"t"),colnumb(r(RESULTS),"p"))
mat T0 = r(T0)
di "The t statistic is $tstat, with p-value $pval"

// Visualizing result
preserve
drop _all 
svmat T0, names(tstat)
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

//  demo the ri_ci command.
ri_ci, numtrials(10) permutations($R) ///
    t1(t, filename(`T0') key(i) ) ///
    teststat(t) ///
    dgp(y ~ t ) ///
    ci0( 10 -10) ///
    : reg y t // estimation command delivering the test statistic.


//  Visualize results
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
    xline(1, lcolor(blue) lpattern(dash)) ytitle("treatment effect")
restore 

capture drop t_*
capture program drop ri_ci
ri_ci, numtrials(10) permutations($R) ///
    t1(t, filename(`T0') key(i) ) ///
    teststat(t) ///
    dgp(y ~ t ) ///
    analytic_initial /// ci0( 10 -10) /// noci ///
    pzero ///
    : reg y t // estimation command delivering the test statistic.
ret list 

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
    xline(1, lcolor(blue) lpattern(dash)) ytitle("treatment effect") ///
    =name(analytic_bounds, replace )
restore 
