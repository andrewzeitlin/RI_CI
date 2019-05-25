

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

global R = 200

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

ri_estimates, permutations($R) key(i) t1(t , filename(`T0')) teststat(t) /* pvalues */ : regress y t
// mat li r(RESULTS)
// global tstat = el(r(RESULTS),rownumb(r(RESULTS),"t"),colnumb(r(RESULTS),"t"))
// global pval  = el(r(RESULTS),rownumb(r(RESULTS),"t"),colnumb(r(RESULTS),"p"))
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


