function [pvalue TEST1 TEST0 y0 ] = ri_ci(DATA, outcome, txvars, tau0 , T0, P, varargin) % model, stat, varargin

	%  Function to conduct RI, including (optionally) confidence intervaals and test of the no-effect null.

	if size(tau0,1) > size(tau0,2) 
		error('null treatment vector tau0 should be (1 x K) not (K x 1)');
	end

	%  Parse inputs
	params = inputParser ; 
	addOptional(params,'Controls',{'unspecified'});  % observational RHS variables
	addOptional(params,'Clusters',{'unspecified'});  % for clustered standard errors in linear model
	addOptional(params,'Model',{'lm'});  % model type
	addOptional(params,'TestType',{'tStat'});
	addOptional(params,'TestSide',{'twosided'}); % test type for the primary test
	parse(params,varargin{:}); 

	g = params.Results.Clusters;
	model = params.Results.Model; 
	TestType = params.Results.TestType; 
	TestSide = params.Results.TestSide;
	xvars = params.Results.Controls;
	if strcmp(xvars{1}, 'unspecified'), xvars = {}; end

	%  Containers for results 
	TEST0 = NaN(P,length(txvars)); % to hold null distribution for test statistic.

	%  Estimate the model using the actual assignment
	sprintf('__RESULTS OF ANALYTICAL MODEL:__');
	if strcmp(model,'lm')
		lm = fitlm(DATA(:,[txvars xvars outcome]))
		TEST1 = table2array(lm.Coefficients([txvars],[TestType]))';
	end

	[ pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0,xvars, model,T0,P,TestType,TestSide,TEST1) ; 

	sprintf('The p-value from Randomization Inference for the hypothesis that tau = %02.2f is %0.2f',tau0, pvalue)
end


%  Subfunction to conduct RI for a particular value
function [pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0,xvars, model, T0, P,TestType,TestSide,TEST1)
	%  p-value for hypothesized sharp null.
	%  Run DGP in reverse to get y0
	y0 = table2array(DATA(:,outcome)) - table2array(DATA(:,txvars)) * tau0' ;

	%  Now, loop over feasible randomizations, impose treatment effect, re-estimate, and extract test statistic
	x = table2array(DATA(:,xvars)); % for speed.
	for pp = 1 : P

		%  Impose hypothesized DGP
		ystar = y0 + T0(:,pp,:) * tau0' ; % accmmodates possibliity of multiple treatment variables

		%  Estimate model and collect test statistic
		if strcmp(model,'lm')
			lm = fitlm([permute(T0(:,pp,:),[1 3 2]) , x], ystar);
			testStat = table2array( ...
				lm.Coefficients(2:1+length(txvars) ... % leaving room for constant term
				, [TestType]) ...
				)';
			TEST0(pp,:) = testStat;
		end
	end

	%  get p-value
	if ~strcmp(TestSide, "righttail") , p_left = mean(TEST0 < repmat(TEST1,P,1)) ; end
	if ~strcmp(TestSide, "lefttail") , p_right = mean(TEST0 > repmat(TEST1,P,1)) ; end
	
	if strcmp(TestSide, "righttail"), pvalue = p_right;
	elseif strcmp(TestSide,"lefttail"), pvalue = p_left; 
	else pvalue = min(2*min(p_left,p_right),1) ;
	end

end

