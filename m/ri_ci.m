function varargout = ri_ci(DATA, outcome, txvars, tau0 , T0, P, varargin) % model, stat, varargin

	%  Function to conduct RI, including (optionally) confidence intervaals and test of the no-effect null.

	%  TODO (1): add functionality to allow estimation commands other than -lm-.
	%	- clusterreg() 
	%   - rereg()
	%   - fitlme()
	%  TODO (2):  allow specification of *inner* boundaries for CI search, so that this can be performed more efficiently/with greater accuracy/over a smaller search region.
	%  TODO (3):  allow stepsize-based stopping criterion in search for CIs.
	

	%  Container for outputs:
	varargout = cell(1, nargout);
	
	if size(tau0,1) > size(tau0,2) 
		error('null treatment vector tau0 should be (1 x K) not (K x 1)');
	end

	%  Parse inputs
	params = inputParser ; 
	addOptional(params,'Controls',{}); % observational RHS variables
	addOptional(params,'Model',{'lm'});  			% model type
	addOptional(params,'TestType',{'tStat'});		%  name of test statistic
	addOptional(params,'TestSide',{'twosided'}); 	% test type for the primary test
	addOptional(params,'FindCI',false);  			%  Switch: locate confidence interval
	addOptional(params,'CIguess',[0 0]); 			% initial guess for confidence interval. 
	addOptional(params,'MaxQueries',10); 			% maximum number of trials to find each end of the confidence interval
	addOptional(params,'MinStepSize',0); 			% minimum step size for stopping rule.
	addOptional(params,'SignificanceLevel',0.05); 	% alpha for CI
	addOptional(params,'ShowMainEstimates',false); 	% to replay results of primary estimate.
	addOptional(params,'TestZero',true); 			% to add a test of the zero null. Defaults on.
	addOptional(params,'CheckBoundaries',true); 	% check boundaries for CI esitmation.
	addOptional(params,'GroupVar',{''}); 			% group variable for random effects or clustered estimates
	addOptional(params,'TheTx',{}); 				% for vector-valued treatments, this cell array contains the name of the treatment of interest.
	parse(params,varargin{:}); 

	model = params.Results.Model; 
	TestType = params.Results.TestType; 
	TestSide = params.Results.TestSide;
	xvars = params.Results.Controls;
	TheTX = params.Results.TheTx ; 
	FindCI = params.Results.FindCI ;
	MaxQueries = params.Results.MaxQueries; 
	MinStepSize = params.Results.MinStepSize; 
	SignificanceLevel = params.Results.SignificanceLevel;  
	groupvar = params.Results.GroupVar; 

	%  If treatment assignment is vector-valued, need to decide which variable will be basis for the test for CI purposes.
	if length(txvars)==1 | length(TheTx) == 0 
		tx = table2array(DATA(:,txvars(1))) ; % if scalar-valued treatment or no specific choice specified, assume this is the first argument.
	else 
		tx = table2array(DATA(:,params.Results.TheTx )); 
	end
	if FindCI & length(txvars)> 1
		error('Confidence interval search currently supports only one-dimensional treatment.')
	end


	%  Containers for results 
	TEST0 = NaN(P,length(txvars)); % to hold null distribution for test statistic.

	%  Estimate the model using the actual assignment. Obtain test statistic.  
	if strcmp(model,'lm')
		lm = fitlm(DATA(:,[txvars xvars outcome])) ;
		TEST1 = table2array(lm.Coefficients([txvars],[TestType]))';

		%  Starting values for search for 95% CI
		if FindCI 
			beta = table2array(lm.Coefficients(2,'Estimate')); 
			se   = table2array(lm.Coefficients(2,'SE')); 
		end
	elseif strcmp(model,'rereg') 
		result = rereg(DATA,outcome,[txvars xvars] ,groupvar )
		TEST1 = table2array(result([txvars],[TestType]));
		if FindCI 
			beta = table2array(result(txvars,{'beta'}));
			se = table2array(result(txvars,{'SE'}));
		end
	elseif strcmp(model,'ks')
		% TODO:  Allow residualization of outcome at this stage, then used residualized outcome for all subsequent analysis.
		if length(xvars)> 0 | length(txvars) > 1
			y = table2array(DATA(:,outcome)) ; 
			rhs = [ones(size(y)),table2array(DATA(:,[~strcmp(txvars,tx),Controls]))];
			bb = rhs\y ;        % regression coefficient
			ydd = y - rhs*bb ;  % residuals 
		else 
			ydd = table2array(DATA(:,outcome));
			ydd = ydd - mean(ydd); 
		end
		[~,~,TEST1] = kstest2(ydd(tx==0),ydd(tx==1)) % two-sample KS stat

		%  replace outcome variable in DATA with residualized version
		if params.Results.TestZero | FindCI 
			DATA(:,outcome) = array2table(ydd);
		end

		if FindCI 
			%  use simple linear regression to initialize search region.
			lm = fitlm(tx,ydd);
			beta = table2array(lm.Coefficients(2,'Estimate')); 
			se   = table2array(lm.Coefficients(2,'SE')); 
		end
	end

	%  Conduct RI on the model, using the point estimmate as a point of comparison. Two-sided test, null hypothesis as specified by user (parameter tau0). 
	if params.Results.TestZero 
		[ pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0 ...
			,xvars, model,T0,P ...
			,TestType,TestSide,TEST1 ...
			,'GroupVar',groupvar ...
			) ; 
	end

	if params.Results.ShowMainEstimates 
		sprintf('__RESULTS OF ANALYTICAL MODEL:__')
		lm 
		sprintf('The p-value from Randomization Inference for the hypothesis that tau = %02.2f is %0.2f',tau0, pvalue)
	end
	%  If requested, find confidence interval
	if FindCI 

		%  Find upper boundary of CI. -------------------------------------%
		lb = beta ;
		ub = beta + 10*1.96*se ; 
		middle = (lb + ub) / 2 ; 

		%  Confirm that p-value at upper bound of search region is below significance threshold
		if params.Results.CheckBoundaries 
			[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,ub,xvars, model,T0,P,TestType,'right',TEST1,'GroupVar',groupvar); 
			if p > SignificanceLevel/2
				error('Initial value for upper bound of CI not big enough.')
			end
		end

		QUERIES_UB = NaN(MaxQueries,2); 
		q = 1 ; 
		stepsize = MinStepSize + 1 ;
		while q <= MaxQueries & stepsize >= MinStepSize % TODO: add step-size constraint.
			% fprintf('This is counter number %i \n', q)
			QUERIES_UB(q,1) = middle;
			[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,middle,xvars, model,T0,P,TestType,'right',TEST1,'GroupVar',groupvar); 
			QUERIES_UB(q,2) = p;

			%  If reject, move left.  Otherwise, move right.
			if q > 1, m0 = middle ; end % store previous tested value if appropriate
			if p <= SignificanceLevel/2
				ub = middle;
				middle = (ub + lb ) / 2;  
			else
				lb = middle;
				middle = (ub + lb ) / 2;  
			end
			if q > 1 
				stepsize = abs(m0 - middle) ;
			else 
				stepsize = MinStepSize + 1 ; % don't let this constraint bind on first try
			end 

			%  Update counter
			q = q+1;
		end
		CI_UB = max(QUERIES_UB(QUERIES_UB(:,2) > SignificanceLevel/2, 1)) ; 

		%  Find lower boundary of CI. -------------------------------------%
		ub = beta ;
		lb = beta - 10*1.96*se ; 
		middle = (lb + ub) / 2 ; 

		%  Confirm that p-value at lower bound of search region is below significance threshold
		if params.Results.CheckBoundaries 
			[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,lb, xvars, model,T0,P,TestType,'left',TEST1,'GroupVar',groupvar); 
			if p > SignificanceLevel/2
				error('Initial value for lower bound of CI not low enough.')
			end
		end

		QUERIES_LB = NaN(MaxQueries,2); 
		q = 1 ; 
		stepsize = MinStepSize + 1 ;
		while q <= MaxQueries & stepsize >= MinStepSize 
			%  fprintf('This is counter number %i \n', q)
			QUERIES_LB(q,1) = middle;
			[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,middle, xvars, model,T0,P,TestType,'left',TEST1,'GroupVar',groupvar); 
			QUERIES_LB(q,2) = p;

			%  If reject, move left.  Otherwise, move right.
			if q > 1, m0 = middle ; end % store previous tested value if appropriate
			if p <= SignificanceLevel/2
				lb = middle;
				middle = (ub + lb ) / 2;  
			else
				ub = middle;
				middle = (ub + lb ) / 2;  
			end

			%  Update counter(s)
			if q > 1 
				stepsize = abs(m0 - middle) ;
			else 
				stepsize = MinStepSize + 1 ; % don't let this constraint bind on first query
			end 
			q = q+1;
		end
		CI_LB = min(QUERIES_LB( QUERIES_LB(:,2) > SignificanceLevel/2, 1)); 

		%  Return CI as a vector
		CI = [CI_LB , CI_UB];
	end


	% Write function outputs
	if nargout >= 1 , varargout{1} = pvalue ; end
	if nargout >= 2 , varargout{2} = TEST1 ; end 
	if nargout >= 3 , varargout{3} = TEST0 ; end 
	if nargout >= 4 , varargout{4} = y0 ; end
	if nargout >= 5 , varargout{5} = CI ; end
	if nargout >= 6 , varargout{6} = QUERIES_UB ; end
	if nargout >= 7 , varargout{7} = QUERIES_LB ; end 
end


%  Subfunction to conduct RI for a particular value
function [pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0,xvars, model, T0, P,TestType,TestSide,TEST1,varargin)
	%  p-value for hypothesized sharp null.

	%  Unpack.
	options = inputParser ; 
	addOptional(options,'GroupVar',{});
	addOptional(options,'TheTx',{}) 
	parse(options,varargin{:}); 
	groupvar = options.Results.GroupVar; 
	theTx = options.Results.TheTx ;

	%  Run DGP in reverse to get y0
	y0 = table2array(DATA(:,outcome)) - table2array(DATA(:,txvars)) * tau0' ;

	%  Now, loop over feasible randomizations, impose treatment effect, re-estimate, and extract test statistic
	if strcmp(model,'lm'), x = table2array(DATA(:,xvars)); end % for speed.

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
		elseif strcmp(model,'rereg') 
			DATA.ystar = ystar ; % rereg() syntax requires this to be part of the table.
			DATA(:,txvars) = array2table(permute(T0(:,pp,:),[1 3 2])) ; 
			result = rereg(DATA,{'ystar'},[txvars xvars],groupvar);
			testStat = table2array(result(txvars, TestType))';
		elseif strcmp(model,'ks') 
			%  FOR THE KS TEST, NEED TO SPECIFY ONE-SIDED TESTS HERE. 
			if length(txvars) > 1 
				tx = T0(:,pp,find(strcmp(txvars,theTx)))
				%  TODO:  For the KS test, when multiple treatment variables, need to residualize outcome to remove this as a potential source of bias here.
			else 
				tx = (T0(:,pp,1)); 
			end
			% if strcmp(TestSide,'right')
			% 	[~,~,testStat ] = kstest2(ystar(tx==1),ystar(tx==0),'Tail','larger') ; % two-sample KS stat
			% elseif strcmp(TestSide,'left')
			% 	[~,~,testStat ] = kstest2(ystar(tx==1),ystar(tx==0),'Tail','smaller') ; % two-sample KS stat
			% else 
				[~,~,testStat ] = kstest2(ystar(tx==1),ystar(tx==0)) ; % two-sample KS stat
			% end
		end
		TEST0(pp,:) = testStat;

	end

	%  get p-value
	if strcmp(model,'ks')
		pvalue = mean(TEST0 > TEST1) ; % For KS test, one-sided testing is handled in the evaluation of the test statistic itself.
	else
		if ~strcmp(TestSide(1), 'right') , p_left = mean(TEST0 < repmat(TEST1,P,1)) ; end
		if ~strcmp(TestSide(1), 'left') ,  p_right = mean(TEST0 > repmat(TEST1,P,1)) ; end
		if strcmp(TestSide(1), 'right'), pvalue = p_right;
		elseif strcmp(TestSide(1), 'left'), pvalue = p_left; 
		else pvalue = min(2*min(p_left,p_right),1) ;
		end
	end

end

