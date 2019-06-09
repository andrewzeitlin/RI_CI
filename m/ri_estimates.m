function [pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0,xvars, model, T0, P,varargin)
	%  Subfunction to conduct RI for a particular value
	%  p-value for hypothesized sharp null.

	%  Unpack.
	options = inputParser ; 
	addOptional(options,'GroupVar',{});
	addOptional(options,'TheTx',{}) ;
	addOptional(options,'TestSide','twosided');
	addOptional(options,'TestValue',{});
	addOptional(options,'TestType', {}); 
	parse(options,varargin{:}); 
	groupvar = options.Results.GroupVar; 
	theTx = options.Results.TheTx ;
	TestSide = options.Results.TestSide; 
	TestValue = options.Results.TestValue; 
	TestType = options.Results.TestType; 

	%  Run DGP in reverse to get y0
	y0 = table2array(DATA(:,outcome)) - table2array(DATA(:,txvars)) * tau0' ;

	%  Now, loop over feasible randomizations, impose treatment effect, re-estimate, and extract test statistic
	if strcmp(model,'lm'), x = table2array(DATA(:,xvars)); end % for speed.

	for pp = 1 : P

		%  Impose hypothesized DGP
		ystar = y0 + T0(:,pp,:) * tau0' ; % accommodates possibliity of multiple treatment variables

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
		pvalue = mean(TEST0 > TestValue) ; % For KS test, one-sided testing is handled in the evaluation of the test statistic itself.
	else
		if ~strcmp(TestSide(1), 'right') , p_left = mean(TEST0 < repmat(TestValue,P,1)) ; end
		if ~strcmp(TestSide(1), 'left') ,  p_right = mean(TEST0 > repmat(TestValue,P,1)) ; end
		if strcmp(TestSide(1), 'right'), pvalue = p_right;
		elseif strcmp(TestSide(1), 'left'), pvalue = p_left; 
		else pvalue = min(2*min(p_left,p_right),1) ;
		end
	end

end

