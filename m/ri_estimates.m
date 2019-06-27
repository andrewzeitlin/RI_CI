function [pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0, model, T0, P,varargin)
	%  Subfunction to conduct RI for a particular value
	%  p-value for hypothesized sharp null.

	%  Unpack.
	options = inputParser ; 
	addOptional(options,'GroupVar',{});
	addOptional(options,'TheTx',{}) ;
	addOptional(options,'TestSide','twosided');
	addOptional(options,'TestValue',{});
	addOptional(options,'TestType', {}); 
	addOptional(options,'Controls',{});
	parse(options,varargin{:}); 
	groupvar = options.Results.GroupVar; 
	theTx = options.Results.TheTx ;
	TestSide = options.Results.TestSide; 
	TestValue = options.Results.TestValue; 
	TestType = options.Results.TestType;
	Controls = options.Results.Controls;  

	%  Run DGP in reverse to get y0
	y0 = table2array(DATA(:,outcome)) - table2array(DATA(:,txvars)) * tau0 ;

	%  If model is ks and TestValue has not been specified, estimate KS statistic on the hypothesized y0
	if strcmp(model,'ks') && length(TestValue) == 0
		if length(txvars) > 1 
			tx = table2array(DATA(:,txvars(find(strcmp(txvars,theTx)))));
		else 
			tx = table2array(DATA(:,txvars)); 
		end
		[~,~,TestValue] = kstest2(y0(tx==1),y0(tx==0));
	end

	%  Now, loop over feasible randomizations, impose treatment effect, re-estimate, and extract test statistic
	if strcmp(model,'lm'), x = table2array(DATA(:,Controls)); end % for speed.
	for pp = 1 : P

		%  For KS stat, RI based on y0
		if strcmp(model, 'ks')
			if length(txvars) > 1 
				tx = T0(:,pp,find(strcmp(txvars,theTx)))
				%  TODO:  For the KS test, when multiple treatment variables, need to residualize outcome to remove this as a potential source of bias here.
			else 
				tx = (T0(:,pp)); 
			end
			[~,~,testStat ] = kstest2(y0(tx==1),y0(tx==0)) ; % two-sample KS stat
		else 
			%  Impose hypothesized DGP
			t0 = permute(T0(:,pp,:),[1 3 2]);  % accommodates possibliity of multiple treatment variables
			ystar = y0 + t0 * tau0 ;

			%  Estimate model and collect test statistic
			if strcmp(model,'lm')
				lm = fitlm([t0 , x ], ystar) ; % 
				testStat = table2array(lm.Coefficients(1+find(strcmp(txvars,theTx)), [TestType]));

			elseif strcmp(model,'re') 
				DATA.ystar = ystar ; % rereg() syntax requires this to be part of the table.
				DATA(:,txvars) = array2table(t0) ; 
				result = rereg(DATA,{'ystar'},[txvars Controls],groupvar);
				testStat = table2array(result(find(strcmp(txvars,theTx)), TestType));
			else
				error('Must specify a valid model')
			end
		end
		TEST0(pp,:) = testStat;
	end

	%  Calculate left and right tails as required 
	if ~strcmp(TestSide(1), 'right') , p_left = mean(TEST0 < repmat(TestValue,P,1)) ; end
	if ~strcmp(TestSide(1), 'left') ,  p_right = mean(TEST0 > repmat(TestValue,P,1)) ; end

	%  Report the appropriate test statistic
	if strcmp(TestSide(1), 'right'), pvalue = p_right;
	elseif strcmp(TestSide(1), 'left'), pvalue = p_left; 
	else pvalue = min(2*min(p_left,p_right),1) ;
	end

end

