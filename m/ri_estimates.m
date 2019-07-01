function [pvalue TEST0 test1 y0 ] = ri_estimates(DATA,outcome,txvars,tau0, model, T0, P,varargin)
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
	addOptional(options,'Support',[-inf,inf]);
	addOptional(options,'RunParallel',false); 
	addOptional(options,'ShowWaitBar',false) ; 
	addOptional(options,'WaitMessage','RI in progress...'); 
	parse(options,varargin{:}); 
	groupvar = options.Results.GroupVar; 
	TheTx = options.Results.TheTx ;
	TestSide = options.Results.TestSide; 
	TestValue = options.Results.TestValue; 
	TestType = options.Results.TestType;
	Controls = options.Results.Controls;  
	Support = sort(options.Results.Support); 
	ShowWaitBar = options.Results.ShowWaitBar; 
	WaitMessage = options.Results.WaitMessage; 
	RunParallel = options.Results.RunParallel; 

	%  Run DGP in reverse to get y0
	y0 = table2array(DATA(:,outcome)) - table2array(DATA(:,txvars)) * tau0 ; 
	y0(~isnan(y0)) = min(max(y0(~isnan(y0)), Support(1)),Support(2)) ;  % impose support on non-missing values of y0

	%  If model is ks and TestValue has not been specified, estimate KS statistic on the hypothesized y0
	if strcmp(model,'ks') && length(TestValue) == 0
		if length(txvars) > 1 
			tx = table2array(DATA(:,txvars(find(strcmp(txvars,TheTx)))));
		else 
			tx = table2array(DATA(:,txvars)); 
		end
		[~,~,TestValue] = kstest2(y0(tx==1),y0(tx==0));
	end

	%  Now, loop over feasible randomizations, impose treatment effect, re-estimate, and extract test statistic

	% For speed, and potential parallelizability: don't swap treatments into data table, but run this as a matrix. 
	%  Extract only once..
	if strcmp(model,'lm')
		x = [ ... %ones(size(DATA,1),1) , 
			table2array(DATA(:,Controls)) ...
			]; 
	else 
		x = []; % creating empty matrix in this case just so it can be passed to the getTestStat function without worrying about model-specific cases.
	end 
	
	if RunParallel
		parfor pp = 1 : P 
			if strcmp(model, 'rereg')
				testStat = getTestStat(y0,tau0,txvars,TheTx,model,TestType,T0,pp,'Controls',Controls,'Support',Support,'x',x,'DATA',DATA); 
			else
				testStat = getTestStat(y0,tau0,txvars,TheTx,model,TestType,T0,pp,'Controls',Controls,'Support',Support,'x',x);  
			end
			TEST0(pp,:) = testStat;
		end
	else 
		if ShowWaitBar, hh = waitbar(0, WaitMessage); end
		for pp = 1 : P
			% Using model-specific calls for rereg to avoid creating extra copies of the dataset otherwise 
			if strcmp(model, 'rereg')
				testStat = getTestStat(y0,tau0,txvars,TheTx,model,TestType,T0,pp,'Controls',Controls,'Support',Support,'x',x,'DATA',DATA); 
			else
				testStat = getTestStat(y0,tau0,txvars,TheTx,model,TestType,T0,pp,'Controls',Controls,'Support',Support,'x',x);  
			end
			TEST0(pp,:) = testStat;
			if ShowWaitBar, waitbar(pp/P); end 
		end
		if ShowWaitBar, close(hh) ; end 
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


function testStat = getTestStat(y0,tau0,txvars,TheTx,model,TestType,T0,pp,varargin)
	%  Parse input 
	parameters = inputParser ; 
	addOptional(parameters,'x',[]); 
	addOptional(parameters,'Controls',{});
	addOptional(parameters,'Support',[-inf,inf]);
	addOptional(parameters,'DATA',{});
	addOptional(parameters,'groupvar',{});
	parse(parameters,varargin{:}); 
	x = parameters.Results.x;
	Controls = parameters.Results.Controls;
	Support = parameters.Results.Support; 
	DATA = parameters.Results.DATA;
	groupvar = parameters.Results.groupvar; 	

	%  For KS stat, RI based on y0
	if strcmp(model, 'ks')
		if length(txvars) > 1 
			tx = T0(:,pp,find(strcmp(txvars,TheTx)))
			%  TODO:  For the KS test, when multiple treatment variables, need to residualize outcome to remove this as a potential source of bias here.
		else 
			tx = (T0(:,pp)); 
		end
		[~,~,testStat ] = kstest2(y0(tx==1),y0(tx==0)) ; % two-sample KS stat
	else 
		%  Impose hypothesized DGP
		t0 = permute(T0(:,pp,:),[1 3 2]);  % accommodates possibliity of multiple treatment variables

		ystar = y0 + t0 * tau0 ; 
		ystar(~isnan(ystar)) = min(Support(2),max(Support(1), ystar(~isnan(ystar))));

		%  Estimate model and collect test statistic
		if strcmp(model,'lm')
			lm = fitlm([t0 , x ], ystar) ; % 
			testStat = table2array(lm.Coefficients(1+find(strcmp(txvars,TheTx)), [TestType]));
		elseif strcmp(model,'re') 
			DATA.ystar = ystar ; % rereg() syntax requires this to be part of the table.
			DATA(:,txvars) = array2table(t0) ; 
			result = rereg(DATA,{'ystar'},[txvars Controls],groupvar);
			testStat = table2array(result(find(strcmp(txvars,TheTx)), TestType));
		else
			error('Must specify a valid model')
		end
	end
end