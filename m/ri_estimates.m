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
	addOptional(options,'formula',{});

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
	formula = options.Results.formula; 

	%  Run DGP in reverse to get y0
	t1 = table2array(DATA(:,txvars)); % actual assignment vector(s)
	y0 = table2array(DATA(:,outcome)) - t1 * tau0 ; 
	y0(~isnan(y0)) = min(max(y0(~isnan(y0)), Support(1)),Support(2)) ;  % impose support on non-missing values of y0

	%%  SETUP 
	if strcmp(model, 'lme') 	%  LME:  add y0 to data table, interpret into formula 
		DATA.y0 = y0 ; 
		if ~strcmp(formula,'')
			%  Translate estimmation parameters into formula for lme 
			rhs = [txvars Controls ] ; 
			formula = [ outcome{1} ' ~ 1 ' ] ; 
			for k = 1:length(rhs) 
				formula = [ formula ' + ' rhs{k}] ; 
			end	
			for k = 1:length(groupvar)
				formula = [ formula ' + (1 | ' groupvar{k} ' ) '] ; 
			end
		end
	elseif strcmp(model,'lm') || strcmp(model,'re') 	% Otherwise, for speed, and potential parallelizability: don't swap treatments into data table, but run this as a matrix. 
		x = table2array(DATA(:,Controls)); 
		if strcmp(model,'re')
			g = table2array(DATA(:,groupvar));
		end
	else 
		x = []; % creating empty matrix in this case just so it can be passed to the getTestStat function without worrying about model-specific cases.
	end 

	%5  EXTRACT TEST STAT
	%  Estimate test statistic on the hypothesized y0 using the *actual* assignment
	if length(TestValue) == 0
		if strcmp(model, 'lme') % function call is different for this one because it operates only on tables
			TestValue = lmeTestStat(DATA,formula,TheTx,TestType)
		elseif strcmp(model, 're') % separating this because it requires an extra argument (the group variable) to be passed as an array.
			TestValue = getTestStat(y0,txvars,TheTx,t1,model,TestType,'Controls',Controls,'Support',Support,'x',x,'g',g); 
		else
			TestValue = getTestStat(y0,txvars,TheTx,t1,model,TestType,'Controls',Controls,'Support',Support,'x',x);  
		end
	end

	%  Now, loop over feasible randomizations, impose treatment effect, re-estimate, and extract test statistic
	if RunParallel
		parfor pp = 1 : P 
			%  Draw treatment 
			t0 = permute(T0(:,pp,:),[1 3 2]);  % accommodates possibliity of multiple treatment variables

			if strcmp(model, 're')
				testStat = getTestStat(y0,txvars,TheTx,t0,model,TestType,'Controls',Controls,'Support',Support,'x',x,'g',g); 
			elseif strcmp(model, 'lme')
				%  TODO:  replace treatment variables in memory with values from this permutation
				%  TODO:  Then estimate




			else
				testStat = getTestStat(y0,txvars,TheTx,t0,model,TestType,'Controls',Controls,'Support',Support,'x',x);  
			end
			TEST0(pp,:) = testStat;
		end
	else 
		if ShowWaitBar, hh = waitbar(0, WaitMessage); end
		for pp = 1 : P
			%  Draw treatment 
			t0 = permute(T0(:,pp,:),[1 3 2]);  % accommodates possibliity of multiple treatment variables

			% Using model-specific calls for rereg to avoid creating extra copies of the dataset otherwise 
			if strcmp(model, 're')
				testStat = getTestStat(y0,txvars,TheTx,t0,model,TestType,'Controls',Controls,'Support',Support,'x',x,'g',g); 
			else
				testStat = getTestStat(y0,txvars,TheTx,t0,model,TestType,'Controls',Controls,'Support',Support,'x',x);  
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

function lmeTestStat = lmeTestStat(DATA,formula,TheTx,TestType) 
	lme = fitlme(DATA,formula); 
	[~,~,stats] = fixedEffects(lme) ;  % extract coefficients ('Estimate') and t-stats ('tStat') . 
	%  Check for lme's tendency to rename variables with _1 
	indicators = find(contains(stats.Name,'_1')) ; 
	for ii = indicators 
		stats.Name(ii) = strrep(stats.Name(ii),'_1','');
	end
	stats.Properties.RowNames = stats.Name; % add row names 
	testStat = stats(TheTx,TestType);
end


function testStat = getTestStat(y0,txvars,TheTx,t0,mod el,TestType,varargin)
	%  Parse input 
	parameters = inputParser ; 
	addOptional(parameters,'x',[]); 
	addOptional(parameters,'Controls',{});
	addOptional(parameters,'Support',[-inf,inf]);
	addOptional(parameters,'g',[]);
	parse(parameters,varargin{:}); 
	x = parameters.Results.x;
	Controls = parameters.Results.Controls;
	Support = parameters.Results.Support; 
	g = parameters.Results.g; 	

	if strcmp(model, 'ks')
		if length(txvars) > 1 
			tx = t0(find(strcmp(txvars,TheTx)));
			[~,~,testStat ] = kstest2(y0(tx==1),y0(tx==0)) ; % two-sample KS stat
		else 
			[~,~,testStat ] = kstest2(y0(t0==1),y0(t0==0)) ; % two-sample KS stat
		end
	elseif strcmp(model,'lm')
		lm = fitlm([t0 , x ], y0) ; % 
		testStat = table2array(lm.Coefficients(1+find(strcmp(txvars,TheTx)), [TestType]));
	elseif strcmp(model,'re') 
		result = rereg_array(y0,[t0, x],g,[txvars Controls]) ; 
		testStat = table2array(result(TheTx, TestType));
	else
		error('Must specify a valid model')
	end
end