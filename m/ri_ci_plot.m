function ri_ci_plot(y0,y1,beta,CI,varargin)
	%fig = 

	%  Unpack
	options = inputParser ;
	addOptional(options,'Colors',{'b' 'r'});
	addOptional(options,'Boundaries',[-inf, inf]); % boundary values for y.
	addOptional(options,'Alpha',0.2) ; % transparency level	
	addOptional(options,'Figure',[]) ; % allow putting this on an existing figure
	addOptional(options,'LegendLabels',{});
	addOptional(options,'LegendLocation','nw');
	parse(options,varargin{:});
	Colors = options.Results.Colors;
	Boundaries = sort(options.Results.Boundaries);
	Alpha = options.Results.Alpha;
	Figure = options.Results.Figure; 
	LegendLabels = options.Results.LegendLabels;
	LegendLocation = options.Results.LegendLocation;


	%  How much to add or subtract
	CI = sort(CI);
	lb = beta - CI(1);
	ub = CI(2) - beta;

	%  Construct figure 
	if Figure > 0 
		figure(Figure);
	else 
		figure; 
	end

	% hold on 
	[F0,x0] = ecdf(y0); 
	[F1,x1] = ecdf(y1);
	plot(x0,F0,'Color',Colors{1})
	plot(x1,F1,'Color',Colors{2})
	x1_lb = max( x1 - lb , ones(size(x1))*Boundaries(1)) ;
	x1_ub = min(ones(size(x1))*Boundaries(2) , x1 + ub)  ;
	%  TODO:  improve interpolation at edges.  
	%  Specifically, don't taper to corners, but follow slope of the primary estimates until they hit the boundaries [0,1].

	patch([x1_lb ; flip(x1_ub)], [F1;flip(F1)],Colors{2},'LineStyle','none')
	alpha(Alpha)

	if length(LegendLabels) > 0 
		legend(LegendLabels,'Location',LegendLocation);
	end

	% hold off
	
end