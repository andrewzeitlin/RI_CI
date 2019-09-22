function OUT = leftjoin(A,B,varargin)
	%  LEFT outer join that preserves sort order in left matrix, A.
	%  Specify Keys and other options as in outerjoin
	%  Assumes 'Type' = 'left' so no need to specify that.
	
	%  Default to MergeKeys == true
	if max(strcmp('MergeKeys',varargin))
		mergekeys = {}; 
	else
		mergekeys = {'MergeKeys',true} ;
	end
	varargin = [varargin,mergekeys] ; 

	%  Merge
	[OUT, ia] = outerjoin(A,B,'Type','left',varargin{:}) ; 

	%  Restore original sort order from A
	OUT.ia = ia ; 
	OUT = sortrows(OUT,'ia');
	OUT.ia = []; 
end