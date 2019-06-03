function [yvar xvars] = parseWilkinson(tbl,mdl)
	%  Function to decode wilkinson notation
	model = split(mdl,'~')
	if length(model) > 2
		error('There should only be one instance of the character ~ in the model.')
	end
	yvar = model{1}
	xvars = model{2} 
end