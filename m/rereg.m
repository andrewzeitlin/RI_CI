function [results,N] = rereg(DATA,yvar,xvars,groupvar)
	%  Function to estimate random-effects model.
	%  Assumes random-effect groups are reported in an (N x 1) vector of group indices.
	%
	%  Andrew Zeitlin
	%  
	%  Notation:  model below assumed to be of the form y_{ig} = x_{ig} + u_g + e_{ig}
	%  And we will define v_{ig} = u_g + e_{ig} 
	%
	%  Estimation of idiosyncratic error component will follow Stata: 
	%  https://www.stata.com/manuals13/xtxtreg.pdf#xtxtregMethodsandformulasxtreg%2Cre
	%  TBD:  simple etimator or Swamy-Arora estimator for etimating \hat{\sigma}^2_{uSA}
    %  Variance estimation follows Inmaculada C. Álvarez, Javier Barbero, José L. Zofío, http://www.paneldatatoolbox.com
    
    %  Check for missing values, and subset the non-missing data that will be used for estimation
    if sum(max(ismissing(DATA(:,[yvar xvars groupvar])),[],2)) > 0 
        DATA = DATA(...
                min(...
                    ~ismissing(DATA(:,[yvar xvars groupvar])) ...
                    ,[],2 ...
                )...
            ,:) ;
    end

    %  Check for categorical variables in xvars, and expand/convert to indicators
    xvars_temp = xvars ; 
    for k = 1:length(xvars)
        xvar = xvars(k);
        if iscategorical(DATA{:,xvar}) 
            Dx = array2table(dummyvar(DATA{:,xvar})); % dummy variables 
            for ll = 1:size(Dx,2)
                Dx.Properties.VariableNames(ll) = strcat(xvar, '_', num2str(ll));
            end
            Dx = Dx(:,2:end);  %  Leave one category out.
            DATA = [DATA,Dx];  %  Append these indicators to the table.

            %  Now remove xvar from Controls and add VariableNames from table Dx instead.
            xvars_temp = xvars_temp(~strcmp(xvars_temp,xvar));
            xvars_temp = [xvars_temp, Dx.Properties.VariableNames];
        end
    end
    xvars = xvars_temp ; 

    %  Extract data as matrices from table DATA.
    y = table2array(DATA(:,yvar)) ;
    x = table2array(DATA(:,xvars));
    g = table2array(DATA(:,groupvar)); 

    %  The computation itself
    [results,N] = rereg_array(y,x,g,xvars);
    
end

