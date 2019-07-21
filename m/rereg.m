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

    %  Extract data as matrices from table DATA.
    y = table2array(DATA(:,yvar)) ;
    x = table2array(DATA(:,xvars));
    g = table2array(DATA(:,groupvar)); 

    %  The computation itself
    [results,N] = rereg_array(y,x,g,xvars);
    
end

