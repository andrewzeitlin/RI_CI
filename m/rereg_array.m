function [results,N] = rereg_array(y,x,g,xnames)
	%  core of rereg(). Takes arguments as arrays.

    %  Tolerance for claim that there is no within-group variation in some
    %  characteristic 
    tol = 1e-9; % ToDo:  Allow specifying this as an option.

	%  0.  Collect some useful ingredients
    %  Check if regressors, x, include a constant term. If so, do nothing.
    %  If not, add a constant.
    %  Scratching this: we make have a full-rank case with no constant but
    %  where we do *NOT* want to add a constant.
        %     if sum(mean(x==1)==1) == 0 
        %         x = [x , ones(size(x,1),1)];  % constant is appended to end.
        %         flag_addedconstant = 1; % flag for whether constant has been added.
        %     else 
        %         flag_addedconstant = 0;
        %     end
    
	%  0.1.  Group means of outcome and RHS variables
    %  If adding a constant term to x increases its rank, do so
    %  ToDo:  replace use of rank() function with something more efficient for large matrices.
    N = length(y); 
    if rank([x,ones(N,1)]) > rank(x) 
        x = [ones(N,1),x];
        xnames = [{'Constant'},xnames];
    end
    K = size(x,2); 
    glist = unique(g); 
    G = length(glist); 

    %  Group stats, demeaned outcomes, mean outcomes
    T_g = grpstats(ones(N,1),g,'sum'); % Count of observations in each group
	ybar_g = grpstats(y,g); %  Group mean of y
    xbar_g = grpstats(x,g); %  Group mean of x
    ybarbar = mean(ybar_g); % grand mean, weighting *panels* equally.
    xbarbar = mean(xbar_g); % grand means, weighting *panels* equally
    yddot = demean(y,g); % -ybarbar; % ; % group-wise demeaned outcome, centered on grand mean
    xddot = demean(x,g); % group-wise demeaned RHS variables, centered on grand mean

    
    %  Remove linearly dependent variables from xddot. NB we will not be
    %  looking at coefficients directly so don't care about ordering.
    %  Case:  *no* within-group variation in RHS variablese.
    if ~nnz(xddot)  % Case: xddot is all zeros. Don't bother with licols, it will return an empty matrix. Instead, replace with a constant. 
        xddot = ones(N,1);
    else
        xddot = licols(xddot); 
    end

    %  Remove linearly dependent columns from dataset of group means, xbar_g.
    if size(xbar_g,2) > 1
        xbar_g_li = licols(xbar_g,tol);
    else 
        xbar_g_li = xbar_g;
    end

    %  add a constant if doing so would increase rank
    % if rank([xbar_g_li,ones(G,1)]) > rank(xbar_g_li)
    %     xbar_g_li = [ones(G,1) xbar_g_li];
    % end
    
    %  1.  Within estimate
    beta_w = [ ... 
                xddot ...  %
            ] \ (yddot) ; 
    sigma2e = sum( (yddot - [ ... ones(N,1),  <-- all variables centered on zero so constant term not necessary
            xddot ...
            ]*beta_w ).^2 ) ...
            /(N - G - size(xddot,2) + 1) ; % using count of x variables that have within-group variation, rather than K, as dof adjustment
    
	%  2.  Between estimate
    beta_b = xbar_g_li \ ybar_g ; 
    SSR_b = sum((ybar_g - xbar_g_li*beta_b).^2); % SSR from between estimate
    Tbar = G / sum( T_g.^(-1)) ; 
    sigma2u = max(0, SSR_b/(G-K) - sigma2e/Tbar); 
    
    %  3.  GLS transform for each observation
    theta_g = 1 - sqrt(sigma2e./(T_g*sigma2u + sigma2e));
    theta_g = repelem(theta_g,T_g,1); 
    
    %  4.  GLS  
    x_tr = x - repmat(theta_g,1,size(xbar_g,2)).*repelem(xbar_g,T_g,1);   % included regressors; there is a constant in x.
    y_tr = y - theta_g .* repelem(ybar_g,T_g,1) ; 
    beta_re = x_tr \ y_tr ; 
    
    %  Compute variance-covariance matrix
    %  Fitteed values 
    yhattr = x_tr * beta_re;  
    yhat   = x * beta_re ; 
    % Residuals
    res = y_tr - yhattr;
    % Inverse of x_tr'x_tr
    invXtrXtr = ((x_tr'*x_tr)\eye(K));
    
    % Residual variance
    resdf = G - 1;
    resvar = (res'*res) ./ resdf;
        
    %  Variance-covariance matrix (cluster-robust)
    
    % Compute total sum
    total_sum = 0;
    for gg=1:G;  
        total_sum = total_sum + x_tr(g==glist(gg),:)'*res(g==glist(gg))*res(g==glist(gg))'*x_tr(g==glist(gg),:);
    end
    varcoef = invXtrXtr * total_sum * invXtrXtr; 
    %  Need degrees of freedom correction?
    varcoef = (G/(G-1)) * ((N-1)/(N-K)) .* varcoef;

    %  Standard errors
    stderr = sqrt(diag(varcoef));
    t = beta_re ./ stderr ; 

    %  Export results as table
    results = array2table([beta_re, stderr, t ] ...
        , 'VariableNames',{'Estimate', 'SE', 'tStat'} ...
        ,'RowNames', xnames);


end