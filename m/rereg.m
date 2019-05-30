function beta_re = rereg(y,x,g)
	%  Function to estimate random-effects model.
	%  Delivers coefficients only.
	%  Assumes random-effect groups are reported in an (N x 1) vector of group indices.
	%
	%  Andrew Zeitlin
	%  This verison:  2018.03.20
	%  
	%  Notation:  model below assumed to be of the form y_{ig} = x_{ig} + u_g + e_{ig}
	%  And we will define v_{ig} = u_g + e_{ig} 
	%
	%  Estimation of idiosyncratic error component will follow Stata: 
	%  https://www.stata.com/manuals13/xtxtreg.pdf#xtxtregMethodsandformulasxtreg%2Cre
	%  TBD:  simple etimator or Swamy-Arora estimator for etimating \hat{\sigma}^2_{uSA}
    
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
    G = length(unique(g)); 
    N = length(y); 
    K = size(x,2); 
    T_g = grpstats(ones(N,1),g,'sum'); % Count of observations in each group
	ybar_g = grpstats(y,g); %  Group mean of y
    xbar_g = grpstats(x,g); %  Group mean of x
    ybarbar = mean(ybar_g); % grand mean, weighting *panels* equally.
    xbarbar = mean(xbar_g); % grand means, weighting *panels* equally
    yhat = demean(y,g); % -ybarbar; % ; % group-wise demeaned outcome, centered on grand mean
    xhat = demean(x,g); % group-wise demeaned RHS variables, centered on grand mean

    
    %  Remove linearly dependent variables from xhat. NB we will not be
    %  looking at coefficients directly so don't care about ordering.
    %  xvaries = ( std(xhat) >= tol );     %  indicator for elements of x that vary within groups
    %     [xmin, xmax] = grpstats(x,g,{'min','max'}); 
    %     xvaries = (sum(xmin ~= xmax) > 0 );
    %     xhat = xhat(:,xvaries); %- repmat(xbarbar(xvaries),N,1) ;
    %  Not convinced we should subtract the above off if we are going to
    %  run a constant-less regression here.  Already demeaned by group.
    %  Case:  *no* within-group variation in RHS variablese.
    if sum(sum(xhat==0)) == size(xhat,1)*size(xhat,2)
        xhat = ones(N,1);
    else
        xhat = licols(xhat,tol); 
        %  Confirm that a constant would not be linearly dependent before adding.
        if rank([ones(N,1),xhat]) > rank(xhat) 
            xhat = [ones(N,1), xhat ];
        end
    end

    %  Remove linearly dependent columns from dataset of group means, xbar_g.
    if size(xbar_g,2) > 1
        xbar_g_li = licols(xbar_g,tol);
    else 
        xbar_g_li = xbar_g;
    end

    %  add a constant if doing so would increase rank
    if rank([xbar_g_li,ones(G,1)]) > rank(xbar_g_li)
        xbar_g_li = [ones(G,1) xbar_g_li];
    end
    
    %  1.  Within estimate
    beta_w = [ ... 
                xhat ...  % constant term  included in this, if appropriate.
            ] \ (yhat) ; 
    sigma2e = sum( (yhat - [ ... ones(N,1),  <-- all variables centered on zero so constant term not necessary
            xhat ...
            ]*beta_w ).^2 ) ...
            /(N - G - size(xhat,2) + 1) ; % using count of x variables that have within-group variation, rather than K, as dof adjustment
    
	%  2.  Between estimate
    beta_b = xbar_g_li \ ybar_g ; 
    SSR_b = sum((ybar_g - xbar_g_li*beta_b).^2); % SSR from between estimate
    Tbar = G / sum( T_g.^(-1)) ; 
    sigma2u = max(0, SSR_b/(G-K) - sigma2e/Tbar); 
    
    %  3.  GLS transform for each observation
    theta_g = 1 - sqrt(sigma2e./(T_g*sigma2u + sigma2e));
    theta_g = repelem(theta_g,T_g,1); 
    
    %  4.  GLS  
    beta_re = [ ... ones(N,1)- theta_g , ...  constant term is now in x.
            x - repmat(theta_g,1,K).*repelem(xbar_g,T_g,1) ...  included regressors
            ] ...
            \ ...
            ( y - theta_g .* repelem(ybar_g,T_g,1)) ; % transformed outcome
    
    %  For ease of conformity in the calling program, if a constant was not supplied, suppress its coefficient from output.
    %     if flag_addedconstant == 1 
    %         beta_re = beta_re(1:end-1); 
    %     end
end

function datahat = demean(data,g) 
    [~, datahat,idx] = aggregate(g, data , @(x) bsxfun(@minus, x, mean(x,1)));
    [~, isrt] = sort(cat(1, idx{:}));
    datahat = cat(1, datahat{:});
    datahat = datahat(isrt,:);    
end

function [Xsub,idx]=licols(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns
% Source:  https://www.mathworks.com/matlabcentral/answers/108835-how-to-get-only-linearly-independent-rows-in-a-matrix-or-to-remove-linear-dependency-b-w-rows-in-a-m

    if ~nnz(X) %X has no non-zeros and hence no independent columns
        Xsub=[]; idx=[];
        return
    end

    if nargin<2, tol=1e-10; end
    [Q, R, E] = qr(X,0); 

    if ~isvector(R)
        diagr = abs(diag(R));
    else
        diagr = R(1);   
    end

    %Rank estimation
    r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation

    idx=sort(E(1:r));

    Xsub=X(:,idx);                      

end

