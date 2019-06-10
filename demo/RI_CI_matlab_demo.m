
rng('default'); % set seed for replicability
addpath('../m'); % assumes we are in the /demo/ folder as pwd

%  Parameters of the simulation
N = 1000; % number of observations
R = 500 ; % number of alternative permutations of the treatment assignment to be used
tau = 1 ; % treatment effect

%  Generate data 
x = randn(N,1);
y0 = x + randn(N,1);
y1 = y0 + tau;
t = (rand(N,1) >= 0.5 ) ; % treatment status
y = y0 + t.*(y1 - y0) ; % switching regression
T0 = (rand(N,R) >= 0.5) ; % set of potential randomizations
DATA = array2table([y,t,x] , 'VariableNames',{'y','t','x'});

clear ri_ci
tau0 = 0;
[pval, ~, t1, t0] = ri_ci(DATA,{'y'},{'t'},tau0, T0, R ,'ShowMainEstimates',true);

figure(1)
clf
hax = axes;
hold on
ksdensity(t0) % distribution of test statistic under the null.
line([t1 t1],get(hax,'YLim'),'Color','red'); % [0 1])
legend('Distribution under null','Estimated test statistic')
hold off

clear ri_ci
leftprob = ri_ci(DATA,{'y'},{'t'},0, T0, 100,'TestSide','lefttail')
rightprob = ri_ci(DATA,{'y'},{'t'},0, T0, 100,'TestSide','righttail')

clear ri_ci
tau0 = 1;
[pval ~ t1 t0] = ri_ci(DATA,{'y'},{'t'},tau0, T0, 100);

figure(2)
clf
hax = axes;
hold on
ksdensity(t0) % distribution of test statistic under the null.
line([t1 t1],get(hax,'YLim'),'Color','red'); % [0 1])
legend('Distribution under null','Estimated test statistic')
hold off

clear ri_ci
tau0 = 0;
tic
[pval ,CI,~,~,Q_UB, Q_LB ] = ri_ci(DATA,{'y'},{'t'},tau0, T0, R ,'FindCI',true);
toc

pval
CI 

[Q_UB, Q_LB]

figure(3)
clf
hax = axes;
hold on
scatter(Q_UB(:,1),Q_UB(:,2))
scatter(Q_LB(:,1),Q_LB(:,2))
line([CI(1) CI(1)],get(hax,'YLim'),'Color','red'); % [0 1])
line([CI(2) CI(2)],get(hax,'YLim'),'Color','red'); % [0 1])
hold off




%--------------------------------------------------------------------------%
%%  DEMO FOR OTHER TYPES OF MODELS 
%--------------------------------------------------------------------------%

R = 100 ; % number of alternative permutations of the treatment assignment to be usedtau = 1 ;

const = 10 ; 
G=200;  % number of clusters
n=5;    % number of observations/cluster
t = (rand(G,1) >= 0.5 ) ; % treatment status
T0 = (rand(G,R) >= 0.5) ; % set of potential randomizations
x_g = randn(G,1) ; % observable at cluster level
e_0g = randn(G,1) ; % error term at cluster level
g = [1:G]';  % group index

%  Expand all group-level datasets to the individual level
t = kron(t,ones(n,1));
T0 = kron(T0,ones(n,1));
e_0g = kron(e_0g,ones(n,1));
x_g = kron(x_g,ones(n,1));
g = kron(g,ones(n,1));

%  constant additive treatment effect
tau = 1; 

%  remainder of DGP at individual level
e_0i = randn(G*n,1);
x_i = randn(G*n,1);
y0 = const + x_g + x_i + e_0g + e_0i ;
y1 = y0 + tau ;
y = y0 + tau * t;

%  Data to table, for passing to rereg.
D = array2table([y,t,x_i,x_g,g],'VariableNames',{'y' 't' 'x_i' 'x_g' 'g'});

clear rereg
lm = fitlm(D,'y~ t + x_i + x_g') % ([t x_g x_i ],y)
%beta_re = rereg(y,[t , x_g, x_i],g)  % <- ols syntax:  assumes data are in arrays rather than tables.
re = rereg(D,{'y'},{'t' 'x_i' 'x_g'},{'g'})


%---  RE regression ---%
clear rereg
clear ri_ci
tic
pval = ri_ci(D,{'y'},{'t'}, T0, R,'Model','rereg','GroupVar',{'g'} )
toc

clear ri_ci
tic
[pval,CI,~,~,Q_UB, Q_LB ] = ri_ci(D,{'y'},{'t'}, T0, R ...
    , 'Model', 'rereg', 'GroupVar',{'g'} ...
    ,'FindCI',true ...
    , 'TestZero', false ... % don't bother with p-value for tau=0
    );
toc

CI
[Q_UB, Q_LB]
figure(4)
clf
hax = axes;
hold on
scatter(Q_UB(:,1),Q_UB(:,2))
scatter(Q_LB(:,1),Q_LB(:,2))
line([CI(1) CI(1)],get(hax,'YLim'),'Color','red'); % [0 1])
line([CI(2) CI(2)],get(hax,'YLim'),'Color','red'); % [0 1])
hold off


%--- KS TEST  ---%
clear ri_ci
tic
[pval ] = ri_ci(D,{'y'},{'t'}, T0, R ...
    , 'Model', 'ks' ...
    , 'TestZero', true ... 
    , 'FindCI', false ...
    )
toc

clear ri_ci
tau0 = 0 ;
tic
[pval, CI,~,~,Q_UB,Q_LB] = ri_ci(D,{'y'},{'t'}, T0, R ...
    , 'Model', 'ks' ...
    , 'TestZero', true ... % don't bother with p-value for tau=0
    , 'FindCI', true ...
    );  % ,~,~,~,CI ,Q_UB, Q_LB 
toc

[Q_UB, Q_LB]

figure(5)
clf
hax = axes;
hold on
scatter(Q_UB(:,1),Q_UB(:,2))
scatter(Q_LB(:,1),Q_LB(:,2))
line([CI(1) CI(1)],get(hax,'YLim'),'Color','red'); % [0 1])
line([CI(2) CI(2)],get(hax,'YLim'),'Color','red'); % [0 1])
hold off



% a little demo to amke sure we are correctly interpreting the KS stat
P = 100;
KS0 = NaN(P,1);
for pp = 1:P 
    [~,~,ks] = kstest2( ...
        y(T0(:,pp)==0) ...
        , ...
        y(T0(:,pp) ==1) ...
        );
    KS0(pp,1) = ks;
end
[~,~,ks] = kstest2(y(t==1),y(t==0))
figure(9)
clf
ksdensity(KS0)

%  Understanding why is the KS stat so large for randomizations that are not the true randomization?
figure(10)
clf 
hold on 
h1 = cdfplot(y(T0(:,11)==1))
h2 = cdfplot(y(T0(:,11)==0))
hold off 

[pval,~,test1,test0] = ri_ci(D,{'y'},{'t'},T0,100,'Model',{'ks'}); 

figure(99)
clf
ksdensity(test0)
%  that seems fine
%  So why are we failing to reject the UB?    

lm = fitlm(D.t,D.y)
beta = table2array(lm.Coefficients(2,'Estimate'))
se = table2array(lm.Coefficients(2,'SE'))
ub = beta + 10*1.96*se ; 

%  This is finding an upper bound
TRIALS = [beta:0.01:2*beta]';
PVALS = NaN(size(TRIALS));
KS0 = NaN(length(TRIALS),100);

[~,~,ks1_ub] = kstest2(table2array(D(D.t==1,{'y'})),table2array(D(D.t==0,{'y'})),'Tail','smaller')

clear ri_estimates
for tt = 1:length(TRIALS)
    [p test0] = ri_estimates(D,{'y'},{'t'},TRIALS(tt),{},{'oks_geq'},T0,100,'TestValue',ks1_ub ... %,'TestSide','right'
        ); 
    PVALS(tt) = p;
    KS0(tt,:) = test0';
end

figure(97)
clf 
hax = axes;
scatter(TRIALS,PVALS)
line([beta beta],get(hax,'YLim'),'Color','red'); % [0 1])


%  This is finding a lower bound
TRIALS_LB = [beta:-0.05:beta - 3*beta]';
PVALS_LB = NaN(size(TRIALS_LB)); 
[~,~,ks1] = kstest2(table2array(D(D.t==1,{'y'})),table2array(D(D.t==0,{'y'})),'Tail','larger')
for tt = 1:length(TRIALS)
    p  = ri_estimates(D,{'y'},{'t'},TRIALS(tt),{},'ks',T0,100,'TestValue',ks1,'TestSide','left'); 
    PVALS_LB(tt) = p;
end

figure(98)
clf 
scatter(TRIALS_LB,PVALS_LB)


%%  Let's just see if we can generate a distribution for the KS statistic under the null
[~,~,ks] = kstest2(y(t==1),y(t==0),'Tail','smaller');;

tau0 = 1.5 ;
y0hat = y - t*tau0; 
KS = NaN(P,1);
for pp = 1:P
    ystar = y0hat + T0(:,pp)*tau0 ; % using the true y0, since this is what we would get if we use the true tau anyway.
    [~,~,KS(pp)] = kstest2(ystar(T0(:,pp)==1),ystar(T0(:,pp)==0),'Tail','smaller');
end 

figure(11)
clf
hax = axes;
ksdensity(KS)
line([ks ks],get(hax,'YLim'),'Color','red'); % [0 1])




