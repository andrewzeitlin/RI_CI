
rng('default'); % set seed for replicability
addpath('../m'); % assumes we are in the /demo/ folder as pwd

%  Parameters of the simulation
N = 1000; % number of observations
R = 100 ; % number of alternative permutations of the treatment assignment to be used
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
[pval t1 t0] = ri_ci(DATA,{'y'},{'t'},tau0, T0, R ,'ShowMainEstimates',true);
%scatter(y0,y0star) % confirm this recovers the true y0

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
[pval t1 t0] = ri_ci(DATA,{'y'},{'t'},tau0, T0, 100);
%scatter(y0,y0star) % confirm this recovers the true y0

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
[pval ,~,~,~,CI ,Q_UB, Q_LB ] = ri_ci(DATA,{'y'},{'t'},tau0, T0, R ,'FindCI',true);
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

R = 500 ; % number of alternative permutations of the treatment assignment to be usedtau = 1 ;
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

D.Properties.VariableNames

clear rereg
clear ri_ci
[pval,t1,t0] = ri_ci(D,{'y'},{'t'},0, T0, 50,'Model','rereg','GroupVar',{'g'} ); % is the problem here that groupvar needs to be passed to rereg as the name of a variable, not as an actual vector/array?

foo = array2table([1,2;3,4],'RowNames',{'a' 'b'},'VariableNames',{'c' 'd'});
therow = {'a'};
thecolumn = {'d'};
foo(therow,thecolumn)


