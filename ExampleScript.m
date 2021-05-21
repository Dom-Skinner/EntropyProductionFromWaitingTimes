%% This script shows how the saved data was generated to produce the figures.
% We include here how to run minimizations over the discrete and continuous
% systems. First we demonstrate how easy it is to obtain a bound given
% waiting time statistics from the saved function
addpath('Code')

tA = 1.;   % mean time spent in A, <t>_A
tB = 2.;   % mean time spent in B, <t>_B
tA2 = 1.4; % second moment of time spent in A, <t^2>_A
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
disp(sig) % that simple

%% Example of running trials
sig_est = zeros(10,4);
sig_est(:,1) = linspace(0.1,2,size(sig_est,1));

ntrials = 20; % may need more than 20 trials for convergence, low number for demonstration
hess=true; % use Hessian info for minimization
gs = false; % Don't use global search, instead descend from random start points
rand_init=true;

for nvars = 2:size(sig_est,2)
    for i = 1:size(sig_est,1)
    sig_est = run_trials(ntrials,sig_est,nvars,i,hess,rand_init,gs);
    end
end

for i = 2:size(sig_est,2)
    hold on
    plot([0;sig_est(:,1)],[2;sig_est(:,i)],'-o');
end
plot([0,2],2 - [0,2]/4,'k--') % This is the asymptotic function
xlabel('Normalized entropy production rate, \sigma')
ylabel('Normalized second moment, \theta')

rmpath('Code')
addpath('ContinuousCode')
%% Example of running trials for continuous system
t2 = linspace(1.2,2,12);
sig = zeros(size(t2));
N = 30;
M = 5;
ctol = 1e-6;
otol = 1e-6;
for i = 1:length(t2)
    % Have used global search here, but could run trials if needed
    [~,~,~,~,~,~,fval,~] = min_cont(N,M, t2(i),ctol,otol);
    sig(i) = fval;
end
plot([sig';0],[t2';2])
xlabel('Normalized entropy production rate, \sigma')
ylabel('Normalized second moment, \theta')
% If you want to explore what the solutions might look like
xs = linspace(0,1,400);
real_space = @(fk,x) 2*real(sum( fk(2:end) .* exp( 1j * 2* pi* (1:length(fk)-1)' * x))) + real(fk(1));
[a,b,Dr,Fr,Tf,pf,fval,eflag] = min_cont(N,M, 1.2,1e-6,1e-6);
% For demonstration we plot a, b, and p.
plot(xs,real_space(a,xs),xs,real_space(b,xs),xs,real_space(pf,xs)) 
xlabel('x')
%% Example of the delta-function minimization

t2 = linspace(1.2,2,12);
sig = zeros(size(t2));
ctol = 1e-6;
otol = 1e-6;
x1_val = 0.8; % This could be optimized over, but seems to make little difference
for i = 1:length(t2)
    % Have used global search here, but could run trials if needed
    [~,fval,~] = min_delta(t2(i),x1_val,ctol,otol) ;
    sig(i) = fval;
end
plot([sig';0],[t2';2])
xlabel('Normalized entropy production rate, \sigma')
ylabel('Normalized second moment, \theta')
