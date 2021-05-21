% This function creates the data needed to compare TUR with sig_T for the
% active sensor example that makes up figure 3. It takes a while to run,
% since many trajectories have to be simulated.

% The control parameter of the system is really Wm, but for plotting equal
% spaced sigma will look best, so quick convert from sigma to Wm values
sig_vals_fine = linspace(0.01,5,100);
Wm = zeros(size(sig_vals_fine));
ent_rate = @(sig,wm) (1-wm) .*abs(log(abs(1 ./ wm))) - sig;
for i = 1:length(sig_vals_fine)
    f = @(wm) ent_rate(sig_vals_fine(i),wm);
    Wm(i) = fzero(f,[eps 1]);
end



% For every one of these Wm values, we will find the empirical median,
% as well as the 2.5 and 97.5 percentile values for 95% confidence 
% intervals for both estimators. One 'trial' is made of a generous
% 200 experiments, each of length t=4000. From one of these trials both
% methods will give a bound. This is then repeated ntrials times to find the
% percentile statistics. The larger the ntrials the better the statistics will
% be, but it will take a while to run. In any case, the trend is obvious.
TUR_med = zeros(size(Wm));
TUR_97 = zeros(size(Wm));
TUR_3 = zeros(size(Wm));
OPTIM_med = zeros(size(Wm));
OPTIM_97 = zeros(size(Wm));
OPTIM_3 = zeros(size(Wm));

T = 4000;
ntrials = 50;
nexp = 100; % 200

for i = 1:length(Wm)

qA_var_est = zeros(ntrials,1);
qA_mean= zeros(ntrials,1);
N_emp_mean = zeros(ntrials,1);
t2_norm = zeros(ntrials,1);
tau = zeros(ntrials,1);

for j  = 1:ntrials    
    % For each trial, traj_sample will sample the relevant statistics, from
    % which we can compute both estimators.
    [q_A,N_emp,TA_emp,TB_emp,TA2_emp] = traj_sample(T,Wm(i),nexp);
    
    % The relevant statistics for TUR are the fraction of time spent in A,
    % qA, as well as the variance of this quantity across nexp experiments,
    % together with N, the number of transitions.
    qA_var_est(j) = var(q_A);
    qA_mean(j) = mean(q_A);
    N_emp_mean(j) = mean(N_emp);
    
    % For sigma_T, each experiment measures <t_A>,<t_B>, <t_A^2>, so just average
    % that across experiments
    t2_norm(j) = mean(TA2_emp)/mean(TA_emp)^2;
    tau(j) = mean( TA_emp+TB_emp)/2;
end

% Compute the entropy estimates for each trial
TUR_ests = 8* (qA_mean.^2) .* (1 - qA_mean).^2 ./ (T.*qA_var_est)  - 2*N_emp_mean./T;
OPTIM_ests = gamma_precomp(t2_norm,true)./tau;

% Derive their statistics
OPTIM_med(i) = median(OPTIM_ests);
OPTIM_3(i) = quantile(OPTIM_ests,0.025);
OPTIM_97(i) = quantile(OPTIM_ests,0.975);
TUR_med(i) = median(TUR_ests);
TUR_3(i) = quantile(TUR_ests,0.025);
TUR_97(i) = quantile(TUR_ests,0.975);

end

TUR_med(TUR_med < 0) = 0; % If the estimate is < 0 for entropy production (not possible)
TUR_3(TUR_3 < 0) = 0;     % simply set it as zero.
TUR_97(TUR_97 < 0) = 0;
clearvars -except Wm OPTIM_med OPTIM_3 OPTIM_97 TUR_med TUR_3 TUR_97 sig_vals_fine
save sensor_empirical

