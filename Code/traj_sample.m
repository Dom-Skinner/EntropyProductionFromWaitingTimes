function [q_A,N_emp,TA_emp,TB_emp,TA2_emp] = traj_sample(T,Wm,nexp)
% Simulate the waiting times of the active sensor for a window T. return
% the empirically observed qA. Do this nexp times


% To simulate, we figure out the distribution of f_A, f_B (this one is just 
% exponential) from what we know analytically, and then we just draw from 
% this repeatedly. In retrospect it might be faster to simulate the full 
% system and then figure out the observed dynamics, but this way works.
n = 4;
vec1 = ones(n,1);
pi_A = ones(n,1)/(n+1);

WA_fun = @(wp,wm) [-(wp+wm), wp, 0, 0; ... 
                  wm, -(wp+wm), wp, 0; ... 
                  0, wm, -(wp+wm), wp; ...
                0, 0, wm, -(wm+wp)];
            
t = linspace(0,40,2000); % This is going to be the t we define the distribution over
f = zeros(length(t),1);
f_cs = zeros(length(t),1);

% Define the partial transition rate matrix WA, and calculate the
% cumalitive probability density from analytic formula
WA = WA_fun(1,Wm);
K = - (pi_A')*WA * vec1;

for i = 1:length(t)
    f(i) = (pi_A')*(WA*WA* expm(WA*t(i))) * vec1/K;
    f_cs(i) = (pi_A')*(WA*expm(WA*t(i)) - WA) * vec1/K;
end

% sample from wait time distribution by sampling uniform dist and transform
f_cs_inv = @(x_in) interp1(f_cs,t,x_in,'pchip'); 

q_A = zeros(nexp,1);
N_emp = zeros(nexp,1);
TA_emp = zeros(nexp,1);
TA2_emp = zeros(nexp,1);
TB_emp = zeros(nexp,1);
for n = 1:nexp
    X = rand(round(T) + 1e4,1);
    TA_rand = f_cs_inv(X);
    TB_rand =  exprnd(1/(1 + Wm),size(X)); % Sample from f_A,f_B many times

    TA_glob = zeros(2*length(TA_rand),1);
    TA_glob(1:2:end) = TA_rand;
    TB_glob = zeros(2*length(TB_rand),1);
    TB_glob(2:2:end) = TB_rand;
    T_tot = cumsum(TA_glob+TB_glob); % Construct the times of all jumps.



    T_start = 1e4; % Start a while later to avoid any initial effects.
    T_end = T_start + T;
    idx_start = find(T_tot > T_start,1,'first'); % Find number of jumps in the time interval
    idx_end = find(T_tot > T_end,1,'first');

    if isempty(idx_end)
        error("Trajectory window too long")
    end

    % For TUR
    q_A(n) = (sum(TA_glob(idx_start+1:idx_end)) + ...
        (TA_glob(idx_start) > 0) * (T_tot(idx_start) - T_start) -...
        (TA_glob(idx_end) >0) * (T_tot(idx_end) - T_end)  )/T; % The exact fraction of time spent in A
    N_emp(n) = idx_end - idx_start; % Number of jumps observed.
    
    % For sig_T
    TA_glob = TA_glob(idx_start+1: idx_end-1);
    TB_glob = TB_glob(idx_start+1: idx_end-1);
    TA_emp(n) = mean(TA_glob(TA_glob>0)); 
    TB_emp(n) = mean(TB_glob(TB_glob>0)); 
    TA2_emp(n) = mean(TA_glob(TA_glob>0).^2); 
end
end