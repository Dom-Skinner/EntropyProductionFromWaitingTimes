% Creates the plot for figure 2.
addpath('Code')
addpath('Data')

hold on

% Results from full numerical minimization
load cont_sys
plot(sig_est(:,13),sig_est(:,1)-1,'-')

% results from delta-function minimization.
load delta_sweep
plot(sig_est(:,2),sig_est(:,1)-1,'-')

% Asymptotics from analytical calculation
s = linspace(12,52,400);
plot(s,1 ./s,'k:' )
plot(s,1 ./s + 4*log(s) ./ (s.^2),'k:')

xlim([15,50])
ylim([0,0.1])
xlabel('(Normalized) Entropy production rate, \sigma')
ylabel('Waiting time variance (Var t_A)')

saveas(gcf,'AsymptoticComparison.pdf')