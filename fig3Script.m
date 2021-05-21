% This script makes plots for the active sensor example
addpath('Code')
addpath('Data')
%% First set up problem
WA_fun = @(wp,wm) [-(wp+wm), wp, 0, 0; ... 
                  wm, -(wp+wm), wp, 0; ... 
                  0, wm, -(wp+wm), wp; ...
                0, 0, wm, -(wm+wp)];
n = 4;
vec1 = ones(n,1);
pi_A = ones(n,1)/(n+1);
ent_rate = @(sig,wm) (1-wm) .*abs(log(abs(1 ./ wm))) - sig;

%% make plots of time spent bound distributions, i.e. fig 3b
sig_vals = [0.,1,2,3,4,5];
Wm = zeros(size(sig_vals));
for i = 1:length(sig_vals)
    f = @(wm) ent_rate(sig_vals(i),wm);
    Wm(i) = fzero(f,[eps 1]);
end

t = linspace(0,80,3200);
f = zeros(length(t),length(Wm));
% know the analytic distribution
for j = 1:length(Wm)
    WA = WA_fun(1,Wm(j));    
    for i = 1:length(t)
        f(i,j) = (pi_A')*(WA*WA* expm(WA*t(i))) * vec1;
    end
    K = - (pi_A')*WA * vec1;
    f(:,j) = f(:,j)/K;
end

subplot(1,2,1)
plot(t,f(:,1),t,f(:,2),t,f(:,3),t,f(:,4),t,f(:,5))
xlim([0,10])
ylim([0,0.5])
ylabel('Time spent bound')
xlabel('Probability density')


%% Make plots 3c. Load the data for the empirical testing as it takes a while to generate.
% To generate the data run script Code/TUR_comp
load sensor_empirical

subplot(1,2,2)
hold on
plot(sig_vals_fine,TUR_med,sig_vals_fine,OPTIM_med,sig_vals_fine,sig_vals_fine)
s2 = [sig_vals_fine, fliplr(sig_vals_fine)];
inBetween = [TUR_3, fliplr(TUR_97)];
fill(s2, inBetween, 'g','FaceAlpha',0.2);


plot(sig_vals_fine,OPTIM_med)
inBetween = [OPTIM_3, fliplr(OPTIM_97)];
fill(s2, inBetween, 'r','FaceAlpha',0.2);
xlabel('Exact entropy production rate')
ylabel('Estimated entropy production rate')

%saveas(gcf,'ActiveSensor.pdf')