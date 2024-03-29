% This script contains all that is needed to recreate the heartbeat figures
% from fig 4, as well as the heartbeat entropy production estimates. 
addpath('Code')
addpath('Data')
%% Get human data
% Data from Goldberger et al. (2000). Not on the github as it is too large,
% but can be obtained from https://physionet.org/content/ltdb/1.0.0/, sample 14046
A = readtable('Data/ecg.txt');

% We do a very rough estimate of the time between beats by thresholding and
% saying X > 2 is a beat. This is only used visualization purposes not for 
% the entropy calculation, which uses the more precise calculation from
% Umetani et al.
T = A.Var1(1:end);
X = A.Var3(1:end);
X = X > 2;

t_A = [];
t0 = 0;

for i = 2:length(T)
    if X(i) ~= X(i-1)
        t = T(i-1) - t0;
        t0 = T(i-1);
        if ~X(i-1)
            t_A(1+end) = t;
        end
    end
end

t_human_scaled = t_A/mean(t_A); % For histograms
T_human = T(520:776) - T(520); % For example ECG trace
X_human = A.Var3(520:776);

%% Make plots of Beat histogram
figure
subplot(1,3,1)
histogram(t_human_scaled,'Normalization','pdf','BinWidth',0.08)
title('Human')
xlim([0,1.7])
ylim([0,4.5])

% Dog and mouse data from PhysioZoo, see Behar et al.

load('Data/peaks_Dog.mat','Data','Fs')
wt = (Data(2:end) - Data(1:end-1))/Fs;
subplot(1,3,2)
histogram(wt/mean(wt),'Normalization','pdf','BinWidth',0.08)
title('Dog')
xlim([0,1.7])
ylim([0,4.5])


load('Data/peaks_Mouse','Data','Fs')
wt = (Data(2:end) - Data(1:end-1))/Fs;
subplot(1,3,3)
histogram(wt/mean(wt),'Normalization','pdf','BinWidth',0.08)
title('Mouse')
xlim([0,1.7])
ylim([0,4.5])


saveas(gcf,'EcgHist.pdf')

%% Make plots of ECG trace

figure

subplot(1,3,1)
plot(T_human,X_human)
title('Human')

load('Data/elctrography_Dog.mat','Data','Fs')
subplot(1,3,2)
plot((1:2*Fs)/(Fs),Data(500+(1:2*Fs)))
title('Dog')

load('Data/elctrography_Mouse','Data','Fs')
subplot(1,3,3)
plot((1:2*Fs)/(Fs),Data(500+(1:2*Fs)))
title('Mouse')

saveas(gcf,'EcgTrace.pdf')

%% Plot fig 4(a) and where the estimates lie on it.
figure
clear
load n_conv_fine
c1 = [217,240,163]/256;
c2 = 0.6*[0,68,27]/256;
for i = 2:19
    hold on
    plot([0;sig_est(:,1)],[1;sig_est(:,i)-1],'-','Color',c1 + (i-2)/17 * (c2-c1));     
end

s = linspace(20,200,400);
plot(s,1 ./s + 4*log(s) ./ (s.^2),':')

s = zeros(65,1);
s1 = zeros(size(s));
s2 = zeros(size(s));
s3 = zeros(size(s));
s4 = zeros(size(s));
for k = 1:size(s,1)
    s(k) = s_extrap(sig_est(k,12:19),12:19);
    s2(k) = s_extrap(sig_est(k,13:19),13:19);
    s3(k) = s_extrap(sig_est(k,12:17),12:17);
    s4(k) = median([s(k)  s2(k) s3(k)]);
end
sig = sig_est(1:size(s,1),1);
hold on
plot(sig,s4-1)



%% Calculate entropy estimates
% Here we use the data directly as measured in Umetani et al. and Behar et
% al. (see main text and SI).

clear

s_asy = @(t2) 1/t2 + 4*log(1/t2);

% Human data directly from table 2 of Umetani et al.
% For 10-19 year old humans
var_A = 81^2; % square of stdev
tA = 1000*60/80; % convert from beats per min to ms
var_A_scaled = var_A/tA^2;
sig = (2/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for 10-19 year old human\n',sig)
hold on
plot([1,s_asy(var_A_scaled)],var_A_scaled *[1 1],'k:')
plot(s_asy(var_A_scaled) *[1 1],[0.001,var_A_scaled ],'k:')

% For 40-49 year old humans
var_A = 60^2; % square of stdev
tA = 1000*60/78; % convert from beats per min to ms
var_A_scaled = var_A/tA^2;
sig = (2/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for 40-49 year old human\n',sig)
hold on
plot([1,s_asy(var_A_scaled)],var_A_scaled *[1 1],'k:')
plot(s_asy(var_A_scaled) *[1 1],[0.001,var_A_scaled ],'k:')

% For dogs
var_A = 69.20^2;
tA = 482.63;
var_A_scaled = var_A/tA^2;
sig = (2/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for dogs\n',sig)
hold on
plot([1,s_asy(var_A_scaled)],var_A_scaled *[1 1],'k:')
plot(s_asy(var_A_scaled) *[1 1],[0.001,var_A_scaled ],'k:')

% For mice
var_A = 10.39^2;
tA = 108.46;
var_A_scaled = var_A/tA^2;
sig = (2/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for mice\n',sig)
hold on
plot([1,s_asy(var_A_scaled)],var_A_scaled *[1 1],'k:')
plot(s_asy(var_A_scaled) *[1 1],[0.001,var_A_scaled ],'k:')

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([1,200])
ylim([0.005,1])
xlabel('Normalized entropy production rate')
ylabel('Waiting time variance')

saveas(gcf,'Fig4a.pdf')
