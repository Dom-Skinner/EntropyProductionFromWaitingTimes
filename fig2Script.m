% Creates the plot for figure 2.
addpath('Code')
addpath('Data')
%% Figure 2(a)
clear
load n_conv_init

subplot(1,3,1)
% Show the first 8 for convergence
for i = 2:8
    hold on
    plot([0;sig_est(:,1)],[2;sig_est(:,i)]-1,'-'); 
end
% Show 18 for the 'converged' value, could show extrapolated value, but
% they amount to essentially the same thing for this scale.
hold on
plot([0;sig_est(:,1)],[2;sig_est(:,9)]-1,'k-'); 
xlim([0,5])

%% Now for the cow data
%% Cow experiment 1
subplot(1,3,3)
probA = readtable('Data/lying_exp_1.txt').Var1;
% probA extracted from fig 1 (b) of Tolkamp et al.
probB = readtable('Data/stand_exp_1.txt').Var1;
% probA extracted from fig 2 (c) of Tolkamp et al. (hence log conversion
% later)
histogram('BinEdges',linspace(0,4,length(probA)+1),'BinCounts',probA,...
    'Normalization','pdf','DisplayStyle', 'stairs')


X_A = linspace(0,4,length(probA)+1);
X_A = 0.5*(X_A(2:end) + X_A(1:end-1));
tA = X_A*probA/sum(probA);
tA2 = (X_A.^2)*probA/sum(probA);

X_B = linspace(0,12,length(probB)+1);
X_B = 0.5*(X_B(2:end) + X_B(1:end-1));
tB = exp(X_B/60)*probB/sum(probB);

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For experiment 1, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)
subplot(1,3,1)
hold on
plot([0,sig*(tA + tB)/2],(tA2/tA^2-1) *[1 1],'k:')
plot(sig*(tA + tB)/2 *[1 1],[tA2/tA^2-1 , 0],'k:')


%% Cow experiment 2
probA = readtable('Data/lying_exp_2.txt').Var1;
probB = readtable('Data/stand_exp_2.txt').Var1;
subplot(1,3,3)
hold on
histogram('BinEdges',linspace(0,4,length(probA)+1),'BinCounts',probA,...
    'Normalization','pdf','DisplayStyle', 'stairs')

X_A = linspace(0,4,length(probA)+1);
X_A = 0.5*(X_A(2:end) + X_A(1:end-1));
tA = X_A*probA/sum(probA);
tA2 = (X_A.^2)*probA/sum(probA);

X_B = linspace(0,12,length(probB)+1);
X_B = 0.5*(X_B(2:end) + X_B(1:end-1));
tB = exp(X_B/60)*probB/sum(probB);

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For experiment 2, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)
subplot(1,3,1)
hold on
plot([0,sig*(tA + tB)/2],(tA2/tA^2 -1) *[1 1],'k:')
plot(sig*(tA + tB)/2 *[1 1],[tA2/tA^2-1 , 0],'k:')

%% Cow experiment 3
probA = readtable('Data/lying_exp_3.txt').Var1;
probB = readtable('Data/stand_exp_3.txt').Var1;
subplot(1,3,3)
hold on
histogram('BinEdges',linspace(0,4,length(probA)+1),'BinCounts',probA,...
    'Normalization','pdf','DisplayStyle', 'stairs')
xlabel('Time spent lying (h)')
ylabel('Probability density')
ylim([0,1.6])
legend({'Indoor beef','Outdoor beef','Dairy'})


X_A = linspace(0,4,length(probA)+1);
X_A = 0.5*(X_A(2:end) + X_A(1:end-1));
tA = X_A*probA/sum(probA);
tA2 = (X_A.^2)*probA/sum(probA);

X_B = linspace(0,12,length(probB)+1);
X_B = 0.5*(X_B(2:end) + X_B(1:end-1));
tB = exp(X_B/60)*probB/sum(probB);

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For experiment 3, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)
subplot(1,3,1)
hold on
plot([0,sig*(tA + tB)/2],(tA2/tA^2 -1) *[1 1],'k:')
plot(sig*(tA + tB)/2 *[1 1],[tA2/tA^2-1 , 0],'k:')

%% And the gene data
%% For GTGlutaminas
X = readtable('Data/GTGlutaminase.txt').Var1;
Y = linspace(0,10,length(X)+1);
Y = 0.5*(Y(2:end) + Y(1:end-1));
subplot(1,3,2)
hold on
histogram('BinEdges',linspace(0,10,length(X)+1),'BinCounts',X,...
    'Normalization','pdf','DisplayStyle', 'stairs')
xlabel('Time spent in off state (h)')
ylabel('Probability density')
tA = Y*X/sum(X);
tA2 = (Y.^2)*X/sum(X);
tB = 5/60; % from median on figure - estimate but not so important since it is so small.
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For GTGlutaminas gene, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)
subplot(1,3,1)
hold on
plot([0,sig*(tA + tB)/2],(tA2/tA^2-1) *[1 1],'k--')
plot(sig*(tA + tB)/2 *[1 1],[tA2/tA^2-1 , 0],'k--')

%% For Bmal1
subplot(1,3,2)
X = readtable('Data/Bmal1a_off.txt').Var1;
Y = linspace(0,10,length(X)+1);
Y = 0.5*(Y(2:end) + Y(1:end-1));
histogram('BinEdges',linspace(0,10,length(X)+1),'BinCounts',X,...
    'Normalization','pdf','DisplayStyle', 'stairs')
legend({'GT:Glutaminase','Bmal1a'})
xlim([0,6])
tA = Y*X/sum(X);
tA2 = (Y.^2)*X/sum(X);
tB = 5/60; % from median on figure - estimate but not so important since it is so small.
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Bmal1 gene, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)
subplot(1,3,1)
hold on
plot([0,sig*(tA + tB)/2],(tA2/tA^2-1) *[1 1],'k--')
plot(sig*(tA + tB)/2 *[1 1],[tA2/tA^2-1 , 0],'k--')
ylim([0.3,1])

saveas(gcf,'fig2.pdf')