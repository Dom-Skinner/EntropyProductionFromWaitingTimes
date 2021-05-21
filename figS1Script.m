% Creates the plot for figure S1.
addpath('Code')
addpath('Data')

clear
load('Data/small_asymp')

for i = 2:5
    hold on
    plot([0;sig_est(:,1)],[2;sig_est(:,i)]);
end
plot([0,2],2 - [0,2]/4,'k--') % This is the asymptotic function
xlabel('Normalized entropy production rate, \sigma')
ylabel('Normalized second moment, \theta')
