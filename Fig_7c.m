clc; clear; close all

Fig7c=figure();

N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p2=1; % embedded pattern: 50ms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Jitter
Results1=zeros(N_mu,1e4);
for p1=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseJit_2/L_NoiseJit_',num2str(p1),'_',num2str(p2),  '.mat NoiseJit']); 
    catch       
    end
    Results1(p1,:)=NoiseJit(2,:)./(1e-16+NoiseJit(3,:));
end

R=mean(Results1);
xax=(1e-4:1e-4:0.06)*1e3;
I50=plot(xax,R(1:numel(xax)),'color',[0 0. 0.],'DisplayName','Jitter','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gaussian
Results2=zeros(N_mu,1e4);
for p1=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseGaus_2/L_NoiseGaus_',num2str(p1),'_',num2str(p2),  '.mat NoiseGaus']);        
    catch       
    end
    Results2(p1,:)=NoiseGaus(2,:)./(1e-16+NoiseGaus(3,:));
end

R=mean(Results2);

I300=plot(xax,R(1:numel(xax)),'color',[0.5 0.5 0.5],'DisplayName','Gaussian','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('\sigma (ms)','FontSize',14)
legend([I50 I300],'Location','southwest','color','none')
title('(c)')
set(gca,'FontSize',14,'color','none')
ylim([0. 1.1])
xlim([1e-4 0.05]*1e3)


print('Fig7c','-depsc2','-r300')

