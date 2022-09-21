clc; clear; close all

Fig6a=figure();

N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 50 ms
Results1=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseRemove/NoiseRemove_',num2str(ic),'_',num2str(1),  '.mat NoiseRemove']); 
    catch       
    end
    Results1(ic,:)=NoiseRemove(2,:)./(1e-16+NoiseRemove(3,:));
end
R=mean(Results1);
xax50=linspace(1,200,200)/2;
I50=plot(xax50,R(1:200),'color',[0 0 0.5],'DisplayName','L_{em} = 50 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 100 ms
Results2=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseRemove/NoiseRemove_',num2str(ic),'_',num2str(2),  '.mat NoiseRemove']); 
    catch       
    end
    Results2(ic,:)=NoiseRemove(2,:)./(1e-16+NoiseRemove(3,:));
end
R=mean(Results2);
xax100=linspace(1,400,400)/4;
I100=plot(xax100,R(1:400),'color',[0 0.5 0],'DisplayName','L_{em} = 100 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 300 ms
Results3=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseRemove/NoiseRemove_',num2str(ic),'_',num2str(3),  '.mat NoiseRemove']); 
    catch       
    end
    Results3(ic,:)=NoiseRemove(2,:)./(1e-16+NoiseRemove(3,:));
end
R=mean(Results3);
xax300=linspace(1,1200,1200)/12;
I300=plot(xax300,R(1:1200),'color',[0.5 0 0],'DisplayName','L_{em} = 300 ms','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('percentage of removed spikes','FontSize',14)
legend([I50 I100 I300 ],'Location','southwest','color','none' )
title('(a)')
set(gca,'FontSize',14,'color','none')
ylim([0. 1.1])
xlim([0 25])

print('Fig6a','-dpng','-r300')

