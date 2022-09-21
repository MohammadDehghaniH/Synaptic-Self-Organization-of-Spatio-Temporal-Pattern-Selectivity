clc; clear; close all

Fig6b=figure();


N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 50 ms
Results1=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseAdd/NoiseAdd_',num2str(ic),'_',num2str(1),  '.mat NoiseAdd']); 
    catch       
    end
    Results1(ic,:)=NoiseAdd(2,:)./(1e-16+NoiseAdd(3,:));
end
R=mean(Results1);
xax50=5+(0.1:0.1:150);
xax50=1.00*(xax50-0.0)/5;
I50=plot(xax50,R(1:numel(xax50)),'color',[0 0 0.5],'DisplayName','L_{em} = 50 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 300 ms
Results2=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseAdd/NoiseAdd_',num2str(ic),'_',num2str(2),  '.mat NoiseAdd']); 
    catch       
    end
    Results2(ic,:)=NoiseAdd(2,:)./(1e-16+NoiseAdd(3,:));
end
R=mean(Results2);
xax300=5+(0.1:0.1:150);
xax300=1.00*(xax300-0.0)/5;
I300=plot(xax300,R(1:numel(xax50)),'color',[0 0.5 0],'DisplayName','L_{em} = 300 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 300 ms
Results3=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseAdd/NoiseAdd_',num2str(ic),'_',num2str(3),  '.mat NoiseAdd']); 
    catch       
    end
    Results3(ic,:)=NoiseAdd(2,:)./(1e-16+NoiseAdd(3,:));
end
R=mean(Results3);
xax500=5+(0.1:0.1:150);
xax500=1.00*(xax500-0.0)/5;
I500=plot(xax500,R(1:numel(xax50)),'color',[0.5 0 0],'DisplayName','L_{em} = 300 ms','Linewidth',2.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('r^{*}','FontSize',14)
legend([I50 I300 I500 ],'Location','southwest','color','none' )
title('(b)')
set(gca,'FontSize',14,'color','none')
ylim([0. 1.1])
xlim([1. 2.5])

print('Fig6b','-dpng','-r300')


