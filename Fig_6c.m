clc; clear; close all

Fig6c=figure();

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
        eval(['load Noise_Robustness/NoiseJit/NoiseJit_',num2str(ic),'_',num2str(1),  '.mat NoiseJit']); 
    catch       
    end
    Results1(ic,:)=NoiseJit(2,:)./(1e-16+NoiseJit(3,:));
end
R=mean(Results1);
xax=(1e-4:1e-4:0.06)*1e3;
I50=plot(xax,R(1:numel(xax)),'color',[0 0 0.5],'DisplayName','L_{em} = 50 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 100 ms
Results2=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseJit/NoiseJit_',num2str(ic),'_',num2str(2),  '.mat NoiseJit']); 
    catch       
    end
    Results2(ic,:)=NoiseJit(2,:)./(1e-16+NoiseJit(3,:));
end
R=mean(Results2);
I300=plot(xax,R(1:numel(xax)),'color',[0 0.5 0],'DisplayName','L_{em} = 100 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_{em} = 300 ms
Results3=zeros(N_mu,1e4);
for ic=1:N_mu
    try
        eval(['load Noise_Robustness/NoiseJit/NoiseJit_',num2str(ic),'_',num2str(3),  '.mat NoiseJit']); 
    catch       
    end
    Results3(ic,:)=NoiseJit(2,:)./(1e-16+NoiseJit(3,:));
end
R=mean(Results3);
I500=plot(xax,R(1:numel(xax)),'color',[0.5 0 0],'DisplayName','L_{em} = 300 ms','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('\sigma (ms)','FontSize',14)
legend([I50 I300 I500 ],'Location','southwest','color','none' )
title('(c)')
set(gca,'FontSize',14,'color','none')
ylim([0. 1.1])
xlim([1e-4 0.025]*1e3)

print('Fig6c','-dpng','-r300')
