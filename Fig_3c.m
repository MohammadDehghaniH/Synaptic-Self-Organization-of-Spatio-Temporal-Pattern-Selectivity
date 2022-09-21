clc; clear; close all

Fig3c = figure();

Ntrials_2 = 10000; 
%Ntrials_2: number of learning cycles to learn the background containing the embedded pattern (step 2).

xax = 1:1:Ntrials_2;
N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%line: R=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Results1=zeros(N_mu,Ntrials_2);
Results2=zeros(N_mu,Ntrials_2);

for p1=1:N_mu
    try
       eval(['load Learning_BP/Data_ns_0/ER_',num2str(p1),'_',num2str(1),  '.mat ns_0']);
       eval(['load Learning_BP/Data_ns_15/ER_',num2str(p1),'_',num2str(1),  '.mat ns_15']);
       eval(['load Learning_BP/Spikes/spike_',num2str(p1),'_',num2str(1),  '.mat Spikess']);
    catch       
    end
    
    Results1(p1,:)=ns_15./(1e-16+Spikess);  
    Results2(p1,:)= ns_0./(1e-16+Spikess);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1=mean(Results1);
I15=plot(xax,R1,'color',[0 0 0.5],'DisplayName','L = 15 ms','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R2=mean(Results2);
I0=plot(xax,R2,'color',[0.85 0 0],'DisplayName','L = 0','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('Learning cycle','FontSize',14)
legend([I15 I0],'Location','southeast','color','none' )
title('(c)')
set(gca,'FontSize',14,'color','none')
ylim([0 1.1])
xlim([1 2000])

print('Fig3c','-dpng','-r300')

