clc; clear; close all

Fig3d=figure();

Ntrials_2 = 10000; 
%Ntrials_2: number of learning cycles to learn the background containing the embedded pattern (step 2).

xax=(1:1:Ntrials_2);
N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%line: R=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Data: 50 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results1=zeros(N_mu,Ntrials_2);
for p1=1:N_mu
    try
       eval(['load Learning_BP/Data_ns_15/ER_',num2str(p1),'_',num2str(1),  '.mat ns_15']);
       eval(['load Learning_BP/Spikes/spike_',num2str(p1),'_',num2str(1),  '.mat Spikess']);
    catch       
    end
    Results1(p1,:)=ns_15./(1e-16+Spikess);
end

R=mean(Results1);
I50=plot(xax,R,'color',[0 0 0.5],'DisplayName','L_{em} = 50 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Data: 100 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results2=zeros(N_mu,Ntrials_2);
for p1=1:N_mu
    try
       eval(['load Learning_BP/Data_ns_15/ER_',num2str(p1),'_',num2str(2),  '.mat ns_15']);
       eval(['load Learning_BP/Spikes/spike_',num2str(p1),'_',num2str(2),  '.mat Spikess']);
    catch       
    end
    Results2(p1,:)=ns_15./(1e-16+Spikess);
end


R=mean(Results2);
I100=plot(xax,R,'color',[0 0.5 0],'DisplayName','L_{em} = 100 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Data: 300 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results3=zeros(N_mu,Ntrials_2);
for p1=1:N_mu
    try
       eval(['load Learning_BP/Data_ns_15/ER_',num2str(p1),'_',num2str(3),  '.mat ns_15']);
       eval(['load Learning_BP/Spikes/spike_',num2str(p1),'_',num2str(3),  '.mat Spikess']);
    catch       
    end
    Results3(p1,:)=ns_15./(1e-16+Spikess);
end


R=mean(Results3);
I300=plot(xax,R,'color',[0.5 0 0],'DisplayName','L_{em} = 300 ms','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('Learning cycle','FontSize',14)
legend([I50 I100 I300 ],'Location','southeast','color','none' )
title('(d)')
set(gca,'FontSize',14,'color','none')
ylim([0 1.1])
xlim([0 2000])

print('Fig3d','-dpng','-r300')

