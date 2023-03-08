clc;clear;close all;

Fig12=figure();

Ntrials_2 = 50000; %Ntrials_2: number of learning cycles
N_mu=100; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%line: R=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Data: f_0 = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN=Ntrials_2/10;
xax=(10:10:Ntrials_2);
Results1=zeros(N_mu,NN);

for ic=1:N_mu
    try
       eval(['load Learning_BP_f0/Data_ns_15/ER_',num2str(ic),'_',num2str(1),'_',num2str(10),  '.mat ns_15']);
       eval(['load Learning_BP_f0/Spikes/spike_',num2str(ic),'_',num2str(1),'_',num2str(10),  '.mat Spikess']);
    catch       
    end
    Results1(ic,:)=ns_15(10:10:end)./(1e-16+Spikess(10:10:end));
end

R=mean(Results1);
I50=plot(xax,R,'.','color',[0 0 0.5],'DisplayName','f_{0} = 0.1','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Data: f_0 = 0.02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=Ntrials_2/50;
xax=(50:50:Ntrials_2);
Results2=zeros(N_mu,NN);

for ic=1:N_mu
    try
       eval(['load Learning_BP_f0/Data_ns_15/ER_',num2str(ic),'_',num2str(2),'_',num2str(50),  '.mat ns_15']);
       eval(['load Learning_BP_f0/Spikes/spike_',num2str(ic),'_',num2str(2),'_',num2str(50),  '.mat Spikess']);
    catch       
    end
    Results2(ic,:)=ns_15(50:50:end)./(1e-16+Spikess(50:50:end));
end


R=mean(Results2);
I100=plot(xax,R,'.','color',[0 0.5 0],'DisplayName','f_{0} = 0.02','Linewidth',2.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Data:: f_0 = 0.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=Ntrials_2/100;
xax=(100:100:Ntrials_2);
Results3=zeros(N_mu,NN);

for ic=1:N_mu
    try
       eval(['load Learning_BP_f0/Data_ns_15/ER_',num2str(ic),'_',num2str(3),'_',num2str(100),  '.mat ns_15']);
       eval(['load Learning_BP_f0/Spikes/spike_',num2str(ic),'_',num2str(3),'_',num2str(100),  '.mat Spikess']);
    catch       
    end
    Results3(ic,:)=ns_15(100:100:end)./(1e-16+Spikess(100:100:end));
end


R=mean(Results3);
I300=plot(xax,R,'.','color',[0.5 0 0],'DisplayName','f_{0} = 0.01','Linewidth',2.5);
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
xlim([0 Ntrials_2])

print('Fig12','-dpng','-r300')

