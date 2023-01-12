clc; clear; close all

Fig5e = figure();

Ntrials_5=10000; 
%Ntrials_2: number of learning cycles to learn the background containing a new embedded pattern (step 5).
%(Note initial conditions in this step come from step2)

N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line y=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Results1=zeros(N_mu,Ntrials_5);
Results2=zeros(N_mu,Ntrials_5);
p2=1;
for p1=1:N_mu
    try
        eval(['load Learning_BPP/Data_ns_15/ER2_',num2str(p1),'_',num2str(p2),  '.mat ns_15_0']);
        eval(['load Learning_BPP/Spikes/spike2_',num2str(p1),'_',num2str(p2),  '.mat Spikess_0']);
    catch       
    end
    
    Results1(p1,:)=ns_15_0./(1e-16+Spikess_0);
  
    
     try
        eval(['load Learning_BPP/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
        eval(['load Learning_BPP/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
    catch       
     end  
    
    Results2(p1,:)=ns_15./(1e-16+Spikess);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xax = 1:1:Ntrials_5;


R=mean(Results1);
I50=plot(xax,R,'color',[0.5 0 0],'DisplayName','First Pattern','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=mean(Results2);
I300=plot(xax,R,'color',[0 0 0.5],'DisplayName','Second Pattern','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('Learning cycle','FontSize',14)
legend([I50 I300  ],'Location','southeast','color','none' )
title('(e)')
set(gca,'FontSize',14,'color','none')
ylim([0 1.1])
xlim([1 4000])

print('Fig5e','-dpng','-r300')


