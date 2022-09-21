clc; clear; close all

Fig4b=figure();

Ntrials_3=20000;
%Ntrials_3: number of learning cycles to learn the background. 
%There is no embedded pattern in afferents.
%(Note initial conditions in this step are from step 2)

N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

Results1=zeros(N_mu,Ntrials_3); % there is no embedded pattern
Results2=zeros(N_mu,Ntrials_3); % there is one embedded pattern

for p1=1:N_mu
    try        
        eval(['load Learning_BPB/Spikes/spike_',num2str(p1),'_',num2str(1),  '.mat Spikess']);
    catch       
    end   
    Results1(p1,:)=Spikess; 
    try        
        eval(['load Learning_BPB/Spikes/spike_p_',num2str(p1),'_',num2str(1),  '.mat Spikess_p']);     
    catch       
    end   
    Results2(p1,:)=Spikess_p;
end

x=1:1:Ntrials_3;
po=50:50:Ntrials_3;
xax=x(po);


R2=mean(Results2);
I1=plot(xax,R2(po),'color',[0 0.5 0],'DisplayName','embedded pattern','Linewidth',2.5);

hold on

R1=mean(Results1);
I2=plot(xax,R1(po),'color',[0.5 0 0],'DisplayName','no  pattern','Linewidth',2.5);

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('Number of Spikes','FontSize',14)
xlabel('Learning cycle','FontSize',14)
legend([I1 I2],'Location','northwest','color','none' )
box off
title('(b)')
set(gca,'FontSize',14,'color','none')
ylim([0 20])
xlim([1 10000])

print('Fig4b','-dpng','-r300')
