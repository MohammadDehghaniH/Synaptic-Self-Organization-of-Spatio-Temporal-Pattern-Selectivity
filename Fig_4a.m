clc; clear; close all

Fig4a=figure();

Ntrials_3=20000;
%Ntrials_3: number of learning cycles to learn the background. 
%There is no embedded pattern in afferents.
%(Note initial conditions in this step are from step 2.)

N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

Results1=zeros(N_mu,Ntrials_3);
for p1=1:N_mu  
     try
        eval(['load Learning_BPB/Data_ns_15/ER_',num2str(p1),'_',num2str(1),  '.mat ns_15']);
        eval(['load Learning_BPB/Spikes/spike_p_',num2str(p1),'_',num2str(1),  '.mat Spikess_p']);     
    catch       
    end   
    Results1(p1,:)=ns_15./(1e-16+Spikess_p);
end

R=mean(Results1);

x=1:1:Ntrials_3;
s_t=50:50:Ntrials_3;

xax=x(s_t);

I1=plot(xax,R(s_t),'color',[0 0.5 0],'Linewidth',2.5);

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('Learning cycle','FontSize',14)
box off
title('(a)')
set(gca,'FontSize',14,'color','none')
xlim([0 10000])
ylim([0. 1.05])


print('Fig4a','-dpng','-r300')


