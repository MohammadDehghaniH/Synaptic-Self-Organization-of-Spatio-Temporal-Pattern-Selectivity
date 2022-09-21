clc; clear; close all;

Fig8d=figure();

Ntrials=20000; %number of learning cycles
xax=1:1:Ntrials;
N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)
Results1=zeros(N_mu,Ntrials); % total number of spikes
Results2=zeros(N_mu,Ntrials); % ns1+ns2  ns1: pattern 1 ns2: pattern 2
Results3=zeros(N_mu,Ntrials); % R
Figa=figure();

for ic=1:N_mu
   
     try      
       eval(['load Learning_BPB/Data_ns_15/2ER_',num2str(ic),'_',num2str(1),  '.mat ns_15']);
       eval(['load Learning_BPB/Spikes/2spike_p_',num2str(ic),'_',num2str(1),  '.mat Spikess_p']);
    catch       
     end     

    Results2(ic,:)=ns_15(:,1)+ns_15(:,2);
    Results1(ic,:)=Spikess_p;
    
    Results3(ic,:)=Results2(ic,:)./(1e-16+Results1(ic,:));
end

y=mean(Results3);
t_s=50:50:Ntrials;

plot(xax(t_s),y(t_s),'color',[0 0.5 0],'Linewidth',2.5);


h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('R','FontSize',14)
xlabel('Learning cycle','FontSize',14)
box off
title('(d)')
set(gca,'FontSize',14,'color','none')
xlim([0 Ntrials])
ylim([0. 1.05])


print('Fig8d','-dpng','-r300')




