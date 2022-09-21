clc;clear;close all;

Fig9c = figure();

N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

N_trials=5000; %number of learning cycles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_psn=7;      %number of post-synaptic neuron

Results1=zeros(N_trials,N_psn);
Results2=zeros(N_trials,N_psn);

for ic=1:N_mu
    try
        eval(['load Fig9b_data/ns15_',num2str(ic),'_',num2str(N_psn),  '.mat ns_15']);
        eval(['load Fig9b_data/Spikess_',num2str(ic),'_',num2str(N_psn),  '.mat Spikess_p']);
    catch       
    end
    Results2=Results2+(ns_15./(1e-16+Spikess_p));
end

R=Results2/N_mu;

imagesc(R')
colormap(gray)
colorbar()
set(gca,'YDir','normal','FontSize',14,'color','none')
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('Post-synaptic neurons','FontSize',14)
xlabel('Learning cycle','FontSize',14)
title('(c)')
ylim([1 N_psn])
xlim([1 5000])
box off

print('Fig9c','-dpng','-r300')
