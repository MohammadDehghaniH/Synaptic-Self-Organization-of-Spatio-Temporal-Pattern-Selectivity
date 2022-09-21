clc; clear; close all

Fig3a=figure();

Ntrials_1=2000; 
%Ntrials_1: number of learning cycles to learn the background.
%There is no embedded pattern in afferents for 2000 learning cycles. (step 1)
N_mu=500; %There are 500 independent simulations (number of simulations)
xax = 1:1:Ntrials_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results1=zeros(N_mu,Ntrials_1);

for p1=1:N_mu
    try
        eval(['load Learning_B/Spikes/spike_',num2str(p1),'_',num2str(1),  '.mat Spikess']);
    catch       
    end
    Results1(p1,:)=Spikess;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=mean(Results1);
I50=plot(xax,R,'color',[0 0 0],'Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('(a)')

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('# spikes','FontSize',14)
xlabel('Learning cycle','FontSize',14)
set(gca,'FontSize',14,'color','none')
ylim([0 3])
xlim([1 2000])
box off
print('Fig3a','-dpng','-r300')

