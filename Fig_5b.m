clc; clear; close all

Fig5b=figure();

Ntrials_3=20000;
%Ntrials_3: number of learning cycles to learn the background. 
%There is no embedded pattern in afferents.
%(Note initial conditions in this step are from step 2)

N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)
Results1=zeros(N_mu,Ntrials_3);

for p1=1:N_mu
    try
        eval(['load Learning_BPB/Data_Cos/Cos_',num2str(p1),'_',num2str(1),  '.mat Cos_W']);
    catch       
    end
    Results1(p1,:)=Cos_W;
end

xax=1:1:Ntrials_3;

R=mean(Results1);
I50=plot(xax,R,'color',[0 0.6 0],'Linewidth',2.5);


h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('cos(\theta)','FontSize',14)
xlabel('Learning cycle','FontSize',14)
title('(b)')
set(gca,'FontSize',14,'color','none')
ylim([0 1])
box off
print('Fig5b','-dpng','-r300')


