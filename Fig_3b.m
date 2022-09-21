clc; clear; close all

Fig3b = figure();

Ntrials_2 = 10000; 
%Ntrials_2: number of learning cycles to learn the background containing the embedded pattern (step 2).

xax=(1:1:Ntrials_2);
N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%line: R=1
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results1=zeros(N_mu,Ntrials_2);

for p1=1:N_mu
    try      
        eval(['load Learning_BP/Data_Cos/Cos_',num2str(p1),'_',num2str(1),  '.mat Cos_W']);
    catch 
    end
    Results1(p1,:)=Cos_W;
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=mean(Results1);
I50=plot(xax,R,'color',[0.7 0 0],'DisplayName','L = 15 ms','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('cos(\theta)','FontSize',14)
xlabel('Learning cycle','FontSize',14)
set(gca,'FontSize',14,'color','none')
ylim([0 1.1])
xlim([1 10000])
title('(b)')


print('Fig3b','-dpng','-r300')
