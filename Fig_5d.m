clc; clear; close all

Fig5d=figure();

NN = 500; %number of afferents
N_mu=20; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

Results1=zeros(N_mu,NN); % w_{\it{BPBP}}
Results2=zeros(N_mu,NN); % w_{\it{BP}}


for p1=1:N_mu
   
    try    
        eval(['load Learning_BPBP/Weights/W_',num2str(p1),'_',num2str(1),  '.mat W_vec']);

    catch       
    end   
    W_vec(1:100)=-1*W_vec(1:100);
    Results1=W_vec; 
    
    try        
       eval(['load Learning_BP/Weights/W_',num2str(p1),'_',num2str(1),  '.mat W_vec']);

    catch       
    end   
    W_vec(1:100)=-1*W_vec(1:100);
    Results2=W_vec;

    
    plot(Results1,Results2,'.','color',[0.3 0.3 0.3],'Markersize',1)
    hold on
    
end




h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;  
xlabel('w_{\it{BPBP}}','FontSize',14)
ylabel('w_{\it{BP}}','FontSize',14)
box off
title('(d)')
set(gca,'FontSize',14,'color','none')

print('Fig5d','-dpng','-r300')

