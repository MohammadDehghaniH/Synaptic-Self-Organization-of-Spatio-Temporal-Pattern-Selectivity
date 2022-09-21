clc; clear; close all

Fig5c=figure();

NN = 500; %number of afferents
N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)


Results1=zeros(N_mu,NN); % w_{\itBP}
Results2=zeros(N_mu,NN); % w_{\itBPB}
Results3=zeros(N_mu,NN); % w_{\itBPBP}
Results4=zeros(N_mu,NN); % w_{\itB}

xax=1:1:NN;

for p1=1:N_mu
   
    try        
       eval(['load Learning_BP/Weights/W_',num2str(p1),'_',num2str(1),  '.mat W_vec']);
    catch       
    end       
    W_vec(1:100)=-1*W_vec(1:100);
    [xxx,yyy]=sort(W_vec); 
    Results1(p1,:)=W_vec(yyy);
    
    try        
        eval(['load Learning_BPB/Weights/W_',num2str(p1),'_',num2str(1),  '.mat W_vec']);
    catch       
    end   
    W_vec(1:100)=-1*W_vec(1:100);   
    Results2(p1,:)=W_vec(yyy); 
    
    try
        eval(['load Learning_BPBP/Weights/W_',num2str(p1),'_',num2str(1),  '.mat W_vec']);
    catch       
    end   
    W_vec(1:100)=-1*W_vec(1:100);
    Results3(p1,:)=W_vec(yyy);     
    
    
    try        
        eval(['load Learning_B/Weights/W_',num2str(p1),'_',num2str(1),  '.mat W_vec']);
    catch       
    end   
    W_vec(1:100)=-1*W_vec(1:100);
    Results4(p1,:)=sort(W_vec);  
    
end

R=mean(Results1);
I1=plot(xax,R,'color',[0 0 0.5],'DisplayName','w_{\itBP}','Linewidth',2);


hold on

R=mean(Results2);
I2=plot(xax,R,'color',[0 0.5 0],'DisplayName','w_{\itBPB}','Linewidth',2);

hold on

R=mean(Results3);
I3=plot(xax,R,'color',[0.5 0 0],'DisplayName','w_{\itBPBP}','Linewidth',2);

hold on

R=mean(Results4);
I4=plot(xax,R,'color',[0.5 0.5 0.5],'DisplayName','w_{\itB}','Linewidth',2);

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
xlabel('Neurons','FontSize',14)
ylabel('<w>','FontSize',14)
box off
title('(c)')
set(gca,'FontSize',14,'color','none')
ylim([-1 1])
legend([I4 I1 I2 I3],'Location','southeast','color','none' )

print('Fig5c','-depsc2','-r300')
