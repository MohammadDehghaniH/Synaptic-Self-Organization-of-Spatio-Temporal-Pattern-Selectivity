clc;clear;close all;

Figb = figure();

Ntrials = 5000; %total number of time steps (learning cycles)
N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

N_psn=7;      %number of post-synaptic neurons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yline(1,'--k','Linewidth',2)
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results1=zeros(N_mu,Ntrials); %Omega: Fixed
Results2=zeros(N_mu,Ntrials); %Omega: Random 
Results3=zeros(N_mu,Ntrials); %Random: number of the embedded patterns in each learning cycle

for p1=1:N_mu
    p1;
    try
       eval(['load Fig9b_data/Omega_',num2str(p1),'_',num2str(N_psn),  '.mat Omega_R']);   
    catch       
    end
    
    Results1(p1,:)=Omega_R;  
    
    try
       eval(['load Fig9b_data/ROmega_',num2str(p1),'_',num2str(N_psn),  '.mat Omega_R']);
       eval(['load Fig9b_data/Rnum_emb_vec_',num2str(p1),'_',num2str(N_psn),  '.mat num_emb_vec']);
    catch       
    end
    
    Results2(p1,:)=Omega_R;
    Results3(p1,:)=num_emb_vec;
end


xax = 1:1:Ntrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 R=mean(Results1,1);
 I50=plot(xax,R,'color',[0 0 0],'DisplayName','Fixed','Linewidth',1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ir=1:Ntrials
    R(1,ir)=sum(Results2(:,ir))./numel(find(Results3(:,ir)>0));
end

I300=plot(xax,R,'color',[0.5 0.5 0.5],'DisplayName','Random','Linewidth',1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('\Omega','FontSize',14)
xlabel('Learning cycle','FontSize',14)
legend([I50, I300],'Location','southeast','color','none')
title('(b)')
set(gca,'FontSize',14,'color','none')
ylim([0 1.1])
xlim([1 5000])

print('Fig9b','-dpng','-r300')
