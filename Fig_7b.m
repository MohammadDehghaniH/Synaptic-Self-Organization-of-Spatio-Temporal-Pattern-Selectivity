clc; clear; close all

Fig7b=figure();


N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)
Ntrials=5000; %number of learning cycles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R=1
yline(1,'--','LineWidth',2,'HandleVisibility','off');
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Results0=zeros(N_mu,Ntrials);
Results1=zeros(N_mu,Ntrials);
Results3=zeros(N_mu,Ntrials);

p2=1;
for p1=1:N_mu

    try
        eval(['load Noise_Robustness/Learn_Gau_2/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
        eval(['load Noise_Robustness/Learn_Gau_2/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
        eval(['load Noise_Robustness/Learn_Gau_2/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
    catch       
    end
    Results0(p1,:)=ns_15(:,1)./(1e-16+Spikess(:,1));
    Results1(p1,:)=ns_15(:,1);
    Results3(p1,:)=Spikess(:,1);
end

xax=(1:1:Ntrials);
R=mean(Results0);
I1=plot(xax,R,'Color',[0 0.5 0],'DisplayName','R','Linewidth',2.5);

hold on

R=zeros(1,Ntrials);
for i=1:Ntrials
    fg=Results3(:,i);
    lk=find(fg>0);
    R(1,i)=sum(Results0(:,i))/(numel(lk));
end

I2=plot(xax,R,'Color',[0.5 0. 0],'DisplayName','R^*','Linewidth',2.5);

hold on

Results0=zeros(N_mu,Ntrials);
Results1=zeros(N_mu,Ntrials);
Results3=zeros(N_mu,Ntrials);
for p1=1:N_mu
    
    try
        eval(['load Noise_Robustness/Learn_Gau_2/Data_ns_0/ER_0_',num2str(p1),'_',num2str(p2),  '.mat ns_0_1']);
        eval(['load Noise_Robustness/Learn_Gau_2/Data_ns_15/ER_1_',num2str(p1),'_',num2str(p2),  '.mat ns_15_1']);
        eval(['load Noise_Robustness/Learn_Gau_2/Spikes/spike_1_',num2str(p1),'_',num2str(p2),  '.mat Spikess_1']);
    catch       
    end
    Results0(p1,:)=ns_15_1(:,1)./(1e-16+Spikess_1(:,1));
    Results1(p1,:)=ns_15_1(:,1);
    Results3(p1,:)=Spikess_1(:,1);
end

R=zeros(1,Ntrials);
for i=1:Ntrials
    fg=Results3(:,i);
    lk=find(fg>0);
    R(1,i)=sum(Results0(:,i))/(numel(lk));
   
end

I3=plot(xax,R,'Color',[0 0 0.5],'DisplayName','R^* for original pattern','Linewidth',2.5);

legend([I1,I2, I3 ],'Location','southeast','color','none' )
ylabel('R')
xlabel('Learning cycle')
title('(b)')
set(gca,'FontSize',14,'color','none' )
ylim([0 1.1])
xlim([1 3000])
set(gca,'FontName','Helvetp1a','FontSize',14,'FontWeight','normal','box','off','LineWidth', 2);

print('Fig7b','-dpng','-r300')


