clc;clear;close all;

Fig9a = figure();

N_mu=50; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

%R1:  with pre-synaptic competition
%R(1,1): 4 post-synaptic neurons
%R(1,2): 5 post-synaptic neurons
%R(1,3): 6 post-synaptic neurons
%R(1,4): 7 post-synaptic neurons
R1=zeros(1,4);
R1e=zeros(1,4);  %errorbar for R1

%R2: without pre-synaptic competition
%R(1,1): 4 post-synaptic neurons
%R(1,2): 5 post-synaptic neurons
%R(1,3): 6 post-synaptic neurons
%R(1,4): 7 post-synaptic neurons
R2=zeros(1,4);
R2e=zeros(1,4); %errorbar for R2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results1=zeros(N_mu,1);
Results2=zeros(N_mu,1);
N_psn=4;      %number of post-synaptic neuron
for ic=1:N_mu

    try
        eval(['load ',num2str(N_psn,'%.i') '/With_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end
    Results1(ic,1)=mean(Omega_R(end-100:end));
    
    try
        eval(['load ',num2str(N_psn,'%.i') '/Without_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end    
    
    Results2(ic,1)=mean(Omega_R(end-100:end));
  
end

R1(1,1)=mean(Results1);
R2(1,1)=mean(Results2);
R1e(1,1)=std(Results1)/sqrt(N_mu);
R2e(1,1)=std(Results2)/sqrt(N_mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results1=zeros(N_mu,1);
Results2=zeros(N_mu,1);
N_psn=5;%number of post-synaptic neuron
for ic=1:N_mu

    try
        eval(['load ',num2str(N_psn,'%.i') '/With_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end
    Results1(ic,1)=Omega_R(end);
    
    try
        eval(['load ',num2str(N_psn,'%.i') '/Without_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end    
    
    Results2(ic,1)=Omega_R(end);
  
end

R1(1,2)=mean(Results1);
R2(1,2)=mean(Results2);
R1e(1,2)=std(Results1)/sqrt(N_mu);
R2e(1,2)=std(Results2)/sqrt(N_mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results1=zeros(N_mu,1);
Results2=zeros(N_mu,1);
N_psn=6;%number of post-synaptic neuron
for ic=1:N_mu
    try
        eval(['load ',num2str(N_psn,'%.i') '/With_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end
    Results1(ic,1)=Omega_R(end);  
    try
        eval(['load ',num2str(N_psn,'%.i') '/Without_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end      
    Results2(ic,1)=Omega_R(end);
  
end

R1(1,3)=mean(Results1);
R2(1,3)=mean(Results2);
R1e(1,3)=std(Results1)/sqrt(N_mu);
R2e(1,3)=std(Results2)/sqrt(N_mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results1=zeros(N_mu,1);
Results2=zeros(N_mu,1);
N_psn=7;%number of post-synaptic neuron
for ic=1:N_mu

    try
        eval(['load ',num2str(N_psn,'%.i') '/With_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end
    Results1(ic,1)=Omega_R(end);
    
    try
        eval(['load ',num2str(N_psn,'%.i') '/Without_Omega_',num2str(ic),'_',num2str(N_psn),  '.mat Omega_R']);
    catch       
    end    
    Results2(ic,1)=Omega_R(end);
end

R1(1,4)=mean(Results1);
R2(1,4)=mean(Results2);
R1e(1,4)=std(Results1)/sqrt(N_mu);
R2e(1,4)=std(Results2)/sqrt(N_mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes1 = axes('Parent',Fig9a,...
    'Position',[0.13 0.14 0.8 0.8]);
hold(axes1,'on');
yline(1,'Parent',axes1,'LineWidth',3);
hold on
a=0.15;
x=[4,5,6,7];
bar1=bar (x-a, R2,0.14,'FaceColor','r','EdgeColor','r','LineWidth',.5)  ;              
hold on
er = errorbar (x-a, R2, R2e, R2e,'LineWidth',1.6); 
er.Color = [0 0.5 0]; 
er.LineStyle = 'none' ; 
hold on
bar2=bar (x+a, R1,0.14,'FaceColor','b','EdgeColor','b','LineWidth',.5)  ;              
hold on
er = errorbar (x+a, R1, R1e, R1e,'LineWidth',1.6);
er.Color = [0 0.5 0]; 
er.LineStyle = 'none' ;
set(bar1,'DisplayName','without pre-synaptic competition',...
    'FaceColor',[0.5 0 0],...
    'EdgeColor',[0.5 0 0]);
set(bar2,'DisplayName','with pre-synaptic competition',...
    'FaceColor',[0 0 0.5],...
    'EdgeColor',[0 0 0.5]);


% Create label
xlabel('number of post-synaptic neurons');
ylabel('\Omega');
title('(a)')
legend([bar1 bar2],{'without pre-synaptic competition','with pre-synaptic competition'},'Location','northwest','color','none')
set(axes1,'FontSize',14,'LineWidth',2,'XTick',[4 5 6 7],'color','none');
legend1 = legend(axes1,'show','color','none');
ylim([0 1.3])
set(legend1,'FontSize',12,'color','none');

print('Fig9a','-dpng','-r300')

