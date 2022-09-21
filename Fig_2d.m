clc; clear; close all

Fig2d = figure();
%xcorr function is used.
%xcorr(INPI,INPE,'normalized');
%INPI: inhibitory input INPE: Excitatory input

N_psn=1;      %number of post-synaptic neuron
N_mu=20; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

NN = 500; %number of afferents
percI = 20; % inhibitory percentages
NI = fix((percI/100)*NN); % number of inhibitory neurons
T = 1000; %total time [all in ms]
delt = 0.1; % integration step.
TT = T/delt; 
times = 1:TT; % discretization
times = times*delt;


% Generating kernels for excitatory and inhibitory inputs
taumem = 15; %membrane time constant
tautraceE =3;%Rise time of exctatoy currents
tautraceI = 5;%Decay time of inhibitory currents
tauriseE =0.5;%Decay time of excitatoy currents
tauriseI = 1;%Rise time of inhibitory currents
traceE = exp(-times/tautraceE) - exp(-times/tauriseE);
ettaE=tautraceE/tauriseE;
VnormE = (ettaE.^(ettaE/(ettaE - 1)))/(ettaE-1);
traceE = traceE*VnormE;
traceI = exp(-times/tautraceI) - exp(-times/tauriseI);
ettaI=tautraceI/tauriseI;
VnormI = (ettaI.^(ettaI/(ettaI - 1)))/(ettaI-1);
traceI = traceI*VnormI;


%rates
rateE = 5; %in [Hz] rate of excitatory neuron
rateI = 20 ;%in [Hz] rate of inhibitory neuron

intimeP_vec=5000; %ms Inserting time of the embedded pattern ()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TP=50/delt; %pattern length (ms/delt= steps)
Results50=zeros(N_mu,2*TP+1);
for i_t=1:N_mu
  
    try
        eval(['load Learning_BP/Data_Pattern/EP_',num2str(i_t),'_',num2str(1),  '.mat PAT']);
        eval(['load Learning_BP/Weights/W_',num2str(i_t),'_',num2str(1),  '.mat W_vec']);
    catch       
    end
    
    PATR(NI+1:NN,:) = rand(400,TT) < (rateE*delt/1000); 
    inall=zeros(NN,2*TT-1);
    PATR(1:NI,:) = (rand(NI,TT) < rateI*delt/1000);
   
    PATP = PATR;
    PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
    WW_vec = W_vec;
    WW_vec(1:NI,:) = -W_vec(1:NI,:);
    for i = 1: NI
        inall(i,:) = conv(PATP(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    INPE=WW_vec(101:500,1)'*INP(101:500,:);
    INPI=WW_vec(1:100,1)'*INP(1:100,:);

    INPI=INPI(intimeP_vec:intimeP_vec+TP)-mean(INPI(intimeP_vec:intimeP_vec+TP));
    INPE=INPE(intimeP_vec:intimeP_vec+TP)-mean(INPE(intimeP_vec:intimeP_vec+TP));
    Results50(i_t,:)=xcorr(INPI,INPE,'normalized');
    
end
xax=(-TP:1:TP)/10;
R=mean(Results50);
hold on
I50=plot(xax,R,'color',[0 0 0.5],'DisplayName','L_{em} = 50 ms ','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TP=100/delt; %pattern length (ms/delt= steps)
Results100=zeros(N_mu,2*TP+1);
for i_t=1:N_mu

    try
        eval(['load Learning_BP/Data_Pattern/EP_',num2str(i_t),'_',num2str(2),  '.mat PAT']);
        eval(['load Learning_BP/Weights/W_',num2str(i_t),'_',num2str(2),  '.mat W_vec']);
    catch       
    end
    
    PATR(NI+1:NN,:) = rand(400,TT) < (rateE*delt/1000); 
    inall=zeros(NN,2*TT-1);
    PATR(1:NI,:) = (rand(NI,TT) < rateI*delt/1000);
   
    PATP = PATR;
    PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
    WW_vec = W_vec;
    WW_vec(1:NI,:) = -W_vec(1:NI,:);
    for i = 1: NI
        inall(i,:) = conv(PATP(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    INPE=WW_vec(101:500,1)'*INP(101:500,:);
    INPI=WW_vec(1:100,1)'*INP(1:100,:);

    INPI=INPI(intimeP_vec:intimeP_vec+TP)-mean(INPI(intimeP_vec:intimeP_vec+TP));
    INPE=INPE(intimeP_vec:intimeP_vec+TP)-mean(INPE(intimeP_vec:intimeP_vec+TP));
    Results100(i_t,:)=xcorr(INPI,INPE,'normalized');
    
end
xax=(-TP:1:TP)/10;
R=mean(Results100);
hold on
I100=plot(xax,R,'color',[0 0.5 0],'DisplayName','L_{em} = 100 ms ','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TP=300/delt; %pattern length (ms/delt= steps)
Results300=zeros(N_mu,2*TP+1);
for i_t=1:N_mu
    
    try
        eval(['load Learning_BP/Data_Pattern/EP_',num2str(i_t),'_',num2str(3),  '.mat PAT']);
        eval(['load Learning_BP/Weights/W_',num2str(i_t),'_',num2str(3),  '.mat W_vec']);
    catch       
    end
    
    PATR(NI+1:NN,:) = rand(400,TT) < (rateE*delt/1000); 
    inall=zeros(NN,2*TT-1);
    PATR(1:NI,:) = (rand(NI,TT) < rateI*delt/1000);
   
    PATP = PATR;
    PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
    WW_vec = W_vec;
    WW_vec(1:NI,:) = -W_vec(1:NI,:);
    for i = 1: NI
        inall(i,:) = conv(PATP(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    INPE=WW_vec(101:500,1)'*INP(101:500,:);
    INPI=WW_vec(1:100,1)'*INP(1:100,:);
    
    INPI=INPI(intimeP_vec:intimeP_vec+TP)-mean(INPI(intimeP_vec:intimeP_vec+TP));
    INPE=INPE(intimeP_vec:intimeP_vec+TP)-mean(INPE(intimeP_vec:intimeP_vec+TP));
    Results300(i_t,:)=xcorr(INPI,INPE,'normalized');
    
end
xax=(-TP:1:TP)/10;
R=mean(Results300);
hold on
I300=plot(xax,R,'color',[0.5 0 0],'DisplayName','L_{em} = 300 ms ','Linewidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xlabel('Lag (ms)','FontSize', 14)
ylabel('Cross Correlation','FontSize', 14)
box off
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;      
set(gca,'FontSize',24,'color','none');
ylim([-1 0.5])
xlim([-200 200])
legend([I50 I100 I300 ],'Location','southwest','color','none')
title('(d)')

print('Fig2d','-dpng','-r300')


