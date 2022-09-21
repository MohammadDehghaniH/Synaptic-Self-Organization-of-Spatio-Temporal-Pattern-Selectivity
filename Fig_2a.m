clc; clear; close all

rng(1101)

Fig2a=figure();

i_t=randi(500); %i'th trial (choose the i'th trial randomly.)

N_psn=1;      %number of post-synaptic neuron
desired_S = 2; % desired number of spikes r0 = 2 Hz

w_max = 1; % maximum amount of synapsis.

resy = 0.0; %resting potential
threshy = 1; % threshold potential

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


inall=zeros(NN,2*TT-1);
tr=200;%transient time
dta=delt/taumem;

intimeP_vec = 5000; %ms Inserting time of the embedded pattern ()

eval(['load Learning_B/Weights/W_',num2str(i_t),'_',num2str(1),  '.mat W_vec']);
W_vec1=W_vec;

eval(['load Learning_BP/Data_Pattern/EP_',num2str(i_t),'_',num2str(1),  '.mat PAT']);
eval(['load Learning_BP/Weights/W_',num2str(i_t),'_',num2str(1),  '.mat W_vec']);

TP = 50/delt; %pattern length (ms/delt= steps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating randomly embedded pattern
PATP=zeros(NN,TT);
PATP(NI+1:NN,:) = rand(400,TT) < (rateE*delt/1000); 
PATP(1:NI,:) = (rand(NI,TT) < rateI*delt/1000);   
PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 

Navg=1; % number of averaging

P_vec=zeros(Navg,TT); % membrane potential
E_vec=zeros(Navg,TT); % Excitatory input
I_vec=zeros(Navg,TT); % Inhibitory input

WW_vec = W_vec;
WW_vec(1:NI,:) = -W_vec(1:NI,:);

for i=1:Navg

    for ti = 1: NI
        inall(ti,:) = conv(PATP(ti,:), traceI);        
    end 
    for ti = NI+1: NN
        inall(ti,:) = conv(PATP(ti,:), traceE);
    end    
    INP = inall(:,1:TT);   
    I_ext(:,1)=WW_vec(:,1)'*INP; 
    INPE=WW_vec(NI+1:NN,1)'*INP(NI+1:NN,:);
    INPI=WW_vec(1:NI,1)'*INP(1:NI,:);
    
    
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
        if V_vec(it-1,:) >= 1
            V_vec(it,1)=resy;
        end
    end   
    
    P_vec(i,:)=V_vec;
    E_vec(i,:)=INPE;
    I_vec(i,:)=INPI;
    
end

xt=(1:1:TT)/10000;
PSP=mean(P_vec,1);
PE=mean(E_vec,1);
PI=mean(I_vec,1);
v=plot(xt,PSP,'g','DisplayName','V(t) ','Linewidth',2);hold on;
i=plot(xt,PE,'r','DisplayName','excitatory input','Linewidth',1.9);hold on;
e=plot(xt,PI,'b','DisplayName','inhibitory input','Linewidth',1.9);
hold on

plot(xt(intimeP_vec:intimeP_vec+TP), PSP(intimeP_vec:intimeP_vec+TP), 'k','Linewidth',2.5);
hold on
plot(xt(intimeP_vec:intimeP_vec+TP), PE(intimeP_vec:intimeP_vec+TP), 'k','Linewidth',2.5);
hold on
plot(xt(intimeP_vec:intimeP_vec+TP), PI(intimeP_vec:intimeP_vec+TP), 'k','Linewidth',2.5);

hleg.String = {'Random Red','Random Blue'};



hold on

xlabel('time (sec)','FontSize', 14)
box off
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;      
set(gca,'FontSize',18, 'color', 'none');


xlim([0.3 0.7])
ylim([min(PI)-0.5 max(PE)+0.5])

legend([v e i ],'Location','northwest' )
title('(a)')
print('Fig2a','-dpng','-r300')



