function Model_Learning_f0(p1,p2)

%p1: seed for random number and p1=1:1:\mu. There are \mu independent simulations
%p2: Pattern frequency. Every f learning cycle, the embedded pattern is shown.
%p2=1: f_0 = 0.1 ; p2=2: f_0 = 0.01 ; p2=3: f_0 = 0.02

rng(10+(p1-1)*1e6);


Ntrials_1=2000; 
%Ntrials_1: number of learning cycles to learn the background.
%There is no embedded pattern in afferents for 2000 learning cycles. (step 1)

Ntrials_2=50000; 
%Ntrials_2: number of learning cycles to learn the background containing the embedded pattern (step 2).
%(Note initial conditions in this step are from step1)
% There is one embedded pattern in afferents every 1/f_0 steps.


%Spikess: vector of output spikes.
%ns_0: vector of output spikes in embedded pater duration
%ns_15:  vector of output spikes in embedded pater duration + 15ms
%em:  em=0: there is no embedded pattern in afferents. em=1: there is an embedded pattern in afferents.
alpha_1=1e-2; %scaling factor 
C_E=0.9e-4;   %excitatory learning rate
bet_a=0.9e-4; %scaling coefficient
AA=1-bet_a;
C_I = 1e-4;%inhibitory learning rate


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
tautraceE =3;%Rise time of excitatory currents
tautraceI = 5;%Decay time of inhibitory currents
tauriseE =0.5;%Decay time of excitatory currents
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

W_vec = 0.01+1e-3*randn(NN,N_psn); % initial weight vector
W_vec(W_vec<0)=0;

Filterd_Dw=zeros(NN,N_psn); %\tilde{\varepsilon} : eligibility after applying hetero-synaptic plasticity


inall=zeros(NN,2*TT-1);
Filterd_Spikes=zeros(1,N_psn); %long-time firing rate
tr=200;%transient time
dta=delt/taumem;

intimeP_vec = 5000; %ms Inserting time of the embedded pattern ()
TP_vec = [50 100 300];  %ms pattern leangth
TP=TP_vec(1)/delt; 

f0_vec=[10 50 100];
f0=f0_vec(p2);
% generating randomly embedded pattern
PAT1 = rand(NN, TP,2);
PAT1(1:NI,:,:) = 1.*(PAT1(1:NI,:,:) < rateI*delt/1000);
PAT1(NI+1:NN,:,:) = 1.*(PAT1(NI+1:NN,:,:) < rateE*delt/1000);
PAT=PAT1(:,:,1);


V0I=0; %modification threshold for inhibitory neurons
V0E=0; %modification threshold for excitatory neurons

%step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spikess=zeros(Ntrials_1,N_psn);
em=0; %if em=0, there is no embedded pattern
for itrial=1:Ntrials_1

    [I_ext,INP]=Model_Input_f0(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em,itrial,f0);                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron Model
    num_spikes_vec=zeros(N_psn,TT);
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
        if V_vec(it-1,:) >= threshy
            V_vec(it,1)=resy;
            num_spikes_vec(1,it)=1;
        end
    end   
    num_spikes=sum(num_spikes_vec(:,tr:TT),2);
    Filterd_Spikes = 0.9*Filterd_Spikes + 0.1*num_spikes;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Spikess(itrial,:)=num_spikes;

    eligibility_vec_I=INP(:,tr:TT)*(V_vec(tr:TT,:)-V0I);
    eligibility_vec_E=INP(:,tr:TT)*((V_vec(tr:TT,:) - V0E).*((V_vec(tr:TT,:) - V0E) > 0));

    Filterd_Dw=0.99*Filterd_Dw+...
        0.01*[eligibility_vec_I(1:NI,:);C_E.*(eligibility_vec_E(NI+1:NN,:)-mean(eligibility_vec_E(NI+1:NN,:)))];
    W_vec(NI+1:NN)=AA*W_vec(NI+1:NN)*exp(alpha_1*(desired_S-Filterd_Spikes));
    W_vec=[W_vec(1:NI,:) + C_I.*Filterd_Dw(1:NI,:);W_vec(NI+1:NN,:) + Filterd_Dw(NI+1:NN,:)];

    W_vec=W_vec.*(W_vec > 0.);
    W_vec(NI+1:NN)=(W_vec(NI+1:NN).*(W_vec(NI+1:NN) <= w_max) + w_max*(W_vec(NI+1:NN) > w_max));

end
%
eval(['save Learning_B_f0/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Learning_B_f0/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['save Learning_B_f0/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['save Learning_B_f0/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['save Learning_B_f0/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['save Learning_B_f0/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);

%step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns_0=zeros(Ntrials_2,N_psn);
ns_15=zeros(Ntrials_2,N_psn);
Spikess=zeros(Ntrials_2,N_psn);
em=1;
for itrial=1:Ntrials_2

    [I_ext,INP]=Model_Input_f0(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em,itrial,f0);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron Model
    num_spikes_vec=zeros(N_psn,TT);
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
        if V_vec(it-1,:) >= threshy
            V_vec(it,1)=resy;
            num_spikes_vec(1,it)=1;
        end
    end   
    num_spikes=sum(num_spikes_vec(:,tr:TT),2);
    Filterd_Spikes = 0.9*Filterd_Spikes + 0.1*num_spikes;    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ns_0(itrial,:) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1: intimeP_vec(1,1) + TP),2);
    ns_15(itrial,:) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1:intimeP_vec(1,1) + TP+150),2);
    Spikess(itrial,:)=num_spikes;

    eligibility_vec_I=INP(:,tr:TT)*(V_vec(tr:TT,:)-V0I);
    eligibility_vec_E=INP(:,tr:TT)*((V_vec(tr:TT,:) - V0E).*((V_vec(tr:TT,:) - V0E) > 0));

    Filterd_Dw=0.99*Filterd_Dw+...
        0.01*[eligibility_vec_I(1:NI,:);C_E.*(eligibility_vec_E(NI+1:NN,:)-mean(eligibility_vec_E(NI+1:NN,:)))];
    W_vec(NI+1:NN)=AA*W_vec(NI+1:NN)*exp(alpha_1*(desired_S-Filterd_Spikes));
    W_vec=[W_vec(1:NI,:) + C_I.*Filterd_Dw(1:NI,:);W_vec(NI+1:NN,:) + Filterd_Dw(NI+1:NN,:)];

    W_vec=W_vec.*(W_vec > 0.);
    W_vec(NI+1:NN)=(W_vec(NI+1:NN).*(W_vec(NI+1:NN) <= w_max) + w_max*(W_vec(NI+1:NN) > w_max));


end
eval(['save Learning_BP_f0/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat PAT']);
eval(['save Learning_BP_f0/Weights/W_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat W_vec']);
eval(['save Learning_BP_f0/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat ns_0']);
eval(['save Learning_BP_f0/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat ns_15']);
eval(['save Learning_BP_f0/Spikes/spike_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat Spikess']);
eval(['save Learning_BP_f0/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat Filterd_Spikes']);
eval(['save Learning_BP_f0/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat eligibility_vec_I']);
eval(['save Learning_BP_f0/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat eligibility_vec_E']);
eval(['save Learning_BP_f0/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),'_',num2str(f0),  '.mat Filterd_Dw']);


end
