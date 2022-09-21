function Model_Learning(p1,p2)


%p1: seed for random number and p1=1:1:\mu. 
%There are \mu independent simulations (500)
%p2: pattern length. p2=1 is for 50 ms;  p2=2 is for 100 ms; p2=3 is for 300 ms

rng(10+(p1-1)*1e6);

%There are 5 steps in this program.  
%%Step 1: There is no embedded pattern in afferents for 2000 learning cycles.

%%Step 2: There is one embedded pattern in afferents for 10000 learning cycles.

%%Step 3: There is no embedded pattern in afferents for 20000 learning cycles.

%%Step 4: There is the embedded pattern (from step 2) in afferents for 20000 learning cycles.

%%step 5: There is a new embedded pattern in afferents for 10000 learning cycles.
%%(initial conditions are from step 2)

%Data from this program is used in other programs to plot figures. 

Ntrials_1=2000; 
%Ntrials_1: number of learning cycles to learn the background.
%There is no embedded pattern in afferents for 2000 learning cycles. (step 1)

Ntrials_2=10000; 
%Ntrials_2: number of learning cycles to learn the background containing the embedded pattern (step 2).
%(Note initial conditions in this step are from step1)

Ntrials_3=20000; 
%Ntrials_3: number of learning cycles to learn the background. 
%There is no embedded pattern in afferents.
%(Note initial conditions in this step are from step 2)

Ntrials_4=20000; 
%Ntrials_4: number of learning cycles to learn the background containing the embedded pattern (step 4).
%(Note initial conditions in this step are from step 3)


Ntrials_5=10000; 
%Ntrials_2: number of learning cycles to learn the background containing a new embedded pattern (step 5).
%(Note initial conditions in this step are from step2)


%Spikess: vector of output spikes.
%ns_0: vector of output spikes in embedded pattern duration
%ns_15:  vector of output spikes in embedded pater duration + 15ms
%Norm_W: the norm of weight vector in each learning cycle
%Cos_W: cosine between weight vector in each learning cycle with initial weight vector
%em:  em=0: there is no embedded pattern in afferents. em=1: there is an embedded pattern in afferents.
%w =ab, w: W_vec, a=1, b: W_vecb, W_vec=W_vecb


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


Filterd_Spikes=zeros(1,N_psn);%long-time firing rate
tr=200;%transient time
dta=delt/taumem;

intimeP_vec = 5000; %ms Inserting time of the embedded pattern
TP_vec = [50 100 300];  %ms pattern leangth
TP=TP_vec(p2)/delt; 

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
inall=zeros(NN,2*TT-1);
Spikess=zeros(Ntrials_1,N_psn);
em=0;
for itrial=1:Ntrials_1

    [I_ext,INP] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em);
                
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
eval(['save Learning_B/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Learning_B/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['save Learning_B/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['save Learning_B/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['save Learning_B/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['save Learning_B/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);

%step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns_0=zeros(Ntrials_2,N_psn);
ns_15=zeros(Ntrials_2,N_psn);
Spikess=zeros(Ntrials_2,N_psn);
Norm_W=zeros(Ntrials_2,N_psn);
Cos_W=zeros(Ntrials_2,N_psn);
W_vec_0=W_vec;
em=1;
for itrial=1:Ntrials_2

    [I_ext,INP] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em);
    
    
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

    Norm_W(itrial,1)=norm(W_vec);
    Cos_W(itrial,1)=sum(W_vec.*W_vec_0)/(norm(W_vec)*norm(W_vec_0));

end
eval(['save Learning_BP/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),  '.mat PAT']);
eval(['save Learning_BP/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Learning_BP/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
eval(['save Learning_BP/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
eval(['save Learning_BP/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['save Learning_BP/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['save Learning_BP/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['save Learning_BP/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['save Learning_BP/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);
eval(['save Learning_BP/Data_Cos/Cos_',num2str(p1),'_',num2str(p2),  '.mat Cos_W']);
eval(['save Learning_BP/Data_Norm/Norm_',num2str(p1),'_',num2str(p2),  '.mat Norm_W']);

%step 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns_0=zeros(Ntrials_3,N_psn);
ns_15=zeros(Ntrials_3,N_psn);
Spikess=zeros(Ntrials_3,N_psn);
Spikess_p=zeros(Ntrials_3,N_psn);
Norm_W=zeros(Ntrials_3,N_psn);
Cos_W=zeros(Ntrials_3,N_psn);
W_vec_0=W_vec;

for itrial=1:Ntrials_3
    em=0;
    [I_ext,INP] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em);
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

    Spikess(itrial,:)=num_spikes;

    eligibility_vec_I=INP(:,tr:TT)*(V_vec(tr:TT,:)-V0I);
    eligibility_vec_E=INP(:,tr:TT)*((V_vec(tr:TT,:) - V0E).*((V_vec(tr:TT,:) - V0E) > 0));


    Filterd_Dw=0.99*Filterd_Dw+...
        0.01*[eligibility_vec_I(1:NI,:);C_E.*(eligibility_vec_E(NI+1:NN,:)-mean(eligibility_vec_E(NI+1:NN,:)))];
    W_vec(NI+1:NN)=AA*W_vec(NI+1:NN)*exp(alpha_1*(desired_S-Filterd_Spikes));
    W_vec=[W_vec(1:NI,:) + C_I.*Filterd_Dw(1:NI,:);W_vec(NI+1:NN,:) + Filterd_Dw(NI+1:NN,:)];


    W_vec=W_vec.*(W_vec > 0.);
    W_vec(NI+1:NN)=(W_vec(NI+1:NN).*(W_vec(NI+1:NN) <= w_max) + w_max*(W_vec(NI+1:NN) > w_max));
    Norm_W(itrial,1)=norm(W_vec);
    Cos_W(itrial,1)=sum(W_vec.*W_vec_0)/(norm(W_vec)*norm(W_vec_0));

    if rem(itrial,50)==0
       em=1;
       [I_ext_o,INP_o] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        %Neuron Model
        num_spikes_vec0=zeros(N_psn,TT);
        V_vec0=zeros(TT,N_psn);
        for it = 2:TT
            V_vec0(it,:) = ((1- dta)*V_vec0(it-1,:)) + I_ext_o(it,:)*dta;
            if V_vec0(it-1,:) >= threshy
                V_vec0(it,1)=resy;
                num_spikes_vec0(1,it)=1;
            end
        end   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ns_0(itrial,:) = sum(num_spikes_vec0(:,intimeP_vec(1,1)+1: intimeP_vec(1,1) + TP),2);
        ns_15(itrial,:) = sum(num_spikes_vec0(:,intimeP_vec(1,1)+1:intimeP_vec(1,1) + TP+150),2);
        Spikess_p(itrial,:)=sum(num_spikes_vec0(:,tr:TT),2);

    end

end


eval(['save Learning_BPB/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Learning_BPB/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['save Learning_BPB/Spikes/spike_p_',num2str(p1),'_',num2str(p2),  '.mat Spikess_p']);
eval(['save Learning_BPB/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['save Learning_BPB/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['save Learning_BPB/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['save Learning_BPB/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);
eval(['save Learning_BPB/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
eval(['save Learning_BPB/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
eval(['save Learning_BPB/Data_Cos/Cos_',num2str(p1),'_',num2str(p2),  '.mat Cos_W']);
eval(['save Learning_BPB/Data_Norm/Norm_',num2str(p1),'_',num2str(p2),  '.mat Norm_W']);

%step 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns_0=zeros(Ntrials_4,N_psn);
ns_15=zeros(Ntrials_4,N_psn);
Spikess=zeros(Ntrials_4,N_psn);
Norm_W=zeros(Ntrials_4,N_psn);
Cos_W=zeros(Ntrials_4,N_psn);
W_vec_0=W_vec;
em=1;
for itrial=1:Ntrials_4

    [I_ext,INP] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em);
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

    Norm_W(itrial,1)=norm(W_vec);
    Cos_W(itrial,1)=sum(W_vec.*W_vec_0)/(norm(W_vec)*norm(W_vec_0));

end

eval(['save Learning_BPBP/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Learning_BPBP/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
eval(['save Learning_BPBP/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
eval(['save Learning_BPBP/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['save Learning_BPBP/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['save Learning_BPBP/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['save Learning_BPBP/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['save Learning_BPBP/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);
eval(['save Learning_BPBP/Data_Cos/Cos_',num2str(p1),'_',num2str(p2),  '.mat Cos_W']);
eval(['save Learning_BPBP/Data_Norm/Norm_',num2str(p1),'_',num2str(p2),  '.mat Norm_W']);

%step 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Learning the second pattern

eval(['load Learning_BP/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),  '.mat PAT']);
eval(['load Learning_BP/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['load Learning_BP/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
eval(['load Learning_BP/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
eval(['load Learning_BP/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['load Learning_BP/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['load Learning_BP/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['load Learning_BP/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['load Learning_BP/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);
eval(['load Learning_BP/Data_Cos/Cos_',num2str(p1),'_',num2str(p2),  '.mat Cos_W']);
eval(['load Learning_BP/Data_Norm/Norm_',num2str(p1),'_',num2str(p2),  '.mat Norm_W']);


PAT2=PAT1(:,:,2);
PAT=PAT1(:,:,1);

ns_0=zeros(Ntrials_5,N_psn);
ns_15=zeros(Ntrials_5,N_psn);
Spikess=zeros(Ntrials_5,N_psn);

ns_0_0=zeros(Ntrials_5,N_psn);
ns_15_0=zeros(Ntrials_5,N_psn);
Spikess_0=zeros(Ntrials_5,N_psn);

Norm_W=zeros(Ntrials_5,N_psn);
Cos_W=zeros(Ntrials_5,N_psn);
W_vec_0=W_vec;
em=1;
for itrial=1:Ntrials_5

    [I_ext,INP] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT2,W_vec,inall,NI,NN,TT,traceI,traceE,em);
    
    
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

    Norm_W(itrial,1)=norm(W_vec);
    Cos_W(itrial,1)=sum(W_vec.*W_vec_0)/(norm(W_vec)*norm(W_vec_0));

    [I_ext_0,INP_0] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron Model
    num_spikes_vec_0=zeros(N_psn,TT);
    V_vec_0=zeros(TT,N_psn);
    for it = 2:TT
        V_vec_0(it,:) = ((1- dta)*V_vec_0(it-1,:)) + I_ext_0(it,:)*dta;
        if V_vec_0(it-1,:) >= threshy
            V_vec_0(it,1)=resy;
            num_spikes_vec_0(1,it)=1;
        end
    end   

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ns_0_0(itrial,:) = sum(num_spikes_vec_0(:,intimeP_vec(1,1)+1: intimeP_vec(1,1) + TP),2);
    ns_15_0(itrial,:) = sum(num_spikes_vec_0(:,intimeP_vec(1,1)+1:intimeP_vec(1,1) + TP+150),2);
    Spikess_0(itrial,:)=sum(num_spikes_vec_0(:,tr:TT),2);

    
end
eval(['save Learning_BPP/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),  '.mat PAT']);
eval(['save Learning_BPP/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Learning_BPP/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
eval(['save Learning_BPP/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
eval(['save Learning_BPP/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);

eval(['save Learning_BPP/Data_ns_0/ER2_',num2str(p1),'_',num2str(p2),  '.mat ns_0_0']);
eval(['save Learning_BPP/Data_ns_15/ER2_',num2str(p1),'_',num2str(p2),  '.mat ns_15_0']);
eval(['save Learning_BPP/Spikes/spike2_',num2str(p1),'_',num2str(p2),  '.mat Spikess_0']);


eval(['save Learning_BPP/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['save Learning_BPP/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['save Learning_BPP/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['save Learning_BPP/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);
eval(['save Learning_BPP/Data_Cos/Cos_',num2str(p1),'_',num2str(p2),  '.mat Cos_W']);
eval(['save Learning_BPP/Data_Norm/Norm_',num2str(p1),'_',num2str(p2),  '.mat Norm_W']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
