function Model_Learning_noise(p1,p2)


%p1: seed for random number and p1=1:1:\mu. There are \mu independent simulations
%p2: pattern length. p2=1 is for 50 ms;  p2=2 is for 100 ms; p2=3 is for 300 ms

rng(10+(p1-1)*1e6);

%There are 3 steps in this program.  

%%Step 1: There is no embedded pattern in afferents for 2000 learning cycles.
%%Step 2: learning with jitter.
%%Step 3: learning with Gaussian.

%Data from this program is used in other programs to plot figures. 

Ntrials_1=2000;
%Ntrials_1: number of learning cycles to learn the background.
%There is no embedded pattern in afferents for 2000 learning cycles. (step 1)


Ntrials_2=5000; 
%Ntrials_2: number of learning cycles to learn the background containing
%the embedded pattern.
%(Note initial conditions in this step are from step1)


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
T = 1000; % background leangth
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

Filterd_Dw=zeros(NN,N_psn); %\tilde{\varepsilon} : eligibility after applying hetero-synaptic plasticity


inall=zeros(NN,2*TT-1);
Filterd_Spikes=zeros(1,N_psn);%long-time firing rate
tr=200;%transient time
dta=delt/taumem;

intimeP_vec = 5000;     %ms Inserting time of the embedded pattern
TP_vec = [50 100 300];  %ms pattern leangth
TP=TP_vec(1)/delt; 


V0I=0; %modification threshold for inhibitory neurons
V0E=0; %modification threshold for excitatory neurons


% generating randomly embedded pattern
PAT1 = rand(NN, TP,2);
PAT1(1:NI,:,:) = 1.*(PAT1(1:NI,:,:) < rateI*delt/1000);
PAT1(NI+1:NN,:,:) = 1.*(PAT1(NI+1:NN,:,:) < rateE*delt/1000);
PAT=PAT1(:,:,1);

W_vec = 0.01+1e-3*randn(NN,N_psn); % initial weight vector
W_vec(W_vec<0)=0;

%step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spikess=zeros(Ntrials_1,N_psn);
em=0;
for itrial=1:Ntrials_1

    [I_ext,INP] = Model_Input(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em);
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron model
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eval(['save Noise_Robustness/Learning_B/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Noise_Robustness/Learning_B/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['save Noise_Robustness/Learning_B/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['save Noise_Robustness/Learning_B/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['save Noise_Robustness/Learning_B/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);


%step 2(jitter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ns_0=zeros(Ntrials_2,N_psn);
ns_15=zeros(Ntrials_2,N_psn);
Spikess=zeros(Ntrials_2,N_psn);
ns_0_1=zeros(Ntrials_2,N_psn);
ns_15_1=zeros(Ntrials_2,N_psn);
Spikess_1=zeros(Ntrials_2,N_psn);

em=1;
fb=200;
for itrial=1:Ntrials_2
    
    %%%%original pattern
    BCV1=rand(NN,TT);
    PATP=zeros(NN,TT);
    PATP(NI+1:NN,:) = BCV1(NI+1:NN,:) < (rateE*delt/1000); 
    PATP(1:NI,:) = BCV1(1:NI,:) < (rateI*delt/1000);  
    if em>0
        PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
    end
    WW_vec = W_vec;
    WW_vec(1:NI,:) = -W_vec(1:NI,:);
    for i = 1: NI
        inall(i,:) = conv(PATP(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    I_ext1(:,1)=WW_vec(:,1)'*INP; 

    %%%%jittering pattern
    ts=3000;
    [y,x]=find(PATP(:,ts+1:end-ts)==1);
    
    nb=round(0+200*randn(1,numel(x)));
    for i=1:numel(x)
        PATP(y(i),ts+x(i))=0;
        PATP(y(i),ts+x(i)+nb(i))=1;
    end    

    for i = 1: NI
        inall(i,:) = conv(PATP(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    I_ext(:,1)=WW_vec(:,1)'*INP; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron model
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron model
    num_spikes_vec=zeros(N_psn,TT);
    V_vec1=zeros(TT,N_psn);
    for it = 2:TT
        V_vec1(it,:) = ((1- dta)*V_vec1(it-1,:)) + I_ext1(it,:)*dta;
        if V_vec1(it-1,:) >= threshy
            V_vec1(it,1)=resy;
            num_spikes_vec(1,it)=1;
        end
    end   
    num_spikes1=sum(num_spikes_vec(:,tr:TT),2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ns_0_1(itrial,:) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb: intimeP_vec(1,1) + TP+fb),2);
    ns_15_1(itrial,:) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb:intimeP_vec(1,1) + TP+150+fb),2);
    Spikess_1(itrial,:)=num_spikes1;

    
    
    eligibility_vec_I=INP(:,tr:TT)*(V_vec(tr:TT,:)-V0I);
    eligibility_vec_E=INP(:,tr:TT)*((V_vec(tr:TT,:) - V0E).*((V_vec(tr:TT,:) - V0E) > 0));

    Filterd_Dw=0.99*Filterd_Dw+...
        0.01*[eligibility_vec_I(1:NI,:);C_E.*(eligibility_vec_E(NI+1:NN,:)-mean(eligibility_vec_E(NI+1:NN,:)))];
    W_vec(NI+1:NN)=AA*W_vec(NI+1:NN)*exp(alpha_1*(desired_S-Filterd_Spikes));
    W_vec=[W_vec(1:NI,:) + C_I.*Filterd_Dw(1:NI,:);W_vec(NI+1:NN,:) + Filterd_Dw(NI+1:NN,:)];

    W_vec=W_vec.*(W_vec > 0.);
    W_vec(NI+1:NN)=(W_vec(NI+1:NN).*(W_vec(NI+1:NN) <= w_max) + w_max*(W_vec(NI+1:NN) > w_max));


end
eval(['save Noise_Robustness/Learn_Jit_2/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),  '.mat PAT']);
eval(['save Noise_Robustness/Learn_Jit_2/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Noise_Robustness/Learn_Jit_2/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
eval(['save Noise_Robustness/Learn_Jit_2/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
eval(['save Noise_Robustness/Learn_Jit_2/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['save Noise_Robustness/Learn_Jit_2/Data_ns_0/ER_0_',num2str(p1),'_',num2str(p2),  '.mat ns_0_1']);
eval(['save Noise_Robustness/Learn_Jit_2/Data_ns_15/ER_1_',num2str(p1),'_',num2str(p2),  '.mat ns_15_1']);
eval(['save Noise_Robustness/Learn_Jit_2/Spikes/spike_1_',num2str(p1),'_',num2str(p2),  '.mat Spikess_1']);



%step 3 (Gaussian)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eval(['load Noise_Robustness/Learning_B/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['load Noise_Robustness/Learning_B/Filterd_Spikes/Filterd_Spikes_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Spikes']);
eval(['load Noise_Robustness/Learning_B/Eligibility_vec/eligibility_vec_I_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_I']);
eval(['load Noise_Robustness/Learning_B/Eligibility_vec/eligibility_vec_E_',num2str(p1),'_',num2str(p2),  '.mat eligibility_vec_E']);
eval(['load Noise_Robustness/Learning_B/Filterd_Dw/Filterd_Dw_',num2str(p1),'_',num2str(p2),  '.mat Filterd_Dw']);

ns_0=zeros(Ntrials_2,N_psn);
ns_15=zeros(Ntrials_2,N_psn);
Spikess=zeros(Ntrials_2,N_psn);
ns_0_1=zeros(Ntrials_2,N_psn);
ns_15_1=zeros(Ntrials_2,N_psn);
Spikess_1=zeros(Ntrials_2,N_psn);

em=1;

sigha=20; %standard deviation
traceG=normpdf(times-500,0,sigha); %Gaussian kernel
inall0=zeros(NN,TT);

for itrial=1:Ntrials_2
    
    %%%%original pattern
    BCV1=rand(NN,TT);
    BCV2=rand(NN,TT);
    PATP=zeros(NN,TT);
    PATP(NI+1:NN,:) = BCV1(NI+1:NN,:) < (rateE*delt/1000); 
    PATP(1:NI,:) = BCV1(1:NI,:) < (rateI*delt/1000);  
    if em>0
        PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
    end
    WW_vec = W_vec;
    WW_vec(1:NI,:) = -W_vec(1:NI,:);
    for i = 1: NI
        inall(i,:) = conv(PATP(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    I_ext1(:,1)=WW_vec(:,1)'*INP; 

    %%%%gaussian
    
    for i=1:NN
       inall0(i,:) =conv(PATP(i,:), traceG,'same');
    end
    PATP2=BCV2<inall0*delt;
    for i = 1: NI
        inall(i,:) = conv(PATP2(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP2(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    I_ext(:,1)=WW_vec(:,1)'*INP;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron model
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
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron model
    num_spikes_vec=zeros(N_psn,TT);
    V_vec1=zeros(TT,N_psn);
    for it = 2:TT
        V_vec1(it,:) = ((1- dta)*V_vec1(it-1,:)) + I_ext1(it,:)*dta;
        if V_vec1(it-1,:) >= threshy
            V_vec1(it,1)=resy;
            num_spikes_vec(1,it)=1;
        end
    end   
    num_spikes1=sum(num_spikes_vec(:,tr:TT),2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ns_0_1(itrial,:) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb: intimeP_vec(1,1) + TP+fb),2);
    ns_15_1(itrial,:) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb:intimeP_vec(1,1) + TP+150+fb),2);
    Spikess_1(itrial,:)=num_spikes1;
    
    eligibility_vec_I=INP(:,tr:TT)*(V_vec(tr:TT,:)-V0I);
    eligibility_vec_E=INP(:,tr:TT)*((V_vec(tr:TT,:) - V0E).*((V_vec(tr:TT,:) - V0E) > 0));

    Filterd_Dw=0.99*Filterd_Dw+...
        0.01*[eligibility_vec_I(1:NI,:);C_E.*(eligibility_vec_E(NI+1:NN,:)-mean(eligibility_vec_E(NI+1:NN,:)))];
    W_vec(NI+1:NN)=AA*W_vec(NI+1:NN)*exp(alpha_1*(desired_S-Filterd_Spikes));
    W_vec=[W_vec(1:NI,:) + C_I.*Filterd_Dw(1:NI,:);W_vec(NI+1:NN,:) + Filterd_Dw(NI+1:NN,:)];

    W_vec=W_vec.*(W_vec > 0.);
    W_vec(NI+1:NN)=(W_vec(NI+1:NN).*(W_vec(NI+1:NN) <= w_max) + w_max*(W_vec(NI+1:NN) > w_max));

end

eval(['save Noise_Robustness/Learn_Gau_2/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),  '.mat PAT']);
eval(['save Noise_Robustness/Learn_Gau_2/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);
eval(['save Noise_Robustness/Learn_Gau_2/Data_ns_0/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_0']);
eval(['save Noise_Robustness/Learn_Gau_2/Data_ns_15/ER_',num2str(p1),'_',num2str(p2),  '.mat ns_15']);
eval(['save Noise_Robustness/Learn_Gau_2/Spikes/spike_',num2str(p1),'_',num2str(p2),  '.mat Spikess']);
eval(['save Noise_Robustness/Learn_Gau_2/Data_ns_0/ER_0_',num2str(p1),'_',num2str(p2),  '.mat ns_0_1']);
eval(['save Noise_Robustness/Learn_Gau_2/Data_ns_15/ER_1_',num2str(p1),'_',num2str(p2),  '.mat ns_15_1']);
eval(['save Noise_Robustness/Learn_Gau_2/Spikes/spike_1_',num2str(p1),'_',num2str(p2),  '.mat Spikess_1']);

end
