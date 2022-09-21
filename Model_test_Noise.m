function Model_test_Noise(p1,p2)


%p1: seed for random number and p1=1:1:\mu. There are \mu independent simulations
%p2: pattern length. p2=1 is for 50 ms;  p2=2 is for 100 ms; p2=3 is for 300 ms

rng(10+(p1-1)*1e7);

%There are 3 steps in this program.  

%%Step 1:  removing spikes

%%Step 2:  jittering spikes

%%Step 3:  adding spikes


%Data from this program is used in other programs to plot figures. 


N_psn=1;      %number of post-synaptic neuron


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

rateE = 5; %in [Hz] rate of excitatory neurons
rateI = 20 ;%in [Hz] rate of inhibitory neurons


tr=200;%transient time
dta=delt/taumem;

intimeP_vec = 5000; %ms Inserting time of the embedded pattern ()
TP_vec = [50 100 300];  %ms pattern leangth
TP=TP_vec(p2)/delt; 


inall=zeros(NN,2*TT-1);

%step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
em=1;

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

PAT=PAT(:,:,1);
NoiseRemove=zeros(3,1e4);
[y,x]=find(PAT==1);
tx=randperm(numel(y));

for jj=1:numel(y)
    PAT(y(tx(jj)),x(tx(jj)))=0;
    
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
    I_ext(:,1)=WW_vec(:,1)'*INP; 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %leaky integrate and fire model
    num_spikes_vec=zeros(N_psn,TT);
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
        if V_vec(it-1,:) >= threshy
            V_vec(it,1)=resy;
            num_spikes_vec(1,it)=1;
        end
    end   

    NoiseRemove(1,jj)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1: intimeP_vec(1,1) + TP),2);
    NoiseRemove(2,jj)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1:intimeP_vec(1,1) + TP+150),2);
    NoiseRemove(3,jj)=sum(num_spikes_vec(:,tr:TT),2);
   
end

eval(['save Noise_Robustness/NoiseRemove/NoiseRemove_',num2str(p1),'_',num2str(p2),  '.mat NoiseRemove']); 

%step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



NoiseJit=zeros(3,1e4);
fbb=1:1:600; % up to 60 ms jitter
PAT=PAT(:,:,1);
for jj=1:1:numel(fbb)
    fb=fbb(1,jj);

    BCV1=rand(NN,TT);
    PATP=zeros(NN,TT);
    PATP(NI+1:NN,:) = BCV1(NI+1:NN,:) < (rateE*delt/1000); 
    PATP(1:NI,:) = BCV1(1:NI,:) < (rateI*delt/1000);  
    if em>0
        PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
    end
    
    
    ts=3000;
    [y,x]=find(PATP(:,ts+1:end-ts)==1);
    nb=round(1+fb*randn(1,numel(x)));
    for i=1:numel(x)
        PATP(y(i),ts+x(i))=0;
        PATP(y(i),ts+x(i)+nb(i))=1;
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
    I_ext(:,1)=WW_vec(:,1)'*INP; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %leaky integrate and fire model
    num_spikes_vec=zeros(N_psn,TT);
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
        if V_vec(it-1,:) >= threshy
            V_vec(it,1)=resy;
            num_spikes_vec(1,it)=1;
        end
    end   

    
    NoiseJit(1,jj)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb:intimeP_vec(1,1) + TP+fb),2);
    NoiseJit(2,jj)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb:intimeP_vec(1,1) + TP+150+fb),2);
    NoiseJit(3,jj)=sum(num_spikes_vec(:,tr:TT),2);
    
end

eval(['save Noise_Robustness/NoiseJit/NoiseJit_',num2str(p1),'_',num2str(p2),  '.mat NoiseJit']); 


%step 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

NoiseAdd=zeros(3,1e4);
jno=0.1:0.1:150; %up to 150 adding spikes
PAT=PAT(:,:,1);
for jj=1:numel(jno)
    jadd=jno(jj);
    BCV1=rand(NN,TT);
    PATP=zeros(NN,TT);
    PATP(NI+1:NN,:) = BCV1(NI+1:NN,:) < (rateE*delt/1000); 
    PATP(1:NI,:) = BCV1(1:NI,:) < (rateI*delt/1000);  
    if em>0
        PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
    end
    PATP(NI+1:NN,:) = PATP(NI+1:NN,:)+ (rand(400,TT) < (jadd*delt/1000));     
    PATP(1:NI,:) = PATP(1:NI,:)+(rand(NI,TT) < (4*jadd*delt/1000));

    WW_vec = W_vec;
    WW_vec(1:NI,:) = -W_vec(1:NI,:);
    for i = 1: NI
        inall(i,:) = conv(PATP(i,:), traceI);        
    end 
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end    
    INP = inall(:,1:TT);   
    I_ext(:,1)=WW_vec(:,1)'*INP; 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %leaky integrate and fire model
    num_spikes_vec=zeros(N_psn,TT);
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
        if V_vec(it-1,:) >= threshy
            V_vec(it,1)=resy;
            num_spikes_vec(1,it)=1;
        end
    end   

    NoiseAdd(1,jj)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1: intimeP_vec(1,1) + TP),2);
    NoiseAdd(2,jj)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1:intimeP_vec(1,1) + TP+150),2);
    NoiseAdd(3,jj)=sum(num_spikes_vec(:,tr:TT),2);
    
end

eval(['save Noise_Robustness/NoiseAdd/NoiseAdd_',num2str(p1),'_',num2str(p2),  '.mat NoiseAdd']); 

end
