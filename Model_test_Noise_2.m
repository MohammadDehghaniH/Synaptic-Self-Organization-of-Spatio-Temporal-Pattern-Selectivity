function Model_test_Noise_2(p1,p2)


%p1: seed for random number and p1=1:1:\mu. There are \mu independent simulations
%p2: pattern length. p2=1 is for 50 ms; 

rng(10+(p1-1)*1e7);

%There are 2 steps in this program.  

%%Step 1: jittering spikes

%%Step 2: gaussian distribution


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
rateE = 5; %in [Hz] rate of excitatory neuron
rateI = 20 ;%in [Hz] rate of inhibitory neuron

W_vec = 0.01+1e-3*randn(NN,N_psn); % initial weight vector
W_vec(W_vec<0)=0;

tr=200;%transient time
dta=delt/taumem;

intimeP_vec = 5000; %ms Inserting time of the embedded pattern ()
TP_vec = [50 100 300];  %ms pattern leangth
TP=TP_vec(p2)/delt; 


%step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
em=1;

eval(['load Noise_Robustness/Learn_Jit_2/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),  '.mat PAT']);
eval(['load Noise_Robustness/Learn_Jit_2/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);


NoiseJit=zeros(3,1e4);
fbb=1:1:600; % up to 60 ms jitter
PAT=PAT(:,:,1);
inall=zeros(NN,2*TT-1);
for j=1:1:numel(fbb)
    fb=fbb(1,j);

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

    
    NoiseJit(1,j)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb:intimeP_vec(1,1) + TP+fb),2);
    NoiseJit(2,j)=sum(num_spikes_vec(:,intimeP_vec(1,1)+1-fb:intimeP_vec(1,1) + TP+150+fb),2);
    NoiseJit(3,j)=sum(num_spikes_vec(:,tr:TT),2);
    
end

eval(['save Noise_Robustness/NoiseJit_2/L_NoiseJit_',num2str(p1),'_',num2str(p2),  '.mat NoiseJit']); 


%step 2: Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['load Noise_Robustness/Learn_Gau_2/Data_Pattern/EP_',num2str(p1),'_',num2str(p2),  '.mat PAT']);
eval(['load Noise_Robustness/Learn_Gau_2/Weights/W_',num2str(p1),'_',num2str(p2),  '.mat W_vec']);


NoiseGaus=zeros(3,1e4);
sigha=0.1:0.1:60; %standard deviation
inall0=zeros(NN,TT);
PAT=PAT(:,:,1);
for j=1:1:numel(sigha)
    fb=sigha(1,j);

    traceG=normpdf(times-500,0,fb); %Gaussian kernel
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
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NoiseGaus(1,j) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1: intimeP_vec(1,1) + TP),2);
    NoiseGaus(2,j) = sum(num_spikes_vec(:,intimeP_vec(1,1)+1:intimeP_vec(1,1) + TP+150),2);
    NoiseGaus(3,j)=num_spikes;    
end

eval(['save Noise_Robustness/NoiseGaus_2/L_NoiseGaus_',num2str(p1),'_',num2str(p2),  '.mat NoiseGaus']); 

end
