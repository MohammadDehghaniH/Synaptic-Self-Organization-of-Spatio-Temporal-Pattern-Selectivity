function Model_Learning_NET_1(p1,p2)
%p1: seed for random number and p1=1:1:\mu. There are \mu independent simulations
%p2: N_psn: number of post-synaptic neurons. (p2=1: N_psn=4);(p2=1: N_psn=5)
%(p2=1: N_psn=6);(p2=1: N_psn=7)
%w =ab w: W_vec a: W_veca b: W_vecb

SDF1=10+(p1-1)*1e6;
rng(SDF1);


Ntrials0=2000;
%%Step 1: There is no embedded pattern in afferents for 2000 learning cycles.


Ntrials1 =5000; 
%%Step 2: There are 4 embedded patterns in afferents for 5000 learning
%%cycles. There is pre-synaptic competition.

Ntrials2 =5000; 
%%Step 3: There is one embedded pattern in afferents for 5000 learning
%%cycles. There is no pre-synaptic competition.


N_psn_vec=[4 5 6 7];
N_psn=N_psn_vec(p2); %number of post-synaptic neurons

em=4; % number of embedded patterns
w_max = 1; % maximum amount of synapsis.
NN = 500;  %number of input neurons
percI = 0.2; %percentage of inhibition
NI = fix(percI*NN); %number of inhibitory neurons



delt = 0.1; %integration steps

TP=50/delt; %pattern length (ms/delt= steps!)

Filterd_Spikes=zeros(1,N_psn);%long-time firing rate

desired_S=2; %desired number of spikes

%initializations: 
del_b_E  =  zeros(NN-NI,N_psn); %\del b_E
del_a_E  =  zeros(NN-NI,N_psn); %\del a_E
del_b_I  = zeros(NI,N_psn);     %\del b_I and \del a_I

alpha_1=1e-2; %scaling factor 
C_E=0.9e-4; %excitatory learning rate
bet_a=0.9e-4;%scaling coefficient
AA=1-bet_a;
C_I = 1e-4;%inhibitory learning rate

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

taumem = 15;   %membrane time constant
tautraceE =3;  %Rise time of excitatory currents
tautraceI = 5; %Decay time of inhibitory currents
tauriseE =0.5; %Decay time of excitatory currents
tauriseI = 1;  %Rise time of inhibitory currents

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

tr=200;%transient time
dta=delt/taumem;

%for putting the pattern into the background:
%starting time (in ms, like everywhere!)
intimeP_vec=(1./delt)*[250 450 650 850];

% initial weight vectors
W_veca = 0.1+1e-2*randn(NN,N_psn);
W_veca(W_veca<0)=0;

W_vecb = 0.1+1e-2*randn(NN,N_psn);
W_vecb(W_vecb<0)=0;

%the patterns
PAT1 = rand(NN, TP,em);
PAT1(1:NI,:,:) = PAT1(1:NI,:,:) < (rateI*delt/1000);
PAT1(NI+1:NN,:,:) = PAT1(NI+1:NN,:,:) < (rateE*delt/1000);
PAT=PAT1;

%step1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emm=0;
for itrial=1:Ntrials0
        
        [I_ext,INP]=Model_Input_NET(rateI,delt,rateE,intimeP_vec,W_veca.*W_vecb,traceI,traceE,NN,TT,N_psn,NI,em,PAT,TP,emm);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        %Neuron Model
        num_spikes=zeros(N_psn,TT);
        V_vec=zeros(TT,N_psn);
        for it = 2:TT
            V_vec(it,:) = (1- dta)*(V_vec(it-1,:)) + I_ext(it-1,:)*dta;
            f_p = find(V_vec(it-1,:) >= threshy);
            if ~isempty(f_p)
                V_vec(it,f_p) =resy; 
                num_spikes(f_p,it) = num_spikes(f_p) + 1; 
            end
        end
        Filterd_Spikes = 0.9*Filterd_Spikes + 0.1*sum(num_spikes(:,tr:TT),2)';   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        del_b_I = 0.99*del_b_I + 0.01*(abs(INP(1:NI,tr:TT))*(V_vec(tr:TT,:)));  

        Filterd_Dw_E= C_E.*((INP(NI+1:NN,tr:TT)*((V_vec(tr:TT,:) ).*((V_vec(tr:TT,:) ) > 0)))...
            -mean((INP(NI+1:NN,tr:TT)*((V_vec(tr:TT,:) ).*((V_vec(tr:TT,:) ) > 0)))));

        del_b_E = 0.99*del_b_E+ 0.01*Filterd_Dw_E;

        W_vecb=[W_vecb(1:NI,:) + C_I*del_b_I;...
            (AA*W_vecb(NI+1:NN,:).*exp(alpha_1*(desired_S-Filterd_Spikes)))+del_b_E];
        SDF=sum((Filterd_Dw_E>0),2);
        del_a_E = 0.99*del_a_E+ 0.01*((Filterd_Dw_E.*(Filterd_Dw_E>0))-((SDF>1).*sum((Filterd_Dw_E.*(Filterd_Dw_E>0)),2)./(SDF+1e-16)));

        W_veca=[(W_veca(1:NI,:) + C_I*del_b_I);...
            (AA*W_veca(NI+1:NN,:).*exp(alpha_1*(desired_S-Filterd_Spikes))) + (del_a_E)];

        W_veca = W_veca.*(W_veca > 0.);
        W_vecb = W_vecb.*(W_vecb > 0.);

        W_veca(NI+1:NN,:)=W_veca(NI+1:NN,:).*(W_veca(NI+1:NN,:) <= w_max) + w_max*(W_veca(NI+1:NN,:) > w_max);
        W_vecb(NI+1:NN,:)=W_vecb(NI+1:NN,:).*(W_vecb(NI+1:NN,:) <= w_max) + w_max*(W_vecb(NI+1:NN,:) > w_max);

end
    
eval(['save ',num2str(N_psn,'%.i') '/N_Wa_',num2str(p1),'_',num2str(N_psn),  '.mat W_veca']);
eval(['save ',num2str(N_psn,'%.i') '/N_Wb_',num2str(p1),'_',num2str(N_psn),  '.mat W_vecb']);
eval(['save ',num2str(N_psn,'%.i') '/N_EP_',num2str(p1),'_',num2str(N_psn),  '.mat PAT']);
eval(['save ',num2str(N_psn,'%.i') '/N_Filterd_Spikes_',num2str(p1),'_',num2str(N_psn),  '.mat Filterd_Spikes']);
eval(['save ',num2str(N_psn,'%.i') '/N_del_a_E_',num2str(p1),'_',num2str(N_psn),  '.mat del_a_E']);
eval(['save ',num2str(N_psn,'%.i') '/N_del_b_I_',num2str(p1),'_',num2str(N_psn),  '.mat del_b_I']);
eval(['save ',num2str(N_psn,'%.i') '/N_del_b_E_',num2str(p1),'_',num2str(N_psn),  '.mat del_b_E']);


%step2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%with competition

Omega_R=zeros(1,Ntrials1);
emm=1;
ns_15= zeros(Ntrials1,N_psn);
Spikess_p = zeros(Ntrials1,N_psn);
for itrial=1:Ntrials1

        [I_ext,INP]=Model_Input_NET(rateI,delt,rateE,intimeP_vec,W_veca.*W_vecb,traceI,traceE,NN,TT,N_psn,NI,em,PAT,TP,emm);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        %Neuron Model
        num_spikes=zeros(N_psn,TT);
        V_vec=zeros(TT,N_psn);
        for it = 2:TT
            V_vec(it,:) = (1- dta)*(V_vec(it-1,:)) + I_ext(it-1,:)*dta;
            f_p = find(V_vec(it-1,:) >= threshy);
            if ~isempty(f_p)
                V_vec(it,f_p) =resy; 
                num_spikes(f_p,it) = num_spikes(f_p) + 1; 
            end
        end
        Filterd_Spikes = 0.9*Filterd_Spikes + 0.1*sum(num_spikes(:,tr:TT),2)';   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ns_15(itrial,:) =sum(num_spikes(:,intimeP_vec(1)+1: intimeP_vec(1)+ TP+150),2)+...
            sum(num_spikes(:,intimeP_vec(2)+1: intimeP_vec(2)+ TP+150),2)+...
            sum(num_spikes(:,intimeP_vec(3)+1: intimeP_vec(3)+ TP+150),2)+...
            sum(num_spikes(:,intimeP_vec(4)+1: intimeP_vec(4)+ TP+150),2);

        Spikess_p(itrial,:)=sum(num_spikes(:,tr:TT),2)';
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Neuron model
        V_vec2=V_vec>=1;
        V_vec1=zeros(em,N_psn);
        for io=1:em
            V_vec1(io,:)=sum(V_vec2(intimeP_vec(io)+1:intimeP_vec(io)+TP+150,:),1);
        end
        Omega_R(1,itrial)=rank(double(V_vec1>0))/(1e-16+em);  
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        del_b_I = 0.99*del_b_I + 0.01*(abs(INP(1:NI,tr:TT))*(V_vec(tr:TT,:)));  


        Filterd_Dw_E= C_E.*((INP(NI+1:NN,tr:TT)*((V_vec(tr:TT,:) ).*((V_vec(tr:TT,:) ) > 0)))...
            -mean((INP(NI+1:NN,tr:TT)*((V_vec(tr:TT,:) ).*((V_vec(tr:TT,:) ) > 0)))));


        del_b_E = 0.99*del_b_E+ 0.01*Filterd_Dw_E;

        W_vecb=[W_vecb(1:NI,:) + C_I*del_b_I;...
            (AA*W_vecb(NI+1:NN,:).*exp(alpha_1*(desired_S-Filterd_Spikes)))+del_b_E];

        SDF=sum((Filterd_Dw_E>0),2);

        del_a_E = 0.99*del_a_E+ 0.01*((Filterd_Dw_E.*(Filterd_Dw_E>0))-((SDF>1).*sum((Filterd_Dw_E.*(Filterd_Dw_E>0)),2)./(SDF+1e-16)));


        W_veca=[(W_veca(1:NI,:) + C_I*del_b_I);...
            (AA*W_veca(NI+1:NN,:).*exp(alpha_1*(desired_S-Filterd_Spikes))) + (del_a_E)];

        W_veca = W_veca.*(W_veca > 0.);
        W_vecb = W_vecb.*(W_vecb > 0.);

        W_veca(NI+1:NN,:)=W_veca(NI+1:NN,:).*(W_veca(NI+1:NN,:) <= w_max) + w_max*(W_veca(NI+1:NN,:) > w_max);
        W_vecb(NI+1:NN,:)=W_vecb(NI+1:NN,:).*(W_vecb(NI+1:NN,:) <= w_max) + w_max*(W_vecb(NI+1:NN,:) > w_max);
end



eval(['save ',num2str(N_psn,'%.i') '/With_Wa_',num2str(p1),'_',num2str(N_psn),  '.mat W_veca']);
eval(['save ',num2str(N_psn,'%.i') '/With_Wb_',num2str(p1),'_',num2str(N_psn),  '.mat W_vecb']);
eval(['save ',num2str(N_psn,'%.i') '/With_EP_',num2str(p1),'_',num2str(N_psn),  '.mat PAT']);
eval(['save ',num2str(N_psn,'%.i') '/With_Omega_',num2str(p1),'_',num2str(N_psn),  '.mat Omega_R']);
eval(['save ',num2str(N_psn,'%.i') '/With_ns15_',num2str(p1),'_',num2str(N_psn),  '.mat ns_15']);
eval(['save ',num2str(N_psn,'%.i') '/With_Spikess_',num2str(p1),'_',num2str(N_psn),  '.mat Spikess_p']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 3
%%%%without competition

eval(['load ',num2str(N_psn,'%.i') '/N_Wa_',num2str(p1),'_',num2str(N_psn),  '.mat W_veca']);
eval(['load ',num2str(N_psn,'%.i') '/N_Wb_',num2str(p1),'_',num2str(N_psn),  '.mat W_vecb']);
eval(['load ',num2str(N_psn,'%.i') '/N_EP_',num2str(p1),'_',num2str(N_psn),  '.mat PAT']);
eval(['load ',num2str(N_psn,'%.i') '/N_Filterd_Spikes_',num2str(p1),'_',num2str(N_psn),  '.mat Filterd_Spikes']);
eval(['load ',num2str(N_psn,'%.i') '/N_del_a_E_',num2str(p1),'_',num2str(N_psn),  '.mat del_a_E']);
eval(['load ',num2str(N_psn,'%.i') '/N_del_b_I_',num2str(p1),'_',num2str(N_psn),  '.mat del_b_I']);
eval(['load ',num2str(N_psn,'%.i') '/N_del_b_E_',num2str(p1),'_',num2str(N_psn),  '.mat del_b_E']);



Omega_R=zeros(1,Ntrials1);
emm=1;
for itrial=1:Ntrials2
    
    [I_ext,INP]=Model_Input_NET(rateI,delt,rateE,intimeP_vec,W_veca.*W_vecb,traceI,traceE,NN,TT,N_psn,NI,em,PAT,TP,emm);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Neuron Model
    num_spikes=zeros(N_psn,TT);
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = (1- dta)*(V_vec(it-1,:)) + I_ext(it-1,:)*dta;
        f_p = find(V_vec(it-1,:) >= threshy);
        if ~isempty(f_p)
            V_vec(it,f_p) =resy; 
            num_spikes(f_p,it) = num_spikes(f_p) + 1; 
        end
    end
    Filterd_Spikes = 0.9*Filterd_Spikes + 0.1*sum(num_spikes(:,tr:TT),2)';   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_vec2=V_vec>=1;
    V_vec1=zeros(em,N_psn);
    for io=1:em
        V_vec1(io,:)=sum(V_vec2(intimeP_vec(io)+1:intimeP_vec(io)+TP+150,:),1);
    end
    Omega_R(1,itrial)=rank(double(V_vec1>0))/(1e-16+em);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

       
        
    eligibility_vec_I=INP(:,tr:TT)*(V_vec(tr:TT,:));  
    del_b_I = 0.99*del_b_I + 0.01*eligibility_vec_I(1:NI,:);  
 
    eligibility_vec_E=INP(:,tr:TT)*((V_vec(tr:TT,:) ).*((V_vec(tr:TT,:)) > 0));  
    W_vecb(NI:NN,:)=AA*W_vecb(NI:NN,:).*exp(alpha_1*(desired_S-Filterd_Spikes));

    
    Filterd_Dw_E= C_E.*(eligibility_vec_E(NI+1:NN,:)-mean(eligibility_vec_E(NI+1:NN,:)));

    
    del_b_E = 0.99*del_b_E+ 0.01*Filterd_Dw_E;
    
    W_vecb=[W_vecb(1:NI,:) + C_I*del_b_I; W_vecb(NI+1:NN,:)+del_b_E];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Filterd_Dw_EX=Filterd_Dw_E.*(Filterd_Dw_E>0);
    Filterd_Dw_Ea= Filterd_Dw_EX;
    
    del_a_E = 0.99*del_a_E+ 0.01*Filterd_Dw_Ea;

    W_veca(NI:NN,:)=AA*W_veca(NI:NN,:).*exp(alpha_1*(desired_S-Filterd_Spikes));

    W_veca=[(W_veca(1:NI,:) + C_I*del_b_I);....
        W_veca(NI+1:NN,:) + (del_a_E)];

    W_veca = W_veca.*(W_veca > 0.);
    W_vecb = W_vecb.*(W_vecb > 0.);

    W_veca=W_veca.*(W_veca <= w_max) + w_max*(W_veca > w_max);
    W_vecb=W_vecb.*(W_vecb <= w_max) + w_max*(W_vecb > w_max);


end

eval(['save ',num2str(N_psn,'%.i') ,'/Without_Wa_',num2str(p1),'_',num2str(N_psn),  '.mat W_veca']);
eval(['save ',num2str(N_psn,'%.i') ,'/Without_Wb_',num2str(p1),'_',num2str(N_psn),  '.mat W_vecb']);
eval(['save ',num2str(N_psn,'%.i') ,'/Without_del_b_I_',num2str(p1),'_',num2str(N_psn),  '.mat del_b_I']);

eval(['save ',num2str(N_psn,'%.i') ,'/Without_del_b_E_',num2str(p1),'_',num2str(N_psn),  '.mat del_b_E']);
eval(['save ',num2str(N_psn,'%.i') ,'/Without_EP_',num2str(p1),'_',num2str(N_psn),  '.mat PAT']);
eval(['save ',num2str(N_psn,'%.i') ,'/Without_Filterd_Spikes_',num2str(p1),'_',num2str(N_psn),  '.mat Filterd_Spikes']);

eval(['save ',num2str(N_psn,'%.i') '/Without_Omega_',num2str(p1),'_',num2str(N_psn),  '.mat Omega_R']);

end


