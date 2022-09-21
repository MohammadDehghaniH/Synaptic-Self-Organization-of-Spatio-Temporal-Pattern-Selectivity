clc; clear; close all

Fig2b=figure();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yline(0,'-k','Linewidth',2)
hold on
xline(0,'-k','Linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters (see Model_Learning.m )
N_mu=500; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

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

intimeP_vec = 5000; %[ms] time of the embedded pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E00V=zeros(N_mu,TT); %excitatory input
I00V=zeros(N_mu,TT); %inhibitory input
S00V=zeros(N_mu,TT); %spike time


TP=50/delt; %pattern length (ms/delt= steps)
for p1=1:N_mu

   try
        eval(['load Learning_BP/Data_Pattern/EP_',num2str(p1),'_',num2str(1),  '.mat PAT']);
        eval(['load Learning_BP/Weights/W_',num2str(p1),'_',num2str(1),  '.mat W_vec']);
    catch       
    end
    
    PATR(NI+1:NN,:) = rand(400,TT) < (rateE*delt/1000); 
    inall=zeros(NN,2*TT-1);
    PATR(1:NI,:) = (rand(NI,TT) < rateI*delt/1000);
   
    PATP = PATR;
    PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT(:,:,1); 
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

     I_ext=(INPI+INPE)';
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S00V(p1,:)=num_spikes_vec;
    E00V(p1,:) = INPE;    
    I00V(p1,:) = INPI;

end
Lem=700;
dfff=(-Lem+1:1:Lem)/10;
DF=sum(S00V(:));
ET=zeros(50*DF,numel(dfff));
IT=zeros(50*DF,numel(dfff));
ghj=0;

for j=1:N_mu
    df=find(S00V(j,:)>0);
    ghj=ghj+numel(df);
    for kj=1:numel(df)
        if df(kj)<20000-Lem
            if df(kj)>Lem
            i=(j-1)*20+kj;
            df(kj);
            ET(i,:)=E00V(j,df(kj)-Lem:df(kj)+Lem-1);
            IT(i,:)=I00V(j,df(kj)-Lem:df(kj)+Lem-1);
            end
        end
    end
end

FGG=sum(ET,2);
gk=numel(find(FGG>0));

plot(dfff,sum(ET,1)/gk,'color',[0 0 0.5],'DisplayName','L_{em} = 50 ms ','Linewidth',2.5);
hold on
I1= plot(dfff,sum(IT,1)/gk,'color',[0 0 0.5],'DisplayName','L_{em} = 50 ms ','Linewidth',2.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E00V=zeros(N_mu,TT);
I00V=zeros(N_mu,TT);
S00V=zeros(N_mu,TT);

TP=100/delt; %pattern length (ms/delt= steps)

for p1=1:N_mu

   try
        eval(['load Learning_BP/Data_Pattern/EP_',num2str(p1),'_',num2str(2),  '.mat PAT']);
        eval(['load Learning_BP/Weights/W_',num2str(p1),'_',num2str(2),  '.mat W_vec']);
    catch       
    end
    
    PATR(NI+1:NN,:) = rand(400,TT) < (rateE*delt/1000); 
    inall=zeros(NN,2*TT-1);
    PATR(1:NI,:) = (rand(NI,TT) < rateI*delt/1000);
   
    PATP = PATR;
    PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT(:,:,1); 
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

     I_ext=(INPI+INPE)';
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S00V(p1,:)=num_spikes_vec;
    E00V(p1,:) = INPE;    
    I00V(p1,:) = INPI;

end
dfff=(-Lem+1:1:Lem)/10;
DF=sum(S00V(:));
ET=zeros(50*DF,numel(dfff));
IT=zeros(50*DF,numel(dfff));
ghj=0;

for j=1:N_mu
    df=find(S00V(j,:)>0);
    ghj=ghj+numel(df);
    for kj=1:numel(df)
        if df(kj)<20000-Lem
            if df(kj)>Lem
            i=(j-1)*20+kj;
            df(kj);
            ET(i,:)=E00V(j,df(kj)-Lem:df(kj)+Lem-1);
            IT(i,:)=I00V(j,df(kj)-Lem:df(kj)+Lem-1);
            end
        end
    end
end

FGG=sum(ET,2);
gk=numel(find(FGG>0));

plot(dfff,sum(ET,1)/gk,'color',[0 0.5 0],'DisplayName','L_{em} = 100 ms ','Linewidth',2.5);
hold on
I2= plot(dfff,sum(IT,1)/gk,'color',[0 0.5 0],'DisplayName','L_{em} = 100 ms ','Linewidth',2.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E00V=zeros(N_mu,TT);
I00V=zeros(N_mu,TT);
S00V=zeros(N_mu,TT);
TP=300/delt; %pattern length (ms/delt= steps)
for p1=1:N_mu

   try
        eval(['load Learning_BP/Data_Pattern/EP_',num2str(p1),'_',num2str(3),  '.mat PAT']);
        eval(['load Learning_BP/Weights/W_',num2str(p1),'_',num2str(3),  '.mat W_vec']);
    catch       
    end
    
    
    PATR(NI+1:NN,:) = rand(400,TT) < (rateE*delt/1000); 
    inall=zeros(NN,2*TT-1);
    PATR(1:NI,:) = (rand(NI,TT) < rateI*delt/1000);
   
    PATP = PATR;
    PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT(:,:,1); 
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

     I_ext=(INPI+INPE)';
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    S00V(p1,:)=num_spikes_vec;
    E00V(p1,:) = INPE;    
    I00V(p1,:) = INPI;

end

dfff=(-Lem+1:1:Lem)/10;

DF=sum(S00V(:));
ET=zeros(50*DF,numel(dfff));
IT=zeros(50*DF,numel(dfff));
ghj=0;

for j=1:N_mu
    df=find(S00V(j,:)>0);
    ghj=ghj+numel(df);
    for kj=1:numel(df)
        if df(kj)<20000-Lem
            if df(kj)>Lem
            i=(j-1)*20+kj;
            df(kj);
            ET(i,:)=E00V(j,df(kj)-Lem:df(kj)+Lem-1);
            IT(i,:)=I00V(j,df(kj)-Lem:df(kj)+Lem-1);
            end
        end
    end
end

FGG=sum(ET,2);
gk=numel(find(FGG>0));

plot(dfff,sum(ET,1)/gk,'color',[0.5 0 0],'DisplayName','L_{em} = 300 ms ','Linewidth',2.5);
hold on
I3= plot(dfff,sum(IT,1)/gk,'color',[0.5 0 0],'DisplayName','L_{em} = 300 ms ','Linewidth',2.5);



ylim([-3 4.6])
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('Inhibitory            excitatory','FontSize',14)
xlabel('Lag (ms)','FontSize',14)

legend([I1 I2 I3 ],'Location','northwest','color','none' )
title('(b)')
set(gca,'FontSize',24,'color','none')

print('Fig2b','-dpng','-r300')
