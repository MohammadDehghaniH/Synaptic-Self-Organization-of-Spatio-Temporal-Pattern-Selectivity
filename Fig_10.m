clc;clear;close all;

rng(2^28)

Fig10 = figure();

N_mu=50; %There are \mu independent embedded patterns, 500 simulations (number of simulations)
N_psn=30; %number of post-synaptic neurons
em=10; % number of embedded patterns

delt = 0.1; %integration steps 
TP=50/delt; %pattern length (ms/delt= steps)
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
inall=zeros(NN,2*TT-1);

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

tr=200;%transient time
dta=delt/taumem;


dis_pat=zeros(N_psn,N_mu); %distribution of spikes

intimeP_vec=(1./delt)*(85:90:895); %ms Inserting time of the embedded patterns

for ic=1:N_mu
    try
        eval(['load ',num2str(N_psn,'%.i') '/Wa_',num2str(ic),'_',num2str(N_psn),  '.mat W_veca']);
        eval(['load ',num2str(N_psn,'%.i') '/Wb_',num2str(ic),'_',num2str(N_psn),  '.mat W_vecb']);
        eval(['load ',num2str(N_psn,'%.i') '/EP_',num2str(ic),'_',num2str(N_psn),  '.mat PAT']);
    catch       
    end
        
    W_vec=W_veca.*W_vecb;
    PATR(1:NI,:) = 1.*(rand(NI,TT) < rateI*delt/1000);
    PATR(NI+1:NN,:) = 1.*(rand(NN-NI,TT) < rateE*delt/1000); 
    PATP = PATR;
    for i=1:em
         PATP(:, intimeP_vec(1,i)+1: intimeP_vec(1,i) + TP) = PAT(:,:,i);  
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
    I_ext=zeros(TT,N_psn);
    for iw = 1:N_psn
        I_ext(:,iw)=WW_vec(:,iw)'*INP;
    end
    %Neuron model
    V_vec=zeros(TT,N_psn);
    for it = 2:TT
        V_vec(it,:) = (1- delt/taumem)*(V_vec(it-1,:)) + I_ext(it-1,:)*delt/taumem;
        f_p = find(V_vec(it-1,:) >= threshy);
        if ~isempty(f_p)
            V_vec(it,f_p) =resy; 
        end
    end
    V_vec=V_vec>=1;
    OMEGA_2=zeros(N_psn,em);
    for il=1:N_psn
       for KL=1:em
            OMEGA_2(il,KL)=sign(sum(V_vec(intimeP_vec(KL):intimeP_vec(KL)+TP+150,il)));
       end
    end
    dis_pat(:,ic)=sort(sum(OMEGA_2,2));
   
end



%distribution
YYY=zeros(10,N_mu);
for i=1:N_mu
    for j=1:em
        YYY(j,i)=numel(find(dis_pat(:,i)==j));
    end
end

yu=mean(YYY,2);
eyu=zeros(1,em);

%errorbar
for i=1:10
    eyu(1,i)=std(YYY(i,:))/sqrt(N_mu);
end
n_x=1:1:em;
bar(n_x,yu,0.33,'EdgeColor',[0 0.5 0],'faceColor',[0 0.5 0])
hold on
errorbar(n_x,yu,eyu,'color',[0.5 0 0],'LineStyle','none', 'LineWidth', 3.3);    


h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel('Neurons','FontSize',14)
xlabel('Number of patterns','FontSize',14)
set(gca,'FontSize',14,'color','none')
ylim([0 20])
box off
