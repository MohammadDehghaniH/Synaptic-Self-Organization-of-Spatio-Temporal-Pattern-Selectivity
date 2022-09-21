clc; clear; close all

Fig8a_b=figure();

rng(2^23)

N_psn=1;     % number of post-synaptic neuron
resy = 0.0;  % resting potential
threshy = 1; % threshold potential
NN = 500;    % number of afferents
i_t=randi(500); %i'th trial (choose the i'th trial randomly.);

p1=randi(20);

percI = 0.2; % inhibitory percentages
NI = fix(percI*NN);  % number of inhibitory neurons
T = 1000; %total time [all in ms]
delt = .1; % integration step.
TT = T/delt;
times = 1:TT;
times = times*delt;

% Generating kernels for excitatory and inhibitory inputs

taumem = 15;  %membrane time constant
tautraceE =3; %Rise time of exctatoy currents
tautraceI = 5;%Decay time of inhibitory currents
tauriseE =0.5;%Decay time of excitatoy currents
tauriseI = 1; %Rise time of inhibitory currents

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

TP = 50/delt; %pattern length (ms/delt= steps)
intimeP_vec = [4000 6000]; %ms Inserting time of the embedded pattern ()

inall=zeros(NN,2*TT-1);
tr=200;%transient time
dta=delt/taumem;


eval(['load Learning_B/Weights/2W_',num2str(i_t),'_',num2str(1),  '.mat W_vec']);
W_vecI=W_vec;

eval(['load Learning_BP/Data_Pattern/2EP_',num2str(i_t),'_',num2str(1),  '.mat PAT']);
eval(['load Learning_BP/Weights/2W_',num2str(i_t),'_',num2str(1),  '.mat W_vec']);


PATa=PAT(:,:,1);
PATb=PAT(:,:,2);
subplot(2,1,1)

%random pattern
BCV1=rand(NN,TT);
PATP=zeros(NN,TT);
PATP(NI+1:NN,:) = BCV1(NI+1:NN,:) < (rateE*delt/1000); 
PATP(1:NI,:) = BCV1(1:NI,:) < (rateI*delt/1000);   
PATP(:, intimeP_vec(1)+1: intimeP_vec(1) + TP) = PATa; 
PATP(:, intimeP_vec(2)+1: intimeP_vec(2) + TP) = PATb; 

[y,x]=find(PATP>0);
plot(x/10000,y,'k.','MarkerSize',2.5)
hold on
[y1,x1]=find(PATb>0);
plot(.6+(x1/10000),y1,'b.','MarkerSize',2.5)

hold on
[y1,x1]=find(PATa>0);
plot(0.4+(x1/10000),y1,'r.','MarkerSize',2.5)
axis off
box off
xlim([0.3 0.7])
ylabel('input','Visible','on','FontSize', 14);
tz=title('(a)');
get(tz,'position')
set(tz,'position',[0.3  510.3125   0]);
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;      
set(gca,'FontSize',14,'color','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BCV1=rand(NN,TT);
WW_vec = W_vecI;
WW_vec(1:NI,:) = -W_vecI(1:NI,:);
for i = 1: NI
    inall(i,:) = conv(PATP(i,:), traceI);        
end 
for i = NI+1: NN
    inall(i,:) = conv(PATP(i,:), traceE);
end    
INP = inall(:,1:TT);   
I_ext(:,1)=WW_vec(:,1)'*INP; 

V_vec1=zeros(TT,N_psn);
for it = 2:TT
    V_vec1(it,:) = ((1- dta)*V_vec1(it-1,:)) + I_ext(it,:)*dta;
    if V_vec1(it-1,:) > threshy
        V_vec1(it,1)=resy;
    end
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xt=(1:1:TT)/10000;

 
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
V_vec=zeros(TT,N_psn);
for it = 2:TT
    V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
    if V_vec(it-1,:) > threshy
        V_vec(it,1)=resy;
    end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

h3=subplot(2,1,2);
restingPotential=zeros(1,numel(V_vec));
threshold=ones(1,numel(V_vec));
plot(xt,restingPotential,'g','DisplayName','resting potential','Linewidth',2.);
hold on
plot(xt,threshold,'r','DisplayName','threshold','Linewidth',2.);
hold on
plot(xt,V_vec1,'b','DisplayName','before learning','Linewidth',2.);
hold on
plot(xt,V_vec,'k','DisplayName','after learning','Linewidth',2.);
xlim([0.3 0.7])
ylabel('V(t)','FontSize', 14)
ylim([-1 1.2])
xlabel('time (sec)','FontSize', 14)
tz=title('(b)');
get(tz,'position')
set(tz,'position',[0.3    1.5         0]);
box off
box off

% Create line
annotation(Fig8a_b,'line',[0.71 0.71],...
    [0.92 0.15]);
annotation(Fig8a_b,'line',[0.81 0.81],...
    [0.92 0.15]);
annotation(Fig8a_b,'line',[0.32 0.32],...
    [0.92 0.15]);
annotation(Fig8a_b,'line',[0.42 0.42],...
    [0.92 0.15]);


h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;      
set(gca,'FontSize',14,'color','none');
ttl=title('(b)','FontSize', 14);

print('Fig8ab','-dpng','-r300')



