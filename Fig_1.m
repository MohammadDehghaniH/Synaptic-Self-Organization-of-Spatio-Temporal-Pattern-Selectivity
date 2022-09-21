clc; clear; close all

rng(1101)

Fig1=figure();

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
taumem = 15;  %membrane time constant
tautraceE =3; %Rise time of excitatoy currents
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

inall=zeros(NN,2*TT-1);
tr=200;%transient time
dta=delt/taumem;

subplot(2,1,1)    
[y,x]=find(PATP>0);

plot(x/10000,y,'k.','MarkerSize',3)

hold on
[y1,x1]=find(PAT>0);
plot(0.5+(x1/10000),y1,'r.','MarkerSize',3)
ylabel('input','Visible','on','FontSize', 14);
tz=title('(a)');
get(tz,'position')
set(tz,'position',[0.3  510.3125   -0.0000]);
box off
axis off
set(gca,'FontSize',14,'color','none');
xlim([0.3 0.7])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%weight vector
WW_vec = W_vec1;
WW_vec(1:NI,:) = -W_vec1(1:NI,:);
for i = 1: NI
    inall(i,:) = conv(PATP(i,:), traceI);        
end 
for i = NI+1: NN
    inall(i,:) = conv(PATP(i,:), traceE);
end    
INP = inall(:,1:TT);   
I_ext(:,1)=WW_vec(:,1)'*INP; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%Neuron Model
V_vec1=zeros(TT,N_psn);
for it = 2:TT
    V_vec1(it,:) = ((1- dta)*V_vec1(it-1,:)) + I_ext(it,:)*dta;
    if V_vec1(it-1,:) >= threshy
        V_vec1(it,1)=resy;
    end
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%Neuron Model
V_vec=zeros(TT,N_psn);
for it = 2:TT
    V_vec(it,:) = ((1- dta)*V_vec(it-1,:)) + I_ext(it,:)*dta;
    if V_vec(it-1,:) >= threshy
        V_vec(it,1)=resy;
    end
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


h3=subplot(2,1,2);
restingPotential=zeros(1,numel(V_vec));
threshold=ones(1,numel(V_vec));
plot(xt,restingPotential,'g','DisplayName','resting potential','Linewidth',2.5);
hold on
plot(xt,threshold,'r','DisplayName','threshold','Linewidth',2.5);
hold on
plot(xt,V_vec1,'b','DisplayName','before learning','Linewidth',2.5);
hold on
plot(xt,V_vec,'k','DisplayName','after learning','Linewidth',2.5);
ylabel('V(t)','FontSize', 14)
xlabel('time (sec)','FontSize', 14)
tz=title('(b)');
get(tz,'position')
set(tz,'position',[0.3    1.5         0]);
box off

% Create line

annotation(Fig1,'line',[0.519 0.519],...
    [0.93 0.155],'LineStyle','-');
annotation(Fig1,'line',[0.614 0.614],...
    [0.93 0.155],'LineStyle','-');
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;      
set(gca,'FontSize',14, 'color', 'none');

ylim([-1 1.2])
xlim([0.3 0.7])

print('Fig1','-dpng','-r300')


