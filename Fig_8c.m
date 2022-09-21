clc; clear; close all

Fig8c=figure();

r0_t = 2:1:9; %r_0 (desired number of spikes)
N_mu=100; %There are \mu independent embedded patterns, 500 simulations (number of simulations)

RESULTS=zeros(4,numel(r0_t));
for p1=1:numel(r0_t)
    SD=r0_t(p1);
    for p2=1:N_mu
        
        eval(['load Learning_BP_4em/S_',num2str(p2),'_',num2str(p1),  '.mat S']);

        if S==1
            RESULTS(1,p1)=RESULTS(1,p1)+1;
        end

        if S==2
            RESULTS(2,p1)=RESULTS(2,p1)+1;
        end

        if S==3
            RESULTS(3,p1)=RESULTS(3,p1)+1;
        end

        if S==4
             RESULTS(4,p1)=RESULTS(4,p1)+1;
        end
    end
end


I1=plot(r0_t,RESULTS(1,:),'-o','color',[0 0 0.5],'DisplayName','One Pattern','Linewidth',2.5);
hold on
I2=plot(r0_t,RESULTS(2,:),'-o','color',[0 0.5 0],'DisplayName','Two Patterns','Linewidth',2.5);
hold on
I3=plot(r0_t,RESULTS(3,:),'-o','color',[0.5 0 0],'DisplayName','Three Patterns','Linewidth',2.5);
hold on
I4=plot(r0_t,RESULTS(4,:),'-o','color',[0.5 0.5 0.5],'DisplayName','Four Patterns','Linewidth',2.5);

box off
h=gca;                   
h.LineWidth=2;            
hLg.LineWidth=1;            
ylabel(["number of simulations", "learns n patterns"],'FontSize',14)
xlabel('r_{0}','FontSize',14)
legend([I1 I2 I3 I4  ],'Location','northeast','color','none' )
title('(c)')
ylim([0 100])
set(gca,'FontSize',14,'color','none')

print('Fig8c','-dpng','-r300')

    


