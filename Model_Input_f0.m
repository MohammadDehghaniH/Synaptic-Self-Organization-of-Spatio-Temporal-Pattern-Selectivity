function [I_ext,INP]=Model_Input_f0(rateI,rateE,delt,intimeP_vec,TP,PAT,W_vec,inall,NI,NN,TT,traceI,traceE,em,itrial,f0)
    
    %I_ext: external input
    %INP: corresponding kernel to external input
    BCV1=rand(NN,TT);
    PATP=zeros(NN,TT);
    PATP(NI+1:NN,:) = BCV1(NI+1:NN,:) < (rateE*delt/1000); 
    PATP(1:NI,:) = BCV1(1:NI,:) < (rateI*delt/1000);  
    if em>0
        if rem(itrial,f0)==0
            PATP(:, intimeP_vec+1: intimeP_vec + TP) = PAT; 
        end
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
end
