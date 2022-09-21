function [I_ext,INP]=Model_Input_NET(rateI,delt,rateE,intimeP_vec,W_vec,traceI,traceE,NN,TT,Npsn,NI,em,PAT,TP,emm);
    %I_ext: external input
    %INP: corresponding kernel to external input
    
    inall=zeros(NN,2*TT-1);
    I_ext=zeros(TT,Npsn);
    PATP=rand(NN,TT);
    PATP(NI+1:NN,:) = PATP(NI+1:NN,:) < (rateE*delt/1000); 
    PATP(1:NI,:) = PATP(1:NI,:) < (rateI*delt/1000);  
    
    if emm>0
        for i=1:em
             PATP(:, intimeP_vec(1,i)+1: intimeP_vec(1,i) + TP) = PAT(:,:,i);  
        end
    end
    for i = 1: NI
        inall(i,:) = -1*conv(PATP(i,:), traceI);        
    end    
    for i = NI+1: NN
        inall(i,:) = conv(PATP(i,:), traceE);
    end
    INP = abs(inall(:,1:TT)); 
    for iw = 1:Npsn
        I_ext(:,iw)=W_vec(:,iw)'*inall(:,1:TT);
    end
end

