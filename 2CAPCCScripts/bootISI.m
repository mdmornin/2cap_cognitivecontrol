function stMtx = bootISI(stMtx,p,refPeriod)
    % Given STMTX of shape TxN (Time by Neurons) adjust based on a value of
    % ISI violation occurence P
    % ISI Violations occur when neurons spike faster than their refractory
    % period
    % Here we are assuming the refractory period is 3 ms
    % A simple threshold will be used to determine if a neuron fires too
    % many times outside its refractory period.
    % is p greater than NTotalViolations / NTotalSpikes 
    TotalSpikes = sum(~isnan(stMtx));
    ISI = diff(stMtx);

    for k = 1:size(ISI,2)
        ISI_Violations(k) = sum(ISI(1:TotalSpikes(k)-1,k) < refPeriod);
    end
    
    Prop_Violations = ISI_Violations./TotalSpikes;
    Determination = Prop_Violations < p;
    
    stMtx = stMtx(:,Determination);

end