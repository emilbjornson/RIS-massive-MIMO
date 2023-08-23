function [SE_MR, SE_RZF, SE_AMMSE, SE_MR_maxmin, SE_RZF_maxmin, SE_AMMSE_maxmin] = ...
    functionComputeSERIS(bb,bHat,Ctilde,K,M,nbrOfRealizations,poww,etaa)

SE_MR = zeros(K,1);
SE_RZF = zeros(K,1);
SE_AMMSE = zeros(K,1);

Cp = zeros(M,M,K);

for k = 1:K

    Cp(:,:,k) = poww(k)*Ctilde(:,:,k);
end


Numerr_MR = zeros(K,1);
Denomm_MR = zeros(K,K);

Numerr_RZF = zeros(K,1);
Denomm_RZF = zeros(K,K);


Numerr_AMMSE = zeros(K,1);
Denomm_AMMSE = zeros(K,K);

for tt = 1:nbrOfRealizations


    V_MR = bHat(:,:,tt);
    V_RZF = (V_MR*diag(poww)*V_MR'+eye(M))\V_MR;
    V_AMMSE = (V_MR*diag(poww)*V_MR'+sum(Cp,3)+eye(M))\V_MR;

    V_MR = V_MR./vecnorm(V_MR);
    V_RZF = V_RZF./vecnorm(V_RZF);
    V_AMMSE = V_AMMSE./vecnorm(V_AMMSE);

    TemporMR = V_MR'*bb(:,:,tt);
    Numerr_MR = Numerr_MR + diag(TemporMR)/nbrOfRealizations;
    Denomm_MR = Denomm_MR + abs(TemporMR).^2/nbrOfRealizations;

    TemporRZF = V_RZF'*bb(:,:,tt);
    Numerr_RZF = Numerr_RZF + diag(TemporRZF)/nbrOfRealizations;
    Denomm_RZF = Denomm_RZF + abs(TemporRZF).^2/nbrOfRealizations;

    TemporAMMSE = V_AMMSE'*bb(:,:,tt);
    Numerr_AMMSE = Numerr_AMMSE + diag(TemporAMMSE)/nbrOfRealizations;
    Denomm_AMMSE = Denomm_AMMSE + abs(TemporAMMSE).^2/nbrOfRealizations;


end

for k = 1:K

    SE_MR(k) = log2(1+poww(k)*abs(Numerr_MR(k))^2/...
        (Denomm_MR(k,:)*poww-poww(k)*abs(Numerr_MR(k))^2+1));
    SE_RZF(k) = log2(1+poww(k)*abs(Numerr_RZF(k))^2/...
        (Denomm_RZF(k,:)*poww-poww(k)*abs(Numerr_RZF(k))^2+1));
    SE_AMMSE(k) = log2(1+poww(k)*abs(Numerr_AMMSE(k))^2/...
        (Denomm_AMMSE(k,:)*poww-poww(k)*abs(Numerr_AMMSE(k))^2+1));
end


iter = 0;
etap = etaa*ones(K,1);
denominator = Denomm_MR*etap-etap.*abs(Numerr_MR).^2+1;
SINR = etap.*abs(Numerr_MR).^2./denominator;

while max(SINR)-min(SINR)>0.01
    iter = iter+1;


    etap = denominator./(abs(Numerr_MR).^2);
    etap = etap*etaa/max(etap);

    denominator = Denomm_MR*etap-etap.*(abs(Numerr_MR).^2)+1;

    SINR = etap.*abs(Numerr_MR).^2./denominator;

end

SE_MR_maxmin = log2(1+SINR);


iter = 0;
etap = etaa*ones(K,1);
denominator = Denomm_RZF*etap-etap.*abs(Numerr_RZF).^2+1;
SINR = etap.*abs(Numerr_RZF).^2./denominator;

while max(SINR)-min(SINR)>0.01
    iter = iter+1;


    etap = denominator./(abs(Numerr_RZF).^2);
    etap = etap*etaa/max(etap);

    denominator = Denomm_RZF*etap-etap.*(abs(Numerr_RZF).^2)+1;

    SINR = etap.*abs(Numerr_RZF).^2./denominator;

end

SE_RZF_maxmin = log2(1+SINR);



iter = 0;
etap = etaa*ones(K,1);
denominator = Denomm_AMMSE*etap-etap.*abs(Numerr_AMMSE).^2+1;
SINR = etap.*abs(Numerr_AMMSE).^2./denominator;

while max(SINR)-min(SINR)>0.01
    iter = iter+1;


    etap = denominator./(abs(Numerr_AMMSE).^2);
    etap = etap*etaa/max(etap);

    denominator = Denomm_AMMSE*etap-etap.*(abs(Numerr_AMMSE).^2)+1;

    SINR = etap.*abs(Numerr_AMMSE).^2./denominator;

end

SE_AMMSE_maxmin = log2(1+SINR);
