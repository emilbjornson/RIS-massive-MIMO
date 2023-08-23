function   [Rh2,HBar2,H,Hhat_Conv,bb1,bHat1,Ctilde1,bb2,bHat2,Ctilde2,bb3,bHat3,Ctilde3,...
    bb4,bHat4,Ctilde4,bb5,bHat5,Ctilde5] = ...
    LS_functionChannelEstimationRIS_individual(Rh,HBar,Rf,fBar,RG_BS,RG_RIS,GBar,...
    channelGaindB_h,channelGaindB_f,channelGaindB_G,...
    probLOS_h,probLOS_f,ricianFactor_h,ricianFactor_f,ricianFactor_G,...
    nbrOfRealizations,L,K,M,etaa,NVer,NHor,NVerUE,NHorUE,RISassignments,tau_p,PorVer,PorHor,...
    SpecNum_h,SpecNum_G,SpecNum_f,LOStoSpecPow_h,LOStoSpecPow_G,LOStoSpecPow_f,...
    PowDistSpec_h,PowDistSpec_f,PowDistSpec_G)


N = NHor*NVer;

NPorHor = round(NHor/PorHor);
NPorVer = round(NVer/PorVer);

nbrOfUEperRISHor = floor(NHor/NHorUE);

Rh2 = zeros(M,M,K);

HBar2 = zeros(M,K,SpecNum_h);

Rf2 = zeros(N,N,K,L);

fBar2 = zeros(N,K,L,SpecNum_f);

RG_BS2 = zeros(M,M,L);

RG_RIS2 = zeros(N,N,L);

GBar2 = zeros(M,N,L,SpecNum_G);


%Prepare to store
channelGain_LOS_h = zeros(K,1);
channelGain_NLOS_h = zeros(K,1);

channelGain_LOS_f = zeros(K,L);
channelGain_NLOS_f = zeros(K,L);

channelGain_LOS_G = zeros(L,1);
channelGain_NLOS_G = zeros(L,1);

Hhat_Conv = zeros(M,K,nbrOfRealizations);

hh = zeros(M,nbrOfRealizations,K);

hHat = zeros(M,nbrOfRealizations,K);

GG = zeros(M,N,nbrOfRealizations,L);

ff = zeros(N,nbrOfRealizations,K,L);

Hhat_RIS = zeros(M,nbrOfRealizations,N,K,L);

%Go through all UEs and apply the channel gains to the spatial
%correlation matrices and mean vectors
for k = 1:K

    if probLOS_h(k)==1 %The LoS Path exists, Rician Factor ~= 0
        channelGain_LOS_h(k) = ricianFactor_h(k)/(ricianFactor_h(k) +1 )*db2pow(channelGaindB_h(k));
        channelGain_NLOS_h(k) = 1/(ricianFactor_h(k) +1 )*db2pow(channelGaindB_h(k));
    else  %Pure NLoS case
        channelGain_LOS_h(k) = 0;
        channelGain_NLOS_h(k) = db2pow(channelGaindB_h(k));
    end
    %Scaling operation
    HBar2(:,k,1) = sqrt(channelGain_LOS_h(k)*LOStoSpecPow_h)*HBar(:,k,1);
    TemporPowTotal = channelGain_LOS_h(k)*(1-LOStoSpecPow_h);
    TemporPow = PowDistSpec_h(:,k)/sum(PowDistSpec_h(:,k))*TemporPowTotal;
    for ss = 2:SpecNum_h
        HBar2(:,k,ss) = sqrt(TemporPow(ss-1))*HBar(:,k,ss);
    end
    Rh2(:,:,k) = channelGain_NLOS_h(k)*Rh(:,:,k);

    for ell = 1:L
        if probLOS_f(k,ell)==1 %The LoS Path exists, Rician Factor ~= 0
            channelGain_LOS_f(k,ell) = ricianFactor_f(k,ell)/(ricianFactor_f(k,ell) +1 )*db2pow(channelGaindB_f(k,ell));
            channelGain_NLOS_f(k,ell) = 1/(ricianFactor_f(k,ell) +1 )*db2pow(channelGaindB_f(k,ell));
        else  %Pure NLoS case
            channelGain_LOS_f(k,ell) = 0;
            channelGain_NLOS_f(k,ell) = db2pow(channelGaindB_f(k,ell));
        end


        Rf2(:,:,k,ell) = channelGain_NLOS_f(k,ell)*Rf(:,:,k,ell);

        fBar2(:,k,ell,1) = sqrt(channelGain_LOS_f(k,ell)*LOStoSpecPow_f)*fBar(:,k,ell,1);

        TemporPowTotal = channelGain_LOS_f(k,ell)*(1-LOStoSpecPow_f);
        TemporPow = PowDistSpec_f(:,k,ell)/sum(PowDistSpec_f(:,k,ell))*TemporPowTotal;
        for ss = 2:SpecNum_f
            fBar2(:,k,ell,ss) = sqrt(TemporPow(ss-1))*fBar(:,k,ell,ss);
        end

    end


end

for ell = 1:L

    channelGain_LOS_G(ell) = ricianFactor_G(ell)/(ricianFactor_G(ell) +1 )*db2pow(channelGaindB_G(ell));
    channelGain_NLOS_G(ell) = 1/(ricianFactor_G(ell) +1 )*db2pow(channelGaindB_G(ell));
    RG_BS2(:,:,ell) = RG_BS(:,:,ell);

    RG_RIS2(:,:,ell) = channelGain_NLOS_G(ell)*RG_RIS(:,:,ell);

    GBar2(:,:,ell,1) = sqrt(channelGain_LOS_G(ell)*LOStoSpecPow_G)*GBar(:,:,ell,1);
    TemporPowTotal = channelGain_LOS_G(ell)*(1-LOStoSpecPow_G);
    TemporPow = PowDistSpec_G(:,ell)/sum(PowDistSpec_G(:,ell))*TemporPowTotal;
    for ss = 2:SpecNum_G
        GBar2(:,:,ell,ss) = sqrt(TemporPow(ss-1))*GBar(:,:,ell,ss);
    end
end
%Generate the channel realizations
%Generate uncorrelated random variables
W_h = sqrt(0.5)*(randn(M,nbrOfRealizations,K)+1i*randn(M,nbrOfRealizations,K));
W_f = sqrt(0.5)*(randn(N,nbrOfRealizations,K,L)+1i*randn(N,nbrOfRealizations,K,L));
W_G = sqrt(0.5)*(randn(M,N,nbrOfRealizations,L)+1i*randn(M,N,nbrOfRealizations,L));



%Go through all UEs and generate channel realizations

for k = 1:K
    Rsqrt = sqrtm(Rh2(:,:,k));
    hh(:,:,k) = Rsqrt*W_h(:,:,k);
    for tt = 1:nbrOfRealizations
        for ss = 1:SpecNum_h
            temporr = exp(1i*2*pi*rand)*HBar2(:,k,ss);
            hh(:,tt,k) = hh(:,tt,k) + temporr;
        end
    end

    for ell = 1:L
        Rsqrt = sqrtm(Rf2(:,:,k,ell));
        ff(:,:,k,ell) = Rsqrt*W_f(:,:,k,ell);
        for tt = 1:nbrOfRealizations
            for ss = 1:SpecNum_f
                temporr = exp(1i*2*pi*rand)*fBar2(:,k,ell,ss);
                ff(:,tt,k,ell) = ff(:,tt,k,ell) + temporr;


            end

        end
    end
end

for ell = 1:L
    RsqrtBS = sqrtm(RG_BS2(:,:,ell));
    RsqrtRIS = sqrtm(RG_RIS2(:,:,ell));

    for tt = 1:nbrOfRealizations

        GG(:,:,tt,ell) = GBar2(:,:,ell,1) + RsqrtBS*W_G(:,:,tt,ell)*RsqrtRIS;
        for ss = 2:SpecNum_G
            temporr = exp(1i*2*pi*rand)*GBar2(:,:,ell,ss);
            GG(:,:,tt,ell) = GG(:,:,tt,ell) + temporr;


        end
    end
end


%%%++++%%%%&&&&!!!!!




pilCoef = sqrt((L*PorVer*PorHor+1)*K*etaa);

ReceivedPilot0 = sqrt(0.5)*(randn(M,nbrOfRealizations,K)...
    +1i*randn(M,nbrOfRealizations,K));
ReceivedPilot0 = ReceivedPilot0 + pilCoef*hh;

for k = 1:K
    hHat(:,:,k) = ReceivedPilot0(:,:,k)/pilCoef;
end

for pver = 1:PorVer
    for phor = 1:PorHor
        ReceivedPilotL = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L)...
            +1i*randn(M,nbrOfRealizations,K,L));
        rangeVer = (pver-1)*NPorVer+1:pver*NPorVer;
        rangeHor = (phor-1)*NPorHor:phor*NPorHor-1;

        posIndexx = vec(repmat(rangeVer.',1,NPorHor)+repmat(NVer*rangeHor,NPorVer,1));


        for ell = 1:L
            for k = 1:K

                for tt = 1:nbrOfRealizations
                    ReceivedPilotL(:,tt,k,ell) = ReceivedPilotL(:,tt,k,ell) + ...
                        pilCoef*GG(:,posIndexx,tt,ell)*ff(posIndexx,tt,k,ell);
                end




                posind = 0;
                for nn = posIndexx.'
                    posind = posind + 1;

                    Hhat_RIS(:,:,nn,k,ell) = ReceivedPilotL(:,:,k,ell)/pilCoef;

                end
            end

        end



    end
end



bb1 = permute(hh,[1 3 2]);
bHat1 = permute(hHat,[1 3 2]);

bb2 = permute(hh,[1 3 2]);
bHat2 = permute(hHat,[1 3 2]);

bb3 = permute(hh,[1 3 2]);
bHat3 = permute(hHat,[1 3 2]);

bb4 = permute(hh,[1 3 2]);
bHat4 = permute(hHat,[1 3 2]);

bb5 = permute(hh,[1 3 2]);
bHat5 = permute(hHat,[1 3 2]);

perRIS = zeros(L,1);

for jj = 1:K

    ell = RISassignments(jj);
    if ell>0
        whichPartition = perRIS(ell);
        perRIS(ell) = perRIS(ell)+1;
        whichHor = mod(whichPartition,nbrOfUEperRISHor);
        whichVer = floor(whichPartition/nbrOfUEperRISHor);
        posIndexx = vec((whichVer*NVerUE+1:(whichVer+1)*NVerUE).'+NVer*(whichHor*NHorUE:(whichHor+1)*NHorUE-1));


        for tt = 1:nbrOfRealizations

            hHatj = hHat(:,tt,jj);
            HHatPrimejj = reshape(Hhat_RIS(:,tt,posIndexx,jj,ell), [M,length(posIndexx)]);

            [anglee1, anglee2, anglee3] = RISPhaseShiftSelection(hHatj,HHatPrimejj,hHatj,HHatPrimejj);
            anglee4 = ones(length(posIndexx),1);
            anglee5 = exp(1i*2*pi*rand(length(posIndexx),1));

            for k = 1:K
                bb1(:,k,tt) = bb1(:,k,tt) + GG(:,posIndexx,tt,ell)*diag(anglee1)*ff(posIndexx,tt,k,ell);
                bb2(:,k,tt) = bb2(:,k,tt) + GG(:,posIndexx,tt,ell)*diag(anglee2)*ff(posIndexx,tt,k,ell);
                bb3(:,k,tt) = bb3(:,k,tt) + GG(:,posIndexx,tt,ell)*diag(anglee3)*ff(posIndexx,tt,k,ell);
                bb4(:,k,tt) = bb4(:,k,tt) + GG(:,posIndexx,tt,ell)*diag(anglee4)*ff(posIndexx,tt,k,ell);
                bb5(:,k,tt) = bb5(:,k,tt) + GG(:,posIndexx,tt,ell)*diag(anglee5)*ff(posIndexx,tt,k,ell);

                HHatPrimekj = reshape(Hhat_RIS(:,tt,posIndexx,k,ell), [M,length(posIndexx)]);

                bHat1(:,k,tt) = bHat1(:,k,tt) + HHatPrimekj*anglee1;
                bHat2(:,k,tt) = bHat2(:,k,tt) + HHatPrimekj*anglee2;
                bHat3(:,k,tt) = bHat3(:,k,tt) + HHatPrimekj*anglee3;
                bHat4(:,k,tt) = bHat4(:,k,tt) + HHatPrimekj*anglee4;
                bHat5(:,k,tt) = bHat5(:,k,tt) + HHatPrimekj*anglee5;

            end

        end
    end
end



Ctilde1 = zeros(M,M,K);
Ctilde2 = zeros(M,M,K);
Ctilde3 = zeros(M,M,K);
Ctilde4 = zeros(M,M,K);
Ctilde5 = zeros(M,M,K);

for k = 1:K
    temporMatr = permute(bb1(:,k,:)-bHat1(:,k,:),[1 3 2]);
    Ctilde1(:,:,k) = temporMatr*temporMatr'/nbrOfRealizations;

    temporMatr = permute(bb2(:,k,:)-bHat2(:,k,:),[1 3 2]);
    Ctilde2(:,:,k) = temporMatr*temporMatr'/nbrOfRealizations;

    temporMatr = permute(bb3(:,k,:)-bHat3(:,k,:),[1 3 2]);
    Ctilde3(:,:,k) = temporMatr*temporMatr'/nbrOfRealizations;

    temporMatr = permute(bb4(:,k,:)-bHat4(:,k,:),[1 3 2]);
    Ctilde4(:,:,k) = temporMatr*temporMatr'/nbrOfRealizations;

    temporMatr = permute(bb5(:,k,:)-bHat5(:,k,:),[1 3 2]);
    Ctilde5(:,:,k) = temporMatr*temporMatr'/nbrOfRealizations;
end



H = permute(hh,[1 3 2]);

for k = 1:K
    noisee = sqrt(0.5)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));
    ReceivedPilot = sqrt(tau_p*etaa)*reshape(H(:,k,:),M,nbrOfRealizations)+noisee;


    Hhat_Conv(:,k,:) = ReceivedPilot/sqrt(tau_p*etaa);

end

