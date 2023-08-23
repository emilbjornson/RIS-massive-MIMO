% This Matlab script generates Figure 8 in the paper:
%
% Özlem Tugfe Demir and Emil Björnson,
% "Is Channel Estimation Necessary to Select Phase-Shifts for
% RIS-Assisted Massive MIMO?,"
% IEEE Transactions on Wireless Communications, vol. 21, no. 11, November
% 2022.
%
% This is version 1.0 (Last edited: 2023-08-23)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.

for scenario = 1:2
    %Number of RISs
    L = 2;

    %Number of vertical and horizontal RIS units per RIS
    NVer = 16;
    NHor = 16;

    PorVer = 4;
    PorHor = 4;

    NVerUE = 16;
    NHorUE = 16;

    %Number of UEs
    K = 8;

    %Number of BS antennas
    M = 100;

    %Select the number of setups with random UE locations
    nbrOfSetups = 600;

    %Number of channel realizations per setup
    nbrOfRealizations = 1000;

    %Total uplink pilot transmit power per UE (W)
    etaa = 0.1;


    %Communication bandwidth
    B = 1e6;

    %Noise figure (in dB)
    noiseFigure = 7;

    %Compute noise power
    noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

    %Noise power (dB)
    noiseVariancedB = noiseVariancedBm-30;


    %Select length of coherence block
    tau_c = 10000;

    %Angular standard deviation per path in the local scattering model (in degrees)
    ASDazimDeg = 15;
    ASDelevDeg = 15;

    %Number of specular components with LOS component, h, G, f
    if scenario ==1
        SpecNum_h = 1;
        SpecNum_G = 1;
        SpecNum_f = 1;
    else
        SpecNum_h = 3;
        SpecNum_G = 3;
        SpecNum_f = 3;
    end

    %Power ratio of LOS to other specular components
    if scenario ==1
        LOStoSpecPow_h = 1;
        LOStoSpecPow_G = 1;
        LOStoSpecPow_f = 1;
    else
        LOStoSpecPow_h = 0.5;
        LOStoSpecPow_G = 0.5;
        LOStoSpecPow_f = 0.5;
    end





    %Direct link loss (dB)
    DirectLoss = 0;

    %Total pilot length for RIS
    tau_p_RIS = (L*PorHor*PorVer+1)*K;

    tau_coef = (L*PorHor*PorVer+1);
    tau_p_RIS2 = tau_coef*K;

    %Total pilot length for conventional operation without RIS
    tau_coef_Conv = 28;
    tau_p_Conv = tau_coef_Conv*K;

    %Prelog factor for RIS
    prelogFactor_RIS1 = (tau_c - tau_p_RIS)/tau_c;
    prelogFactor_RIS2 = (tau_c - tau_p_RIS2)/tau_c;
    prelogFactor_RIS3 = (tau_c - tau_p_RIS - tau_p_RIS2)/tau_c;

    %Prelog factor for conventional operation without RIS
    prelogFactor_Conv = (tau_c - tau_p_Conv)/tau_c;

    SE_RIS_MR1 = zeros(K,nbrOfSetups,5);
    SE_RIS_MR2b = zeros(K,nbrOfSetups,5);

    SE_Conv_MR = zeros(K,nbrOfSetups);

    SE_RIS_RZF1 = zeros(K,nbrOfSetups,5);
    SE_RIS_RZF2b = zeros(K,nbrOfSetups,5);

    SE_Conv_RZF = zeros(K,nbrOfSetups);

    SE_RIS_AMMSE1 = zeros(K,nbrOfSetups,5);
    SE_RIS_AMMSE2b = zeros(K,nbrOfSetups,5);

    SE_Conv_AMMSE = zeros(K,nbrOfSetups);

    SE_RIS_MR_maxmin1 = zeros(K,nbrOfSetups,5);
    SE_RIS_MR_maxmin2b = zeros(K,nbrOfSetups,5);

    SE_Conv_MR_maxmin = zeros(K,nbrOfSetups);

    SE_RIS_RZF_maxmin1 = zeros(K,nbrOfSetups,5);
    SE_RIS_RZF_maxmin2b = zeros(K,nbrOfSetups,5);

    SE_Conv_RZF_maxmin = zeros(K,nbrOfSetups);

    SE_RIS_AMMSE_maxmin1 = zeros(K,nbrOfSetups,5);
    SE_RIS_AMMSE_maxmin2b = zeros(K,nbrOfSetups,5);

    SE_Conv_AMMSE_maxmin = zeros(K,nbrOfSetups);


    %%%%%%%%%%%%%%%%


    RISpositions = zeros(L,1);
    RISpositions(1) = 200+50*1i;
    RISpositions(2) = 200-50*1i;


    %Base UE distance to BS
    BaseDistHor = 200;

    %UE dropping area
    AreaHorSize = 100;
    AreaVerSize = 100;

    %LOS probability (it is 0.5 if it will be determined by the formula)
    probLOSbinary_f = 0.5;
    probLOSbinary_h = 0.5;

    %Rician factor loss for RIS-UE channels compared to the original formula
    %(dB)
    ricianFactorLoss_f = 0;


    %Azimuth and elevation angle deviations for the specular components
    %(radians)
    AzimSpecDev = 60/180*pi;
    ElevSpecDev = 15/180*pi;

    %% Go through all setups
    for setupp = 1:nbrOfSetups
        disp(setupp)

        %Power distribution among specular components other than LOS
        PowDistSpec_h = rand(SpecNum_h-1,K);
        PowDistSpec_f = rand(SpecNum_f-1,K,L);
        PowDistSpec_G = rand(SpecNum_G-1,L);

        %R and HMean is normalized
        [Rh,HBar,Rf,fBar,RG_BS,RG_RIS,GBar,...
            channelGaindB_h,channelGaindB_f,channelGaindB_G,...
            probLOS_h,probLOS_f,ricianFactor_h,ricianFactor_f,ricianFactor_G,...
            Rf_b,fBar_b,RG_RIS_b,GBar_b,...
            channelGaindB_f_b,channelGaindB_G_b] =...
            function_Setup_specular(L,K,M,NHor,NVer,ASDazimDeg,ASDelevDeg,...
            SpecNum_h,SpecNum_G,SpecNum_f,DirectLoss,RISpositions,BaseDistHor,AreaHorSize,AreaVerSize,...
            probLOSbinary_h,probLOSbinary_f,ricianFactorLoss_f,AzimSpecDev,ElevSpecDev);

        channelGaindB_h = channelGaindB_h - noiseVariancedB;
        channelGaindB_f = channelGaindB_f - noiseVariancedB;
        channelGaindB_f_b = channelGaindB_f_b - noiseVariancedB;


        poww = etaa*ones(K,1);




        RISassignments = functionRISassignment(channelGaindB_h,channelGaindB_f,channelGaindB_G,...
            L,K,NVer,NHor,NVerUE,NHorUE);

        [Rh2,HBar2,H,Hhat_Conv,bb1,bHat1,Ctilde1,bb2,bHat2,Ctilde2,bb3,bHat3,Ctilde3,...
            bb4,bHat4,Ctilde4,bb5,bHat5,Ctilde5] = ...
            functionChannelEstimationRIS_individual(Rh,HBar,Rf,fBar,RG_BS,RG_RIS,GBar,...
            channelGaindB_h,channelGaindB_f,channelGaindB_G,...
            probLOS_h,probLOS_f,ricianFactor_h,ricianFactor_f,ricianFactor_G,...
            nbrOfRealizations,L,K,M,etaa,NVer,NHor,NVerUE,NHorUE,RISassignments,tau_p_Conv,PorVer,PorHor,...
            SpecNum_h,SpecNum_G,SpecNum_f,LOStoSpecPow_h,LOStoSpecPow_G,LOStoSpecPow_f,...
            PowDistSpec_h,PowDistSpec_f,PowDistSpec_G);


        [SE_MR0, SE_RZF0, SE_AMMSE0, SE_MR_maxmin0, SE_RZF_maxmin0, SE_AMMSE_maxmin0 ] = ...
            functionComputeSEConv(Rh2,HBar2,H,Hhat_Conv,K,M,tau_p_Conv,etaa,nbrOfRealizations,poww,SpecNum_h);

        SE_Conv_MR(:,setupp) = prelogFactor_Conv*SE_MR0;
        SE_Conv_RZF(:,setupp) = prelogFactor_Conv*SE_RZF0;
        SE_Conv_AMMSE(:,setupp) = prelogFactor_Conv*SE_AMMSE0;
        SE_Conv_MR_maxmin(:,setupp) = prelogFactor_Conv*SE_MR_maxmin0;
        SE_Conv_RZF_maxmin(:,setupp) = prelogFactor_Conv*SE_RZF_maxmin0;
        SE_Conv_AMMSE_maxmin(:,setupp) = prelogFactor_Conv*SE_AMMSE_maxmin0;



        [SE_MR1, SE_RZF1, SE_AMMSE1, SE_MR_maxmin1, SE_RZF_maxmin1, SE_AMMSE_maxmin1 ] = functionComputeSERIS(bb1,bHat1,Ctilde1,K,M,nbrOfRealizations,poww,etaa);

        SE_RIS_MR1(:,setupp,1) = prelogFactor_RIS1*SE_MR1;
        SE_RIS_RZF1(:,setupp,1) = prelogFactor_RIS1*SE_RZF1;
        SE_RIS_AMMSE1(:,setupp,1) = prelogFactor_RIS1*SE_AMMSE1;
        SE_RIS_MR_maxmin1(:,setupp,1) = prelogFactor_RIS1*SE_MR_maxmin1;
        SE_RIS_RZF_maxmin1(:,setupp,1) = prelogFactor_RIS1*SE_RZF_maxmin1;
        SE_RIS_AMMSE_maxmin1(:,setupp,1) = prelogFactor_RIS1*SE_AMMSE_maxmin1;




        %%%%%%%%%%%%%%%

        [bb1,bHat1,Ctilde1,bb2,bHat2,Ctilde2,bb3,bHat3,Ctilde3,...
            bb4,bHat4,Ctilde4,bb5,bHat5,Ctilde5]  =  ...
            functionChannelEstimationRIS_cascadedDominant(Rh,HBar,Rf,fBar,RG_BS,RG_RIS,GBar,...
            channelGaindB_h,channelGaindB_f,channelGaindB_G,...
            probLOS_h,probLOS_f,ricianFactor_h,ricianFactor_f,ricianFactor_G,...
            nbrOfRealizations,L,K,M,etaa,NVer,NHor,NVerUE,NHorUE,RISassignments,tau_p_RIS2,...
            SpecNum_h,SpecNum_G,SpecNum_f,LOStoSpecPow_h,LOStoSpecPow_G,LOStoSpecPow_f,...
            PowDistSpec_h,PowDistSpec_f,PowDistSpec_G);

        [SE_MR2, SE_RZF2, SE_AMMSE2, SE_MR_maxmin2, SE_RZF_maxmin2, SE_AMMSE_maxmin2 ] = functionComputeSERIS(bb1,bHat1,Ctilde1,K,M,nbrOfRealizations,poww,etaa);

        SE_RIS_MR2b(:,setupp,1) = prelogFactor_RIS2*SE_MR2;
        SE_RIS_RZF2b(:,setupp,1) = prelogFactor_RIS2*SE_RZF2;
        SE_RIS_AMMSE2b(:,setupp,1) = prelogFactor_RIS2*SE_AMMSE2;
        SE_RIS_MR_maxmin2b(:,setupp,1) = prelogFactor_RIS2*SE_MR_maxmin2;
        SE_RIS_RZF_maxmin2b(:,setupp,1) = prelogFactor_RIS2*SE_RZF_maxmin2;
        SE_RIS_AMMSE_maxmin2b(:,setupp,1) = prelogFactor_RIS2*SE_AMMSE_maxmin2;


        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%



    end
    if scenario==1
        ConvA =  vec(SE_Conv_RZF_maxmin);
        RIS1A = vec(SE_RIS_RZF_maxmin1(:,:,1));
        RIS2bA = vec(SE_RIS_RZF_maxmin2b(:,:,1));
    else
        ConvB =  vec(SE_Conv_RZF_maxmin);
        RIS1B = vec(SE_RIS_RZF_maxmin1(:,:,1));
        RIS2bB = vec(SE_RIS_RZF_maxmin2b(:,:,1));
    end
end


%% Figure 8

nbrOfPoints = length(ConvA);

figure;
hold on; box on; grid on;
set(gca,'fontsize',16);

ppp1 = plot(sort(ConvA),linspace(0,1,nbrOfPoints),'k:','LineWidth',2);
ppp2 = plot(sort(RIS1A),linspace(0,1,nbrOfPoints),'k:o','LineWidth',2);
ppp3 = plot(sort(RIS2bA),linspace(0,1,nbrOfPoints),'k-','LineWidth',2);

ppp4 = plot(sort(ConvB),linspace(0,1,nbrOfPoints),'b:','LineWidth',3);
ppp5 = plot(sort(RIS1B),linspace(0,1,nbrOfPoints),'b:o','LineWidth',3);
ppp6 = plot(sort(RIS2bB),linspace(0,1,nbrOfPoints),'b-','LineWidth',3);

ppp2.MarkerSize = 6;
ppp2.MarkerIndices = 1:ceil(nbrOfPoints/7):nbrOfPoints;


ppp5.MarkerSize = 6;
ppp5.MarkerIndices = 1:ceil(nbrOfPoints/7):nbrOfPoints;




xlabel('SE per UE [b/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Conv-mMIMO $S=1$', 'RIS-mMIMO, Short-Term, $S=1$', 'RIS-mMIMO, Long-Term, $S=1$',  'Conv-mMIMO $S=3$', 'RIS-mMIMO, Short-Term, $S=3$', 'RIS-mMIMO, Long-Term, $S=3$' },'Interpreter','Latex','Location','SouthEast');
xlim([0 6]);
