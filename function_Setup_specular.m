function [Rh,HBar,Rf,fBar,RG_BS,RG_RIS,GBar,...
    channelGaindB_h,channelGaindB_f,channelGaindB_G,...
    probLOS_h,probLOS_f,ricianFactor_h,ricianFactor_f,ricianFactor_G,...
    Rf_b,fBar_b,RG_RIS_b,GBar_b,...
    channelGaindB_f_b,channelGaindB_G_b] =...
    function_Setup_specular(L,K,M,NHor,NVer,ASDazimDeg,ASDelevDeg,...
    SpecNum_h,SpecNum_G,SpecNum_f,DirectLoss,RISpositions,BaseDistHor,AreaHorSize,AreaVerSize,...
    probLOSbinary_h,probLOSbinary_f,ricianFactorLoss_f,AzimSpecDev,ElevSpecDev)


N = NHor*NVer;
%Standard deviation of shadow fading in dB
sigma_sf_NLOS = 10; %for NLOS
sigma_sf_LOS = 4;   %for LOS

verticalDistance = 10;
%Minimum 3D distance between BSs and UEs
minDistance = 20;

%Define the antenna spacing (in number of wavelengths) for BS
antennaSpacingBS = 1/2; %Half wavelength distance

%Define the antenna spacing (in number of wavelengths) for RIS
antennaSpacingRIS = 1/4;
antennaSpacingRIS_b = 1/4;

ASDazimRad = ASDazimDeg/180*pi;
ASDelevRad = ASDelevDeg/180*pi;

posRISY = vec(repmat(0:NHor-1,NVer,1)*antennaSpacingRIS);
posRISZ = vec(repmat((0:NVer-1).',1,NHor)*antennaSpacingRIS);

posRISY_b = vec(repmat(0:NHor-1,NVer,1)*antennaSpacingRIS_b);
posRISZ_b = vec(repmat((0:NVer-1).',1,NHor)*antennaSpacingRIS_b);

BSposition = 0;

UEpositions = zeros(K,1);

%Prepare to store normalized spatial covariance matrices and array steering vectors
%and phase shifts
Rh = zeros(M,M,K);
HBar = zeros(M,K,SpecNum_h);

Rf = zeros(N,N,K,L);
Rf_b = zeros(N,N,K,L);

fBar = zeros(N,K,L,SpecNum_f);
fBar_b = zeros(N,K,L,SpecNum_f);

RG_BS = zeros(M,M,L);

RG_RIS = zeros(N,N,L);
RG_RIS_b = zeros(N,N,L);

GBar = zeros(M,N,L,SpecNum_G);
GBar_b = zeros(M,N,L,SpecNum_G);

%The maximum distance for LOS
maxdistLOS = 300;

%Prepare to store channel gain (in dB), Rician Factor \kappa and
%probabaility of LoS for each UE

channelGaindB_h = zeros(K,1);
channelGaindB_f = zeros(K,L);

probLOS_h = zeros(K,1);
probLOS_f = zeros(K,L);

ricianFactor_h = zeros(K,1);
ricianFactor_f = zeros(K,L);

distancesBS_RIS = abs(BSposition - RISpositions);
ricianFactor_G = db2pow(13-0.03*distancesBS_RIS);
channelGaindB_G =  -30.18-26*log10(distancesBS_RIS) + sigma_sf_LOS*randn(L,1);

%Put out UEs in the cells randomly
UEcounter = 0;
while UEcounter < K

    posXY = BaseDistHor + AreaHorSize*rand - AreaVerSize/2*1i + AreaVerSize*1i*rand;
    distancesBS_UE = sqrt(abs(posXY - BSposition)^2 + verticalDistance^2);
    distancesRIS_UE = sqrt(abs(posXY - RISpositions).^2 + verticalDistance^2);


    if min(distancesBS_UE, min(distancesRIS_UE)) >= minDistance

        UEcounter = UEcounter + 1;
        UEpositions(UEcounter) = posXY;

        if probLOSbinary_h == 1
            probLOS_h(UEcounter) = 1;
        elseif probLOSbinary_h == 0
            probLOS_h(UEcounter) = 0;
        else
            probLOS_h(UEcounter) = rand<((maxdistLOS-distancesBS_UE)/maxdistLOS);
        end

        if probLOSbinary_f == 1
            probLOS_f(UEcounter,:) = 1;
        elseif probLOSbinary_f == 0
            probLOS_f(UEcounter,:) = 0;
        else
            probLOS_f(UEcounter,:) = rand(1,L)<((maxdistLOS-distancesRIS_UE).'/maxdistLOS);
        end


        ricianFactor_h(UEcounter) = db2pow(13-0.03*distancesBS_UE);
        ricianFactor_f(UEcounter,:) = db2pow(13-0.03*distancesRIS_UE-ricianFactorLoss_f);

        if probLOS_h(UEcounter)==1
            channelGaindB_h(UEcounter) =  -30.18-26*log10(distancesBS_UE) + sigma_sf_LOS*randn;
        else
            channelGaindB_h(UEcounter) = -34.53-38*log10(distancesBS_UE) + sigma_sf_NLOS*randn;

        end

        for ell = 1:L

            if probLOS_f(UEcounter,ell)==1
                channelGaindB_f(UEcounter,ell) =  -30.18-26*log10(distancesRIS_UE(ell)) + sigma_sf_LOS*randn;
            else
                channelGaindB_f(UEcounter,ell) = -34.53-38*log10(distancesRIS_UE(ell)) + sigma_sf_NLOS*randn;
            end
        end

        azimBS_UE = angle(UEpositions(UEcounter)-BSposition);
        elevBS_UE = -acos(abs(UEpositions(UEcounter)-BSposition)...
            /sqrt(verticalDistance^2+abs(UEpositions(UEcounter)-BSposition)^2));
        for ss = 1:SpecNum_h
            if (probLOS_h(UEcounter) == 1)&&(ss==1)
                HBar(:,UEcounter,ss) = exp(1i*2*pi.*(0:(M-1))*cos(elevBS_UE)*sin(azimBS_UE)*antennaSpacingBS);
            else
                elevAng = elevBS_UE-ElevSpecDev+2*ElevSpecDev*rand;
                azimAng = azimBS_UE-AzimSpecDev+2*AzimSpecDev*rand;
                HBar(:,UEcounter,ss) = exp(1i*2*pi.*(0:(M-1))*cos(elevAng)*sin(azimAng)*antennaSpacingBS);
            end
        end

        Rh(:,:,UEcounter) = functionRlocalscatteringBS(M,azimBS_UE,elevBS_UE,ASDazimRad,ASDelevRad,antennaSpacingBS);

        for ell = 1:L

            azimRIS_UE = angle(UEpositions(UEcounter)-RISpositions(ell))-pi/2;
            elevRIS_UE = -acos(abs(UEpositions(UEcounter)-RISpositions(ell))...
                /sqrt(verticalDistance^2+abs(UEpositions(UEcounter)-RISpositions(ell))^2));
            for ss = 1:SpecNum_f
                if (probLOS_f(UEcounter,ell) == 1)&&(ss==1)
                    fBar(:,UEcounter,ell,ss) = exp(1i*2*pi*(posRISY*cos(elevRIS_UE)*sin(azimRIS_UE)+posRISZ*sin(elevRIS_UE)));
                    fBar_b(:,UEcounter,ell,ss) = exp(1i*2*pi*(posRISY_b*cos(elevRIS_UE)*sin(azimRIS_UE)+posRISZ_b*sin(elevRIS_UE)));

                else
                    elevAng = elevRIS_UE-ElevSpecDev+2*ElevSpecDev*rand;
                    azimAng = azimRIS_UE-AzimSpecDev+2*AzimSpecDev*rand;
                    fBar(:,UEcounter,ell,ss) = exp(1i*2*pi*(posRISY*cos(elevAng)*sin(azimAng)+posRISZ*sin(elevAng)));
                    fBar_b(:,UEcounter,ell,ss) = exp(1i*2*pi*(posRISY_b*cos(elevAng)*sin(azimAng)+posRISZ_b*sin(elevAng)));

                end
            end


            Rf(:,:,UEcounter,ell) = functionRlocalscatteringRIS(NHor,NVer,azimRIS_UE,elevRIS_UE,ASDazimRad,ASDelevRad,antennaSpacingRIS);
            Rf_b(:,:,UEcounter,ell) = functionRlocalscatteringRIS(NHor,NVer,azimRIS_UE,elevRIS_UE,ASDazimRad,ASDelevRad,antennaSpacingRIS_b);

            azimRIS_BS = angle(RISpositions(ell)-BSposition);
            azimRIS_BS2 = azimRIS_BS + pi/2;
            elevRIS_BS = 0;
            GBar(:,:,ell,1) = exp(1i*2*pi*(0:(M-1)).'*cos(elevRIS_BS)*sin(azimRIS_BS)*antennaSpacingBS)...
                *(exp(1i*2*pi*(posRISY*cos(elevRIS_BS)*sin(azimRIS_BS2)+posRISZ*sin(elevRIS_BS)))).';

            GBar_b(:,:,ell,1) = exp(1i*2*pi*(0:(M-1)).'*cos(elevRIS_BS)*sin(azimRIS_BS)*antennaSpacingBS)...
                *(exp(1i*2*pi*(posRISY_b*cos(elevRIS_BS)*sin(azimRIS_BS2)+posRISZ_b*sin(elevRIS_BS)))).';

            for ss = 2:SpecNum_G
                elevSign = randn;
                azimSign = randn;
                if elevSign > 0
                    elevAng = elevRIS_BS + ElevSpecDev*rand;
                    elevAng2 = elevRIS_BS + ElevSpecDev*rand;
                else
                    elevAng = elevRIS_BS - ElevSpecDev*rand;
                    elevAng2 = elevRIS_BS - ElevSpecDev*rand;
                end

                if azimSign > 0
                    azimAng = azimRIS_BS + AzimSpecDev*rand;
                    azimAng2 = azimRIS_BS2 - AzimSpecDev*rand;
                else
                    azimAng = azimRIS_BS - AzimSpecDev*rand;
                    azimAng2 = azimRIS_BS2 + AzimSpecDev*rand;
                end

                GBar(:,:,ell,ss) = exp(1i*2*pi*(0:(M-1)).'*cos(elevAng)*sin(azimAng)*antennaSpacingBS)...
                    *(exp(1i*2*pi*(posRISY*cos(elevAng2)*sin(azimAng2)+posRISZ*sin(elevAng2)))).';

                GBar_b(:,:,ell,ss) = exp(1i*2*pi*(0:(M-1)).'*cos(elevAng)*sin(azimAng)*antennaSpacingBS)...
                    *(exp(1i*2*pi*(posRISY_b*cos(elevAng2)*sin(azimAng2)+posRISZ_b*sin(elevAng2)))).';

            end

            RG_BS(:,:,ell) = functionRlocalscatteringBS(M,azimRIS_BS,elevRIS_BS,ASDazimRad,ASDelevRad,antennaSpacingBS);

            RG_RIS(:,:,ell) = functionRlocalscatteringRIS(NHor,NVer,azimRIS_BS2,elevRIS_BS,ASDazimRad,ASDelevRad,antennaSpacingRIS);
            RG_RIS_b(:,:,ell) = functionRlocalscatteringRIS(NHor,NVer,azimRIS_BS2,elevRIS_BS,ASDazimRad,ASDelevRad,antennaSpacingRIS_b);

        end

    end
end
channelGaindB_f_b = channelGaindB_f + 10*log10(antennaSpacingRIS_b^2*4*pi);
channelGaindB_G_b = channelGaindB_G + 10*log10(antennaSpacingRIS_b^2*4*pi);

channelGaindB_f = channelGaindB_f + 10*log10(antennaSpacingRIS^2*4*pi);
channelGaindB_G = channelGaindB_G + 10*log10(antennaSpacingRIS^2*4*pi);

channelGaindB_h = channelGaindB_h - DirectLoss;