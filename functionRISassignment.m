function RISassignments = functionRISassignment(channelGaindB_h,channelGaindB_f,channelGaindB_G,...
    L,K,NVer,NHor,NVerUE,NHorUE)


nbrOfUEperRISVer = floor(NVer/NVerUE);
nbrOfUEperRISHor = floor(NHor/NHorUE);

nbrOfUEperRIS = nbrOfUEperRISVer*nbrOfUEperRISHor;
RISassignments = zeros(K,1);
perRIS = zeros(L,1);

[~, sortedUEs] = sort(channelGaindB_h,'ascend');
for ii = 1:K
    k = sortedUEs(ii);
    [~, sortedRIS] = sort(channelGaindB_f(k,:).'+channelGaindB_G,'descend');
    for jj = 1:L
        ell = sortedRIS(jj);
        if perRIS(ell)<nbrOfUEperRIS
            RISassignments(k) = ell;
            perRIS(ell) = perRIS(ell)+1;
            break
        end
    end
end