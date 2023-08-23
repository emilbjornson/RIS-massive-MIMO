function [anglee1, anglee2, anglee3] = RISPhaseShiftSelection(hHatj,HHatPrimejj,hHatj2,HHatPrimej2)

[~,Nj] = size(HHatPrimejj);

anglee1 = zeros(Nj,1);

A = (HHatPrimejj'*HHatPrimejj);

[U,Lambda] = eig(A);

[maxx,indexx] = max(diag(Lambda));

if sum(abs(hHatj))==0

    anglee1 = exp(1i*angle(U(:,indexx)));


else

    temporr = zeros(Nj,1);
    for d = 1:Nj
        temporr(d) = U(:,d)'*HHatPrimejj'*hHatj;
    end

    gammaLow = maxx;
    gammaUpp = inf;
    gamma = gammaLow+1;

    diff = 100;
    while diff>0.00001
        summ = sum((abs(temporr).^2)./(gamma-diag(Lambda)).^2);
        if summ>Nj
            gammaLow = gamma;
        else
            gammaUpp = gamma;
        end

        if gammaUpp<inf
            gamma = (gammaLow+gammaUpp)/2;
        else
            gamma = 2*gammaLow;
        end

        diff = abs(gammaLow-gammaUpp);
    end


    for d = 1:Nj
        anglee1 = anglee1 + temporr(d)/(gamma-Lambda(d,d))*U(:,d);
    end
    anglee1 = exp(1i*angle(anglee1));




end


%%%%%%%%%%%%%

[~,Nj] = size(HHatPrimejj);

anglee2 = zeros(Nj,1);
maxEigvalue1 = max(eig(HHatPrimejj'*HHatPrimejj));
maxEigvalue2 = max(eig(HHatPrimej2'*HHatPrimej2));

epsi = maxEigvalue1/maxEigvalue2*0.2;
A = (HHatPrimejj'*HHatPrimejj)-epsi*(HHatPrimej2'*HHatPrimej2);

[U,Lambda] = eig(A);

[maxx,indexx] = max(diag(Lambda));

if sum(abs((HHatPrimejj'*hHatj-epsi*HHatPrimej2'*hHatj2)))==0

    anglee2 = exp(1i*angle(U(:,indexx)));

else

    temporr = zeros(Nj,1);
    for d = 1:Nj
        temporr(d) = U(:,d)'*(HHatPrimejj'*hHatj-epsi*HHatPrimej2'*hHatj2);
    end

    gammaLow = maxx;
    gammaUpp = inf;
    gamma = gammaLow+1;

    diff = 100;
    while diff>0.00001
        summ = sum((abs(temporr).^2)./(gamma-diag(Lambda)).^2);
        if summ>Nj
            gammaLow = gamma;
        else
            gammaUpp = gamma;
        end

        if gammaUpp<inf
            gamma = (gammaLow+gammaUpp)/2;
        else
            gamma = 2*gammaLow;
        end

        diff = abs(gammaLow-gammaUpp);
    end


    for d = 1:Nj
        anglee2 = anglee2 + temporr(d)/(gamma-Lambda(d,d))*U(:,d);
    end
    anglee2 = exp(1i*angle(anglee2));




end

%%%%%%%%%%%%%%%%%555

anglee3 = exp(1i*angle(HHatPrimejj'*hHatj));

