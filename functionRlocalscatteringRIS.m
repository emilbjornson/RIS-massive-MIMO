function R = functionRlocalscatteringRIS(Nh,Nv,varphi,theta,ASD_varphi,ASD_theta,antennaSpacing)

N = Nh*Nv;
posY = vec(repmat(0:Nh-1,Nv,1)*antennaSpacing);
posZ = vec(repmat((0:Nv-1).',1,Nh)*antennaSpacing);

%Prepare the matrix
R = zeros(N,N);

%% Go through all the columns of the first row
for row = 1:N
    for column = row+1:N



        distanceY = posY(row)-posY(column);
        distanceZ = posZ(row)-posZ(column);

        A = exp(1i*2*pi*(distanceY*sin(varphi)*cos(theta)+distanceZ*sin(theta)));
        B = 2*pi*distanceY*cos(varphi)*cos(theta);
        C = -2*pi*distanceY*cos(varphi)*sin(theta);
        D = -2*pi*distanceY*sin(varphi)*sin(theta)+2*pi*distanceZ*cos(theta);

        sigma2tilde = ASD_varphi^2/(1+C^2*ASD_varphi^2*ASD_theta^2);
        R(row,column) = A*sqrt(sigma2tilde)/ASD_varphi*...
            exp(D^2*ASD_theta^2*(C^2*ASD_theta^2*sigma2tilde -1)/2)*...
            exp(-1i*B*C*D*ASD_theta^2*sigma2tilde)*...
            exp(-B^2*sigma2tilde/2);

    end
end

R = (R+R')+eye(N);
