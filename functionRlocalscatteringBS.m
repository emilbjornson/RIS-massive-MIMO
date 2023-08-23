function R = functionRlocalscatteringBS(N,varphi,theta,ASD_varphi,ASD_theta,antennaSpacing)
%Prepare to compute the first row of the matrix
firstRow  = zeros(N,1);
firstRow(1) = 1;

%% Go through all the columns of the first row
for column = 2:N

    %Distance from the first antenna
    distanceY = antennaSpacing*(1-column);

    A = exp(1i*2*pi*(distanceY*sin(varphi)*cos(theta)));
    B = 2*pi*distanceY*cos(varphi)*cos(theta);
    C = -2*pi*distanceY*cos(varphi)*sin(theta);
    D = -2*pi*distanceY*sin(varphi)*sin(theta);

    sigma2tilde = ASD_varphi^2/(1+C^2*ASD_varphi^2*ASD_theta^2);
    firstRow(column) = A*sqrt(sigma2tilde)/ASD_varphi*...
        exp(D^2*ASD_theta^2*(C^2*ASD_theta^2*sigma2tilde -1)/2)*...
        exp(-1i*B*C*D*ASD_theta^2*sigma2tilde)*...
        exp(-B^2*sigma2tilde/2);


end

%Compute the spatial correlation matrix by utilizing the Toeplitz structure
%and scale it properly
R = toeplitz(firstRow);
R = R*(N/trace(R));
