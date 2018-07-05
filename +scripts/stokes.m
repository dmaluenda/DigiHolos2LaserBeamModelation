%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%                  David Maluenda Niubo - dmaluendn@gmail.com            
%
%      CC: by, NC, SA                                    2012-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[s0,s1,s2,s3,DOP]=stokes(Ex,Ey)

    A0   =     [1  0;  0  0];
    A90  =     [0  0;  0  1];
    A45  = 0.5*[1  1;  1  1];
    A135 = 0.5*[1 -1; -1  1];
    ACD  = 0.5*[1  1; 1i 1i];
    ACL  = 0.5*[1 -1;-1i 1i];

    NN = size(Ex);
    st=4;
    E(:,:,1) = Ex(1:st:NN(1), 1:st:NN(2));
    E(:,:,2) = Ey(1:st:NN(1), 1:st:NN(2));

    N = size(E); 
    I0a = zeros(N);
    I135a = I0a ; I45a = I0a ; I90a = I45a ; ICDa = I0a ; ICLa = I0a ;
    
    JV = [0;0];
    for ii=1:N(1)
        tic
        for jj=1:N(2)
            JV(1) = E(ii,jj,1);
            JV(2) = E(ii,jj,2);
            E0 = A0*JV;
            I0a(ii,jj) = abs(E0(1))^2 + abs(E0(2))^2;
            
            E45=A45*JV;
            I45a(ii,jj) = abs(E45(1))^2 + abs(E45(2))^2;
            
            E90=A90*JV;
            I90a(ii,jj) = abs(E90(1))^2 + abs(E90(2))^2;
            
            E135=A135*JV;
            I135a(ii,jj) = abs(E135(1))^2 + abs(E135(2))^2;
            
            ECD=ACD*JV;
            ICDa(ii,jj)= abs(ECD(1))^2 + abs(ECD(2))^2;
            
            ECL=ACL*JV;
            ICLa(ii,jj)= abs(ECL(1))^2 + abs(ECL(2))^2;
        end 
        disp([toc 's: ' num2str(ii) '/' num2str(N(1))])
    end

    I0  = I0a(:,:,1);
    I45 = I45a(:,:,1);
    I90 = I90a(:,:,1);
    I135= I135a(:,:,1);
    ICD = ICDa(:,:,1);
    ICL = ICLa(:,:,1);

    figure;
    subplot(2,3,1); imagesc( I0' ,[0 1]); title 'I0'
    subplot(2,3,2); imagesc( I45',[0 1]); title 'I45'
    subplot(2,3,3); imagesc( I90',[0 1]); title 'I90'
    subplot(2,3,4); imagesc(I135',[0 1]); title 'I135'
    subplot(2,3,5); imagesc( ICD',[0 1]); title 'ICD'
    subplot(2,3,6); imagesc( ICL',[0 1]); title 'ICL'

    s0 = I0+I90;
    s1 = I0-I90;
    s2 = I45-I135;
    s3 = ICD-ICL;
    DOP=(s1.^2+s2.^2+s3.^2)./s0.^2;
