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


A0   =     [1  0;  0  0];
A90  =     [0  0;  0  1];
A45  = 0.5*[1  1;  1  1];
A135 = 0.5*[1 -1; -1  1];
ACD  = 0.5*[1  1; 1i 1i];
ACL  = 0.5*[1 -1;-1i 1i];

NN = size(Ex);
st = 10;

E(:,:,1) = Ex(1:st:NN(1), 1:st:NN(2));
E(:,:,2) = Ey(1:st:NN(1), 1:st:NN(2));

N = size(E);

E45  = zeros(N);
E135 = E45 ; E0 = E45 ; E90 = E45 ; ECD = E45 ; ECL = E45;

I45  = zeros(N(1),N(2));
I135 = I45;

JV=[0;0];
for ii=1:N(1)
    tic
    for jj=1:N(2)
        JV(:) = E(ii,jj,:);
        E0(ii,jj,:)   = A0*JV;
        E45(ii,jj,:)  = A45*JV;
        E90(ii,jj,:)  = A90*JV;
        E135(ii,jj,:) = A135*JV;
        ECD(ii,jj,:)  = ACD*JV;
        ECL(ii,jj,:)  = ACL*JV;
    end 
    disp([toc 's: ' num2str(ii) '/' num2str(N(1))])
end

I0(:,:)   = abs(  E0(:,:,1).^2 +   E0(:,:,2).^2);
I45(:,:)  = abs( E45(:,:,1).^2 +  E45(:,:,2).^2);
I90(:,:)  = abs( E90(:,:,1).^2 +  E90(:,:,2).^2);
I135(:,:) = abs(E135(:,:,1).^2 + E135(:,:,2).^2);
ICD(:,:)  = abs( ECD(:,:,1).^2 +  ECD(:,:,2).^2);
ICL(:,:)  = abs( ECL(:,:,1).^2 +  ECL(:,:,2).^2);

figure; imagesc(I0');  title 'I0'
figure; imagesc(I45'); title 'I45'
figure; imagesc(I90'); title 'I90'
figure; imagesc(I135');title 'I135'
figure; imagesc(ICD'); title 'ICD'
figure; imagesc(ICL'); title 'ICL'
