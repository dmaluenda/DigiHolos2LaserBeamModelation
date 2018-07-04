%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubó - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% David Maluenda (2013)
%
% Entren les componets complexes (x,y) del camp. Està pensat per 16x16.
% Selecciona els 8x8 de la part central.
% Np és el nombre de punts en que dibuixa cada el·lipse.


function[] = delipse(S0,S1,S2,S3)

    Np  = 30;
    N   = 9;
    Nim = max(size(S0));
    s0  = zeros(N,N);
    s1 = s0 ; s2 = s0 ; s3 = s0 ;

    for i=1:N
        for j=1:N
            s0(i,j) = mean(mean(S0((i-1)*Nim/N+1:i*Nim/N, (j-1)*Nim/N+1:j*Nim/N)));
            s1(i,j) = mean(mean(S1((i-1)*Nim/N+1:i*Nim/N, (j-1)*Nim/N+1:j*Nim/N)));
            s2(i,j) = mean(mean(S2((i-1)*Nim/N+1:i*Nim/N, (j-1)*Nim/N+1:j*Nim/N)));
            s3(i,j) = mean(mean(S3((i-1)*Nim/N+1:i*Nim/N, (j-1)*Nim/N+1:j*Nim/N)));
        end
    end
         
    % L=sqrt(s1.^2+s2.^2);
    A0 = (s0+s1); %sqrt((s0+L)/2);
    B0 = (s0-s1); %sqrt((s0-L)/2);
    theta = atan2(s2, s1)/2;

    M = max(max(max(A0)), max(max(B0)));
    A = A0/max(max(M))/2.25;
    B = B0/max(max(M))/2.25;

    Cph = s2./(sqrt(s0.^2-s1.^2-s2.^2));
    Sph = s3./sqrt(s0.^2-s1.^2-s2.^2);
    ph  = acos(Cph);
    phs = asin(Sph);

    figure;
    hold on;
    for ii=1:N;
        for jj=1:N;
            for t1=1:Np;
                t = t1*2*pi/Np;

                x = ii + A(ii,jj).*cos(t);
                y = jj + B(ii,jj).*cos(t+ph(ii,jj));

                plot(x,y,'.k','MarkerSize',8);
            end
        end
    end
    axis([0 N+1 0 N+1]);
    hold off;
