%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Complex values. (2012)
%David Maluenda Niubó - Applied Physics and Optics (UB)

clear variables; close all
%clear all;close all
%k=2;
date=clock;
N1=256;N2=175;
A_max=1;%just for ploting
%usefull values and SemiCercle
phi1_0=48*pi/180;  %    <-Rotate
phi2_0=60*pi/180;  %    <-Rotate
A1_maxCM = 0.30;   %  <-Trim1 in Complex Modulation
A2_maxCM = 0.25;   %  <-Trim2 in Complex Modulation
A1_maxAM = 0.65;   %  <-Trim1 in Amplitude Modulation
A2_maxAM = 0.5;   %  <-Trim2 in Amplitude Modulation
errA = 0.1;
errP = 10/360;
% A_max = A1_max;
mapa1=load(['valors_p' num2str(1) '.txt']);
ph1=mapa1(1:N1,3);A1_0=mapa1(1:N1,2);
mapa2=load(['valors_p' num2str(2) '.txt']);
ph2=mapa2(1:N2,3);A2_0=mapa2(1:N2,2);

A1=A1_0.*exp(1i*ph1);
A2=A2_0.*exp(1i*ph2);

count = 1; count1 = 1;
for i=1:N1
    for j=i:N1
        aux=(A1(i)+A1(j))/2;
        B1(count)=aux;
        aux1_i(count)=i;
        aux1_j(count)=j;
        if abs(aux)<=A1_maxCM
            pC1(count1) = count;
            count1 = count1 + 1; 
        end
            
    auxANG = 180 + angle(A1(i)) - angle(A1(j));    
    errAmpA = ( abs2( abs(A1(i)) - abs(A1(j))*cos( auxANG ) ) + ...
                abs2( abs(A1(j)) - abs(A1(i))*cos( auxANG ) ) ) * errA^2;
    errAmpP = abs2( 2*abs(A1(i))*abs(A1(j)) * sin( auxANG ) ) * errP^2;
    factorERR = ( abs2(A1(j)) - abs2(A1(i)) - abs2(2*aux) )/(2*abs(A1(i))*abs(2*aux));
    errPh   = errP^2 + abs2( 2*abs(A1(j))/sqrt(1+factorERR^2))*errA^2;    
    errB1(count) = sqrt(errAmpA+errAmpP)/4/abs(aux)*exp(1i*sqrt(errPh));
            
       count=count+1;
    end
end
count = 1; count2 = 1;
for i=1:N2
    for j=i:N2
        aux=(A2(i)+A2(j))/2;
        B2(count)=aux;
        aux2_i(count)=i;
        aux2_j(count)=j;
        if abs(aux)<=A2_maxCM
            pC2(count2) = count;
            count2 = count2 + 1;
        end
        count = count + 1;
    end
end

% Amplitude only modulation
N_At   = 300;
Ats    = linspace(-1,1,N_At);
count1 = 1 ; count2 = 1 ;
aux1   = 0 ; aux2   = 0 ;
for iA = 1:N_At

    At1 = Ats(iA)*exp(1i*phi1_0)*A1_maxAM;
    At2 = Ats(iA)*exp(1i*phi2_0)*A2_maxAM;

    [~,pA] = min(abs( At1 - B1 ));
    [~,pB] = min(abs( At2 - B2 ));
    
    if pA~=aux1
        Am1(count1) = B1(pA);
        pA1(count1) = pA;
        aux1   = pA;
        count1 = count1+1;
    end
    
    if pB~=aux2
        Am2(count2) = B2(pB);
        pA2(count2) = pB;
        aux2   = pB;
        count2 = count2+1;
    end
end
% Contrast1 = (max(abs2(Am1))-min(abs2(Am1)))/(max(abs2(Am1))+min(abs2(Am1)));
% ExtintionRatio1 = min(abs2(Am1))/max(abs2(Am1));

h=figure;
plotPMsigma((1:256)',abs(A1),ones(256,1)*0.02);
axis([1 256 0 1.05])

h=figure;
plotPMsigma((1:256)',mod(angle(A1)*180/pi,360),ones(256,1)*18);
axis([1 256 0 360])

A1_C = A1;
% ploting SLM1
h = figure;
polar( angle(A1_C) , abs(A1_C) , '-r' ) %Real curve
hold on
polar( angle(B1) , abs(B1) , '+g' ) %Possibles values
polar( angle(B1(pC1)) , abs(B1(pC1)) , '+b' ) %Complex Modulation
hp=polar(angle(B1(pA1)),abs(B1(pA1)),'ok');set(hp,'MarkerSize',10); % Amplitude Mod
title 'SLM1'
legend 'SLM curve' 'Possible coding values' 'Complex Modulation' 'Amplitude only modulation'
hold off
cd plots
print(h,'-depsc',[num2str(date(3),'%02.0f') num2str(date(2),'%02.0f') ...
    num2str(date(1)-2000) '_polar_SLM1.eps']);
cd ..

% % Creating data for SLM1
% Complex Modulation
phi1    = mod(angle(B1)-phi1_0,2*pi);
dataCM1 = [ abs(B1(pC1))/A1_maxCM ; phi1(pC1) ; aux1_i(pC1) ; aux1_j(pC1) ];
fid     = fopen( 'ComplexValues_SLM1.txt' , 'wt' );
fprintf( fid , '%10.20f  %10.20f  %3.0f  %3.0f\n' , dataCM1 );
fclose(fid);
% Amplitude Modulation
dataAM1 = [ abs(B1(pA1))/A1_maxAM ; phi1(pA1) ; aux1_i(pA1) ; aux1_j(pA1) ];
fid     = fopen( 'AmplitudeValues_SLM1.txt' , 'wt' );
fprintf( fid , '%10.20f  %10.20f  %3.0f  %3.0f\n' , dataAM1 );
fclose(fid);


%ploting SLM2
h = figure;
polar( angle(A2) , abs(A2) , '-r' ) %Real curve
hold on
polar( angle(B2) , abs(B2) , '+g' ) %Possibles values
polar( angle(B2(pC2)) , abs(B2(pC2)) , '+b' ) %Complex Modulation
hp=polar(angle(B2(pA2)),abs(B2(pA2)),'ok');set(hp,'MarkerSize',10); % Amplitude Mod
title 'SLM2'
legend 'SLM curve' 'Possible coding values' 'Complex Modulation' 'Amplitude only modulation'
hold off
cd plots
print(h,'-depsc',[num2str(date(3),'%02.0f') num2str(date(2),'%02.0f') ...
    num2str(date(1)-2000) '_polar_SLM2.eps']);
cd ..


% % % Creating data for SLM2
% Complex Modulation
phi2    = mod(angle(B2)-phi2_0,2*pi);
dataCM2 = [ abs(B2(pC2))/A2_maxCM ; phi2(pC2) ; aux2_i(pC2) ; aux2_j(pC2) ];
fid     = fopen( 'ComplexValues_SLM2.txt' , 'wt' );
fprintf( fid , '%10.20f  %10.20f  %3.0f  %3.0f\n' , dataCM2);
fclose(fid);
% Amplitude Modulation
dataAM2 = [ abs(B2(pA2))/A2_maxAM ; phi2(pA2) ; aux2_i(pA2) ; aux2_j(pA2) ];
fid     = fopen( 'AmplitudeValues_SLM2.txt' , 'wt' );
fprintf( fid , '%10.20f  %10.20f  %3.0f  %3.0f\n' , dataAM2 );
fclose(fid);




