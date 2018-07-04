%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all


data = load('Holos/T1_float.txt');

Amax = 0.57;
A = data(:,1)/100*Amax;
P = data(:,2).*10.^(-round(data(:,3)));
PB= 2.548*10^(-3);
P = (P)/PB;

figure;
Atest=linspace(-Amax,Amax,100);
plot(A,P,'xk',Atest,Atest.^2,'-k');

figure;
title 'Relative Error'
plot(A,abs(sqrt(P)-abs(A))./sqrt(P))

figure;
title 'Absolute Error'
plot(A,abs(sqrt(P)-abs(A)))

absErr = mean(sqrt(P)-abs(A)) % sqrt(sum(abs2((sqrt(P)-abs(A))./sqrt(P))))


figure;
Ptest=linspace(0,Amax^2,5);
plot(A( 1: 50).^2,P( 1: 50),'xk', ...
     A(52:101).^2,P(52:101),'+k', ...
     A(  51  ).^2,P(  51  ),'ok', ...
         Ptest   ,  Ptest  ,'-k')  ;

figure;
plot(A,sqrt(P).*sign(A))


Contrast = ((max(P))-(min(P)))/((max(P))+(min(P)))
ExtintionRatio = (min(P)/max(P))


