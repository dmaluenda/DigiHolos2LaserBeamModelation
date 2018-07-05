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



% xcorr demo
%
% The signals
t = [1:451];
P = 30*sqrt(2);
s1 = sin(2*pi*t/P);
s2 = sin(2*pi*t/P+90*pi/180); % s1 lags s2 by 0.35s
figure
subplot(2,1,1);
plot(t,s1,'r',t,s2,'b');
grid
title('signals')
%
% Now cross-correlate the two signals
%
x = xcorr(s1,s2,'coeff');
tx = [-450:450];
subplot(2,1,2)
plot(tx,x)
grid
%
% Determine the lag
%
[mx,ix] = max(x);
lag = tx(ix)*360/P;
hold on
tm = [lag,lag];
mm = [-1,1];
plot(tm,mm,'k')
hold off
% Note that the lag is only as close as the time resolution.
% i.e. actual = -0.35, calculated = -0.34
S = sprintf('Lag = %5.2f',lag);
title(S)