%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%R = 70;
%rho = 0:0.005:1;
theta = 0:0.01:2*pi;    
k = 1/2;
l = 1;
theta_amp = 0;
theta_ph = 0;

Ph = mod(l*theta+theta_ph,2*pi);
Ph_x = zeros(size(theta));

px = find(Ph<pi);
Ph_x(px) = Ph(px);

Ph_y = zeros(size(theta));
py = find(Ph>=pi);

Ph_y(py)=mod(-Ph(py),2*pi);