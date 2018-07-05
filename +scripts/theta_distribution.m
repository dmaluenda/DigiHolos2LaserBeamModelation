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

clear all
%R = 70;
%rho  = 0:0.005:1;
theta = 0:0.01:2*pi;

%A_0 = sqrt(rho.^2.*exp(-rho));
%A_0 = A_0/max(max(A_0));
k = 1/2;
l = 1/2;
theta_amp = 0;
theta_ph = 0;

E_x = abs(cos(theta*k + theta_amp));
E_y = abs(sin(theta*k + theta_amp));
% Ph_x = angle(E_x);
% Ph_y = angle(E_y);
% E_x = A_0.*abs(E_x);
% E_y = A_0.*abs(E_y);

Ph_x = -l*theta+pi/2;

px = find(l*theta>pi/4);
Ph_x(px) = l*theta(px);

px = find(l*theta>3*pi/4);
Ph_x(px) = -l*theta(px)-pi/2;

px = find(l*theta>5*pi/4);
Ph_x(px) = l*theta(px)-pi;

px = find(l*theta>7*pi/4);
Ph_x(px) = -l*theta(px)+pi/2;

Ph_x = mod(Ph_x, 2*pi);



Ph_y = l*theta+pi/2;

py = find(l*theta>pi/4);
Ph_y(py) = -l*theta(py)-pi;

py = find(l*theta>3*pi/4);
Ph_y(py) = l*theta(py)-pi/2;

py = find(l*theta>5*pi/4);
Ph_y(py) = -l*theta(py);

py = find(l*theta>7*pi/4);
Ph_y(py) = l*theta(py)+pi/2;

Ph_y = mod(Ph_y, 2*pi);
