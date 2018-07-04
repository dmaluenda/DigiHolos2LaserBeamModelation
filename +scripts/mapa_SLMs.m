%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubó - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[SLM1,SLM2]=mapa_SLMs(phase1,phase2)

SLM1 = phase1; %to initiate
data = load('curve_SLM1.txt');
mapa = data(:,2); %[0 255]
phase1(:) = floor(mod(phase1(:),2*pi)*1000+1); 
SLM1(:) = mapa(phase1(:));
SLM1 = SLM1/255;

SLM2 = phase2; %to initiate
data = load('curve_SLM2.txt');
mapa = data(:,2);
phase2(:) = floor(mod(phase2(:),2*pi)*1000+1);
SLM2(:) = mapa(phase2(:));
SLM2 = SLM2/255;