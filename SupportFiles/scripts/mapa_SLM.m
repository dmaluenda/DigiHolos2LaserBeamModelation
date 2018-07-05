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

function[SLM]=mapa_SLM(phase,SLM_number)
k=SLM_number;
SLM=phase; %to initiate
data=load(['curve_SLM' num2str(k) '.txt']);
mapa=data(:,2);
phase(:)=floor(mod(phase(:),2*pi)*1000+1);
SLM(:)=mapa(phase(:));

SLM=SLM/255; %1-?????