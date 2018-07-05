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

data1=load(['ComplexValues_SLM' num2str(SLM_number) '.txt']);
T_slm1=data1(:,1);
ph_slm1=mod(data1(:,2),2*pi);
mapa1_1=data1(:,3);
mapa2_1=data1(:,4);


U1=(T_slm1<=A_max&ph_slm1>=20*pi/180&ph_slm1<=40*pi/180+pi);
pU1=find(U1)';
T_SLM1=T_slm1(pU1);
ph_SLM1=ph_slm1(pU1);
Mapa1_1=mapa1_1(pU1);
Mapa2_1=mapa2_1(pU1);
M1=max(size(T_SLM1));
