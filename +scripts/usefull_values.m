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


U1  = ( T_slm1<=A_max & ph_slm1>=30*pi/180 & ph_slm1<=30*pi/180+pi );
pU1 = find(U1)';

T_SLM1  = T_slm1(pU1);
ph_SLM1 = ph_slm1(pU1);
Mapa1_1 = mapa1_1(pU1);
Mapa2_1 = mapa2_1(pU1);

M1 = max(size(T_SLM1));



U2  = ( T_slm2<=A_max & ph_slm2>=30*pi/180 & ph_slm2<=30*pi/180+pi );
pU2 = find(U2)';

T_SLM2  = T_slm2(pU2);
ph_SLM2 = ph_slm2(pU2);
Mapa1_2 = mapa1_2(pU2);
Mapa2_2 = mapa2_2(pU2);
M2 = max(size(T_SLM2));
