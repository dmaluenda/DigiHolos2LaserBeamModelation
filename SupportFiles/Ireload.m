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



mapa1=load(['valors_p' num2str(1) '.txt']);
ph1=mapa1(:,3);A1_0=mapa1(:,2);

mapa2=load(['valors_p' num2str(2) '.txt']);
ph2=mapa2(:,3);A2_0=mapa2(:,2);

T1=sqrt(I1/max(I1));
T2=sqrt(I2/max(I2));

curve1=[(1:256)' T1 ph1];
curve2=[(1:256)' T2 ph2];

fid2 = fopen(['valors_p' num2str(1) '.txt'],'wt');
fprintf(fid2,'%3.0f %6.4f% 6.4f\n',curve1');
fclose(fid2);

fid2 = fopen(['valors_p' num2str(2) '.txt'],'wt');
fprintf(fid2,'%3.0f %6.4f% 6.4f\n',curve2');
fclose(fid2);