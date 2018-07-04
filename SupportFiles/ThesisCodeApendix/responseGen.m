%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



load('ExperimentalDATA.mat');

header=['gl' 'T1' 'ph1' 'T2' 'ph2'];

data=[(1:256)' T1' ph1' T2' ph2'];

delete('SLMresponse.txt');
fid = fopen('SLMresponse.txt', 'at');

fprintf(fid, '%s\n'   , ' gray  |______SLM1_________|______SLM2_______|');

fprintf(fid, '%s\n\n' , '_level_|__T1___|__ph1______|__T2___|__ph2____|');

fprintf(fid, '%.0f\t\t%.4f\t%.4f\t\t%.4f\t%.4f\n', data.' );

fclose(fid);