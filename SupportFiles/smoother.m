%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Smooths the curve phi, can be used several times as smoothed as you want
in=phi2;  %Ak phik_45 phik_135
% out=in;
% for i=2:255
%     out(i)=mean([in(i-1) in(i-1) in(i) in(i+1) in(i+1)]);
% end

out=smooth(in);

plot(out)

phi2=out;