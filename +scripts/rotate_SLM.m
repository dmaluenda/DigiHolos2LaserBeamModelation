%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare the images to be displayed on the SLMs
%   ---> SLM1 is rotated and flipped <---
%   --->    SLM2 is just flipped     <---
function[OUT_1,OUT_2]=rotate_SLM(IN_1,IN_2)

	N = size(IN_1);
	n = min(N);
	m = max(N);

	out_1 = IN_1';
	OUT_1 = zeros(N);
	OUT_1(m/2-n/2:m/2+n/2-1,1:n) = out_1(n:-1:1, m/2+n/2:-1:m/2-n/2+1);

	OUT_2(N(1):-1:1,N(2):-1:1)=IN_2(1:N(1),1:N(2));
