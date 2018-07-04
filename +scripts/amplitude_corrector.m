%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubó - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns the amplitudes to achieve a flat profile
function[Trans1,Trans2] = amplitude_corrector(Amp1, Amp2)

	Amp1 = CircularPupil(Amp1,5);
	Amp2 = CircularPupil(Amp2,5);

	A_0 = min( [max(max(Amp1)) , max(max(Amp2))] )/exp(1);
	
	N1 = size(Amp1);
	N2 = size(Amp2);

	P1 = Amp1>A_0;
	P2 = Amp2>A_0;
	p1 = find(P1)';
	p2 = find(P2)';

	Trans1 = ones(N1);
	Trans2 = ones(N2);

	Trans1(p1) = A_0./Amp1(p1);
	Trans2(p2) = A_0./Amp2(p2);

	% figure
	% imagesc(Trans1)
	% figure
	% imagesc(Trans2)