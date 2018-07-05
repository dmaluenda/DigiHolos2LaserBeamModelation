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

%Returns the pointer which points to the earliest complex value on the SLM
%C is a SCALAR anf T_SLM and ph_SLM are ARRAYs of the possibles candidates.
function[m,p] = min_dist(C, T_slm, ph_slm, M)

	for k=1:M
	        dist = (abs(C-T_slm(k).*exp(1i*ph_slm(k))));
	end
	[m,p] = min(dist);