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


function[im_out] = CircularPupil(im_in,R)
	N = size(im_in);

	Im = fftshift(fft2(im_in));

	[x,y] = meshgrid(-N(2)/2:N(2)/2-1,-N(1)/2:N(1)/2-1);

	circ = (x.^2+y.^2<R^2).*1;
	IM = Im.*circ;
	
	im_out = ifft2(fftshift(IM));
	im_out = abs(im_out)/max(max(im_out));