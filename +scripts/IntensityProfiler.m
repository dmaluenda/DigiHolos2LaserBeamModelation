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

function[im_out] = IntensityProfiler(im_in,Q)

	Im = fftshift(fft2(im_in));

	[x,y] = meshgrid(-1024/2:1024/2-1,-768/2:768/2-1);

	circ = (x.^2+y.^2<Q^2).*1;
	IM   = Im.*circ;

	figure;
	imshow(log(abs(IM)/max(max(IM))))

	im_out=ifft2(fftshift(IM));