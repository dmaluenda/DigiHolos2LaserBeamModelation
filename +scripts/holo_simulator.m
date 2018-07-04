%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubó - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Retuns the ElectricField distribution, as an images for each complex
% component, at the conjugate plane of C_x and C_y, whose are determined by
% the pointres gl1 and gl2 on SLMs
% [C_img,C_img] = f([img,[img])  -->  [a+bi,c+di] = f([1 256],[1 256])

function[E_x,E_y]=holo_simulator(slm1,slm2)

	slm1 = uint8(slm1*255)+1;
	slm2 = uint8(slm2*255)+1;
	
	N = size(slm1);
	C_SLM1 = zeros(N);
	C_SLM2 = zeros(N);
	% figure
	% imagesc(slm1)
	% figure
	% imagesc(slm2)

	values1 = load([pwd '\SupportFiles\Valors_p1.txt']);
	mapa1   = values1(:,1);
	T_slm1  = values1(:,2);
	ph_slm1 = mod(values1(:,3),2*pi);

	values2 = load([pwd '\SupportFiles\Valors_p2.txt']);
	mapa2   = values2(:,1);
	T_slm2  = values2(:,2);
	ph_slm2 = mod(values2(:,3),2*pi);


	C_SLM1(:) = T_slm1(slm1(:)) .* exp(1i* ph_slm1(slm1(:)) );
	C_SLM2(:) = T_slm2(slm2(:)) .* exp(1i* ph_slm2(slm2(:)) );

	% figure
	% imagesc(angle(C_SLM1))
	% figure
	% imagesc(abs(C_SLM1))
	% figure
	% imagesc(angle(C_SLM2))
	% figure
	% imagesc(abs(C_SLM2))

	FC_x = fftshift(fft2(C_SLM1));
	FC_y = fftshift(fft2(C_SLM2));

	spatial_filter = zeros(N);
	spatial_filter(3*N(1)/8 : 5*N(1)/8-1 , 3*N(2)/8 : 5*N(2)/8-1) = ones(N/4);
	
	FC_x = FC_x.*spatial_filter;
	FC_y = FC_y.*spatial_filter;

	E_x = ifft2(FC_x);
	E_y = ifft2(FC_y);

	% figure
	% imagesc(abs(E_x))
	% figure
	% imagesc(abs(E_y))
	% figure
	% imagesc(angle(E_x))
	% figure
	% imagesc(angle(E_y))
	