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

% Fills a hologram with a chees distribution values (Arizon procedure)
%  in order to get a certain Transmision and Phase shift

function[SLM1] = mapa_Holo3(Trans,Phase,N,SLM_number)
	if SLM_number==1
	    ph0 = 0*pi/180;
	    Phase = Phase*pi/180+ph0;
	else
	    ph0 = 0*pi/180;
	    Phase = Phase*pi/180+ph0;
	end

	SLM1 = zeros(N);

	Cpath = cd;
	cd ..
	PATH = cd;
	cd(Cpath)
	data1 = load([PATH '/ComplexValues_SLM' num2str(SLM_number) '.txt']);

	T_SLM1  = data1(:,1);
	ph_SLM1 = mod(data1(:,2),2*pi);
	Mapa1_1 = data1(:,3);
	Mapa2_1 = data1(:,4);

	% Desirable value
	C = Trans*exp(1i*Phase);

	[~,p]=min(abs(C-T_SLM1.*exp(1i*ph_SLM1)));

	SLM1(1:2:N(1), 1:2:N(2)) = Mapa1_1(p); % M_L (top-left)
	SLM1(2:2:N(1), 2:2:N(2)) = Mapa1_1(p); % M_L (botom-right)
	SLM1(2:2:N(1), 1:2:N(2)) = Mapa2_1(p); % M_L (top-right)
	SLM1(1:2:N(1), 2:2:N(2)) = Mapa2_1(p); % M_L (botom-left)

	SLM1=(SLM1-1)/255;

	% imshow(SLM1')
	% figure
	% imshow(SLM2')
	% figure
	% imagesc(m1)
	% figure
	% imagesc(m2)
