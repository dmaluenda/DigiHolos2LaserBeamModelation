%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a serial of holograms to encode T=0:Tstep:1 and Ph=0:PhStep:PhFinal
%  in a half of the image remaining the other side at T=1 and Ph=0

function[SLM1]=mapa_Holo2(Trans,Phase,N,SLM_number)
	if SLM_number==1
	    X0 = [553 640];
	    Y0 = [270 344];
	    X1 = [385 468];
	    Y1 = [270 344];
	    ph0 = 0*pi/180;
	    Phase = Phase*pi/180+ph0;
	else
	    D = (N(1)-N(2))/2;
	    X1 = [553-D 640-D];
	    Y1 = [270+D 350+D];
	    X0 = [385-D 468-D];
	    Y0 = [270+D 350+D];
	    ph0 = 0*pi/180;
	    Phase = Phase*pi/180+ph0;
	end

	X = [1 N(2)];
	Y = [round(N(1)/2) N(1)];

	SLM1 = zeros(N);

	Amp_max1 = 1;
	Amp_max2 = 1;
	A_max = min([Amp_max1 Amp_max2]);

	Cpath = cd;
	cd ..
	PATH = cd;
	cd(Cpath)
	data1 = load([PATH '/ComplexValues_SLM' num2str(SLM_number) '.txt']);

	T_SLM1=data1(:,1);
	ph_SLM1=mod(data1(:,2),2*pi);
	Mapa1_1=data1(:,3);
	Mapa2_1=data1(:,4);

	Trans=Trans*A_max;

	C=Trans*exp(1i*Phase);

	% REFERENCE SIDE : T=1 ; Ph=0
	[~,p]=min(abs(1*exp(1i*ph0+pi/2)-T_SLM1.*exp(1i*ph_SLM1)));
	SLM1(1:2:N(1), 1:2:N(2)) = Mapa1_1(p); % M_L (top-left)
	SLM1(2:2:N(1), 2:2:N(2)) = Mapa1_1(p); % M_L (botom-right)
	SLM1(2:2:N(1), 1:2:N(2)) = Mapa2_1(p); % M_L (top-right)
	SLM1(1:2:N(1), 2:2:N(2)) = Mapa2_1(p); % M_L (botom-left)

	% TARGET SIDE
	[~,px]=min(abs(C-T_SLM1.*exp(1i*ph_SLM1)));
	SLM1(Y(1)  :2:Y(2), X(1)  :2:X(2)) = Mapa1_1(px); % M_L (top-left)
	SLM1(Y(1)+1:2:Y(2), X(1)+1:2:X(2)) = Mapa1_1(px); % M_L (botom-right)
	SLM1(Y(1)+1:2:Y(2), X(1)  :2:X(2)) = Mapa2_1(px); % M_L (top-right)
	SLM1(Y(1)  :2:Y(2), X(1)+1:2:X(2)) = Mapa2_1(px); % M_L (botom-left)

	SLM1 = (SLM1-1)/255;

	% imshow(SLM1')
	% figure
	% imshow(SLM2')
	% figure
	% imagesc(m1)
	% figure
	% imagesc(m2)
