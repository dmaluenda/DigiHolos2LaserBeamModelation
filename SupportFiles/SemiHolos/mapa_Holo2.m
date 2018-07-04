%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[SLM1]=mapa_Holo2(Trans,Phase,N,SLM_number)
if SLM_number==1
    X0=[553 640];Y0=[270 344];
    X1=[385 468];Y1=[270 344];
    ph0=0*pi/180;
    Phase=Phase*pi/180+ph0;
else
    D=(N(1)-N(2))/2;
    X1=[553-D 640-D];Y1=[270+D 350+D];
    X0=[385-D 468-D];Y0=[270+D 350+D];
    ph0=0*pi/180;
    Phase=Phase*pi/180+ph0;
end
X=[1 N(2)];Y=[round(N(1)/2) N(1)];

SLM1=zeros(N);

Amp_max1=1;
Amp_max2=1;
A_max=min([Amp_max1 Amp_max2]);

Cpath=cd;
cd ..
path=cd;
cd(Cpath)
data1=load([path '/ComplexValues_SLM' num2str(SLM_number) '.txt']);

T_SLM1=data1(:,1);
ph_SLM1=mod(data1(:,2),2*pi);
Mapa1_1=data1(:,3);
Mapa2_1=data1(:,4);


%usefull_values; %load the usefull values in T_SLMi, ph_SLMi and Mapaj_i

Trans=Trans*A_max;

C=Trans*exp(1i*Phase);

[~,p]=min(abs(1*exp(1i*ph0+pi/2)-T_SLM1.*exp(1i*ph_SLM1)));
SLM1(1:2:N(1),1:2:N(2)) = Mapa1_1(p);
SLM1(1+1:2:N(1),1+1:2:N(2)) = Mapa1_1(p); % \
SLM1(1+1:2:N(1),1:2:N(2)) = Mapa2_1(p);
SLM1(1:2:N(1),1+1:2:N(2))= Mapa2_1(p); % /

% [~,p0]=min(abs(0-T_SLM1.*exp(1i*ph_SLM1)));
% SLM1(Y0(1):2:Y0(2),X0(1):2:X0(2)) = Mapa1_1(p0);
% SLM1(Y0(1)+1:2:Y0(2),X0(1)+1:2:X0(2)) = Mapa1_1(p0); % \
% SLM1(Y0(1)+1:2:Y0(2),X0(1):2:X0(2)) = Mapa2_1(p0);
% SLM1(Y0(1):2:Y0(2),X0(1)+1:2:X0(2))= Mapa2_1(p0); % /
% 
% [~,p1]=min(abs(C-T_SLM1.*exp(1i*ph_SLM1)));
% SLM1(Y1(1):2:Y1(2),X1(1):2:X1(2)) = Mapa1_1(p1);
% SLM1(Y1(1)+1:2:Y1(2),X1(1)+1:2:X1(2)) = Mapa1_1(p1); % \
% SLM1(Y1(1)+1:2:Y1(2),X1(1):2:X1(2)) = Mapa2_1(p1);
% SLM1(Y1(1):2:Y1(2),X1(1)+1:2:X1(2))= Mapa2_1(p1); % /

[~,px]=min(abs(C-T_SLM1.*exp(1i*ph_SLM1)));
SLM1(Y(1):2:Y(2),X(1):2:X(2)) = Mapa1_1(px);
SLM1(Y(1)+1:2:Y(2),X(1)+1:2:X(2)) = Mapa1_1(px); % \
SLM1(Y(1)+1:2:Y(2),X(1):2:X(2)) = Mapa2_1(px);
SLM1(Y(1):2:Y(2),X(1)+1:2:X(2))= Mapa2_1(px); % /

SLM1=(SLM1-1)/255;

% 
% imshow(SLM1')
% figure
% imshow(SLM2')
% figure
% imagesc(m1)
% figure
% imagesc(m2)
