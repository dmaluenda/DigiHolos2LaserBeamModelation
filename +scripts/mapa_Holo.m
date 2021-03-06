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

% Hologram generator
%
%  grayL = f(Trans_1, Trans_2, Phase )
%  [0 1] = f( [0 1] ,  [0 1] ,[0 2pi])

function[SLM1,SLM2]=mapa_Holo(Trans1,Trans2,Phase1,Phase2,ModulationType)

% --- Just for scripting to debug ---
% clear all;
% [Trans1,Trans2,Phase1,Phase2]=scripts.beam_design([1024 768],1);
% -----------------------------------

if ~strcmp(ModulationType,'Amplitude') && ~strcmp(ModulationType,'Complex')
    error('ModulationType must be either "Amplitude" or "Complex"')
end
PATH = [pwd '/SupportFiles/' ModulationType 'Values_SLM'];

% sizes
N = size(Phase1);

% ROI
X = [1 768];
Y = [1 1024];

% loading maps of codificable values
Amp_max1 = 1;
Amp_max2 = 1;
A_max    = min( [Amp_max1 Amp_max2] );

data1 =load([PATH '1.txt']);
T_SLM1  = data1(:,1);
ph_SLM1 = mod(data1(:,2),2*pi);
Mapa1_1 = data1(:,3);
Mapa2_1 = data1(:,4);

data2 = load([PATH '2.txt']);
T_SLM2  = data2(:,1);
ph_SLM2 = mod(data2(:,2),2*pi);
Mapa1_2 = data2(:,3);
Mapa2_2 = data2(:,4);

% accesible values
C_SLM1 = T_SLM1.*exp(1i*ph_SLM1);
C_SLM2 = T_SLM2.*exp(1i*ph_SLM2);

% desirable values
C1 = Trans1.*exp(1i*Phase1).*A_max; 
C2 = Trans2.*exp(1i*Phase2).*A_max;

%---resizing-for-macropixel-procedure----------------
C1_mean(1:(Y(2)-Y(1)+1)/2,1:(X(2)-X(1)+1)/2) = ... 
    ( C1(Y(1)  :2:Y(2), X(1)  :2:X(2)) + ...
      C1(Y(1)+1:2:Y(2), X(1)  :2:X(2)) + ...
      C1(Y(1)  :2:Y(2), X(1)+1:2:X(2)) + ...
      C1(Y(1)+1:2:Y(2), X(1)+1:2:X(2))  )/4 ;

C2_mean(1:(Y(2)-Y(1)+1)/2,1:(X(2)-X(1)+1)/2) = ...
    ( C2(Y(1)  :2:Y(2), X(1)  :2:X(2)) + ...
      C2(Y(1)+1:2:Y(2), X(1)  :2:X(2)) + ...
      C2(Y(1)  :2:Y(2), X(1)+1:2:X(2)) + ... 
      C2(Y(1)+1:2:Y(2), X(1)+1:2:X(2))  )/4 ;

SLM1=zeros(N);
SLM2=zeros(N);

m1 = ones(size(C1_mean)); 
m2 = ones(size(C2_mean));
p1 = ones(size(C1_mean));
p2 = ones(size(C2_mean));
aux = 0;
time0 = clock;
for i=1:(Y(2)-Y(1)+1)/2
   time1 = clock;
   for j=1:(X(2)-X(1)+1)/2
        [m1(i,j), p1(i,j)] = min(abs(C1_mean(i,j)-C_SLM1));
        [m2(i,j), p2(i,j)] = min(abs(C2_mean(i,j)-C_SLM2));
   end

   dis=round(i/(Y(2)-Y(1)+1)*2*100);
   if i==1
       time2 = clock;
       interval1 = time2-time1;
       TIME = interval1*(Y(2)-Y(1)+1)/2;
   end
   if dis~=aux && mod(dis,5)==0
       aux = dis;
       time3 = clock;
       interval2 = time3-time1;
       TIME = interval2*((Y(2)-Y(1)+1)/2-i);
   end
end


% M_L (top-left) for SLM1
SLM1(  Y(1) :2:Y(2) , X(1):2:X(2) ) = ...
    Mapa1_1( p1(1:(Y(2)-Y(1)+1)/2,1:(X(2)-X(1)+1)/2) );

% M_L (botom-right) for SLM1
SLM1( Y(1)+1:2:Y(2) , X(1)+1:2:X(2) ) = ...
    Mapa1_1( p1(1:(Y(2)-Y(1)+1)/2, 1:(X(2)-X(1)+1)/2)); % \

% M_R (top-right) for SLM1
SLM1( Y(1)+1:2:Y(2) ,  X(1) :2:X(2) ) = ...
    Mapa2_1( p1(1:(Y(2)-Y(1)+1)/2, 1:(X(2)-X(1)+1)/2));

% M_R (botom-left) for SLM1
SLM1(  Y(1) :2:Y(2) , X(1)+1:2:X(2) ) = ...
    Mapa2_1( p1(1:(Y(2)-Y(1)+1)/2, 1:(X(2)-X(1)+1)/2)); % /


% M_L (top-left) for SLM2
SLM2(  Y(1) :2:Y(2) ,  X(1) :2:X(2) ) = ...
    Mapa1_2( p2(1:(Y(2)-Y(1)+1)/2, 1:(X(2)-X(1)+1)/2));

% M_L (botom-right) for SLM2   
SLM2( Y(1)+1:2:Y(2) , X(1)+1:2:X(2) ) = ...
    Mapa1_2( p2(1:(Y(2)-Y(1)+1)/2, 1:(X(2)-X(1)+1)/2)); % \

% M_R (top-right) for SLM2    
SLM2( Y(1)+1:2:Y(2) ,  X(1) :2:X(2) ) = ...
    Mapa2_2( p2(1:(Y(2)-Y(1)+1)/2, 1:(X(2)-X(1)+1)/2));

% M_R (botom-left) for SLM1
SLM2(  Y(1) :2:Y(2) , X(1)+1:2:X(2) ) = ...
    Mapa2_2( p2(1:(Y(2)-Y(1)+1)/2, 1:(X(2)-X(1)+1)/2)); % /  


SLM1 = (SLM1-1)/255;
SLM2 = (SLM2-1)/255;

[SLM1,SLM2] = scripts.rotate_SLM(SLM1,SLM2);

% imshow(SLM1')
% figure
% imshow(SLM2')
% figure
% imagesc(m1)
% figure
% imagesc(m2)
