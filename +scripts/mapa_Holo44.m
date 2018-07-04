%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubó - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hologram generator
% grayL = f(Trans_1, Trans_2, Phase )
% [0 1] = f( [0 1] ,  [0 1] ,[0 2pi])
function[SLM1,SLM2]=mapa_Holo44(Trans1,Trans2,Phase1,Phase2)
disp('Calculating the holograms, this may take some time. Take it easy!');


% --- Just for script to debug ---
% clear all;
% [Trans1,Trans2,Phase1,Phase2]=scripts.beam_design([1024 768],23,0);
% --------------------------------

% sizes
N = size(Phase1);
X = [1 768];
Y = [1 1024];

% loading maps of codificable values
PATH = [pwd '\SupportFiles\ComplexValues4_10_SLM'];

Amp_max1 = 0.5;
Amp_max2 = 0.5;
A_max = min([Amp_max1 Amp_max2]);

data1 = load([PATH num2str(1) '.txt']);
T_SLM1  = data1(:,1)/10000;
ph_SLM1 = mod(data1(:,2)/10000,2*pi);
MapaI_1 = data1(:,3);
MapaJ_1 = data1(:,4);
MapaK_1 = data1(:,5);
MapaL_1 = data1(:,6);

data2 = load([PATH num2str(2) '.txt']);
T_SLM2  = data2(:,1)/10000;
ph_SLM2 = mod(data2(:,2)/10000,2*pi);
MapaI_2 = data2(:,3);
MapaJ_2 = data2(:,4);
MapaK_2 = data2(:,5);
MapaL_2 = data2(:,6);

% accesible values
C_SLM1=T_SLM1.*exp(1i*ph_SLM1);
C_SLM2=T_SLM2.*exp(1i*ph_SLM2);

% dasirable values
C1 = Trans1.*exp(1i*Phase1).*A_max;
C2 = Trans2.*exp(1i*Phase2).*A_max;

%---resizing-for-macropixel-procedure----------------
C1_mean(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4) = ...
    (C1(Y(1):4:Y(2),X(1):4:X(2)) + ...
    C1(Y(1):4:Y(2),X(1)+1:4:X(2)) + ...
    C1(Y(1):4:Y(2),X(1)+2:4:X(2)) + ...
    C1(Y(1):4:Y(2),X(1)+3:4:X(2)) + ... %11
    C1(Y(1)+1:4:Y(2),X(1):4:X(2)) + ...
    C1(Y(1)+1:4:Y(2),X(1)+1:4:X(2)) + ...
    C1(Y(1)+1:4:Y(2),X(1)+2:4:X(2)) + ...
    C1(Y(1)+1:4:Y(2),X(1)+3:4:X(2)) + ... %12
    C1(Y(1)+2:4:Y(2),X(1):4:X(2)) + ...
    C1(Y(1)+2:4:Y(2),X(1)+1:4:X(2)) + ...
    C1(Y(1)+2:4:Y(2),X(1)+2:4:X(2)) + ...
    C1(Y(1)+2:4:Y(2),X(1)+3:4:X(2)) + ... %21
    C1(Y(1)+3:4:Y(2),X(1):4:X(2)) + ...
    C1(Y(1)+3:4:Y(2),X(1)+1:4:X(2)) + ...
    C1(Y(1)+3:4:Y(2),X(1)+2:4:X(2)) + ...
    C1(Y(1)+3:4:Y(2),X(1)+3:4:X(2)))/16 ; %22
C2_mean(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4) = ...
    (C2(Y(1):4:Y(2),X(1):4:X(2)) + ...
    C2(Y(1):4:Y(2),X(1)+1:4:X(2)) + ...
    C2(Y(1):4:Y(2),X(1)+2:4:X(2)) + ...
    C2(Y(1):4:Y(2),X(1)+3:4:X(2)) + ... %11
    C2(Y(1)+1:4:Y(2),X(1):4:X(2)) + ...
    C2(Y(1)+1:4:Y(2),X(1)+1:4:X(2)) + ...
    C2(Y(1)+1:4:Y(2),X(1)+2:4:X(2)) + ...
    C2(Y(1)+1:4:Y(2),X(1)+3:4:X(2)) + ... %12
    C2(Y(1)+2:4:Y(2),X(1):4:X(2)) + ...
    C2(Y(1)+2:4:Y(2),X(1)+1:4:X(2)) + ...
    C2(Y(1)+2:4:Y(2),X(1)+2:4:X(2)) + ...
    C2(Y(1)+2:4:Y(2),X(1)+3:4:X(2)) + ... %21
    C2(Y(1)+3:4:Y(2),X(1):4:X(2)) + ...
    C2(Y(1)+3:4:Y(2),X(1)+1:4:X(2)) + ...
    C2(Y(1)+3:4:Y(2),X(1)+2:4:X(2)) + ...
    C2(Y(1)+3:4:Y(2),X(1)+3:4:X(2)))/16 ; %22

SLM1 = zeros(N);
SLM2 = zeros(N);

m1 = zeros(size(C1_mean)); 
m2 = zeros(size(C2_mean));
p1 = zeros(size(C1_mean));
p2 = zeros(size(C2_mean));
aux=0;
time0 = clock;
for i=1:(Y(2)-Y(1)+1)/4
   time1 = clock;
   for j = 1:(X(2)-X(1)+1)/4
        [m1(i,j), p1(i,j)] = min(abs(C1_mean(i,j)-C_SLM1));
        [m2(i,j), p2(i,j)] = min(abs(C2_mean(i,j)-C_SLM2));
   end
   dis=round(i/(Y(2)-Y(1)+1)*2*100);
   if i==1
       time2 = clock;
       interval1 = time2-time1;
       TIME = interval1*(Y(2)-Y(1)+1)/2;
       disp(['Estimated time: ' num2str(TIME(5)+floor(TIME(6)/60)) 'min ' ...
           num2str(floor(mod(TIME(6),60))) 'sec']);
   end
   if dis~=aux && mod(dis,5)==0
       aux = dis;
       time3 = clock;
       interval2 = time3-time1;
       TIME = interval2*((Y(2)-Y(1)+1)/2-i);
       if dis==100
          TIME = clock-time0;
          disp(['Done! Total time spent: ' num2str(TIME(5)+floor(TIME(6)/60)) 'min ' ...
              num2str(floor(mod(TIME(6),60))) 'sec left']); 
       else
          disp([num2str(dis) '% done: ' num2str(TIME(5)+floor(TIME(6)/60)) 'min ' ...
              num2str(floor(mod(TIME(6),60))) 'sec left']);
       end
   end
end



SLM1(Y(1):4:Y(2),X(1):4:X(2)) = ...
    MapaI_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));   %I
SLM1(Y(1)+1:4:Y(2),X(1)+1:4:X(2)) = ...
    MapaI_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM1(Y(1)+2:4:Y(2),X(1)+2:4:X(2)) = ...
    MapaI_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM1(Y(1)+3:4:Y(2),X(1)+3:4:X(2)) = ...
    MapaI_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));

SLM1(Y(1)+1:4:Y(2),X(1):4:X(2)) = ...
    MapaJ_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));  %J
SLM1(Y(1)+2:4:Y(2),X(1)+1:4:X(2)) = ...
    MapaJ_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4)); 
SLM1(Y(1)+3:4:Y(2),X(1)+2:4:X(2)) = ...
    MapaJ_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4)); 
SLM1(Y(1):4:Y(2),X(1)+3:4:X(2)) = ...
    MapaJ_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4)); 

SLM1(Y(1)+2:4:Y(2),X(1):4:X(2)) = ...
    MapaK_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));  %K
SLM1(Y(1)+3:4:Y(2),X(1)+1:4:X(2)) = ...
    MapaK_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM1(Y(1):4:Y(2),X(1)+2:4:X(2)) = ...
    MapaK_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM1(Y(1)+1:4:Y(2),X(1)+3:4:X(2)) = ...
    MapaK_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));

SLM1(Y(1)+3:4:Y(2),X(1):4:X(2))= ...
    MapaL_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));  %L
SLM1(Y(1):4:Y(2),X(1)+1:4:X(2))= ...
    MapaL_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM1(Y(1)+1:4:Y(2),X(1)+2:4:X(2))= ...
    MapaL_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM1(Y(1)+2:4:Y(2),X(1)+3:4:X(2))= ...
    MapaL_1(p1(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));

%SLM2

SLM2(Y(1):4:Y(2),X(1):4:X(2)) = ...
    MapaI_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));   %I
SLM2(Y(1)+1:4:Y(2),X(1)+1:4:X(2)) = ...
    MapaI_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM2(Y(1)+2:4:Y(2),X(1)+2:4:X(2)) = ...
    MapaI_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM2(Y(1)+3:4:Y(2),X(1)+3:4:X(2)) = ...
    MapaI_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));

SLM2(Y(1)+1:4:Y(2),X(1):4:X(2)) = ...
    MapaJ_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));  %J
SLM2(Y(1)+2:4:Y(2),X(1)+1:4:X(2)) = ...
    MapaJ_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4)); 
SLM2(Y(1)+3:4:Y(2),X(1)+2:4:X(2)) = ...
    MapaJ_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4)); 
SLM2(Y(1):4:Y(2),X(1)+3:4:X(2)) = ...
    MapaJ_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4)); 

SLM2(Y(1)+2:4:Y(2),X(1):4:X(2)) = ...
    MapaK_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));  %K
SLM2(Y(1)+3:4:Y(2),X(1)+1:4:X(2)) = ...
    MapaK_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM2(Y(1):4:Y(2),X(1)+2:4:X(2)) = ...
    MapaK_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM2(Y(1)+1:4:Y(2),X(1)+3:4:X(2)) = ...
    MapaK_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));

SLM2(Y(1)+3:4:Y(2),X(1):4:X(2))= ...
    MapaL_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));  %L
SLM2(Y(1):4:Y(2),X(1)+1:4:X(2))= ...
    MapaL_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM2(Y(1)+1:4:Y(2),X(1)+2:4:X(2))= ...
    MapaL_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));
SLM2(Y(1)+2:4:Y(2),X(1)+3:4:X(2))= ...
    MapaL_2(p2(1:(Y(2)-Y(1)+1)/4,1:(X(2)-X(1)+1)/4));




SLM1 = (SLM1)/255;
SLM2 = (SLM2)/255;

[SLM1,SLM2] = scripts.rotate_SLM(SLM1,SLM2);

% imshow(SLM1')
% figure
% imshow(SLM2')
% figure
% imagesc(m1)
% figure
% imagesc(m2)
