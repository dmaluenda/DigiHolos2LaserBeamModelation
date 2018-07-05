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

%AD=[0 1;1 0];

%% CCD key points [px] 
% A = [263 -23];
% B = [450 182];
% C = [266 377];
% D = [72 177]; 

%% SLM key points [mm]
a = [421 324];
b = [308 215];
c = [210 332];
d = [328 424]; 

% O_CCD = round([(B(1)+D(1))/2 (C(2)+A(2))/2]);
O_SLM = round([(a(1)+c(1))/2 (d(2)+b(2))/2]);
% D_CCD = (norm(B-D)+norm(C-A))/2;D_SLM=(norm(b-d)+norm(c-a))/2;
% step = D_SLM/D_CCD; %[mm]/[px]

[xx,yy] = meshgrid(-O_SLM(1):640-O_SLM(1)-1, -O_SLM(2):480-O_SLM(2)-1);

%Oa = a+A*AD*step - CCD_size*AD*step;
%Ob = b+B*AD*step - CCD_size*AD*step;
%Oc = c+C*AD*step - CCD_size*AD*step;
%Od = d+D*AD*step - CCD_size*AD*step;
%O = round((Oa+Ob+Oc+Od)/4)*AD;
%O_CCD_SLM = round(O+O_CCD*AD*step+offset*AD*step);

%aa = round(O+A*AD*step);
%bb = round(O+B*AD*step);
%cc = round(O+C*AD*step);
%dd = round(O+D*AD*step);

SLM_size = [640 480];
SLM = zeros(SLM_size)';

theta = atan(yy./xx);

% for i=1:SLM_size(1)
%     for j=1:SLM_size(2)
%         SLM(i,j) = theta(i,j)    
%     end
% end

imshow(abs(theta)/max(max(theta)))%,'colormap',colormap)
imwrite(abs(theta)/max(max(theta)),'RadialV.png','png')
