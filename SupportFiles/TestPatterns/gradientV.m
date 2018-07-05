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


%A=[263 -23];B=[450 182];C=[266 377];D=[72 177]; %CCD key points [px] 
a=[421 324];b=[308 215];c=[210 332];d=[328 424]; %SLM key points [mm]

%O_CCD=round([(B(1)+D(1))/2 (C(2)+A(2))/2]);
O_SLM=round([(a(1)+c(1))/2 (d(2)+b(2))/2]);
%D_CCD=(norm(B-D)+norm(C-A))/2;D_SLM=(norm(b-d)+norm(c-a))/2;
%step=D_SLM/D_CCD; %[mm]/[px]

%Oa=a+A*AD*step-CCD_size*AD*step;
%Ob=b+B*AD*step-CCD_size*AD*step;
%Oc=c+C*AD*step-CCD_size*AD*step;
%Od=d+D*AD*step-CCD_size*AD*step;
%O=round((Oa+Ob+Oc+Od)/4)*AD;
%O_CCD_SLM=round(O+O_CCD*AD*step+offset*AD*step);

%aa=round(O+A*AD*step);
%bb=round(O+B*AD*step);
%cc=round(O+C*AD*step);
%dd=round(O+D*AD*step);

SLM_size=[640 480];
SLM=zeros(SLM_size)';

% for i=1:SLM_size(1)
%     for j=c(1):a(1)
%         SLM(j,i)=(j-a(1))/(c(1)-a(1));
%     end
%     for j=1:c(1)
%         SLM(j,i)=1;
%     end
% end

for i=1:SLM_size(1)
    for j=c(1):floor((2*c(1)+a(1))/3)
        SLM(j,i)=3*(j-c(1))/(a(1)-c(1));
    end
    for j=ceil((2*c(1)+a(1))/3):a(1)
        SLM(j,i)=-3*(j-a(1))/(a(1)-c(1))/2;
    end
end

imshow(SLM)%,'colormap',colormap)
imwrite(SLM,'terçV.png','png')
