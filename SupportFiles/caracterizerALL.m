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



clear all;
close all

st=4;
% HWP=60;
% QWP=[165:5:175 167 172];

slm=1;
labelDIR='160214';
BAULdir=['C:\Users\Lab541\Documents\BAUL - TRASH\' labelDIR '_SLM' num2str(slm)];

SFpath=cd;
cd(BAULdir)
files = dir;
files=files(3:end,:);
cd(SFpath);

for file = files'
    file.name(end+1:10)='_';
    clear hwp qwp;
    if file.isdir==1 && file.name(2)=='_' && file.name(4)=='_'       %X_X_
        hwp=str2num(file.name(1));
        qwp=str2num(file.name(3));
    elseif file.isdir==1 && file.name(2)=='_' && file.name(5)=='_'   %X_XX_
        hwp=str2num(file.name(1));
        qwp=str2num([file.name(3) file.name(4)]);
    elseif file.isdir==1 && file.name(2)=='_' && file.name(6)=='_'   %X_XXX_
        hwp=str2num(file.name(1));
        qwp=str2num([file.name(3) file.name(4) file.name(5)]);
    elseif file.isdir==1 && file.name(3)=='_' && file.name(5)=='_'   %XX_X_
        hwp=str2num([file.name(1) file.name(2)]);
        qwp=str2num(file.name(4));
    elseif file.isdir==1 && file.name(3)=='_' && file.name(6)=='_'   %XX_XX_
        hwp=str2num([file.name(1) file.name(2)]);
        qwp=str2num([file.name(4) file.name(5)]);
    elseif file.isdir==1 && file.name(3)=='_' && file.name(7)=='_'   %XX_XXX_
        hwp=str2num([file.name(1) file.name(2)]);
        qwp=str2num([file.name(4) file.name(5) file.name(6)]);
    elseif file.isdir==1 && file.name(4)=='_' && file.name(6)=='_'   %XXX_X_
        hwp=str2num([file.name(1) file.name(2) file.name(3)]);
        qwp=str2num(file.name(5));
    elseif file.isdir==1 && file.name(4)=='_' && file.name(7)=='_'   %XXX_XX_
        hwp=str2num([file.name(1) file.name(2) file.name(3)]);
        qwp=str2num([file.name(5) file.name(6)]);
    elseif file.isdir==1 && file.name(4)=='_' && file.name(8)=='_'   %XXX_XXX_
        hwp=str2num([file.name(1) file.name(2) file.name(3)]);
        qwp=str2num([file.name(5) file.name(6) file.name(7)]);
    else
        display(['ERROR! The file `' file.name '´ do NOT match the format expected'])
    end

% for hwp=HWP
%     for qwp=QWP
%     close all

%Transmisions-----------------------------------------
% if (hwp<99&&mod(hwp,10)==0&&mod(qwp,15)==0&&~(hwp==40&&qwp==120))
%     label=[num2str(hwp) '_' num2str(qwp)];
%     path=['C:\Users\Lab541\Documents\BAUL - TRASH\310114_ALL\' label '_I\I1_'];
%     im=double(imread([path '0.png']));
%     Xi_1=300;Xf_1=600;
%     Y1_1=320;Y2_1=370;
%     Y3_1=430;Y4_1=480;
% else
    label=[num2str(hwp) '_' num2str(qwp)];
%     path=['C:\Users\Lab541\Documents\BAUL - TRASH\070214_SLM2\' label
%     '\M_'];
im=double(imread([BAULdir '\' label '\M_0.png']));
    Xi_1=150;Xf_1=350;
    Y1_1=170;Y2_1=370;
    Y3_1=430;Y4_1=630;
% end
Xi_2=Xi_1;Xf_2=Xf_1;
Y1_2=Y1_1;Y2_2=Y2_1;
Y3_2=Y3_1;Y4_2=Y4_1;

up1=im(Y1_1:Y2_1,Xi_1:Xf_1);
down1=im(Y3_1:Y4_1,Xi_1:Xf_1);
fact=sum(sum(up1))/sum(sum(down1));
I1=ones(1,256);


for i=1:st:256
    filename=[BAULdir '\' label '\M_' num2str(i-1) '.png'];
    im=double(imread(filename));
    up=fact*im(Y1_1:Y2_1,Xi_1:Xf_1);
    down=im(Y3_1:Y4_1,Xi_1:Xf_1);
    I1(i)=sum(sum(down))/sum(sum(up));
end

% if slm~=1
% I1=1./I1;%for flip ref to signal
% end

A1=sqrt(I1(1:st:256)/max(I1(1:st:256)));%from intensity to amplitude
h1=figure;
plot(1:st:256,A1,'LineWidth',2)
title(['\lambda/2=' num2str(hwp) 'º ; \lambda/4=' num2str(qwp) 'º'])
axis([1 256 0 1.05])
cd plots
[~,~,~]=mkdir([labelDIR '_SLM' num2str(slm)]);
cd([labelDIR '_SLM' num2str(slm)]);
print(h1,'-dpng',[num2str(hwp) '_' num2str(qwp) '_T_SLM' num2str(slm) '.png']);
cd(SFpath)
close(h1)
% h=figure;
% plot(I2)
% axis([1 256 0 1.1])
% cd plots
% print(h,'-depsc','amplitude_SLM2.eps');
% cd ..


%Phases-----------------------------------------------
% if (hwp<99&&mod(hwp,10)==0&&mod(qwp,15)==0&&~(hwp==40&&qwp==120))
%     Xi=Xi_1;Xf=Xf_1;
%     Y1=Y1_1;Y2=Y2_1;
%     Y3=Y3_1;Y4=Y4_1;
%     
%     label=[num2str(hwp) '_' num2str(qwp)];
%     path=['C:\Users\Lab541\Documents\BAUL - TRASH\310114_ALL\' label '_45\ph1_'];
%     im=double(imread([path '32.png']));
%     down_a=im(Y3:Y4,Xi:Xf);
%     down=mean(down_a,1);
%     DOWN=abs(fft(down));
%     [~,Pf]=max(DOWN(5:100));%input('Type the frequency of the first peak:\n  ');
%     Q=1;%input('and the tolerance:\n  ');
%     Pf=Pf+4;
% 
% else
    Xi=450;Xf=900;
    Y1=170;Y2=370;
    Y3=430;Y4=630;
    label=[num2str(hwp) '_' num2str(qwp)];
    im=double(imread([BAULdir '\' label '\M_0.png']));
    down_a=im(Y3:Y4,Xi:Xf);
    down=mean(down_a,1);
    DOWN=abs(fft(down));
    [~,Pf]=max(DOWN(7:100));%input('Type the frequency of the first peak:\n  ');
    Q=1;%input('and the tolerance:\n  ');
    Pf=Pf+6;
% end



im=zeros(1024,768);
N=zeros(1,256);
M=zeros(1,256);
P=zeros(1,256);
phi=zeros(1,256);
Dx=Xf-Xi+1;
delta_phi=zeros(1,256);
aux=zeros(Dx,256);
aux2=zeros(Dx,256);

for i=1:st:256
    filename=[BAULdir '\' label '\M_' num2str(i-1) '.png'];
    im=double(imread(filename));
    up_a=im(Y1:Y2,Xi:Xf);
    up=mean(up_a,1);
    down_a=im(Y3:Y4,Xi:Xf);
    down=mean(down_a,1);
%      figure;plot(up)
    UP=fft(up);DOWN=fft(down);
    aux(:,i)=UP;aux2(:,i)=DOWN;
%      figure;plot(abs(UP))
    if Pf~=0
        UP(2:Pf-Q)=0;DOWN(2:Pf-Q)=0;
        %UP(Pf+Q:Dx)=0;DOWN(Pf+Q:Dx)=0;
        UP(Pf+Q:Dx-Pf+2-Q)=0;DOWN(Pf+Q:Dx-Pf+2-Q)=0;
        UP(Dx-Pf+2+Q:Dx)=0;DOWN(Dx-Pf+2+Q:Dx)=0;
    end
    up2=abs(ifft(UP));down2=abs(ifft(DOWN));
%         if(i==109)
%             figure;plot(up2);
%             hold on;plot(down2);hold off
%         end
    [peakU,peakUP]=findpeaks(up2);
    [peakD,peakDW]=findpeaks(down2);
    P(i)=peakUP(5)-peakUP(4);
    n=size(peakUP);
    m=size(peakDW);
    N(i)=n(2);M(i)=m(2);
    if N(i)>M(i)
        D=peakUP(1:M(i))-peakDW(1:M(i));
    elseif M(i)>N(i)
        D=peakUP(1:N(i))-peakDW(1:N(i));
    else
        D=peakUP-peakDW;
    end

        Dmean=mean(D);

    Dstd=std(D);
    phi(i)=360*Dmean/P(i);
    delta_phi(i)=Dstd*360/P(i);
end

phi_0=min(phi(1:25));
phi=phi-phi_0;
phi=mod(phi,360);

%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi1_45=phi;
            delta_phi1_45=delta_phi;



phi1_45=mod(phi1_45-phi1_45(1),360);
phi1=(phi1_45);%+phi1_135)/2;
% phi2=(phi2_45+phi2_135)/2;


h2=figure;
plot(1:st:256,phi1_45(1:st:256),'LineWidth',2)%,':b',1:256,phi1_135,':r',1:256,phi1,'-k')
title(['\lambda/2=' num2str(hwp) 'º ; \lambda/4=' num2str(qwp) 'º'])
axis([1 256 0 360])
cd plots
[~,~,~]=mkdir([labelDIR '_SLM' num2str(slm)]);
cd([labelDIR '_SLM' num2str(slm)]);
print(h2,'-dpng',[num2str(hwp) '_' num2str(qwp) '_ph_SLM' num2str(slm) '.png']);
cd(SFpath);
close(h2)
h3=figure;
polar(phi1(1:st:256)*pi/180,A1/max(A1))
title(['\lambda/2=' num2str(hwp) 'º ; \lambda/4=' num2str(qwp) 'º'])
set(findobj(gca, 'Type', 'line'),'linewidth',2);
cd plots
[~,~,~]=mkdir([labelDIR '_SLM' num2str(slm)]);
cd([labelDIR '_SLM' num2str(slm)]);
print(h3,'-dpng',[num2str(hwp) '_' num2str(qwp) '_polar_SLM' num2str(slm) '.png']);
close(h3)

phi_def=phi1_45(1:st:256)*pi/180;
T_def=A1/max(A1);
curve=[(1:st:256)' T_def' phi_def'];

fid = fopen([num2str(hwp) '_' num2str(qwp) '_curve2.txt'],'wt');
fprintf(fid,'%3.0f %6.4f% 6.4f\n',curve');
fclose(fid);
cd(SFpath);

% h=figure;
% plot(1:256,phi2_45,':b',1:256,phi2_135,':r',1:256,phi2,'-k')
% cd plots
% print(h,'-depsc','phase_SLM2.eps');
% cd ..
% h=figure;
% polar(phi2*pi/180,I2/max(I2))
% cd plots
% print(h,'-depsc','polar_SLM2.eps');
% cd ..
%     end
end
