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

%Transmisions-----------------------------------------
Xi_1=360;Xf_1=540;
Y1_1=225;Y2_1=325;
Y3_1=400;Y4_1=500;
Xi_2=360;Xf_2=540;
Y1_2=225;Y2_2=325;
Y3_2=400;Y4_2=500;

im=double(imread('I1_0.png'));
up=fact1*im(Y1_1:Y2_1,Xi_1:Xf_1);
down=im(Y3_1:Y4_1,Xi_1:Xf_1);
fact1=sum(sum(up))/sum(sum(down));
im=double(imread('I2_0.png'));
up=fact1*im(Y1_1:Y2_1,Xi_1:Xf_1);
down=im(Y3_1:Y4_1,Xi_1:Xf_1);
fact2=sum(sum(up))/sum(sum(down));
I1(1:256)=0;I2(1:256)=0;

for i=1:256
    filename=['I1_' num2str(i-1) '.png'];
    im=double(imread(filename));
    up=fact1*im(Y1_1:Y2_1,Xi_1:Xf_1);
    down=im(Y3_1:Y4_1,Xi_1:Xf_1);
    I1(i)=sum(sum(down))/sum(sum(up));
end
for i=1:256
    filename=['I2_' num2str(i-1) '.png'];
    im=double(imread(filename));
    up=fact2*im(Y1_2:Y2_2,Xi_2:Xf_2);
    down=im(Y3_2:Y4_2,Xi_2:Xf_2);
    I2(i)=sum(sum(down))/sum(sum(up));
end

h=figure;
plot(I1/max(I1))
axis([1 256 0 1.1])
print(h,'-depsc','amplitude_SLM1.eps');
print(h,'-dpng','amplitude_SLM1.png');
h=figure;
plot(I2/max(I2))
axis([1 256 0 1.1])
print(h,'-depsc','amplitude_SLM2.eps');
print(h,'-dpng','amplitude_SLM2.png');

%Phases-----------------------------------------------
Xi=360;Xf=540;
Y1=225;Y2=325;
Y3=400;Y4=500;

im=zeros(1024,768);
N=zeros(1,256);
M=zeros(1,256);
P=zeros(1,256);
phi=zeros(1,256);
Dx=Xf-Xi+1;
delta_phi=zeros(1,256);
aux=zeros(Dx,256);
aux2=zeros(Dx,256);

for j=1:2
    im=double(imread(['ph' num2str(j) '_30.png']));
    plot(abs(fft(im)));
    Pf=input('Type the frequency of the first peak:\n  ');
    Q=input('and the tolerance:\n  ');
    for i=1:256
        filename=['ph' num2str(j) '_' num2str(i-1) '.png'];
        im=double(imread(filename));
        %im=imread('img0.png');
        up_a(:,1)=im(Y1,Xi:Xf);
        up_a(:,2)=im(Y2,Xi:Xf);
        up=mean(up_a,2);
        down_a(:,1)=im(Y3,Xi:Xf);
        down_a(:,2)=im(Y4,Xi:Xf);
        down=mean(down_a,2);
        UP=fft(up);DOWN=fft(down);
        aux(:,i)=UP;aux2(:,i)=DOWN;
        if Pf~=0
            UP(2:Pf-Q)=0;DOWN(2:Pf-Q)=0;
            UP(Pf+Q:Dx)=0;DOWN(Pf+Q:Dx)=0;
        end
        up2=abs(ifft(UP));down2=abs(ifft(DOWN));
        %figure;plot(1:401,up2,1:401,down2)

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
        Dmean=mean(D);%--------------!!!!!
        Dstd=std(D);
        phi(i)=360*Dmean/P(i);
        delta_phi(i)=Dstd*360/P(i);
    end
    phi_0=min(phi(1:20));
    phi=phi-phi_0;
    phi=mod(phi,360);
    if j==1
        phi1=phi;
        delta_phi1=delta_phi;
    else
        phi2=phi;
        delta_phi2=delta_phi;
    end
end

h=figure;
plot(1:256,phi1,1:256,delta_phi1);
print(h,'-depsc','phase_SLM1.eps');
print(h,'-dpng','phase_SLM1.png');
h=figure;
polar(phi1*pi/180,I1/max(I1))
print(h,'-depsc','polar_SLM1.eps');
print(h,'-dpng','polar_SLM1.png');

h=figure;
plot(1:256,phi2,1:256,delta_phi2);
print(h,'-depsc','phase_SLM2.eps');
print(h,'-dpng','phase_SLM2.png');
h=figure;
polar(phi2*pi/180,I2/max(I2))
print(h,'-depsc','polar_SLM2.eps');
print(h,'-dpng','polar_SLM2.png');