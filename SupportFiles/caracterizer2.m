%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear variables;
close all

%Transmisions-----------------------------------------
Xi_1 = 250  ;  Xf_1 = 800 ;
Y1_1 = 320  ;  Y2_1 = 370;
Y3_1 = 430  ;  Y4_1 = 480;
Xi_2 = Xi_1 ;  Xf_2 = Xf_1;
Y1_2 = Y1_1 ;  Y2_2 = Y2_1;
Y3_2 = Y3_1 ;  Y4_2 = Y4_1;

ZERO = sum(sum(double(imread('ZEROcal.png'))));

label = '15122016' ; %input('   Label of the calibration?\n','s');
path  = [pwd '/' label '_I/'];


I1(1:256)=1;I2(1:256)=1;

SLM=1:2;

for j = SLM
    for i = 1:256
        filename = [path 'I' num2str(j) '_' num2str(i-1) '.png'];
        im       = double(imread(filename));
        
        if j==1
            I1(i)=sum(sum(im)) - ZERO;
        else
            I2(i)=sum(sum(im)) - ZERO;
        end
    end
end


A1=sqrt(I1/max(I1));A2=sqrt(I2/max(I2)); %from intensity to amplitude
h=figure;
plot(A1)
axis([1 256 0 1.1])
cd plots
print(h,'-depsc',[label 'amplitude_SLM1.eps']);
cd ..
h=figure;
plot(A2)
axis([1 256 0 1.1])
cd plots
print(h,'-depsc',[label 'amplitude_SLM2.eps']);
cd ..
%Phases-----------------------------------------------
Xi=Xi_1;Xf=Xf_1;
Y1=Y1_1;Y2=Y2_1;
Y3=Y3_1;Y4=Y4_1;

im=zeros(1024,768);
N=zeros(1,256);
M=zeros(1,256);
P=zeros(1,256);
phi=zeros(1,256);
Dx=Xf-Xi+1;
delta_phi=zeros(1,256);
aux=zeros(Dx,256);
aux2=zeros(Dx,256);


for k=1:2
    if k==1
        path=[pwd '/' label '_45/'];
        im=double(imread([path 'ph' num2str(2) '_30.png']));
        down_a(:,:)=im(Y3:Y4,Xi:Xf);
        down=mean(down_a,1);
        figure
        plot(abs(fft(down)));
        Pf=31;%input('Type the frequency of the first peak:\n  ');
        Q=1;%input('and the tolerance:\n  ');
    else
        path=[pwd '/' label '_135/'];
    end
for j=SLM
    for i=1:256
        filename=[path 'ph' num2str(j) '_' num2str(i-1) '.png'];
        im=double(imread(filename));
        %im=imread('img0.png');
        up_a(:,:)=im(Y1:Y2,Xi:Xf);
        %up_a(:,2)=im(Y2,Xi:Xf);
        up=mean(up_a,1);
        down_a(:,:)=im(Y3:Y4,Xi:Xf);
        %down_a(:,2)=im(Y4,Xi:Xf);
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
%         if(i==110)
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
        if j==1
            Dmean=mean(D);
        else
            Dmean=mean(D);
        end%--------------!!!!!
        Dstd=std(D);
        phi(i)=360*Dmean/P(i);
        delta_phi(i)=Dstd*360/P(i);
    end
    phi_0=min(phi(1:25));
    phi=phi-phi_0;
    phi=mod(phi,360);
    if k==1
        if j==1
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi1_45=phi;
            delta_phi1_45=delta_phi;
        else
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi2_45=phi;
            delta_phi2_45=delta_phi;
        end
    else
        if j==1
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi1_135=phi;
            delta_phi1_135=delta_phi;
        else
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi2_135=phi;
            delta_phi2_135=delta_phi;
        end
    end
end
end

phi1=(phi1_45+phi1_135)/2;
phi2=(phi2_45+phi2_135)/2;

h=figure;
plot(1:256,phi1_45,':b',1:256,phi1_135,':r',1:256,phi1,'-k')
cd plots
print(h,'-depsc',[label 'phase_SLM1.eps']);
cd ..
h=figure;
polar(phi1*pi/180,A1/max(A1))
cd plots
print(h,'-depsc',[label 'polar_SLM1.eps']);
cd ..


h=figure;
plot(1:256,phi2_45,':b',1:256,phi2_135,':r',1:256,phi2,'-k')
cd plots
print(h,'-depsc',[label 'phase_SLM2.eps']);
cd ..
h=figure;
polar(phi2*pi/180,A2/max(A2))
cd plots
print(h,'-depsc',[label 'polar_SLM2.eps']);
cd ..
