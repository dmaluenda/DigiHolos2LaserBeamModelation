%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Phases-----------------------------------------------
close all
Xi_1=400;Xf_1=800;
Y1_1=250;Y2_1=325;
Y3_1=425;Y4_1=550;
Xi=Xi_1;Xf=Xf_1;
Y1=Y1_1;Y2=Y2_1;
Y3=Y3_1;Y4=Y4_1;
Analizer=45;

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
    down_a(:,1)=im(Y3,Xi:Xf);
    down_a(:,2)=im(Y4,Xi:Xf);
    down=mean(down_a,2);
    figure
    plot(abs(fft(down)));
    Pf=14;%input('Type the frequency of the first peak:\n  ');
    Q=1;%input('and the tolerance:\n  ');
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
%             UP(Pf+Q:391-Q)=0;DOWN(Pf+Q:391-Q)=0;
%             UP(391+Q:Dx)=0;UP(391+Q:Dx)=0;
        end
        %plot(abs(UP))
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
        if j==1
            Dmean=-mean(-D);
        else
            Dmean=mean(-D);
        end%--------------!!!!!
        Dstd=std(D);
        phi(i)=360*Dmean/P(i);
        delta_phi(i)=Dstd*360/P(i);
    end
    phi=mod(phi,360);
    phi_0=min(phi(1:10));
    phi=phi-phi_0;
    phi=mod(phi,360);
    if j==1
%         phi_0=min(phi(1:35));
%         phi=phi-phi_0;
%         phi=mod(phi,360);
        phi1=phi;
        delta_phi1=delta_phi;
        if Analizer==45
            phi1_45=phi1;
        else
            phi1_135=phi1;
        end
    else
%         phi_0=min(phi(1:35));
%         phi=phi-phi_0;
%         phi=mod(phi,360);
        phi2=phi;
        delta_phi2=delta_phi;
        if Analizer==45
            phi2_45=phi2;
        else
            phi2_135=phi2;
        end
    end
end

h=figure;
plot(1:256,phi1,1:256,delta_phi1);
print(h,'-depsc',['phase_SLM1_' num2str(Analizer) '.eps']);
print(h,'-dpng',['phase_SLM1_' num2str(Analizer) '.png']);
h=figure;
polar(phi1*pi/180,I1/max(I1))
print(h,'-depsc',['polar_SLM1_' num2str(Analizer) '.eps']);
print(h,'-dpng',['polar_SLM1_' num2str(Analizer) '.png']);

h=figure;
plot(1:256,phi2,1:256,delta_phi2);
print(h,'-depsc',['phase_SLM2_' num2str(Analizer) '.eps']);
print(h,'-dpng',['phase_SLM2_' num2str(Analizer) '.png']);
h=figure;
polar(phi2*pi/180,I2/max(I2))
print(h,'-depsc',['polar_SLM2_' num2str(Analizer) '.eps']);
print(h,'-dpng',['polar_SLM2_' num2str(Analizer) '.png']);