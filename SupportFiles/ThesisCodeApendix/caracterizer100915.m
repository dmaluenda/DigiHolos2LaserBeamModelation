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

close all
firstPATH=pwd;
if ~exist('LV','var'), LV = 0; end

if LV==0
    path=['C:\Users\Lab541\Documents\BAUL - TRASH' '\260314_135'];
    label='Thesis';
    SLM=1:2; %SLM1 / SLM2
    AMP=0; %1 for amp / 0 for phase
    ANALIZER=1; %0 for 45º / 1 for 135º
end


for slm=SLM
for amp=AMP
for analizer=ANALIZER


%Transmisions-----------------------------------------
Xi_1=300;Xf_1=600;
Y1_1=300;Y2_1=350;
Y3_1=450;Y4_1=500;
Xi_2=Xi_1;Xf_2=Xf_1;
Y1_2=Y1_1;Y2_2=Y2_1;
Y3_2=Y3_1;Y4_2=Y4_1;

if amp==1
    if analizer==0
    
%path=[pwd '\' label '_I'];
im=double(imread([path '\I' num2str(slm) '_0.png']));
up1=im(Y1_1:Y2_1,Xi_1:Xf_1);
down1=im(Y3_1:Y4_1,Xi_1:Xf_1);
fact=sum(sum(up1))/sum(sum(down1));
I1(1:256)=0;



for i=1:256
    filename=[path '\I' num2str(slm) '_' num2str(i-1) '.png'];
    im=double(imread(filename));
    up=fact*im(Y1_1:Y2_1,Xi_1:Xf_1);
    down=im(Y3_1:Y4_1,Xi_1:Xf_1);
    I1(i)=sum(sum(up))/sum(sum(down));
end
 I1=1./I1; %flips signal for ref

I1=sqrt(I1./max(I1));%from intensity to amplitude
out=I1;
h=figure;
plot(out)
axis([1 256 0 1.1])
cd(path)
print(h,'-depsc',['amplitude_SLM' num2str(slm) '.eps']);
cd ..
    end

else
%Phases-----------------------------------------------
Xi=Xi_1;Xf=Xf_1;
Y1=Y1_1;Y2=Y2_1;
Y3=Y3_1;Y4=Y4_1;

N=zeros(1,256);
M=zeros(1,256);
P=zeros(1,256);
phi=zeros(1,256);
Dx=Xf-Xi+1;
delta_phi=zeros(1,256);
aux=zeros(Dx,256);
aux2=zeros(Dx,256);


k=-analizer+1;
% if k==1
%     path=[pwd '\' label '_45'];
% else
%     path=[pwd '\' label '_135'];
% end
im=double(imread([path '\ph' num2str(slm) '_30.png']));
down_a(:,:)=im(Y3:Y4,Xi:Xf);
down=mean(down_a,1);
DOWN=abs(fft(down));
[~,Pf]=max(DOWN(7:100));%input('Type the frequency of the first peak:\n  ');
Q=1;%input('and the tolerance:\n  ');
Pf=Pf+6;
% DOWN2=DOWN;DOWN2(2:Pf-Q)=0;DOWN2(Pf+Q:Dx-Pf+2-Q)=0;DOWN2(Dx-Pf+2+Q:Dx)=0;
% down2=abs(ifft(DOWN2));
% subplot(2,2,1);plot(down)
% subplot(2,2,2);plot(abs(DOWN))
% subplot(2,2,3);plot(down2)
% subplot(2,2,4);plot(abs(DOWN2))



for i=1:256
    filename=[path '\ph' num2str(slm) '_' num2str(i-1) '.png'];
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
%     if(i==110)
%         figure;plot(up2);
%         hold on;plot(down2);hold off
%     end
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
    if slm==1
        Dmean=-mean(D);
    else
        Dmean=-mean(D);
    end%--------------!!!!!
    Dstd=std(D);
    phi(i)=360*Dmean/P(i);
    delta_phi(i)=Dstd*360/P(i);
end
phi_0=min(phi(1:10));
phi=phi-phi_0;
out=mod(phi,360);

h=figure;
plot(out);hold on
plot(out+delta_phi,':b')
plot(out-delta_phi,':b')
hold off

axis([1 256 0 360])
cd(path)
print(h,'-depsc',['phase_SLM' num2str(slm) '.eps']);
cd ..

end


out_aux=1:max(size(out));

if slm==1
    if amp==0
        if analizer==0
            phi1_45=out;
        else
            phi1_135=out;
        end
    else
        A1=out;
    end
else
    if amp==0
        if analizer==0
            phi2_45=out;
        else
            phi2_135=out;
        end
    else
        A2=out;
    end
end


end
end
end

cd(firstPATH);