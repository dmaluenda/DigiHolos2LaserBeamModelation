%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Complex values. (2012)
%David Maluenda Niubó - Applied Physics and Optics (UB)

clear B1 aux1_i aux1_j B2 aux2_i aux2_j; close all
%clear all;close all
%k=2;
date=clock;
N1=255;N2=255;
A_max=1;%just for ploting
%usefull values and SemiCercle
phi1_0=45*pi/180;  %    <-Rotate
phi2_0=45*pi/180;  %    <-Rotate
A1_max=0.35; %         %  <-Trim1
A2_max=0.35;  %         %  <-Trim2

mapa1=load(['valors_p' num2str(1) '.txt']);
ph1=mapa1(1:N1,3);A1_0=mapa1(1:N1,2);
mapa2=load(['valors_p' num2str(2) '.txt']);
ph2=mapa2(1:N2,3);A2_0=mapa2(1:N2,2);

figure;polar(ph2,A2_0)
A1=A1_0.*exp(1i*ph1);
A2=A2_0.*exp(1i*ph2);

% B=zeros(1,(N*(N-1))/2);
% aux_i=zeros(1,(N*(N-1))/2);
% aux_j=zeros(1,(N*(N-1))/2);

count=1;
for i=1:N1
    for j=i:N1
        aux=(A1(i)+A1(j))/2;
        if abs(aux)<=A1_max
            B1(count)=aux;
            aux1_i(count)=i;
            aux1_j(count)=j;
            count=count+1;
        end
    end
end
count=1;
for i=1:N2
    for j=i:N2
        aux=(A2(i)+A2(j))/2;
        if abs(aux)<=A2_max
            B2(count)=aux;
            aux2_i(count)=i;
            aux2_j(count)=j;
            count=count+1;
        end
    end
end

% M=(N*(N-1))/2;
% C=zeros(1,(M*(M-1))/20);
% count=1;
% for i=1:5:M
%     for j=i+1:5:M
%         C(count)=(B(i)+B(j))/2;
%         count=count+1;
%     end
% end

%Cercles to plot
phi1_SC=phi1_0:0.05:pi+phi1_0;
A1_SC=A_max*ones(size(phi1_SC));
phi2_SC=phi2_0:0.05:pi+phi2_0;
A2_SC=A_max*ones(size(phi2_SC));

A_bottom=A_max:-0.01:-A_max;
phi1_bottom=(phi1_0+angle(A_bottom)).*ones(size(A_bottom));
phi2_bottom=(phi2_0+angle(A_bottom)).*ones(size(A_bottom));

%ploting SLM1
h=figure;
polar(angle(A1),abs(A1),'-r') %Real curve
hold on
polar(angle(B1),abs(B1),'+b') %Possibles values
polar(phi1_SC,A1_SC,'-k') %Semicircle
polar(pi+phi1_SC,A1_SC,'-k') %Semicircle
polar(phi1_bottom,abs(A_bottom),'-k'); %bottom SC
plot(0,0,'*k','markersize',10) %Center
title (['SLM ' num2str(1)])
legend 'SLM curve' 'Possible coding values' 'Usefull zone'
hold off
cd plots
print(h,'-depsc',[num2str(date(3)) num2str(date(2)) ...
    num2str(date(1)-2000) '_polar_SLM1.eps']);
cd ..
% figure
% polar(angle(C),abs(C),'x')

%Creating data for SLM1
rho1=abs(B1)/A1_max;phi1=mod(angle(B1)-phi1_0,2*pi);
data=[rho1;phi1;aux1_i;aux1_j];
fid = fopen(['ComplexValues_SLM' num2str(1) '.txt'],'wt');
fprintf(fid,'%10.20f  %10.20f  %3.0f  %3.0f\n',data);
fclose(fid);




%ploting SLM2
h=figure;
polar(angle(A2),abs(A2),'-r') %Real curve
hold on
polar(angle(B2),abs(B2),'+b') %Possibles values
polar(phi2_SC,A2_SC,'-k') %Semicircle
polar(pi+phi2_SC,A2_SC,'-k') %Semicircle
polar(phi2_bottom,abs(A_bottom),'-k'); %bottom SC
plot(0,0,'*k','markersize',10) %Center
title (['SLM ' num2str(2)])
legend 'SLM curve' 'Possible coding values' 'Usefull zone'
hold off
cd plots
print(h,'-depsc',[num2str(date(3)) num2str(date(2)) ...
    num2str(date(1)-2000) '_polar_SLM2.eps']);
cd ..


%Creating data for SLM2
rho2=abs(B2)/A2_max;phi2=mod(angle(B2)-phi2_0,2*pi);
data=[rho2;phi2;aux2_i;aux2_j];
fid = fopen(['ComplexValues_SLM' num2str(2) '.txt'],'wt');
fprintf(fid,'%10.20f  %10.20f  %3.0f  %3.0f\n',data);
fclose(fid);