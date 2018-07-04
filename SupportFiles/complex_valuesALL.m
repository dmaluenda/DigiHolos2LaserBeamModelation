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


cd plots
cd 160214_SLM1;

st=1;st0=4;

HWP=[62];
QWP=[165:5:175 167 172];

% for hwp=HWP
%     for qwp=QWP

files = dir('*.txt');
for file = files'
    clear hwp qwp;
    if file.name(2)=='_' && file.name(4)=='_'       %X_X_
        hwp=str2num(file.name(1));
        qwp=str2num(file.name(3));
    elseif file.name(2)=='_' && file.name(5)=='_'   %X_XX_
        hwp=str2num(file.name(1));
        qwp=str2num([file.name(3) file.name(4)]);
    elseif file.name(2)=='_' && file.name(6)=='_'   %X_XXX_
        hwp=str2num(file.name(1));
        qwp=str2num([file.name(3) file.name(4) file.name(5)]);
    elseif file.name(3)=='_' && file.name(5)=='_'   %XX_X_
        hwp=str2num([file.name(1) file.name(2)]);
        qwp=str2num(file.name(4));
    elseif file.name(3)=='_' && file.name(6)=='_'   %XX_XX_
        hwp=str2num([file.name(1) file.name(2)]);
        qwp=str2num([file.name(4) file.name(5)]);
    elseif file.name(3)=='_' && file.name(7)=='_'   %XX_XXX_
        hwp=str2num([file.name(1) file.name(2)]);
        qwp=str2num([file.name(4) file.name(5) file.name(6)]);
    elseif file.name(4)=='_' && file.name(6)=='_'   %XXX_X_
        hwp=str2num([file.name(1) file.name(2) file.name(3)]);
        qwp=str2num(file.name(5));
    elseif file.name(4)=='_' && file.name(7)=='_'   %XXX_XX_
        hwp=str2num([file.name(1) file.name(2) file.name(3)]);
        qwp=str2num([file.name(5) file.name(6)]);
    elseif file.name(4)=='_' && file.name(8)=='_'   %XXX_XXX_
        hwp=str2num([file.name(1) file.name(2) file.name(3)]);
        qwp=str2num([file.name(5) file.name(6) file.name(7)]);
    else
        display(['ERROR! The file `' file.name '´ do NOT match the format expected'])
    end
    

clear B1 aux1_i aux1_j B2 aux2_i aux2_j; close all
%clear all;close all
%k=2;
N1=256/st0;%N2=255;
A_max=1;%just for ploting
%usefull values and SemiCercle
phi1_0=0;%37*pi/180;  %    <-Rotate
% phi2_0=28*pi/180;  %    <-Rotate
A0_SC=1; %         %  <-Trim

mapa1=load([num2str(hwp) '_' num2str(qwp) '_curve2.txt']);
ph1=mapa1(1:N1,3);A1_0=mapa1(1:N1,2);
% mapa2=load(['valors_p' num2str(2) '.txt']);
% ph2=mapa2(1:N2,3);A2_0=mapa2(1:N2,2);

% figure
% polar(ph1,A1_0)

% A1_0=sqrt(A1_0);
A1=A1_0.*exp(1i*ph1);
% A2=A2_0.*exp(1i*ph2);

% B=zeros(1,(N*(N-1))/2);
% aux_i=zeros(1,(N*(N-1))/2);
% aux_j=zeros(1,(N*(N-1))/2);

count=1;
for i=1:st:N1
    for j=i:st:N1
        aux=(A1(i)+A1(j))/2;
        if abs(aux)<=A0_SC
            B1(count)=aux;
            aux1_i(count)=i;
            aux1_j(count)=j;
            count=count+1;
        end
    end
end
% count=1;
% for i=1:N2
%     for j=i:N2
%         aux=(A2(i)+A2(j))/2;
%         if abs(aux)<=A0_SC
%             B2(count)=aux;
%             aux2_i(count)=i;
%             aux2_j(count)=j;
%             count=count+1;
%         end
%     end
% end

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
A1_SC=A0_SC*ones(size(phi1_SC));
% phi2_SC=phi2_0:0.05:pi+phi1_0;
% A2_SC=A0_SC*ones(size(phi2_SC));

A_bottom=A0_SC:-0.01:-A0_SC;
phi1_bottom=(phi1_0+angle(A_bottom)).*ones(size(A_bottom));
% phi2_bottom=(phi2_0+angle(A_bottom)).*ones(size(A_bottom));

%ploting SLM1
h=figure;
polar(angle(A1),abs(A1),'xk') %Real curve
hold on
polar(angle(B1),abs(B1),'+b') %Possibles values
%polar(angle(C),abs(C),'xg') %Possibles values
polar(phi1_SC,A1_SC,'-k') %Semicircle
polar(pi+phi1_SC,A1_SC,'-k') %Semicircle
polar(phi1_bottom,abs(A_bottom),'-k'); %bottom SC
plot(0,0,'*k','markersize',10) %Center
title (['\lambda/2=' num2str(hwp) 'º ; \lambda/4=' num2str(qwp) 'º'])
% legend 'SLM curve' 'Possible coding values' 'Usefull zone'
hold off
print(h,'-dpng',[num2str(hwp) '_' num2str(qwp) '_polar_coded2.png']);
close(h)

end %for files

%     end %for QWPs
% end %for HWP

cd ..
cd ..
% figure
% polar(angle(C),abs(C),'x')
% 
% %Creating data for SLM1
% rho1=abs(B1)/A0_SC;phi1=mod(angle(B1)-phi1_0,2*pi);
% data=[rho1;phi1;aux1_i;aux1_j];
% fid = fopen(['ComplexValues_SLM' num2str(1) '.txt'],'wt');
% fprintf(fid,'%10.20f  %10.20f  %3.0f  %3.0f\n',data);
% fclose(fid);



% 
% %ploting SLM2
% h=figure;
% polar(angle(A2),abs(A2),'-r') %Real curve
% hold on
% polar(angle(B2),abs(B2),'+b') %Possibles values
% %polar(angle(C),abs(C),'xg') %Possibles values
% polar(phi2_SC,A2_SC,'-k') %Semicircle
% polar(pi+phi2_SC,A2_SC,'-k') %Semicircle
% polar(phi2_bottom,abs(A_bottom),'-k'); %bottom SC
% plot(0,0,'*k','markersize',10) %Center
% title (['SLM ' num2str(2)])
% legend 'SLM curve' 'Possible coding values' 'Usefull zone'
% hold off
% cd plots
% print(h,'-depsc','polar_SLM2.eps');
% cd ..
% % figure
% % polar(angle(C),abs(C),'x')
% 
% %Creating data for SLM2
% rho2=abs(B2)/A0_SC;phi2=mod(angle(B2)-phi2_0,2*pi);
% data=[rho2;phi2;aux2_i;aux2_j];
% fid = fopen(['ComplexValues_SLM' num2str(2) '.txt'],'wt');
% fprintf(fid,'%10.20f  %10.20f  %3.0f  %3.0f\n',data);
% fclose(fid);
