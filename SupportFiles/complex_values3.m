%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Generates the posible complex values achievables by SLMs using
%the macropixel Arizon's codification procedure. (2013)
%David Maluenda Niubó - Applied Physics and Optics (UB)

clear all;%B1 aux1_i aux1_j aux1_k aux1_l B2 aux2_i aux2_j aux2_k aux2_l;
close all
tic

%--Parameters-of-usefull-values-and-SemiCercle--
N1=256;N2=256; %Usable values from SLMs
phi1_0=37*pi/180;%  <-Rotate (initial phase in the output file)
phi2_0=28*pi/180;%  <-Rotate (initial phase in the output file)
A0_SC=1; %          <-Trim (Normalizing factor in the output file)
A_max=1;%just for ploting

%---getting-data-------->  complex(Nx1): [ A1 ; A2 ]
mapa1=load(['valors_p' num2str(1) '.txt']);
ph1=mapa1(1:N1,3);A1_0=mapa1(1:N1,2);
A1=A1_0.*exp(1i*ph1);
mapa2=load(['valors_p' num2str(2) '.txt']);
ph2=mapa2(1:N2,3);A2_0=mapa2(1:N2,2);
A2=A2_0.*exp(1i*ph2);


%---step to scan the values
st=2;

for slm=1:2

if slm==1
    A=A1;
    phi_0=phi1_0;
else
    A=A2;
    phi_0=phi2_0;
end
B1(1,1:round((N1/st*2)^3/2))=0;   %initializing
aux1_i=B1;aux1_j=B1;aux1_k=B1;aux1_l=B1;


%--------------SLM1-------------------------
%---in-pairs--------------------------------
count=1;
% for i=1:N1
%     for j=i:N1
%         aux=(A(i)+A(j))/2;
%         alpha=mod(angle(aux),2*pi);         
%         if (abs(aux)<=A0_SC && alpha>phi_0 && alpha<phi_0+pi)
%             B1(count)=aux;
%             aux1_i(count)=i;
%             aux1_j(count)=j;
%             aux1_k(count)=i;
%             aux1_l(count)=j;
%             count=count+1;
%         end
%     end
% end
% clear i j k l
%---4-single-------------------------------
for i=1:st:N1
    st1=randi([1 st],1);
    for j=i:st1:N1
        st2=randi([1 2*st],1);
        for k=j:st2:N1
%             st3=randi([1 2*st],1);
%             for l=k:st3:N1
                aux=(A(i)+A(j)+A(k))/3;
                alpha=mod(angle(aux),2*pi);
%                 if (abs(aux)<=A0_SC && (alpha<phi_0 || alpha>phi_0+pi))
                    B1(count)=aux;
                    aux1_i(count)=i;
                    aux1_j(count)=j;
                    aux1_k(count)=k;
%                     aux1_l(count)=l;
                    count=count+1;
%                 end
%             end
        end
    end
end
B1=resizer(B1);   %fit in non-zero array
aux1_i=resizer(aux1_i);
aux1_j=resizer(aux1_j);
aux1_k=resizer(aux1_k);
% aux1_l=resizer(aux1_l);

%---Ploting-data--------------------------------------------------
%------Cercles-to-plot--------------------------------------------
phi1_SC=phi1_0:0.05:pi+phi1_0;
A1_SC=A0_SC*ones(size(phi1_SC));
phi2_SC=phi2_0:0.05:pi+phi1_0;
A2_SC=A0_SC*ones(size(phi2_SC));
A_bottom=A0_SC:-0.01:-A0_SC;
phi1_bottom=(phi1_0+angle(A_bottom)).*ones(size(A_bottom));
phi2_bottom=(phi2_0+angle(A_bottom)).*ones(size(A_bottom));


if slm==1

%---Writing-data-in-file------------------------------------------
%------SLM1-------------------------------------------------------
rho1=uint16(abs(B1)/A0_SC*10000);phi1=uint16(mod(angle(B1)-phi1_0,2*pi)*10000);
data=[rho1;phi1;aux1_i;aux1_j;aux1_k];
fid = fopen(['ComplexValues3_' num2str(st) '_SLM' num2str(1) '.txt'],'wt');
fprintf(fid,'%hu  %hu  %hu  %hu  %hu\n',data);
fclose(fid);

%---ploting
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
print(h,'-depsc','polar4_SLM1.eps');
cd ..

else


%----write----------------------------------------------------
factorUINT=10000;
rho1=uint16(abs(B1)/A0_SC*factorUINT);phi1=uint16(mod(angle(B1)-phi2_0,2*pi)*factorUINT);
data=[rho1;phi1;aux1_i;aux1_j;aux1_k];
fid = fopen(['ComplexValues3_' num2str(st) '_SLM' num2str(2) '.txt'],'wt');
fprintf(fid,'%hu  %hu  %hu  %hu  %hu\n',data);
fclose(fid);


%----plot----------------------------------------------------
h=figure;
polar(angle(A2),abs(A2),'-r') %Real curve
hold on
polar(angle(B1),abs(B1),'+b') %Possibles values
%polar(angle(C),abs(C),'xg') %Possibles values
polar(phi2_SC,A2_SC,'-k') %Semicircle
polar(pi+phi2_SC,A2_SC,'-k') %Semicircle
polar(phi2_bottom,abs(A_bottom),'-k'); %bottom SC
plot(0,0,'*k','markersize',10) %Center
title (['SLM ' num2str(2)])
legend 'SLM curve' 'Possible coding values' 'Usefull zone'
hold off
cd plots
print(h,'-depsc','polar4_SLM2.eps');
cd ..
% figure
% polar(angle(C),abs(C),'x')

end
%clear B1 aux1_i aux1_j aux1_k aux1_l;

end

toc
