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



%Complex values. (2012)
%David Maluenda Niubó - Applied Physics and Optics (UB)
close all

N1=255;%max useful GrayLevel
N2=255;
phi1_0=0*pi/180;  %    <-Rotate
phi2_0=45*pi/180;  %    <-Rotate
A0_SC=0.32; %         %  <-Trim



if ~exist('LV','var'), LV = 0; end

if LV==0
    Ks=2;
    for k=Ks
        if k==1
            A1_0=A1;ph1=(phi1_45+phi1_135)/2;
        elseif k==2
            A2_0=A2;ph2=(phi2_45+phi2_135)/2;
        end
    end
else 
    Ks=1:2;
end

%Generates a file with the LCD's response.
for k=Ks
%     if k==1
%         phi_def=ph1*pi/180;
%         T_def=A1_0/max(A1_0);
%     else
%         phi_def=ph2*pi/180;
%         T_def=A2_0/max(A2_0);
%     end
%     
% %Generates the response of every GrayLevel (256x3)
%     curve=[(1:256)' T_def' phi_def'];
%     fid = fopen(['valors_p' num2str(k) '.txt'],'wt');
%     fprintf(fid,'%3.0f %6.4f% 6.4f\n',curve');
%     fclose(fid);
    
%Generates a phase map, ie. GrayLevel VS Phase (6285x2)
    mapa=load(['valors_p' num2str(k) '.txt']);
    
    if k==1
        A1_0=mapa(:,2);
        ph1=mapa(:,3)*180/pi;
    elseif k==2
        A2_0=mapa(:,2);
        ph2=mapa(:,3)*180/pi;
    end
    
    mapa=mapa(:,3);
    mapa=round(mapa*1000);
    gl=zeros(1,6285);
    for i=1:255
        gl(mapa(i)+1:mapa(i+1)+1)=i;
    end

    gl(max(mapa):round((6285+max(mapa))/2))=gl(max(mapa));
    data=[(0:0.001:2*pi+0.001);gl];
    fid = fopen(['curve_SLM' num2str(k) '.txt'],'wt');
    fprintf(fid,'%6.4f  %3.0f\n',data);
    fclose(fid);

if k==1
M1=A1_0.*exp(1i*ph1*pi/180);
elseif k==2
M1=A2_0.*exp(1i*ph2*pi/180);
end

C1=zeros(1,(N1*(N1-1))/2);
aux1_i=zeros(1,(N1*(N1-1))/2);
aux1_j=zeros(1,(N1*(N1-1))/2);
% C2=zeros(1,(N2*(N2-1))/2);
% aux2_i=zeros(1,(N2*(N2-1))/2);
% aux2_j=zeros(1,(N2*(N2-1))/2);

count=1;
for i=1:N1
    for j=i:N1
        aux=(M1(i)+M1(j))/2;
        if abs(aux)<=A0_SC
            C1(count)=aux;
            aux1_i(count)=i;
            aux1_j(count)=j;
            count=count+1;
        end
    end
end
C1=resizer(C1);   %fits in non-zero array
aux1_i=resizer(aux1_i);
aux1_j=resizer(aux1_j);
% count=1;
% for i=1:N2
%     for j=i:N2
%         aux=(M2(i)+M2(j))/2;
%         if abs(aux)<=A0_SC
%             C2(count)=aux;
%             aux2_i(count)=i;
%             aux2_j(count)=j;
%             count=count+1;
%         end
%     end
% end
% C2=resizer(C2);   %fits in non-zero array
% aux2_i=resizer(aux2_i);
% aux2_j=resizer(aux2_j);

date=clock;

%Creating data
rho1=abs(C1)/A0_SC;phi1=mod(angle(C1)-phi1_0,2*pi);
data=[rho1;phi1;aux1_i;aux1_j];
fid = fopen(['ComplexValues_SLM' num2str(k) '.txt'],'wt');
fprintf(fid,'%10.20f  %10.20f  %3.0f  %3.0f\n',data);
fclose(fid);

%Ploting data
phi1_SC=phi1_0:0.05:pi+phi1_0;
A1_SC=A0_SC*ones(size(phi1_SC));
A_bottom=A0_SC:-0.01:-A0_SC;
phi1_bottom=(phi1_0+angle(A_bottom)).*ones(size(A_bottom));

h=figure;
polar(angle(M1),abs(M1),'xk') %Real curve
hold on
kkkkk=polar(angle(C1),abs(C1),'+b'); %Possibles values
set(kkkkk,'markersize',0.5);
polar(phi1_SC,A1_SC,'-k') %Semicircle
polar(pi+phi1_SC,A1_SC,'-k') %Semicircle
polar(phi1_bottom,abs(A_bottom),'-k'); %bottom SC
plot(0,0,'*k','markersize',10) %Center
legend 'SLM curve' 'Possible coding values'
hold off
print(h,'-dpng',['plots\' num2str(date(3)) num2str(date(2)) ...
    num2str(date(1)-2000) 'polar.png']);
%close(h)

end


% st=10;
% for slm=1:2
% if slm==1
%     M=M1;phi_0=phi1_0;N=N1;
% else
%     M=M2;phi_0=phi2_0;N=N2;
% end
% C=zeros(1,1:round((N/st*2)^4/2));   %initializing
% aux_i=C;aux_j=C;aux_k=C;aux_l=C;
% 
% 
% %--------------SLM1-------------------------
% %---in-pairs--------------------------------
% count=1;
% for i=1:N
%     for j=i:N
%         aux=(M(i)+M(j))/2;
%         alpha=mod(angle(aux),2*pi);         
%         if (abs(aux)<=A0_SC && alpha>phi_0 && alpha<phi_0+pi)
%             C(count)=aux;
%             aux_i(count)=i;
%             aux_j(count)=j;
%             aux_k(count)=i;
%             aux_l(count)=j;
%             count=count+1;
%         end
%     end
% end
% clear i j k l
% %---4-single-------------------------------
% for i=1:st:N
%     st1=randi([1 st],1);
%     for j=i:st1:N
%         st2=randi([1 2*st],1);
%         for k=j:st2:N
%             st3=randi([1 2*st],1);
%             for l=k:st3:N
%                 aux=(M(i)+M(j)+M(k)+M(l))/4;
%                 alpha=mod(angle(aux),2*pi);
%                 if (abs(aux)<=A0_SC && (alpha<phi_0 || alpha>phi_0+pi))
%                     C(count)=aux;
%                     aux_i(count)=i;
%                     aux_j(count)=j;
%                     aux_k(count)=k;
%                     aux_l(count)=l;
%                     count=count+1;
%                 end
%             end
%         end
%     end
% end
% C=resizer(C);   %fits in non-zero array
% aux_i=resizer(aux_i);
% aux_j=resizer(aux_j);
% aux_k=resizer(aux_k);
% aux_l=resizer(aux_l);
% 
% %---Writing-data-in-file------------------------------------------
% factorUINT=10000;
% rho=uint16(abs(C)/A0_SC*factorUINT);phi=uint16(mod(angle(C)-phi_0,2*pi)*factorUINT);
% data=[rho;phi;aux_i;aux_j;aux_k;aux_l];
% fid = fopen(['ComplexValues4_' num2str(st) '_SLM' num2str(slm) '.txt'],'wt');
% fprintf(fid,'%hu  %hu  %hu  %hu  %hu  %hu\n',data);
% fclose(fid);
% 
% clear M C aux_i aux_j aux_k aux_l;
% end