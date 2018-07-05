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


% phi2=phi2-min(phi2);
% phi2(gl_max+1:256)=phi2(gl_max);%gl_max is the gray level whose phase is
% the maximum one

SLM=2; %SLM's label that is caracterazing 


for j=SLM
    if j==1
        I_1=A1;
        %PHI1_45=phi1_45;PHI1_135=phi1_45;
        %PHI1_45=PHI1_45-min(PHI1_45);
        %PHI1_135=PHI1_135-min(PHI1_135);
        PHI1=phi1;%mod( (PHI1_45+PHI1_135)/2 ,360);
        %PHI1=PHI1-min(PHI1);
        phi_def=PHI1*pi/180;
        T_def=I_1/max(I_1);
        GLmax=256;
    else
        I_2=A2;
        %PHI2_45=phi2_45;PHI2_135=phi2_45;
        %PHI2_45=PHI2_45-min(PHI2_45);
        %PHI2_135=PHI2_135-min(PHI2_135);
        PHI2=phi2;%mod( (PHI2_45+PHI2_135)/2 ,360);
        %PHI2=PHI2-min(PHI2);
        phi_def=PHI2*pi/180;
        T_def=I_2/max(I_2);
        GLmax=256;
    end
    
    if size( T_def ,1)<size( T_def ,2), T_def = T_def.' ;end
    if size(phi_def,1)<size(phi_def,2),phi_def=phi_def.';end
        
    curve=[(1:GLmax).' T_def(1:GLmax) phi_def(1:GLmax)];

    fid = fopen(['valors_p' num2str(j) '.txt'],'wt');
    fprintf(fid,'%3.0f %6.4f% 6.4f\n',curve');
    fclose(fid);
    figure
    polar(phi_def(1:GLmax), T_def(1:GLmax))
    title (['SLM ' num2str(j)])

    mapa=round(phi_def*1000);
    gl=zeros(1,6285);
    for i=1:255
        gl(mapa(i)+1:mapa(i+1)+1)=i;
    end

    gl(max(mapa):round((6285+max(mapa))/2))=gl(max(mapa));
    data=[(0:0.001:2*pi+0.001);gl];
    figure
    plot(data(1,:),data(2,:))
    title (['SLM ' num2str(j)])
    fid = fopen(['curve_SLM' num2str(j) '.txt'],'wt');
    fprintf(fid,'%6.4f  %3.0f\n',data);
    fclose(fid);
end