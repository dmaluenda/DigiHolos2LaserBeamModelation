%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Map_generator.m (2012)
% David Maluenda NiubÛ - Applied Physics and Optics (UB)
% phi2=phi2-min(phi2);
% phi2(gl_max+1:256)=phi2(gl_max);%gl_max is the gray level whose phase is
% the maximum one
clear all
%%

DATA=dlmread('SLMresponse.txt','',3,0);
I_1=DATA(:,2);PHI1=DATA(:,3);GLM1=find(I_1,1,'last')-1;
I_2=DATA(:,4);PHI2=DATA(:,5);GLM2=find(I_2,1,'last')-1;



%%

SLM=1:2; %SLM's label that is caracterazing 

for k=SLM
    if k==1
        phi_def=PHI1;
        T_def=I_1/max(I_1);
        GLmax=GLM1;
    else
        phi_def=PHI2;
        T_def=I_2/max(I_2);
        GLmax=GLM2;
    end
    curve=[(1:GLmax)' T_def(1:GLmax)' phi_def(1:GLmax)'];

    fid = fopen(['valors_p' num2str(k) '.txt'],'wt');
    fprintf(fid,'%3.0f %6.4f% 6.4f\n',curve');
    fclose(fid);
    figure
    polar(phi_def(1:GLmax), T_def(1:GLmax))
    title (['SLM ' num2str(k)])

    mapa=round(phi_def*1000);
    gl=zeros(1,6285);
    for i=1:255
        gl(mapa(i)+1:mapa(i+1)+1)=i;
    end

    gl(max(mapa):round((6285+max(mapa))/2))=gl(max(mapa));
    data=[(0:0.001:2*pi+0.001);gl];
    figure
    plot(data(1,:),data(2,:))
    title (['SLM ' num2str(k)])
    fid = fopen(['curve_SLM' num2str(k) '.txt'],'wt');
    fprintf(fid,'%6.4f  %3.0f\n',data);
    fclose(fid);
end