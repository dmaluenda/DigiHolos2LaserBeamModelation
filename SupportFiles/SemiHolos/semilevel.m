%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=[768 1024];
SLM_number=1;

aux=1;
for Trans=0:Tstep:100
    for Phase=0:PHstep:PHf
        [SLM1]= mapa_Holo2(Trans/100,Phase,N,SLM_number);
        imwrite(SLM1,['SLM1_T' num2str(Trans) '_ph' num2str(Phase) '.png'],'PNG');
    end
    if mod(Trans,Cstep)<=aux
        disp(['SemiHolo: ' num2str(Trans) '% done']);
        aux=mod(Trans,Cstep);
    end
end
%display('SLM1 is done')
% imshow(SLM1)