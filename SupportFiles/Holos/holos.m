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
aux=1;

for SLM_number=1:2
for Trans=0:Tstep:100
    for Phase=0:PHstep:PHf
        [SLM]= mapa_Holo3(Trans/100,Phase,N,SLM_number);%Holo=f(T,ph[º],N,SLM)
        imwrite(SLM,['SLM' num2str(SLM_number) '_T' num2str(Trans) '_ph' num2str(Phase) '_uint8.png'],'PNG');
    end
    if mod(Trans,Cstep)<=aux
        disp(['Holo for display ' num2str(SLM_number) ': ' num2str(Trans) '% done']);
        aux=mod(Trans,Cstep);
    end
end

end
% imshow(SLM1)