%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a serial of holograms to encode T=0:Tstep:1 and Ph=0:PhStep:PhFinal
%  in a half of the image remaining the other side at T=1 and Ph=0 (Vertically)

N = [1024 768];
SLM_number = 2;

aux = 1;
for Trans=0:Tstep:100
    for Phase=0:PHstep:PHf
        [SLM2] = mapa_Holo2(Trans/100, Phase, N, SLM_number);
        imwrite(SLM2', ['SLM2_T' num2str(Trans) '_ph' num2str(Phase) '.png'], 'PNG');
    end
    if mod(Trans,Cstep)<=aux
        disp(['SemiHoloV: ' num2str(Trans) '% done']);
        aux=mod(Trans, Cstep);
    end
end
display('SLM2 is done')

% imshow(SLM2')
        