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


fid = fopen('+scripts\design_kinds.txt');

tline = fgetl(fid);
i=1;
while ischar(tline)
    disp(tline)
    names[i]=tline;
    tline = fgetl(fid);
    i=i+1;
end

fclose(fid);