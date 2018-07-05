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

function[DR,ratio] = dynamic_range(imageIN)

N = size(imageIN);

count0   = 0;
count255 =0;

DR = [min(min(imageIN)) max(max(imageIN))];

for i=1:N(1)
    for j=1:N(2)
        if imageIN(i,j) <= DR(1)
            count0 = count0+1;
        elseif imageIN(i,j)>=DR(2)
            count255 = count255+1;
        end
    end
end

ratio0   = count0/N(1)/N(2);
ratio255 = count255/N(1)/N(2);

ratio = [ratio0 ratio255];