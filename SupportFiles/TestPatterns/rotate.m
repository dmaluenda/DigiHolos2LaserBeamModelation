%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[imout] = rotate(imin)
    N = size(imin);
    imin = normalize_2D(imin);

    if N(1)>N(2)
        imin(1:N(1),N(2):N(1)) = 0;
        AD = zeros(N(1),N(1));
        M = N(1);
    else
        imin(N(1):N(2),1:N(2)) = 0;
        AD = zeros(N(2),N(2));
        M = N(2);
    end

    for i=1:M
        for j=1:M
            if i+j==M
                AD(i,j) = 1;
            end
        end
    end



    imout = (imin*AD)'*AD;