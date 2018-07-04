%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns a (square) 90º rotated for SLM_1
%    --->   SLM1 is rotated in the setup   <---
function[im_out]=rotate(im_in)
    N = size(im_in);

    % Padding to make an square
    if N(1)==N(2) 
        M = N(1);
    elseif N(1)>N(2)
        M = N(1);
        im_in(1:N(1), N(2):N(1)) = 1;
    else
        M = N(2);
        im_in(N(1):N(2), 1:N(2)) = 1;
    end

    % Rotation matrix
    AD = zeros(M,M);
    for i=1:M
        for j=1:M
            if i+j==M
                AD(i,j) = 1;
            end
        end
    end

    im_out = (im_in*AD)'*AD;
    