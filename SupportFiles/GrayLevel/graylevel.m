%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fills this directory with a plain images of every gray level [0:255]

im = ones(768, 1024);
for i=0:255
    im(:) = i;
    imwrite(im/255, [num2str(i) 'level.png'], 'PNG');
end
