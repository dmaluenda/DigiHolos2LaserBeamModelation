%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fills this directory with a semi-plain images of every gray level [0:255]
%  semi-plain is understood as half-image with a certain value and
%  the other with 0value 

im=zeros(768, 1024);
for i=0:255
    im(1:1024*768/2) = i;
    imwrite(im/255, [num2str(i) 'semilevel.png'], 'PNG');
end
