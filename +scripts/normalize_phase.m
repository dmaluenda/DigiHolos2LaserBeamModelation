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

%Return the phase shift between the arms of the interferometer
function[phase]=normalize_phase(interference1, interference2, intensity)
	N = size(interference1);

	phase = zeros(N);
	sign  = ones(N);

	data = load('response.txt');
	mapa = data(:,1);
	for j=1:N(2)
	    for i=1:N(1)
	        if (interference1(i,j)==0 || intensity(i,j)==0)
	            phase(i,j) = 0;
	        elseif (interference1(i,j) > intensity(i,j))
	            phase(i,j) = 1;
	        else
	            phase(i,j) = interference1(i,j)/intensity(i,j);
	%             if interference2(i,j)/intensity(i,j)>0.5
	%                 sign(i,j)=-1;
	%             end
	        end
	    end
	end

	phase(:) = mapa(round(phase(:)*255+1));
	phase = phase.*sign;
