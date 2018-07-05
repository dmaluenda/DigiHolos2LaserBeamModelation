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


if ~exist('SLM_size','var'), SLM_size = [1024 768]; end
if ~exist('beam_type','var'), beam_type = 45; end
if ~exist('infile','var'), infile = 0; end
if ~exist('draw','var'), draw = 1; end 
if ~exist('stokes','var'), stokes = 0; end
if ~exist('gauss_correction','var'), gauss_correction = 0; end
if ~exist('prefix','var'), prefix = clock; end	

[y,x] = meshgrid( SLM_size(2)/2:-1:-SLM_size(2)/2+1, ...
	             -SLM_size(1)/2+1:SLM_size(1)/2);
theta = mod(atan2(y,x),2*pi); %[0 2pi]
rho   = sqrt(x.^2+y.^2);

E_x = zeros(SLM_size);
E_y = E_x ; Ph_x = E_x ; Ph_y = E_x ;

R = 200;%/sqrt(2);
%f0 = 2;
NA = 0.85;
rho = sqrt((x.^2+y.^2)/R^2);
phi = theta;
theta2 = asin(rho);

%g = exp(-(rho/f0).^2);        

E_x = -1i*cos(theta2).*sin(phi)+cos(phi);%.*g;
E_y =  1i*cos(theta2).*cos(phi)+sin(phi);%.*g;

pp = rho>NA;
E_x(pp) = 0;
E_y(pp) = 0;

Ph_x = angle(E_x);
Ph_y = angle(E_y);
E_x = abs(E_x);
E_y = abs(E_y); 

data1 = load('SupportFiles/Trans_SLM1.txt');
data2 = load('SupportFiles/Trans_SLM2.txt');

T1 = data1(:,1);
T2 = data2(:,1);

N = size(E_x);
holo1 = zeros(N);
holo2 = zeros(N);

for ii=1:N(1)
    for jj=1:N(2)
        [~,holo1(ii,jj)] = min(abs(E_x(ii,jj)-T1));
        [~,holo2(ii,jj)] = min(abs(E_y(ii,jj)-T2));
    end
end

[SLM1,SLM2] = scripts.rotate_SLM(holo1/255,holo2/255);
        
imwrite(SLM1', 'Holograms/' prefix '_holo1.bmp');
imwrite(SLM2', 'Holograms/' prefix '_holo2.bmp');  
