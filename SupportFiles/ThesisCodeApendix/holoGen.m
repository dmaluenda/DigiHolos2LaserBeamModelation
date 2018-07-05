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



% HoloGenerator (2013)
% David Maluenda Niubo - Applied Physics and Optics (UB)
% 
% This function take Amplitudes and Phases and return the
% macro pixel holograms following the accessibleValues.dat
%
%       output                        input
% -------------------|--------------------------------------
% (  SLM1 ,  SLM2 ) = f( Amp_1 , Amp_2 ,  Phase1 ,  Phase2 )
% ( [0 1] , [0 1] ) = f( [0 1] , [0 1] , [0 2pi] , [0 2pi] )
%       dynamic range of in- and out-put matrices

function[SLM1,SLM2]=holoGen(Amp1,Amp2,Phase1,Phase2)
clear all; close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FOR TESTING THE FUNTION USE THE FOLLOWING BEAM DESIGN
%
% LG01 beam radially (azimuthally) polarized for the inner (outer) part 
%
%
%         [y,x] = meshgrid( linspace( -768/2+1, 768/2, 768) ,...
%                           linspace(-1024/2+1,1024/2,1024)  );
%         theta = mod( atan2(y,x) , 2*pi ); %[0 2pi]
%         rho   = sqrt(x.^2+y.^2);
% 
%         % Radial part
%         R_x   = (abs(cos(theta)));
%         R_y   = (abs(sin(theta)));
%         R_Phx = angle(cos(theta));
%         R_Phy = angle(sin(theta));
%         % Azimuthal part
%         A_x   = (abs(sin(theta)));
%         A_y   = (abs(cos(theta)));
%         A_Phx = angle(sin(theta));
%         A_Phy = angle(cos(theta))+pi;
% 
%         R    = 70;
%         prof = abs((-2*rho.^2/R^2+1).*exp(-rho.^2/R^2));
%         core = (rho.^2<50).*1;
%         clad = 1-core;
% 
%         Amp1   = prof.*( core.*R_x + clad.*A_x );
%         Amp2   = prof.*( core.*R_y + clad.*A_y );
%         Phase1 = core.*R_Phx + clad.*A_Phx;
%         Phase2 = core.*R_Phy + clad.*A_Phy;
% 
%         Phase1 = mod(Phase1,2*pi);
%         Phase2 = mod(Phase2,2*pi);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% initializing

% size of matrices
N = size(Phase1);

%zone of interest
X = [1  N(2)] ; Y = [1  N(1)];
%X= [221 570] ; Y = [317 740];

% desirable values from input matrices
C1 = Amp1.*exp(1i*Phase1); 
C2 = Amp2.*exp(1i*Phase2);

%% loading maps of accessible values values
data = dlmread('AccessValues.dat','',3,0); % All data together

Map1_1  = data(:,1); % index of M1 for SLM1
Map2_1  = data(:,2); % index of M2 for SLM2
T_SLM1  = data(:,3); % amplitude of coded value for SLM1 
ph_SLM1 = data(:,4); % phase of coded value for SLM1

Map1_2  = data(:,5); % index of M1
Map2_2  = data(:,6); % index of M2
T_SLM2  = data(:,7); % amplitude of coded value for SLM1 
ph_SLM2 = data(:,8); % phase of coded value for SLM1

% accesible complex values
C_SLM1  = T_SLM1.*exp(1i*ph_SLM1); 
C_SLM2  = T_SLM2.*exp(1i*ph_SLM2);



%% resizing desired modulation for macropixel procedure

C1_mean( 1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2 ) = ... 
     (   C1( Y(1)  :2:Y(2) , X(1)  :2:X(2) ) + ...
         C1( Y(1)+1:2:Y(2) , X(1)  :2:X(2) ) + ...
         C1( Y(1)  :2:Y(2) , X(1)+1:2:X(2) ) + ...
         C1( Y(1)+1:2:Y(2) , X(1)+1:2:X(2) )     )/4;

C2_mean( 1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2 ) = ... 
     (   C2( Y(1)  :2:Y(2) , X(1)  :2:X(2) ) + ...
         C2( Y(1)+1:2:Y(2) , X(1)  :2:X(2) ) + ...
         C2( Y(1)  :2:Y(2) , X(1)+1:2:X(2) ) + ...
         C2( Y(1)+1:2:Y(2) , X(1)+1:2:X(2) )     )/4;

%% Arizon's double phase cell holographic technique

% initialing the holograms
SLM1 = zeros(N);
SLM2 = zeros(N);


% seeking the minimal euclidian distance between
%       the desired value and coded values

m1 = ones(size(C1_mean)); % auxiliar values with the
m2 = ones(size(C2_mean)); % minimum euclidian distance
p1 = ones(size(C1_mean)); % index with the nearest
p2 = ones(size(C2_mean)); % accessible value to the desired
for i=1:(Y(2)-Y(1)+1)/2
   for j=1:(X(2)-X(1)+1)/2
        [m1(i,j),p1(i,j)] = min(abs( C1_mean(i,j)-C_SLM1 ));
        [m2(i,j),p2(i,j)] = min(abs( C2_mean(i,j)-C_SLM2 ));
   end
end


% filling the SLM1 with the indices p1 for the main diagonal
%         and p2 for the other diagonal

% the upper left value in the macropixel
SLM1( Y(1)  :2:Y(2) , X(1)  :2:X(2) ) = ...  
    Map1_1( p1(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) );
% the upper right value in the macropixel
SLM1( Y(1)+1:2:Y(2) , X(1)+1:2:X(2)) = ...
    Map1_1( p1(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) );
% the downer left value in the macropixel
SLM1( Y(1)+1:2:Y(2) , X(1)  :2:X(2)) = ...
    Map2_1( p1(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) );
% the downer right value in the macropixel
SLM1( Y(1)  :2:Y(2) , X(1)+1:2:X(2))= ...
    Map2_1( p1(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) );


% filling the SLM2 with the indices p1 for the main diagonal
%         and p2 for the other diagonal

% the upper left value in the macropixel
SLM2( Y(1)  :2:Y(2) , X(1)  :2:X(2) ) = ...  
    Map1_2( p2(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) );
% the upper right value in the macropixel
SLM2( Y(1)+1:2:Y(2) , X(1)+1:2:X(2)) = ...
    Map1_2( p2(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) );
% the downer left value in the macropixel
SLM2( Y(1)+1:2:Y(2) , X(1)  :2:X(2)) = ...
    Map2_2( p2(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) );
% the downer right value in the macropixel
SLM2( Y(1)  :2:Y(2) , X(1)+1:2:X(2))= ...
    Map2_2( p2(1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2) ); 


% from double to integer image
SLM1 = SLM1/255;
SLM2 = SLM2/255;

% plotting the holograms
figure;imagesc(SLM1',[0 1]);colormap('gray');
axis equal;axis([0 N(1) 0 N(2)]);title SLM1
figure;imagesc(SLM2',[0 1]);colormap('gray');
axis equal;axis([0 N(1) 0 N(2)]);title SLM2

imwrite(SLM1','Hologram1.bmp')
imwrite(SLM2','Hologram2.bmp')
