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


% HoloGenerator. PROGRAM (2012)
% David Maluenda Niubo - Applied Physics and Optics (UB)
% 
% This script take Amplitudes and Phases and return the
% macro pixel holograms following the accessibleValues.dat
%       output                      input
% (  SLM1 ,  SLM2 )=f( Amp_1 , Amp_2 ,  Phase1 ,  Phase2 )
%       dynamic range of in- and out-put matrices
% ( [0 1] , [0 1] )=f( [0 1] , [0 1] , [0 2pi] , [0 2pi] ) 
function[SLM1,SLM2]=HoloGen(Amp1,Amp2,Phase1,Phase2)

% size of matrices
N = size(Phase1);

% loading maps of accessible values values
data = load('accessValues.dat'); % All data together

Map1_1  = data(:,1); % index of M1
Map2_1  = data(:,2); % index of M2
T_SLM1  = data(:,3); % amplitude of  SLM1 
ph_SLM1 = data(:,4); % 

Map1_2  = data(:,5); % index of M1
Map2_2  = data(:,6); % index of M2
T_SLM2  = data(:,7); % modulation of SLM2
ph_SLM2 = data(:,8);


% accesible complex values
C_SLM1  = T_SLM1.*exp(1i*ph_SLM1); 
C_SLM2  = T_SLM2.*exp(1i*ph_SLM2);


% desirable values from input matrices
C1 = Amp1.*exp(1i*Phase1); 
C2 = Amp2.*exp(1i*Phase2);

% resizing for macropixel procedure
C1_mean( 1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2 ) = ... 
      (  C1( Y(1)  :2:Y(2) , X(1)  :2:X(2) ) + ...
         C1( Y(1)+1:2:Y(2) , X(1)  :2:X(2) ) + ...
         C1( Y(1)  :2:Y(2) , X(1)+1:2:X(2) ) + ...
         C1( Y(1)+1:2:Y(2) , X(1)+1:2:X(2) )  )/4;

C2_mean( 1:(Y(2)-Y(1)+1)/2 , 1:(X(2)-X(1)+1)/2 ) = ... 
      (  C2( Y(1)  :2:Y(2) , X(1)  :2:X(2) ) + ...
         C2( Y(1)+1:2:Y(2) , X(1)  :2:X(2) ) + ...
         C2( Y(1)  :2:Y(2) , X(1)+1:2:X(2) ) + ...
         C2( Y(1)+1:2:Y(2) , X(1)+1:2:X(2) )  )/4;

% initialing the holograms
SLM1 = zeros(N);
SLM2 = zeros(N);

m1 = ones(size(C1_mean)); % auxiliar values with the
m2 = ones(size(C2_mean)); % minimum euclidian distance
p1 = ones(size(C1_mean)); % index with the nearest
p2 = ones(size(C2_mean)); % accessible value to the desired
aux=0;

for i=1:(Y(2)-Y(1)+1)/2
   for j=1:(X(2)-X(1)+1)/2
        [m1(i,j),p1(i,j)]=min(abs(C1_mean(i,j)-C_SLM1));
        [m2(i,j),p2(i,j)]=min(abs(C2_mean(i,j)-C_SLM2));
   end
end

% filling the SLM1 with the indices p1 for the main diagonal
%   and p2 for the other diagonal
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
%   and p2 for the other diagonal
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
SLM1=(SLM1)/255;
SLM2=(SLM2)/255;

