%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% mapGen.m (2013)
% David Maluenda Niubó - Applied Physics and Optics (UB)
%
% Generates the 'AccessValues.dat' file with the asseccible values
% by means of Arizon's Macropixel Double Cell Hologram technique
% of a couple of modulator with responses saved in 'SLMresponse.dat' file
clear all;close all

%% Fixing some parameters

phi1_0 = 0;   % Fixing the relative phase  <- Rotate SLM1
phi2_0 = 0;   % Fixing the relative phase  <- Rotate SLM1
A1_max = 0.3; % renormalization of SLM2    <- Trim SLM1
A2_max = 0.3; % renormalization of SLM2    <- Trim SLM2
A_max  = 1;   % just for ploting


%% Loading the response of modulators

DATA = dlmread('SLMresponse.dat','',3,0);
A1_0 = DATA(:,2) ; ph1=DATA(:,3) ; N1=find(A1_0,1,'last');
A2_0 = DATA(:,4) ; ph2=DATA(:,5) ; N2=find(A2_0,1,'last');

% figure;polar(ph1,A1_0)
% figure;polar(ph2,A2_0)

A1 = A1_0.*exp( 1i*ph1 ); % Complex modulation response of SLM1
A2 = A2_0.*exp( 1i*ph2 ); % Complex modulation response of SLM1


%% Exploring all the accessible values via Arizon's Cell technique

% SLM1
count=1;count2=1;
for i=1:N1
    for j=i:N1 
        aux = ( A1(i)+A1(j) )/2;  % accessible value
        if abs(aux)<=A1_max
            C1(count) = aux;
            aux1_i(count) = i;    % M1
            aux1_j(count) = j;    % M2
            count = count+1;
        else
            B1(count2) = aux;     % non-useful accessible value
            count2 = count2+1;
        end
    end
end

% SLM2
count=1;count2=1;
for i=1:N2
    for j=i:N2
        aux = ( A2(i)+A2(j) )/2;  % accessible value
        if abs(aux)<=A2_max
            C2(count) = aux;
            aux2_i(count) = i;    % M1
            aux2_j(count) = j;    % M2
            count = count+1;
        else
            B2(count2) = aux;     % non-useful accessible value
            count2 = count2+1;
        end
    end
end


%% Plotting data

% Circles of renormalization
phi1_SC = 0:0.05:2*pi;
A1_SC   = A1_max*ones(size(phi1_SC));
phi2_SC = 0:0.05:2*pi;
A2_SC   = A2_max*ones(size(phi2_SC));

% SLM1
figure;
h = polar( angle(A1) , abs(A1) , '.r' );    % Real curve
set(gca, 'fontsize',14 );
hold on

h = polar( angle(B1) , abs(B1) , '.b' );    % Accesible values
set(h, 'color',[0.5 0.5 0.5] , 'markersize',0.1 );

h = polar( angle(C1) , abs(C1) , '.b' );    % Renormalized values
set(h, 'color',[0 0 0.5] , 'markersize',0.1 );

h = polar( phi1_SC , A1_SC , '-.k' );       % Renormalization Circle
set(h, 'linewidth',2 );

plot(0,0, '*k' , 'markersize',10 )           % Center

title (['SLM ' num2str(1)])
legend({'Modulation response','Accessible values', ...
    'Renormalized values','Renormalizated zone'} , ...
    'FontSize',14 , 'FontWeight','bold' , 'location','southeast' )
hold off


% SLM2
figure;
h = polar( angle(A2) , abs(A2) , '.r' );    % Real curve
set(gca, 'fontsize',14 );
hold on

h = polar( angle(B2) , abs(B2) , '.b' );    % Accesible values
set(h, 'color',[0.5 0.5 0.5] , 'markersize',0.1 );

h = polar( angle(C2) , abs(C2) , '.b' );    % Renormalized values
set(h, 'color',[0 0 0.5] , 'markersize',0.1 );

h = polar( phi2_SC , A2_SC , '-.k' );       % Renormalization Circle
set(h, 'linewidth',2 );
plot(0,0, '*k' , 'markersize',10)           % Center

title (['SLM ' num2str(2)])
legend({'Modulation response','Accessible values',...
    'Renormalized values','Renormalizated zone'} ,...
    'FontSize',14 , 'FontWeight','bold' , 'location','southeast')
hold off


%% Saving data to AccessValues.dat file

Cx0 = abs(C1)/A1_max;                 % Amplitude of coded value for SLM1
Phx = mod( angle(C1)-phi1_0 , 2*pi ); % Phase of coded value for SLM1
Cy0 = abs(C2)/A2_max;                 % Amplitude of coded value for SLM2
Phy = mod( angle(C2)-phi2_0 , 2*pi ); % Phase of coded value for SLM2

% filling the holes with NaN
if max(size(Cx0))>max(size(Cy0))
    Cy0   ( end+1:max(size(Cx0)) ) = NaN;
    Phy   ( end+1:max(size(Cx0)) ) = NaN;
    aux2_i( end+1:max(size(Cx0)) ) = NaN;
    aux2_j( end+1:max(size(Cx0)) ) = NaN;
elseif max(size(Cx0))<max(size(Cy0))
    Cx0   ( end+1:max(size(Cy0)) ) = NaN;
    Phx   ( end+1:max(size(Cy0)) ) = NaN;
    aux1_i( end+1:max(size(Cy0)) ) = NaN;
    aux1_j( end+1:max(size(Cy0)) ) = NaN;
end
data = [aux1_i' aux1_j' Cx0' Phx' aux2_i' aux2_j' Cy0' Phy'];

delete('AccessValues.dat');
fid = fopen('AccessValues.dat', 'at');

fprintf(fid, '%s\n'   , '_______SLM1______________|_________SLM2___________|');
fprintf(fid, '%s\n\n' , 'M1_|M2_|_Cx0___|_Phx_____|__M1_|M2_|_Cy0___|_Phy__|');
fprintf(fid, '%.0f\t%.0f\t%.4f\t%.4f\t\t%.0f\t%.0f\t%.4f\t%.4f\n', data.' );

fclose(fid);

