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



% Evaluate the complex modulation response via evaluating the phase modulation
% and the amplitude modulation.
%
% The phase modulation is evaluated via the shift of vertical-interferences
% fringe-pattern produced when different gray level is displayed on two
% certain ROIs according to
%
%       Xi                  Xf
%   Y1   --------------------
%       |      ROI_up        | vertical interferences are inside of ROI_up
%   Y2   --------------------
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ fringes-pattern breaks here 
%   Y3   --------------------
%       |      ROI_dw        | vertical interferences are inside of ROI_dw
%   Y4   --------------------
%
% 
% To evaluate the amplitude modulation two options are considered. One is
% directly from the sum of the pixels values of each images. The other, is 
% comparing the sum of the pixels value in each ROI -descrived above- for
% every image (i.e. gray level).
% 


clear variables ; close all

%% Parameters

% ROIs definition
Xi = 250 ; Xf = 700 ;
Y1 = 320 ; Y2 = 370 ;
Y3 = 430 ; Y4 = 480 ;

% SLMs to be determined
SLM = 1; % can be [1 2]

% Positions of Pol: 1 = P@45º  ;  2 = P@135º
POLs = 1:2; % can be scalar, array or 0 to avoid evaluation

% to cheack the FFT peak 
check_peak = 0; % =1 to check the FFT peak
Pf = 24;
Q  = 1;

% Transmitions to evaluate
Ts = 1; % 0 to avoid evaluation 

% mode to determine the transmittion modulation response
whole = 0; % 0: ROIs is used ; 1: whole sum is used

% label used to read and write info
label = '310117';
path  = [pwd '/' label '_I/'];


%% Amplittude modulation

if Ts ~= 0

I0   = [0 0]; % whole energy for im0.png
fact = [0 0]; % correction factor when ROIs mode is used
I1(1:256) = 0; % array with the intensity modulation
I2(1:256) = 0; % array with the intensity modulation

for j=SLM
    
    im      = double(imread([path 'I' num2str(j) '_0.png']));
    up      = im(Y1:Y2,Xi:Xf);
    down    = im(Y3:Y4,Xi:Xf);
    fact(j) = sum(sum(up))/sum(sum(down));
    I0(j)   = sum(sum(im));
    
    for i=1:256

        filename = [path 'I' num2str(j) '_' num2str(i-1) '.png'];
        im       = double(imread(filename));

        if whole == 0
            up   = fact(j)*im(Y1:Y2,Xi:Xf);
            down = im(Y3:Y4,Xi:Xf);
        else
            up   = I0(j);
            down = sum(sum(im));
        end

        if j==1
            
            I1(i) = sum(sum(down))/sum(sum(up));
            
            if i==256
            % To flip ref to signal
            % I1 = 1./I1;
            
            % from intensity to amplitude
            A1 = sqrt(I1/max(I1));
            
            % plot and save the info
            h = figure;
            plot(A1)
            axis([1 256 0 1.1])
            cd plots
            print(h,'-depsc',[label 'amplitude_SLM1.eps']);
            cd ..
            
            end
            
        else
            
            I2(i) = sum(sum(down))/sum(sum(up));
            
            if i==256
            
            % To flip ref to signal
            % I2 = 1./I2; 
            
            % from intensity to amplitude
            A2 = sqrt(I2/max(I2));
            
            % plot and save the info
            h=figure;
            plot(A2)
            axis([1 256 0 1.1])
            cd plots
            print(h,'-depsc',[label 'amplitude_SLM2.eps']);
            cd ..
            
            end
            
        end

    end % loop over gray levels
end % loop over SLMs


end % MAIN sitch




%% Phase modulation

if POLs ~= 0

N   = zeros(1,256);
M   = zeros(1,256);
P   = zeros(1,256);
phi = zeros(1,256);
Dx  = Xf-Xi+1;

delta_phi = zeros(1,256);
FFTup  = zeros(Dx,256);
FFTdw = zeros(Dx,256);


phi1_45  = zeros(1,256);
phi1_135 = zeros(1,256);
phi2_45  = zeros(1,256);
phi2_135 = zeros(1,256);


for k=POLs % 1 => P@45º and 2 => P@135º
    if k==1 , k_str = '45'; else k_str = '135'; end
    path=[pwd '/' label '_' k_str '/'];
    
for j=SLM
    
    if check_peak == 1 % plot FFT(im) to evaluate where is the peak
        
        im   = double(imread([path 'ph' num2str(j) '_30.png']));
        down = mean(im(Y3:Y4,Xi:Xf),1);
        
        figure; plot(abs(fft(down)));
        
        Pf = input('Type the frequency of the first peak:\n  ');
        Q  = 1;%input('and the tolerance:\n  ');
        check_peak = 0;
        
    end
    
    for i=1:256
        filename = [path 'ph' num2str(j) '_' num2str(i-1) '.png'];
        im = double(imread(filename));

        up   = mean( im(Y1:Y2,Xi:Xf) ,1);
        down = mean( im(Y3:Y4,Xi:Xf) ,1);
        
        UP = fft(up)  ; FFTup(:,i) = UP;
        DW = fft(down); FFTdw(:,i) = DW;
   
        if Pf~=0
            UP(2:Pf-Q) = 0 ; UP(Pf+Q:Dx-Pf+2-Q) = 0 ; UP(Dx-Pf+2+Q:Dx) = 0;
            DW(2:Pf-Q) = 0 ; DW(Pf+Q:Dx-Pf+2-Q) = 0 ; DW(Dx-Pf+2+Q:Dx) = 0;
        end
        
        up2 = abs(ifft(UP));
        dw2 = abs(ifft(DW));
     
        if i==31 && k==min(POLs)
            
            up_draw = up2-min(up2);
            up_draw = up_draw/max(up_draw)*2-1;
            up_draw = up_draw*(Y1/2-Y2/2)+Y1/2+Y2/2;
            
            dw_draw = dw2-min(dw2);
            dw_draw = dw_draw/max(dw_draw)*2-1;
            dw_draw = dw_draw*(Y3/2-Y4/2)+Y3/2+Y4/2;
            
            figure;
            imagesc(im)
            hold on
            plot(Xi:Xf,up_draw,'LineWidth',3)
            plot(Xi:Xf,dw_draw,'LineWidth',3)
            title(['ph' num2str(j) '\_' num2str(i-1) '\_45.png'])
            hold off
        end
        
        [peakU,peakUP] = findpeaks(up2);
        [peakD,peakDW] = findpeaks(dw2);
        
        P(i) = peakUP(5)-peakUP(4); % Period
        
        N(i) = size(peakUP,2);
        M(i) = size(peakDW,2);
        if N(i)>M(i)
            D = peakUP(1:M(i))-peakDW(1:M(i));
        elseif M(i)>N(i)
            D = peakUP(1:N(i))-peakDW(1:N(i));
        else
            D = peakUP-peakDW;
        end
        
        if j==1   % --- The sign change decreasing to increasing phi(i)
            Dmean = mean(D);
        else
            Dmean =-mean(D);
        end
        phi(i) = 360*Dmean/P(i);
        
        delta_phi(i) = std(D)*360/P(i); % deviation of the mean
        
    end % loop over the 256 gray values
        
    [phi_N,pFirst] = antimod(phi,360);
    phi_N = phi_N - min(phi_N(1:pFirst));
    
    if k==1 % P@45º
        if j==1 % SLM1
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi1_45=phi_N;
            delta_phi1_45=delta_phi;
        else % SLM2
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi2_45=phi_N;
            delta_phi2_45=delta_phi;
        end
    else % P@135º
        if j==1 % SLM1
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi1_135=phi_N;
            delta_phi1_135=delta_phi;
        else % SLM2
%           phi_0=min(phi(1:35));
%           phi=phi-phi_0;
%           phi=mod(phi,360);
            phi2_135=phi_N;
            delta_phi2_135=delta_phi;
        end
    end
end % loop over SLMs
end % loop over P@Xº

if length(POLs)>1
    phi1 = mod(  ( phi1_45 + phi1_135 )/2  ,360);
    phi2 = mod(  ( phi2_45 + phi2_135 )/2  ,360);
elseif any( POLs==1 )
    phi1 = mod(phi1_45,360);
    phi2 = mod(phi2_45,360);
elseif any( POLs==2 )
    phi1 = mod(phi1_135,360);
    phi2 = mod(phi2_135,360);
else
    phi1 = zeros(1,256);
    phi2 = zeros(1,256);
end

% plot and save info
if any( SLM==1 )
h = figure; % phase for SLM1
plot( 1:256,mod(phi1_45,360),':b' , 1:256,mod(phi1_135,360),':r' , 1:256,phi1,'-k' )
cd plots
print(h,'-depsc',[label 'phase_SLM1.eps']);
cd ..
end
if any( SLM==2 )
h = figure; % phase for SLM2
plot( 1:256,mod(phi2_45,360),':b' , 1:256,mod(phi2_135,360),':r' , 1:256,phi2,'-k' )
cd plots
print(h,'-depsc',[label 'phase_SLM2.eps']);
cd ..
end


end


%% Plot and save polar info


if all([Ts SLM]) 

if any( SLM==1 )
h = figure; % polar for SLM1
polar( phi1*pi/180 , A1/max(A1) )
cd plots
print(h,'-depsc',[label 'polar_SLM1.eps']);
cd ..
end
if any( SLM==2 )
h = figure; % polar for SLM2
polar( phi2*pi/180 , A2/max(A2) )
cd plots
print(h,'-depsc',[label 'polar_SLM2.eps']);
cd ..
end


end






