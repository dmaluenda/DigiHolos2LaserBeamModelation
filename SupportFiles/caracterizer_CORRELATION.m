%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
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
Xi    = 175; % left
Xf    = 600; % rigth
ROI_C = 370; % vertical center
ROI_S = 100; % heigth
ROI_D = 100; % separation
Y1    = ROI_C-ROI_D/2-ROI_S;
Y2    = ROI_C-ROI_D/2;
Y3    = ROI_C+ROI_D/2;
Y4    = ROI_C+ROI_D/2+ROI_S ;

% SLMs to be determined
SLM = 1:2; % can be [1 2]

% Positions of Pol: 1 = P@45º  ;  2 = P@135º
POLs = 1:2; % can be scalar, array or 0 to avoid phase evaluation

% to cheack the FFT peak 
check_peak = 0; % =1 to check the FFT peak
Pf = 5;
Q  = 1;

% Transmitions to evaluate
Ts = 1; % 0 to avoid amplitude evaluation 

% mode to determine the transmittion modulation response
whole   = 3; % 0: ROIs is used ; 1: Whole sum is used ; 3: PowerMeter eval.
Ioffset = 1.7; % some Intensity offset

% label used to read and write info
label = '150317';



%% Amplittude modulation

if Ts ~= 0

path  = [pwd '/' label '_I'];
I0   = zeros(1,length(SLM)); % whole energy for im0.png
fact = zeros(1,length(SLM)); % correction factor when ROIs mode is used
Is   = zeros(256,length(SLM)); % array with the intensity modulation

j = 0;
for slm=SLM
j=j+1; % slm index

    if whole == 3 % PowerMeter evaluation
        
        filename = [path '/T' num2str(slm) '.txt'];
        Idata    = dlmread( filename );
        
        units   = round(Idata(:,3)); % 1: mW ; 3: uW ; 6: nW ; 9:pW
        Is(:,j) = Idata(:,2).*10.^(-units);
        
    else   % CCD evaluation

        im      = double(imread([path '/I' num2str(slm) '_0.png']));
        up      = im(Y1:Y2,Xi:Xf);
        dw      = im(Y3:Y4,Xi:Xf);
        fact(j) = sum(sum(up))/sum(sum(dw));
        I0(j)   = sum(sum(im));
        
        for i=2:256 % loop over all images
            
            filename = [path '/I' num2str(slm) '_' num2str(i-1) '.png'];
            im       = double(imread(filename));

            if whole == 0 % using ROI 
                up = fact(j)*im(Y1:Y2,Xi:Xf);
                dw = im(Y3:Y4,Xi:Xf);
            else   % whole CCD
                up = I0(j);
                dw = sum(sum(im));
            end

            Is(i,j) = sum(sum(dw))/sum(sum(up));

        end % loop over the 256 gl
        
    end % eval. method switch
    
    if slm==1 % SLM1
        
        I1 = Is(:,j).';


        % To flip ref to signal on ROI
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
        
    elseif slm==2 %SLM2
        
        I2 = Is(:,j).';
            
        % To flip ref to signal on ROI
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
        
    end % SLM switch
    
end % LOOP over SLM


end % MAIN switch




%% Phase modulation

if POLs ~= 0

N   = zeros(1,256);
M   = zeros(1,256);
P   = zeros(1,256);
D   = zeros(1,256);
phi = zeros(1,256);
Dx  = Xf-Xi+1;

delta_phi = zeros(1,256);

phi1_45    = zeros(1,256);
phi1_135   = zeros(1,256);
phi2_45    = zeros(1,256);
phi2_135   = zeros(1,256);
phi1_45CC  = zeros(1,256);
phi1_135CC = zeros(1,256);
phi2_45CC  = zeros(1,256);
phi2_135CC = zeros(1,256);


for k=POLs % 1 => P@45º and 2 => P@135º
    if k==1 , k_str = '45'; else k_str = '135'; end
    path=[pwd '/' label '_' k_str '/'];
    
for j=SLM  % 1 => SLM1  and 2 => SLM2
    
    if check_peak == 1 % plot FFT(im) to evaluate where is the peak
        
        im   = double(imread([path 'ph' num2str(j) '_30.png']));
        dw = mean(im(Y3:Y4,Xi:Xf),1);
        
        figure; plot(abs(fft(dw)));
        
        Pf = input('Type the frequency of the first peak:\n  ');
        %check_peak = 0;
        
    end
    
    for i=1:256
        filename = [path 'ph' num2str(j) '_' num2str(i-1) '.png'];
        im = double(imread(filename));

        up   = mean( im(Y1:Y2,Xi:Xf) ,1);
        dw = mean( im(Y3:Y4,Xi:Xf) ,1);
        
        if Q~=0 % Band pass filtter centred in Pf and Q width 
        
            UP = fft(up)  ;
            DW = fft(dw);
        
            UP(2:Pf-Q) = 0 ; UP(Pf+Q:Dx-Pf+2-Q) = 0 ; UP(Dx-Pf+2+Q:Dx) = 0;
            DW(2:Pf-Q) = 0 ; DW(Pf+Q:Dx-Pf+2-Q) = 0 ; DW(Dx-Pf+2+Q:Dx) = 0;

            up = abs(ifft(UP));
            dw = abs(ifft(DW));
            
        end
        
        if i==61  % Draw one case to check that all is correct
            
            up_draw = up-min(up);
            up_draw = up_draw/max(up_draw)*2-1;
            up_draw = up_draw*(Y1/2-Y2/2)+Y1/2+Y2/2;
            
            dw_draw = dw-min(dw);
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
        
        CrossCorr = xcorr(up-min(up),dw-min(dw)); % ,'coeff'
        [~,pCC]   = max(CrossCorr);
        Dcc(i)    = -abs(pCC - length(CrossCorr));
        
        [~,peakCC] = findpeaks(CrossCorr);
        
        % Period
        Pcc(i) = abs(mean(peakCC(2:end)-peakCC(1:end-1))) ;
        
        % find maximums to evaluate the displacement and the period
        [peakU,peakUP] = findpeaks(up);
        [peakD,peakDW] = findpeaks(dw);
        
        % Period
        P(i) = (  abs(mean(peakUP(2:end)-peakUP(1:end-1)))  +   ...
                  abs(mean(peakDW(2:end)-peakDW(1:end-1)))  ) / 2 ;
        
        % Displacement
        N(i) = length(peakUP);
        M(i) = length(peakDW);
        if N(i)>M(i)
            Di(1:M(i),i) = peakUP(1:M(i))-peakDW(1:M(i));
        elseif M(i)>N(i)
            Di(1:N(i),i) = peakUP(1:N(i))-peakDW(1:N(i));
        else
            Di(1:N(i),i) = peakUP-peakDW;
        end
        
        if j==1 % --- The sign changes from a decreasing to a increasing phi
            D(i) =mean(Di(:,i));
        else
            D(i) =mean(Di(:,i));
        end

        delta_D(i) = std(Di(:,i)); % deviation of the mean
        
    end % loop over the 256 gray values
    
    % from pixels to rad
    phi       =    D    * 360 ./ mean(P);
    phiCC     =   Dcc   * 360 ./ mean(P);
    delta_phi = delta_D * 360 ./ mean(P);
    
    % to avoid jumps of 360º
    [phi_N,pFirst] = antimod(phi,360);
    if pFirst == 0,pFirst=20;end
    phi_N = phi_N - min(phi_N(1:pFirst));
    
    % to avoid jumps of 360º
    [phiCC,pFirstCC] = antimod(phiCC,360);
    if pFirstCC == 0,pFirstCC=20;end
    phiCC = phiCC - min(phiCC(1:pFirst));
    
    if k==1 % P@45º
        if j==1 % SLM1
            phi1_45 = phi_N;
            phi1_45CC = phiCC;
            delta_phi1_45 = delta_phi;
        else % SLM2
            phi2_45 = phi_N;
            phi2_45CC = phiCC;
            delta_phi2_45 = delta_phi; 
        end
    else % P@135º
        if j==1 % SLM1
            phi1_135 = phi_N;
            phi1_135CC = phiCC;
            delta_phi1_135 = delta_phi;

        else % SLM2
            phi2_135 = phi_N;
            phi2_135CC = phiCC;
            delta_phi2_135 = delta_phi;
        end
    end
end % loop over SLMs
end % loop over P@Xº

 phi1_45(150:end)  = phi1_135(150:end);
 phi2_45(1:100)    = phi2_45CC(1:100);
 phi2_135(1:100)   = phi2_135CC(1:100);
 phi2_45CC(129:end)= phi2_45(129:end);
% phi1_135(1:60)   = phi1_45(1:60);
% phi1_135CC       = phi1_45CC;

if length(POLs)>1
    phi1   = mod(  ( phi1_45 + phi1_135 )/2  ,360);
    phi2   = mod(  ( phi2_45 + phi2_135 )/2  ,360);
    phi1CC = mod(  (phi1_45CC+phi1_135CC)/2  ,360);
    phi2CC = mod(  (phi2_45CC+phi2_135CC)/2  ,360);
elseif any( POLs==1 )
    phi1   = mod( (phi1_45+phi1_45CC)/2 ,360);
    phi2   = mod( (phi2_45+phi2_45CC)/2 ,360);
    phi1CC = mod(phi1_45CC,360);
    phi2CC = mod(phi2_45CC,360);
elseif any( POLs==2 )
    phi1   = mod( phi1_135 ,360);
    phi2   = mod( phi2_135 ,360);
    phi1CC = mod(phi1_135CC,360);
    phi2CC = mod(phi2_135CC,360);
else
    return;
end

% plot and save info
if any( SLM==1 )
h = figure; % phase for SLM1
plot( 1:256,mod(phi1_45,360),':b' , 1:256,mod(phi1_135,360),':r' , ...
      1:256,phi1,'-k' , 1:256,phi1_45CC,'xk' , 1:256,phi1_135CC,'+k')
hold on;plot(1:256,(phi1+phi1CC)/2,'LineWidth',3)
cd plots
print(h,'-depsc',[label 'phase_SLM1.eps']);
cd ..

h=figure;
plotPMsigma((1:256)',A1,ones(256,1)'*0.02)
end
if any( SLM==2 )
h = figure; % phase for SLM2
plot( 1:256,mod(phi2_45,360),':b' , 1:256,mod(phi2_135,360),':r' , ...
      1:256,phi2,'-k' , 1:256,phi2_45CC,'xk' , 1:256,phi2_135CC,'+k')
hold on;plot(1:256,(phi2+phi2CC)/2,'LineWidth',3)
cd plots
print(h,'-depsc',[label 'phase_SLM2.eps']);
cd ..
phi2 = (phi2+phi2CC)/2;

h=figure;
plotPMsigma((1:256)',phi2,abs(phi2_45CC-phi2_135CC));
end


end


%% Plot and save polar info


if all([Ts POLs]) 

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




%% PLOTING ALL INFO  [FIGURE 3 in the paper]

% parameters to set up figure
plotheight = 20;
plotwidth  = 16;
subplotsx  = 2;
subplotsy  = 1;   
leftedge   = 1.5;
rightedge  = 1.5;   
topedge    = 2;
bottomedge = 3;
spacex     = 0.5;
spacey     = 0.2;
fontsize   = 15;    
 
sub_pos = subplot_pos( plotwidth,plotheight , ...
                       leftedge,rightedge,bottomedge,topedge , ...
                       subplotsx,subplotsy , spacex,spacey );
 
%setting up the figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'Position', [500 300 900 600]);
set(gcf, 'Name','FIGURE 3 on the PAPER')
set(gcf, 'Color',[1 1 1]);



% loop creating axes and displaying information
for i  = 1:subplotsx; % axes from left to right
for ii = 1:subplotsy; % axes from down to up

% crate the current axes
ax = axes( 'position',sub_pos{i,ii} , ...
            'XMinorGrid','off' , ...
           'FontSize',fontsize , 'Box','on' , 'Layer','top' );


if i==1
    [h1,h3] = plotPMsigma((1:256)',A1,ones(256,1)'*0.02);
    axis([1 256 0 1.1]);
    ax.YTick = 0:0.25:1;
    ax.YLabel.String = 'Amplitude';
else
    phErr = abs(phi2_45CC-phi2_135CC)/sqrt(2);
    phErr(phErr<10)=10;
    [h1,h3] = plotPMsigma((1:256)',phi2,phErr);
    axis([1 256 0 360]);
    ax.YTick = 0:90:360;
    ax.YLabel.String = 'Phase';
end
hold on





ax.FontSize      = 20;
ax.XTick         = 0:64:256;
% ax.XTickLabels   = {'2','1','0','1','2'}; % r>0 allways
ax.XLabel.String = 'Gray level applied';
ax.LineWidth     = 2;

if i==1 % left column
    legend off
end
if i==2 % central column
    ax.YAxisLocation='right';
    hl=legend([h3 h1],...
        'SLM response','Estimateted Error');
    hl.Position=[0.4364 0.52 0.1359 0.1620];
end

ax.GridLineStyle = ':'
grid on

hold off


end
end



