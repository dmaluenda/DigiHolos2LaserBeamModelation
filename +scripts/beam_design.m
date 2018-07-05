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


% David Maluenda Niubó - FAO (UB) [2014]
%
% [Ex,Ey,Px,Py] = beam_design(SLMsize,BeamType,infile,draw,stokes,gauss_correction)
%
% Returns tha amplitude of X and Y components and the phase between both. 
% Beam_type is the label of the designed beam
%
%   Example:  SLMsize = [1024 768]
%             BeamType = last
%             infile = 0
%             draw = 0
%             stokes = 0
%             gauss_correction=0
% --------------------------------------           

function[E_x,E_y,Ph_x,Ph_y] = beam_design(SLM_size, beam_type, infile, draw, ...
                                          stokes, gauss_correction)

    %% --- TO USE IT AS SCRIPT -------
    if ~exist('SLM_size','var'), clear all, SLM_size = [1024 768]; end
    if ~exist('beam_type','var'), beam_type = 45; end
    if ~exist('infile','var'), infile = 0; end
    if ~exist('draw','var'), draw = 1; end 
    if ~exist('stokes','var'), stokes = 0; end
    if ~exist('gauss_correction','var'), gauss_correction = 0; end
    % --------------------------------


    PATH = cd;
    designPATH = [PATH '\Designs\design'];


    [y,x] = meshgrid( SLM_size(2)/2:-1:-SLM_size(2)/2+1, ...
                     -SLM_size(1)/2+ 1: SLM_size(1)/2);
    theta = mod(atan2(y,x),2*pi); %[0 2pi]
    rho   = sqrt(x.^2+y.^2);
    E_x   = zeros(SLM_size) ;
    E_y = E_x ; Ph_x = E_x ; Ph_y = E_x ;


    % --- MAIN SWICH ------
    % Add new beam types here
    switch beam_type

        case 86
        
            N    = 8;
            m    = 8;
            As   = NaN;
            NA   = 0.6;
            rhoM = 165/2;

            beamNAME = ['Vec.Needle N=' num2str(N) ' m=' num2str(m) ...
                        ' As=' num2str(As)];

            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*cos(theta);
            Ey = h.*sin(theta);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        

        case 85
        
            N    = 4;
            m    = 4;
            As   = NaN;
            NA   = 0.6;
            rhoM = 165/2;

            beamNAME = ['Vec.Needle N=' num2str(N) ' m=' num2str(m) ...
                        ' As=' num2str(As)];

            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*cos(theta);
            Ey = h.*sin(theta);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        

        case 84
        
            N    = 8;
            m    = 8;
            As   = NaN;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Vec.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*cos(theta);
            Ey = h.*sin(theta);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
 

         case 83
        
            N    = 8;
            m    = 8;
            As   = -100;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = -h;
            Ey = -h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
            
            
        case 82
        
            N    = 8;
            m    = 8;
            As   = 14;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);


        case 81
      
            N    = 8;
            m    = 8;
            As   = 6;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 80
        
            N    = 8;
            m    = 8;
            As   = 12;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);


         case 79
      
            N    = 8;
            m    = 8;
            As   = 8;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        

        case 78          
      
            N    = 8;
            m    = 8;
            As   = -20;
            NA   = 0.6;
            rhoM = 165/2;
            
    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];
        
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        

        case 77
      
            N    = 8;
            m    = 8;
            As   = -10;
            NA   = 0.6;
            rhoM = 165/2;
            
    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);


        case 76

            N    = 8;
            m    = 8;
            As   = 20;
            NA   = 0.6;
            rhoM = 165/2;
            
    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        

        case 75  
      
            N    = 8;
            m    = 8;
            As   = 10;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        

         case 74    
      
            N    = 8;
            m    = 8;
            As   = -2;
            NA   = 0.6;
            rhoM = 165/2;
            
    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 73   
      
            N    = 8;
            m    = 8;
            As   = -1.5;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 72
            
            N    = 8;
            m    = 8;
            As   = -1;
            NA   = 0.6;
            rhoM = 165/2;
            
    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 71   
      
            N    = 8;
            m    = 8;
            As   = -0.5;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        

        case 70
           
            N    = 8;
            m    = 8;
            As   = 2;
            NA   = 0.6;
            rhoM = 165/2;
           
    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 69 
      
            N    = 8;
            m    = 8;
            As   = 1.5;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 68     
      
            N    = 8;
            m    = 8;
            As   = 1;
            NA   = 0.6;
            rhoM = 165/2;
      
    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
        
        case 67   
      
            N    = 8;
            m    = 8;
            As   = 0.5;
            NA   = 0.6;
            rhoM = 165/2;

    beamNAME=['Esc.Needle N=' num2str(N) ' m=' num2str(m) ' As=' num2str(As)];

            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h.*exp(1i*2*pi*As*sin(teta).^4);
            Ey = h.*exp(1i*2*pi*As*sin(teta).^4);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 66
            beamNAME='Esc.Needle N=5 m=5';
        
      
            N    = 5;
            m    = 5;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 65
            beamNAME='Esc.Needle N=5 m=5';
        
      
            N    = 5;
            m    = 5;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 64
            beamNAME='Esc.Needle N=8 m=7';
        
      
            N    = 8;
            m    = 7;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
        
        case 63
            beamNAME='Esc.Needle N=10 m=5';
        
      
            N    = 10;
            m    = 5;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
        
        case 62
            beamNAME='Esc.Needle N=10 m=8';
        
      
            N    = 10;
            m    = 8;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        case 61
            beamNAME='Esc.Needle N=6 m=6';
        
      
            N    = 6;
            m    = 6;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
        
        case 60
            beamNAME='Esc.Needle N=3 m=3';
        
      
            N    = 3;
            m    = 3;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
        
        case 59
            beamNAME='Esc.Needle N=2 m=2';
        
      
            N    = 2;
            m    = 2;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
          case 58
            beamNAME='Esc.Needle N=5 m=5';
        
      
            N    = 5;
            m    = 5;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
          case 57
            beamNAME='Esc.Needle N=4 m=4';
        
      
            N    = 4;
            m    = 4;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
         case 56
            beamNAME='Esc.Needle N=8 m=8';
        
      
            N    = 8;
            m    = 8;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);


         case 55
            beamNAME='Esc.Needle N=13 m=7';
        
      
            N    = 13;
            m    = 7;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        
        case 54
            beamNAME='Esc.Needle N=13 m=13';
        
      
            N    = 13;
            m    = 13;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
        
        case 53
            beamNAME='Esc.Needle N=15 m=2';
        
      
            N    = 15;
            m    = 2;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
        
        case 52
            beamNAME='Esc.Needle N=8 m=2';
        
      
            N    = 8;
            m    = 2;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
              
        case 51
            beamNAME='Esc.Needle N=8 m=4';
        
      
            N    = 8;
            m    = 4;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
         case 50
            beamNAME='Esc.Needle N=10 m=2';
        
      
            N    = 10;
            m    = 2;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
  
        
            case 49
            beamNAME='Esc.Needle N=10 m=4';
        
      
            N    = 10;
            m    = 4;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
              
        
        case 48
            beamNAME='Esc.Needle N=10 m=10';
        
      
            N    = 10;
            m    = 10;
            NA   = 0.6;
            rhoM = 150/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);

        
        
         case 47
            beamNAME='Esc.Needle N=10 m=10';
        
      
            N    = 10;
            m    = 10;
            NA   = 0.6;
            rhoM = 165/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;
            Ey = h;
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
              
        
        case 46
            beamNAME='Esc.Needle N=10 m=10';
        
      
            N    = 10;
            m    = 10;
            NA   = 0.6;
            rhoM = 180/2;
            
            mask = (rho<=rhoM).*1;
            teta = asin(rho/rhoM*NA).*mask;
            
            alpha = cos(teta);
            alpha0= cos(max(teta(:)));
            alpha1= cos(min(teta(:)));
            alphaB= m/2/N*(1-alpha0) + alpha0 ;
         
            h = N*sinc( 2*N * (alpha-alphaB)./(alpha1-alpha0) ) .*mask ;
            
            h = h/max(h(:));
            
            Ex = h;%.*cos(theta);
            Ey = h;%.*sin(theta);
            
            E_x  = abs(Ex);
            Ph_x = angle(Ex);
            E_y  = abs(Ey);
            Ph_y = angle(Ey);
        
              
        
        case 45
            beamNAME='Encryption3 OnlyWhiteNoise64';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=scripts.normalize_2D(fftshift(fft2(fftshift(im))));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=64; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(Cmask1))+0.1;
            Ey=scripts.normalize_2D(abs(Cmask2))+0.1;
            Ex=Ex/max(Ex(:));Ey=Ey/max(Ey(:));
            Phx=angle(Cmask1);
            Phy=angle(Cmask2);
            Cmask1=Ex.*exp(1i*Phx);
            Cmask2=Ey.*exp(1i*Phy);
            
            N=10;
            Cmask1p=Cmask1;
            for i=1:N
            aux=randi(512*384);
            Cmask1p(aux)=Cmask1(aux)+double(randi(200)/100)-1;
            end
            Cmask1p(Cmask1p>1)=Cmask1p(Cmask1p>1)./abs(Cmask1p(Cmask1p>1));
            
            %q=2;
            %Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            %Ex=scripts.normalize_2D(Ex);
            %TH=0.04;
            %Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            focalP=IM.*Cmask1;
            CCD=scripts.normalize_2D(ifftshift(ifft2(ifftshift(focalP))));
            N=1000;
            CCDp=CCD;
            i=1;
            for i=1:N
            aux=round(randi(512*384));
            AUX=CCD(aux);
            CDDp(aux)=(AUX+double(randi(200)/100)-1)*exp(1i*2*pi*double(randi(100)/100));
            end
            out=ifftshift(ifft2(ifftshift( fftshift(fft2(fftshift(CCDp)))./Cmask1 )));
            figure;imagesc(abs(focalP))
            figure;imagesc(abs(CCD));
            figure;imagesc(abs(out)) 
        
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
            
            dlmwrite('Designs/design43_Ex.txt',E_x);
            dlmwrite('Designs/design43_Phx.txt',Ph_x);
        
        
        case 44
            beamNAME='Encryption3 OnlyWhiteNoise32';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=scripts.normalize_2D(fftshift(fft2(fftshift(im))));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=32; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(Cmask1))+0.1;
            Ey=scripts.normalize_2D(abs(Cmask2))+0.1;
            Ex=Ex/max(Ex(:));Ey=Ey/max(Ey(:));
            Phx=angle(Cmask1);
            Phy=angle(Cmask2);
            Cmask1=Ex.*exp(1i*Phx);
            Cmask2=Ey.*exp(1i*Phy);
            
            N=10;
            Cmask1p=Cmask1;
            for i=1:N
            aux=randi(512*384);
            Cmask1p(aux)=Cmask1(aux)+double(randi(200)/100)-1;
            end
            Cmask1p(Cmask1p>1)=Cmask1p(Cmask1p>1)./abs(Cmask1p(Cmask1p>1));
            
            %q=2;
            %Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            %Ex=scripts.normalize_2D(Ex);
            %TH=0.04;
            %Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            focalP=IM.*Cmask1;
            CCD=scripts.normalize_2D(ifftshift(ifft2(ifftshift(focalP))));
            N=1000;
            CCDp=CCD;
            i=1;
            for i=1:N
            aux=round(randi(512*384));
            AUX=CCD(aux);
            CDDp(aux)=(AUX+double(randi(200)/100)-1)*exp(1i*2*pi*double(randi(100)/100));
            end
            out=ifftshift(ifft2(ifftshift( fftshift(fft2(fftshift(CCDp)))./Cmask1 )));
            figure;imagesc(abs(focalP))
            figure;imagesc(abs(CCD));
            figure;imagesc(abs(out)) 
        
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
            
            dlmwrite('Designs/design43_Ex.txt',E_x);
            dlmwrite('Designs/design43_Phx.txt',Ph_x);
        
        
       
        case 43
            beamNAME='Encryption3 OnlyWhiteNoise4';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=scripts.normalize_2D(fftshift(fft2(fftshift(im))));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=4; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(Cmask1))+0.1;
            Ey=scripts.normalize_2D(abs(Cmask2))+0.1;
            Ex=Ex/max(Ex(:));Ey=Ey/max(Ey(:));
            Phx=angle(Cmask1);
            Phy=angle(Cmask2);
            Cmask1=Ex.*exp(1i*Phx);
            Cmask2=Ey.*exp(1i*Phy);
            
            N=10;
            Cmask1p=Cmask1;
            for i=1:N
            aux=randi(512*384);
            Cmask1p(aux)=Cmask1(aux)+double(randi(200)/100)-1;
            end
            Cmask1p(Cmask1p>1)=Cmask1p(Cmask1p>1)./abs(Cmask1p(Cmask1p>1));
            
            %q=2;
            %Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            %Ex=scripts.normalize_2D(Ex);
            %TH=0.04;
            %Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            focalP=IM.*Cmask1;
            CCD=scripts.normalize_2D(ifftshift(ifft2(ifftshift(focalP))));
            N=1000;
            CCDp=CCD;
            i=1;
            for i=1:N
            aux=round(randi(512*384));
            AUX=CCD(aux);
            CDDp(aux)=(AUX+double(randi(200)/100)-1)*exp(1i*2*pi*double(randi(100)/100));
            end
            out=ifftshift(ifft2(ifftshift( fftshift(fft2(fftshift(CCDp)))./Cmask1 )));
    %         figure;imagesc(abs(focalP))
    %         figure;imagesc(abs(CCD));
    %         figure;imagesc(abs(out)) 
        
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
            
            dlmwrite('Designs/design43_Ex.txt',E_x);
            dlmwrite('Designs/design43_Phx.txt',Ph_x);
        
        case 42
            beamNAME='Encryption3 OnlyWhiteNoise8';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=scripts.normalize_2D(fftshift(fft2(fftshift(im))));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=8; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(Cmask1))+0.1;
            Ey=scripts.normalize_2D(abs(Cmask2))+0.1;
            Ex=Ex/max(Ex(:));Ey=Ey/max(Ey(:));
            Phx=angle(Cmask1);
            Phy=angle(Cmask2);
            Cmask1=Ex.*exp(1i*Phx);
            Cmask2=Ey.*exp(1i*Phy);
            
            
            %q=2;
            %Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            %Ex=scripts.normalize_2D(Ex);
            %TH=0.04;
            %Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            focalP=IM.*Cmask1;
            CCD=ifftshift(ifft2(ifftshift(focalP)));
            out=ifftshift(ifft2(ifftshift( fftshift(fft2(fftshift(CCD)))./Cmask1 )));
    %         figure;imagesc(abs(focalP))
    %         figure;imagesc(abs(CCD));
    %         figure;imagesc(abs(out)) 
        
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;

            dlmwrite('Designs/design43_Ex.txt',E_x);
            dlmwrite('Designs/design43_Phx.txt',Ph_x);
        
        case 41
            beamNAME='Encryption3 OnlyWhiteNoise16';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=fftshift(fft2(fftshift(im)));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=16; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(Cmask1))+0.1;
            Ey=scripts.normalize_2D(abs(Cmask2))+0.1;
            Ex=Ex/max(Ex(:));Ey=Ey/max(Ey(:));
            Phx=angle(Cmask1);
            Phy=angle(Cmask2);
            Cmask1=Ex.*exp(1i*Phx);
            Cmask2=Ey.*exp(1i*Phy);
            
            
            %q=2;
            %Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            %Ex=scripts.normalize_2D(Ex);
            %TH=0.04;
            %Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            focalP=IM.*Cmask1;
            CCD=ifftshift(ifft2(ifftshift(focalP)));
            out=ifftshift(ifft2(ifftshift( fftshift(fft2(fftshift(CCD)))./Cmask1 )));
    %         figure;imagesc(abs(focalP))
    %         figure;imagesc(abs(CCD));
    %         figure;imagesc(abs(out)) 
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
        
            dlmwrite('Designs/design43_Ex.txt',E_x);
            dlmwrite('Designs/design43_Phx.txt',Ph_x);
        
            case 40
            beamNAME='Encryption2 PerlinNoise';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=scripts.normalize_2D(ifftshift(ifft2(ifftshift(im))));
            
            mask1=zeros(SLM_size);
            mask2=zeros(SLM_size);
            mask3=zeros(SLM_size);
            mask4=zeros(SLM_size);
     
            
            %RandStream.setGlobalStream(RandStream('mt19937ar','seed',0));
            key = zeros(SLM_size(1),SLM_size(2),4);    % output image
            
            for mk=1:4
                w = max(SLM_size);  % width of current layer
                i = 0;       % iterations
                while w >3
                    i = i + 1;
                    d = interp2(randn(w), i-1, 'splines');
                    key(:,:,mk) = key(:,:,mk) + i * d(1:SLM_size(1), 1:SLM_size(2));
                    w = w - ceil(w/2 - 1);
                end
                key(:,:,mk)=scripts.normalize_2D(key(:,:,mk));
            end
            
            key1=key(:,:,1);key2=key(:,:,2);key3=key(:,:,3);key4=key(:,:,4);
            
            key1=key1+0.1;
            key1=key1/max(key1(:));
            
            Cmask1=key1.*exp(1i*2*pi*key2);
            Cmask2=key3.*exp(1i*2*pi*key4);
            
            fCmask1=fft2(Cmask1);
            fCmask1(1,1)=0;
            Cmask1=ifft2(fCmask1);
            
            
            q=3;
            IM(513-q:513+q,385-q:385+q)=0;IM(512:514,384:386)=0;
            E_x=scripts.normalize_2D(abs(IM.*Cmask1));
            E_y=abs(Cmask2);
            Ph_x=angle(IM.*Cmask1);
            Ph_y=angle(Cmask2);
            
            E_x=scripts.normalize_2D(E_x);
    %         TH=0.07;
    %         E_x(E_x<TH)=TH;Ph_x(E_x<TH)=rand*2*pi;
         
            
            CCD=fftshift(fft2(fftshift(E_x.*exp(1i*Ph_x))));
            out=ifftshift(ifft2(ifftshift(CCD)))./Cmask1;
            figure;imagesc(E_x)
            figure;imagesc(abs(CCD));
            figure;imagesc(abs(fftshift(fft2(fftshift(out))))) 
        
        
        case 39
            beamNAME='Encryption2 Uncrypted';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=scripts.normalize_2D(ifftshift(ifft2(ifftshift(im))));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=8; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk))+0.1;
            key1=key1/max(key1(:));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk))+0.1;
            key3=key3/max(key3(:));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(IM));
            Ey=abs(Cmask2);
            Phx=angle(IM.*Cmask1);
            Phy=angle(Cmask2);
            
            q=1;
            Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            Ex=scripts.normalize_2D(Ex);
    %         TH=0.01;
    %         Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            
            CCD=fftshift(fft2(fftshift(Ex.*exp(1i*Phx))));
            out=ifftshift(ifft2(ifftshift(CCD)))./Cmask1;
    %         figure;imagesc(Ex)
    %         figure;imagesc(abs(CCD));
    %         figure;imagesc(abs(fftshift(fft2(fftshift(out))))) 
        
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
            E_y=E_x;
            Ph_y=Ph_x;
            
            
          case 38
            beamNAME='Encryption2 WhiteNoise8';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=scripts.normalize_2D(ifftshift(ifft2(ifftshift(im))));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=8; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk))+0.1;
            key1=key1/max(key1(:));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk))+0.1;
            key3=key3/max(key3(:));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(IM.*Cmask1));
            Ey=abs(Cmask2);
            Phx=angle(IM.*Cmask1);
            Phy=angle(Cmask2);
            
            q=2;
            Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            Ex=scripts.normalize_2D(Ex);
    %         TH=0.01;
    %         Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            
            CCD=fftshift(fft2(fftshift(Ex.*exp(1i*Phx))));
            out=ifftshift(ifft2(ifftshift(CCD)))./Cmask1;
    %         figure;imagesc(Ex)
    %         figure;imagesc(abs(CCD));
    %         figure;imagesc(abs(fftshift(fft2(fftshift(out))))) 
        
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
            E_y=E_x;
            Ph_y=Ph_x;
        
        case 37
            beamNAME='Encryption2 WhiteNoise16';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=ifftshift(ifft2(ifftshift(im)));
            
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=16; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk))+0.1;
            key1=key1/max(key1(:));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(IM.*Cmask1));
            Ey=abs(Cmask2);
            Phx=angle(IM.*Cmask1);
            Phy=angle(Cmask2);
            
            
            q=2;
            Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            Ex=scripts.normalize_2D(Ex);
            %TH=0.04;
            %Ex(Ex<TH)=TH;Phx(Ex<TH)=rand*2*pi;
         
            
            CCD=fftshift(fft2(fftshift(Ex.*exp(1i*Phx))));
            out=ifftshift(ifft2(ifftshift(CCD)))./Cmask1;
    %         figure;imagesc(Ex)
    %         figure;imagesc(abs(CCD));
    %         figure;imagesc(abs(fftshift(fft2(fftshift(out))))) 
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
        
            E_y=E_x;
            Ph_y=Ph_x;
            
            
            
        case 36
            beamNAME='Encryption2 WhiteNoise32';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            IM=ifftshift(ifft2(ifftshift(im)));
                    
            mask1=zeros(SLM_size/2);
            mask2=zeros(SLM_size/2);
            mask3=zeros(SLM_size/2);
            mask4=zeros(SLM_size/2);
            
    %         RandStream.setStandardStream(RandStream('mt19937ar','seed',0)); %seed
            stk=32; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk))+0.1;
            key1=key1/max(key1(:));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk))+0.1;
            key3=key3/max(key3(:));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk/2-1
                for j=0:SLM_size(2)/stk/2-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            Ex=scripts.normalize_2D(abs(IM.*Cmask1));
            Ey=abs(Cmask2);
            Phx=angle(IM.*Cmask1);
            Phy=angle(Cmask2);
            
            
            q=1;
            Ex(257-q:257+q,193-q:193+q)=0;Phx(257-q:257+q,193-q:193+q)=0;
            Ex=scripts.normalize_2D(Ex);
            %TH=0.01;
            %E_x(E_x<TH)=TH;Ph_x(E_x<TH)=rand*2*pi;
         
            
            CCD=fftshift(fft2(fftshift(Ex.*exp(1i*Phx))));
            out=ifftshift(ifft2(ifftshift(CCD)))./Cmask1;
    %         figure;imagesc(Ex)
    %         figure;imagesc(abs(CCD));
    %         figure;imagesc(abs(fftshift(fft2(fftshift(out))))) 
            
            E_x(1:2:end-1,1:2:end-1)=Ex;E_x(2:2:end,1:2:end-1)=Ex;
            E_y(1:2:end-1,1:2:end-1)=Ey;E_y(2:2:end,1:2:end-1)=Ey;
            Ph_x(1:2:end-1,1:2:end-1)=Phx;Ph_x(2:2:end,1:2:end-1)=Phx;
            Ph_y(1:2:end-1,1:2:end-1)=Phy;Ph_y(2:2:end,1:2:end-1)=Phy;
            E_x(1:2:end-1,2:2:end)=Ex;E_x(2:2:end,2:2:end)=Ex;
            E_y(1:2:end-1,2:2:end)=Ey;E_y(2:2:end,2:2:end)=Ey;
            Ph_x(1:2:end-1,2:2:end)=Phx;Ph_x(2:2:end,2:2:end)=Phx;
            Ph_y(1:2:end-1,2:2:end)=Phy;Ph_y(2:2:end,2:2:end)=Phy;
            
            E_y=E_x;
            Ph_y=Ph_x;
            
        case 35
            beamNAME='Encryption PerlinNoise2';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',0));
            stk=1;
            m=SLM_size/stk;  
            key = zeros(m(1),m(2),4);    % output image

            
            for mk=1:4
                w = max(m);  % width of current layer
                i = 0;       % iterations
                while i < 4
                    i = i + 1;
                    d = interp2(randn(w), i-1, 'splines');
                    key(:,:,mk) = key(:,:,mk) + i * d(1:m(1), 1:m(2));
                    w = w - ceil(w/2 - 1);
                end
                key(:,:,mk)=scripts.normalize_2D(key(:,:,mk));
            end
            
            for i=0:SLM_size(1)/stk-1
                for j=0:SLM_size(2)/stk-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,2);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,3);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,4);
                end
            end
            
            mask2=mod(mask2*2*pi,2*pi);
            mask4=mod(mask4*2*pi,2*pi);
            
            Cmask1=mask1.*exp(1i*mask2);
            Cmask2=mask3.*exp(1i*mask4);
            
            E_y=mask1.*im;
            E_x=mask3;
            Ph_y=zeros(SLM_size);%mask2;
            Ph_x=mask4;
            
            
        
        
        case 34
            beamNAME='Encryption WhiteNoise2';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            
            mask1=zeros(SLM_size);
            mask2=zeros(SLM_size);
            mask3=zeros(SLM_size);
            mask4=zeros(SLM_size);
            
            RandStream.setDefaultStream(RandStream('mt19937ar','seed',0)); %seed
            stk=32; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk-1
                for j=0:SLM_size(2)/stk-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            E_y=mask1.*im;
            E_x=mask3;
            Ph_y=zeros(SLM_size);%mask2*2*pi;
            Ph_x=mask4*2*pi;
        
        
        case 33
            beamNAME='Encryption PerlinNoise';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            
            RandStream.setDefaultStream(RandStream('mt19937ar','seed',0));
            stk=1;
            m=SLM_size/stk;  
            key = zeros(m(1),m(2),4);    % output image

            
            for mk=1:4
                w = max(m);  % width of current layer
                i = 0;       % iterations
                while i < 4
                    i = i + 1;
                    d = interp2(randn(w), i-1, 'splines');
                    key(:,:,mk) = key(:,:,mk) + i * d(1:m(1), 1:m(2));
                    w = w - ceil(w/2 - 1);
                end
                key(:,:,mk)=scripts.normalize_2D(key(:,:,mk));
            end
            
            for i=0:SLM_size(1)/stk-1
                for j=0:SLM_size(2)/stk-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,2);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,3);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key(i+1,j+1,4);
                end
            end
            
            mask2=mod(mask2*2*pi,2*pi);
            mask4=mod(mask4*2*pi,2*pi);
            
            Cmask1=mask1.*exp(1i*mask2);
            Cmask2=mask3.*exp(1i*mask4);
            
            E_y=mask1.*im;
            E_x=mask3;
            Ph_y=mask2;
            Ph_x=mask4;
            
            
        
        
        case 32
            beamNAME='Encryption WhiteNoise';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            
            mask1=zeros(SLM_size);
            mask2=zeros(SLM_size);
            mask3=zeros(SLM_size);
            mask4=zeros(SLM_size);
            
            RandStream.setDefaultStream(RandStream('mt19937ar','seed',0)); %seed
            stk=32; %pixels for rand number  (power of two)
            
            key1=scripts.normalize_2D(rand(SLM_size/stk));
            key2=scripts.normalize_2D(rand(SLM_size/stk));
            key3=scripts.normalize_2D(rand(SLM_size/stk));
            key4=scripts.normalize_2D(rand(SLM_size/stk));
            
            for i=0:SLM_size(1)/stk-1
                for j=0:SLM_size(2)/stk-1
                    mask1(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key1(i+1,j+1);
                    mask2(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key2(i+1,j+1);
                    mask3(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key3(i+1,j+1);
                    mask4(i*stk+1:i*stk+stk,j*stk+1:j*stk+stk)=key4(i+1,j+1);
                end
            end
            
            Cmask1=mask1.*exp(1i*2*pi*mask2);
            Cmask2=mask3.*exp(1i*2*pi*mask4);
            
            E_y=mask1.*im;
            E_x=mask3;
            Ph_y=mask2*2*pi;
            Ph_x=mask4*2*pi;
                    
        case 31
            beamNAME='Encryption Clear';
            currentPATH=cd;
            cd ..
            phdPATH=cd;
            cd(currentPATH);
            imgPATH=[phdPATH '/WorkInProgress/OpticalEncryption/IMAGE.png'];
            im=im2double(imread(imgPATH))';
            E_y=im;
            E_x=ones(SLM_size);
            Ph_x=zeros(SLM_size);
            Ph_y=Ph_x;
        
        case 30
            beamNAME='HammerCircular_FPB';
            
            fx=60/sqrt(2)/1.5;fy=fx;
            
            aux=x;
            x=y/fx;
            y=aux/fy;
            
            S0=(sqrt(x.^2/8+y.^2/2)<=1).*1;
            elip=(sqrt(x.^2/8+y.^2/2)>1).*10;

            chi0=0;
            psi0=pi/2;
            
            z=sqrt(1-(x/4).^2-(y/2).^2).*S0;
            chi2 =( chi0 + asin( z.*y ) ).*S0;
            psi2 =( psi0 + 2*atan2(z.*x,2*(2*z.^2-1)) ).*S0;
            
            S1=cos(psi2).*cos(chi2).*S0;
            S2=sin(psi2).*cos(chi2).*S0;
            S3=sin(chi2).*S0;
            
            E_x=sqrt((S0+S1).*S0/2);
            E_y=sqrt((S0-S1).*S0/2);
            Ph_y=atan2(S3,S2)+pi/2;
            Ph_x=zeros(SLM_size);
            
        
        case 29
            beamNAME='HammerLineal_FPB';
        
            fx=60/sqrt(2)/1.5;fy=fx;
            
            aux=x;
            x=y/fx;
            y=aux/fy;
            
            S0=(sqrt(x.^2/8+y.^2/2)<=1).*1;
            elip=(sqrt(x.^2/8+y.^2/2)>1).*10;

            chi0=0;
            psi0=0;
            
            z=sqrt(1-(x/4).^2-(y/2).^2).*S0;
            chi2 =( chi0 + asin( z.*y ) ).*S0;
            psi2 =( psi0 + 2*atan2(z.*x,2*(2*z.^2-1)) ).*S0;
            
            S1=cos(psi2).*cos(chi2).*S0;
            S2=sin(psi2).*cos(chi2).*S0;
            S3=sin(chi2).*S0;
            
            E_x=sqrt((S0+S1).*S0/2);
            E_y=sqrt((S0-S1).*S0/2);
            Ph_y=atan2(S3,S2);
            Ph_x=zeros(SLM_size);
        
        case 28
            beamNAME='MollweideProjection_FPB';
            
            fx=60/sqrt(2)/1.5;fy=fx;
            
            aux=x;
            x=y/fx;
            y=aux/fy;
            
            psi0=0;
            S0=(sqrt(x.^2/8+y.^2/2)<=1).*1;
            elip=(sqrt(x.^2/8+y.^2/2)>1).*10;

            thetaM=asin(y/sqrt(2));
            chi2 = asin( (2*thetaM+sin(2*thetaM))/pi ).*S0;
            psi2 = psi0 + pi*x./2/sqrt(2)./cos(thetaM).*S0;
            
            S1=cos(psi2).*cos(chi2).*S0;
            S2=sin(psi2).*cos(chi2).*S0;
            S3=sin(chi2).*S0;
            
            E_x=sqrt((S0+S1).*S0/2);
            E_y=sqrt((S0-S1).*S0/2);
            Ph_y=atan2(S3,S2);
            Ph_x=zeros(SLM_size);
            
        case 27
            beamNAME='AlonsoMillorat_FPB';
            
            R=60;
            circ=(x.^2+y.^2<R^2).*1;
            rho=sqrt(x.^2+y.^2);
            
            Ex=cos(rho/R*pi/2).*circ;
            Ey=sin(rho/R*pi/2).*exp(1i*theta).*circ;
            
            E_x=abs(Ex);E_y=abs(Ey);
            Ph_x=angle(Ex);Ph_y=angle(Ey);
            
        case 26
            
            beamNAME='AlonsoCircular_FPB';
            
            elip=zeros(size(x));

            z=0;
            lambda=1;
            k=2*pi/lambda;
            A=1;
            W0=60;
            Zr=pi/lambda*W0^2;
            G=1+1i*z/Zr;

            U00=A/G*exp(1i*k*z-rho.^2/W0^2/G);
            U01=rho.*(sqrt(2)/W0/G).*exp(1i*theta).*U00;

            Ex=(U00+U01)/2;
            Ey=(U01-U00)*1i/2;

            As=-10;
            qua=1;%exp(1i*As*rho.^2);

            Ph_x=angle(Ex)+angle(qua);Ph_y=angle(Ey)+angle(qua);
            E_x=abs(Ex);E_y=abs(Ey);
        
        case 25
            beamNAME='AlonsoLineal_FPB';
            
            rho=sqrt(x.^2+y.^2);
            z=0;
            lambda=1;
            k=2*pi/lambda;
            A=1;
            W0=50;
            Zr=pi/lambda*W0^2;
            G=1+1i*z/Zr;
            
            U00=A/G*exp(1i*k*z-rho.^2/W0^2/G);
            U01=rho.*(sqrt(2)/W0/G).*exp(1i*theta).*U00;
            
            Ex=U00;
            Ey=U01;
            
            
            
            E_x=abs(Ex);E_y=abs(Ey);
            Ph_x=angle(Ex);Ph_y=angle(Ey);

        case 24
            beamNAME='DonutF';
            
            R=100;%/sqrt(2);
            %f0=2;
            NA=0.85;
            rho=sqrt((x.^2+y.^2)/R^2);
            phi=theta;
            theta2=asin(rho);
            %g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi);%.*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi);%.*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y); 
            
        
        case 23
            beamNAME='DonutF';
            
            R=100/sqrt(2);f0=2;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
            g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi).*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi).*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y); 
        
        case 22
            beamNAME='DonutF';
            
            R=100/sqrt(2);f0=1.5;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
            g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi).*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi).*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y); 
        
        case 21
            beamNAME='DonutF';
            
            R=100/sqrt(2);f0=1;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
            g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi).*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi).*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y); 
        
        
        case 20
            beamNAME='DonutF';
            
            R=100/sqrt(2);f0=0.8;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
            g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi).*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi).*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y); 
        
        case 19
            beamNAME='DonutF';
            
            R=100/sqrt(2);f0=0.6;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
            g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi).*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi).*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y);
        
        case 18
            beamNAME='DonutF';
            
            R=100/sqrt(2);f0=0.4;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
            g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi).*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi).*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y);
        
        case 17
            beamNAME='DonutF';
            
            R=100/sqrt(2);f0=0.2;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
            g=exp(-(rho/f0).^2);        
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi).*g;
            E_y=1i*cos(theta2).*cos(phi)+sin(phi).*g;
            pp=rho>1;
            E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y);
            
        case 16
            beamNAME='DonutF';
            
            R=100/sqrt(2);
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            theta2=asin(rho);
                
            E_x=-1i*cos(theta2).*sin(phi)+cos(phi);
            E_y=1i*cos(theta2).*cos(phi)+sin(phi);
     %       pp=rho>1;
      %      E_x(pp)=0;E_y(pp)=0;
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y);
            
    %-------------------------------------        
            
        
        case 4 %LG01
            beamNAME='LG01';
            R_x=(abs(cos(theta)));
            R_y=(abs(sin(theta)));
            R_Phx=angle(cos(theta));
            R_Phy=angle(sin(theta));
    %       R_Ph=((theta>pi/2&theta<pi)|(theta<0&theta>-pi/2))*pi;%zeros(SLM_size);
    %       case 'A'  %Azimutal-polarized beam
            A_x=(abs(sin(theta)));
            A_y=(abs(cos(theta)));
            A_Phx=angle(sin(theta));
            A_Phy=angle(cos(theta))+pi;

            R=70;
            rho=2*(x.^2+y.^2)/R^2;
    %       Gaussian=exp(-rho/2);
    %       D=rho.*cos(theta).^2.*exp(-rho) + ... 
    %         rho.^2.*cos(2*(theta+2*pi/4)).^2.*exp(-rho);% + ...
           % rho.^3.*cos(3*(theta+2*pi/12)).^2.*exp(-rho);
            D=abs((-rho+1).*exp(-rho/2));
            core=((x).^2+(y).^2<(512-462)^2).*1;
            clad=1-core;


            E_x=D.*core.*R_x+D.*clad.*A_x*2.5;%Gaussian;
            E_y=D.*core.*R_y+D.*clad.*A_y*2.5;%Gaussian;
            Ph_x=core.*R_Phx+clad.*A_Phx;
            Ph_y=core.*R_Phy+clad.*A_Phy;

        
        case 5 %O1
            beamNAME='O1';
            R=70;
            rho=2*(x.^2+y.^2)/R^2;
            A_0=sqrt(rho.^2.*exp(-rho));
            A_0=A_0/max(max(A_0));
            k=0;l=1;theta_amp=pi/4;theta_ph=7*pi/4;
            E_x=cos(theta*k+theta_amp);
            E_y=sin(theta*k+theta_amp);
            %Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=A_0.*abs(E_x);E_y=A_0.*abs(E_y);
        
            Ph=mod(l*theta+theta_ph,2*pi);
            Ph_x=zeros(SLM_size);
            px=find(Ph>pi);
            Ph_x(px)=-Ph(px);
            Ph_y=zeros(SLM_size);
            py=find(Ph<=pi);
            Ph_y(py)=Ph(py);
        
        case 6 %O2
            beamNAME='O2';
            R=70;
            rho=2*(x.^2+y.^2)/R^2;
            A_0=sqrt(rho.^2.*exp(-rho));
            A_0=A_0/max(max(A_0));
            k=1;l=1;theta_amp=0;theta_ph=0;
        
            E_x=cos(theta*k+theta_amp);
            E_y=sin(theta*k+theta_amp);
            %Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=A_0.*abs(E_x);E_y=A_0.*abs(E_y);
        
            %Ph=mod(l*theta+theta_ph,2*pi);
            Ph_x=theta+pi/2;
            px=find(theta>pi/4);
            Ph_x(px)=pi-theta(px);
            px=find(theta>3*pi/4);
            Ph_x(px)=-pi/2+theta(px);
            px=find(theta>5*pi/4);
            Ph_x(px)=-theta(px);
            px=find(theta>7*pi/4);
            Ph_x(px)=pi/2+theta(px);
        
            Ph_y=pi/2-theta;
            py=find(theta>pi/4);
            Ph_y(py)=theta(py);
            py=find(theta>3*pi/4);
            Ph_y(py)=-pi/2-theta(py);
            py=find(theta>5*pi/4);
            Ph_y(py)=pi+theta(py);
            py=find(theta>7*pi/4);
            Ph_y(py)=pi/2-theta(py);
        
        case 7 %O1 well done
            beamNAME='O1 well done';
            R=70;
            rho=2*(x.^2+y.^2)/R^2;
            A_0=sqrt(rho.^2.*exp(-rho));
            A_0=A_0/max(max(A_0));
            k=1/2;l=1;theta_amp=0;theta_ph=0;
        
            E_x=abs(sin(theta*k+theta_amp));
            E_y=abs(cos(theta*k+theta_amp));
            %Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=A_0.*abs(E_x);E_y=A_0.*abs(E_y);
        
            Ph_x=zeros(size(theta));
            Ph_x(l*theta<pi)=l*theta(l*theta<pi);
            Ph_x(l*theta>pi)=pi;
        
            Ph_y=zeros(size(theta));
            Ph_y(l*theta>pi)=l*theta(l*theta>pi)+pi;
            Ph_y(l*theta<pi)=pi;
        
        case 8 %O2 wrong done
            beamNAME='O2 wrong done';
            R=70;
            rho=2*(x.^2+y.^2)/R^2;
            A_0=sqrt(rho.^2.*exp(-rho));
            A_0=A_0/max(max(A_0));
            k=1;l=1;theta_amp=0;theta_ph=0;
        
            E_x=cos(theta*k+theta_amp);
            E_y=sin(theta*k+theta_amp);
            %Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=A_0.*abs(E_x);E_y=A_0.*abs(E_y);
        
            %Ph=mod(l*theta+theta_ph,2*pi);
            Ph_x=zeros(size(theta));
            px=find(theta>pi/2);
            Ph_x(px)=2*theta(px)-pi;
            px=find(theta>3*pi/4);
            Ph_x(px)=-2*theta(px);
            px=find(theta>pi);
            Ph_x(px)=0;
            px=find(theta>3*pi/2);
            Ph_x(px)=2*theta(px)-pi;
            px=find(theta>7*pi/4);
            Ph_x(px)=-2*theta(px);
        
            Ph_y=2*theta;
            py=find(theta>pi/4);
            Ph_y(py)=-2*theta(py)+pi;
            py=find(theta>pi/2);
            Ph_y(py)=0;
            py=find(theta>pi);
            Ph_y(py)=2*theta(py);
            py=find(theta>5*pi/4);
            Ph_y(py)=-2*theta(py)+pi;
            py=find(theta>3*pi/2);
            Ph_y(py)=0;
        
        case 9 %O3 wrong done
            beamNAME='O3 (wrong done)';
            R=70;
            rho=2*(x.^2+y.^2)/R^2;
            A_0=sqrt(rho.^2.*exp(-rho));
            A_0=A_0/max(max(A_0));
            k=1;l=1;theta_amp=0;theta_ph=0;
        
            E_x=1;
            E_y=1;
            %Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=A_0.*abs(E_x);E_y=A_0.*abs(E_y);
        
            %Ph=mod(l*theta+theta_ph,2*pi);
            Ph_x=2*theta+pi/2;
            px=find(theta>pi/4);
            Ph_x(px)=0;
            px=find(theta>3*pi/4);
            Ph_x(px)=2*theta(px)+5*pi/2;
            px=find(theta>5*pi/4);
            Ph_x(px)=0;
            px=find(theta>7*pi/4);
            Ph_x(px)=2*theta(px)+pi/2;
        
            Ph_y=zeros(size(theta));
            py=find(theta>pi/4);
            Ph_y(py)=-2*theta(py)+pi/2;
            py=find(theta>3*pi/4);
            Ph_y(py)=0;
            py=find(theta>5*pi/4);
            Ph_y(py)=-2*theta(py)+pi/2;
            py=find(theta>7*pi/4);
            Ph_y(py)=0;
        
        case 10 %O3 Well done
            beamNAME='O3 (well done)';
        
            R=70;
            rho=2*(x.^2+y.^2)/R^2;
            A_0=sqrt(rho.^2.*exp(-rho));
            A_0=A_0/max(max(A_0));
            k=1;l=1;theta_amp=0;theta_ph=0;
        
            E_x=1;
            E_y=1;
            %Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=A_0.*abs(E_x);E_y=A_0.*abs(E_y);
        
            %Ph=mod(l*theta+theta_ph,2*pi);
            Ph_x=-theta+3*pi/4;
            px=find(theta>3*pi/4);
            Ph_x(px)=theta(px)-3*pi/4;
            px=find(theta>7*pi/4);
            Ph_x(px)=-theta(px)+3*pi/4;
        
            Ph_y=theta+pi/4;
            py=find(theta>3*pi/4);
            Ph_y(py)=-theta(py)-pi/4;
            py=find(theta>7*pi/4);
            Ph_y(py)=theta(py)+pi/4;
        
        case 11 %Elliptical
            beamNAME='Elliptical';
            a1=100;b1=50;
            a2=80*2;b2=20*2;R=200;
            k=1;l=0;theta_amp=0;theta_ph=0;
        
            A_0=((x/a1).^2+(y/b1).^2<1)*1.;
            A_0=A_0.*exp(-(x/a2).^2-(y/b2).^2);
        
            E_x=cos(theta*k+theta_amp);
            E_y=sin(theta*k+theta_amp);
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y);
        
            rho=sqrt(x.^2+y.^2);
            gaussian=exp(-rho.^2/R.^2);
        
            E_x=A_0.*E_x.*gaussian;E_y=A_0.*E_y.*gaussian;
        
        
        
        case 12 %Hermite_Gauss 21
            beamNAME='Hermite-Gauss 21';
            theta=mod(theta,2*pi);
            R=60;
            rho=sqrt(x.^2+y.^2);
            E_x=(8*x.^2/R^2-2).*sqrt(8)./R.*y.*exp(-rho.^2/R^2);
            E_x=E_x/max(max(E_x));
            E_y=(8*x.^2/R^2-2).*sqrt(8)./R.*y.*exp(-rho.^2/R^2);
            E_y=E_y/max(max(E_y));
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y);
        
            Z1=zeros(SLM_size);Z2=zeros(SLM_size);
            Z3=zeros(SLM_size);Z4=zeros(SLM_size);
            Z5=zeros(SLM_size);Z6=zeros(SLM_size);
            Z1(1:477,1:385)=1;Z2(1:477,386:end)=1;
            Z3(478:547,1:385)=1;Z4(478:547,386:end)=1;
            Z5(548:end,1:385)=1;Z6(548:end,386:end)=1;
            Ph_x=zeros(SLM_size);Ph_y=zeros(SLM_size);
            Ph_x=-Z3*pi/2-Z1*pi;Ph_y=Z6*pi-Z4*pi/2;
            
        case 13
            beamNAME='AnellNestor2';
            rho=sqrt(x.^2+y.^2);
            factor=251.25/4;
            Di=factor*1;%mida en milimetres
            Gr=2;
            Rout=round(Di/2+Gr);Rin=round(Di/2-Gr);
            Circ1=(rho<Rout).*1;
            Circ2=(rho<Rin).*1;
            Anell=Circ1-Circ2;
            E_x=Anell;
            E_y=zeros(SLM_size);
            Ph_x=E_y;Ph_y=E_y;
            
        case 14
            beamNAME='AnellNestor';
            rho=sqrt(x.^2+y.^2);
            factor=251.25/4;
            Di=factor*1;%mida en milimetres
            Gr=2;
            Rout=round(Di/2+Gr);Rin=round(Di/2-Gr);
            Circ1=(rho<Rout).*1;
            Circ2=(rho<Rin).*1;
            Anell=Circ1-Circ2;
            E_x=Anell;
            E_y=zeros(SLM_size);
            Ph_x=E_y;Ph_y=E_y;
            
            
        case 15
            beamNAME='DonutF';
            R=100;
            rho=(x.^2+y.^2)/R^2;
            phi=theta;
            R0=R/300;
                
            E_x=rho.*cos(phi).*exp(-rho.^2/R0^2);
            E_y=rho.*sin(phi).*exp(-rho.^2/R0^2);
            Ph_x=angle(E_x);Ph_y=angle(E_y);
            E_x=abs(E_x);E_y=abs(E_y);
            
          
        otherwise
            if beam_type==1 %Radial
                beamNAME='Radial';
                k=1;l=0;theta_amp=0;theta_ph=0;
            elseif beam_type==2 %Azimuthal
                beamNAME='Azimuthal';
                k=1;l=0;theta_amp=-pi/2;theta_ph=0;
            elseif beam_type==3 %Star-like
                beamNAME='Star-like';
                k=4;l=0;theta_amp=0;theta_ph=0;
            else
                beamNAME='Horizontal-Homogenous';
                k=0;l=0;theta_amp=0;theta_ph=0;
            end
        
        E_x=cos(theta*k+theta_amp);
        E_y=sin(theta*k+theta_amp);
        Ph_x=angle(E_x);Ph_y=angle(E_y);
    %     p=find(xor(E_x<0,E_y<0));
        E_x=abs(E_x);E_y=abs(E_y);
    %     Ph=l*theta+theta_ph;
    %     Ph(p)=Ph(p)+pi;

    end

    if gauss_correction~=0
        R = 5;
        f0 = 1;
        igauss = 1.05-exp(-1/f0^2*(x.^2+y.^2).^1/R^2)*0.1;
        % figure;imagesc(igauss,[0 1])
        E_x = E_x.*igauss;
        E_y = E_y.*igauss;
    end
    Ph_x = mod(Ph_x,2*pi);
    Ph_y = mod(Ph_y,2*pi);



    if infile~=0 
        % imwrite() normalizes the images !!! 
        disp(['beam_type=' beamNAME]);

        fNAME=[designPATH num2str(beam_type) '_Ex'];
        dlmwrite([fNAME '.dat'], E_x);
        imwrite(scripts.normalize_2D(E_x), [fNAME '.png'], 'png')

        fNAME=[designPATH num2str(beam_type) '_Phx'];
        dlmwrite([fNAME '.dat'], Ph_x);
        imwrite(scripts.normalize_2D(Ph_x), [fNAME '.png'], 'png')

        fNAME=[designPATH num2str(beam_type) '_Ey'];
        dlmwrite([fNAME '.dat'], E_y);
        imwrite(scripts.normalize_2D(E_y), [fNAME '.png'], 'png')

        fNAME=[designPATH num2str(beam_type) '_Phy'];
        dlmwrite([fNAME '.dat'], Ph_y);
        imwrite(scripts.normalize_2D(Ph_y), [fNAME '.png'], 'png')
    end

    if draw~=0
        I = (abs(E_x)).^2+(abs(E_y)).^2;
        I = I/max(max(I));

        figure
        imagesc(I',[0 1])
        title 'Design-I'

        figure
        imagesc(E_x')
        title 'Design-Ex'

        figure
        imagesc(E_y')
        title 'Design-Ey'

        figure
        imagesc(Ph_x') 
        title 'Design-Ph_x'

        figure
        imagesc(Ph_y') 
        title 'Design-Ph_y'

        figure
        Ph=mod(Ph_y-Ph_x,2*pi);
        imagesc(Ph',[0 2*pi]) 
        title 'Design-Ph'
    end

    fid = fopen('+scripts/design_kinds.txt');
    tline = fgetl(fid);
    beams_in_file = 0;
    while ischar(tline)
        tline = fgetl(fid);
        beams_in_file = beams_in_file+1;
    end
    fclose(fid);

    % add the name in the list 'beam_kinds.txt' if it's not there
    if beam_type>beams_in_file
        dlmwrite('+scripts/design_kinds.txt', [num2str(beam_type) ': ' beamNAME], ...
                 '-append', 'delimiter', '', 'newline', 'pc')
    else
        disp('Designing the beam:');
        fid = fopen('+scripts/design_kinds.txt');
        beamN = textscan(fid, '%s', 1, 'delimiter', ',', 'headerlines', beam_type-1);
        disp(beamN{1});
    end


    if stokes~=0
        [s0,s1,s2,s3,DOP] = scripts.stokes(E_x.*exp(1i*Ph_x), E_y.*exp(1i*Ph_y));
        figure
        imagesc(s0',[-1 1]);title 'S0'
        figure
        imagesc(s1',[-1 1]);title 'S1'
        figure
        imagesc(s2',[-1 1]);title 'S2'
        figure
        imagesc(s3',[-1 1]);title 'S3'
        figure
        imagesc(DOP',[-1 1]);title 'DOP'
    end
