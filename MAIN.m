%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Hologram generator for a Mach-Zehnder interferometer with
% orthogonal-polarized arms. (2012)


clear variables; close all
PATH=cd;

%%%TO MAKE THE GENERAL PATTERNS for the LabView PROGRAMS%%%
% cd SupportFiles
% TestPaternsGenerator
% cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Modulation mode -------------------- 
Modulation = 'Amplitude';


%---Label of the beam (1:Radial,2:Athimutal,3-Star,...)
beam_types = 86; % it can be a vector
for beam_type=beam_types

    %---Aditional-options: (0)NO - (else)YES  --------- 
    makeHolo = 1;

    %[2 4 44] % (2) for a simple Arizon's procedure and (4) for amplied 
    for macropixel = 2;

        makeDesign   = 1;  % makes the design instead just taking a given file
        designINfile = 0;  % writes the design in a file
        simulation   = 0;  % simulates the reconstructed paraxial beam
        stokes       = 0;  % gets the stokes images of the simulated beam
        corrected    = 0;  % adds a phase map to correct the beam (Not aviable yet!!)


        %---Parameters----------------------
        SLM_size = [1024 768];
        D_SLM    = 250;  % Spot's Diameter at SLM [px]
        T_max    = 0.5;  % Max transmittance 

        A = [500 14]  ; B = [892 378];  % CCD Key Points [px]
        C = [480 770] ; D = [98 392];   % Just for a corrected holograms
        
        tic
        options = '';  % empty sting to initialize or any other option to add
        
        [~,count] = max(beam_types==beam_type);
        disp(['The HoloGenerator has been started (' num2str(count) ...
              '/' num2str(max(size(beam_types))) ')'])

        % --- Creating/Importing the beam's design -------
        if(makeDesign==1)
            [h_T1,h_T2,h_ph1,h_ph2] = scripts.beam_design(SLM_size, beam_type, ...
                                                          designINfile, 0, 0);
        else
            h_T1  = dlmread(['Designs/design' num2str(beam_type) '_Ex.dat']);
            h_T2  = dlmread(['Designs/design' num2str(beam_type) '_Ey.dat']);
            h_ph1 = dlmread(['Designs/design' num2str(beam_type) '_Phx.dat']);
            h_ph2 = dlmread(['Designs/design' num2str(beam_type) '_Phy.dat']);
        end

        %---GENERATE-the-HOLOGRAM------------------------------------------------
        if corrected~=0
        % Option: ----correcting-beam--- (NOT aviable yet!)
            % scripts.beam_caracterizer; %Produce: ARM_1, ARM_2 [0 255] and PHASE [0 2pi]
            % scripts.CCD2SLM; %Produce: A_1, A_2 [0 255] and intrinsic_phase [0 2pi]
            % Phase_c=mod(-intrinsic_phase,2*pi); %[0 2pi]
            % [c_T1,c_T2]=scripts.amplitude_corrector(A_1,A_2); %[0 1] %
            % Amp1_c=ones(SLM_size);Amp2_c=ones(SLM_size);
            c_T1  = 1;
            c_T2  = 1;
            c_ph1 = 1;
            c_ph2 = 1;
            h_T1  = h_T1.*c_T1;
            h_T2  = h_T2.*c_T2;
            h_ph1 = h_ph1.*c_ph1;
            h_ph2 = h_ph2.*c_ph2;

            options = strcat(options,'_c');
        end


        %---Makes-Shows-and-SAVES-the-holograms------------------------------------------
        if makeHolo~=0
            if macropixel==2
                [hologram_SLM1,hologram_SLM2] = scripts.mapa_Holo(h_T1, h_T2, ...
                                                                  h_ph1, h_ph2, ...
                                                                  Modulation);
            elseif macropixel==4
                [hologram_SLM1,hologram_SLM2] = scripts.mapa_Holo4(h_T1, h_T2, ...
                                                                   h_ph1,h_ph2);
                options = strcat(options,'_M4');
            elseif macropixel==44
                [hologram_SLM1,hologram_SLM2] = scripts.mapa_Holo44(h_T1, h_T2, ...
                                                                    h_ph1, h_ph2);
                options=strcat(options,'_M44');
            end
            %just for check matrix
            % [DRn_Hol1,DRt_Hol1] = scripts.dynamic_range(hologram_SLM1);
            % [DRn_Hol2,DRt_Hol2] = scripts.dynamic_range(hologram_SLM2);
            % figure;
            % imshow(hologram_SLM1') %,'colormap',colormap)
            % figure;
            % imshow(hologram_SLM2') %,'colormap',colormap)
            switch Modulation
                case 'Amplitude'
                    imwrite(hologram_SLM1', [PATH '/Holograms_AmplitudeOnly/Hologram1_' ...
                                             num2str(beam_type) options '.bmp'],'bmp')
                    imwrite(hologram_SLM2', [PATH '/Holograms_AmplitudeOnly/Hologram2_' ...
                                             num2str(beam_type) options '.bmp'],'bmp')
                otherwise
                    imwrite(hologram_SLM1', [PATH '/Holograms/Hologram1_' ...
                                             num2str(beam_type) options '.bmp'],'bmp')
                    imwrite(hologram_SLM2', [PATH '/Holograms/Hologram2_' ...
                                             num2str(beam_type) options '.bmp'],'bmp')
            end
        end
        %--------------------------------------------------------------------------


        %Option:---Simulates-the-reconstructed-paraxial-beam---
        if simulation==1
            SIM_PATH = [PATH '\Simulations']
            [Esim_x,Esim_y] = scripts.holo_simulator(hologram_SLM1, hologram_SLM2);

            imwrite(scripts.normalize_2D(abs(Esim_x))',...
                [SIM_PATH '\Sim_Ex_' num2str(beam_type) options '.png'],'png')
            imwrite(scripts.normalize_2D(abs(Esim_y))',...
                [SIM_PATH '\Sim_Ey_' num2str(beam_type) options '.png'],'png')
            imwrite(scripts.normalize_2D(mod(angle(Esim_x),2*pi))',...
                [SIM_PATH '\Sim_Phx_' num2str(beam_type) options '.png'],'png')
            imwrite(scripts.normalize_2D(mod(angle(Esim_y),2*pi))',...
                [SIM_PATH '\\Sim_Phy_' num2str(beam_type) options '.png'],'png')

            sim_ph = mod(angle(Esim_y), 2*pi) - mod(angle(Esim_x), 2*pi);
            imwrite(scripts.normalize_2D(mod(sim_ph,2*pi))',...
                [SIM_PATH '\Sim_Ph_' num2str(beam_type) options '.png'],'png')
        
        %Option:---gets-the-stokes-images-of-the-simulated beam---
            if stokes==1
                [s0,s1,s2,s3] = scripts.stokes(Esim_x, Esim_y);

                imwrite(normalize_2D(s0)',...
                    [SIM_PATH '/S0_' num2str(beam_type) options '.png'],'png')
                imwrite(normalize_2D(s1)',...
                    [SIM_PATH '/S1_' num2str(beam_type) options '.png'],'png')
                imwrite(normalize_2D(s2)',...
                    [SIM_PATH '/S2_' num2str(beam_type) options '.png'],'png')
                imwrite(normalize_2D(s3)',...
                    [SIM_PATH '/S3_' num2str(beam_type) options '.png'],'png')
            end
        end
            
        toc
        disp(' ');
    end
end