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

%% Carecterize the beam after the Mach-Zehnder interferometer

% 1/4wave + Analizer
INTERFERENCE_1  = double(imread('REF_PHASE1_CCD.png'))'; 
[DRn_IF,DRt_IF] = dynamic_range(INTERFERENCE_1);

% 1/4wave + Analizer(+45º)
INTERFERENCE_2    =double(imread('REF_PHASE2_CCD.png'))';
[DRn_IF2,DRt_IF2] = dynamic_range(INTERFERENCE_2);

% Just the first ARM
ARM_1 = double(imread('REF_ARM1_CCD.png'))'; 
[DRn_A1,DRt_A1] = dynamic_range(ARM_1);

% Just the secon ARM
ARM_2 = double(imread('REF_ARM2_CCD.png'))'; 
[DRn_A2,DRt_A2] = dynamic_range(ARM_2);

% INTENSITY profile
INTENSITY = (ARM_1+ARM_2)*2;  
[DRn_ITY,DRt_ITY] = dynamic_range(INTENSITY);

% Phase definition from the images
PHASE = normalize_phase(INTERFERENCE_1,INTERFERENCE_2,INTENSITY); 
[DRn_PH,DRt_PH] = dynamic_range(PHASE);

% Returns a (square) 90º rotated for SLM_1 (and NOT flipped image) [0 2pi]
PHASE = rotate_CCD(PHASE); 
ARM_1 = rotate_CCD(ARM_1);
ARM_2 = rotate_CCD(ARM_2);

% New size of image
CCD_size = size(PHASE); 
