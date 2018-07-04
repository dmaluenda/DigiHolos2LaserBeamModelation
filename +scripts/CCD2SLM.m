%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubó - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% David Maluenda (2011-2012)
%
% Resizes the information (ARM_1, ARM_2 and PHASE) from CCD plane to
% the SLMs plane keeping the total size (#pixels at CCD = #pixels at SLM)
% because the spot size on CCD is different than the spot size on SLMs


C_CCD = round([(B(1)+D(1))/2 (C(2)+A(2))/2]); %Spot's Center at CCD
C_SLM = round(SLM_size/2); %Spot's Center at SLM (SLMs can scroll through the plane)
D_CCD = (norm(B-D)+norm(C-A))/2; %Spot's Diameter at CCD
step  = D_SLM/D_CCD; %Conversion factor [mm]/[px]
O     = round(C_SLM+C_CCD*AD*step-CCD_size*AD*step);%Main corner of CCD in SLM as Origen (top-left in SLM)

intrinsic_phase = zeros(SLM_size);
A_1 = zeros(SLM_size);
A_2 = zeros(SLM_size);

i = 1:CCD_size(1);
j = 1:CCD_size(2);

intrinsic_phase(round(O(1)+i*step), round(O(2)+j*step)) = PHASE(i,j);
A_1(round(O(1)+i*step), round(O(2)+j*step)) = ARM_1(i,j);
A_2(round(O(1)+i*step), round(O(2)+j*step)) = ARM_2(i,j);

[DRn_iPh,DRt_iPh]=dynamic_range(intrinsic_phase); %just a control