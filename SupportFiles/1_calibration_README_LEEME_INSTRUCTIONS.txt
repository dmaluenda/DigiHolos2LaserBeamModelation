%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Claibration method and system set up:

1: Take the images from Mach-Zehnder interferometer by using the VI_LabView 'SLM_Calibrator_AUTO.vi' following the indications.
2: Analize the images by using the matlab program 'caracterizer.m'.
3: Smooth the curve with 'smoother.m' progrma as times as needed.
4: Generate the map files by using 'map_generator.m' program.
5: Generate the files, whose have the posibles complex values and how it can be obtained by using 'complex_values.m' program.

Now, ComplexValues_SLMx.txt have been made and they must be in the CurveCaract folder.

6: Run the 'TestPaternsGenerator.m' to generate the general holograms like 'SwitchOff.png' and plain complex values like 'SLM1_Tx_phX.png'.