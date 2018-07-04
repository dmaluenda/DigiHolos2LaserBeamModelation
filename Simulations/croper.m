%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for type=1:23

fileNAME1=['Sim_Ex_' num2str(type) '.png'];
fileNAME2=['Sim_Ey_' num2str(type) '.png'];
fileNAME3=['Sim_Ph_' num2str(type) '.png'];
fileNAME4=['Sim_Phx_' num2str(type) '.png'];
fileNAME5=['Sim_Phy_' num2str(type) '.png'];

fileNAME14=['Sim_Ex_' num2str(type) '_M4.png'];
fileNAME24=['Sim_Ey_' num2str(type) '_M4.png'];
fileNAME34=['Sim_Ph_' num2str(type) '_M4.png'];
fileNAME44=['Sim_Phx_' num2str(type) '_M4.png'];
fileNAME54=['Sim_Phy_' num2str(type) '_M4.png'];

fileNAME144=['Sim_Ex_' num2str(type) '_M44.png'];
fileNAME244=['Sim_Ey_' num2str(type) '_M44.png'];
fileNAME344=['Sim_Ph_' num2str(type) '_M44.png'];
fileNAME444=['Sim_Phx_' num2str(type) '_M44.png'];
fileNAME544=['Sim_Phy_' num2str(type) '_M44.png'];

im1=imread(fileNAME1);
im2=imread(fileNAME2);
im3=imread(fileNAME3);
im4=imread(fileNAME4);
im5=imread(fileNAME5);

im14=imread(fileNAME1);
im24=imread(fileNAME2);
im34=imread(fileNAME3);
im44=imread(fileNAME4);
im54=imread(fileNAME5);

im144=imread(fileNAME144);
im244=imread(fileNAME244);
im344=imread(fileNAME344);
im444=imread(fileNAME444);
im544=imread(fileNAME544);

IM1=im1(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM2=im2(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM3=im3(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM4=im4(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM5=im5(768/2-128:768/2+128,1024/2-128:1024/2+128);

IM14=im14(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM24=im24(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM34=im34(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM44=im44(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM54=im54(768/2-128:768/2+128,1024/2-128:1024/2+128);

IM144=im144(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM244=im244(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM344=im344(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM444=im444(768/2-128:768/2+128,1024/2-128:1024/2+128);
IM544=im544(768/2-128:768/2+128,1024/2-128:1024/2+128);

fileNAME1=['cSim_Ex_' num2str(type) '.png'];
fileNAME2=['cSim_Ey_' num2str(type) '.png'];
fileNAME3=['cSim_Ph_' num2str(type) '.png'];
fileNAME4=['cSim_Phx_' num2str(type) '.png'];
fileNAME5=['cSim_Phy_' num2str(type) '.png'];

fileNAME14=['cSim_Ex_' num2str(type) '_M4.png'];
fileNAME24=['cSim_Ey_' num2str(type) '_M4.png'];
fileNAME34=['cSim_Ph_' num2str(type) '_M4.png'];
fileNAME44=['cSim_Phx_' num2str(type) '_M4.png'];
fileNAME54=['cSim_Phy_' num2str(type) '_M4.png'];

fileNAME144=['cSim_Ex_' num2str(type) '_M44.png'];
fileNAME244=['cSim_Ey_' num2str(type) '_M44.png'];
fileNAME344=['cSim_Ph_' num2str(type) '_M44.png'];
fileNAME444=['cSim_Phx_' num2str(type) '_M44.png'];
fileNAME544=['cSim_Phy_' num2str(type) '_M44.png'];

imwrite(fileNAME1);
imwrite(fileNAME2);
imwrite(fileNAME3);
imwrite(fileNAME4);
imwrite(fileNAME5);

im14=imwrite(fileNAME1);
im24=imwrite(fileNAME2);
im34=imwrite(fileNAME3);
im44=imwrite(fileNAME4);
im54=imwrite(fileNAME5);

im144=imwrite(fileNAME144);
im244=imwrite(fileNAME244);
im344=imwrite(fileNAME344);
im444=imwrite(fileNAME444);
im544=imwrite(fileNAME544);


end