%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%        https://github.com/dmaluenda/DigiHolos2LaserBeamModelation
%
%    David Maluenda Niubo - dmaluendn@gmail.com            CC: by, NC, SA 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transmisions-----------------------------------------
Xi_1=400;Xf_1=800;
Y1_1=250;Y2_1=325;
Y3_1=425;Y4_1=550;
Xi_2=Xi_1;Xf_2=Xf_1;
Y1_2=Y1_1;Y2_2=Y2_1;
Y3_2=Y3_1;Y4_2=Y4_1;

im=double(imread('I1_0.png'));
up=im(Y1_1:Y2_1,Xi_1:Xf_1);
down=im(Y3_1:Y4_1,Xi_1:Xf_1);
fact1=sum(sum(down))/sum(sum(up));
im=double(imread('I2_0.png'));
up=im(Y1_2:Y2_2,Xi_2:Xf_2);
down=im(Y3_2:Y4_2,Xi_2:Xf_2);
fact2=sum(sum(down))/sum(sum(up));
I1(1:256)=0;
I2(1:256)=0;

for i=1:256
    filename=['I1_' num2str(i-1) '.png'];
    im=double(imread(filename));
    up=fact1*im(Y1_1:Y2_1,Xi_1:Xf_1);
    down=im(Y3_1:Y4_1,Xi_1:Xf_1);
    I1(i)=sum(sum(up))/sum(sum(down));
end
for i=1:256
    filename=['I2_' num2str(i-1) '.png'];
    im=double(imread(filename));
    up=fact2*im(Y1_2:Y2_2,Xi_2:Xf_2);
    down=im(Y3_2:Y4_2,Xi_2:Xf_2);
    I2(i)=sum(sum(up))/sum(sum(down));
end

I1=sqrt(I1);
I2=sqrt(I2);
h=figure;
plot(I1)
axis([1 256 0 1.1])
print(h,'-depsc','amplitude_SLM1.eps');
print(h,'-dpng','amplitude_SLM1.png');
h=figure;
plot(I2)
axis([1 256 0 1.1])
print(h,'-depsc','amplitude_SLM2.eps');
print(h,'-dpng','amplitude_SLM2.png');