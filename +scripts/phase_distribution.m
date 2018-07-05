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


%To distribute a phase mask ph0 in two masks ph1 and ph2
function[ph1,ph2,T1,T2]=phase_distribution(ph,T1,T2)
	%ph_ref = pi;%4.0310;
	cal1 = 0*pi/180; %ph0_1
	cal2 = 0*pi/180; %ph0_2

	N = size(ph);
	% d = max(N)-min(N);
	% ph10 = zeros(max(N),max(N));
	% ph20 = zeros(max(N),max(N));
	% % T2o = zeros(max(N),max(N));
	% ph10 = zeros(N);
	% ph20 = zeros(N);

	ph = mod(ph,2*pi);
	P = ph>=pi;
	p = find(P);

	phP = pi/2+ph/2;
	phM = pi/2-ph/2;
	ph1 = phM;
	ph2 = phP;

	ph1(p) = 2*pi-phP(p);
	ph2(p) = -phM(p);

	ph1 = ph1+cal1;
	ph2 = ph2+cal2;

	% for i=1:N(1)
	%     for j=1:N(2)
	% %         T2o(i,j+d/2)=T2(i,j);
	%         if ph0(i,j)<pi
	%             ph1(i,j)=pi/2-ph(i,j)/2+cal1;%3.8990;
	%             ph2(i,j)=pi/2+ph(i,j)/2+cal2; %(i,j+d/2)
	%         else
	%             ph1(i,j)=-pi/2-ph(i,j)/2+cal1;%3.8990-ph(i,j)+ph_ref;
	%             ph2(i,j)=-pi/2+ph(i,j)/2+cal2;
	%         end
	%     end
	% end

	% ph1=ph10;
	% ph2=ph20;


	% ph1=ph10(1:N(1),1:N(2));
	% ph1(:,1:768)=ph1(:,768:-1:1);
	% T1(:,1:768)=T1(:,768:-1:1);
	% 
	% ph20=rotate(ph20);
	% T2=rotate(T2o);
	% ph2=ph20(1:N(1),d/2:d/2+N(2)-1);
	% T2=T2(1:N(1),d/2:d/2+N(2)-1);
	% ph2(:,1:768)=ph2(:,768:-1:1);
	% T2(:,1:768)=T2(:,768:-1:1);
