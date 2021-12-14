%% 2D animáció
clear all;

for alfa1 = -90:30:90


ss=0.05; %términtavételi felbontás [m]

L=10;   %pályamáret
[x,y]=meshgrid(0:ss:L); %x és y is -3-tól 3-ig megy 1-esével


% alfa1=0;    %irányszög [deg]
alfa1=30;
xs1=-500*ss*L*sind(alfa1);   %hangforrás pozíciója
ys1=-500*ss*L*cosd(alfa1);

alfa2=-60;
xs2=-500*ss*L*sind(alfa2);   %hangforrás pozíciója
ys2=-500*ss*L*cosd(alfa2);

rLim=0.5; %távlimit
r1=sqrt((x-xs1).^2+(y-ys1).^2);
r1(r1<rLim)=rLim;
r2=sqrt((x-xs2).^2+(y-ys2).^2);
r2(r2<rLim)=rLim;

c=2;    %hullámsebesség
micNum= 5;  %MICNUM


% jeltulajdonságok
Amp=20;

Amp1=1;
lambda=2;
k=2*pi/lambda;
w=2*pi/(lambda/c);

Amp2=0.8;
lambda2=3.141592;
k2=2*pi/lambda2;
w2=2*pi/(lambda2/c);

Amp3=0;
lambda3=2;
k3=2*pi/lambda3;
w3=2*pi/(lambda3/c);



d=1; % Antennasor távolsága [m] (1 tizedes jegy!)


ts=0.001; %mintavételi idő  [sec]
T=20;    %időablak szélessége [sec]

rec=[]; %vett jelmátrix (oszlopok a jelsorok)

    
for t=0:ts:T

    %MINTAVÉTEL
    recTemp=[];
    for m=1:micNum
        rRed=sqrt((ys1-L)^2+(xs1-(m-1)*d)^2);
        rRed2=sqrt((ys2-L)^2+(xs2-(m-1)*d)^2);
        fRed=Amp*real(Amp1*exp(1i*(k.*rRed-t*w))./rRed + Amp2*exp(1i*(k2.*rRed-t*w2))./rRed + Amp3*exp(1i*(k3.*rRed2-t*w3))./rRed2); %jelforrás
        recTemp=[recTemp fRed];    %mintavétel
    end
    
    rec=[rec; recTemp];

end

% for m=1:micNum
%     figure(22)
%     hold on;
%     plot(rec(:,m));
% end

% DELAY AND SUM
t=0:ts:T;
x=0:1:T/ts;

o=-1;
delSumGraph=[];


Theta=-pi/2:0.05:pi/2;
for k=1:length(Theta)
    
    %ThetaD=-40;
    %Theta=pi*ThetaD/180;
    %floor(sin(Theta)*2/v/0.05)
    

    delVal=[];  %mikrofonok késleltetési vektora
    for m=1:micNum
        delVal=[delVal -floor(sin(Theta(k))*d*(m-1)/c/ts)];
    end
    
    o=delVal(2);
    
    del=[]; %késleletetett jelmátrix
    %if flag~=o
        for m=1:micNum
            delTemp=circshift(rec(:,m)',delVal(m));
            %delTemp=[zeros(1, delVal(m)+3*length(rec(:,m))) rec(:,m)' zeros(1, -delVal(m)+3*length(rec(:,m)))];
            del=[del; delTemp];
        end
        
        
    delSum=sum(del);
    
    delSumGraph=[delSumGraph delSum*delSum'];
    ref=rec(:,1)'*rec(:,1);
    o=delVal(1);

%         polarplot(Theta(1:k)+pi/2, 2*sqrt(delSumGraph)/T);
%         axis([0 180 0 max(2*sqrt(delSumGraph)/T)])
% %         polarplot(Theta(1:k)+pi/2, log(0.01+delSumGraph/max(delSumGraph)));
% %         axis([0 180 -2 max(log(0.01+delSumGraph/max(delSumGraph)))])
        
        pause(0.01);
        
       
end
        figure(8);
        polarplot(Theta(1:k)+pi/2, sqrt(delSumGraph/ref));
        axis([0 180 0 max(sqrt(delSumGraph/ref))])
        hold on
end