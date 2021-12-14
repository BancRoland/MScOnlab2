%% 2D animáció

clear all;
close all;

TempMap=[];
for i=0:0.01:1
   TempMap=[TempMap; 0+i 0+i 1];
end
for i=0.01:0.01:1
   TempMap=[TempMap; 1 1-i 1-i];
end


%avimaker
v = VideoWriter('xxx');
open(v);


ss=0.1; %def: 0.05 términtavételi felbontás [m]

L=20;   %pályamáret
[x,y]=meshgrid(0:ss:L); %x és y is -3-tól 3-ig megy 1-esével

% alfa1=90;    %irányszög [deg]
% xs1=5*ss*L*sind(alfa1);   %hangforrás pozíciója
% ys1=5*ss*L*cosd(alfa1);

xs1=0;
ys1=10;
% xs2=7.5;
% ys2=10;
vx=0.2;
vy=0;

% alfa2=-60;
% xs2=-100*ss*L*sind(alfa2);   %hangforrás pozíciója
% ys2=-100*ss*L*cosd(alfa2);

rLim=0.2; %távlimit

% r2=sqrt((x-xs2).^2+(y-ys2).^2);
% r2(r2<rLim)=rLim;

c=1;    %hullámsebesség
micNum= 3;  %MICNUM


% jeltulajdonságok
Amp=1;

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



d=3; % Antennasor távolsága [m] (1 tizedes jegy!)


ts=0.1; %0.05 az alap   mintavételi idő  [sec]
T=15;    %időablak szélessége [sec]

rec=[]; %vett jelmátrix (oszlopok a jelsorok)

figure(1);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1])



t=0:ts:T;
signal=sin(w*t);
impulse=zeros(length(x(1,:)));
impulse(:,:,200)=0;

%pályaleíráshoz
vc=2;
rc=2;
wc=vc/rc;
% xs1v=[ 1:0.03:2 2:0.06:4 4:0.09:6 6:0.12:20];
a=0;
for t=0:ts:T

%     %mozgás útvonala van megadva 0 sebességekkel!!!
%     xs1=10+rc*cos(wc*t);
%     ys1=10+rc*sin(wc*t);

    xs1=5+vx*t;
    ys1=10+vy*t;

% 
    
%     xs1=xs1v(floor(t/ts+1))+1;
    
%     vx=a*t;

    clf
    r1=sqrt((x-xs1-t*vx).^2+(y-ys1-t*vy).^2);
    r1(r1<rLim)=rLim;
    del1=floor(r1/c/ts);
     
%     r2=sqrt((x-xs2-t*vx).^2+(y-ys2-t*vy).^2);
%     r2(r2<rLim)=rLim;
%     del2=floor(r2/c/ts);

    for xTemp= 1:L/ss
        %xTemp
        for yTemp= 1:L/ss
            if del1(xTemp,yTemp)>length(impulse(xTemp,yTemp,:))
                impulse(xTemp,yTemp,del1(xTemp,yTemp))=signal(floor(t/ts+1))/r1(xTemp, yTemp);
            else
                impulse(xTemp,yTemp,del1(xTemp,yTemp))=impulse(xTemp,yTemp,del1(xTemp,yTemp))+signal(floor(t/ts+1))/r1(xTemp, yTemp);
            end
        end
    end
    
%     for xTemp= 1:L/ss
%         %xTemp
%         for yTemp= 1:L/ss
%             if del2(xTemp,yTemp)>length(impulse(xTemp,yTemp,:))
%                 impulse(xTemp,yTemp,del2(xTemp,yTemp))=signal(floor(t/ts+1))/r2(xTemp, yTemp);
%             else
%                 impulse(xTemp,yTemp,del2(xTemp,yTemp))=impulse(xTemp,yTemp,del2(xTemp,yTemp))+signal(floor(t/ts+1))/r2(xTemp, yTemp);
%             end
%         end
%     end

    
    



%     %MINTAVÉTEL
%     recP=[];
%     for m=1:micNum
% %         valTemp=(impulse( 7/ss+1,((m-1)*d+3)/ss+1,1)  +   impulse( 7/ss+1,((m-1)*d+3)/ss+1 ,2)   +   impulse( 7/ss+1,((m-1)*d+3)/ss+1,3)  )/3;
%         valTemp=impulse( 7/ss+1,((m-1)*d+3)/ss+1,1);
%         recP=[recP valTemp];    %mintavétel
%     end    
%     rec=[rec; recP];
%     


% HANGTÉR KIJELZÉSE

%     subplot(micNum,5,picPos);

%     imagesc(impulse(:,:,1));
%     imagesc((impulse(:,:,1)+impulse(:,:,2)+impulse(:,:,3))/3);
    imagesc((impulse(:,:,1)+impulse(:,:,2))/2);
    colormap(TempMap);

    
    title('soundfield t=', t);
    hold on;
    

    
%     %PONTOK KIÍRÁSA
%     %itt valami fura, meg kellett fordítani az x meg az y koordinátákat, különben tükrözve ment azt egész
%     for m=1:micNum
%         plot(((m-1)*d+3)/ss+1, 7/ss+1, 'k.','MarkerSize',10);
%     end
%         
%     plot((xs1+L)/ss+1, (ys1+L)/ss+1,'r.');

    
%     %IDŐFÜGGVÉNYEK KIÍRÁSA
%      for m=1:micNum
%          subplot(micNum,5,[(m-1)*5+4 (m-1)*5+5]); 
%          plot(rec(:,m));
%          axis([0 T/ts -1 1]);
%      end

    %GIF
%       plot(x,y)
      %drawnow
      zlim([-1 1]);
      caxis([-1 1]);
      frame = getframe(1);
      writeVideo(v,frame);
      writeVideo(v,frame);
      
    impulse=impulse(:,:,2:end);
    pause(0.001);

end
close(v);


% for m=1:micNum
%     figure(22)
%     hold on;
%     plot(rec(:,m));
% end

%% DELAY AND SUM
t=0:ts:T;
x=0:1:T/ts;
figure(7);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
o=-1;
delSumGraph=[];


Theta=-pi/2:0.1:pi/2;
for k=1:length(Theta)
    
    clf
    
    %ThetaD=-40;
    %Theta=pi*ThetaD/180;
    %floor(sin(Theta)*2/v/0.05)
    

    delVal=[];  %mikrofonok késleltetési vektora
    for m=1:micNum
        delVal=[delVal -floor(sin(Theta(k))*d*(m-1)/c/ts)];
    end
    
    o=delVal(2);
    
    del1=[]; %késleletetett jelmátrix
    
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

        subplot(5, 1, 1);
        
        for m=1:micNum
            hold on;
            plot(del(m,:));
        end
        
        title('theta= ', Theta/pi*180)
        %plot(A(n):end);
        axis([0 length(rec(:,1)) -1 1])
        
        subplot(5,1,2);
        
        plot(delSum);
        title('B');
        axis([0 length(rec(:,1)) -3 3])
        
        subplot(5,1,3);
        plot(Theta(1:k), delSumGraph);
        axis([-pi/2 pi/2 0 max(delSumGraph)])
        subplot(5,1,[4 5]);
        polarplot(Theta(1:k)+pi/2, sqrt(delSumGraph/ref));
        axis([0 180 0 max(sqrt(delSumGraph/ref))])
%         polarplot(Theta(1:k)+pi/2, 2*sqrt(delSumGraph)/T);
%         axis([0 180 0 max(2*sqrt(delSumGraph)/T)])
% %         polarplot(Theta(1:k)+pi/2, log(0.01+delSumGraph/max(delSumGraph)));
% %         axis([0 180 -2 max(log(0.01+delSumGraph/max(delSumGraph)))])
        
        pause(0.01);
        
       
end
