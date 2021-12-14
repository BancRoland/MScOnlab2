%% antennasor 3 antennával
clear all;

TempMap=[];
for i=0:0.01:1
   TempMap=[TempMap; 0+i 0+i 1];
end
for i=0.01:0.01:1
   TempMap=[TempMap; 1 1-i 1-i];
end

ss=0.1; %términtavételi felbontás [m]

L=10;   %pályamáret
[x,y]=meshgrid(-L:ss:L); %x és y is -3-tól 3-ig megy 1-esével

alfa=60;    %irányszög [deg]
xs=-1*L*sind(alfa);   %hangforrás pozíciója
ys=-1*L*cosd(alfa);

rl=0.5; %távlimit
r=sqrt((x-xs).^2+(y-ys).^2);
r(r<rl)=rl;

c=2;    %hullámsebesség


% jeltulajdonságok
Amp=5;

Amp1=1;
lambda=2;
k=2*pi/lambda;
w=2*pi/(lambda/c);

Amp2=0.5;
lambda2=3.141592;
k2=2*pi/lambda2;
w2=2*pi/(lambda2/c);

Amp3=0.2;
lambda3=1;
k3=2*pi/lambda3;
w3=2*pi/(lambda3/c);

d=1; % Antennasor távolsága [1m] (1 tizedes jegy!)

Scx=L; %mikrofontömb középpontja
Scy=0;

Rx=0;
Gx=0;
Bx=0;
Ry=0;
Gy=d;
By=-d;

%mikrofonok által vett jelsorok
R=[];
G=[];
B=[];

ts=0.01; %mintavételi idő  [sec]
T=5;    %időablak szélessége [sec]


figure(1);
set(gcf,'units','normalized','outerposition',[0 0 1 1])

for t=0:ts:T

    clf
    f=Amp*real(Amp1*exp(1i*(k.*r-t*w))./r + Amp2*exp(1i*(k2.*r-t*w2))./r + Amp3*exp(1i*(k3.*r-t*w3))./r); %jelforrás
    
    R=[R f((Scx+Rx+L)/ss+1, (Scy+Ry+L)/ss+1)];    %mintavétel
    G=[G f((Scx+Gx+L)/ss+1, (Scy+Gy+L)/ss+1)];
    B=[B f((Scx+Bx+L)/ss+1, (Scy+By+L)/ss+1)];
    
    subplot(3,5,[1 2 3 6 7 8 11 12 13]);
    %imagesc(f, title="R");
    imagesc(f);
    colormap(TempMap);
    title('soundfield t=', t);
    hold on;
    %itt valami fura, meg kellett fordítani az x meg az y koordinátákat, különben tükrözve ment azt egész
%     plot((Scy+Ry+L)/ss+1, (Scx+Rx+L)/ss+1, '.','filled','r');
%     plot((Scy+Gy+L)/ss+1, (Scx+Gx+L)/ss+1, '.','filled','g');
%     plot((Scy+By+L)/ss+1, (Scx+Bx+L)/ss+1, '.','filled','b');

    scatter((Scy+Ry+L)/ss+1, (Scx+Rx+L)/ss+1, 'filled', 'ro');
    scatter((Scy+Gy+L)/ss+1, (Scx+Gx+L)/ss+1, 'filled', 'go');
    scatter((Scy+By+L)/ss+1, (Scx+Bx+L)/ss+1, 'filled', 'bo');
    scatter((xs+L)/ss+1, (ys+L)/ss+1, 'filled', 'ko');
    zlim([-1 1]);
    caxis([-1 1]);
    
    subplot(3,5,[4 5]);
    plot(B,'-b','LineWidth',2);
    title('B');
    axis([0 T/ts -1 1]);
    
    subplot(3,5,[9 10]); 
    plot(R,'-r','LineWidth',2);
    title('R');
    axis([0 T/ts -1 1]);
    
    subplot(3,5,[14 15]);  
    plot(G,'-g','LineWidth',2);
    title('G');
    axis([0 T/ts -1 1]);
    
    pause(ts/200); 

end

hold off;
%%
figure(2)
plot(R);
hold on
plot(G);
plot(B);

% %%
% figure(3)
% plot(angle(fft(R)));
% 
% figure(4)
% plot(angle(fft(G)));
% 
% figure(5)
% plot(angle(fft(B)));
% 
% figure(6)
% plot(abs(fft(G)));

%% delay
t=0:ts:T;
x=0:1:T/ts;
figure(7);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
o=0;
delSumGraph=[];

delSumGraph=[];
Theta=-pi/2:0.1:pi/2;
for k=1:length(Theta)
    
    clf
    %ThetaD=-40;
    %Theta=pi*ThetaD/180;
    %floor(sin(Theta)*2/v/0.05)
    
    nR=(floor(sin(Theta(k))*0/c/ts));
    nG=(floor(sin(Theta(k))*-d/c/ts));
    nB=(floor(sin(Theta(k))*d/c/ts));
    
    o=1;
    if nR~=o
        delR=[zeros(1, nR+length(R)) R zeros(1, -nR+length(R))];
        delG=[zeros(1, nG+length(G)) G zeros(1, -nG+length(G))];
        delB=[zeros(1, nB+length(B)) B zeros(1, -nB+length(B))];

        delSum=delR + delG + delB;

    end
    
    delSumGraph=[delSumGraph delSum*delSum'];
    o=nR;



        subplot(5, 1, 1);
        hold on;
        title('Késleltetett jelek');
        plot(delR, 'r');
        plot(delG, 'g');
        plot(delB, 'b');
%         plot(x-nR,R);
%         hold on;
%         plot(x-nG,G);
%         plot(x-nB,B);
        title('theta= ', Theta/pi*180)
        %plot(A(n):end);
        axis([0 3*length(R) -1 1])
        
        subplot(5,1,2);
        hold on;
        title('Összegzett késleltetett jelek');
        delSum=delR + delG + delB;
        plot(delSum);
        title('B');
        axis([0 3*length(R) -3 3])
        
        subplot(5,1,3);
        hold on;
        title('Iránykarakterisztika');
        plot(Theta(1:k), delSumGraph);
        axis([-pi/2 pi/2 0 15])
        subplot(5,1,[4 5]);
        polarplot(Theta(1:k)+pi/2, delSumGraph);
        hold on;
        title('Iránykarakterisztika');
        axis([0 180 0 max(delSumGraph)])
        
        pause(0.1);
        
        

end




% for Theta=-pi/2:0.1:pi/2
%     
%     clf
%     %ThetaD=-40;
%     %Theta=pi*ThetaD/180;
%     %floor(sin(Theta)*2/v/0.05)
%     
%     nR=(floor(sin(Theta)*0/c/ts));
%     nG=(floor(sin(Theta)*-d/c/ts));
%     nB=(floor(sin(Theta)*d/c/ts));
%     
%     o=1;
%     if nR~=o
%         subplot(3, 1, 1);
% 
%         delR=[zeros(1, nR+length(R)) R zeros(1, -nR+length(R))];
%         delG=[zeros(1, nG+length(G)) G zeros(1, -nG+length(G))];
%         delB=[zeros(1, nB+length(B)) B zeros(1, -nB+length(B))];
%         plot(delR);
%         hold on;
%         plot(delG);
%         plot(delB);
% %         plot(x-nR,R);
% %         hold on;
% %         plot(x-nG,G);
% %         plot(x-nB,B);
%         title('theta= ', Theta/pi*180)
%         %plot(A(n):end);
%         axis([0 3*length(R) -1 1])
%         
%         
%         subplot(3,1,2);
%         delSum=delR + delG + delB;
%         plot(delSum);
%         title('B');
%         axis([0 3*length(R) -3 3])
%         
%         subplot(3,1,3);
%         
%         plot(delSumGraph);
%         title('B');
%         %axis([0 pi/2 -1 2]);
%         
%         pause(0.1);
%     end
%     
%     delSumGraph=[delSumGraph delSum*delSum'];
%     o=nR;
%     figure(7);
%     plot(delSumGraph);
% end