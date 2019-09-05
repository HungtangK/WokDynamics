clc
close all
% Two link model

global s1 s2 c1 c2 L d D t;
global X Y VX VY AX AY deltat PhiIndex;


deltat=0.01;
t=(0:deltat:4.25)';
xstove=930;
ystove=675;
ratio=0.765;
omega=3.047;    
d=130;
D=390;
L=1.625*d;
PHI0=1.176;


%  s1(t)
phis1=deg2rad(0);
As1=0.20;
s10=0.24820843;
s1=s10+As1*cos(2*pi*omega*t+phis1);

% c1(t)
c1=(1-s1.^2).^0.5;


%  s2(t)
phis2=deg2rad(8.178000808-(-48.16109633));
% phis2=phis1;
As2=0.20;
s20=0.031723578;
s2=s20+As2*cos(2*pi*omega*t+phis2);

% c2(t)
c2=(1-s2.^2).^0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
plot(t,s1,t,s2);
legend("sin(\theta_1)","sin(\theta_2)")


figure 
hold all
axis equal
grid on
xlabel("x/mm")
ylabel("y/mm")
title("Trajectories at different wok locations")
for i=-1:0.25:1
    [xtemp,ytemp]=rigid(PHI0*i,t);
    plot(xtemp,ytemp)
end

load expdata
plot(xL_exp,yL_exp,'.')
plot(xR_exp,yR_exp,'.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% horizontal velocity

[xtemp,ytemp]=rigid(0,t);
v=zeros(size(xtemp));
v(1:end-1)=xtemp(2:end)-xtemp(1:end-1);
v=v/deltat/1000;
disp("The largest horizontal velocity is: "+num2str(min(v))+" m/s")

figure(3)
plot(t(1:end-1),v(1:end-1))
title("horizontal velocity at the bottom(m/s)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vertical acceleration

[xtemp1,ytemp]=rigid(-PHI0,t);
a=zeros(size(ytemp));
a(2:end-1)=ytemp(1:end-2)-2*ytemp(2:end-1)+ytemp(3:end);
a=a/deltat^2/1000;
disp("The largest vertical acceleration is: "+num2str(min(a))+" m/s^2")

figure(4)
plot(t(2:end-1),a(2:end-1))
title("vertical accleration on the left end(m/s^2)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trajectory

i_to=find(a==min(a));

PhiIndex=PHI0*[-1:0.1:1];

[X,Y]=rigid(PhiIndex,t);

VX=zeros(size(X));
VY=zeros(size(X));
AX=zeros(size(X));
AY=zeros(size(X));

VX(2:end-1,:)=(X(3:end,:)-X(1:end-2,:))/2/deltat;
VY(2:end-1,:)=(Y(3:end,:)-Y(1:end-2,:))/2/deltat;
AX(2:end-1,:)=(X(3:end,:)-2*X(2:end-1,:)+X(1:end-2,:))/deltat^2;
AY(2:end-1,:)=(Y(3:end,:)-2*Y(2:end-1,:)+Y(1:end-2,:))/deltat^2;

XP=zeros(size(X));
YP=zeros(size(X));
i_takeoff=zeros(size(PhiIndex));
i_land=zeros(size(PhiIndex));

[XP,YP,i_takeoff,i_land]=projectile();

figure(3)
hold all
plot(t,VX(:,(end+1)/2)/1000)
figure(4)
hold all
plot(t,AY(:,1)/1000)


figure(5)
% hold all
axis equal  
clear F
grid on

for i=1:0.33/deltat*3
    
    plot([X(i,:),X(i,1)],[Y(i,:),Y(i,1)]);
    hold on
    quiver(X(i,:),Y(i,:),VX(i,:),VY(i,:));
    for j=1:length(PhiIndex)
        if AY(i,j)<-9800
            plot(X(i,j),Y(i,j),'ok')
        end
        if i>=i_takeoff(j) && i<=i_land(j)
            plot(XP(i_takeoff(j):i,j),YP(i_takeoff(j):i,j),'-o')
        end
    end
    grid on
    axis equal
    text(100,-50,['t=',num2str(t(i))])
    hold off
    xlim([-300,200])
    ylim([-150,350])
    F(i)=getframe;
end


vidObj = VideoWriter('movie.avi');
vidObj.FrameRate=1/deltat;
vidObj.FrameRate=10;
open(vidObj);
writeVideo(vidObj,F);
% movie(F,20)
close(vidObj);
  




function [xt,yt]=rigid(PHI,t_temp)
    global s1 s2 c1 L t;
    index=find(t==t_temp);
    xt=-L*(s1-sin(asin(s2)+PHI));
    yt=L*(c1-cos(asin(s2)+PHI));
    xt=xt(index,:);
    yt=yt(index,:);
end

function [xt,yt]=projectile_single(j_PHI,i_to)
    global X Y VX VY AY deltat L t;
    global s1 c1;
    if AY(i_to,j_PHI)<-9800
        xt=X(i_to,j_PHI);
        yt=Y(i_to,j_PHI);
        Vxt=VX(i_to,j_PHI);
        Vyt=VY(i_to,j_PHI);
        tt=t(i_to);
        
        xt=[xt;xt(end)+Vxt*deltat];
        yt=[yt;yt(end)+Vyt*deltat];
        Vyt=Vyt-9800*deltat;
        tt=tt+deltat;
        i_time=i_to;
        
        while norm([xt(end)-(-L*s1(i_time+1)),yt(end)-(L*c1(i_time+1))])<1*L || Vyt>0
            xt=[xt;xt(end)+Vxt*deltat];
            yt=[yt;yt(end)+Vyt*deltat];
            Vyt=Vyt-9800*deltat;
            tt=tt+deltat;
            i_time=i_time+1;
            
        end
    end
end

function [XP,YP,i_takeoff,i_land]=projectile()
    global PhiIndex t AY;
    flag=zeros(size(PhiIndex));
    for i=1:length(t)
        for j=1:length(PhiIndex)
            if AY(i,j)<-10000 && flag(j)==0
                flag(j)=1;
                i_takeoff(j)=i;
                [xt,yt]=projectile_single(j,i);
                i_land(j)=i_takeoff(j)-1+length(xt);
                XP(i_takeoff(j):i_land(j),j)=xt;
                YP(i_takeoff(j):i_land(j),j)=yt;
            end
        end
    end
end
