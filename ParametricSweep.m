clc
close all
clear
% Two link model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter definitions
global s1 s2 c1 c2 L d D t;
global X Y VX VY AX AY deltat PhiIndex PHI0;
% Stove rim location 
% xstove=930;
% ystove=675;
% ratio=0.765;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Parameter Changes    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega=-1;
A=-1;
phis=-1;
% load datamaster.mat
% Parameters
datamaster=[];
% for phis=deg2rad(0:5:360)
% for phis=deg2rad(0:5:180)
% for phis=[0,deg2rad(8.178000808-(-48.16109633)),1.5,2,pi,4,4.625,5,6]
% for omega=1:0.1:5
% for omega=1:0.25:5
% for A=0.1:0.1:1
for q=0
for p=0
% for phis=4.625
    q=1;
    p=1;
    
    close all
    
    flag_movie=1;
    flag_load=0;
    flag_figs=0;
    flag_longoutput=0;
    flag_moreSite=0;
    total_cycles=3;
    
    if omega==-1
        omega=3.047; 
    end
    if phis==-1
        phis=deg2rad(8.178000808-(-48.16109633));
    end
    if A==-1
        A=0.2;
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

As1=A;
As2=A;

% Time 
deltat=0.002;
t=(0:deltat:total_cycles/omega)';

% Wok dimension
d=130;
D=390;
L=1.625*d;
PHI0=1.176;

% s1(t)
s10=0.24820843;
s1=s10+As1*cos(2*pi*omega*t);
s1=sin(s1);

% c1(t)
c1=(1-s1.^2).^0.5;
% s2(t)
s20=0.031723578;
s2=s20+As2*cos(2*pi*omega*t+phis);
s2=sin(s2);
% c2(t)
c2=(1-s2.^2).^0.5;


% Fig 1 Trajectory at different locations
if flag_figs
    figure (100) 
    hold all
    axis equal
    xlabel("x/mm")
    ylabel("y/mm")
%     title("Trajectories at different wok locations")
    for i=-1:0.25:1
        [xtemp,ytemp]=rigid(PHI0*i,t);
        plot(xtemp,ytemp,'color','[0.31,0.31,0.31]','LineWidth',4)
    end
    [xtemp,ytemp]=rigid(PHI0*[-1:0.05:1],0);
    plot(xtemp,ytemp,'k','LineWidth',8)

    if flag_load
        load expdata
        xL_exp=0.765*(xL_exp-930);
        yL_exp=-0.765*(yL_exp-675);
        xR_exp=0.765*(xR_exp-930);
        yR_exp=-0.765*(yR_exp-675);
        xC_exp=0.765*(xC_exp-930);
        yC_exp=-0.765*(yC_exp-675);
        plot(xL_exp,yL_exp,'.','color',[0.6,0.2,0]')
        plot(xR_exp,yR_exp,'.','color',[0.6,0.2,0]')
        plot(xC_exp,yC_exp,'g.')
    end
    axis off
    saveas(gcf,"WokTraj_Omega="+num2str(omega,'%.2f')+"_Phid="+num2str(phis,'%.2f')+"_A="+num2str(A,'%.2f')+".fig")
    saveas(gcf,"WokTraj_Omega="+num2str(omega,'%.2f')+"_Phid="+num2str(phis,'%.2f')+"_A="+num2str(A,'%.2f')+".svg")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jumping site kinematics analysis

% The Phi we are investigating as jumping sites
if flag_moreSite
    PhiIndex=PHI0*[-1:0.01:1];
else
    PhiIndex=PHI0*[-1:0.1:1];
end
% The trajectoryies of jumping sites
% First index of X and so on is time and second location
[X,Y]=rigid(PhiIndex,t);

% The velocity of the jumping sites
VX=zeros(size(X));
VY=zeros(size(X));
VX(2:end-1,:)=(X(3:end,:)-X(1:end-2,:))/2/deltat;
VY(2:end-1,:)=(Y(3:end,:)-Y(1:end-2,:))/2/deltat;
VX(1,:)=(-X(3,:)+4*X(2,:)-3*X(1,:))/2/deltat;
VY(1,:)=(-Y(3,:)+4*Y(2,:)-3*Y(1,:))/2/deltat;
VX(end,:)=(3*X(end,:)-4*X(end-1,:)+X(end-2,:))/2/deltat;
VY(end,:)=(3*Y(end,:)-4*Y(end-1,:)+Y(end-2,:))/2/deltat;

% The acceleration of the jumping sites
AX=zeros(size(X));
AY=zeros(size(X));
AX(2:end-1,:)=(X(3:end,:)-2*X(2:end-1,:)+X(1:end-2,:))/deltat^2;
AY(2:end-1,:)=(Y(3:end,:)-2*Y(2:end-1,:)+Y(1:end-2,:))/deltat^2;
AX(1,:)=(-VX(3,:)+4*VX(2,:)-3*VX(1,:))/2/deltat;
AY(1,:)=(-VY(3,:)+4*VY(2,:)-3*VY(1,:))/2/deltat;
AX(end,:)=(3*VX(end,:)-4*VX(end-1,:)+VX(end-2,:))/2/deltat;
AY(end,:)=(3*VY(end,:)-4*VY(end-1,:)+VY(end-2,:))/2/deltat;

% Calculate all the projectile trajectories
% First index is time and second location
global XP YP i_takeoff i_land LP PhiP;
XP=zeros(size(X));
YP=zeros(size(X));
i_takeoff=zeros(size(PhiIndex));
i_land=zeros(size(PhiIndex));
LP=zeros(size(X));
PhiP=zeros(size(X));
    
% calculate XP,YP,i_takeoff,i_land
projectile();
% calculate LP,PhiP
WokCoordinate(XP,YP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Work Analaysis

% Mass/kg
M=1.7;

% Center of mass
Dc=L/4*((1-cos(2*PHI0))/(1-cos(PHI0)));

% Moment of Inertia
Ic=M*(L^2/6*(4-3*cos(PHI0)-cos(PHI0)^3)/(1-cos(PHI0))/1000^2-Dc^2);
% Ic=0;

% Translational accleration Angular accleration of the center of mass
% axC in m/s^2 and alpha in rad/s^2
[axC,ayC,xC,yC,alphaC]=CenterKinematics(Dc);

% Inertial force
FxC=zeros(size(t));
FyC=zeros(size(t));
MC=zeros(size(t));
FxC=-M*axC;
FyC=-M*ayC;
MC=-Ic*alphaC;

% Gravity
FG=-M*9.8;

% Handle position
[xH,yH]=rigid(PHI0+deg2rad(5),t);

% handle velocity
VxH=zeros(size(xH));
VyH=VxH;
VxH(2:end-1)=(xH(3:end)-xH(1:end-2))/2/deltat;
VyH(2:end-1)=(yH(3:end)-yH(1:end-2))/2/deltat;

% moment arm
armN=norm([xC,yC])*sin(atan2(-xC,yC)-asin(s1));
armH1=norm([xC-xH,yC-yH])*sin(asin(s1)-atan2(yH-yC,xH-xC));
armH2=norm([xC-xH,yC-yH])*cos(asin(s1)-atan2(yH-yC,xH-xC));

% Calculate force
FH1=-(FxC.*c1+(FyC+FG).*s1);
N=-(-FxC.*s1+(FyC+FG).*c1);
M=-(FH1.*armH1+N.*armN+MC);

FHx=FH1.*c1;
FHy=FH1.*s1;

save('handle.mat','t','xH','yH','FHx','FHy')

% Plot Forces
if flag_figs
    PlotForces(FxC,FyC,FG,FH1,M,N,armH1,armN,MC)
    saveas(gcf,"Forces="+num2str(omega,'%.2f')+"_Phid="+num2str(phis,'%.2f')+"_A="+num2str(A,'%.2f')+".fig")
    saveas(gcf,"Forces="+num2str(omega,'%.2f')+"_Phid="+num2str(phis,'%.2f')+"_A="+num2str(A,'%.2f')+".svg")
end
%     figure(8)
%     hold all
%     plot(t,FH1)
%     plot(t,(FH1.^2+N.^2).^0.5)
%     legend("Force with support","Force without support")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output



% figure(3)
% plot(t,VX(:,(end+1)/2)/1000)
% title("horizontal velocity at the bottom(m/s)")
% 
% figure(4)
% plot(t,AY(:,1)/1000)
% title("vertical accleration on the left end(m/s^2)")

if flag_figs
    figure(2)
    hold all
    for j=1:length(PhiIndex)
        if i_takeoff(j)
            plot(PhiP(i_takeoff(j):i_land(j),j),LP(i_takeoff(j):i_land(j),j));
        end
    end
    xlabel("\Phi_P")
    ylabel("L_P")
    % ylim([0,0.5])
end

i_tcoarse=1:20:size(t);
% saveas(gcf,"ProjectileTrace_Omega="+num2str(omega,'%.2f')+"_Phid="+num2str(phis,'%.2f')+"_A="+num2str(A,'%.2f')+".jpeg")

if flag_figs
    figure(1)
    hold all
    ittemp=40;
    for j=1:length(PhiIndex)
        if i_takeoff(j)
            xtemp=sin(PHI0*PhiP(i_takeoff(j):i_land(j),j)).*(1-LP(i_takeoff(j):i_land(j),j))*L;
            ytemp=-cos(PHI0*PhiP(i_takeoff(j):i_land(j),j)).*(1-LP(i_takeoff(j):i_land(j),j))*L;
            xtemp=xtemp-L*s1(ittemp);
            ytemp=ytemp+L*c1(ittemp);
            plot(xtemp,ytemp,'--b','LineWidth',4);
    %         index=intersect(i_tcoarse,i_takeoff(j):i_land(j));
    %         xtemp=sin(PHI0*PhiP(index,j)).*(1-LP(index,j))*L;
    %         ytemp=-cos(PHI0*PhiP(index,j)).*(1-LP(index,j))*L;
    %         plot(xtemp,ytemp,'o');
         end
    end
    xtemp=sin(PHI0*[-1:0.05:1])*L;
    ytemp=-cos(PHI0*[-1:0.05:1])*L;
    xtemp=xtemp-L*s1(ittemp);
    ytemp=ytemp+L*c1(ittemp);
    h=plot(xtemp,ytemp,'k');
    h.LineWidth=8;
    axis equal
    axis off
    saveas(gcf,"RiceTraj_Omega="+num2str(omega,'%.2f')+"_Phid="+num2str(phis,'%.2f')+"_A="+num2str(A,'%.2f')+".fig")
    saveas(gcf,"RiceTraj_Omega="+num2str(omega,'%.2f')+"_Phid="+num2str(phis,'%.2f')+"_A="+num2str(A,'%.2f')+".svg")
end

% display and data saving
disp("**********************************************************************")
disp("omega="+num2str(omega)+" Phis="+num2str(phis)+" A="+num2str(A))
disp("The largest vertical acceleration is: "+num2str(min(AY(:,1)/1000))+" m/s^2")
disp("The largest horizontal velocity is: "+num2str(min(VX(:,(end+1)/2)/1000))+" m/s")
if i_land(1)
disp("The time it takes for rice grain on the most left fall is: "+num2str((i_land(1)-i_takeoff(1))*deltat)+" s")
disp("The farthest landing location of the rice:"+num2str(PhiP(i_land(1)))+" rad")
disp("The range where the location satisfy the takeoff condition is: "+num2str(sum(i_takeoff(1:(end+1)/2)~=0))+" no unit")
end
disp("The largest force needed at the handle : "+num2str(max(FH1))+" N")
disp("The average force needed at the handle : "+num2str(mean(FH1(3:3+floor(1/omega/deltat))))+" N")

if i_land(1)
    tosave=[
        omega;
        phis;
        A;
        min(AY(:,1)/1000);
        min(VX(:,(end+1)/2)/1000);
        PhiP(i_land(1));
        (i_land(1)-i_takeoff(1))*deltat;
        sum(i_takeoff(1:(end+1)/2)~=0);
        max(FH1);
        mean(FH1(3:3+floor(1/omega/deltat)))
        max(LP(i_takeoff(1):i_land(1),1))];
else 
    tosave=[
        omega;
        phis;
        A;
        min(AY(:,1)/1000);
        min(VX(:,(end+1)/2)/1000);
        0;
        (i_land(1)-i_takeoff(1))*deltat;
        sum(i_takeoff(1:(end+1)/2)~=0);
        max(FH1);
        mean(FH1(3:3+floor(1/omega/deltat)))
        0];
end
if flag_longoutput
   tosave=[
        phis;
        min(AY/1000,[],1)'
        ];
end

datamaster=[datamaster,tosave];

end
end

save("datamaster.mat","datamaster")

%%
T=1/omega/deltat;
playSpeed=0.2;
% Movie
if flag_movie
    f=figure;
    f.Visible='off';
    F=[];
    ax1 = axes('Position',[0.1 0.1 0.8 0.8]);
%     ax2 = axes('Position',[0.55 0.72 0.3 0.15]);
    ax2 = axes('Position',[0.4 0.72 0.43 0.15]);
    plot(t,s1,t,s2)
%     text(t(end)-0.2,0.35,'\theta_1','Color','b')
%     text(t(end)-0.1,-0.08,'\theta_2','Color','r')
    text(t(end)+0.05,0.30,'\theta_1','Color','b')
    text(t(end)+0.05,-0.00,'\theta_2','Color','r')
    hold on

    
    h2=plot(ax2,[0,0],[-100,100]);
    ax2.YLim=[-0.2,0.5];
    ax2.XLim=[0,t(end)];

    hold(ax1,'on')
    ax1.DataAspectRatio=[1 1 1];
    xlim(ax1,[-300,200])
%     ylim(ax1,[-100,250])
    ylim(ax1,[-50,350])
    
    i=1;
    plot(ax1,0,-10,'^','MarkerSize',10,'MarkerFaceColor',[1,140/255,0],'MarkerEdgeColor',[1,140/255,0])
%     t1=text(ax1,-250,-60,['t=',num2str(t(i))]);
%     t2=text(ax1,-250,-30,[num2str(playSpeed),'x']);
    t1=text(ax1,-250,310,['t=',num2str(t(i)),' s']);
    t2=text(ax1,-250,280,[num2str(playSpeed),'x']);
    h1wok=plot(ax1,[X(i,:),X(i,1)],[Y(i,:),Y(i,1)],'k','LineWidth',2);
    h1lever=plot(ax1,[0,-L*s1(i),-L*(s1(i)-s2(i))],[0,L*c1(i),L*(c1(i)-c2(i))],'.-','MarkerSize',30,'LineWidth',5,'Color',[1,140/255,0]);
    for j=1:length(PhiIndex)
        h1_p(j)=plot(ax1,0,0,'--b','LineWidth',1);
        h1_p(j).Visible='off';
    end
        
    for i=1:length(t)
        if mod(i,50)==1
            disp(['t=',num2str(t(i))])
        end
        t1.String=['t=',num2str(t(i)),' s'];
        h1wok.XData=[X(i,:),X(i,1)];
        h1wok.YData=[Y(i,:),Y(i,1)];
        h1lever.XData=[0,-L*s1(i),-L*(s1(i)-s2(i))];
        h1lever.YData=[0,L*c1(i),L*(c1(i)-c2(i))];
        
        for j=1:length(PhiIndex)
            if mod(i,T)>=i_takeoff(j) && mod(i,T)<=i_land(j)
                h1_p(j).Visible='on';
                Base=floor(i-mod(i,T));
                h1_p(j).XData=XP(i_takeoff(j):i-Base,j);
                h1_p(j).YData=YP(i_takeoff(j):i-Base,j);
%                 plot(ax1,XP(i_takeoff(j):i-Base,j),YP(i_takeoff(j):i-Base,j),'--b','LineWidth',1);
            else
                h1_p(j).Visible='off';
            end
        end
        h2.XData=[t(i),t(i)];

%         F(i)=getframe(gcf);
        if t(i)<3 || t(i)>5
            F=[F,getframe(gcf)];
            t2.String=[num2str(playSpeed),'x'];
        else
            F=[F,getframe(gcf),getframe(gcf),getframe(gcf),getframe(gcf),getframe(gcf)];
            t2.String=[num2str(playSpeed/5),'x'];
        end
    end

 
    % video output

    vidObj = VideoWriter('movie');
    vidObj.FrameRate=playSpeed/deltat;
    open(vidObj);
    writeVideo(vidObj,F);
    close(vidObj);
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%      Functions      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Rigid takes in the angular position of the wok and the time in sec and
% return the xy coordinate of the wok at the time instant
function [xt,yt]=rigid(PHI,t_temp)
    global s1 s2 c1 L t;
    index=find(t==t_temp);
    xt=-L*(s1-sin(asin(s2)+PHI));
    yt=L*(c1-cos(asin(s2)+PHI));
    xt=xt(index,:);
    yt=yt(index,:);
end

% Projectile_single takes in the angular position INDEX and the time
% INDEX and returns the trajectory of the projectile
function [xt,yt]=projectile_single(j_PHI,i_to)
    global X Y VX VY deltat L t;
    global s1 c1;

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

    while norm([xt(end)-(-L*s1(i_time+1)),yt(end)-(L*c1(i_time+1))])<L
        xt=[xt;xt(end)+Vxt*deltat];
        yt=[yt;yt(end)+Vyt*deltat];
        Vyt=Vyt-9800*deltat;
        tt=tt+deltat;
        i_time=i_time+1;     
    end
    
end
%%
% Projectile returns the trajectories of all the projectiles as well as
% their takeoff and landing time INDEX
function projectile()
    global PhiIndex t AY;
    global XP YP i_takeoff i_land;
    global AX s2 deltat;
    flag=zeros(size(PhiIndex));
    for i=1:length(t)
        for j=1:(length(PhiIndex)+1)/2
          if  -AX(i,j)*sin(-PhiIndex(j)-asin(s2(i)))+(-9800-AY(i,j))*cos(-PhiIndex(j)-asin(s2(i)))>0 && flag(j)==0
                flag(j)=1;
                i_takeoff(j)=i;
                [xt,yt]=projectile_single(j,i);
                i_land(j)=i_takeoff(j)-1+length(xt);
                XP(i_takeoff(j):i_land(j),j)=xt;
                YP(i_takeoff(j):i_land(j),j)=yt;
            end
        end
    end
    for i=0.1/deltat:length(t)
        for j=(length(PhiIndex)+1)/2:length(PhiIndex)
          if  -AX(i,j)*sin(-PhiIndex(j)-asin(s2(i)))+(-9800-AY(i,j))*cos(-PhiIndex(j)-asin(s2(i)))>0 && flag(j)==0
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
%%
% WokCoordinate returns the XY cordinate in the stationary wok coordinate
% system
function WokCoordinate(XP,YP)
    global t PhiIndex L s1 c1 s2;
    global LP PhiP PHI0;

    for i=1:length(t)
        xCenter=-s1(i)*L;
        yCenter=c1(i)*L;
        for j=1:length(PhiIndex)
            if XP(i,j)~=0
                LP(i,j)=norm([XP(i,j),YP(i,j)]-[xCenter,yCenter]);
                LP(i,j)=1-LP(i,j)/L;
                PhiP(i,j)=atan2(xCenter-XP(i,j),yCenter-YP(i,j))+asin(s2(i));
                PhiP(i,j)=-PhiP(i,j);
                PhiP(i,j)=PhiP(i,j)/PHI0;
            end           
        end
    end 
end
%%
function [axC,ayC,xC,yC,alphaC]=CenterKinematics(Dc)
    global s1 s2 c1 c2 L deltat t;
    xC=-L*s1+Dc*s2;
    yC=L*c1-Dc*c2;
    axC=zeros(size(xC));
    ayC=zeros(size(xC));
    alphaC=zeros(size(xC));
    axC(2:end-1)=(xC(1:end-2)-2*xC(2:end-1)+xC(3:end))/deltat^2/1000;
    ayC(2:end-1)=(yC(1:end-2)-2*yC(2:end-1)+yC(3:end))/deltat^2/1000;
    alphaC(2:end-1)=(asin(s2(1:end-2))-2*asin(s2(2:end-1))+asin(s2(3:end)))/deltat^2;
%     figure(300)
%     hold all
%     plot(t,axC)
%     plot(t,2*L/1000*2*(pi*3.05)^2*0.2*(cos(2*pi*3.05*t)-0.75*cos(2*pi*3.05*t+1.176)))

%     figure(400)
%     hold all
%     plot(t,ayC)
%     ayyC=(sin(2*pi*3.05*t).^2-0.7*sin(2*pi*3.05*t+1.176).^2);
%     ayyC=ayyC+s1.^2-0.7*s2.^2;
%     ayyC=ayyC*2*L/1000*2*(pi*3.05)^2*0.2^2;
% %     plot(t,ayyC)
%     plot(t,yC)
%     ayyC= -c1.*s1.^2+s1.*c1;
%     ayyC=ayyC+0.7* (-c2.*s2.^2+s2.*c2);
%     ayyC=ayyC*2*L/1000*2*(pi*3.05)^2*0.2;
%     plot(t,ayyC)
end
%%
function PlotForces(FxC,FyC,FG,FH1,M,N,armH1,armN,MC)
    global t c1 s1 L;
    figure(9)
    subplot(2,1,1)
    hold all
    plot(t,FH1.*c1,'k','LineWidth',2);
    plot(t,FxC,'r','LineWidth',2);
    plot(t,N.*-s1,'b','LineWidth',2)
%     plot(t,FH1.*c1+FxC+N.*-s1);sum
    A=0.2;
    omega=3.047;
    Dc=0.6923;
    phis=0.9833;
%     plot(t,1.7*L/1000*(2*pi*omega)^2*A*(cos(2*pi*omega*t)-Dc*cos(2*pi*omega*t+phis)),'k--','LineWidth',1)

    legend('Handle','Inertia','Support')
    ylabel("x direction")
    set(gca,'fontsize', 18)
    set(gca,'Position', [0.13,0.55,0.775,0.35])
    set(gca,'XTick', [])
    xlim([0,0.6])
    ylim([-35,35])
    
    subplot(2,1,2)
    hold all
    plot(t,FH1.*s1,'k','LineWidth',2);
    plot(t,FyC,'r','LineWidth',2);
    plot(t,N.*c1,'b','LineWidth',2);
    plot(t,FG.*ones(size(t)),'g','LineWidth',2);
    
%     plot(t,FH1.*s1+FyC+FG+N.*c1)sum
    legend('Handle','Inertia','Support','Gravity')
    ylabel("y direction")
    xlabel("t/s")
    set(gca,'fontsize', 18)
    set(gca,'Position', [0.13,0.175,0.775,0.35])
    ylim([-35,35])
    xlim([0,0.6])

%     subplot(3,1,3)
%     hold all
%     plot(t,FH1.*armH1)
%     plot(t,M)
%     plot(t,N.*armN)
%     plot(t,MC)
%     plot(t,FH1.*armH1+M+N.*armN-MC)
end 