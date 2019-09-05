clear
load("datamaster_omega.mat")
close all



% hold all
% plot(datamaster(1,:),datamaster(8,:))
% plot(datamaster(1,:),datamaster(9,:))
% 
% plot([3.047,3.047],[0,150],'k--')
% saveas(gcf,"force_omega.fig")
% saveas(gcf,"JPEGs\force_omega.jpg")
% xlabel("Frequency/Hz")
% ylabel("Handle force/N")
% legend("Average","Maximum")


% figure(1)
% 
% hold all
% PHI0=1.176;
% % hr=plot(datamaster(1,:),datamaster(4,:)/-10);
% % hr=plot(datamaster(1,:),datamaster(5,:));
% subplot(2,1,1)
% hold all
% datamaster(6,1:4)=-PHI0;
% plot(datamaster(1,:),datamaster(6,:)/PHI0,'LineWidth',2);
% plot([3.047,3.047],[-1,1.50],'k--','LineWidth',2)
% ylabel("\zeta/\zeta_{max}")
% ylim([-1,0.5])
% set(gca,'fontsize', 18)
% 
% % subplot(,1,2)
% % hold all
% % plot(datamaster(1,:),datamaster(7,:),'LineWidth',2);
% % plot([3.047,3.047],[-1,1.50],'k--','LineWidth',2)
% % plot([1,5],[1,1],'k:','LineWidth',2)
% % ylim([0,0.3])
% % % xlabel("Frequency/Hz")
% % ylabel("t_{travel}")
% % % h=legend("Landing location","Landing time");
% % % h.Location="southeast";
% % set(gca,'fontsize', 18)
% 
% subplot(2,1,2)
% hold all
% plot(datamaster(1,:),datamaster(7,:).*datamaster(1,:),'LineWidth',2);
% plot([3.047,3.047],[-1,1.50],'k--','LineWidth',2)
% plot([1,5],[1,1],'k:','LineWidth',2)
% ylim([0,1.5])
% xlabel("Frequency/Hz")
% ylabel("t_{travel}/T_{cycle}")
% % h=legend("Landing location","Landing time");
% % h.Location="southeast";
% set(gca,'fontsize', 18)
% saveas(gcf,"omega.fig")
% saveas(gcf,"JPEGs\omega.jpg")
% saveas(gcf,"SVGs\omega.svg")
% 
% load("datamaster_phis.mat")
% % figure(2)
% % hold all
% % plot(datamaster(2,:),-datamaster(4,:)/9.8)
% % plot(datamaster(2,:),-datamaster(5,:))
% % plot([deg2rad(8.178000808-(-48.16109633)),deg2rad(8.178000808-(-48.16109633))],[0,2],'k--')
% % xlabel("Phase difference/rad")
% % ylabel("g or m/s^2")
% % legend("Maximum Vertical Accleration at the left","Maximum horizontatl velocity at the bottom")
% % set(gca,'fontsize', 18)
% % saveas(gcf,"velocity_phasediff.fig")
% % saveas(gcf,"JPEGs\velocity_phasediff.jpg")
% % saveas(gcf,"SVGs\velocity_phasediff.svg")
% 
% figure(2)
% hold all
% % index=1:
% hr=plot(datamaster(6,:),datamaster(2,:),'k','LineWidth',2);
% plot([-1,5],[deg2rad(8.178000808-(-48.16109633)),deg2rad(8.178000808-(-48.16109633))],'k--','LineWidth',2)
% plot([-1,5],[4.625,4.625],'--','LineWidth',2,'Color',[0.75,0.75,0.75])
% % text(datamaster(6,end),datamaster(2,end),"\omega_{r}")
% ylabel("Phase Difference/rad")
% xlabel("\zeta/\zeta_{max}")
% 
% figure(3)
% hold all
% hr=plot(datamaster(11,:),datamaster(2,:),'k','LineWidth',2);
% plot([-1,5],[deg2rad(8.178000808-(-48.16109633)),deg2rad(8.178000808-(-48.16109633))],'k--','LineWidth',2)
% plot([-1,5],[4.625,4.625],'--','LineWidth',2,'Color',[0.75,0.75,0.75])
% ylabel("Phase Difference/rad")
% xlabel("h_{max}/L")
% 
% % figure(4)
% % hold all
% % hr=plot(datamaster(7,:).*3.047,datamaster(2,:));
% % plot([0,1.3],[deg2rad(8.178000808-(-48.16109633)),deg2rad(8.178000808-(-48.16109633))],'k--')
% % ylabel("Phase Difference/rad")
% % xlabel("Normalized falling time")
% % set(hr,'LineWidth',2)
% 
% figure(2)
% load("datamaster_omegaphis.mat")
% hold all
% len_phis=length(1:0.25:5);
% len_omega=length(0:5:360);
% datamaster(6,:)=datamaster(6,:)-(datamaster(6,:)==0);
% for j=1:1:len_phis
%     plot(datamaster(6,j:len_phis:end),datamaster(2,j:len_phis:end),'LineWidth',1)
% end
% for j=3:2:len_phis
%     text(datamaster(6,j+len_phis*(len_omega-1))-0.02,datamaster(2,j+len_phis*(len_omega-1))+0.1,num2str(datamaster(1,j)),'FontSize',12)
% end
% j=1;
% text(datamaster(6,j+len_phis*(len_omega-1))-0.02,datamaster(2,j+len_phis*(len_omega-1))+0.1,['f=',num2str(datamaster(1,j))],'FontSize',12)
% 
% xlim([-1,0.5])
% set(gca,'fontsize', 18)
% saveas(gcf,"landloc_phasediff.fig")
% saveas(gcf,"JPEGs\landloc_phasediff.jpg")
% saveas(gcf,"SVGs\landloc_phasediff.svg")
% 
% figure(3)
% hold all
% for j=1:1:len_phis
%     plot(datamaster(11,j:len_phis:end),datamaster(2,j:len_phis:end),'LineWidth',1)
% end
% for j=3:2:len_phis
%     text(datamaster(11,j+len_phis*(len_omega-1))-0.02,datamaster(2,j+len_phis*(len_omega-1))+0.1,num2str(datamaster(1,j)),'FontSize',12)
% end
% j=1;
% text(datamaster(11,j+len_phis*(len_omega-1))-0.02,datamaster(2,j+len_phis*(len_omega-1))+0.1,['f=',num2str(datamaster(1,j))],'FontSize',12)
% 
% xlim([0,1])
% set(gca,'fontsize', 18)
% saveas(gcf,"height_phasediff.fig")
% saveas(gcf,"JPEGs\height_phasediff.jpg")
% saveas(gcf,"SVGs\height_phasediff.svg")

% figure(4)
% for j=1:1:len_phis
%     plot(datamaster(7,j:len_phis:end).*datamaster(1,j:len_phis:end),datamaster(2,j:len_phis:end))
% end
% for j=1:2:len_phis
%     text(datamaster(7,j+len_phis*(len_omega-1)).*datamaster(1,j+len_phis*(len_omega-1))-0.02,datamaster(2,j+len_phis*(len_omega-1))+0.1,num2str(datamaster(1,j)))
% end


% figure(2)
% load("datamaster_phis_long.mat")
% for i=2:size(datamaster,2)/2
%     plot(datamaster(1,:),-datamaster(i,:)/9.8)
% end
% legend("")

%% Calculate the boundary of the no flight region
load("datamaster_omegaphis_fine.mat")
[X,Y] = meshgrid(deg2rad(0:5:360),1:0.1:5);  
[len_omega,len_phis]=size(X);
Z = reshape(datamaster(6,:),len_omega,len_phis) ;
phis=deg2rad(0:5:360);
omega=1:0.1:5;
lb=zeros(1,len_phis);
ub=zeros(1,len_phis);
phi=deg2rad(8.178000808-(-48.16109633));
f=3.047;
for i=1:len_phis
    index=max(find(Z(:,i)==0));
    ub(i)=omega(index);
    lb(i)=1;
end


%%
figure(6)
load("datamaster_omegaphis_fine.mat")
datamaster(6,find(datamaster(6,:)==0))=-1;
Z = reshape(datamaster(6,:),len_omega,len_phis) ;
Z=(Z+1).*1.176;

contourf(X,Y,Z,20);
colorbar
hold all
fill([phis,phis(end:-1:1)],[lb,ub(end:-1:1)],[0.5,0.5,0.5],'FaceAlpha',1)
text(pi,1.5,'No Flight','HorizontalAlignment','Center','FontSize',15)
plot([0,2,2,0],[f,f,5,5],'k--','LineWidth',2)
% plot([phi,phi],[1,f],'k--','LineWidth',2)
plot(phi,f,'ko','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(phi,4,'Optimal','HorizontalAlignment','Center','FontSize',15)
xlabel("Phase Difference/Rad")
ylabel("Frequency/Hz")
set(gca,'fontsize', 18)

% Plot max line
ztemp=zeros(size(Z));
max1x=[];
max1y=[];
max2x=[];
max2y=[];
max3x=[];
max3y=[];
max4x=[];
max4y=[];
for i=1:41
    for j=1:73
        if (i>2&&i<42)&&(j>2&&j<72)&&(Z(i,j)>Z(i,j-1)&&Z(i,j)>Z(i,j+1))&&(Z(i,j)>Z(i,j-2)&&Z(i,j)>Z(i,j+2))
            if X(i,j)<1 && Y(i,j)>4
                max1x=[max1x,X(i,j)];
                max1y=[max1y,Y(i,j)];
            else if X(i,j)<3
                max2x=[max2x,X(i,j)];
                max2y=[max2y,Y(i,j)];
                else if X(i,j)>5 && Y(i,j)>4
                        max3x=[max3x,X(i,j)];
                        max3y=[max3y,Y(i,j)];
                    else if Y(i,j)>2.2
                        max4x=[max4x,X(i,j)];
                        max4y=[max4y,Y(i,j)];
                        end
                    end
                end
            end
            ztemp(i,j)=1;
        end
    end
end
% plot(max1x,max1y,'k:','LineWidth',2)
% plot(max2x,max2y,'k:','LineWidth',2)
% plot(max3x,max3y,'k:','LineWidth',2)
% plot(max4x,max4y,'k:','LineWidth',2)

saveas(gcf,"LandLocPD.fig")
saveas(gcf,"JPEGs\LandLocPD.jpg")
saveas(gcf,"SVGs\LandLocPD.svg")
title("Flight distance\zeta(rad)")

%%
figure(7)
Z = reshape(datamaster(11,:),len_omega,len_phis) ;
Z=Z.*211.25;
contourf(X,Y,Z,20);
colorbar
hold all
fill([phis,phis(end:-1:1)],[lb,ub(end:-1:1)],[0.5,0.5,0.5],'FaceAlpha',1)
text(pi,1.5,'No Flight','HorizontalAlignment','Center','FontSize',15)
plot([0,2,2,0],[f,f,5,5],'k--','LineWidth',2)
% plot([phi,phi],[1,f],'k--','LineWidth',2)
plot(phi,f,'ko','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
text(phi,4,'Optimal','HorizontalAlignment','Center','FontSize',15)
xlabel("Phase Difference/Rad")
ylabel("Frequency/Hz")
set(gca,'fontsize', 18)

% Plot max line
ztemp=zeros(size(Z));
max1x=[];
max1y=[];
max2x=[];
max2y=[];
max3x=[];
max3y=[];
max4x=[];
max4y=[];
for i=1:41
    for j=1:73
        if (i>2&&i<42)&&(j>2&&j<72)&&(Z(i,j)>Z(i,j-1)&&Z(i,j)>Z(i,j+1))%&&(Z(i,j)>Z(i,j-2)&&Z(i,j)>Z(i,j+2))
            if X(i,j)<=0.6 && Y(i,j)>4
                max1x=[max1x,X(i,j)];
                max1y=[max1y,Y(i,j)];
            else if X(i,j)<2
                max2x=[max2x,X(i,j)];
                max2y=[max2y,Y(i,j)];
                else if Y(i,j)>3
                        max3x=[max3x,X(i,j)];
                        max3y=[max3y,Y(i,j)];
%                     else if Y(i,j)>2.2
%                         max4x=[max4x,X(i,j)];
%                         max4y=[max4y,Y(i,j)];
%                         end
                    end
                end
            end
            ztemp(i,j)=1;
        end
    end
end
% plot(max1x,max1y,'k:','LineWidth',2)
% plot(max2x,max2y,'k:','LineWidth',2)
% plot(max3x,max3y,'k:','LineWidth',2)
% plot(max4x,max4y,'k:','LineWidth',2)

saveas(gcf,"heightPD.fig")
saveas(gcf,"JPEGs\heightPD.jpg")
saveas(gcf,"SVGs\heightPD.svg")
title("Flight height(cm)")



%%
% figure(9)
% Z = reshape(datamaster(7,:),len_omega,len_phis) ;
% contourf(X,Y,Z,20);
% colorbar
% hold all
% fill([phis,phis(end:-1:1)],[lb,ub(end:-1:1)],[0.5,0.5,0.5],'FaceAlpha',1)
% text(pi,1.5,'No Flight Zone','HorizontalAlignment','Center','FontSize',15)
% plot([0,2,2,0],[f,f,5,5],'k--','LineWidth',2)
% plot([phi,phi],[1,f],'k--','LineWidth',2)
% plot(phi,f,'ko','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
% xlabel("Phase Difference/Rad")
% ylabel("Frequency/Hz")
% set(gca,'fontsize', 18)
% saveas(gcf,"FlightTimePD.fig")
% saveas(gcf,"JPEGs\FlightTimePD.jpg")
% saveas(gcf,"SVGs\FlightTimePD.svg")
% title('Flight time(s)')

%%
figure(10)
Z = reshape(datamaster(8,:)/201,len_omega,len_phis) ;
contourf(X,Y,Z,20);
hold all
fill([phis,phis(end:-1:1)],[lb,ub(end:-1:1)],[0.5,0.5,0.5],'FaceAlpha',1)
text(pi,1.5,'No Flight','HorizontalAlignment','Center','FontSize',15)
% plot([0,2,2,0],[f,f,5,5],'k--','LineWidth',2)
% plot([phi,phi],[1,f],'k--','LineWidth',2)
plot(phi,f,'kp','MarkerSize',8,'MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0])
% plot(phi,f,'ko','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','r','FaceAlpha',0.5)
% text(phi,4,'Optimal','HorizontalAlignment','Center','FontSize',15)
xlabel("Phase Difference/Rad")
ylabel("Frequency/Hz")
set(gca,'fontsize', 18)
colorbar

% Plot max line
ztemp=zeros(size(Z));
max1x=[];
max1y=[];
max2x=[];
max2y=[];
max3x=[];
max3y=[];
max4x=[];
max4y=[];
for i=1:41
    for j=1:73
        if (i>2&&i<42)&&(j>2&&j<72)&&(Z(i,j)>Z(i,j-1)&&Z(i,j)>Z(i,j+1))%&&(Z(i,j)>Z(i,j-2)&&Z(i,j)>Z(i,j+2))
            if X(i,j)<=3
                max1x=[max1x,X(i,j)];
                max1y=[max1y,Y(i,j)];
            else 
                max2x=[max2x,X(i,j)];
                max2y=[max2y,Y(i,j)];
            end
            ztemp(i,j)=1;
        end
    end
end
plot(max1x,max1y,'k:','LineWidth',2)
plot(max2x,max2y,'k:','LineWidth',2)
% plot(max3x,max3y,'k:','LineWidth',2)
% plot(max4x,max4y,'k:','LineWidth',2)

saveas(gcf,"DetachCondPD.fig")
saveas(gcf,"JPEGs\DetachCondPD.jpg")
saveas(gcf,"SVGs\DetachCondPD.svg")
title("Porportion Satisfies the Detachment Condition")
