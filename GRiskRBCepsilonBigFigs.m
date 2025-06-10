clear; clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline: alpha=1/3, beta=.99, rho=.9, delta=.025


fontS=14;

load GRiskCESRCBigPCnonexpB_2024_02_13_21_19-8k51.mat;

sigmavect'
alpha
beta
delta
rho
A

S_SDZt1=S_SDZ;
S_CSDt1=S_CSD;
S_KSDt1=S_KSD;
S_ISDt1=S_ISD;
S_YSDt1=S_YSD;
S_V1t1=S_V1;
S_V3t1=S_V3;
S_CBARt1=S_CBAR;
S_KBARt1=S_KBAR;
S_IBARt1=S_IBAR;
S_YBARt1=S_YBAR;

clear S_SDZ S_CSD S_KSD S_ISD S_YSD S_V1 S_V3 S_CBAR S_KBAR S_IBAR S_YBAR ;


load GRiskCESRCBigPCnonexpB_2024_02_21_19_15-8k51b.mat;

sigmavect'
alpha
beta
delta
rho
A


S_SDZt2=S_SDZ;
S_CSDt2=S_CSD;
S_KSDt2=S_KSD;
S_ISDt2=S_ISD;
S_YSDt2=S_YSD;
S_V1t2=S_V1;
S_V3t2=S_V3;
S_CBARt2=S_CBAR;
S_KBARt2=S_KBAR;
S_IBARt2=S_IBAR;
S_YBARt2=S_YBAR;

clear S_SDZ S_CSD S_KSD S_ISD S_YSD S_V1 S_V3 S_CBAR S_KBAR S_IBAR S_YBAR ;


S_SDZ=nan(1,5,6);
S_CSD=nan(1,5,6);
S_KSD=nan(1,5,6);
S_ISD=nan(1,5,6);
S_YSD=nan(1,5,6);
S_V1=nan(1,5,6);
S_V3=nan(1,5,6);
S_CBAR=nan(1,5,6);
S_KBAR=nan(1,5,6);
S_IBAR=nan(1,5,6);
S_YBAR=nan(1,5,6);


for i=1:6
    for j=3:5 % Only up to .1
        S_SDZ(:,2*(j-1)+1-4,i)=S_SDZt1(:,j,i);
        S_CSD(:,2*(j-1)+1-4,i)=S_CSDt1(:,j,i);
        S_KSD(:,2*(j-1)+1-4,i)=S_KSDt1(:,j,i);
        S_ISD(:,2*(j-1)+1-4,i)=S_ISDt1(:,j,i);
        S_YSD(:,2*(j-1)+1-4,i)=S_YSDt1(:,j,i);
        S_V1(:,2*(j-1)+1-4,i)=S_V1t1(:,j,i);
        S_V3(:,2*(j-1)+1-4,i)=S_V3t1(:,j,i);
        S_CBAR(:,2*(j-1)+1-4,i)=S_CBARt1(:,j,i);
        S_KBAR(:,2*(j-1)+1-4,i)=S_KBARt1(:,j,i);
        S_IBAR(:,2*(j-1)+1-4,i)=S_IBARt1(:,j,i);        
        S_YBAR(:,2*(j-1)+1-4,i)=S_YBARt1(:,j,i);
    end

    for j=3:4 % Only up to .1
        S_SDZ(:,2*j-4,i)=S_SDZt2(:,j,i);
        S_CSD(:,2*j-4,i)=S_CSDt2(:,j,i);
        S_KSD(:,2*j-4,i)=S_KSDt2(:,j,i);
        S_ISD(:,2*j-4,i)=S_ISDt2(:,j,i);
        S_YSD(:,2*j-4,i)=S_YSDt2(:,j,i);
        S_V1(:,2*j-4,i)=S_V1t2(:,j,i);
        S_V3(:,2*j-4,i)=S_V3t2(:,j,i);
        S_CBAR(:,2*j-4,i)=S_CBARt2(:,j,i);
        S_KBAR(:,2*j-4,i)=S_KBARt2(:,j,i);
        S_IBAR(:,2*j-4,i)=S_IBARt2(:,j,i);        
        S_YBAR(:,2*j-4,i)=S_YBARt2(:,j,i);
    end
end


lsdzvect = size(S_SDZ,2);
lsigmavect = length(sigmavect);



figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
tiledlayout(2,2,'TileSpacing','compact');

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_CSD(:,:,i)./S_CBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)./S_CBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.1 8.7]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
legend(labels, 'Position', [.42 .6 .07 .3/2],'FontSize',fontS+2)
title('Consumption SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_KSD(:,:,i)./S_KBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, S_KSD(:,1,i)./S_KBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.2 18]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Location','northwest','FontSize',fontS+2)
title('Capital SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_YSD(:,:,i)./S_YBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, S_YSD(:,1,i)./S_YBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.15 13.5]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Location','northwest','FontSize',fontS+2)
title('GDP  SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_ISD(:,:,i)./S_IBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    % text(S_SDZ(:,1,i)*1.01, S_ISD(:,1,i)./S_IBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.4 35]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Position', [.42 .2 .07 .3/2],'FontSize',fontS+2)
title('Investment SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

exportgraphics(gcf,"GRiskRBCepsilonBigrSD.eps","Resolution",300)





figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
tiledlayout(2,2,'TileSpacing','compact');

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_CSD(:,:,i)/cssTh*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)/cssTh*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
text(.025, 16, labels{i,1},'color',colors{i,1},'FontSize',12);
axis([-.0015 .1 -.1 8.7]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
legend(labels, 'Position', [.42 .6 .07 .3/2],'FontSize',fontS+2)
title('Consumption SD/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_KSD(:,:,i)./S_KBAR(:,lsdzvect,i)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, S_KSD(:,1,i)./S_KBAR(:,lsdzvect,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.2 18]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Location','northwest','FontSize',fontS+2)
title('Capital SD/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_YSD(:,:,i)/(A*kssTh^alpha)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, S_YSD(:,1,i)/(A*kssTh^alpha)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.15 13.5]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Location','northwest','FontSize',fontS+2)
title('GDP SD/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_ISD(:,:,i)./S_IBAR(:,lsdzvect,i)*100, 'color',colors{i,1},'linewidth',1.5);
    % text(S_SDZ(:,1,i)*1.01, S_ISD(:,1,i)./S_IBAR(:,lsdzvect,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.4 35]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
title('Investment SD/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;


exportgraphics(gcf,"GRiskRBCepsilonBigSD.eps","Resolution",300)





figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
tiledlayout(2,2,'TileSpacing','compact');

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_CBAR(:,:,i)/S_CBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, (S_CBAR(:,1,i)/S_CBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.005 .4]);
ytickformat('%.2f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
title('Consumption AVG/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_KBAR(:,:,i)/S_KBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, (S_KBAR(:,1,i)/S_KBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.05 4.99]);
ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
legend(labels, 'Position', [.57 .65 .07 .3/2],'FontSize',fontS+2)
title('Capital and investment AVG/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_YBAR(:,:,i)/S_YBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, (S_YBAR(:,1,i)/S_YBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.001 .1 -.01 1.55]);
ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Position', [.2 .55 .07 .3],'FontSize',fontS+2)
title('GDP AVG/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

S_ZBAR=S_YBAR./(A*S_KBAR.^alpha);
nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_ZBAR(:,:,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, (S_ZBAR(:,1,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.001 .1 -.16 .2]);
ytickformat('%.2f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Position', [.2 .55 .07 .3],'FontSize',fontS+2)
title('TFP AVG/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

exportgraphics(gcf,"GRiskRBCepsilonBigM.eps","Resolution",300)




% figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_ISD(:,:,i)./S_IBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
%     % text(S_SDZ(:,1,i)*1.01, S_ISD(:,1,i)./S_IBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.4 35]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
% title('Investment SD/AVG','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_KSD(:,:,i)./S_KBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_KSD(:,1,i)./S_KBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.2 18]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% %legend(labels, 'Location','northwest','FontSize',fontS+2)
% title('Capital SD/AVG','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% exportgraphics(gcf,"GRiskRBCepsilonBigrSDyi.eps","Resolution",300)






% figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_YBAR(:,:,i)/S_YBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     %text(S_SDZ(:,1,i)*1.01, (S_YBAR(:,1,i)/S_YBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.001 .1 -.01 1.55]);
% ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Position', [.2 .55 .07 .3],'FontSize',fontS+2)
% title('GDP AVG/SS','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% S_ZBAR=S_YBAR./(A*S_KBAR.^alpha);
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_ZBAR(:,:,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, (S_ZBAR(:,1,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.001 .1 -.16 .2]);
% ytickformat('%.2f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Position', [.2 .55 .07 .3],'FontSize',fontS+2)
% title('TFP AVG/SS','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% exportgraphics(gcf,"GRiskRBCepsilonBigMy.eps","Resolution",300)




figure;
hold on;
for i=1:lsigmavect
    plot(S_CSD(:,:,i), S_CBAR(:,:,i)-S_CBAR(:,lsdzvect,i), 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)./S_CBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.01 .7 -.0005 .095]);
ytickformat('%.2f')
xlabel('SD','FontSize',fontS+2);
ylabel('AVG','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
title('Consumption AVG vs SD','Interpreter','Latex','FontSize',fontS+4);
hold off;


% CSD=nan(lsdzvect,lsigmavect);
% CBAR=nan(lsdzvect,lsigmavect);
% for i=1:lsigmavect
%     for j=1:lsdzvect
%         CSD(j,i)=S_CSD(1,j,i);
%         CBAR(j,i)=S_CBAR(1,j,i);
%     end
% end
% figure;
% hold on;
% for j=1:lsdzvect-1
%     plot(CSD(j,:)/cssTh*100, (CBAR(j,:)/cssTh-1)*100, 'color',colors{j,1},'linewidth',1.5);
%     %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)./S_CBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% %axis([-.01 .7 -.0005 .095]);
% %ytickformat('%.2f')
% xlabel('SD','FontSize',fontS+2);
% ylabel('AVG','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% %centeraxes(gca);
% %legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
% title('Consumption AVG vs SD','Interpreter','Latex','FontSize',fontS+4);
% hold off;




% figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% nexttile;
% hold on;
% cssGainPercentVkss=nan(lsigmavect,lsdzvect);
% for i=1:lsigmavect
%     for j=1:lsdzvect
%         cssGainPercentVkss(i,j)=((1+(1-beta)*(1-sigmavect(i))*(S_V1(1,j,i)-S_V1(1,lsdzvect,i))*cssTh^(sigmavect(i)-1))^(1/(1-sigmavect(i)))-1)*100;
%     end
% end
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), cssGainPercentVkss(i,:), 'color',colors{i,1},'linewidth',1.5);
% end
% for i=[1:4 lsigmavect]
%     text(S_SDZ(:,1,i)*1.01, cssGainPercentVkss(i,1), labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% i=5;
% text(S_SDZ(:,1,i), cssGainPercentVkss(i,1)-.1, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% axis([-.003 .2 -3.8 1.2]);
% ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% title('$x\%$ based on $v(k_{ss},E(z))$','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% nexttile;
% hold on;
% cssGainPercentVkbar=nan(lsigmavect,lsdzvect);
% for i=1:lsigmavect
%     for j=1:lsdzvect
%         cssGainPercentVkbar(i,j)=((1+(1-beta)*(1-sigmavect(i))*(S_V3(1,j,i)-S_V3(1,lsdzvect,i))*cssTh^(sigmavect(i)-1))^(1/(1-sigmavect(i)))-1)*100;
%     end
% end
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), cssGainPercentVkbar(i,:), 'color',colors{i,1},'linewidth',1.5);
% end
% for i=[3 4 lsigmavect]
%     text(S_SDZ(:,1,i)*1.01, cssGainPercentVkbar(i,1), labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% i=1;
% text(S_SDZ(:,1,i)*1.01, cssGainPercentVkbar(i,1)-.05, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% i=2;
% text(S_SDZ(:,1,i)*1.01, cssGainPercentVkbar(i,1)+.04, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% i=5;
% text(S_SDZ(:,1,i)*1.03, cssGainPercentVkbar(i,1)+.05, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% axis([-.003 .2 -3.8 1.2]);
% ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% title(['$x\%$ based on $v(\bar{k},E(z))$'],'Interpreter','Latex','FontSize',fontS+4);
% hold off;
% exportgraphics(gcf,"GRiskRBCepsilonvksskbar.eps","Resolution",300)






% clearvars -except fontS;
% load GRiskCESRCBigHPPC_2024_03_03_12_25-4k51hp.mat;
% 
% figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_CSDCOMPDATA(:,:,i), 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_CSDCOMPDATA(:,1,i), labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% text(.01,3.5, labels{i,1},'color',colors{i,1},'FontSize',12);
% plot(S_SDZ(:,:,i), 1.3856*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle',':');
% axis([-.0015 .2 -.03 3.5]);
% ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% title('Consumption SD HP-filtered','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_YSDCOMPDATA(:,:,i), 'color',colors{i,1},'linewidth',1.5);
%     %text(S_SDZ(:,1,i)*1.01, S_YSDCOMPDATA(:,1,i), labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% plot(S_SDZ(:,:,i), 1.6443*ones(lsdzvect),'color',[.5 .5 .5],'LineStyle',':');
% axis([-.0015 .2 -.1 10.8]);
% % ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Location','northwest','FontSize',fontS+2)
% title('GDP SD HP-filtered','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% exportgraphics(gcf,"GRiskRBCepsilonBigHPSDcy.eps","Resolution",300)









clearvars -except fontS;
 
load GRiskCESRCBigPCB_2025_02_25_23_17-8k51a.mat;


figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
tiledlayout(2,2,'TileSpacing','compact');

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_CSD(:,:,i)./S_CBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)./S_CBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.25 15.3]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
legend(labels, 'Position', [.42 .6 .07 .3/2],'FontSize',fontS+2)
title('Consumption SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_KSD(:,:,i)./S_KBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, S_KSD(:,1,i)./S_KBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.4 22]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
title('Capital SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_YSD(:,:,i)./S_YBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, S_YSD(:,1,i)./S_YBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.25 26]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Location','northwest','FontSize',fontS+2)
title('GDP  SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), S_ISD(:,:,i)./S_IBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
    % text(S_SDZ(:,1,i)*1.01, S_ISD(:,1,i)./S_IBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.4 26]);
%ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
title('Investment SD/AVG','Interpreter','Latex','FontSize',fontS+4);
hold off;

exportgraphics(gcf,"GRiskRBCepsilonBigarSD.eps","Resolution",300)






% figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_CSD(:,:,i)./S_CBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
%     %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)./S_CBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.25 15.3]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
% title('Consumption SD/AVG','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_YSD(:,:,i)./S_YBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_YSD(:,1,i)./S_YBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.25 26]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% %legend(labels, 'Location','northwest','FontSize',fontS+2)
% title('GDP  SD/AVG','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% exportgraphics(gcf,"GRiskRBCepsilonBigarSDcy.eps","Resolution",300)





% figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_CSD(:,:,i)/cssTh*100, 'color',colors{i,1},'linewidth',1.5);
%     %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)/cssTh*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.25 17.5]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
% title('Consumption SD/SS','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_YSD(:,:,i)/(A*kssTh^alpha)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_YSD(:,1,i)/(A*kssTh^alpha)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.25 26]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% title('GDP SD/SS','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% exportgraphics(gcf,"GRiskRBCepsilonBigaSDcy.eps","Resolution",300)





figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
tiledlayout(2,2,'TileSpacing','compact');

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_CBAR(:,:,i)/S_CBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, (S_CBAR(:,1,i)/S_CBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.05 4.9]);
ytickformat('%.1f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
title('Consumption AVG/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_KBAR(:,:,i)/S_KBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, (S_KBAR(:,1,i)/S_KBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.0015 .1 -.15 10.2]);
%ytickformat('%.0f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
legend(labels, 'Position', [.57 .65 .07 .3/2],'FontSize',fontS+2)
title('Capital and investment AVG/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_YBAR(:,:,i)/S_YBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, (S_YBAR(:,1,i)/S_YBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.001 .1 -.1 8]);
%ytickformat('%.0f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
title('GDP AVG/SS','Interpreter','Latex','FontSize',fontS+4);
hold off;

S_ZBAR=S_YBAR./(A*S_KBAR.^alpha);
nexttile;
hold on;
for i=1:lsigmavect
    plot(S_SDZ(:,:,i), (S_ZBAR(:,:,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
    text(S_SDZ(:,1,i)*1.01, (S_ZBAR(:,1,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.001 .1 -.16 .2]);
ytickformat('%.2f')
xlabel('\sigma','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
title('TFP AVG/SS','Interpreter','Latex','FontSize',fontS+4);

exportgraphics(gcf,"GRiskRBCepsilonBigaM.eps","Resolution",300)





% figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_ISD(:,:,i)./S_IBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
%     % text(S_SDZ(:,1,i)*1.01, S_ISD(:,1,i)./S_IBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.4 26]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
% title('Investment SD/AVG','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), S_KSD(:,:,i)./S_KBAR(:,:,i)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, S_KSD(:,1,i)./S_KBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.0015 .1 -.4 22]);
% %ytickformat('%.1f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% title('Capital SD/AVG','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% exportgraphics(gcf,"GRiskRBCepsilonBigarSDyi.eps","Resolution",300)






% figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
% tiledlayout(1,2,'TileSpacing','compact');
% 
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_YBAR(:,:,i)/S_YBAR(:,lsdzvect,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, (S_YBAR(:,1,i)/S_YBAR(:,lsdzvect,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.001 .1 -.1 8]);
% %ytickformat('%.0f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% %legend(labels, 'Location','northwest','FontSize',fontS+2)
% title('GDP AVG/SS','Interpreter','Latex','FontSize',fontS+4);
% hold off;
% 
% S_ZBAR=S_YBAR./(A*S_KBAR.^alpha);
% nexttile;
% hold on;
% for i=1:lsigmavect
%     plot(S_SDZ(:,:,i), (S_ZBAR(:,:,i)-1)*100, 'color',colors{i,1},'linewidth',1.5);
%     text(S_SDZ(:,1,i)*1.01, (S_ZBAR(:,1,i)-1)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
% end
% axis([-.001 .1 -.16 .2]);
% ytickformat('%.2f')
% xlabel('\sigma','FontSize',fontS+2);
% fontsize(gca,fontS,"points")
% centeraxes(gca);
% legend(labels, 'Position', [.2 .55 .07 .3],'FontSize',fontS+2)
% title('TFP AVG/SS','Interpreter','Latex','FontSize',fontS+4);
% 
% exportgraphics(gcf,"GRiskRBCepsilonBigaMy.eps","Resolution",300)




figure;
hold on;
for i=1:lsigmavect
    plot(S_CSD(:,:,i), S_CBAR(:,:,i)-S_CBAR(:,lsdzvect,i), 'color',colors{i,1},'linewidth',1.5);
    %text(S_SDZ(:,1,i)*1.01, S_CSD(:,1,i)./S_CBAR(:,1,i)*100, labels{i,1},'color',colors{i,1},'FontSize',fontS);
end
axis([-.01 .7 -.0005 .095]);
ytickformat('%.2f')
xlabel('SD','FontSize',fontS+2);
ylabel('AVG','FontSize',fontS+2);
fontsize(gca,fontS,"points")
centeraxes(gca);
%legend(labels, 'Position', [.42 .2 .07 .3],'FontSize',fontS+2)
title('Consumption AVG vs SD','Interpreter','Latex','FontSize',fontS+4);
hold off;


