clear all; close all;

%% figure 4A-1C'

load('data_exp1a.mat') ;
load('data_exp1b.mat') ;
load('data_exp1c.mat') ;

bl_w = 1:50; comp_w= 51:60; wo_w = 171:180;
blb5 = 41:50;

% parameters for analytical solution of optimal aim point
a = 70; % gain start point(mm)
b = 300; % gain end point(mm)
c = 100; % maximum reward
alpha = c/(b-a); % slope
beta = -c*a/(b-a); %
mu_seq = linspace(220,300,200); % aimpoint

figure(1)
for condition = 1:3
    % grab data
    if condition == 1 % Exp 1a
        endpoint_bl = endpoint_opt(:,bl_w); % endpoints in the baseline
        endpoint_comp = endpoint_opt(:,comp_w); % endpoints in the competition
        endpoint_com = endpoint_com_opt(:, comp_w); % computer's endpoints
        facealpha = 0.2; barcolor = [0.8 0.8 0.8];
    elseif condition == 2 % Exp 1b
        endpoint_bl = endpoint_ave(:,bl_w); endpoint_comp = endpoint_ave(:,comp_w); endpoint_com = endpoint_com_ave(:, comp_w);
        facealpha = 0.2; barcolor = [0.8 0.8 0.8];
    elseif condition == 3 % Exp 1c
        endpoint_bl = endpoint_indv(:,bl_w); endpoint_comp = endpoint_indv(:,comp_w);
        facealpha = 0.5; barcolor = [0.5 0.5 0.5];
    end
    meanendpoint_bl_indv = mean(endpoint_bl')';
    meanendpoint_bl = mean(meanendpoint_bl_indv);
    mean_comp = mean(endpoint_comp(:,1:5)')' ;
    
    % analytical solution of optimal aim point
    OptAim_indv = [];
    [N T] = size(endpoint_bl) ;
    std_endpoint_bl = std(endpoint_bl(:,:)')';
    for subi = 1:N
        sigma  = std_endpoint_bl(subi);
        clear EG
        for tmp = 1:length(mu_seq)
            mu = mu_seq(tmp);
            qfunca = 0.5*erfc(((a-mu)/sigma)/sqrt(2));
            qfuncb = 0.5*erfc(((b-mu)/sigma)/sqrt(2));
            EG(tmp) = (alpha*mu+beta)*(qfunca - qfuncb)...
                -alpha*sigma/(sqrt(2*pi))*( exp(-0.5./(sigma^2).*(b-mu)^2) - exp(-0.5./(sigma^2).*(a-mu)^2) );
        end
        [optPoint, optAimd] = max(EG);
        OptAim_indv(subi,1) = mu_seq(optAimd);
    end
    
    % time-series
    subplot(3,100,condition*100-99:condition*100-45)
    seshade(endpoint_bl/10, 0.5, 'k', 'k-', bl_w); hold on
    lineplot(meanendpoint_bl/10, 'h', 'k-', 'linewidth', 1.5);
    seshade(endpoint_comp/10, facealpha, 'k', 'k-', comp_w);
    xlim([0 61]); ylim([25 30]);
    yticks(25:1:30); xticks(0:10:60);
    ylabel('Endpoint [cm]'); xlabel('Trials');
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);
    
    % bar graph
    subplot(3,100,condition*100-35:condition*100)
    bar(1,mean(OptAim_indv)/10, 0.8, 'FaceColor',[0.5 0.5 0.5]); hold on
    bar(2,mean(meanendpoint_bl_indv)/10, 0.8,  'FaceColor',[0.5 0.5 0.5]);
    bar(3,mean(mean_comp)/10, 0.8,  'FaceColor',barcolor);
    
    % individual plot
    ms = 5;
    plot(1, OptAim_indv/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);
    plot(2, meanendpoint_bl_indv/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);
    plot(3, mean_comp/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);
    
    errorbar(1, mean(OptAim_indv)/10, std(OptAim_indv/10)/sqrt(N), 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',6) ;
    errorbar(2, mean(meanendpoint_bl_indv)/10, std(meanendpoint_bl_indv/10)/sqrt(N), 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',6) ;
    errorbar(3, mean(mean_comp)/10, std(mean_comp/10)/sqrt(N), 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',6) ;
    
    ylim([25 30]); xlim([0 4]);
    yticks(25:1:30); xticks(0:1:4);
    xticklabels(char('', 'Opt.', 'Obs.','1st-5th'));
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);
    
end
figure(1)
pos(3) = 550; pos(4) = 950;
set(gcf, 'Position', pos);

%% figure 4D-F'

load('data_exp1c.mat') ;
load('data_exp2.mat') ;
load('data_exp3.mat') ;

figure(2)
for condition = 1:4
    % grab data
    if condition == 1 % Exp 1c
        bl_w = 171:220; comp_w = 221:230;
        endpoint_bl = endpoint_indv(:,bl_w); endpoint_comp = endpoint_indv(:,comp_w);
        barcolor1 = [0.5 0.5 0.5];  barcolor2 = [0.8 0.8 0.8];
    elseif condition == 2 % Exp 2
        bl_w = 1:50; comp_w = 51:60;
        endpoint_bl = endpoint_presen(:,bl_w); endpoint_comp = endpoint_presen(:,comp_w);
        endpoint_computer = endpoint_com_presen(:, 241:250);
        barcolor1 = [0.5 0.5 0.5];  barcolor2 = [0.8 0.8 0.8];
    elseif condition == 3 % Exp 3 (Preceding subject)
        bl_w = 1:50; comp_w= 51:60;
        endpoint_bl = [endpoint_P1_2p(:, bl_w)] ; endpoint_comp = [endpoint_P1_2p(:, comp_w)] ;
        barcolor1 = [0.5 0 0];  barcolor2 = [0.8 0 0];
    elseif condition == 4 % Exp 3 (Following subject)
        bl_w = 1:50; comp_w= 51:60;
        endpoint_bl = [endpoint_P2_2p(:, bl_w)] ; endpoint_comp = [endpoint_P2_2p(:, comp_w)] ;
    end
    meanendpoint_bl_indv = mean(endpoint_bl')';
    meanendpoint_bl = mean(meanendpoint_bl_indv);
    mean_comp = mean(endpoint_comp(:,1:5)')' ;
    
    % analytical solution of optimal aim point
    OptAim_indv = [];
    [N T] = size(endpoint_bl) ;
    std_endpoint_bl = std(endpoint_bl(:,:)')';
    for subi = 1:N
        sigma  = std_endpoint_bl(subi);
        clear EG
        for tmp = 1:length(mu_seq)
            mu = mu_seq(tmp);
            qfunca = 0.5*erfc(((a-mu)/sigma)/sqrt(2));
            qfuncb = 0.5*erfc(((b-mu)/sigma)/sqrt(2));
            EG(tmp) = (alpha*mu+beta)*(qfunca - qfuncb)...
                -alpha*sigma/(sqrt(2*pi))*( exp(-0.5./(sigma^2).*(b-mu)^2) - exp(-0.5./(sigma^2).*(a-mu)^2) );
        end
        [optPoint, optAimd] = max(EG);
        OptAim_indv(subi,1) = mu_seq(optAimd);
    end
    
    % time-series
    if condition == 1
        subplot(3,100,condition*100-99:condition*100-45)
        seshade(endpoint_bl/10, 0.5, 'k', 'k-', bl_w); hold on
        lineplot(meanendpoint_bl/10, 'h', 'k-', 'linewidth', 1.5);
        seshade(endpoint_comp/10, 0.2, 'k', 'k-', comp_w);
        xlim([170 231]); ylim([25 30]);
        yticks(25:1:30); xticks(170:10:230);
    elseif condition == 2
        subplot(3,100,condition*100-99:condition*100-45)
        seshade(endpoint_bl/10, 0.5, 'k', 'k-', bl_w); hold on
        lineplot(meanendpoint_bl/10, 'h', 'k-', 'linewidth', 1.5);
        seshade(endpoint_computer/10, 0.5, 'r', 'r-', 51:60);
        seshade(endpoint_comp/10, 0.2, 'k', 'k-', 61:70);
        xlim([0 71]); ylim([24 30]);
        yticks(24:1:30); xticks(0:10:70);
    elseif condition == 3
        subplot(3,100,condition*100-99:condition*100-45)
        seshade(endpoint_bl/10, 0.5, 'b', 'b-', bl_w); hold on
        lineplot(meanendpoint_bl/10, 'h', 'b-', 'linewidth', 1.5);
        seshade(endpoint_comp/10, 0.2, 'b', 'b-', comp_w);
    elseif condition == 4
        subplot(3,100,3*100-99:3*100-45)
        seshade(endpoint_bl/10, 0.5, 'r', 'r-', bl_w); hold on
        lineplot(meanendpoint_bl/10, 'h', 'r-', 'linewidth', 1.5);
        seshade(endpoint_comp/10, 0.2, 'r', 'r-', comp_w);
        xlim([0 61]); ylim([25 30]);
        yticks(25:1:30); xticks(0:10:60);
    end
    ylabel('Endpoint [cm]'); xlabel('Trials');
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);
    
    if condition <=3
        if condition <= 3
            subplot(3,100,condition*100-35:condition*100)
        else
            subplot(3,100,3*100-35:3*100)
        end
        % bar graph
        bar(1,mean(OptAim_indv)/10, 0.8, 'FaceColor',barcolor1); hold on
        bar(2,mean(meanendpoint_bl_indv)/10, 0.8, 'FaceColor',barcolor1);
        bar(3,mean(mean_comp)/10, 0.8, 'FaceColor',barcolor2);
        
        % individual plot
        ms = 5;
        plot(1, OptAim_indv/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);
        plot(2, meanendpoint_bl_indv/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);
        plot(3, mean_comp/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);
        
        errorbar(1, mean(OptAim_indv)/10, std(OptAim_indv/10)/sqrt(N), 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',6) ;
        errorbar(2, mean(meanendpoint_bl_indv)/10, std(meanendpoint_bl_indv/10)/sqrt(N), 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',6) ;
        errorbar(3, mean(mean_comp)/10, std(mean_comp/10)/sqrt(N), 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',6) ;
        
        if condition == 2
            ylim([24 30]);
        else
            ylim([25 30]);
        end
        xlim([0 4]);
        yticks(24:1:30); xticks(0:1:5);
        xticklabels(char('', 'Opt.', 'Obs.','1st-5th'));
        set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);
    end
end
figure(2)
pos(3) = 550; pos(4) = 950;
set(gcf, 'Position', pos);

