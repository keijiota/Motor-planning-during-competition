clear all; close all;

load('data_exp1a.mat') ;
load('data_exp1b.mat') ;
load('data_exp1c.mat') ;

% parameters for analytical solution of optimal aim point
a = 70; % gain start point(mm)
b = 300; % gain end point(mm)
c = 100; % maximum reward
alpha = c/(b-a); % slope
beta = -c*a/(b-a); %
mu_seq = linspace(220,300,200); % aimpoint

figure(1)
bl_w = 1:50; comp_w = 51:170;
d = 0.5; ms = 8;

for condition = 1:3
    % grab data
    Std_subdata_indv = []; Std_subdata = [];
    Mean_subdata_indv = []; Mean_subdata = [];
    OptAim_indv = []; OptAim = [];
    RS_indv = []; RS = []; Mean_RS_indv = [];
    if condition == 1 % Exp 1a
        subdata_indv = endpoint_opt(:, bl_w); subdata =  endpoint_opt(:, comp_w);
    elseif condition == 2 % Exp 1b
        subdata_indv = endpoint_ave(:, bl_w); subdata =  endpoint_ave(:, comp_w);
    elseif condition == 3 % Exp 1c
        subdata_indv = endpoint_indv(:, bl_w); subdata =  endpoint_indv(:, comp_w);
    end
    [N T] = size(subdata) ; B = T/10;
    
    Std_subdata_indv = std(subdata_indv(:, :)')' ;
    Mean_subdata_indv = mean(subdata_indv(:, :)')' ;    
    for i = 1:B
        Std_subdata(:, i) = std(subdata(:, i*10-9:i*10)')' ;
        Mean_subdata(:, i) = mean(subdata(:, i*10-9:i*10)')' ;
    end
    
    % analytical solution of optimal aim point    
    % baseline
    for subi = 1:N
        sigma  = Std_subdata_indv(subi,1);
        mu_obs = Mean_subdata_indv(subi,1);
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
        RS_indv(subi, 1) = mu_obs - mu_seq(optAimd);
    end
    
    % competitive or individual task
    for subi = 1:N
        for i = 1:B
            sigma  = Std_subdata(subi,i);
            mu_obs = Mean_subdata(subi,i);
            clear EG
            for tmp = 1:length(mu_seq)
                mu = mu_seq(tmp);
                qfunca = 0.5*erfc(((a-mu)/sigma)/sqrt(2));
                qfuncb = 0.5*erfc(((b-mu)/sigma)/sqrt(2));
                EG(tmp) = (alpha*mu+beta)*(qfunca - qfuncb)...
                    -alpha*sigma/(sqrt(2*pi))*( exp(-0.5./(sigma^2).*(b-mu)^2) - exp(-0.5./(sigma^2).*(a-mu)^2) );
            end
            [optPoint, optAimd] = max(EG);
            OptAim(subi,i)      = mu_seq(optAimd);
            RS(subi, i) = mu_obs - mu_seq(optAimd);
        end
    end
    
    % data plot
    Mean_RS_indv = RS_indv;
    if condition == 1 % Exp 1a
        errorbar(1.5, mean(Mean_RS_indv/10), std(Mean_RS_indv/10)/sqrt(N), 'ro', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(RS(:,1:B)/10, d, 'g', 'g-', 4:B+3) ;         
    elseif condition == 2 % Exp 1b
        errorbar(2, mean(Mean_RS_indv/10), std(Mean_RS_indv/10)/sqrt(N), 'bo', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(RS(:,1:B)/10, d, 'b', 'b-', 4:B+3) ;         
    elseif condition == 3 % Exp 1c
        errorbar(2.5, mean(Mean_RS_indv/10), std(Mean_RS_indv/10)/sqrt(N), 'go', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(RS(:,1:B)/10, d, 'r', 'r-', 4:B+3) ;
    end
end

figure(1)
lineplot(0, 'h', 'k--');
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1);
xlim([0 15.5]); ylim([-1 2.2]); yticks(-1:0.5:2);
xticks([2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
xticklabels(char('BL', '1','','','4','','','','8','','','','12'));
xlabel('Blocks'); ylabel('Risk-sensitivity [cm]'); 
pos(3) = 500; pos(4) = 350;
set(gcf, 'Position', pos);


