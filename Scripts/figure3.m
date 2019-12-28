clear all; close all;

load('data_exp1a.mat') ;
load('data_exp1b.mat') ;
load('data_exp1c.mat') ;

figure(1)
d = 0.5; ms = 8;

bl_w = 1:50; comp_w = 51:170;
for condition = 1:3
    % grab data
    Mean_subdata_indv = []; MMean_subdata_indv = []; Mean_subdata = []; Mean_comdata = [];
    if condition == 1 % Exp. 1a
        subdata_indv = endpoint_opt(:, bl_w); subdata =  endpoint_opt(:, comp_w); comdata = endpoint_com_opt(:, comp_w);
    elseif condition == 2 % Exp. 1b
        subdata_indv = endpoint_ave(:, bl_w); subdata =  endpoint_ave(:, comp_w); comdata = endpoint_com_ave(:, comp_w);
    elseif condition == 3 % Exp. 1c
        subdata_indv = endpoint_indv(:, bl_w); subdata =  endpoint_indv(:, comp_w);
    end
    [N T] = size(subdata) ; B = T/10;
    
    % average endpoint in each block
    for i = 1:B
        if i <= 5
            Mean_subdata_indv(:, i) = mean(subdata_indv(:, i*10-9:i*10)')' ;
        end
        Mean_subdata(:, i) = mean(subdata(:, i*10-9:i*10)')' ;
        Mean_comdata(:, i) = mean(comdata(:, i*10-9:i*10)')' ;
    end
    MMean_subdata_indv = mean(Mean_subdata_indv')';
    
    % plot data
    if condition == 1 % Exp. 1a
        subplot(1,3,condition)
        errorbar(1.5, mean(MMean_subdata_indv)/10, std(MMean_subdata_indv/10)/sqrt(N), 'bo', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(Mean_subdata(:,1:B)/10, d, 'b', 'b-', 3:B+2) ;
        seshade(Mean_comdata(:,1:B)/10, d, 'r', 'r-', 3:B+2) ;        
    elseif condition == 2 % Exp. 1b
        subplot(1,3,condition)
        errorbar(1.5, mean(MMean_subdata_indv)/10, std(MMean_subdata_indv/10)/sqrt(N), 'bo', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(Mean_subdata(:,1:B)/10, d, 'b', 'b-', 3:B+2) ;
        seshade(Mean_comdata(:,1:B)/10, d, 'r', 'r-', 3:B+2) ;        
    elseif condition == 3 % Exp. 1c
        subplot(1,3,condition)
        errorbar(1.5, mean(MMean_subdata_indv)/10, std(MMean_subdata_indv/10)/sqrt(N), 'bo', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(Mean_subdata(:,1:B)/10, d, 'b', 'b-', 3:B+2) ;
    end
    xlim([0 15]); ylim([25 29]); yticks(25:1:29); 
    xticks([1.5, 3, 6, 10, 14]); xticklabels(char('BL','1', '4', '8', '12'));
    xlabel('Blocks'); ylabel('Aim point [cm]');
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);
    
end

figure(1)
pos(3) = 1000; pos(4) = 300;
set(gcf, 'Position', pos);







