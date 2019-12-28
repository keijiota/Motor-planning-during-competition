clear all; close all;

load('data_exp5a');
load('data_exp5b');

bl_w = 1:50; comp_w = 51:170;
figure(1)
d = 0.5; ms = 8;

for condition = 1:2
    % grab data
    Mean_subdata_indv = []; MMean_subdata_indv = []; Mean_subdata = []; Mean_comdata = []; 
    if condition == 1 % Exp. 5a
        subdata_indv = endpoint_observation(:, bl_w); 
        subdata = endpoint_observation(:, comp_w);
        comdata = endpoint_com_observation(:, comp_w);
    elseif condition == 2  % Exp. 5b
        subdata_indv = endpoint_threshold(:, bl_w); 
        subdata = endpoint_threshold(:, comp_w);
    end
    
    % average endpoint in each block    
    [N T] = size(subdata) ; B = T/10;
    for i = 1:B
        if i <= 5
            Mean_subdata_indv(:, i) = mean(subdata_indv(:, i*10-9:i*10)')' ;
        end
        Mean_subdata(:, i) = mean(subdata(:, i*10-9:i*10)')' ;
        Mean_comdata(:, i) = mean(comdata(:, i*10-9:i*10)')' ;
    end    
    MMean_subdata_indv = mean(Mean_subdata_indv')';
    
    % plot data
    if condition == 1  % Exp. 5a
        subplot(1,2,1)
        errorbar(1.5, mean(MMean_subdata_indv)/10, std(MMean_subdata_indv/10)/sqrt(N), 'bo', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(Mean_subdata(:,1:B)/10, d, 'b', 'b-', 3:B+2) ;
        seshade(Mean_comdata(:,1:B)/10, d, 'r', 'r-', 3:B+2) ;
        lineplot(mean(MMean_subdata_indv)/10, 'h', 'k--');
        
    elseif condition == 2  % Exp. 5b
        subplot(1,2,2)
        errorbar(1.5, mean(MMean_subdata_indv)/10, std(MMean_subdata_indv/10)/sqrt(N), 'bo', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
        seshade(Mean_subdata(:,1:B)/10, d, 'b', 'b-', 3:B+2) ;     
        lineplot(mean(MMean_subdata_indv)/10, 'h', 'k--'); 
    end
        set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);
        ylabel('Aim point [cm]'); xlabel('Blocks');         
        xlim([0 15]); ylim([24.8 29]); yticks(24:1:29); xticks([1.5, 3, 6, 10, 14]); xticklabels(char('BL','1', '4', '8', '12'));    
end

figure(1)
pos(3) = 625; pos(4) = 320;
set(gcf, 'Position', pos);




