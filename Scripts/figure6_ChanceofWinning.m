clear all; close all;

load('data_exp1a.mat') ;
load('data_exp1b.mat') ;

bl_w = 1:50; comp_w = 51:170;
for condition = 1:2
    % grab data
    Mean_subdata_indv = []; Mean_subdata = []; Mean_comdata = [];
    sigma_subdata = [];
    
    if condition == 1
        subdata_indv =  endpoint_opt(:, bl_w); subdata =  endpoint_opt(:, comp_w); comdata = endpoint_com_opt(:, comp_w);
        sigma = sigma_com_in_opt(:, 6:17);
    elseif condition == 2
        subdata_indv =  endpoint_ave(:, bl_w); subdata =  endpoint_ave(:, comp_w); comdata = endpoint_com_ave(:, comp_w);
        sigma = sigma_com_in_ave(:, 6:17);
    end
    [N T] = size(subdata) ; B = T/10;
    
    % average endpoint in each block    
    for ib = 1:B
        if ib <= 5
            Mean_subdata_indv(:, ib) = mean(subdata_indv(:, ib*10-9:ib*10)')' ;
        end
        Mean_subdata(:, ib) = mean(subdata(:, ib*10-9:ib*10)')' ;
        Mean_comdata(:, ib) = mean(comdata(:, ib*10-9:ib*10)')' ;
        sigma_subdata(:, ib) = std(subdata(:, ib*10-9:ib*10)')' ;
    end
    Mean_subdata_bl = mean(Mean_subdata_indv')' ;
    
    if condition == 1
        Mean_subdata_bl_opt  = Mean_subdata_bl ;
        Mean_subdata_opt = Mean_subdata;
        Mean_comdata_opt = Mean_comdata;
        Sigma_opt = sigma;
    elseif condition == 2
        Mean_subdata_bl_ave  = Mean_subdata_bl ;
        Mean_subdata_ave = Mean_subdata;
        Mean_comdata_ave = Mean_comdata;
        Sigma_ave = sigma;
    end
end

% Ntrial = 100000; % the number of simulations
Ntrial = 100; % the number of simulations

subii = 0; subN = [9 8];
for condition = 1:2
    for subi = 1:subN(condition)
        subii = subii + 1;
        win_rate = [];
        for ib = 2:12
            if condition == 1
                aim_sub = Mean_subdata_opt(subi, ib);
                aim_com = Mean_comdata_opt(subi, ib);
                sigma_sub = Sigma_opt(subi, ib); % this sigma was subject's sigma for previous 40 trials and was for calculation of risk-sensitivity
                sigma_com = Sigma_opt(subi, ib); % this sigma was also used for opponent input sigma
                aim_bl = Mean_subdata_bl_opt(subi);
                diff_com_bl = aim_com - aim_bl ;
                diff_sub_bl = aim_sub - aim_bl ;
            elseif condition == 2
                aim_sub = Mean_subdata_ave(subi, ib);
                aim_com = Mean_comdata_ave(subi, ib);
                sigma_sub = Sigma_ave(subi, ib);
                sigma_com = Sigma_ave(subi, ib);
                aim_bl = Mean_subdata_bl_ave(subi);
                diff_com_bl = aim_com - aim_bl ;
                diff_sub_bl = aim_sub - aim_bl ;
            end
            
            aimi = 0;
            for diff_aim = -40:0.5:20 % simulational aim point
                aim_sub = diff_aim + aim_bl;
                aimi = aimi + 1;
               % simulation of competitive game for N repetitions               
                for itrial = 1:Ntrial
                    endpoint_sub = aim_sub + randn(10,1)*sigma_sub ; % simulational aim point with movement noise
                    endpoint_com = aim_com + randn(10,1)*sigma_com ; % computer's aim point with movement noise
                    % calculate score given endpoint
                    score_sub = (endpoint_sub - 70) / 230 * 100 ;
                    score_sub(score_sub>100) = 0 ;
                    score_com = (endpoint_com - 70) / 230 * 100 ;
                    score_com(score_com>100) = 0 ;
                    % total score in 10 trials
                    TS_sub(itrial,1) = sum(score_sub) ; 
                    TS_com(itrial,1) = sum(score_com) ;
                end
                % chance of winning in each aim point
                win_rate(aimi, ib-1) = sum(TS_sub > TS_com) / Ntrial ;
            end
            Diff_com_bl(ib-1) = diff_com_bl ;
        end
        
        [sDiff_com_bl, tmp] = sort(Diff_com_bl) ;
        swin_rate = win_rate(:,tmp) ;
        
        % store data across subjects and blocks
        mWin_rate(:,:,subii) = swin_rate;
        mDiff_com_bl(:,:,subii) = sDiff_com_bl;
    end
end

[c r z] = size(mDiff_com_bl) ;
bin = -40:2.5:20; % in steps of 2.5 mm
[n N] = size(bin);
[wc wr wz] = size(mWin_rate) ;

ddiff_com_bl = NaN(1, 50, N-1);
wwin_rate = NaN(wc, 50, N-1);

% rearrange the data in each bin
for bini = 1:N-1
    ii = 0;
    for ri = 1:r
        for zi = 1:z
            if mDiff_com_bl(1, ri, zi) > bin(bini) && mDiff_com_bl(1, ri, zi) <= bin(bini+1)
                ii = ii +1;
                ddiff_com_bl(1, ii, bini) = mDiff_com_bl(1, ri, zi);
                wwin_rate(:, ii, bini) = mWin_rate(:, ri, zi);
            end
        end
    end
end

% average chance of winning in each bin
for bini = 1:N-1
    mean_win_rate(:,bini) = nanmean(wwin_rate(:,:,bini)')' ;
end

mean_win_rate = flipud(mean_win_rate);

figure
imagesc(bin/10, -4:1:2, mean_win_rate); hold on
yticklabels({'2','1','0','-1','-2','-3','-4'});
xticklabels({'-4','-3','-2','-1','0','1','2'});

Yopt = Mean_subdata_opt(:, 2:end) - Mean_subdata_bl_opt;
Xopt = Mean_comdata_opt(:, 2:end) - Mean_subdata_bl_opt;
Yave = Mean_subdata_ave(:, 2:end) - Mean_subdata_bl_ave;
Xave = Mean_comdata_ave(:, 2:end) - Mean_subdata_bl_ave;

ms = 6; min = 20; max = -40;
Yopt = min - abs(-Yopt + max);
Yave = min - abs(-Yave + max);
plot(Xopt/10, Yopt/10, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', ms);
plot(Xave/10, Yave/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);
xlim([-4 2]); ylim([-4 2]);
lineplot(0, 'v', 'k--'); lineplot(min/10 - abs(-0 + max/10), 'h', 'k--');
ylabel('Subject relative aim point [cm]'); xlabel('Opponent relative aim point [cm]');         

h = colorbar; 
ylabel(h,'Chance of winning'); 
axis('square'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 11, 'linewidth', 1);

pos(3) = 900; pos(4) = 350;
set(gcf, 'Position', pos);


