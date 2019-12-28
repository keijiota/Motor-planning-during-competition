clear all; close all;

load('data_exp1a.mat') ;
load('data_exp1b.mat') ;

bl_w = 1:50; comp_w = 51:170; 
X = []; Y = [];
for condition = 1:2
    % grab data
    Mean_subdata_indv = []; Mean_subdata = []; Mean_comdata = [];
    if condition == 1
        subdata_indv =  endpoint_opt(:, bl_w); subdata =  endpoint_opt(:, comp_w); comdata = endpoint_com_opt(:, comp_w);
    elseif condition == 2
        subdata_indv =  endpoint_ave(:, bl_w); subdata =  endpoint_ave(:, comp_w); comdata = endpoint_com_ave(:, comp_w);
    end
    [N T] = size(subdata) ; B = T/10;
    BB(condition) = B; ; NN(condition) = N;
    
    % average endpoint in each block
    for ib = 1:B
        if ib <= 5
            Mean_subdata_indv(:, ib) = mean(subdata_indv(:, ib*10-9:ib*10)')' ;
        end
        Mean_subdata(:, ib) = mean(subdata(:, ib*10-9:ib*10)')' ;
        Mean_comdata(:, ib) = mean(comdata(:, ib*10-9:ib*10)')' ;
    end
    Mean_subdata_bl = mean(Mean_subdata_indv')' ;
    
    if condition == 1
        Mean_subdata_bl_opt  = Mean_subdata_bl ;
        Mean_subdata_opt = Mean_subdata;
        Mean_comdata_opt = Mean_comdata;
    elseif condition == 2
        Mean_subdata_bl_ave  = Mean_subdata_bl ;
        Mean_subdata_ave = Mean_subdata;
        Mean_comdata_ave = Mean_comdata;
    end
end

Xo = []; Yo = []; Xa = []; Ya = [];
for condition = 1:2
    for in = 1:NN(condition)
        for ib = 1:BB(condition)
            if condition == 1
                if ib >= 2 % excludes first block
                    Xo = [Xo; Mean_comdata_opt(in,ib) - Mean_subdata_bl_opt(in)]; % Opponent's aim point - baseline level
                    Yo = [Yo; Mean_subdata_opt(in,ib) - Mean_subdata_bl_opt(in)]; % Subject's aim point - baseline level
                end
            elseif condition == 2
                if ib >= 2
                    Xa = [Xa; Mean_comdata_ave(in,ib) - Mean_subdata_bl_ave(in)]; % Opponent's aim point - baseline level
                    Ya = [Ya; Mean_subdata_ave(in,ib) - Mean_subdata_bl_ave(in)]; % Subject's aim point - baseline level
                end
            end
        end
    end
end
Xoa = [Xo; Xa]; Yoa = [Yo; Ya];

% plot data
ms = 6;
figure(1)
subplot(1,2,1)
plot(Xo/10, Yo/10, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', ms); hold on
plot(Xa/10, Ya/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms);


%% bootstrap and permutation test
originaldataX_left = []; originaldataY_left = [];
originaldataX_right = []; originaldataY_right = [];

bw = 0.05;

A = []; Estimates_left = []; Estimates_right = [];
originaldataX = Xoa; originaldataY = Yoa; thres = mean(Xoa);
[N S] = size(originaldataX);

% dividing data to left half and right half
for i = 1:N
    if originaldataX(i) < thres
        originaldataX_left = [originaldataX_left; originaldataX(i)]; originaldataY_left = [originaldataY_left; originaldataY(i)];
    else
        originaldataX_right = [originaldataX_right; originaldataX(i)]; originaldataY_right = [originaldataY_right; originaldataY(i)];
    end
end

% calculaitng regression lines of original data
estimates_orgleft = polyfit(originaldataX_left, originaldataY_left,1);
estimates_orgright = polyfit(originaldataX_right, originaldataY_right,1);
Estimates_orgleft(1,:) = estimates_orgleft;
Estimates_orgright(1,:) = estimates_orgright;

% take the difference of original data
diff_of_originaldata = estimates_orgleft(1) - estimates_orgright(1);

% bootstraping
bootN = 50000; % times
bootN = 1000; 

for nrep = 1:bootN % repeating N times
    resampledataX = []; resampledataY = [];
    resampledataX_left = []; resampledataY_left = [];
    resampledataX_right = []; resampledataY_right = [];
    
    % produce random number
    indx = fix(rand(N,1) * N + 1);
    % resampling data
    resampledataX = originaldataX(indx) ;
    resampledataY = originaldataY(indx) ;
    
    % dividing resampled data to left half and right half
    for i = 1:N
        if resampledataX(i) < thres
            resampledataX_left = [resampledataX_left; resampledataX(i)]; resampledataY_left = [resampledataY_left; resampledataY(i)];
        else
            resampledataX_right = [resampledataX_right; resampledataX(i)]; resampledataY_right = [resampledataY_right; resampledataY(i)];
        end
    end
    
    % calculaitng regression lines of resampled data
    estimates_left = polyfit(resampledataX_left, resampledataY_left,1);
    estimates_right = polyfit(resampledataX_right, resampledataY_right,1);
    
    % storing coefficient
    Estimates_left(nrep,:) = estimates_left;
    Estimates_right(nrep,:) = estimates_right;
    
    % take the difference of resampled data
    diff_of_randomized_data = estimates_left(1) - estimates_right(1);
    
    % check that this difference is difference from 0
    % (assumption is that there is no difference between the coefficiences of regression lines)
    if diff_of_originaldata < 0
        if diff_of_randomized_data > 0
            A(nrep) = 1;
        else
            A(nrep) = 0;
        end
    elseif diff_of_originaldata > 0
        if diff_of_randomized_data < 0
            A(nrep) = 1;
        else
            A(nrep) = 0;
        end
    end
    Diff_of_randomized_data(nrep) = diff_of_randomized_data;
end

% count how many times this mean difference is larger than the origina mean difference
sumofdiff = sum(A);

% calculate P value for permutation test
Pvalue = sumofdiff/bootN ; 

% average and confidence interval for histograms
mean(Estimates_left(:,1)) ; 
sEstimates_left = sort(Estimates_left(:,1));
CI95_left = [sEstimates_left(bootN*0.025), sEstimates_left(bootN*0.975)]; 
mean(Estimates_right(:,1)); 
sEstimates_right = sort(Estimates_right(:,1));
CI95_right = [sEstimates_right(bootN*0.025), sEstimates_right(bootN*0.975)] ; 

% plot bootstraping samples
subplot(1,2,2)
histogram(Estimates_left(:,1), 'FaceAlpha', 0.5, 'FaceColor','b','Normalization', 'probability', 'BinWidth', bw); hold on
histogram(Estimates_right(:,1), 'FaceAlpha', 0.5, 'FaceColor','r','Normalization', 'probability', 'BinWidth', bw);
lineplot(mean(Estimates_left(:,1)), 'v', 'b', 'linewidth', 2);
lineplot(mean(Estimates_right(:,1)), 'v', 'r', 'linewidth', 2);
xlim([-0.2 1.2]); xticks(-0.2:0.2:1.2); ylim([0 0.3]); yticks(0:0.1:0.3);
ylabel('Normalized histogram'); xlabel('Slope');         
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 11, 'linewidth', 1);


% plot regression lines and a quadratic fit
subplot(1,2,1)
thres = mean(Xoa);
yhut = polyval([Estimates_orgleft(1,1), Estimates_orgleft(1,2)/10], -4:0.01:thres/10);
plot(-4:0.01:thres/10, yhut, 'b-', 'linewidth', 2);
yhut = polyval([Estimates_orgright(1,1), Estimates_orgright(1,2)/10], thres/10:0.01:2);
plot(thres/10:0.01:2, yhut, 'r-', 'linewidth', 2); 

x = -5:0.001:2;
mdl = fitglm(Xoa/10,Yoa/10,'y ~ 1 + x1 + x1^2')  ; 
yhut_q = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2).*x + mdl.Coefficients.Estimate(3).*x.^2 ;
yhut_q = yhut_q' ; 
plot(x,yhut_q, 'g--', 'linewidth', 2); 

xlim([-4 2]); ylim([-4 2]); yticks(-4:1:2); xticks(-4:1:2);
lineplot(0, 'h', 'k--'); lineplot(0, 'v', 'k--');
axis('square'); 
ylabel('Subject relative aim point [cm]'); xlabel('Opponent relative aim point [cm]');         
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 11, 'linewidth', 1);

pos(3) = 900; pos(4) = 350;
set(gcf, 'Position', pos);






