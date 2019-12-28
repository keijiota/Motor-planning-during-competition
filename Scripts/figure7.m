clear all; close all;

load('data_exp4.mat') ;

%% average risk-sensitivity

% parameters for analytical solution of optimal aim point
a = 70; % gain start point(mm)
b = 300; % gain end point(mm)
c = 100; % maximum reward
alpha = c/(b-a); % slope
beta = -c*a/(b-a); %
mu_seq = linspace(220,300,200); % aimpoint

% grab data
bl_w = 1:50; comp_w = 51:170;
Std_subdata_indv = []; Std_subdata = [];
Mean_subdata_indv = []; Mean_subdata = [];
OptAim_indv = []; OptAim = [];
RS_indv = []; RS = []; 

subdata_indv = endpoint_highlyave(:, bl_w); subdata =  endpoint_highlyave(:, comp_w); comdata = endpoint_com_highlyave(:, comp_w);
[N T] = size(subdata) ; B = T/10;

Std_subdata_indv = std(subdata_indv(:, :)')' ;
Mean_subdata_indv = mean(subdata_indv(:, :)')' ;

% analytical solution of optimal aim point
for i = 1:B
    Std_subdata(:, i) = std(subdata(:, i*10-9:i*10)')' ;
    Mean_subdata(:, i) = mean(subdata(:, i*10-9:i*10)')' ;
end

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

% plot data
figure
d = 0.5; ms = 8;
subplot(1,12,1:4)
errorbar(1.5, mean(RS_indv/10), std(RS_indv/10)/sqrt(N), 'co', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
seshade(RS(:,1:B)/10, d, 'c', 'c-', 3:B+2) ;
lineplot(0, 'h', 'k--');
xlim([0 15]); ylim([-1 2]); yticks(-1:0.5:2);
xticks([1.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]); xticklabels(char('BL', '1','','','4','','','','8','','','','12'));
ylabel('Risk-sensitivity [cm]'); xlabel('Blocks');         
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);



%% average endpoint 
% grab data
Mean_subdata_indv = []; MMean_subdata_indv = []; Mean_subdata = []; Mean_comdata = [];

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
subplot(1,12,5.5:7.5)
errorbar(1.5, mean(MMean_subdata_indv)/10, std(MMean_subdata_indv/10)/sqrt(N), 'bo', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', ms) ; hold on
seshade(Mean_subdata(:,1:B)/10, d, 'b', 'b-', 3:B+2) ;
seshade(Mean_comdata(:,1:B)/10, d, 'r', 'r-', 3:B+2) ;
xlim([0 15]); ylim([24 29]); yticks(24:1:29); xticks([1.5, 3, 6, 10, 14]); xticklabels(char('BL', '1', '4', '8', '12'));
ylabel('Aim point [cm]'); xlabel('Blocks');         
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);


%% subject's aim point vs opponent's aim point 

X = []; Y = [];
for in = 1:N
    for ib = 1:B
        if ib >= 2 % excludes the first block
            X = [X; Mean_comdata(in,ib) - MMean_subdata_indv(in)]; % Opponent's aim point - baseline level
            Y = [Y; Mean_subdata(in,ib) - MMean_subdata_indv(in)]; % Subject's aim point - baseline level
        end
    end
end

% plot data
ms = 6;
subplot(1,12,9:12)
plot(X/10, Y/10, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', ms); hold on
xlim([-4 2]); ylim([-4 2]); yticks(-4:1:2); xticks(-4:1:2);
plot(-4:2, -4:2, 'k--');
lineplot(0, 'h', 'k--'); lineplot(0, 'v', 'k--');
axis('square'); 
ylabel('Subject relative aim point [cm]'); xlabel('Opponent relative aim point [cm]');         
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);

pos(3) = 1250; pos(4) = 350;
set(gcf, 'Position', pos);



