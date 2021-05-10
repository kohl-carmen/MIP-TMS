% Kohl C, Wong MXM, Rushworth MFS & Chau BKH: Intraparietal stimulation
% disrupts negative distractor effects in hyman multi-alternative
% decision-making

%% Behavioural Analysis: GLM2
% Script to generate results of Figure 2d
% Written by Carmen Kohl, 2020.
% github.com/kohl-carmen/MIP-TMS

% Applies GLM2 to predict accuracy in non-tms trials 
%       GLM2:	Step 1, β0 + β1 z(HV-LV) + β2 z(HV+LV) + ε1
%               Step 2, β3 + β4 z(D-HV) + ε2
% Reports ttests of each predictor's beta weights & plots Figure 2d 

clearvars 
% set directory
dir = fileparts(which('GLM1.m'));
cd(dir)
 
% GLM2
nPartic = 31; % nr of participants
modelStep1 = '[hv(idx) + lv(idx), hv(idx) - lv(idx)]';
regressor_str_Step2 = {'D-HV'};
modelStep2 = '[d(idx) - hv(idx)]'; 

for iPartic = 1:nPartic
    %% Prepare data
    partic_str = sprintf('%02d', iPartic);
    load(strcat('Data',filesep,partic_str))   
    % select variables of interest (see data.Key)
    tms = [data.MIP(:, 14); data.MT(:, 14)];  % 1=TMS, 0=NonTMS
    d = [data.MIP(:, 4); data.MT(:, 4)]; % distractor value
    lv = [data.MIP(:, 3); data.MT(:, 3)]; % low value
    hv = [data.MIP(:, 2); data.MT(:, 2)];% high value
    accuracy = [data.MIP(:, 18); data.MT(:, 18)]; % 1=high value chosen, 
                                                  % 0=low value chosen, 
                                                  % nan=distractor/empty 
                                                  % quadrant chosen
   
    % exclude trials in which the distractor/empty quadrant was chosen
    rmv =(isnan(accuracy));   
    tms(rmv) = [];        
    accuracy(rmv) = [];                
    d(rmv) = [];        
    lv(rmv) = [];        
    hv(rmv ) = [];           
           
    % median split(HV-LV)
    mediansplit = zeros(size(hv));
    mediansplit(hv-lv >= median(hv-lv)) = 1;
    split = {'low','high'};
    
    for median_split=[0 1]
        %% GLM
        % select all trials in the condition of interest 
        idx = (tms==0 & mediansplit==median_split);
        % Step1: define regressors and criterion
        regressorsStep1 = eval(modelStep1);
        regressorsStep1 = normalise(regressorsStep1);
        criterionStep1 = accuracy(idx);
        % Step1: fit GLM
        [~,~,statsStep1] = glmfit(regressorsStep1,...
                                  criterionStep1,'binomial');
        % Step2: define regressors and criterion 
        regressorsStep2 = eval(modelStep2);
        regressorsStep2 = normalise(regressorsStep2);
        criterionStep2 = statsStep1.resid;
        % Step2: fit GLM
        betatmp = glmfit(regressorsStep2,criterionStep2);
        [betasStep2.(split{median_split + 1})(iPartic,:)] = betatmp;
                                    
    end
end


%% Ttest
% one-sample t-test for D-HV (for high/low HV-LV)
fprintf('\nTtest (GLM2 - NonTMS): %s\n',regressor_str_Step2{:})
% ttest (easy)
fprintf('\nD-HV (high HV-LV): ')
[~, p, ~, stats] = ttest(betasStep2.high(:,end));
cohensD = (mean(betasStep2.high(:,end))-0) / std(betasStep2.high(:,end));
% print output (t-statistic, p-value, Cohen's D)
fprintf('t(%d) = %2.2f, p = %2.3f, d = %2.2f\n', stats.df, ...
        stats.tstat, p, cohensD)
    
% ttest (hard)
fprintf('\nD-HV (low HV-LV): ')
[~, p, ~, stats] = ttest(betasStep2.low(:,end));
cohensD = (mean(betasStep2.low(:,end))-0) / std(betasStep2.low(:,end));
% print output (t-statistic, p-value, Cohen's D)
fprintf('t(%d) = %2.2f, p = %2.3f, d = %2.2f\n', stats.df, ...
        stats.tstat, p, cohensD)


%% Plot Figure 2 d
figure
bar([mean(betasStep2.low(:,end));mean(betasStep2.high(:,end))]);
hold on
% standard error
x = betasStep2.high(:, end);
errors_high = std(x) / sqrt(length(x));
x = betasStep2.low(:, end);
errors_low = std(x) / sqrt(length(x));    
err = errorbar(1:2, [mean(betasStep2.low(:, end)); ...
             mean(betasStep2.high(:, end))], [errors_low, errors_high]);
err.Color = [0 0 0];                            
err.LineStyle = 'none';  
set(gca, 'xticklabel', {'hard', 'easy'})
title('GLM2 - Figure 2d')
ylabel('Effect size on accuracy')


function x = normalise(x)
    % removes mean from data x
    dim = 1;
    dims = size(x);
    dimsize = size(x, dim);
    dimrep = ones(1, length(dims));
    dimrep(dim) = dimsize;

    x = (x - repmat(nanmean(x, dim), dimrep)) ./ ...
         repmat(nanstd(x, 0, dim), dimrep);
end
