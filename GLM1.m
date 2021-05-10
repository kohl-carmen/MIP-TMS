% Kohl C, Wong MXM, Rushworth MFS & Chau BKH: Intraparietal stimulation
% disrupts negative distractor effects in hyman multi-alternative
% decision-making

%% Behavioural Analysis: GLM1
% Script to generate results of Figure 2c
% Written by Carmen Kohl, 2020.
% github.com/kohl-carmen/MIP-TMS

% Applies GLM1 to predict accuracy in non-tms trials 
%       GLM1:	β0 + β1 z(HV-LV) + β2 z(HV+LV) + β3 z(D-HV) 
%               + β4 z(HV-LV) z(D-HV) + ε
% Reports ttests of each predictor's beta weights & plots Figure 2c 

clearvars 
% set directory
dir = fileparts(which('GLM1.m'));
cd(dir)
 
% GLM1
nPartic = 31; % nr of participants
regressor_str = {'HV-LV', 'HV+LV', 'D-HV', '(HV-LV)*(D-HV)'};
model = ['[hv(idx) - lv(idx), hv(idx) + lv(idx), d(idx) - hv(idx),'...
        'normalise(hv(idx) - lv(idx)).*normalise(d(idx) - hv(idx))]'];

beta = nan(nPartic, length(regressor_str) + 1);
for iPartic = 1:nPartic
    %% Prepare data
    partic_str = sprintf('%02d', iPartic);
    load(strcat('Data',filesep, partic_str))
    
    % select variables of interest (see data.Key)
    tms = [data.MIP(:, 14); data.MT(:, 14)]; % 1=TMS, 0=NonTMS
    d = [data.MIP(:, 4); data.MT(:, 4)]; % distractor value
    lv = [data.MIP(:, 3); data.MT(:, 3)]; % low value 
    hv = [data.MIP(:, 2); data.MT(:, 2)]; % high value
    accuracy = [data.MIP(:, 18); data.MT(:, 18)]; % 1=high value chosen, 
                                                  % 0=low value chosen, 
                                                  % nan=distractor/empty 
                                                  % quadrant chosen

   
    % exclude trials in which the distractor/empty quadrant was chosen
    rmv = (isnan(accuracy));   
    tms(rmv) = [];        
    accuracy(rmv) = [];                
    d(rmv) = [];        
    lv(rmv) = [];        
    hv(rmv ) = [];        
   
    %% GLM
    % select all trials in the condition of interest 
    idx = (tms==0);
    % define regressors and criterion
    regressors = eval(model);
    regressors = normalise(regressors);
    criterion = accuracy(idx);
    % fit GLM
    [beta(iPartic, :)] = glmfit(regressors,criterion,'binomial');    
end

%% Ttest
% one-sample t-test for each regressor
for iRegs = 2:length(regressor_str) + 1
    % ttest
    [~, p, ~, stats] = ttest(beta(:, iRegs));
    cohensD = (mean(beta(:, iRegs))-0) / std(beta(:, iRegs));
    % print output (t-statistic, p-value, Cohen's D)
    fprintf('\nTtest (GLM1 - NonTMS): %s \n', regressor_str{iRegs-1})
    fprintf('t(%d) = %2.2f, p = %2.3f, d = %2.2f\n', stats.df,...
            stats.tstat, p, cohensD) 
end


%% Plot Figure 2 c
figure
bar(mean(beta(:, 2:end)));
hold on
% standard error
errors = nan(1, size(beta, 2));
for iRegs = 1:length(regressor_str) + 1
    x = (beta(:, iRegs ));
    errors(iRegs) = std(x) / sqrt(length(x));
end
err = errorbar(1:4, mean(beta(:, 2:end)), errors(2:end));
err.Color = [0 0 0];                            
err.LineStyle = 'none';  
set(gca, 'xticklabel', regressor_str)
title('GLM1 - Figure 2c')
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