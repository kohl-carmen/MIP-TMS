% Kohl C, Wong MXM, Rushworth MFS & Chau BKH: Intraparietal stimulation
% disrupts negative distractor effects in hyman multi-alternative
% decision-making

%% Behavioural Analysis: GLM3
% Script to generate results of Figure3a-c
% Written by Carmen Kohl, 2020.
% github.com/kohl-carmen/MIP-TMS

% Applies GLM3 to predict accuracy in each condition 
%       GLM3:	β0 + β1 z(HV-LV) + β2 z(HV+LV) + β3 z(D-HV) + ε
%       Conditions: MIP/MT x TMS/NonTMS x Ipsilateral/Contralateral D
% Reports Session(MIP/MT) x Stimulation (TMS/NonTMS) x Distractor Location 
% (Ipsilateral/Contralateral)ANOVAs for each GLM3 regressor, reports ttests
% of each regressor's beta coefficients & plots Figure 3a-c  

clearvars 
% set directory
dir = fileparts(which('GLM1.m'));
cd(dir)

nPartic = 31; % nr of participants
% define conditions 
Conditions = {'MIP0ipsi', 'MIP1ipsi', 'MT0ipsi', 'MT1ipsi',...
              'MIP0contra', 'MIP1contra', 'MT0contra', 'MT1contra'};
CondSession = [1, 1, 0, 0, 1, 1, 0, 0];
CondTMS = [0, 1, 0, 1, 0, 1, 0, 1];
CondDLoc = [0, 0, 0, 0, 1, 1, 1, 1];
% GLM3
regressor_str = {'HV+LV' 'HV-LV' 'D-HV' };
model = '[hv(idx) + lv(idx), hv(idx) - lv(idx), d(idx) - hv(idx)]';

beta = nan(nPartic,length(regressor_str) + 1,length(Conditions));
for iPartic = 1:nPartic
    %% Prepare data
    partic_str = sprintf('%02d', iPartic);
    load(strcat('Data\',partic_str))  
    % select variables of interest (see data.Key)
    session = [ones(size(data.MIP, 1),1); zeros(size(data.MT, 1),1)];
    tms = [data.MIP(:,14); data.MT(:,14)]; % 1=TMS, 0=NonTMS 
    d = [data.MIP(:, 4); data.MT(:, 4)]; % distractor value
    lv = [data.MIP(:, 3); data.MT(:, 3)]; % low value
    hv = [data.MIP(:, 2); data.MT(:, 2)]; % high value
    d_loc = [data.MIP(:, 13); data.MT(:, 13)]; %distractor location
    d_loc_binary = zeros(size(d_loc));
    d_loc_binary(d_loc==2 | d_loc==4) = 1; % 0=contralateral, 1=ipsilateral
    accuracy = [data.MIP(:, 18); data.MT(:, 18)]; % 1=high value chosen, 
                                                  % 0=low value chosen, 
                                                  % nan=distractor/empty 
                                                  % quadrant chosen
   
    % exclude trials in which the distractor/empty quadrant was chosen
    rmv = (isnan(accuracy));   
    session(rmv) = [];
    tms(rmv) = [];        
    accuracy(rmv) = [];                
    d(rmv) = [];        
    lv(rmv) = [];        
    hv(rmv) = [];     
    d_loc_binary(rmv) = [];
       
    %% GLM
    for iCondition = 1:length(Conditions)
        % select all trials in the condition of interest 
        idx = (session==CondSession(iCondition) & ...
               tms==CondTMS(iCondition) & ...
               d_loc_binary==CondDLoc(iCondition));
        % define regressors and criterion
        regressors = eval(model);
        regressors = normalise(regressors);
        criterion = accuracy(idx);
        % fit GLM
        beta(iPartic, :, iCondition) = glmfit(regressors,...
                                             criterion, 'binomial');
    end
end


%% ANOVA
% Site x Stimulation x Distractor Location ANOVA for each regressor
fprintf('\nANOVA (GLM3):')
iEffect = 3:2:15;
% set up ANOVA
Effect = {'Session', 'Stim', 'DLoc', 'Session*Stim', 'Session*DLoc*', ...
        'Stim*DLoc', 'Session*Stim*DLoc'};
within1 = categorical(CondSession)';
within2 = categorical(CondTMS)';
within3 = categorical(CondDLoc)';
within = table(within1, within2, within3, 'variablenames',...
               {'Session', 'Stim', 'DLoc'});
factors = 'Session*Stim*DLoc';
allvars = strcat(Conditions{1}, '-', Conditions{end}, ' ~1');
for regressor = 2:length(regressor_str) + 1
    % ANOVA
    tab = array2table(squeeze(beta(:, regressor,:)),...
                    'variablenames', Conditions);   
    rm = fitrm(tab, allvars, 'WithinDesign', within);
    ranovatbl = ranova(rm, 'withinmodel', factors);   
    % print output (F-statistic, p-value, partial eta squared)
    fprintf('\n- - %s - -\n', regressor_str{regressor-1})
    ranovatbl = table2array(ranovatbl);
    for i = 1:length(iEffect)
         fprintf('%s : F(%d, %d) = %4.2f, p = %1.3f, n2 = %4.3f \n',...
                 Effect{i},  ranovatbl(1, 2), ranovatbl(2, 2), ...
                 ranovatbl (iEffect(i), 4), ranovatbl(iEffect(i), 5), ...
                 ranovatbl(iEffect(i), 1) / (ranovatbl(iEffect(i), 1) + ...
                 ranovatbl(iEffect(i) + 1, 1))) 
    end
end         

%% Ttests
% one-sample t-test for each regressor
fprintf('\nTtest (GLM3):')
mean_betas = squeeze(mean(beta, 3));
for i = 2:size(mean_betas,2)
    % ttest
    [~, p, ~, stats] = ttest(mean_betas(:, i));
    cohensd = (mean(mean_betas(:, i))) / std(mean_betas(:, i));
    % print output (t-statistic, p-value, Cohen's D)
    fprintf('\n%s', regressor_str{i-1})
    fprintf('\nt(%i) = %2.3f, p = %2.3f, d = %2.3f\n', stats.df, ...
            stats.tstat, p, cohensd)
end
    
%% Plot Figure 3 a-c
xax = [1, 2; 3, 4];
clr.s1 = {[ 1, 0.75, 0.5],[0.875, 0.438, 0]};
clr.s0 = {[ 0.5, 0.5, 0.5], 'k'};
linetype = {'--', '-'};
for rfig = 2: length(regressor_str) + 1
    figure(rfig-1)
    hold on
    count = 0;
    leg = nan(1, length(unique(session))*length(unique(tms)));
    for sessioni = 0:1
        for tmsi = 0:1
            count = count + 1;
            conds_oi = find(CondSession==sessioni & CondTMS==tmsi);
            to_plot = squeeze(mean(beta(:, rfig, conds_oi)));
            leg(count) = plot(xax(sessioni + 1,:), to_plot, ...
                         linetype{tmsi + 1}, 'Linewidth', 2, 'Color', ...
                         clr.(strcat('s', num2str(sessioni))){tmsi + 1});
            %standard error
            x = beta(:, rfig, conds_oi(1));
            errors(1) = std(x) / sqrt(length(x));       
            x = beta(:, rfig, conds_oi(2));
            errors(2) = std(x) / sqrt(length(x));          
            errorbar(xax(sessioni + 1, :), to_plot, errors, ...
                     'Linewidth',2, 'Linestyle', 'none', 'Color', ...
                     clr.(strcat('s',num2str(sessioni))){tmsi + 1});
        end
    end
    legend(leg,{'MT0', 'MT1', 'MIP0', 'MIP1'})
    xlim([0 5])
    set(gca,'XTick',1:4)
    set(gca,'xticklabel', {'Ipsi' 'Contra' 'Ipsi' 'Contra'})
    title(regressor_str{rfig-1})        
    if rfig <4, ylim([-0.6 1.2]); else, ylim([-0.3 0.6]); end
end
     


function x = normalise(x)
    % removes mean from data x
    dim = 1;
    dims = size(x);
    dimsize = size(x,dim);
    dimrep = ones(1, length(dims));
    dimrep(dim) = dimsize;

    x = (x - repmat(nanmean(x, dim), dimrep)) ./ ...
         repmat(nanstd(x, 0, dim), dimrep);
end

