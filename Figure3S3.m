% Kohl C, Wong MXM, Rushworth MFS & Chau BKH: Intraparietal stimulation
% disrupts negative distractor effects in hyman multi-alternative
% decision-making

%% Behavioural Analysis: GLM3
% Script to generate results of Figure3S3

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
    load(strcat('Data',filesep, partic_str))  
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
    rt = [data.MIP(:, 16); data.MT(:, 16)];
    rt = rt/1000;
   
    % exclude trials in which the distractor/empty quadrant was chosen
    rmv = (isnan(accuracy));   
    session(rmv) = [];
    tms(rmv) = [];        
    accuracy(rmv) = [];                
    d(rmv) = [];        
    lv(rmv) = [];        
    hv(rmv) = [];     
    d_loc_binary(rmv) = [];
    rt(rmv) = [];
       
    %% GLM
    for iCondition = 1:length(Conditions)
        % select all trials in the condition of interest 
        idx = (accuracy==1 & session==CondSession(iCondition) & ...
                tms==CondTMS(iCondition) & ...
                d_loc_binary==CondDLoc(iCondition));
        % define regressors and criterion
        regressors = eval(model);
        regressors = normalise(regressors);
        criterion = rt(idx);
        % fit GLM
        beta(iPartic, :, iCondition) = glmfit(regressors,...
                                             criterion);
    end
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
%     if rfig <4, ylim([-0.6 1.2]); else, ylim([-0.3 0.6]); end
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

