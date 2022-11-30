%% Eyetracking Figures
clear
% runs GLM to see the impact of different predictors on different gaze

% Figure 5 a
% Which_GLM = 1;
% Which criterion = 1;
%
% Figure 5 b
% Which_GLM = 1;
% Which criterion = 2;
%
% Figure 5 c
% Which_GLM = 1;
% Which criterion = 3;
%
% Figure 5 d
% Which_GLM = 3;
%
% Figure 5 S1 a
% Which_GLM = 2;
% Which criterion = 1;
%
% Figure 5 S1 b
% Which_GLM = 2;
% Which criterion = 2;
%
% Figure 5 S1 c
% Which_GLM = 2;
% Which criterion = 3;

Which_GLM = 2;
% 1: HV+LV, HV-LV, D-HV  => these values to predict the number of gaze shifts in different conds (TMSxSess)
% 2: HV, LV, D  => these values to predict the number of gaze shifts in different conds (TMSxSess)
% 3: hv_lv, lv_hv, lv_dv, dv_lv, hv_dv, dv_hv  => these gaze shifts to predict accuracy across all non-tms trials

Which_criterion = 3;
% If GLM is 1 or 2 (predicting gaze shifts), which gaze shifts
% 1: bidirectional D to HV and HV to D
% 2: D to HV
% 3: HV to D

if Which_GLM == 1
    Regs = {'HV+LV' 'HV-LV' 'D-HV'};
    model = '[hv_v(idx)+lv_v(idx), hv_v(idx)-lv_v(idx),d_v(idx)-hv_v(idx)]';
    Testing = 2;
    ANOVA = 1;
elseif Which_GLM == 2
    Regs = {'HV' 'LV' 'DV'};
    model = '[hv_v(idx), lv_v(idx), d_v(idx)]';
    Testing = 2;
    ANOVA = 1;
elseif Which_GLM == 3
    Regs = {'hv_lv' 'lv_hv' 'lv_dv' 'dv_lv' 'hv_dv' 'dv_hv'};
    model = '[hv_lv(idx),lv_hv(idx),lv_dv(idx), dv_lv(idx), hv_dv(idx), dv_hv(idx)]';
    Testing = 1;
    ANOVA = 10;
end

if Which_criterion == 1
    criterion_name = ['bi_dv_hv'];
elseif Which_criterion == 2
    criterion_name = ['dv_hv'];
elseif Which_criterion == 3
    criterion_name = ['hv_dv'];
end

Partic = 1:31;

Sessions = {'MT' 'MIP'};

if ANOVA == 1
    Conds = {'MIP0' 'MIP1' 'MT0' 'MT1'};
elseif ANOVA == 2
    Conds = {'MIP0ipsi' 'MIP1ipsi' 'MT0ipsi' 'MT1ipsi' 'MIP0contra' 'MIP1contra' 'MT0contra' 'MT1contra'};
elseif ANOVA == 10
    Conds = {'nontms'};
end

betas = [];

for partic = 1:length(Partic)
    nr_fix = [];
    dur_fix = [];
    rel_dur_fix = [];
    bad_trials = [];

    for sess = 1:length(Sessions)
        subj_code = 7000 + Partic(partic) * 10 + sess - 1;
        % load fix structure
        load(strcat('Data',filesep, 'Tobii',filesep, num2str(subj_code), filesep, 'FIXATIONS')) 

        temp = length(nr_fix);

        for i = 1:length(FIXATIONS.nr_fix)
            nr_fix(temp + i, :) = FIXATIONS.nr_fix{i};
            dur_fix(temp + i, :) = FIXATIONS.total_dur{i};
            rel_dur_fix(temp + i, :) = FIXATIONS.total_dur{i} ./ sum(FIXATIONS.total_dur{i});

            %gaze
            gaze{temp + i} = [];

            if sum(FIXATIONS.nr_fix{i}) > 0
                cont = FIXATIONS.timeseries{i};
                latest_fix = 0;

                for j = 2:length(cont)

                    for k = 1:4

                        if cont(j, k) == 1 & (cont(j - 1, k) == 0 | j == 2)
                            gaze{temp + i} = [gaze{temp + i}; latest_fix, k];
                            latest_fix = k;
                        end

                    end

                end

            end

        end

    end

    %match trials to beh data
    load(strcat('Data',filesep, sprintf('%02d', partic))) 
    sess_i = [ones(size(data.MT, 1), 1); ones(size(data.MIP, 1), 1) .* 2];
    tms = [data.MT(:, 14); data.MIP(:, 14)]; % 1=TMS, 0=NonTMS
    d_v = [data.MT(:, 4); data.MIP(:, 4)]; % distractor value
    lv_v = [data.MT(:, 3); data.MIP(:, 3)]; % low value
    hv_v = [data.MT(:, 2); data.MIP(:, 2)]; % high value
    d_loc = [data.MT(:, 13); data.MIP(:, 13)]; %distractor location
    dloc = zeros(size(d_loc));
    dloc(d_loc == 2 | d_loc == 4) = 1; % 0=contralateral, 1=ipsilateral
    acc = [data.MT(:, 18); data.MIP(:, 18)]; % 1=high value chosen,
    % 0=low value chosen,
    % nan=distractor/empty
    % quadrant chosen
    rt = [data.MT(:, 16); data.MIP(:, 16)];

    hv_lv = zeros(540, 1);
    hv_dv = zeros(540, 1);
    hv_all = zeros(540, 1);

    lv_hv = zeros(540, 1);
    lv_dv = zeros(540, 1);
    lv_all = zeros(540, 1);

    dv_hv = zeros(540, 1);
    dv_lv = zeros(540, 1);
    dv_all = zeros(540, 1);

    all_dv = zeros(540, 1);
    all_hv = zeros(540, 1);
    all_lv = zeros(540, 1);

    %bidirectional
    every_hv = zeros(540, 1);
    every_lv = zeros(540, 1);
    every_dv = zeros(540, 1);

    bi_dv_hv = zeros(540, 1);
    bi_dv_lv = zeros(540, 1);
    bi_hv_lv = zeros(540, 1);

    for trial = 1:540

        if ~isempty(gaze{trial})

            for i = 2:size(gaze{trial}, 1)

                if gaze{trial}(i, 1) == 1
                    hv_all(trial) = hv_all(trial) + 1;
                    every_hv(trial) = every_hv(trial) + 1;

                    if gaze{trial}(i, 2) == 2
                        hv_lv(trial) = hv_lv(trial) + 1;
                        bi_hv_lv(trial) = bi_hv_lv(trial) + 1;
                    elseif gaze{trial}(i, 2) == 3
                        hv_dv(trial) = hv_dv(trial) + 1;
                        bi_dv_hv(trial) = bi_dv_hv(trial) + 1;
                    end

                elseif gaze{trial}(i, 1) == 2
                    lv_all(trial) = lv_all(trial) + 1;
                    every_lv(trial) = every_lv(trial) + 1;

                    if gaze{trial}(i, 2) == 1
                        lv_hv(trial) = lv_hv(trial) + 1;
                        bi_hv_lv(trial) = bi_hv_lv(trial) + 1;
                    elseif gaze{trial}(i, 2) == 3
                        lv_dv(trial) = lv_dv(trial) + 1;
                        bi_dv_lv(trial) = bi_dv_lv(trial) + 1;
                    end

                elseif gaze{trial}(i, 1) == 3
                    dv_all(trial) = dv_all(trial) + 1;
                    every_dv(trial) = every_dv(trial) + 1;

                    if gaze{trial}(i, 2) == 1
                        dv_hv(trial) = dv_hv(trial) + 1;
                        bi_dv_hv(trial) = bi_dv_hv(trial) + 1;
                    elseif gaze{trial}(i, 2) == 2
                        dv_lv(trial) = dv_lv(trial) + 1;
                        bi_dv_lv(trial) = bi_dv_lv(trial) + 1;
                    end

                end

                if gaze{trial}(i, 2) == 1
                    all_hv(trial) = all_hv(trial) + 1;
                    every_hv(trial) = every_hv(trial) + 1;
                elseif gaze{trial}(i, 2) == 3
                    all_dv(trial) = all_dv(trial) + 1;
                    every_dv(trial) = every_dv(trial) + 1;
                elseif gaze{trial}(i, 2) == 2
                    all_lv(trial) = all_lv(trial) + 1;
                    every_lv(trial) = every_lv(trial) + 1;
                end

            end

        end

    end

    hv_e = nr_fix(:, 1);
    lv_e = nr_fix(:, 2);
    d_e = nr_fix(:, 3);
    x_e = nr_fix(:, 4);

    for run = 1:length(Conds)

        if length(Conds) == 1
            idx = tms == 0;
        else

            if any(run == [1 2 5 6]) % MIP
                this_session = 2;
            elseif any(run == [3 4 7 8]) % MT
                this_session = 1;
            end

            if any(run == [1 3 5 7]) % no stim
                this_tms = 0; %no stim
            elseif any(run == [2 4 6 8]) %  stim
                this_tms = 1;
            end

            if any(run == [1 2 3 4])
                d_loc_oi = 0;

            elseif any(run == [5 6 7 8])
                d_loc_oi = 1;

            end

            if ANOVA == 1
                idx = sess_i == this_session & tms == this_tms & ~isnan(acc);
            elseif ANOVA == 2
                idx = sess_i == this_session & tms == this_tms & dloc == d_loc_oi & ~isnan(acc);
            end

        end

        regressors = [eval(model)];

        regressors = [normalise(regressors)];

        if Testing == 1
            criterion = [acc(idx)];

            [betas(partic, :, run), dev, stats] = glmfit(regressors, criterion, 'binomial');
            [warnMsg, warnId] = lastwarn;

        else
            criterion = eval(criterion_name);
            criterion = [criterion(idx)];
            [betas(partic, :, run), dev, stats] = glmfit(regressors, criterion);

        end

    end

end

%% plot betas
clf

for c = 1:length(Conds)

    for r = 1:length(Regs) + 1
        x = (betas(~isnan(betas(:, r, c)), r, c));
        SEM = std(x) / sqrt(length(x)); % Standard Error
        ts = tinv([0.025 0.975], length(x) - 1); % T-Score
        CI = mean(x) + ts * SEM; % Confidence Intervals
        errors.(strcat('Reg', num2str(r))).(Conds{c}) = SEM;
    end

end

to_plot = [];

for r = 1:length(Regs) + 1

    for c = 1:length(Conds)
        to_plot(r, c) = ([nanmean(betas(:, r, c))]);
    end

end

if ANOVA == 10
    to_plot = to_plot(2:end);
end

h = bar(to_plot, .9);
hold on
%errors
for r = 1:length(Regs) + 1

    for c = 1:length(Conds)
        error_to_plot(r, c) = ([errors.(strcat('Reg', num2str(r))).(Conds{c})]);
    end

end

if ANOVA == 10
    error_to_plot = error_to_plot(2:end);
end

[numgroups, numbars] = size(to_plot);
groupwidth = min(0.8, numbars / (numbars + 1.5));

for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange: https://uk.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab-7-13-r2011b
    x = (1:numgroups) - groupwidth / 2 + (2 * i - 1) * groupwidth / (2 * numbars); % Aligning error bar with individual bar
    errorbar(x, to_plot(:, i), error_to_plot(:, i), 'k', 'linestyle', 'none');
end

hold on
set(gca, 'XTick', 1:length(Regs) + 1)
set(gca, 'xticklabel', {'INTC', Regs{:}})
xlabel('Regressors')
ylabel('Beta Value')
legend(Conds, 'Location', 'southwest')

if ANOVA == 10
    ylim([-0.25 0.3])
    set(gca, 'YTick', [- .25 0 .25])

else

    title(strcat('Nr Fixations (mean/trial)'))

end

if ANOVA < 10
    close all

    for rfig = 2:length(Regs) + 1 %one figure per regressor (except intercept)
        figure(rfig - 1)
        to_plot = [];
        error_to_plot = [];

        for MTMIP = 1:2

            if MTMIP == 1
                these_conds = [3 4];
            else
                these_conds = [1 2];
            end

            for cond = 1:length(Conds) / 2 % ipsi0 contra0 ipsi1 contra1
                to_plot(MTMIP, cond) = mean(betas(:, rfig, these_conds(cond)));
                error_to_plot(MTMIP, cond) = ([errors.(strcat('Reg', num2str(rfig))).(Conds{these_conds(cond)})]);

            end

        end

        h = bar(to_plot, .9);
        colours = {[1 .75 .5] [.875 .438 0]};

        for i = 1:length(Conds) / 2
            h(i).FaceColor = colours{i};
        end

        hold on

        [numgroups, numbars] = size(to_plot);
        groupwidth = min(0.8, numbars / (numbars + 1.5));

        for i = 1:numbars
            % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange: https://uk.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab-7-13-r2011b
            x = (1:numgroups) - groupwidth / 2 + (2 * i - 1) * groupwidth / (2 * numbars); % Aligning error bar with individual bar
            errorbar(x, to_plot(:, i), error_to_plot(:, i), 'k', 'linestyle', 'none');
        end

        if rfig == 2 %| rfig==4
            temp1 = sprintf('MT Non-TMS');
            temp2 = sprintf('MT TMS');
            temp3 = sprintf('MIP Non-TMS');
            temp4 = sprintf('MIP TMS');
            %         legend({'Non-TMS  D-ipsilateral'  'Non-TMS  D-contralateral' 'TMS  D-ipsilateral' 'TMS D-contralateral'})
            legend(temp1, temp2, temp3, temp4)
        end

        if criterion_name(1) == 'b'
            %         ylim([-.11 .07])
            ylim([- .04 .1])
            set(gca, 'YTick', [- .1 - .05 0 .05 .1])
        elseif criterion_name(1) == 'd'
            %         ylim([-.07 .05])
            ylim([- .04 .07])
            set(gca, 'YTick', [- .04 - .02 0 .02 .04 .06])
        elseif criterion_name(1) == 'h'
            %         ylim([-.05 .025])
            ylim([- .02 .04])
            set(gca, 'YTick', [- .02 0 .02 .04])
        end

        set(gca, 'xticklabel', {'MT' 'MIP'})
        ticks = get(gca, 'YTick');
        title(Regs{rfig - 1})

    end

end

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
