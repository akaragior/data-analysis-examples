% Representational Similarity Analysis (RSA) on Time-Frequency MEG data.
%
% Searchlight RSA across time, frequency and sensors, on MEG data that were
% preprocessed and epoched in FieldTrip.
%
% Decoding, for each participant (runDecoding):
%   1. computes log power in sliding time-frequency (TF) windows,
%   2. assigns every trial to one of 6 classes (3 conditions x 2 sequences),
%   3. estimates the noise covariance between channels and derives a whitening
%      matrix for every sensor searchlight and for the whole sensor array,
%   4. computes the cross-validated Mahalanobis (crossnobis) distance between
%      every pair of classes in every TF bin, averaged over repeated stratified
%      k-fold partitions of the trials.
%   Output, one file per participant: a 6 x 6 representational dissimilarity
%   matrix (RDM) for every time x frequency bin, for each searchlight and for
%   the whole array, plus the settings used.
%
% Group analysis (runStatistics, makeFigures):
%   5. reduces every RDM to contrasts: sequence specificity within each
%      condition, condition specificity, and generalization of the sequence
%      code across conditions,
%   6. tests them with cluster-based permutation tests over time and frequency
%      (whole array) and over sensors (searchlight), in the preparation window,
%   7. plots model and mean RDMs, time-frequency maps with the significant clusters, and
%      searchlight topographies.
%
% Requires FieldTrip (tested with fieldtrip-20240515) and MATLAB R2020a or
% later (exportgraphics).


clear

%% Paths (machine-specific)
addpath('C:\Toolboxes\fieldtrip-20240515');
ft_defaults;

baseDir    = 'D:\MEG\4backup\';
megDataDir = [baseDir 'PortfolioRSA\data\simulated\'];   % where simulate_data.m writes

%% What to run
runDecoding   = true;    % per-participant crossnobis RDMs (slow)
runStatistics = true;    % group cluster statistics; if false, load the saved ones
makeFigures   = true;
saveFigures   = true;    % PNG copies of the figures, next to the results

%% Analysis settings
band           = 'low';    % TF band preset, see bandPresets below
epochLatency   = [-2.6 1.2];   % s, relative to the first key press
regularization = 0.01;         % covariance shrinkage, as a fraction of the mean channel variance
nFolds         = 5;            % cross-validation folds
nReps          = 50;           % repetitions of the fold partition, each with a fresh shuffle
randomSeed     = 42;           % participant k uses randomSeed + k
nSubjects      = 23;

outputDir = [baseDir 'PortfolioRSA\' band '\'];
if ~isfolder(outputDir)
    mkdir(outputDir);
end

% Classes: 3 conditions x 2 sequences, numbered condition by condition.
% conditionOrder lists the condition codes used at epoching, in class
% order: 2 = silent production, 1 = production with tones, 3 = passive
% listening to tones.
conditionOrder = [2 1 3];
conditionNames = {'MotorSilent' 'MotorTones' 'PassiveTones'};
classNames     = {'Sil1' 'Sil2' 'MTo1' 'MTo2' 'PTo1' 'PTo2'};
nSequences     = 2;
nClasses       = numel(classNames);
classPairs     = nchoosek(1:nClasses, 2);   % all 15 pairs of classes

%% Time-frequency settings
% Band presets. The window length sets the frequency resolution; the lowest
% frequency is 3/windowLength, so every window spans at least 3 cycles.
bandPresets.verylow = struct('windowLength', 2,    'freqStep', 0.25, 'freqMax', 7);    % 1.5-7 Hz
bandPresets.low     = struct('windowLength', 0.75, 'freqStep', 0.5,  'freqMax', 25);   % 4-25 Hz
bandPresets.high    = struct('windowLength', 0.3,  'freqStep', 2,    'freqMax', 60);   % 10-60 Hz
preset       = bandPresets.(band);
windowLength = preset.windowLength;   % s
lowestFreq   = 3 / windowLength;      % Hz, 3 cycles per window

% Configuration struct for FieldTrip TF analysis.
tfCfg            = [];
tfCfg.method     = 'mtmconvol';       % sliding-window FFT
tfCfg.taper      = 'hanning';
tfCfg.output     = 'pow';
tfCfg.keeptrials = 'yes';
tfCfg.pad        = 'nextpow2';
tfCfg.foi        = lowestFreq:preset.freqStep:preset.freqMax;
tfCfg.t_ftimwin  = windowLength * ones(size(tfCfg.foi));   % same window length at every frequency
tfStep           = 0.02;              % s, spacing of window centers

% A window centered at t spans [t - windowLength/2, t + windowLength/2], so
% this is the last center whose window ends at the press or tone onset.
% This is important so that we avoid bleed of post press/tone onset into
% the preparation phase.
lastPrePressCenter = -windowLength / 2;

%% Sensor searchlights
% Searchlight i is channel i and its nNeighbours nearest sensors, by
% Euclidean distance between sensor positions. The searchlights are computed
% for each participant from the sensor definition in the data (grad), so
% they always follow the channel order of the data (see sensorNeighbours).
nNeighbours = 6;   % searchlights of 7 sensors

%% Participant loop
if runDecoding
    for iSubject = 1:nSubjects
        fprintf('Participant %d\n', iSubject);

        % --- 1. Load epochs --------------------------------------------------
        % Saved after ft_timelockanalysis with kept trials: tmlck_test holds
        % the trials, trialinfo_test their condition and sequence codes.
        subjectData = load([megDataDir 'epocheddata_S' num2str(iSubject)]);
        trialInfo   = subjectData.trialinfo_test;

        selectCfg         = [];
        selectCfg.latency = epochLatency;
        timelock = ft_selectdata(selectCfg, subjectData.tmlck_test);
        clear subjectData;

        if any(~isfinite(timelock.trial(:)))
            warning('Participant %d has non-finite samples, which break the TF analysis. Skipped.', iSubject);
            continue
        end

        % Channel sets: one searchlight per channel, in the channel order of
        % the data, plus the whole array as the last set
        neighbours  = sensorNeighbours(timelock.grad, timelock.label, nNeighbours);
        channelSets = [num2cell(neighbours, 2); {true(1, numel(timelock.label))}];
        nChanSets   = numel(channelSets);

        % --- 2. Time-frequency decomposition ---------------------------------
        % Window centers must stay windowLength/2 inside the epoch, or the
        % window runs past the data (FieldTrip returns NaN there). Centers are
        % snapped to the tfStep grid, with one extra step of margin at each
        % end.
        epochEdges  = timelock.time([1 end]);
        firstCenter = ceil((epochEdges(1) + windowLength / 2) / tfStep) * tfStep + tfStep;
        lastCenter  = floor((epochEdges(2) - windowLength / 2) / tfStep) * tfStep - tfStep;
        tfCfg.toi   = firstCenter:tfStep:lastCenter;

        tfData = ft_freqanalysis(tfCfg, timelock);   % powspctrm: trials x channels x freqs x times
        clear timelock;

        % --- 3. Features and class labels ------------------------------------
        % Class 1-6 = Sil1 Sil2 MTo1 MTo2 PTo1 PTo2. Trials of other conditions
        % or sequences get conditionIndex 0 and are dropped.
        conditionCode = trialInfo.labelCond(:);
        sequenceCode  = trialInfo.labelSeq(:);
        [~, conditionIndex] = ismember(conditionCode, conditionOrder);
        keepTrials = conditionIndex > 0 & ismember(sequenceCode, 1:nSequences);
        classLabel = nSequences * (conditionIndex(keepTrials) - 1) + sequenceCode(keepTrials);

        % Features: log power, which is closer to Gaussian than raw power and
        % turns multiplicative power changes into additive ones
        logPower = log(tfData.powspctrm(keepTrials, :, :, :));
        logPower = permute(logPower, [1 2 4 3]);   % trials x channels x times x freqs
        tfTime        = tfData.time;               % s, window centers
        tfFreq        = tfData.freq;               % Hz
        channelLabels = tfData.label;
        clear tfData;

        [nTrials, nChannels, nTimeBins, nFreqs] = size(logPower);

        % Merge time and frequency into one axis of TF bins. Each bin is
        % analysed independently, so this only lets one computation cover all
        % bins at once. The reshape after the distance computation undoes this
        % exactly.
        nTFBins  = nTimeBins * nFreqs;
        logPower = reshape(logPower, nTrials, nChannels, nTFBins);

        % --- 4. Cross-validation folds ---------------------------------------
        % Repeated stratified k-fold. In each repetition, every class's trials
        % are shuffled and split into folds, so each fold holds about 1/nFolds
        % of every class. foldAssignment(trial, rep) = fold of that trial.
        % Random seed per participant, right before the only random step, so
        % each participant's folds are reproducible on their own.
        rng(randomSeed + iSubject);
        foldAssignment = zeros(nTrials, nReps);
        for iRep = 1:nReps
            for iClass = 1:nClasses
                classTrials = find(classLabel == iClass);
                classTrials = classTrials(randperm(numel(classTrials)));
                foldAssignment(classTrials, iRep) = mod(0:numel(classTrials) - 1, nFolds) + 1;
            end
        end

        % --- 5. Noise covariance ---------------------------------------------
        % How channels co-vary from trial to trial once the class effect is
        % removed. Within one class and TF bin, all trials share the same true
        % pattern, so subtracting the class mean leaves only noise. The
        % residual cross-products R'*R (channels x channels) are summed over
        % classes and bins, then divided by the degrees of freedom: each
        % subtracted mean costs one, so n - 1 per class and bin. Pooling over
        % classes assumes they share one noise structure; pooling over bins
        % assumes it is stable across the band and the epoch.
        noiseCov = zeros(nChannels);
        noiseDof = 0;
        for iBin = 1:nTFBins
            for iClass = 1:nClasses
                classData = logPower(classLabel == iClass, :, iBin);   % trials x channels
                residuals = classData - mean(classData, 1);
                noiseCov  = noiseCov + residuals' * residuals;
                noiseDof  = noiseDof + size(classData, 1) - 1;
            end
        end
        noiseCov = noiseCov / noiseDof;   % noiseDof = nTFBins * (nTrials - nClasses)

        % --- 6. Whitening matrices -------------------------------------------
        % Whitening transforms the channels so that their noise is uncorrelated
        % and equally large in every direction, like a z-score for many channels
        % at once. Distances between whitened patterns are then in units of
        % noise (Mahalanobis distances): a difference along a direction in which
        % trials vary a lot counts little, one along a direction in which they
        % barely vary counts a lot.
        % For a channel set with noise covariance C, W = inv(chol(C, 'lower'))
        % does this: W*C*W' = I. Before that, a small fraction of the mean
        % channel variance is added to every channel (regularization), so that
        % directions that look almost noise-free by chance are not blown up.
        % Each searchlight is whitened with the noise covariance of its own
        % channels only.
        whiteners = cell(nChanSets, 1);
        for iSet = 1:nChanSets
            inChanSet     = channelSets{iSet};
            nFeatures = sum(inChanSet);
            setCov    = noiseCov(inChanSet, inChanSet);
            setCov    = setCov + eye(nFeatures) * regularization * trace(setCov) / nFeatures;
            whiteners{iSet} = chol(setCov, 'lower') \ eye(nFeatures);
        end

        % --- 7. Crossnobis distances -----------------------------------------
        % For classes A and B, fold f gives the whitened difference pattern
        %   a_f = W * (mean_A,f - mean_B,f) = delta + noise_f,
        % where delta is the true difference and noise_f is independent between
        % folds. Squaring one fold's estimate is biased upward by its own
        % noise:
        %   E[a_f' * a_f] = |delta|^2 + E|noise_f|^2
        % Multiplying two different folds' estimates is not, since their noise
        % is independent:
        %   E[a_f' * a_g] = |delta|^2   for f ~= g
        % The crossnobis distance averages a_f' * a_g over all fold pairs with
        % f ~= g.
        distanceSum = zeros(nClasses, nClasses, nTFBins, nChanSets);
        for iRep = 1:nReps

            % Mean pattern of every class in every fold
            foldClassMeans = zeros(nClasses, nFolds, nChannels, nTFBins);
            for iClass = 1:nClasses
                for iFold = 1:nFolds
                    inCell = classLabel == iClass & foldAssignment(:, iRep) == iFold;
                    foldClassMeans(iClass, iFold, :, :) = mean(logPower(inCell, :, :), 1);
                end
            end

            repDistance = zeros(nClasses, nClasses, nTFBins, nChanSets);   % RDMs from this repetition
            for iSet = 1:nChanSets
                inChanSet = channelSets{iSet};
                nFeatures = sum(inChanSet);

                % Whiten all fold means at once: channels to the first dimension
                % (channels x classes x folds x bins), everything else flattened,
                % one matrix product, shape restored
                setMeans      = permute(foldClassMeans(:, :, inChanSet, :), [3 1 2 4]);
                whitenedMeans = whiteners{iSet} * reshape(setMeans, nFeatures, []);
                whitenedMeans = reshape(whitenedMeans, nFeatures, nClasses, nFolds, nTFBins);

                for iPair = 1:size(classPairs, 1)
                    classA = classPairs(iPair, 1);
                    classB = classPairs(iPair, 2);

                    % a_f for every fold: channels x folds x bins
                    foldDiffs = reshape(whitenedMeans(:, classA, :, :) - whitenedMeans(:, classB, :, :), ...
                                        nFeatures, nFolds, nTFBins);

                    % Sum of a_f' * a_g over f ~= g, without a double loop:
                    % |sum_f a_f|^2 contains all nFolds^2 pairs; subtracting the
                    % nFolds same-fold terms |a_f|^2 leaves the nFolds*(nFolds-1) others
                    sumOverFolds  = sum(foldDiffs, 2);                    % sum_f a_f
                    allFoldPairs  = squeeze(sum(sumOverFolds.^2, 1));     % |sum_f a_f|^2
                    sameFoldPairs = squeeze(sum(foldDiffs.^2, [1 2]));    % sum_f |a_f|^2
                    pairDistance  = (allFoldPairs - sameFoldPairs) / (nFolds * (nFolds - 1));

                    % Per channel, so searchlights and the whole array share a scale
                    pairDistance = pairDistance / nFeatures;

                    repDistance(classA, classB, :, iSet) = pairDistance;
                    repDistance(classB, classA, :, iSet) = pairDistance;
                end
            end

            distanceSum = distanceSum + repDistance;
            fprintf('.');   % one dot per repetition
        end
        fprintf('\n');

        % Each repetition is unbiased; averaging reduces the variance that
        % comes from any single partition of the trials. Then reshape back
        % to original dimensions again.
        crossnobis = distanceSum / nReps;
        crossnobis = reshape(crossnobis, nClasses, nClasses, nTimeBins, nFreqs, nChanSets);
        % Searchlights: class x class x time x freq x searchlight center
        % Whole array:  class x class x time x freq
        crossnobisSearchlight = crossnobis(:, :, :, :, 1:nChanSets - 1);
        crossnobisWholeArray  = crossnobis(:, :, :, :, nChanSets);

        % --- 8. Save ---------------------------------------------------------------
        analysisSettings = struct();
        analysisSettings.band               = band;
        analysisSettings.tfCfg              = tfCfg;
        analysisSettings.time               = tfTime;
        analysisSettings.freq               = tfFreq;
        analysisSettings.lastPrePressCenter = lastPrePressCenter;
        analysisSettings.epochLatency       = epochLatency;
        analysisSettings.regularization     = regularization;
        analysisSettings.nFolds             = nFolds;
        analysisSettings.nReps              = nReps;
        analysisSettings.randomSeed         = randomSeed + iSubject;   % seed used for this participant
        analysisSettings.classNames         = classNames;
        analysisSettings.conditionNames     = conditionNames;
        analysisSettings.nNeighbours        = nNeighbours;
        analysisSettings.channelLabels      = channelLabels;           % also the order of searchlight centers
        analysisSettings.dimordSearchlight  = 'class_class_time_freq_chan';
        analysisSettings.dimordWholeArray   = 'class_class_time_freq';

        save(sprintf('%sdecodingresults_S%d.mat', outputDir, iSubject), ...
             'crossnobisSearchlight', 'crossnobisWholeArray', 'analysisSettings', '-v7.3');
    end

    fprintf('Decoding finished.\n');
end

%% Group analysis settings
clusterFormingAlpha = 0.05;   % sample-level threshold for joining a cluster
significanceLevel   = 0.05;   % cluster-level threshold
nRandomizations     = 1000;   % permutations per test
statisticsFile      = [outputDir 'groupStatistics.mat'];
layoutFile          = 'neuromag306mag_helmet.mat';
neighbourTemplate   = 'neuromag306mag_neighb.mat';   % sensor neighbours for spatial clusters

%% Load the RDMs and compute the metrics
if runStatistics || makeFigures

    % --- 1. Load the participants' RDMs ----------------------------------------
    % rdm: participants x class x class x time x freq x channel set. As in the
    % decoding, the searchlights come first and the whole array last.
    clear rdm
    for iSubject = nSubjects:-1:1   % backwards, for preallocation with 1st entry
        result = load([outputDir 'decodingresults_S' num2str(iSubject)]);
        rdm(iSubject, :, :, :, :, :) = cat(5, result.crossnobisSearchlight, result.crossnobisWholeArray);
        % searchlight k must be centered on the same sensor in every participant
        if iSubject < nSubjects
            assert(isequal(result.analysisSettings.channelLabels, channelLabels), ...
                   'Participant %d: channel order differs from participant %d.', iSubject, nSubjects);
        end
        channelLabels = result.analysisSettings.channelLabels;
    end
    analysisSettings = result.analysisSettings;
    clear result
    time          = analysisSettings.time;
    freq          = analysisSettings.freq;
    channelLabels = analysisSettings.channelLabels;
    

    % --- 2. RDM metrics --------------------------------------------------------
    % Each metric condenses the 6 x 6 RDM into one value per participant, time,
    % frequency and channel set (participants x time x freq x channel set).
    % Classes: 1-2 = Sil1 Sil2, 3-4 = MTo1 MTo2, 5-6 = PTo1 PTo2.

    metrics = struct();
    % Sequence specificity: distance between the two sequences of a condition
    metrics.seq_Sil = squeeze(rdm(:, 1, 2, :, :, :));
    metrics.seq_MTo = squeeze(rdm(:, 3, 4, :, :, :));
    metrics.seq_PTo = squeeze(rdm(:, 5, 6, :, :, :));

    % Condition specificity: distance between two conditions for the same
    % sequence, averaged over the two sequences (pairs with different
    % sequences would add sequence information)
    metrics.cnd_Sil_MTo = squeeze((rdm(:, 1, 3, :, :, :) + rdm(:, 2, 4, :, :, :)) / 2);
    metrics.cnd_Sil_PTo = squeeze((rdm(:, 1, 5, :, :, :) + rdm(:, 2, 6, :, :, :)) / 2);
    metrics.cnd_MTo_PTo = squeeze((rdm(:, 3, 5, :, :, :) + rdm(:, 4, 6, :, :, :)) / 2);

    % Sequence generalization between two conditions: distance between
    % different sequences minus distance between the same sequence, across the
    % two conditions. Above zero when they share a sequence code.
    metrics.gen_SilxMTo = squeeze((rdm(:, 1, 4, :, :, :) + rdm(:, 2, 3, :, :, :)) / 2 ...
                             - (rdm(:, 1, 3, :, :, :) + rdm(:, 2, 4, :, :, :)) / 2);
    metrics.gen_SilxPTo = squeeze((rdm(:, 1, 6, :, :, :) + rdm(:, 2, 5, :, :, :)) / 2 ...
                             - (rdm(:, 1, 5, :, :, :) + rdm(:, 2, 6, :, :, :)) / 2);
    metrics.gen_MToxPTo = squeeze((rdm(:, 3, 6, :, :, :) + rdm(:, 4, 5, :, :, :)) / 2 ...
                             - (rdm(:, 3, 5, :, :, :) + rdm(:, 4, 6, :, :, :)) / 2);

    metricNames     = fieldnames(metrics);
    sequenceMetrics = {'seq_Sil' 'seq_MTo' 'seq_PTo'};   % in the order of conditionNames

    % Only the whole-array RDMs are needed from here on, for the RDM figure
    rdmWholeArray = squeeze(rdm(:, :, :, :, :, end));   % participants x class x class x time x freq
    clear rdm
end

%% Group statistics
% Cluster-based permutation tests (Maris & Oostenveld, 2007), restricted to
% window centers up to lastPrePressCenter: every tested window ends before the
% first press or tone, so preparation is tested without post-onset activity.
% Each test compares two sets of maps, paired by participant: a metric against
% zero, or sequence specificity in one condition against another.
if runStatistics
    statCfg                  = [];
    statCfg.parameter        = 'powspctrm';
    statCfg.method           = 'montecarlo';
    statCfg.statistic        = 'depsamplesT';
    statCfg.correctm         = 'cluster';
    statCfg.clusterstatistic = 'maxsum';
    statCfg.clusteralpha     = clusterFormingAlpha;
    statCfg.alpha            = significanceLevel;
    statCfg.numrandomization = nRandomizations;
    statCfg.latency          = [time(1) lastPrePressCenter];
    statCfg.randomseed       = 1;
    statCfg.design           = [1:nSubjects, 1:nSubjects; ones(1, nSubjects), 2 * ones(1, nSubjects)];
    statCfg.uvar             = 1;   % design row 1: participant
    statCfg.ivar             = 2;   % design row 2: first or second set of maps

    % Metrics against zero: one-sided, since distances and generalization are
    % expected to be positive
    oneSidedCfg             = statCfg;
    oneSidedCfg.tail        = 1;
    oneSidedCfg.clustertail = 1;

    % Differences between conditions: two-sided
    twoSidedCfg             = statCfg;
    twoSidedCfg.tail        = 0;
    twoSidedCfg.clustertail = 0;
    twoSidedCfg.correcttail = 'prob';

    % Whole array: clusters over time and frequency
    stats    = struct();
    zeroMaps = zeros(size(metrics.seq_Sil));
    zeroWholeArrayFT = toFieldTrip(zeroMaps(:, :, :, end), {'wholearray'}, time, freq);
    for iMetric = 1:numel(metricNames)
        name = metricNames{iMetric};
        wholeArrayDataFT.(name) = toFieldTrip(metrics.(name)(:, :, :, end), {'wholearray'}, time, freq);
        stats.wholeArray.(name) = ft_freqstatistics(oneSidedCfg, wholeArrayDataFT.(name){:}, zeroWholeArrayFT{:});
    end
    stats.wholeArray.seq_SilVsMTo = ft_freqstatistics(twoSidedCfg, wholeArrayDataFT.seq_Sil{:}, ...
                                                     wholeArrayDataFT.seq_MTo{:});
    stats.wholeArray.seq_SilVsPTo = ft_freqstatistics(twoSidedCfg, wholeArrayDataFT.seq_Sil{:}, ...
                                                     wholeArrayDataFT.seq_PTo{:});
    stats.wholeArray.seq_MToVsPTo = ft_freqstatistics(twoSidedCfg, wholeArrayDataFT.seq_MTo{:}, ...
                                                     wholeArrayDataFT.seq_PTo{:});

    % Searchlight: clusters over sensors, time and frequency, for sequence
    % specificity. The topography is what separates sources of sequence
    % information, e.g. visual responses to the sequence cue from planning.
    zeroSearchlightFT = toFieldTrip(zeroMaps(:, :, :, 1:end - 1), channelLabels, time, freq);

    neighbourCfg              = [];
    neighbourCfg.method       = 'template';
    neighbourCfg.template     = neighbourTemplate;
    searchlightCfg            = oneSidedCfg;
    searchlightCfg.neighbours = ft_prepare_neighbours(neighbourCfg, zeroSearchlightFT{1});
    for iMetric = 1:numel(sequenceMetrics)
        name = sequenceMetrics{iMetric};
        searchlightDataFT        = toFieldTrip(metrics.(name)(:, :, :, 1:end - 1), channelLabels, time, freq);
        stats.searchlight.(name) = ft_freqstatistics(searchlightCfg, searchlightDataFT{:}, zeroSearchlightFT{:});
    end
    clear zeroMaps wholeArrayDataFT zeroSearchlightFT searchlightDataFT

    save(statisticsFile, 'stats', '-v7.3');

elseif makeFigures
    load(statisticsFile, 'stats');
end

%% Figures
% Time-frequency maps show the whole epoch, with color scales fitted to the
% tested (preparation) window, so stronger responses after onset may saturate.
if makeFigures
    figures       = struct('handle', {}, 'fileName', {});    
     
    % Group mean of every metric, whole array: time x freq
    for iMetric = 1:numel(metricNames)
        name = metricNames{iMetric};
        groupMean.(name) = squeeze(mean(metrics.(name)(:, :, :, end), 1));
    end

    % --- Model RDMs and mean RDMs, whole array, all frequencies ----------------
    % Top row: two model RDMs, for comparison with the data (the metrics above
    % are not fitted to them). Sequence model: 1 where two classes differ in
    % sequence, the same sequence code in every condition. Condition model: 1
    % where a motor condition meets PassiveTones, MotorTones like MotorSilent
    % (as expected before the first press; around the presses and tones,
    % MotorTones would sit in between).
    % Bottom row: the data. One color scale per window: distances after onset
    % are much larger.
    classSequence  = repmat(1:nSequences, 1, numel(conditionNames));   % 1 2 1 2 1 2
    classIsPassive = repelem(strcmp(conditionNames, 'PassiveTones'), nSequences);
    modelRDMs      = {double(classSequence' ~= classSequence), double(classIsPassive' ~= classIsPassive)};
    modelTitles    = {'Model: sequence', 'Model: motor vs passive'};

    inPreparation = time <= lastPrePressCenter;
    rdmWindows = {inPreparation, ~inPreparation};
    rdmTitles  = {{'Data: preparation', sprintf('window centers up to %.2f s', lastPrePressCenter)}, ...
                  {'Data: presses or tones', 'windows that include them'}};
    handle = figure('Name', 'Mean RDMs', 'Color', 'w', 'Position', [100 100 950 820]);
    for iPanel = 1:4
        subplot(2, 2, iPanel);
        if iPanel <= 2
            panelRDM = modelRDMs{iPanel};
            imagesc(panelRDM);
            title(modelTitles{iPanel});
            colorLabel = 'model dissimilarity';
        else
            panelRDM = squeeze(mean(rdmWholeArray(:, :, :, rdmWindows{iPanel - 2}, :), [1 4 5]));
            imagesc(panelRDM);
            title(rdmTitles{iPanel - 2});
            colorLabel = 'crossnobis distance';
        end
        colorbarHandle = colorbar;
        colorbarHandle.Label.String = colorLabel;
        axis square;
        colormap(gca, cmapWhiteRed(256));
        caxis([0 max(panelRDM(:))]);
        set(gca, 'XTick', 1:nClasses, 'XTickLabel', classNames, ...
                 'YTick', 1:nClasses, 'YTickLabel', classNames);
    end
    sgtitle(sprintf('Model RDMs and mean crossnobis RDMs, whole array, %.1f-%.1f Hz (n = %d)', ...
                    freq(1), freq(end), nSubjects));
    figures(end + 1) = struct('handle', handle, 'fileName', 'rdm');

    % --- Condition specificity -------------------------------------------------
    conditionMetrics = {'cnd_Sil_MTo' 'cnd_Sil_PTo' 'cnd_MTo_PTo'};  
    limits = symmetricLimits({groupMean.cnd_Sil_MTo, groupMean.cnd_Sil_PTo, groupMean.cnd_MTo_PTo}, ...
                             inPreparation);
    pairTitlesCondition = {'MotorSilent vs MotorTones' 'MotorSilent vs PassiveTones' 'MotorTones vs PassiveTones'};
    handle = figure('Name', 'Condition specificity', 'Color', 'w', 'Position', [100 100 1260 380]);
    for iPair = 1:3
        name = conditionMetrics{iPair};
        subplot(1, 3, iPair);
        drawTFMap(time, freq, groupMean.(name), stats.wholeArray.(name), limits, lastPrePressCenter, ...
                  'crossnobis distance');
        title(pairTitlesCondition{iPair});
    end
    sgtitle('Condition specificity');
    figures(end + 1) = struct('handle', handle, 'fileName', 'condition_specificity');

    % --- Sequence specificity: per condition, and between conditions ---------
    % 3 x 3 grid. Diagonal: each condition, with the clusters of its test
    % against zero. Upper triangle: row condition minus column condition, with
    % the clusters of the two-sided test of that same difference.
    differences      = {groupMean.seq_Sil - groupMean.seq_MTo, ...
                        groupMean.seq_Sil - groupMean.seq_PTo, ...
                        groupMean.seq_MTo - groupMean.seq_PTo};
    differenceTests  = {'seq_SilVsMTo' 'seq_SilVsPTo' 'seq_MToVsPTo'};
    differenceTitles = {'MotorSilent - MotorTones' 'MotorSilent - PassiveTones' ...
                        'MotorTones - PassiveTones'};
    sequenceLimits   = symmetricLimits({groupMean.seq_Sil, groupMean.seq_MTo, groupMean.seq_PTo}, ...
                                       inPreparation);
    differenceLimits = symmetricLimits(differences, inPreparation);
    diagonalPanels   = [1 5 9];
    upperPanels      = [2 3 6];

    handle = figure('Name', 'Sequence specificity', 'Color', 'w', 'Position', [100 100 1250 1000]);
    for iCondition = 1:3
        name = sequenceMetrics{iCondition};
        subplot(3, 3, diagonalPanels(iCondition));
        drawTFMap(time, freq, groupMean.(name), stats.wholeArray.(name), sequenceLimits, lastPrePressCenter, ...
                  'crossnobis distance');
        title(conditionNames{iCondition});
    end
    for iPair = 1:3
        subplot(3, 3, upperPanels(iPair));
        drawTFMap(time, freq, differences{iPair}, stats.wholeArray.(differenceTests{iPair}), ...
                  differenceLimits, lastPrePressCenter, 'difference in crossnobis distance');
        title(differenceTitles{iPair});
    end
    sgtitle('Sequence specificity (diagonal) and its differences between conditions');
    figures(end + 1) = struct('handle', handle, 'fileName', 'sequence_specificity');

    % --- Generalization of the sequence code across conditions --------------
    generalizationMetrics = {'gen_SilxMTo' 'gen_SilxPTo' 'gen_MToxPTo'};
    limits = symmetricLimits({groupMean.gen_SilxMTo, groupMean.gen_SilxPTo, groupMean.gen_MToxPTo}, ...
                             inPreparation);
    pairTitlesGeneralization = {'MotorSilent X MotorTones' 'MotorSilent X PassiveTones' 'MotorTones X PassiveTones'};
    handle = figure('Name', 'Sequence generalization', 'Color', 'w', 'Position', [100 100 1260 380]);
    for iPair = 1:3
        name = generalizationMetrics{iPair};
        subplot(1, 3, iPair);
        drawTFMap(time, freq, groupMean.(name), stats.wholeArray.(name), limits, lastPrePressCenter, ...
                  'distance: different - same sequence');
        title(pairTitlesGeneralization{iPair});
    end
    sgtitle('Sequence generalization across conditions');
    figures(end + 1) = struct('handle', handle, 'fileName', 'sequence_generalization');

    % --- Searchlight topographies ---------------------------------------------
    % One figure per metric family (sequence, condition specificity). Top row:
    % preparation; bottom row: windows that include the presses or tones. Each
    % map is the group mean over its window, in one frequency band per family
    % and window, chosen for what it shows (a display choice; the statistics
    % do not depend on it). Searchlight tests exist only for sequence
    % specificity in the preparation window: there, sensors are marked where a
    % significant cluster covers at least half of the plotted time-frequency
    % bins (a cluster can spread through the whole array, so marking every
    % sensor it touches says little about location).
    familyBands = {[4 8],   [10 16]; ...   % sequence: theta and cue responses | second press
                   [17 23], [4 8]};        % condition: beta desynchronization | tone responses

    topoCfg                 = [];
    topoCfg.layout          = layoutFile;
    topoCfg.parameter       = 'powspctrm';
    topoCfg.colormap        = cmapBlueWhiteRed(256);
    topoCfg.style           = 'straight';   % colors only, no contour lines
    topoCfg.comment         = 'no';
    topoCfg.figure          = 'gca';
    topoCfg.highlightsymbol = 'o';
    topoCfg.highlightcolor  = [0 0 0];

    familyMetrics = {sequenceMetrics, conditionMetrics};
    familyTitles  = {conditionNames, pairTitlesCondition};
    familyNames   = {'Sequence specificity', 'Condition specificity'};
    familyFiles   = {'sequence_searchlight', 'condition_searchlight'};
    windowTimes   = {inPreparation, ~inPreparation};
    windowLabels  = {'preparation', 'presses or tones'};

    for iFamily = 1:2
        handle = figure('Name', [familyNames{iFamily} ', searchlight'], ...
            'Color', 'w', 'Position', [100 100 1200 800]);
        for iWindow = 1:2
            % Topographies (sensors x 1), one per metric of the family
            band   = familyBands{iFamily, iWindow};
            inBand = freq >= band(1) & freq <= band(2);
            maps   = cell(1, 3);
            for iMetric = 1:3
                name = familyMetrics{iFamily}{iMetric};
                maps{iMetric} = squeeze(mean(metrics.(name)(:, windowTimes{iWindow}, inBand, 1:end - 1), [1 2 3]));
            end
            rowValues    = [maps{:}];
            topoCfg.zlim = [-1 1] * max(abs(rowValues(:)));   % one color scale per row
            for iMetric = 1:3
                if iFamily == 1 && iWindow == 1
                    % sensors where a significant cluster covers at least half of
                    % the plotted bins (the test covered all frequencies, so
                    % stat.freq equals freq)
                    stat     = stats.searchlight.(familyMetrics{iFamily}{iMetric});
                    coverage = mean(mean(stat.mask(:, inBand, :), 2), 3);
                    topoCfg.highlight        = 'on';
                    topoCfg.highlightchannel = stat.label(coverage >= 0.5);
                else
                    topoCfg.highlight = 'off';
                end
                % the map as a FieldTrip structure with a single time and frequency
                topoData = toFieldTrip(reshape(maps{iMetric}, 1, 1, 1, []), channelLabels, ...
                    mean(time(windowTimes{iWindow})), mean(freq(inBand)));
                subplot(2, 3, (iWindow - 1) * 3 + iMetric);
                ft_topoplotTFR(topoCfg, topoData{1});
                colorbarHandle = colorbar;
                colorbarHandle.Label.String = 'crossnobis distance';
                title({familyTitles{iFamily}{iMetric}, ...
                    sprintf('%s, %g-%g Hz', windowLabels{iWindow}, band(1), band(2))});
            end
        end
        sgtitle(sprintf('%s per searchlight (top: preparation, window centers up to %.2f s)', ...
                        familyNames{iFamily}, lastPrePressCenter));
        figures(end + 1) = struct('handle', handle, 'fileName', familyFiles{iFamily});
    end

    if saveFigures
        for iFigure = 1:numel(figures)
            exportgraphics(figures(iFigure).handle, [outputDir figures(iFigure).fileName '.png'], ...
                           'Resolution', 200);
        end
    end
end

%% Local functions
function data = toFieldTrip(values, labels, time, freq)
    % One FieldTrip freq structure per participant, from values arranged as
    % participants x time x freq x channels
    data = cell(1, size(values, 1));
    for iParticipant = 1:numel(data)
        data{iParticipant}.label     = labels;
        data{iParticipant}.time      = time;
        data{iParticipant}.freq      = freq;
        data{iParticipant}.dimord    = 'chan_freq_time';
        data{iParticipant}.powspctrm = permute(values(iParticipant, :, :, :), [4 3 2 1]);
    end
end

function drawTFMap(time, freq, values, stat, limits, lastPrePressCenter, colorLabel)
    % Time-frequency map of values (time x freq) with the significant clusters
    % of stat outlined: solid where positive, dashed where negative.
    % colorLabel: label of the color bar.
    imagesc(time, freq, values');
    set(gca, 'YDir', 'normal');
    colormap(gca, cmapBlueWhiteRed(256));
    caxis(limits);
    colorbarHandle = colorbar;
    colorbarHandle.Label.String = colorLabel;
    hold on;

    % The test covers only the preparation window, so its mask is placed on
    % the full freq x time grid
    testedTimes = nearest(time, stat.time(1)):nearest(time, stat.time(end));
    positive    = false(numel(freq), numel(time));
    negative    = false(numel(freq), numel(time));
    positive(:, testedTimes) = squeeze(stat.mask & stat.stat > 0);
    negative(:, testedTimes) = squeeze(stat.mask & stat.stat < 0);
    if any(positive(:))
        contour(time, freq, double(positive), [0.5 0.5], 'LineColor', 'k', 'LineWidth', 1.5);
    end
    if any(negative(:))
        contour(time, freq, double(negative), [0.5 0.5], 'LineColor', 'k', 'LineWidth', 1.5, ...
                'LineStyle', '--');
    end

    xline(0, 'k');                      % first press or tone
    xline(lastPrePressCenter, 'k--');   % last window that ends before it
    hold off;
    xlabel('window center (s)');
    ylabel('frequency (Hz)');
end

function limits = symmetricLimits(maps, rows)
    % Color limits symmetric around zero that cover the given rows (time bins)
    % of every map in a cell array
    peak = 0;
    for iMap = 1:numel(maps)
        values = maps{iMap}(rows, :);
        peak   = max(peak, max(abs(values(:))));
    end
    limits = [-1 1] * max(peak, eps);
end

function map = cmapBlueWhiteRed(nColors)
    % Diverging colormap: blue below zero, white at zero, red above
    blue = [0.13 0.40 0.67];
    red  = [0.70 0.09 0.17];
    ramp = linspace(0, 1, nColors / 2)';
    map  = [blue + ramp .* (1 - blue); 1 - ramp .* (1 - red)];
end

function map = cmapWhiteRed(nColors)
    % Sequential colormap: white at zero, red at the maximum
    red  = [0.70 0.09 0.17];
    ramp = linspace(0, 1, nColors)';
    map  = 1 - ramp .* (1 - red);
end

function neighbours = sensorNeighbours(grad, labels, nNeighbours)
    % Searchlight definition: row i marks channel i and its nNeighbours
    % nearest sensors, by Euclidean distance between sensor positions. Rows
    % and columns follow the order of labels.
    [found, gradIndex] = ismember(labels, grad.label);
    assert(all(found), 'Some channels are missing from the sensor definition (grad).');
    position = grad.chanpos(gradIndex, :);   % channels x 3
    sensorDistance = sqrt((position(:, 1) - position(:, 1)').^2 + ...
                          (position(:, 2) - position(:, 2)').^2 + ...
                          (position(:, 3) - position(:, 3)').^2);
    [~, byDistance] = sort(sensorDistance, 2);   % each row starts with the channel itself
    neighbours = false(numel(labels));
    for iChan = 1:numel(labels)
        neighbours(iChan, byDistance(iChan, 1:nNeighbours + 1)) = true;
    end
end