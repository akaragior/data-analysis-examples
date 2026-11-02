% SIMULATE_DATA  Simulated MEG epochs in the format of the real data.
%
% Writes one file per participant, epocheddata_S<k>.mat, holding the same two
% variables as the real epoched data, so run_decoding.m reads them unchanged:
%   tmlck_test      FieldTrip timelock structure with single trials
%                   (102 magnetometers; trial = trials x channels x samples, in T)
%                   and the sensor definition (grad)
%   trialinfo_test  per-trial condition code (labelCond) and sequence (labelSeq)
% plus simulation_truth.mat, which records what was injected.
%
% Trial timeline (0 = first press or tone):
%   -2.5 s  the sequence cue (a fractal) appears
%   -0.5 s  the cue disappears: go signal
%    0   s  press/tone 1, the same finger and tone in both sequences
%   +0.5 s  press/tone 2, different between the sequences
% The sequence goes on, but the epochs end at 1.2 s, before the third press or
% tone, as in the real data.
%
% Simulated activity (left hand):
%   Visual    occipital responses to cue onset and offset. Each sequence has its
%             own fractal in every condition, so these differ between sequences:
%             a confound that only its topography separates from the rest.
%   Theta     frontal 4-7 Hz power from shortly after cue onset to the end of
%             the epoch, with a sequence-specific topography, the same in all
%             conditions but at half amplitude in PassiveTones.
%   Mu, beta  desynchronization of the sensorimotor rhythms, right > left: mu
%             (8-12 Hz) over the somatosensory and beta (17-23 Hz) over the
%             motor hand areas (Salmelin et al., 1995). In MotorSilent and
%             MotorTones, most of it from shortly after cue onset (planning),
%             full from the go signal into the first press (execution). In
%             PassiveTones, weaker and tied to the tones: beta from the first
%             tone on (Fujioka et al., 2012), mu already from the go signal, in
%             anticipation of tones from practised sequences (Bedford et al.,
%             2025, preprint).
%   Ramps     readiness potential before the first press: a midline part
%             (SMA) from 1.5 s before, and a part over the hand area
%             contralateral to the moving hand from the go signal on (Brunia et
%             al., 2012). The midline part is 2.5 times larger before
%             sound-producing presses, the lateralized part is the same (Reznik
%             et al., 2018). In PassiveTones, a stimulus-preceding negativity
%             instead: the tones are informative, since the sequence is judged
%             against the cue afterwards. It is right-dominant fronto-temporal,
%             with an auditory part (Brunia et al., 2012). All ramps hold after
%             the first press or tone, since the sequence goes on.
%   Presses   motor command over the motor hand area, with a finger-specific
%             topography, rising from ~130 ms before each registered press and
%             peaking 50 ms before it: the motor potential, which starts ~80 ms
%             before movement onset (Brunia et al., 2012), and the button
%             registers the press a little after the movement starts. Then a
%             somatosensory evoked field ~30 ms after each press, with a
%             finger-specific topography slightly behind the motor hand area.
%             From thumb to little finger, both move medially and backwards
%             along the central sulcus (Dechent & Frahm, 2003; van Westen et
%             al., 2004). MotorSilent and MotorTones only.
%   Tones     N1m ~150 ms after each tone, with a pitch-specific topography:
%             higher tones lie more medially, i.e. deeper, in auditory cortex
%             (Pantev et al., 1996). Attenuated for self-generated tones
%             (MotorTones; Martikainen et al., 2005).
%   Noise     spatially correlated 1/f activity, posterior alpha that varies
%             from trial to trial, and sensor noise.
% Effect strengths and source positions vary between participants.
%
% Each source's sensor pattern is a simple cluster: weights that fall off with
% distance from the source position on the 2-D sensor layout. This is not a
% physical forward model; the data only need effects that are spatially
% organized in a plausible way, so that the decoding can be checked against a
% known ground truth.
%
% Requires FieldTrip (template sensor array and plotting layout).

%% Paths (machine-specific)
addpath('C:\Toolboxes\fieldtrip-20240515');
ft_defaults;

baseDir    = 'D:\MEG\4backup\';
outputDir = [baseDir 'PortfolioRSA\data\simulated\'];
if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Simulation settings
nSubjects       = 23;
nTrialsPerClass = 60;
fsample         = 100;          % Hz; covers the verylow and low bands (up to 25 Hz)
epochTime       = [-4.0 1.2];   % s; from early enough for a pre-cue TF baseline to before the 3rd press
nNoiseSources   = 30;           % sources of background 1/f activity
randomSeed      = 42;            % participant k uses randomSeed + k

% Events (s, relative to the first press or tone)
cueOnset  = -2.5;               % sequence cue onset
cueOffset = -0.5;               % sequence cue offset = go signal
itemTimes = [0 0.5];            % press/tone 1 and 2

% Conditions, with the codes of the real data. One row per condition:
%                 MotorTones  MotorSilent  PassiveTones
conditionCodes  = [1;         2;           3];
conditionNames  = {'MotorTones'; 'MotorSilent'; 'PassiveTones'};
hasPresses      = [true;      true;        false];
toneGain        = [0.6;       0;           1];     % self-generated tones are attenuated
thetaGain       = [1;         1;           0.5];   % sequence-specific theta, weaker without movement
motorDesync     = [0.5;       0.5;         0];     % mu/beta amplitude drop for planning and pressing
listenDesync    = [0;         0;           0.2];   % mu/beta amplitude drop for listening to the tones
rampLateralGain = [1;         1;           0];
rampMidlineGain = [2.5;       1;           0];     % larger before sound-producing presses
spnGain         = [0;         0;           1];     % stimulus-preceding negativity before informative tones

% Left-hand fingers: 1 = thumb ... 5 = little. Both sequences start with the
% index finger and differ at the second press. Each finger plays its own tone.
sequenceFingers = [2 3      % sequence 1: index, middle
                   2 5];    % sequence 2: index, little

% Peak amplitude (evoked responses, ramps) or RMS (oscillations, noise) in fT
% at the strongest sensor, before participant scaling.
% These are plausible guesses, not calibrated to a dataset: round numbers in a
% typical range for magnetometers. They are ordered so that the sensory
% responses are the strongest, the motor fields and slow potentials weaker, and
% the effects of interest (sequence-specific theta, mu during listening) close
% to the noise on single trials. The background (1/f activity, posterior alpha,
% sensor noise) makes single trials noisy, as in real recordings: the strongest
% evoked responses are about as large as the background, weaker effects only
% emerge across trials and participants. Only the ratios matter for the
% decoding: with log power and whitening, scaling everything by the same
% factor leaves every distance unchanged.
amp.visual      = 250;
amp.motor       = 150;
amp.sef         = 150;
amp.n1m         = 300;
amp.rampLateral = 100;
amp.rampMidline = 60;
amp.spn         = 60;
amp.theta       = 35;
amp.mu          = 100;
amp.beta        = 80;
amp.alpha       = 120;
amp.brainNoise  = 150;
amp.sensorNoise = 15;

%% Sensors
% The 102 magnetometers of FieldTrip's template Neuromag/Elekta sensor array,
% in its channel order, which is also the order of the recorded data. The
% simulated data carry this sensor definition (grad), from which
% run_decoding.m computes the searchlights, as it does for real data.
template      = load('neuromag306.mat');      % FieldTrip, template/gradiometer
grad          = template.neuromag306;
isMag         = endsWith(grad.label, '1');    % magnetometers end in 1, gradiometers in 2 or 3
channelLabels = grad.label(isMag);
nChannels     = numel(channelLabels);

% Where to place the simulated sources: 2-D positions from the plotting
% layout, centered and scaled to the unit disk (x: left to right, y: back to
% front)
layoutCfg        = [];
layoutCfg.layout = 'neuromag306mag_helmet.mat';
layout           = ft_prepare_layout(layoutCfg);
[found, layoutIndex] = ismember(channelLabels, layout.label);
assert(all(found), 'Some channels are missing from the layout.');
sensorPos = layout.pos(layoutIndex, :);
sensorPos = sensorPos - mean(sensorPos, 1);
sensorPos = sensorPos / max(sqrt(sum(sensorPos.^2, 2)));

% Source positions in the same coordinates, before participant shifts
src.handRight       = [0.40 0.10];    % right motor hand area (left hand)
src.somatoRight     = [0.40 0.02];    % right somatosensory hand area, behind the motor one
src.handLeft        = [-0.40 0.10];
src.somatoLeft      = [-0.40 0.02];
src.sma             = [0.00 0.25];
src.frontLeft       = [-0.25 0.55];
src.frontRight      = [0.25 0.55];
src.frontoTempRight = [0.55 0.30];    % right fronto-temporal (insula, auditory)
src.audLeft         = [-0.75 0.00];
src.audRight        = [0.75 0.00];
src.visual          = [0.00 -0.75];
src.parietal        = [0.00 -0.50];
fingerAxis          = [-0.8 -0.6];    % thumb to little finger: medially and backwards (right hemisphere)
toneAxisLeft        = [1 0];          % low to high tones: medially (deeper) in auditory cortex
toneAxisRight       = [-1 0];

normalizePeak = @(pattern) pattern / max(abs(pattern));

%% Time axis, waveforms and envelopes
time     = (round(epochTime(1) * fsample):round(epochTime(2) * fsample)) / fsample;
nSamples = numel(time);
fprintf('About %.0f MB per participant, %.1f GB in total\n', ...
        6 * nTrialsPerClass * nChannels * nSamples * 8 / 1e6, ...
        6 * nTrialsPerClass * nChannels * nSamples * 8 * nSubjects / 1e9);

gauss   = @(peak, width) exp(-(time - peak).^2 / (2 * width^2));
sigmoid = @(center, slope) 1 ./ (1 + exp(-(time - center) / slope));

% Evoked responses to an event at time onset
visualShape = @(onset) gauss(onset + 0.10, 0.025) - 0.6 * gauss(onset + 0.17, 0.03);
motorShape  = @(onset) gauss(onset - 0.05, 0.035);   % motor command, before the press
sefShape    = @(onset) gauss(onset + 0.03, 0.012) - 0.5 * gauss(onset + 0.065, 0.02);
n1mShape    = @(onset) gauss(onset + 0.15, 0.03) - 0.5 * gauss(onset + 0.26, 0.045);
visualWave  = visualShape(cueOnset) + visualShape(cueOffset);

% Slow potentials: rise linearly into the first press or tone, then hold to
% the end of the epoch. The early ones (midline readiness potential, SPN)
% start 1.5 s before, the lateralized readiness potential at the go signal.
earlyRamp = min(max((time + 1.5) / 1.5, 0), 1);
lateRamp  = min(max((time - cueOffset) / -cueOffset, 0), 1);

% Envelopes of the induced activity
thetaEnvelope = sigmoid(cueOnset + 0.3, 0.1);                         % until the end of the epoch

% Desynchronization of the sensorimotor rhythms, from 0 (none) to 1 (full);
% motorDesync and listenDesync set its depth per condition. Planning and
% pressing: 70% from shortly after cue onset, full from the go signal into
% the first press. Listening: from the first tone; for mu, partly already from
% the go signal, in anticipation of the tones.
motorDesyncShape    = 0.7 * sigmoid(cueOnset + 0.4, 0.15) + 0.3 * sigmoid(itemTimes(1) - 0.2, 0.08);
listenDesyncShape   = sigmoid(itemTimes(1) + 0.1, 0.05);
listenDesyncShapeMu = 0.4 * sigmoid(cueOffset + 0.25, 0.08) + 0.6 * listenDesyncShape;

%% Participants
participantTruth = struct([]);
for iSubject = 1:nSubjects
    rng(randomSeed + iSubject);

    % --- Participant variability ---------------------------------------------
    headShift  = 0.03 * randn(1, 2);       % moves every source a little
    gainOf     = @() exp(0.3 * randn);     % effect strength, median 1
    effectGain = struct('visual', gainOf(), 'theta', gainOf(), 'motor', gainOf(), 'mu', gainOf(), ...
                        'beta', gainOf(), 'ramp', gainOf(), 'sef', gainOf(), 'n1m', gainOf());
    noiseScale = exp(0.15 * randn);
    at         = @(name) src.(name) + headShift;

    % --- Sensor patterns -----------------------------------------------------
    % Fingers are ordered along the central sulcus and tones along auditory
    % cortex; neighbouring fingers or tones shift the cluster by a small step.
    % Motor commands sit over the motor hand area, press responses over the
    % somatosensory hand area behind it.
    motorPattern = cell(1, 5);
    pressPattern = cell(1, 5);
    tonePattern  = cell(1, 5);
    for iFinger = 1:5
        stepsFromMiddle = iFinger - 3;
        motorPattern{iFinger} = clusterPattern(sensorPos, ...
            at('handRight') + stepsFromMiddle * 0.07 * fingerAxis, 0.15);
        pressPattern{iFinger} = clusterPattern(sensorPos, ...
            at('somatoRight') + stepsFromMiddle * 0.07 * fingerAxis, 0.15);
        tonePattern{iFinger} = normalizePeak( ...
            clusterPattern(sensorPos, at('audLeft') + stepsFromMiddle * 0.06 * toneAxisLeft, 0.15) + ...
            clusterPattern(sensorPos, at('audRight') + stepsFromMiddle * 0.06 * toneAxisRight, 0.15));
    end

    visualPattern = {clusterPattern(sensorPos, at('visual') + [-0.06 0], 0.2), ...
                     0.9 * clusterPattern(sensorPos, at('visual') + [0.06 0.015], 0.2)};
    frontLeft     = clusterPattern(sensorPos, at('frontLeft'), 0.25);
    frontRight    = clusterPattern(sensorPos, at('frontRight'), 0.25);
    thetaPattern  = {normalizePeak(frontLeft + 0.35 * frontRight), ...
                     normalizePeak(0.35 * frontLeft + frontRight)};

    rampLateralPattern = clusterPattern(sensorPos, at('handRight'), 0.2);
    rampMidlinePattern = clusterPattern(sensorPos, at('sma'), 0.25);
    spnPattern         = normalizePeak(clusterPattern(sensorPos, at('frontoTempRight'), 0.25) + ...
                                       0.5 * clusterPattern(sensorPos, at('audLeft'), 0.2));
    muPattern    = {clusterPattern(sensorPos, at('somatoRight'), 0.2), ...
                    0.6 * clusterPattern(sensorPos, at('somatoLeft'), 0.2)};
    betaPattern  = {clusterPattern(sensorPos, at('handRight'), 0.2), ...
                    0.6 * clusterPattern(sensorPos, at('handLeft'), 0.2)};
    alphaPattern = clusterPattern(sensorPos, at('parietal'), 0.3);

    % --- Trials --------------------------------------------------------------
    [conditionGrid, sequenceGrid] = ndgrid(1:numel(conditionCodes), 1:2);
    classCondition = repelem(conditionGrid(:), nTrialsPerClass);
    classSequence  = repelem(sequenceGrid(:), nTrialsPerClass);
    trialOrder     = randperm(numel(classCondition));
    trialCond      = classCondition(trialOrder);   % row of the condition table
    trialSeq       = classSequence(trialOrder);    % sequence 1 or 2
    nTrials        = numel(trialCond);

    jitter    = @() exp(0.25 * randn(nTrials, 1));   % trial-to-trial amplitude
    bandNoise = @(band) filteredNoise(nSamples, nTrials, fsample, band, 0).';

    % --- Activity of interest ------------------------------------------------
    data = zeros(nTrials, nChannels, nSamples);   % fT
    for iSeq = 1:2
        inSeq   = double(trialSeq == iSeq);
        fingers = sequenceFingers(iSeq, :);

        % Visual responses to the cue (the same fractal in every condition)
        data = addComponent(data, amp.visual * effectGain.visual * visualPattern{iSeq}, ...
                            (inSeq .* jitter()) .* visualWave);

        % Sequence-specific frontal theta, weaker in PassiveTones
        data = addComponent(data, amp.theta * effectGain.theta * thetaPattern{iSeq}, ...
                            (inSeq .* thetaGain(trialCond)) .* thetaEnvelope .* bandNoise([4 7]));

        % The two presses (motor command, then somatosensory response) and tones
        pressed = inSeq .* hasPresses(trialCond);
        heard   = inSeq .* toneGain(trialCond);
        for iItem = 1:2
            data = addComponent(data, amp.motor * effectGain.motor * motorPattern{fingers(iItem)}, ...
                                (pressed .* jitter()) .* motorShape(itemTimes(iItem)));
            data = addComponent(data, amp.sef * effectGain.sef * pressPattern{fingers(iItem)}, ...
                                (pressed .* jitter()) .* sefShape(itemTimes(iItem)));
            data = addComponent(data, amp.n1m * effectGain.n1m * tonePattern{fingers(iItem)}, ...
                                (heard .* jitter()) .* n1mShape(itemTimes(iItem)));
        end
    end

    % Mu and beta desynchronization over both hemispheres (not sequence-specific)
    muEnvelope   = 1 - motorDesync(trialCond) .* motorDesyncShape ...
                     - listenDesync(trialCond) .* listenDesyncShapeMu;
    betaEnvelope = 1 - motorDesync(trialCond) .* motorDesyncShape ...
                     - listenDesync(trialCond) .* listenDesyncShape;
    for iSide = 1:2
        data = addComponent(data, amp.mu * effectGain.mu * muPattern{iSide}, ...
                            muEnvelope .* bandNoise([8 12]));
        data = addComponent(data, amp.beta * effectGain.beta * betaPattern{iSide}, ...
                            betaEnvelope .* bandNoise([17 23]));
    end

    % Slow potentials (not sequence-specific): readiness potential before the
    % presses, stimulus-preceding negativity before the tones of the tone task
    data = addComponent(data, amp.rampLateral * effectGain.ramp * rampLateralPattern, ...
                        (rampLateralGain(trialCond) .* jitter()) .* lateRamp);
    data = addComponent(data, amp.rampMidline * effectGain.ramp * rampMidlinePattern, ...
                        (rampMidlineGain(trialCond) .* jitter()) .* earlyRamp);
    data = addComponent(data, amp.spn * effectGain.ramp * spnPattern, ...
                        (spnGain(trialCond) .* jitter()) .* earlyRamp);

    % --- Background ----------------------------------------------------------
    % Posterior alpha whose amplitude changes from trial to trial
    data = addComponent(data, amp.alpha * alphaPattern, exp(0.4 * randn(nTrials, 1)) .* bandNoise([8 12]));

    % Spatially correlated 1/f activity from randomly placed sources
    brainNoise = zeros(nTrials, nChannels, nSamples);
    for iSource = 1:nNoiseSources
        radius      = 0.9 * sqrt(rand);
        direction   = 360 * rand;
        noiseCenter = radius * [cosd(direction) sind(direction)];
        noisePattern = clusterPattern(sensorPos, noiseCenter, 0.15 + 0.2 * rand);
        brainNoise  = addComponent(brainNoise, noisePattern, ...
                                   filteredNoise(nSamples, nTrials, fsample, [0.5 45], 1).');
    end
    brainNoise = brainNoise * (amp.brainNoise * noiseScale / sqrt(mean(brainNoise(:).^2)));
    data = data + brainNoise + amp.sensorNoise * randn(size(data));
    clear brainNoise;

    % --- Save in the format of the real data ---------------------------------
    tmlck_test            = [];
    tmlck_test.label      = channelLabels;
    tmlck_test.time       = time;
    tmlck_test.trial      = data * 1e-15;                 % fT to T
    tmlck_test.dimord     = 'rpt_chan_time';
    tmlck_test.sampleinfo = [(0:nTrials - 1)' * nSamples + 1, (1:nTrials)' * nSamples];
    tmlck_test.grad       = grad;                         % sensor definition, all 306 channels

    trialinfo_test           = struct();
    trialinfo_test.labelCond = conditionCodes(trialCond);
    trialinfo_test.labelSeq  = trialSeq;

    save(fullfile(outputDir, sprintf('epocheddata_S%d.mat', iSubject)), ...
         'tmlck_test', 'trialinfo_test', '-v7');
    clear data tmlck_test;

    participantTruth(iSubject).headShift     = headShift;
    participantTruth(iSubject).effectGain    = effectGain;
    participantTruth(iSubject).noiseScale    = noiseScale;
    participantTruth(iSubject).motorPattern  = motorPattern;
    participantTruth(iSubject).pressPattern  = pressPattern;
    participantTruth(iSubject).tonePattern   = tonePattern;
    participantTruth(iSubject).visualPattern = visualPattern;
    participantTruth(iSubject).thetaPattern  = thetaPattern;
    fprintf('Participant %d: %d trials written\n', iSubject, nTrials);
end

%% Ground truth
truth                 = struct();
truth.events          = struct('cueOnset', cueOnset, 'cueOffset', cueOffset, 'itemTimes', itemTimes);
truth.conditionCodes  = conditionCodes;
truth.conditionNames  = conditionNames;
truth.sequenceFingers = sequenceFingers;
truth.amplitudes      = amp;
truth.conditionGains  = struct('hasPresses', hasPresses, 'toneGain', toneGain, 'thetaGain', thetaGain, ...
                               'motorDesync', motorDesync, 'listenDesync', listenDesync, ...
                               'rampLateralGain', rampLateralGain, ...
                               'rampMidlineGain', rampMidlineGain, 'spnGain', spnGain);
truth.channelLabels   = channelLabels;
truth.sensorPos       = sensorPos;
truth.participants    = participantTruth;
save(fullfile(outputDir, 'simulation_truth.mat'), 'truth');
fprintf('Simulation finished: %s\n', outputDir);

%% Local functions
function pattern = clusterPattern(sensorPos, center, width)
    % Sensor pattern of a source: a cluster of sensors around its position,
    % with weights that fall off with distance (Gaussian of the given width).
    % Peak 1.
    distanceSquared = sum((sensorPos - center).^2, 2);
    pattern = exp(-distanceSquared / (2 * width^2));
    pattern = pattern / max(pattern);
end

function noise = filteredNoise(nSamples, nTrials, fsample, band, exponent)
    % Gaussian noise (samples x trials) whose power falls as 1/f^exponent
    % inside band (Hz) and is zero outside it, scaled to unit RMS.
    % exponent = 0 gives flat band-limited noise.
    freq      = (0:nSamples - 1)' * fsample / nSamples;
    freq      = min(freq, fsample - freq);                     % two-sided spectrum
    gainCurve = (freq >= band(1) & freq <= band(2)) ./ max(freq, eps).^(exponent / 2);
    noise     = real(ifft(fft(randn(nSamples, nTrials)) .* gainCurve));
    noise     = noise / sqrt(mean(noise(:).^2));
end

function data = addComponent(data, pattern, timecourse)
    % Add a source with a fixed sensor pattern (channels x 1) and one time
    % course per trial (trials x samples) to data (trials x channels x samples).
    data = data + reshape(timecourse, size(timecourse, 1), 1, []) .* reshape(pattern, 1, [], 1);
end
