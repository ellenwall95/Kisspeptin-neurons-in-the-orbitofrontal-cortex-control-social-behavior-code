% code for reading .doric file continous mode
% code for analysis of ORB photometry signal

clear all
close all
clc

%% Please adjust parameters
FilePath = ['/Users/ellenwall/Documents/PostDoc/Photometry/OFCphotometry'];%i.e. where is the file
Filename = {'5314_female_0006_104247294_incomplete.doric'}; %name of the file

excludeTimeStart = 0; % Define the time to exclude in seconds (adjust as needed)
excludeTimeEnd = 0; % Define the time to exclude at the end in seconds (adjust as needed)

%% extract data
currFile = fullfile(FilePath,Filename{1});
Data_Acquired = ExtractDataAcquisition(currFile);
dataTable = [Data_Acquired(4).Data(2).Data,Data_Acquired(3).Data(1).Data,Data_Acquired(4).Data(1).Data];
recordingEndTime = max(dataTable(:,1)) - excludeTimeEnd;
dataTable = dataTable(dataTable(:,1) >= excludeTimeStart & dataTable(:,1) <= recordingEndTime, :); % Filter the data to exclude rows before this time

Time = dataTable(:,1); %Time
GCaMP = dataTable(:,2); %GCaMP
Background = dataTable (:,3); %405 background

%% remove artifacts (e.g packet loss)
% Compute first derivatives (rate of change)
GCaMP_diff = diff(GCaMP);
Background_diff = diff(Background);

% Threshold multipliers used 
GCaMP_threshold_multiplier = 2;
Background_threshold_multiplier = 2;

% Set a threshold for artifact detection (adjust as needed)
threshold_GCaMP = GCaMP_threshold_multiplier * std(GCaMP_diff);
threshold_Background = Background_threshold_multiplier * std(Background_diff);

% Find artifact indices where the change is above the threshold
artifact_idx_GCaMP = find(abs(GCaMP_diff) > threshold_GCaMP);
artifact_idx_Background = find(abs(Background_diff) > threshold_Background);

% Identify artifacts present in BOTH signals
artifact_idx = intersect(artifact_idx_GCaMP, artifact_idx_Background);

% Expand artifact indices to include surrounding points (optional)
expand_range = 2; % Adjust if needed
artifact_idx_expanded = [];
for i = 1:length(artifact_idx)
    artifact_idx_expanded = [artifact_idx_expanded, (artifact_idx(i)-expand_range):(artifact_idx(i)+expand_range)];
end
artifact_idx_expanded = unique(artifact_idx_expanded);
artifact_idx_expanded(artifact_idx_expanded < 1 | artifact_idx_expanded > length(Time)) = []; % Keep valid indices

% Convert indices to time
artifact_times = Time(artifact_idx_expanded);

% Create cleaned versions of signals
GCaMP_cleaned = GCaMP;
Background_cleaned = Background;

% Set artifact regions to NaN for proper interpolation
GCaMP_cleaned(artifact_idx_expanded) = NaN;
Background_cleaned(artifact_idx_expanded) = NaN;

% Interpolate over artifact regions (linear interpolation)
GCaMP_cleaned = fillmissing(GCaMP_cleaned, 'linear');
Background_cleaned = fillmissing(Background_cleaned, 'linear');

%Plot signals with artifacts highlighted and cleaned signal
figure;

% GCaMP with artifacts
subplot(2,1,1);
plot(Time, GCaMP, 'b'); hold on;
plot(Time(artifact_idx_expanded), GCaMP(artifact_idx_expanded), 'ro', 'MarkerSize', 4, 'DisplayName', 'Artifacts');
plot(Time, GCaMP_cleaned, 'g');
legend('Raw GCaMP', 'Artifacts', 'Cleaned GCaMP');
title('GCaMP Signal - Artifact Detection and Cleaning');
xlabel('Time (s)');
ylabel('Signal');

% Background with artifacts
subplot(2,1,2);
plot(Time, Background, 'k'); hold on;
plot(Time(artifact_idx_expanded), Background(artifact_idx_expanded), 'mo', 'MarkerSize', 4, 'DisplayName', 'Artifacts');
plot(Time, Background_cleaned, 'g');
legend('Raw 405', 'Artifacts', 'Cleaned 405');
title('Background Signal - Artifact Detection and Cleaning');
xlabel('Time (s)');
ylabel('Signal');


%% Low band pass filter 
Fs = 1 / mean(diff(Time)); % Sampling rate in Hz
GCaMP_filt = lowpass(GCaMP_cleaned, 3, Fs, 'Steepness', 0.95); %low-pass filter on 465
Background_filt = lowpass(Background_cleaned, 3, Fs, 'Steepness', 0.95); %low-pass filter on 405

% Trim start and end to avoid edge artifacts
trimSec = 2; % Adjust if needed
trimSamples = round(trimSec * Fs);
validIdx = (1+trimSamples):(length(Time)-trimSamples);

Time_trimmed = Time(validIdx);
GCaMP_filt_trimmed = GCaMP_filt(validIdx);
Background_filt_trimmed = Background_filt(validIdx);

figure;
subplot(2,1,1);
plot(Time, GCaMP_cleaned, 'b'); hold on;
plot(Time_trimmed, GCaMP_filt_trimmed, 'r');
legend('Raw GCaMP','Filtered GCaMP');
title('GCaMP Signal');

subplot(2,1,2);
plot(Time, Background_cleaned, 'k'); hold on;
plot(Time_trimmed, Background_filt_trimmed, 'g');
legend('Raw 405','Filtered 405');
title('Background Signal');

%% IRLS regression
IRLS_constant = 3;
[~, fitted_iso] = IRLS_dFF(GCaMP_filt_trimmed, Background_filt_trimmed, IRLS_constant);

% Plot overlay
figure;
plot(Time_trimmed, GCaMP_filt_trimmed, 'b', 'DisplayName', 'Filtered 465 (GCaMP)');
hold on;
plot(Time_trimmed, fitted_iso, 'r', 'DisplayName', 'Fitted 405 (IRLS Background)');
legend;
xlabel('Time (s)');
ylabel('Signal (a.u.)');
title('Filtered 465 Signal and IRLS-Fitted 405 Overlay');

[dFF, fitted_iso] = IRLS_dFF(GCaMP_filt_trimmed, Background_filt_trimmed, IRLS_constant);

GCaMP_minus_IRLS = GCaMP_filt_trimmed - fitted_iso;
zScore_GCaMP_minus_IRLS = (GCaMP_minus_IRLS - mean(GCaMP_minus_IRLS)) / std(GCaMP_minus_IRLS);
zScore_dFF = (dFF - mean(dFF)) / std(dFF);

figure;

% Subplot 1: dF/F
subplot(3,1,1);
plot(Time_trimmed, dFF, 'm');
xlabel('Time (s)');
ylabel('dF/F');
title('IRLS-Corrected ΔF/F');

% Subplot 2: z-score of 465 - fitted 405
subplot(3,1,2);
plot(Time_trimmed, zScore_GCaMP_minus_IRLS, 'b');
xlabel('Time (s)');
ylabel('zScore');
title('Z-Score of (465 - Fitted 405)');

% Subplot 3: z-score of dF/F
subplot(3,1,3);
plot(Time_trimmed, zScore_dFF, 'g');
xlabel('Time (s)');
ylabel('zScored dF/F');
title('Z-Score of dF/F');

sgtitle('Photometry Signal Comparison (IRLS-Corrected)');

%% Table of data

GCaMP_threshold_column = repmat(GCaMP_threshold_multiplier, length(Time_trimmed), 1);
Background_threshold_column = repmat(Background_threshold_multiplier, length(Time_trimmed), 1);

% Create table
data_table_photometry = table( ...
    Time_trimmed, ...
    GCaMP_filt_trimmed, ...
    Background_filt_trimmed, ...
    fitted_iso, ...
    GCaMP_minus_IRLS, ...
    dFF, ...
    zScore_GCaMP_minus_IRLS, ...
    zScore_dFF, ...
    GCaMP_threshold_column, ...
    Background_threshold_column, ...
    'VariableNames', {'Time_s', 'GCaMP_465', 'Background_405', 'Fitted_405_IRLS','GCaMPminusfitted405', 'dFF', 'zscoreGCaMPminusIRLS405', 'zscoredFF', 'GCaMP_Threshold_Multiplier', 'Background_Threshold_Multiplier'}...
);

writetable(data_table_photometry, [Filename{1}(1:end-6) '_IRLS_zscoredFF.xlsx']);