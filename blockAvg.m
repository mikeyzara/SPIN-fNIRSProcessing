function blkAvg = blockAvg(dc4, stim, SDpairs)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

%% Get the indices for the conditions of interest (cond1 - cond4)
cnt = 1; %A counting variable
for j = 1:size(stim,2)
    if contains(stim(j).name, ['condition',num2str(cnt)]) %Check the name of stim to see if it's a condition
        idx(cnt) = j; %Store the index
        cnt = cnt + 1; %Increment the counting variable
    end
end

%% Get the stim marker timings
Hb = dc4.GetDataTimeSeries; %The Hb data
hbTime = dc4.GetTime; %The time vector
cnt2 = 1; %counting variable

Fs = 10; %Sampling Freq
w1 = 10; %Window of analysis
baseWindow = 5; %Window for baseline
cnt2 = 1; %A counting variable

for jj = idx %Loop through all four conditions
    stimTime = find(ismember(hbTime,stim(jj).data(:,1))); % Find location of stim markers in the time vector
    for kk = 1:size(stimTime,1) %Loop through all 25 stim markers in current condition
        sig{cnt2,kk} = Hb(stimTime(kk):stimTime(kk)+w1*Fs,:); %Get the signal from t=markerStim to 15 sec after
        baseline{cnt2,kk} = mean(Hb(stimTime(kk)-baseWindow*Fs:stimTime(kk),:),1); %Get the baseline period (t=-5sec to t=markerStim)
        sig_baseCorr{cnt2,kk} = sig{cnt2,kk}-baseline{cnt2,kk}; %Baseline correction
    end
    cnt2=cnt2+1; %increment the counting variable
end

%% Get HbO and HbR data from sig_baseCorr
for m = 1:size(sig_baseCorr,1) %4 rows -- represents Conditions
    for mm = 1:size(sig_baseCorr,2) %25 columns -- represents each trial
        HbO{m,mm} = sig_baseCorr{m,mm}(:,1:3:end); %The HbO signal
        HbR{m,mm} = sig_baseCorr{m,mm}(:,2:3:end); %The HbR signal
    end
end

%% Get the average response across blocks and exclude outliers
% *New:
meanHbOCh = cellfun(@mean,HbO,'UniformOutput',false); %Calculate the mean HbO per channel for every trial for each condition
for cond = 1:size(HbO,1)
    for trial = 1:size(HbO,2)
        block(:,:,trial) = HbO{cond,trial}; %Store the trials as a 3D matrix -- size [151 18 25] corresponds to [timePoints channels trials]
        
        %Detect trial outliers -- Step 1
%         tempTrials(:,:,trial) = meanHbOCh{cond,trial}; %Temporary variable to hold the mean HbO per ch for every trial

    end
    blkAvg{cond} = mean(block,3); %Calculate the grand average mean of the block

    % Detect Trial outliers -- Step 2
%     meanHbO_Trials = mean(tempTrials,3); %Calculate the mean HbO across trials for each channel
%     stdevHbO = std(tempTrials,0,3); %Calculate the std across trials for each channel

%     outliers = isoutlier(tempTrials,'mean',3)

%     trialOutlier{cond} = trialO;
end

end